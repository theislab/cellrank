from typing import Any, Union, Callable, Optional, Sequence
from typing_extensions import Literal

from abc import ABC

from anndata import AnnData
from cellrank import logging as logg
from cellrank.tl._enum import _DEFAULT_BACKEND, Backend_t
from cellrank.ul._docs import d, inject_docs
from cellrank.tl.kernels.utils import MonteCarlo, Stochastic, Deterministic
from cellrank.tl.kernels._mixins import ConnectivityMixin
from scvelo.preprocessing.moments import get_moments
from cellrank.tl.kernels._base_kernel import BidirectionalKernel
from cellrank.tl.kernels.utils._velocity_model import BackwardMode, VelocityModel
from cellrank.tl.kernels.utils._similarity_scheme import Scheme, Similarity

import numpy as np
from scipy.sparse import issparse

__all__ = ("DisplacementKernel",)


# TODO: docs
class DisplacementKernel(ConnectivityMixin, BidirectionalKernel, ABC):
    def __init__(
        self,
        adata: AnnData,
        backward: bool = False,
        xkey: Optional[str] = "Ms",
        vkey: Optional[str] = "velocity",
        **kwargs: Any,
    ):
        super().__init__(
            adata,
            backward=backward,
            xkey=xkey,
            vkey=vkey,
            **kwargs,
        )

    def _read_from_adata(
        self,
        xkey: str = "Ms",
        vkey: str = "velocity",
        subset: Optional[Union[str, Sequence[str]]] = None,
        **kwargs: Any,
    ) -> None:
        super()._read_from_adata(**kwargs)

        if subset is None and f"{vkey}_genes" in self.adata.var:
            subset = self.adata.var[f"{vkey}_genes"]

        self._xdata = self._extract_layer(xkey, subset=subset)
        self._vdata = self._extract_layer(vkey, subset=subset)

        # TODO(Marius1311): why not nans = np.any(np.isnan(self._vdata), axis=0)?
        nans = np.isnan(np.sum(self._vdata, axis=0))
        if np.any(nans):
            self._xdata = self._xdata[:, ~nans]
            self._vdata = self._vdata[:, ~nans]
        # fmt: off
        self._vexp = get_moments(self.adata, self._vdata, second_order=False).astype(np.float64, copy=False)
        self._vvar = get_moments(self.adata, self._vdata, second_order=True).astype(np.float64, copy=False)
        # fmt: on

    # TODO(michalk8): check callable signature
    @inject_docs(m=VelocityModel, b=BackwardMode, s=Scheme)  # don't swap the order
    @d.dedent
    def compute_transition_matrix(
        self,
        model: Literal[
            "deterministic", "stochastic", "monte_carlo"
        ] = VelocityModel.DETERMINISTIC,
        backward_mode: Literal["transpose", "negate"] = BackwardMode.TRANSPOSE,
        similarity: Union[
            Literal["correlation", "cosine", "dot_product"],
            Callable[[np.ndarray], np.ndarray],
        ] = Scheme.CORRELATION,
        softmax_scale: Optional[float] = None,
        n_samples: int = 1000,
        seed: Optional[int] = None,
        **kwargs: Any,
    ) -> "DisplacementKernel":
        """
        Compute transition matrix based on velocity directions on the local manifold.

        For each cell, infer transition probabilities based on the cell's velocity-extrapolated cell state and the
        cell states of its *K* nearest neighbors.

        Parameters
        ----------
        %(velocity_mode)s
        %(velocity_backward_mode)s
        %(softmax_scale)s
        %(velocity_scheme)s
        n_samples
            Number of bootstrap samples when ``mode = {m.MONTE_CARLO!r}``.
        seed
            Random seed when ``mode = {m.MONTE_CARLO!r}``.
        %(parallel)s

        Returns
        -------
        Self and updates :attr:`transition_matrix` and :attr:`params`.
        """
        start = logg.info(f"Computing transition matrix using `{model!r}` model")

        # fmt: off
        params = {"model": model, "similarity": str(similarity), "softmax_scale": softmax_scale}
        if self.backward:
            params["bwd_mode"] = str(backward_mode)
        if VelocityModel(model) == VelocityModel.MONTE_CARLO:
            params["n_samples"] = n_samples
            params["seed"] = seed
        if self._reuse_cache(params, time=start):
            return self

        if softmax_scale is None:
            softmax_scale = self._estimate_softmax_scale(backward_mode=backward_mode, similarity=similarity)
            logg.info(f"Using `softmax_scale={softmax_scale:.4f}`")
            params["softmax_scale"] = softmax_scale
        # fmt: on

        model = self._create_model(
            model,
            backward_mode=backward_mode,
            similarity=similarity,
            softmax_scale=softmax_scale,
            n_samples=n_samples,
            seed=seed,
        )
        if isinstance(model, Stochastic):
            kwargs["backend"] = _DEFAULT_BACKEND
        self.transition_matrix, _ = model(**kwargs)

        logg.info("    Finish", time=start)

        return self

    def _create_model(
        self,
        model: Union[str, VelocityModel],
        backward_mode: Literal["negate", "transpose"],
        similarity: Union[str, Scheme, Callable],
        n_samples: int = 1000,
        seed: Optional[int] = None,
        **kwargs: Any,
    ):
        model = VelocityModel(model)
        backward_mode = BackwardMode(backward_mode) if self.backward else None

        if isinstance(similarity, str):
            similarity = Similarity.create(Scheme(similarity))
        elif not callable(similarity):
            raise TypeError(
                f"Expected `scheme` to be a function, found `{type(similarity).__name__}`."
            )

        if self.backward and model != VelocityModel.DETERMINISTIC:
            logg.warning(
                f"Mode `{model!r}` is currently not supported for backward process. "
                f"Using to `mode={VelocityModel.DETERMINISTIC!r}`"
            )
            model = VelocityModel.DETERMINISTIC

        if model == VelocityModel.STOCHASTIC and not hasattr(similarity, "hessian"):
            model = VelocityModel.MONTE_CARLO
            logg.warning(
                f"Unable to detect a method for Hessian computation. If using one of the "
                f"predefined similarity functions, consider installing `jax` as "
                f"`pip install jax jaxlib`. Using `mode={model!r}` and `n_samples={n_samples}`"
            )

        if model == VelocityModel.DETERMINISTIC:
            return Deterministic(
                self._conn,
                self._xdata,
                self._vdata,
                similarity=similarity,
                backward_mode=backward_mode,
                **kwargs,
            )
        if model == VelocityModel.STOCHASTIC:
            return Stochastic(
                self._conn,
                self._xdata,
                self._vexp,
                self._vvar,
                similarity=similarity,
                backward_mode=backward_mode,
                **kwargs,
            )
        if model == VelocityModel.MONTE_CARLO:
            return MonteCarlo(
                self._conn,
                self._xdata,
                self._vexp,
                self._vvar,
                similarity=similarity,
                backward_mode=backward_mode,
                n_samples=n_samples,
                seed=seed,
                **kwargs,
            )

        raise NotImplementedError("TODO")

    def _estimate_softmax_scale(
        self,
        n_jobs: Optional[int] = None,
        backend: Backend_t = _DEFAULT_BACKEND,
        **kwargs,
    ) -> float:
        model = self._create_model(
            VelocityModel.DETERMINISTIC, softmax_scale=1.0, **kwargs
        )
        _, logits = model(n_jobs, backend)
        return 1.0 / np.median(np.abs(logits.data))

    def _extract_layer(
        self,
        key: Optional[str] = None,
        subset: Optional[Union[str, Sequence[str]]] = None,
        dtype: np.dtype = np.float64,
    ) -> np.ndarray:
        if key in (None, "X"):
            data = self.adata.X
        elif key in self.adata.layers:
            data = self.adata.layers[key]
        else:
            raise KeyError("TODO.")

        if isinstance(subset, str):
            subset = self.adata.var[subset]
        subset = np.asarray(subset)
        if np.issubdtype(subset.dtype, bool) and subset.shape == (data.shape[1],):
            data = data[:, subset]
        else:
            data = data[:, self.adata.var_names.isin(subset)]

        data = data.astype(dtype, copy=False)
        return data.toarray() if issparse(data) else data

    def __invert__(self) -> "DisplacementKernel":
        # fmt: off
        dk = self.copy()
        dk._backward = not self.backward
        dk._params = {}
        dk._transition_matrix = None
        return dk
        # fmt: on
