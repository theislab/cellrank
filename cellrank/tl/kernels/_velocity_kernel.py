from typing import Any, Tuple, Union, Callable, Optional, Sequence
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
from cellrank.tl.kernels.utils._similarity import Similarity, SimilarityABC
from cellrank.tl.kernels.utils._velocity_model import BackwardMode, VelocityModel

import numpy as np
from scipy.sparse import issparse

__all__ = ("VelocityKernel",)


class VelocityKernel(ConnectivityMixin, BidirectionalKernel, ABC):
    """
    Kernel which computes a transition matrix based on RNA velocity.

    This borrows ideas from both :cite:`manno:18` and :cite:`bergen:20`. In short, for each cell *i*, we compute
    transition probabilities :math:`p_{i, j}` to each cell *j* in the neighborhood of *i*. The transition probabilities
    are computed as a multinomial logistic regression where the weights :math:`w_j` (for all *j*) are given
    by the vector that connects cell *i* with cell *j* in gene expression space, and the features :math:`x_i` are given
    by the velocity vector :math:`v_i` of cell *i*.

    Parameters
    ----------
    %(adata)s
    %(backward)s
    xkey
        Key in :attr:`anndata.AnnData.layers` where expected gene expression counts are stored.
    vkey
        Key in :attr:`anndata.AnnData.layers` where velocities are stored.
    gene_subset
        List of genes to be used to compute transition probabilities.
        If not specified, genes from :attr:`anndata.AnnData.var` ``['{vkey}_genes']`` are used.
    kwargs
        Keyword arguments for the parent class.
    """

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
        self._logits: Optional[np.ndarray] = None

    def _read_from_adata(
        self,
        xkey: Optional[str] = "Ms",
        vkey: Optional[str] = "velocity",
        gene_subset: Optional[Union[str, Sequence[str]]] = None,
        **kwargs: Any,
    ) -> None:
        super()._read_from_adata(**kwargs)

        if gene_subset is None and f"{vkey}_genes" in self.adata.var:
            gene_subset = self.adata.var[f"{vkey}_genes"]

        self._xdata = self._extract_layer(xkey, subset=gene_subset)
        self._vdata = self._extract_layer(vkey, subset=gene_subset)

        nans = np.isnan(np.sum(self._vdata, axis=0))
        if np.any(nans):
            self._xdata = self._xdata[:, ~nans]
            self._vdata = self._vdata[:, ~nans]
        # fmt: off
        self._vexp = get_moments(self.adata, self._vdata, second_order=False).astype(np.float64, copy=False)
        self._vvar = get_moments(self.adata, self._vdata, second_order=True).astype(np.float64, copy=False)
        # fmt: on

    # TODO(michalk8): remove the docrep in 2.0
    @inject_docs(m=VelocityModel, b=BackwardMode, s=Similarity)  # don't swap the order
    @d.dedent
    def compute_transition_matrix(
        self,
        model: Literal[
            "deterministic", "stochastic", "monte_carlo"
        ] = VelocityModel.DETERMINISTIC,
        backward_mode: Literal["transpose", "negate"] = BackwardMode.TRANSPOSE,
        similarity: Union[
            Literal["correlation", "cosine", "dot_product"],
            Callable[[np.ndarray, np.ndarray, float], Tuple[np.ndarray, np.ndarray]],
        ] = Similarity.CORRELATION,
        softmax_scale: Optional[float] = None,
        n_samples: int = 1000,
        seed: Optional[int] = None,
        **kwargs: Any,
    ) -> "VelocityKernel":
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
            Number of samples when ``mode = {m.MONTE_CARLO!r}``.
        seed
            Random seed when ``mode = {m.MONTE_CARLO!r}``.
        %(parallel)s
        kwargs
            Keyword arguments for the underlying ``model``.

        Returns
        -------
        Self and updates :attr:`transition_matrix`, :attr:`logits` and :attr:`params`.
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
        self.transition_matrix, self._logits = model(**kwargs)

        logg.info("    Finish", time=start)

        return self

    def _create_model(
        self,
        model: Union[str, VelocityModel],
        backward_mode: Literal["negate", "transpose"],
        similarity: Union[str, Similarity, Callable],
        n_samples: int = 1000,
        seed: Optional[int] = None,
        **kwargs: Any,
    ):
        model = VelocityModel(model)
        backward_mode = BackwardMode(backward_mode) if self.backward else None

        if isinstance(similarity, str):
            similarity = SimilarityABC.create(Similarity(similarity))
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

        raise NotImplementedError(f"Model `{model}` is not yet implemented.")

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
            raise KeyError(f"Unable to find data in `adata.layers[{key}]`.")

        if isinstance(subset, str):
            subset = self.adata.var[subset]
        subset = np.asarray(subset)
        if np.issubdtype(subset.dtype, bool) and subset.shape == (data.shape[1],):
            data = data[:, subset]
        else:
            data = data[:, np.isin(self.adata.var_names, subset)]

        data = data.astype(dtype, copy=False)
        return data.toarray() if issparse(data) else data

    @property
    def logits(self) -> Optional[np.ndarray]:
        """Array of shape ``(n_cells, n_cells)`` containing unnormalized transition matrix."""
        return self._logits

    def __invert__(self) -> "VelocityKernel":
        dk = self._copy_ignore("_transition_matrix", "_logits")
        dk._backward = not self.backward
        dk._params = {}
        return dk
