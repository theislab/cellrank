from typing import Any, Callable, Literal, Optional, Sequence, Tuple, Union

from scvelo.preprocessing.moments import get_moments

import numpy as np
import scipy.sparse as sp

from anndata import AnnData

from cellrank import logging as logg
from cellrank._utils._docs import d, inject_docs
from cellrank._utils._enum import DEFAULT_BACKEND, Backend_t
from cellrank.kernels._base_kernel import BidirectionalKernel
from cellrank.kernels.mixins import ConnectivityMixin
from cellrank.kernels.utils import Deterministic, MonteCarlo, Stochastic
from cellrank.kernels.utils._similarity import Similarity, SimilarityABC
from cellrank.kernels.utils._velocity_model import BackwardMode, VelocityModel

__all__ = ["VelocityKernel"]


@d.dedent
class VelocityKernel(ConnectivityMixin, BidirectionalKernel):
    """Kernel which computes a transition matrix based on RNA velocity.

    .. seealso::
        - See :doc:`../../../notebooks/tutorials/kernels/200_rna_velocity` on how to
          compute the :attr:`~cellrank.kernels.VelocityKernel.transition_matrix` based on RNA velocity.

    This borrows ideas from both :cite:`manno:18` and :cite:`bergen:20`. In short, for each cell :math:`i`, we compute
    transition probabilities :math:`T_{i, j}` to each cell :math:`j` in the neighborhood of :math:`i`. We quantify
    how much the velocity vector :math:`v_i` of cell :math:`i` points towards each of its nearst neighbors. For
    this comparison, we support various schemes including cosine similarity and pearson correlation.

    Parameters
    ----------
    %(adata)s
    %(backward)s
    attr
        Attribute of :class:`~anndata.AnnData` to read from.
    xkey
        Key in :attr:`~anndata.AnnData.layers` or :attr:`~anndata.AnnData.obsm`
        where expected gene expression counts are stored.
    vkey
        Key in :attr:`~anndata.AnnData.layers` or :attr:`~anndata.AnnData.obsm` where velocities are stored.
    gene_subset
        List of genes to be used to compute transition probabilities.
        If not specified, genes from :attr:`adata.var['{vkey}_genes'] <anndata.AnnData.var>` are used.
        This feature is only available when reading from :attr:`anndata.AnnData.layers` and will be ignored otherwise.
    kwargs
        Keyword arguments for the :class:`~cellrank.kernels.Kernel`.
    """

    def __init__(
        self,
        adata: AnnData,
        backward: bool = False,
        attr: Optional[Literal["layers", "obsm"]] = "layers",
        xkey: Optional[str] = "Ms",
        vkey: Optional[str] = "velocity",
        **kwargs: Any,
    ):
        super().__init__(
            adata,
            backward=backward,
            xkey=xkey,
            vkey=vkey,
            attr=attr,
            **kwargs,
        )
        self._logits: Optional[np.ndarray] = None

    def _read_from_adata(
        self,
        xkey: Optional[str] = "Ms",
        vkey: Optional[str] = "velocity",
        attr: Optional[Literal["layers", "obsm"]] = "layers",
        gene_subset: Optional[Union[str, Sequence[str]]] = None,
        **kwargs: Any,
    ) -> None:
        super()._read_from_adata(**kwargs)

        if attr == "layers" and gene_subset is None and f"{vkey}_genes" in self.adata.var:
            gene_subset = self.adata.var[f"{vkey}_genes"]
        elif attr == "obsm" and gene_subset is not None:
            logg.warning(
                f"Found `gene_subset != None`, but it is not supported for `adata.{attr}`. Using `gene_subset = None`."
            )
            gene_subset = None

        self._xdata = self._extract_data(key=xkey, attr=attr, subset=gene_subset)
        self._vdata = self._extract_data(key=vkey, attr=attr, subset=gene_subset)
        np.testing.assert_array_equal(
            x=self._xdata.shape,
            y=self._vdata.shape,
            err_msg=f"Shape mismatch: {self._xdata.shape} vs {self._vdata.shape}",
        )

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
        model: Literal["deterministic", "stochastic", "monte_carlo"] = VelocityModel.DETERMINISTIC,
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
        """Compute transition matrix based on velocity directions on the local manifold.

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
        Returns self and updates :attr:`transition_matrix`, :attr:`logits` and :attr:`params`.
        """
        start = logg.info(f"Computing transition matrix using `{model!r}` model")

        # fmt: off
        params = {"model": str(model), "similarity": str(similarity), "softmax_scale": softmax_scale}
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
            kwargs["backend"] = DEFAULT_BACKEND
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
            raise TypeError(f"Expected `scheme` to be a function, found `{type(similarity).__name__}`.")

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
                f"`pip install jax`. Using `mode={model!r}` and `n_samples={n_samples}`"
            )

        if model == VelocityModel.DETERMINISTIC:
            return Deterministic(
                self.connectivities,
                self._xdata,
                self._vdata,
                similarity=similarity,
                backward_mode=backward_mode,
                **kwargs,
            )
        if model == VelocityModel.STOCHASTIC:
            return Stochastic(
                self.connectivities,
                self._xdata,
                self._vexp,
                self._vvar,
                similarity=similarity,
                backward_mode=backward_mode,
                **kwargs,
            )
        if model == VelocityModel.MONTE_CARLO:
            return MonteCarlo(
                self.connectivities,
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
        backend: Backend_t = DEFAULT_BACKEND,
        **kwargs,
    ) -> float:
        model = self._create_model(VelocityModel.DETERMINISTIC, softmax_scale=1.0, **kwargs)
        _, logits = model(n_jobs, backend)
        return 1.0 / np.median(np.abs(logits.data))

    def _extract_data(
        self,
        key: Optional[str] = None,
        attr: Optional[Union[None, Literal["layers", "obsm"]]] = None,
        subset: Optional[Union[str, Sequence[str]]] = None,
        dtype: np.dtype = np.float64,
    ) -> np.ndarray:
        if key in (None, "X"):
            data = self.adata.X
        elif attr == "layers" and key in self.adata.layers:
            data = self.adata.layers[key]
        elif attr == "obsm" and key in self.adata.obsm:
            data = np.asarray(self.adata.obsm[key])
        else:
            raise KeyError(f"Unable to find data in `adata.{attr}[{key!r}]`.")

        if subset is not None:
            if isinstance(subset, str):
                subset = self.adata.var[subset]
            subset = np.asarray(subset)
            if np.issubdtype(subset.dtype, bool) and subset.shape == (data.shape[1],):
                data = data[:, subset]
            else:
                data = data[:, np.isin(self.adata.var_names, subset)]

        data = data.astype(dtype, copy=False)
        return data.toarray() if sp.issparse(data) else data

    @property
    def logits(self) -> Optional[np.ndarray]:
        """Array of shape ``(n_cells, n_cells)`` containing the unnormalized transition matrix."""
        return self._logits

    def __invert__(self) -> "VelocityKernel":
        dk = self._copy_ignore("_transition_matrix", "_logits")
        dk._backward = not self.backward
        dk._params = {}
        return dk
