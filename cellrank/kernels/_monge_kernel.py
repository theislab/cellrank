from typing import Any, Tuple, Union, Literal, Callable, Optional, Sequence

from moscot.problems.time import TemporalNeuralProblem

from anndata import AnnData
from cellrank import logging as logg
from cellrank.kernels import VelocityKernel
from cellrank._utils._docs import d, inject_docs
from cellrank._utils._enum import _DEFAULT_BACKEND, Backend_t
from cellrank.kernels.utils import Deterministic
from cellrank.kernels.mixins import ConnectivityMixin
from scvelo.preprocessing.moments import get_moments
from cellrank.kernels._base_kernel import BidirectionalKernel
from cellrank.kernels.utils._similarity import Similarity, SimilarityABC
from cellrank.kernels.utils._velocity_model import BackwardMode, VelocityModel

import numpy as np
from scipy.sparse import issparse

__all__ = ["MongeKernel"]

_ATTR_OPTIONS = Literal["obsm", "layers"]


@d.dedent
class MongeKernel(VelocityKernel):
    """
    Summary.

    Description.

    Args:
        VelocityKernel (_type_): _description_
    """

    def __init__(
        self,
        tnp: TemporalNeuralProblem,
        backward: bool = False,
        attr: _ATTR_OPTIONS = "obsm",
        xkey: Optional[str] = "X_pca",
        vkey: Optional[str] = "velocity",
        **kwargs: Any,
    ):
        super().__init__(
            tnp.adata,
            backward=backward,
            xkey=xkey,
            vkey=vkey,
            attr=attr,  # this will be passed to _read_from_adata()
            **kwargs,
        )
        self._logits: Optional[np.ndarray] = None

    def _read_from_adata(
        self,
        attr: _ATTR_OPTIONS = "obsm",
        xkey: Optional[str] = "X_pca",
        vkey: Optional[str] = "MongeVelocities",
        gene_subset: Optional[Union[str, Sequence[str]]] = None,
        **kwargs: Any,
    ) -> None:
        super(VelocityKernel, self)._read_from_adata(
            **kwargs
        )  # TODO: Is there a better way?

        if attr == "layers":
            self._xdata = self._extract_layer(xkey, subset=gene_subset)
            self._vdata = self._extract_layer(vkey, subset=gene_subset)
        elif attr == "obsm":
            self._xdata = self.adata.obsm[xkey]
            self._vdata = np.array(
                self.adata.obsm[vkey]
            )  # TODO: I dont think we need self._vdata to be of type
            # 'jaxlib.xla_extension.ArrayImpl', therefore I am converting to numpy array here. But not sure if this is
            # the best way of doing so? Could velocities also be sparse? (jaxlib.xla_extension.ArrayImpl type causes
            # errors in compute_transition_matrix() -> ModelABC.__init__() as jax' astype does not have keyword `copy`:
            # Signature (dtype: Union[Any, str, numpy.dtype, jax._src.SupportsDType]) -> jax.Array
        else:
            raise ValueError(f"Invalid attr {attr}. Valid choices are {_ATTR_OPTIONS}")

        nans = np.isnan(np.sum(self._vdata, axis=0))
        if np.any(nans):
            self._xdata = self._xdata[:, ~nans]
            self._vdata = self._vdata[:, ~nans]

        # fmt: off
        self._vexp = get_moments(self.adata, self._vdata, second_order=False).astype(np.float64, copy=False)
        self._vvar = get_moments(self.adata, self._vdata, second_order=True).astype(np.float64, copy=False)
        # fmt: on

    def compute_transition_matrix(
        self,
        backward_mode: Literal["transpose", "negate"] = BackwardMode.TRANSPOSE,
        similarity: Union[
            Literal["correlation", "cosine", "dot_product"],
            Callable[[np.ndarray, np.ndarray, float], Tuple[np.ndarray, np.ndarray]],
        ] = Similarity.CORRELATION,
        softmax_scale: Optional[float] = None,
        seed: Optional[
            int
        ] = None,  # TODO: does the seed impact results for the deterministic model?
        **kwargs: Any,
    ) -> "VelocityKernel":
        """
        Compute transition matrix based on velocity directions on the local manifold.

        For each cell, infer transition probabilities based on the cell's velocity-extrapolated cell state and the
        cell states of its *K* nearest neighbors.

        Parameters
        ----------
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
        start = logg.info("Computing transition matrix.")

        # fmt: off
        params = {"model": "deterministic", "similarity": str(similarity), "softmax_scale": softmax_scale}
        if self.backward:
            params["bwd_mode"] = str(backward_mode)

        if self._reuse_cache(params, time=start):
            return self

        if softmax_scale is None:
            softmax_scale = self._estimate_softmax_scale(backward_mode=backward_mode, similarity=similarity)
            logg.info(f"Using `softmax_scale={softmax_scale:.4f}`")
            params["softmax_scale"] = softmax_scale
        # fmt: on

        model = self._create_model(
            params["model"],
            backward_mode=backward_mode,
            similarity=similarity,
            softmax_scale=softmax_scale,
            n_samples=1000,  # will be ignored for Deterministic model
            seed=seed,
        )

        self.transition_matrix, self._logits = model(**kwargs)

        logg.info("    Finish", time=start)

        return self
