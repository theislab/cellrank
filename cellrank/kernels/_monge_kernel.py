from typing import Any, Tuple, Union, Literal, Callable, Optional, Sequence

from moscot.problems.time import TemporalNeuralProblem

from cellrank._utils._docs import d, inject_docs
from scvelo.preprocessing.moments import get_moments
from cellrank.kernels._velocity_kernel import VelocityKernel
from cellrank.kernels.utils._similarity import Similarity
from cellrank.kernels.utils._velocity_model import BackwardMode, VelocityModel

import numpy as np

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
        vkey: Optional[str] = "Monge_velocities",
        **kwargs: Any,
    ):
        super().__init__(
            tnp.adata,
            backward=backward,
            xkey=xkey,
            vkey=vkey,
            attr=attr,
            **kwargs,
        )
        self._logits: Optional[np.ndarray] = None

    def _read_from_adata(
        self,
        attr: _ATTR_OPTIONS = "obsm",
        xkey: Optional[str] = "X_pca",
        vkey: Optional[str] = "Monge_velocities",
        gene_subset: Optional[Union[str, Sequence[str]]] = None,
        **kwargs: Any,
    ) -> None:
        super(VelocityKernel, self)._read_from_adata(**kwargs)

        if attr == "layers":
            self._xdata = self._extract_layer(xkey, subset=gene_subset)
            self._vdata = self._extract_layer(vkey, subset=gene_subset)
        elif attr == "obsm":
            self._xdata = self.adata.obsm[xkey]
            self._vdata = np.array(self.adata.obsm[vkey])
            # NOTE: jaxlib.xla_extension.ArrayImpl type causes errors in
            # compute_transition_matrix() -> ModelABC.__init__() as jax' astype does not have keyword `copy`.
            # Converting to numpy array here to avoid this issue.
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
        **kwargs: Any,
    ) -> "MongeKernel":
        """
        Compute transition matrix based on velocity directions on the local manifold.

        For each cell, infer transition probabilities based on the cell's velocity-extrapolated cell state and the
        cell states of its *K* nearest neighbors.

        Parameters
        ----------
        %(velocity_backward_mode)s
        %(softmax_scale)s
        %(velocity_scheme)s
        %(parallel)s
        kwargs
            Keyword arguments for the underlying the Deterministic VelocityModel.

        Returns
        -------
        Self and updates :attr:`transition_matrix`, :attr:`logits` and :attr:`params`.
        """
        super().compute_transition_matrix(
            model=VelocityModel.DETERMINISTIC,
            backward_mode=backward_mode,
            similarity=similarity,
            softmax_scale=softmax_scale,
            **kwargs,
        )
        return self
