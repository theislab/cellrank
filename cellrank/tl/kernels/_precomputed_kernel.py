from typing import Any, Dict, Union, Optional

from anndata import AnnData
from cellrank import logging as logg
from cellrank._key import Key
from cellrank.ul._docs import d
from cellrank.ul._utils import _read_graph_data
from cellrank.tl.kernels._base_kernel import (
    _RTOL,
    KernelExpression,
    UnidirectionalKernel,
)

import numpy as np
from scipy.sparse import spmatrix, csr_matrix

__all__ = ("PrecomputedKernel",)


@d.dedent
class PrecomputedKernel(UnidirectionalKernel):
    """
    Kernel which contains a precomputed transition matrix.

    Parameters
    ----------
    transition_matrix
        Row-normalized transition matrix or a key in :attr:`anndata.AnnData.obsp`.
        or a :class:`cellrank.tl.kernels.KernelExpression` with a precomputed transition matrix.
    %(adata)s
        If `None`, a temporary placeholder object is created.
    %(backward)s
    kwargs
        Keyword arguments for :class:`cellrank.kernels.Kernel`.
    """

    # TODO(michalk8): adata first?
    def __init__(
        self,
        transition_matrix: Optional[
            Union[str, np.ndarray, spmatrix, KernelExpression]
        ] = None,
        adata: Optional[AnnData] = None,
        **kwargs: Any,
    ):
        origin = "'array'"
        params: Dict[str, Any] = {}

        if transition_matrix is None:
            transition_matrix = Key.uns.kernel()
            logg.info(f"Accessing `adata.obsp[{transition_matrix!r}]`")

        if isinstance(transition_matrix, str):
            # TODO(michalk8): update message
            if adata is None:
                raise ValueError(
                    "When `transition_matrix` specifies a key to `adata.obsp`, `adata` cannot be None."
                )
            origin = f"adata.obsp[{transition_matrix!r}]"
            transition_matrix = _read_graph_data(adata, transition_matrix)
        elif isinstance(transition_matrix, KernelExpression):
            # TODO(michalk8): use .transition_matrix
            if transition_matrix._transition_matrix is None:
                raise ValueError(
                    "Compute transition matrix first as `.compute_transition_matrix()`."
                )
            if adata is not None and adata is not transition_matrix.adata:
                logg.warning(
                    "Ignoring supplied `adata` object because it differs from the kernel's `adata` object."
                )

            # use `str` rather than `repr` because it captures the parameters
            origin = str(transition_matrix).strip("~<>")
            params = transition_matrix.params.copy()
            adata = transition_matrix.adata
            transition_matrix = transition_matrix.transition_matrix

        if not isinstance(transition_matrix, (np.ndarray, spmatrix)):
            raise TypeError(
                f"Expected transition matrix to be of type `numpy.ndarray` or `scipy.sparse.spmatrix`, "
                f"found `{type(transition_matrix).__name__}`."
            )

        if transition_matrix.shape[0] != transition_matrix.shape[1]:
            raise ValueError(
                f"Expected transition matrix to be square, found `{transition_matrix.shape}`."
            )

        if not np.allclose(np.sum(transition_matrix, axis=1), 1.0, rtol=_RTOL):
            # TODO(michalk8): update message
            raise ValueError("Not a valid transition matrix, not all rows sum to 1.")

        if adata is None:
            logg.warning("Creating empty `AnnData` object")
            adata = AnnData(
                csr_matrix((transition_matrix.shape[0], 1), dtype=np.float64)
            )

        super().__init__(adata, **kwargs)

        self._transition_matrix = csr_matrix(transition_matrix)
        self._params = params
        self._origin = origin

    def compute_transition_matrix(self, *_: Any, **__: Any) -> "PrecomputedKernel":
        """Return self."""
        return self

    def __repr__(self):
        return f"<{self.__class__.__name__}[origin={self._origin!r}]>"

    def __str__(self):
        return repr(self)
