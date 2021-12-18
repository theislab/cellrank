from typing import Any, Union, Optional

from anndata import AnnData
from cellrank import logging as logg
from cellrank._key import Key
from cellrank.ul._utils import _read_graph_data
from cellrank.tl.kernels._base_kernel import KernelExpression, UnidirectionalKernel

import numpy as np
from scipy.sparse import spmatrix, csr_matrix

__all__ = ("PrecomputedKernel",)


class PrecomputedKernel(UnidirectionalKernel):
    """
    TODO.

    Parameters
    ----------
    object
        Can be one of the following types:

            - :class:`anndata.AnnData`
            - :class:`scipy.sparse.ndarray`, :class:`numpy.ndarray`
            - :class:`cellrank.kernels.KernelExpression`
            - :class:`str`
            - :class:`bool`
    kwargs
        TODO.

    """

    def __init__(
        self,
        object: Union[str, bool, np.ndarray, spmatrix, AnnData, KernelExpression],
        copy: bool = True,
        **kwargs: Any,
    ):
        if isinstance(object, AnnData):
            self._from_adata(object, copy=copy, **kwargs)
        elif isinstance(object, KernelExpression):
            self._from_kernel(object, copy=copy)
        elif isinstance(object, (np.ndarray, spmatrix)):
            self._from_matrix(object, copy=copy, **kwargs)
        else:
            adata = kwargs.pop("adata", None)
            if isinstance(adata, AnnData):
                if isinstance(object, str):
                    self._from_adata(adata, obsp_key=object, copy=copy)
                if object is None or isinstance(object, bool):
                    self._from_adata(adata, backward=object, copy=copy)
                return
            raise TypeError(
                f"Expected object to be either `str`, `bool`, `numpy.ndarray`, "
                f"`scipy.sparse.spmatrix`, `AnnData` or `KernelExpression`, "
                f"found `{type(object).__name__}`"
            )

    def _from_adata(
        self,
        adata: AnnData,
        obsp_key: Optional[str] = None,
        backward: Optional[bool] = None,
        copy: bool = True,
    ) -> None:
        if obsp_key is None:
            obsp_key = Key.uns.kernel(backward)
        tmat = _read_graph_data(adata, obsp_key)
        self._from_matrix(tmat, adata=adata, copy=copy)
        if obsp_key == Key.uns.kernel(bwd=False):
            self._backward = False
        elif obsp_key == Key.uns.kernel(bwd=True):
            self._backward = True
        self.params["origin"] = f"adata.obsp[{obsp_key!r}]"

    def _from_kernel(self, kernel: KernelExpression, copy: bool = True) -> None:
        if kernel.transition_matrix is None:
            raise RuntimeError(
                "Compute transition matrix first as `.compute_transition_matrix()`."
            )
        self._from_matrix(kernel.transition_matrix, adata=kernel.adata, copy=copy)
        self._params = kernel.params.copy()
        self._backward = kernel.backward
        self.params["origin"] = repr(kernel)

    def _from_matrix(
        self,
        matrix: Union[np.ndarray, spmatrix],
        adata: Optional[AnnData] = None,
        copy: bool = True,
    ) -> None:
        # fmt: off
        if adata is None:
            logg.warning(f"Creating empty `AnnData` object of shape `{matrix.shape[0], 1}`")
            adata = AnnData(csr_matrix((matrix.shape[0], 1), dtype=np.float64))
        super().__init__(adata)
        self._backward: Optional[bool] = None
        self.transition_matrix = matrix.copy() if copy else matrix
        self.params['origin'] = "array"
        # fmt: on

    def compute_transition_matrix(self, *_: Any, **__: Any) -> "PrecomputedKernel":
        """Do nothing and return self."""
        return self

    @property
    def backward(self) -> Optional[bool]:
        return self._backward
