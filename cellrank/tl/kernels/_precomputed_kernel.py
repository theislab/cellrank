from typing import Any, Union, Optional

from anndata import AnnData
from cellrank import logging as logg
from cellrank._key import Key
from cellrank.ul._docs import d
from cellrank.ul._utils import _read_graph_data
from cellrank.tl.kernels._base_kernel import KernelExpression, UnidirectionalKernel

import numpy as np
from scipy.sparse import spmatrix, csr_matrix

__all__ = ("PrecomputedKernel",)


@d.dedent
class PrecomputedKernel(UnidirectionalKernel):
    """
    Kernel which contains a precomputed transition matrix.

    Parameters
    ----------
    object
        Can be one of the following types:

            - :class:`anndata.AnnData` - annotated data object.
            - :class:`scipy.sparse.ndarray`, :class:`numpy.ndarray` - row-normalized transition matrix.
            - :class:`cellrank.kernels.KernelExpression` - kernel expression.
            - :class:`str` - key in :attr:`anndata.AnnData.obsp` where the transition matrix is stored.
              ``adata`` must be provided in this case.
            - :class:`bool` - directionality of the transition matrix that will be used to infer its storage location.
              If `None`, the directionality will be determined automatically. ``adata`` must be provided in this case.
    %(adata)s
        Must be provided when ``object`` is :class:`str` or :class:`bool`.
    obsp_key
        Key in :attr:`anndata.AnnData.obsp` where the transition matrix is stored.
        If `None`, it will be determined automatically. Only used when ``object`` is :class:`anndata.AnnData`.
    copy
        Whether or not to copy the stored transition matrix.
    backward
        Hint whether this is a forward, backward or a unidirectional kernel. Only used when ``objects`` is
        :class:`anndata.AnnData`.
    """

    _SENTINEL = object()

    def __init__(
        self,
        object: Union[str, bool, np.ndarray, spmatrix, AnnData, KernelExpression],
        adata: Optional[AnnData] = None,
        obsp_key: Optional[str] = None,
        **kwargs: Any,
    ):
        if isinstance(object, AnnData):
            self._from_adata(object, obsp_key=obsp_key, **kwargs)
        elif isinstance(object, KernelExpression):
            self._from_kernel(object, copy=kwargs.get("copy", False))
        elif isinstance(object, (np.ndarray, spmatrix)):
            self._from_matrix(object, adata=adata, **kwargs)
        elif isinstance(adata, AnnData):
            if isinstance(object, str):
                self._from_adata(adata, obsp_key=object, **kwargs)
            elif object is None or isinstance(object, bool):
                kwargs["backward"] = object
                self._from_adata(adata, **kwargs)
            else:
                raise ValueError("Unable to interpret the data.")
        else:
            raise TypeError(
                f"Expected `object` to be either `str`, `bool`, `numpy.ndarray`, "
                f"`scipy.sparse.spmatrix`, `AnnData` or `KernelExpression`, "
                f"found `{type(object).__name__}`"
            )

    def _from_adata(
        self,
        adata: AnnData,
        obsp_key: Optional[str] = None,
        backward: Optional[bool] = _SENTINEL,
        copy: bool = False,
    ) -> None:
        if obsp_key is None:
            obsp_key = Key.uns.kernel(backward)
        tmat = _read_graph_data(adata, obsp_key)
        if backward is PrecomputedKernel._SENTINEL:
            # not ideal, assume we wan't to estimate
            if obsp_key == Key.uns.kernel(bwd=None):
                backward = None
            elif obsp_key == Key.uns.kernel(bwd=False):
                backward = False
            elif obsp_key == Key.uns.kernel(bwd=True):
                backward = True
            else:
                backward = None
        self._from_matrix(tmat, adata=adata, backward=backward, copy=copy)
        self.params["origin"] = f"adata.obsp[{obsp_key!r}]"

    def _from_kernel(self, kernel: KernelExpression, copy: bool = False) -> None:
        if kernel.transition_matrix is None:
            raise RuntimeError(
                "Compute transition matrix first as `.compute_transition_matrix()`."
            )
        self._from_matrix(
            kernel.transition_matrix,
            backward=kernel.backward,
            adata=kernel.adata,
            copy=copy,
        )
        self._params = kernel.params.copy()
        self.params["origin"] = repr(kernel)

    def _from_matrix(
        self,
        matrix: Union[np.ndarray, spmatrix],
        adata: Optional[AnnData] = None,
        backward: Optional[bool] = None,
        copy: bool = False,
    ) -> None:
        # fmt: off
        if adata is None:
            logg.warning(f"Creating empty `AnnData` object of shape `{matrix.shape[0], 1}`")
            adata = AnnData(csr_matrix((matrix.shape[0], 1), dtype=np.float64))
        super().__init__(adata)
        self._backward: Optional[bool] = backward
        self.transition_matrix = matrix.copy() if copy else matrix
        self.params['origin'] = "array"
        # fmt: on

    def compute_transition_matrix(self, *_: Any, **__: Any) -> "PrecomputedKernel":
        """Do nothing and return self."""
        return self

    @property
    def backward(self) -> Optional[bool]:
        """Direction of the process."""
        return self._backward
