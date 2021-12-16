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


class PrecomputedKernel(UnidirectionalKernel):
    """
    Kernel which contains a precomputed transition matrix.

    Parameters
    ----------
    args
        Positional arguments for :class:`cellrank.kernels.UnidirectionalKernel`.
    kwargs
        Keyword arguments for :class:`cellrank.kernels.UnidirectionalKernel`.
    """

    def __init__(self, *args: Any, **kwargs: Any):
        super().__init__(*args, **kwargs)
        self._origin: Optional[str] = None

    @classmethod
    def from_object(
        cls,
        object: Union[str, bool, np.ndarray, spmatrix, AnnData, KernelExpression],
        **kwargs: Any,
    ):
        if isinstance(object, AnnData):
            return cls.from_adata(object, **kwargs)
        if isinstance(object, KernelExpression):
            return cls.from_kernel(object)
        if isinstance(object, (np.ndarray, spmatrix)):
            return cls.from_matrix(object, **kwargs)
        adata = kwargs.pop("adata", None)
        if isinstance(adata, AnnData):
            if isinstance(object, str):
                return cls.from_adata(adata, obsp_key=object)
            if object is None or isinstance(object, bool):
                return cls.from_adata(adata, backward=object)
        raise TypeError("TODO")

    @classmethod
    @d.dedent
    def from_adata(
        cls,
        adata: AnnData,
        obsp_key: Optional[str] = None,
        backward: Optional[bool] = None,
    ) -> "PrecomputedKernel":
        """
        TODO.

        Parameters
        ----------
        %(adats)s
        """
        if obsp_key is None:
            obsp_key = Key.uns.kernel(backward)
        tmat = _read_graph_data(adata, obsp_key)
        kernel = cls.from_matrix(tmat, adata=adata)
        kernel._origin = f"adata.obsp[{obsp_key!r}]"
        return kernel

    @classmethod
    def from_kernel(cls, kernel: KernelExpression) -> "PrecomputedKernel":
        if isinstance(kernel.transition_matrix, None):
            raise RuntimeError()
        kernel = cls.from_matrix(kernel.transition_matrix, adata=kernel.adata)
        kernel._origin = str(kernel).strip("~<>")
        kernel._params = kernel.params.copy()
        return kernel

    @classmethod
    def from_matrix(
        cls, matrix: Union[np.ndarray, spmatrix], adata: Optional[AnnData] = None
    ) -> "PrecomputedKernel":
        if adata is None:
            logg.warning("Creating empty `AnnData` object")
            adata = AnnData(csr_matrix((matrix.shape[0], 1), dtype=np.float64))
        kernel = super().__init__(adata)
        kernel.transition_matrix = matrix
        kernel._origin = "array"
        return kernel

    def compute_transition_matrix(self, *_: Any, **__: Any) -> "PrecomputedKernel":
        """Return self."""
        return self

    def __repr__(self):
        return f"<{self.__class__.__name__}[origin={self._origin!r}]>"

    def __str__(self):
        return repr(self)
