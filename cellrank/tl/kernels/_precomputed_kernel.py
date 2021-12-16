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
    @classmethod
    def from_object(
        cls,
        object: Union[str, bool, np.ndarray, spmatrix, AnnData, KernelExpression],
        copy: bool = True,
        **kwargs: Any,
    ):
        """
        Parameters
        ----------
        object
            Can be one of the following types:

                - :class:`anndata.AnnData`
                - :class:`scipy.sparse.ndarray`, :class:`numpy.ndarray`
                - :class:`cellrank.kernels.KernelExpression`
                - :class:`str`
                - :class:`kwargs`
        kwargs
            TODO.

        Returns
        -------
        Kernel with supplied transition matrix.
        """
        if isinstance(object, AnnData):
            return cls.from_adata(object, copy=copy, **kwargs)
        if isinstance(object, KernelExpression):
            return cls.from_kernel(object, copy=copy)
        if isinstance(object, (np.ndarray, spmatrix)):
            return cls.from_matrix(object, copy=copy, **kwargs)
        adata = kwargs.pop("adata", None)
        if isinstance(adata, AnnData):
            if isinstance(object, str):
                return cls.from_adata(adata, obsp_key=object, copy=copy)
            if object is None or isinstance(object, bool):
                return cls.from_adata(adata, backward=object, copy=copy)
        raise TypeError("TODO")

    @classmethod
    def from_adata(
        cls,
        adata: AnnData,
        obsp_key: Optional[str] = None,
        backward: Optional[bool] = None,
        copy: bool = True,
    ) -> "PrecomputedKernel":
        if obsp_key is None:
            obsp_key = Key.uns.kernel(backward)
        tmat = _read_graph_data(adata, obsp_key)
        kernel = cls.from_matrix(tmat, adata=adata, copy=copy)
        kernel.params["origin"] = f"adata.obsp[{obsp_key!r}]"
        return kernel

    @classmethod
    def from_kernel(
        cls, kernel: KernelExpression, copy: bool = True
    ) -> "PrecomputedKernel":
        if kernel.transition_matrix is None:
            raise RuntimeError(
                f"Compute transition matrix first as `.compute_transition_matrix()`."
            )
        origin = repr(kernel).strip("~<>")
        params = kernel.params.copy()
        kernel = cls.from_matrix(
            kernel.transition_matrix, adata=kernel.adata, copy=copy
        )
        kernel._params = params
        kernel.params["origin"] = origin
        return kernel

    @classmethod
    def from_matrix(
        cls,
        matrix: Union[np.ndarray, spmatrix],
        adata: Optional[AnnData] = None,
        copy: bool = True,
    ) -> "PrecomputedKernel":
        # fmt: off
        if adata is None:
            logg.warning(f"Creating empty `AnnData` object of shape `{matrix.shape[0], 1}`")
            adata = AnnData(csr_matrix((matrix.shape[0], 1), dtype=np.float64))
        kernel = cls(adata)
        kernel.transition_matrix = matrix.copy() if copy else matrix
        kernel.params['origin'] = "array"
        # fmt: on
        return kernel

    def compute_transition_matrix(self, *_: Any, **__: Any) -> "PrecomputedKernel":
        """Do nothing and return self."""
        return self
