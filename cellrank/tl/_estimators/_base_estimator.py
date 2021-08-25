from typing import Any, Dict, Union, Optional

from abc import ABC, abstractmethod

from anndata import AnnData
from cellrank.tl.kernels import PrecomputedKernel
from cellrank.tl._estimators.mixins import IOMixin, KernelMixin, AnnDataMixin
from cellrank.tl.kernels._base_kernel import KernelExpression

import numpy as np
from scipy.sparse import spmatrix


class BaseEstimator(IOMixin, AnnDataMixin, KernelMixin, ABC):
    def __init__(
        self,
        obj: Union[AnnData, np.ndarray, spmatrix, KernelExpression],
        key: Optional[str] = None,
        obsp_key: Optional[str] = None,
        **kwargs: Any,
    ):
        if isinstance(obj, KernelExpression):
            kernel = obj
        elif isinstance(obj, (np.ndarray, spmatrix)):
            kernel = PrecomputedKernel(obj)
        elif isinstance(obj, AnnData):
            if obsp_key is None:
                raise ValueError(
                    "Specify `obsp_key=...` when supplying an `AnnData` object."
                )
            elif obsp_key not in obj.obsp.keys():
                raise KeyError(
                    f"Transition matrix not found in `adata.obsp[{obsp_key!r}]`."
                )
            kernel = PrecomputedKernel(obj.obsp[obsp_key], adata=obj)
        else:
            raise TypeError(
                f"Expected an object of type `KernelExpression`, `numpy.ndarray`, `scipy.sparse.spmatrix` "
                f"or `anndata.AnnData`, got `{type(obj).__name__}`."
            )

        if kernel._transition_matrix is None:
            # TODO: make sure tests pass
            raise ValueError(
                "Compute transition matrix first as `.compute_transition_matrix()`."
            )

        super().__init__(adata=kernel.adata, kernel=kernel, **kwargs)

        self._params: Dict[str, Any] = {}
        # TODO: make sure tests pass
        # TODO: use this key for params/adata serialization
        kernel.write_to_adata(key=key)

    def __init_subclass__(cls, **kwargs: Any):
        super().__init_subclass__()

    def to_adata(self) -> None:
        self.adata.uns["TODO_kernel"] = self.params.copy()

    @classmethod
    def from_adata(cls, adata: AnnData) -> "BaseEstimator":
        return NotImplemented

    @property
    def params(self) -> Dict[str, Any]:
        """TODO."""
        return self._params

    @abstractmethod
    def fit(self, *args: Any, **kwargs: Any) -> None:
        """TODO."""
