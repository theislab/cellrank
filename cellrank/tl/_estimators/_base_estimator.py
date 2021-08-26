from typing import Any, Dict, Union, Mapping, Optional

from abc import ABC, abstractmethod

from anndata import AnnData
from cellrank.tl import Lineage
from cellrank.tl.kernels import PrecomputedKernel
from cellrank.tl._estimators.mixins import IOMixin, KernelMixin, AnnDataMixin
from cellrank.tl.kernels._base_kernel import KernelExpression
from cellrank.tl._estimators.mixins._constants import Key

import numpy as np
import pandas as pd
from scipy.sparse import spmatrix


class BaseEstimator(IOMixin, AnnDataMixin, KernelMixin, ABC):
    def __init__(
        self,
        obj: Union[AnnData, np.ndarray, spmatrix, KernelExpression],
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
        self._shadow_adata = AnnData(X=self.adata.X)

    def __init_subclass__(cls, **kwargs: Any):
        super().__init_subclass__()

    def _set(
        self,
        attr: Optional[str] = None,
        obj: Optional[Union[pd.DataFrame, Mapping[str, Any]]] = None,
        key: Optional[str] = None,
        value: Optional[
            Union[np.ndarray, pd.Series, pd.DataFrame, Lineage, AnnData]
        ] = None,
        copy: bool = True,
    ):
        if attr is not None:
            if not hasattr(self, attr):
                raise AttributeError(attr)
            setattr(self, attr, value)

        if key is not None:
            if value is None:
                try:
                    if isinstance(obj, pd.DataFrame):
                        obj.drop(key, axis="columns", inplace=True)
                    else:
                        del obj[key]
                except KeyError:
                    pass
            else:
                obj[key] = value.copy() if copy else value

    def to_adata(self) -> AnnData:
        return self._shadow_adata.copy()
        # TODO: control which attrs are serialized
        adata = AnnData(X=self.adata.X)

        try:
            self.kernel.adata = adata
            self.kernel.write_to_adata(None)
        finally:
            self.kernel.adata = self.adata
        adata.uns[
            Key.uns.estimator(None, self.backward) + "_params"
        ] = self.params.copy()

        return adata

    @classmethod
    def from_adata(cls, adata: AnnData, key: str) -> "BaseEstimator":
        return NotImplemented

    @property
    def params(self) -> Dict[str, Any]:
        """TODO."""
        return self._params

    @abstractmethod
    def fit(self, *args: Any, **kwargs: Any) -> None:
        """TODO."""
