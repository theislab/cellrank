from typing import Any, Dict, Tuple, Union, Mapping, Optional
from typing_extensions import Literal

from abc import ABC, abstractmethod
from copy import copy as copy_
from pathlib import Path
from contextlib import contextmanager

from anndata import AnnData
from cellrank.tl import Lineage
from cellrank.tl.kernels import PrecomputedKernel
from cellrank.tl._estimators.mixins import IOMixin, KernelMixin, AnnDataMixin
from cellrank.tl.kernels._base_kernel import KernelExpression
from cellrank.tl._estimators.mixins._constants import Key

import numpy as np
import pandas as pd
from scipy.sparse import spmatrix, csr_matrix


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
                    f"Unable to find transition matrix in `adata.obsp[{obsp_key!r}]`."
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
        # TODO: allow non-skeleton?
        self._shadow_adata = AnnData(
            X=csr_matrix(self.adata.shape, dtype=self.adata.X.dtype),
            obs=self.adata.obs[[]].copy(),
            var=self.adata.var[[]].copy(),
            raw=self.adata.raw.to_adata(),  # TODO: better way?
        )

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
        shadow_only: bool = False,
    ):
        if not self._in_shadow:
            if attr is not None:
                if not hasattr(self, attr):
                    raise AttributeError(attr)
                setattr(self, attr, value)
            if shadow_only:
                return

        if key is not None:
            if value is None:
                try:
                    # TODO: necessary?
                    if isinstance(obj, pd.DataFrame):
                        obj.drop(key, axis="columns", inplace=True)
                    else:
                        del obj[key]
                except KeyError:
                    pass
            else:
                obj[key] = copy_(value) if copy else value

    def _get(
        self,
        attr: str,
        obj: Optional[Union[pd.DataFrame, Mapping[str, Any]]],
        key: str,
        where: Optional[Literal["obs", "obsm", "var", "varm", "uns"]] = None,
        dtype: Optional[Union[type, Tuple[type, ...]]] = None,
        copy: bool = True,
        allow_missing: bool = False,
    ):
        if not hasattr(self, attr):
            raise AttributeError(attr)

        try:
            data = obj[key]
            if dtype is not None and not isinstance(data, dtype):
                raise TypeError(
                    f"Expected `.{attr}` to be of type `{dtype}`, found `{type(data).__name__}`."
                )
            if copy:
                data = copy_(data)
            setattr(self, attr, data)
            if where is not None:
                getattr(self._shadow_adata, where)[key] = data
        except KeyError:
            if not allow_missing:
                raise

    @property
    @contextmanager
    def _shadow(self) -> None:
        if self._in_shadow:
            yield
        else:
            adata = self.adata
            try:
                self._adata = self._shadow_adata
                yield
            finally:
                self._adata = adata

    @property
    def _in_shadow(self) -> bool:
        return self.adata is self._shadow_adata

    @property
    @contextmanager
    def _remove_adata(self) -> None:
        # TODO: handle shadow?
        with super()._remove_adata:
            adata = self.kernel.adata
            try:
                self.kernel.adata = None
                yield
            finally:
                self.kernel.adata = adata

    @staticmethod
    def read(
        fname: Union[str, Path], adata: Optional[AnnData] = None, copy: bool = False
    ) -> "BaseEstimator":
        obj: BaseEstimator = IOMixin.read(fname, adata, copy=copy)
        obj.kernel.adata = adata

        return obj

    def to_adata(self) -> AnnData:
        adata = self._shadow_adata.copy()
        try:
            self.kernel.adata = adata
            self.kernel.write_to_adata(None)
        finally:
            self.kernel.adata = self.adata
        self._shadow_adata.uns[
            Key.uns.estimator(None, self.backward) + "_params"
        ] = self.params.copy()
        return adata

    @classmethod
    def from_adata(cls, adata: AnnData, obsp_key: str) -> "BaseEstimator":
        obj = cls(adata, obsp_key=obsp_key)
        obj._read_from_adata(adata)

        return obj

    @property
    def params(self) -> Dict[str, Any]:
        """TODO."""
        return self._params

    @abstractmethod
    def fit(self, *args: Any, **kwargs: Any) -> None:
        """TODO."""
