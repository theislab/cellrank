from typing import Any, Dict, Tuple, Union, Mapping, Callable, Optional, Sequence
from typing_extensions import Literal

from abc import ABC, abstractmethod
from copy import copy
from copy import copy as copy_
from copy import deepcopy
from inspect import Parameter, signature, getmembers, currentframe
from contextlib import contextmanager

from anndata import AnnData
from cellrank._key import Key
from cellrank.tl.kernels import PrecomputedKernel
from cellrank.tl._lineage import Lineage
from cellrank.tl.estimators.mixins import IOMixin, KernelMixin, AnnDataMixin
from cellrank.tl.kernels._base_kernel import Kernel, KernelExpression

import numpy as np
import pandas as pd
from scipy.sparse import spmatrix, csr_matrix


class BaseEstimator(IOMixin, KernelMixin, AnnDataMixin, ABC):
    def __init__(
        self,
        obj: Union[AnnData, np.ndarray, spmatrix, KernelExpression],
        obsp_key: Optional[str] = None,
        **kwargs: Any,
    ):
        if isinstance(obj, Kernel):
            if obj._transition_matrix is None:
                raise ValueError(
                    "Compute transition matrix first as `.compute_transition_matrix()`."
                )
            kernel = obj
        elif isinstance(obj, KernelExpression):
            # this will fail if not all kernels have transition matrix computed
            kernel = obj.compute_transition_matrix()
        elif isinstance(obj, (np.ndarray, spmatrix)):
            kernel = PrecomputedKernel(obj)
        elif isinstance(obj, AnnData):
            if obsp_key is None:
                raise ValueError(
                    "Specify `obsp_key=...` when supplying an `AnnData` object."
                )
            elif obsp_key not in obj.obsp:
                raise KeyError(
                    f"Unable to find transition matrix in `adata.obsp[{obsp_key!r}]`."
                )
            kernel = PrecomputedKernel(obsp_key, adata=obj)
        else:
            raise TypeError(
                f"Expected an object of type `KernelExpression`, `numpy.ndarray`, `scipy.sparse.spmatrix` "
                f"or `anndata.AnnData`, got `{type(obj).__name__}`."
            )

        super().__init__(kernel=kernel, **kwargs)

        self._params: Dict[str, Any] = {}
        self._shadow_adata = AnnData(
            X=csr_matrix(self.adata.shape, dtype=self.adata.X.dtype),
            obs=self.adata.obs[[]].copy(),
            var=self.adata.var[[]].copy(),
            raw=None if self.adata.raw is None else self.adata.raw.to_adata(),
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
        """TODO: docs."""
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
        """TODO: docs."""
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
                self.adata = self._shadow_adata
                yield
            finally:
                self.adata = adata

    @property
    def _in_shadow(self) -> bool:
        return self.adata is self._shadow_adata

    @property
    @contextmanager
    def _remove_adata(self) -> None:
        # shadow is kept (should be very small)
        with super()._remove_adata:
            adata = self.kernel.adata
            try:
                self.kernel.adata = None
                yield
            finally:
                self.kernel.adata = adata

    def _create_params(
        self,
        locs: Optional[Mapping[str, Any]] = None,
        func: Optional[Callable] = None,
        remove: Sequence[str] = (),
    ) -> Dict[str, Any]:
        frame = currentframe()
        try:
            if locs is None:
                locs = frame.f_back.f_locals
            if func is None or True:
                name = frame.f_back.f_code.co_name
                func = dict(getmembers(self)).get(name, None)
            if not callable(func):
                raise TypeError("TODO")

            params = {}
            for name, param in signature(func).parameters.items():
                if param.kind in (Parameter.VAR_POSITIONAL, Parameter.VAR_KEYWORD):
                    continue
                if name in remove:
                    continue
                if name in locs:
                    params[name] = locs[name]

            return params
        except AttributeError:
            # frame can be None
            return {}
        except TypeError:
            # TODO: parse error and log
            return {}
        finally:
            del frame

    def _read_params(self, key: str) -> Dict[str, Any]:
        ekey = Key.uns.estimator(self.backward) + "_params"
        return dict(self.adata.uns.get(ekey, {}).get(key, {}))

    def to_adata(self) -> AnnData:
        """TODO: docrep."""
        adata = self._shadow_adata.copy()
        try:
            # kernel and estimator share the adata
            self.adata = adata
            self.kernel.write_to_adata()
        finally:
            self.adata = self.adata
        key = Key.uns.estimator(self.backward) + "_params"
        adata.uns[key] = deepcopy(self.params)
        return adata

    @classmethod
    def from_adata(cls, adata: AnnData, obsp_key: str) -> "BaseEstimator":
        """TODO: docrep."""
        return super().from_adata(adata, obsp_key=obsp_key)

    def copy(self, *, deep: bool = False) -> "BaseEstimator":
        k = deepcopy(self.kernel) if deep else copy(self.kernel)
        res = type(self)(k)
        for k, v in self.__dict__.items():
            if isinstance(v, Mapping):
                res.__dict__[k] = deepcopy(v)
            elif k != "_kernel":
                res.__dict__[k] = deepcopy(v) if deep else copy(v)

        return res

    def __copy__(self) -> "BaseEstimator":
        return self.copy(deep=False)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}[n={len(self)}, kernel={repr(self.kernel)}]"

    def __str__(self) -> str:
        return f"{self.__class__.__name__}[n={len(self)}, kernel={str(self.kernel)}]"

    @property
    def params(self) -> Dict[str, Any]:
        """TODO."""
        return self._params

    @abstractmethod
    def fit(self, *args: Any, **kwargs: Any) -> None:
        """TODO."""
