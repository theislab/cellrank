import copy
import enum
from typing import (
    Any,
    Callable,
    Dict,
    Iterable,
    Literal,
    Mapping,
    NamedTuple,
    Optional,
    Protocol,
    Sequence,
    Tuple,
    TypeVar,
    Union,
)

from wrapt import decorator

import numpy as np
import pandas as pd

from anndata import AnnData

from cellrank import logging as logg
from cellrank._utils._enum import ModeEnum
from cellrank._utils._lineage import Lineage

__all__ = [
    "BaseProtocol",
    "logger",
    "shadow",
    "SafeGetter",
    "StatesHolder",
]


class StatesHolder(NamedTuple):  # noqa: D101
    assignment: Optional[pd.Series] = None
    probs: Optional[pd.Series] = None
    colors: Optional[np.ndarray] = None
    memberships: Optional[Lineage] = None

    def set(self, **kwargs: Any) -> "StatesHolder":  # noqa: D102
        return self._replace(**kwargs)


class PlotMode(ModeEnum):
    EMBEDDING = enum.auto()
    TIME = enum.auto()


class BaseProtocol(Protocol):  # noqa: D101
    @property
    def adata(self) -> AnnData:  # noqa: D102
        ...

    @property
    def backward(self) -> bool:  # noqa: D102
        ...

    @property
    def params(self) -> Dict[str, Any]:  # noqa: D102
        ...

    def _set(
        self,
        attr: Optional[str] = None,
        obj: Optional[Union[pd.DataFrame, Mapping[str, Any]]] = None,
        key: Optional[str] = None,
        value: Optional[Union[np.ndarray, pd.Series, pd.DataFrame, Lineage, AnnData, Dict[str, Any]]] = None,
        copy: bool = True,
        shadow_only: bool = False,
    ) -> None:
        ...

    def _get(
        self,
        *,
        obj: Union[pd.DataFrame, Mapping[str, Any]],
        key: str,
        shadow_attr: Optional[Literal["obs", "obsm", "var", "varm", "uns"]] = None,
        dtype: Optional[Union[type, Tuple[type, ...]]] = None,
        copy: bool = True,
        allow_missing: bool = False,
    ) -> Any:
        ...

    def _create_params(
        self,
        locs: Optional[Mapping[str, Any]] = None,
        func: Optional[Callable] = None,
        remove: Sequence[str] = (),
    ) -> Dict[str, Any]:
        ...

    def _read_params(self, key: str) -> Dict[str, Any]:
        ...


@decorator()
def logger(wrapped: Callable[..., str], instance: Any, args: Any, kwargs: Dict[str, Any]) -> str:
    """Handle logging for :class:`~anndata.AnnData` writing functions of :class:`~cellrank.estimators.BaseEstimator`."""
    log, time = kwargs.pop("log", True), kwargs.pop("time", None)
    msg = wrapped(*args, **kwargs)

    if log:
        logg.info(msg, time=time)

    return msg


@decorator()
def shadow(wrapped: Callable[..., str], instance: Any, args: Any, kwargs: Mapping[str, Any]) -> str:
    """Duplicate function call with shadow :class:`~anndata.AnnData` object.

    Used to create :class:`~anndata.AnnData` serialization object in the :class:`~cellrank.estimators.BaseEstimator`.
    """
    res = wrapped(*args, **kwargs)
    with instance._shadow:
        try:
            # don't care what the "shadowed" function returns
            _ = wrapped(*args, **kwargs)
        except Exception as e:  # noqa: BLE001
            logg.error(f"Unable to duplicate function call using shadow `anndata.AnnData` object. Reason: `{e}`")

    return res


class SafeGetter:
    """Context manager that does a simple rollback on :attr:`object.__dict__`.

    Used in :class:`~cellrank.estimators.BaseEstimator` to maintain consistency
    during deserialization from :class:`~anndata.AnnData`.

    Parameters
    ----------
    obj
        Object to watch. Shallow copy of :attr:`obj.__dict__` is made during initialization,
        except `_shadow_adata` (deep copy), if present.
    allowed
        Whether to allow these exception to occur with the context.

    Notes
    -----
    :attr:`ok` is :obj:`False` if and only if any exception, regardless of whether it is in ``allowed``,
    has been encountered.
    """

    def __init__(self, obj: Any, allowed: Union[type, Sequence[type]] = ()):
        self._exc_type: Optional[TypeVar[Exception]] = None
        self._obj = obj
        if not isinstance(allowed, Iterable):
            allowed = (allowed,)
        self._allowed = allowed

        self._dict = copy.copy(obj.__dict__)
        if "_shadow_adata" in self._dict:
            # deepcopy `_shadow_adata` (small cost) to easily revert changes
            # do not run deepcopy on the full `__dict__` (possibly large cost)
            # silent assumption is that within the context `adata` is not modified
            self._dict["_shadow_adata"] = copy.deepcopy(self._dict["_shadow_adata"])

    def __enter__(self) -> "SafeGetter":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> bool:
        self._exc_type = exc_type
        # TODO(michalk8): log properly
        if not self.ok:
            logg.debug(f"The estimator will not be completely initialized, reason: {exc_type(exc_val)}")
            self._obj.__dict__ = self._dict
            return self._exc_type in self._allowed

        return True

    @property
    def ok(self) -> bool:
        """Return :obj:`True` if no exception has been raised."""
        return self._exc_type is None
