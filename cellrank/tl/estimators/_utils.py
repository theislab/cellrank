from typing import Any, Union, TypeVar, Optional, Sequence

from copy import copy, deepcopy
from collections import Iterable

from scanpy import logging as logg


class SafeGetter:
    """
    Context manager that does a simple rollback on :attr:`object.__dict__`.

    Used in :class:`cellrank.estimators.BaseEstimator` to maintain consistency
    during deserialization from :class:`anndata.AnnData`.

    Parameters
    ----------
    obj
        Object to watch. Shallow copy of ``obj.__dict__`` is made during initialization, with the exception
        of the `_shadow_adata` (deep copy), if present.
    allowed
        Allow these exception to occur with the context.

    Notes
    -----
    :attr:`ok` is `False` iff any exception, regardless of whether it is in ``allowed``, has been encountered.
    """

    def __init__(self, obj: Any, allowed: Union[type, Sequence[type]] = ()):
        self._exc_type: Optional[TypeVar[Exception]] = None
        self._obj = obj
        if not isinstance(allowed, Iterable):
            allowed = (allowed,)
        self._allowed = allowed

        self._dict = copy(obj.__dict__)
        if "_shadow_adata" in self._dict:
            # deepcopy `_shadow_adata` (small cost) to easily revert changes
            # do not run deepcopy on the full `__dict__` (possibly large cost)
            # silent assumption is that within the context `adata` is not modified
            self._dict["_shadow_adata"] = deepcopy(self._dict["_shadow_adata"])

    def __enter__(self) -> "SafeGetter":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> bool:
        self._exc_type = exc_type
        # TODO(michalk8): log properly
        if not self.ok:
            logg.debug(
                f"The estimator will not be completely initialized, reason: {exc_type(exc_val)}"
            )
            self._obj.__dict__ = self._dict
            return self._exc_type in self._allowed

        return True

    @property
    def ok(self) -> bool:
        """Return `True` if no exception has been raised."""
        return self._exc_type is None
