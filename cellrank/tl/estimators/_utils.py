from typing import Any, Union, Optional, Sequence

from copy import deepcopy
from collections import Iterable

from scanpy import logging as logg


class SafeGetter:
    def __init__(self, obj: Any, allowed: Union[type, Sequence[type]] = ()):
        self._exc: Optional[Exception] = None
        self._obj = obj
        if not isinstance(allowed, Iterable):
            allowed = (allowed,)
        self._allowed = allowed

        # TODO: better way?
        # deepcopy to be safe (e.g. `_shadow_adata`)
        self._dict = deepcopy(obj.__dict__)

    def __enter__(self) -> "SafeGetter":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> bool:
        # TODO: refactor using custom exc instead of passing flags?
        self._exc = exc_type
        # TODO: log properly
        if not self.ok:
            logg.warning(
                f"TODO: lower verbosity. The estimator will be incompletely initialized, reason: {exc_type(exc_val)}"
            )
            self._obj.__dict__ = self._dict
            return self._exc in self._allowed

        return True

    @property
    def ok(self) -> bool:
        return self._exc is None
