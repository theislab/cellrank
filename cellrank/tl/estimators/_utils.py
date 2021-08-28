from typing import Any, Union, Optional, Sequence

from copy import copy, deepcopy
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
        self._dict = copy(obj.__dict__)
        if "_shadow_adata" in self._dict:
            # deepcopy `_shadow_adata` to easily revert changes
            # do not run deepcopy on the full `__dict__`: if an
            # incomplete read occurs, Estimator.from_adata(adata).adata
            # will not be supplied adata
            # silently assumes readers do not modify the adata
            self._dict["_shadow_adata"] = deepcopy(self._dict["_shadow_adata"])

    def __enter__(self) -> "SafeGetter":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> bool:
        # TODO: refactor using custom exc instead of passing flags?
        self._exc = exc_type
        # TODO: log properly
        if not self.ok:
            logg.debug(
                f"The estimator will be incompletely initialized, reason: {exc_type(exc_val)}"
            )
            self._obj.__dict__ = self._dict
            return self._exc in self._allowed

        return True

    @property
    def ok(self) -> bool:
        return self._exc is None
