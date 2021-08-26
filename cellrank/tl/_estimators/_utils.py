from typing import Any, Union, Optional, Sequence

from collections import Iterable


class SafeGetter:
    def __init__(self, obj: Any, allowed: Union[type, Sequence[type]] = ()):
        self._exc: Optional[Exception] = None
        self._obj = obj
        if not isinstance(allowed, Iterable):
            allowed = (allowed,)
        self._allowed = allowed

        # TODO: is this the "best" way?
        self._dict = obj.__dict__.copy()

    def __enter__(self) -> "SafeGetter":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> bool:
        self._exc = exc_type
        if not self.ok:
            self._obj.__dict__ = self._dict

        return exc_type in self._allowed

    @property
    def ok(self) -> bool:
        return self._exc is None
