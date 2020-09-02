# -*- coding: utf-8 -*-
"""Module containing CellRank contants."""

from abc import ABC, ABCMeta
from enum import Enum, EnumMeta
from typing import Any, Union, Callable
from functools import wraps

from cellrank.ul._docs import _initial, _terminal

_DEFAULT_BACKEND = "loky"


class PrettyEnum(Enum):
    """Enum with a pretty __str__ and __repr__."""

    def __repr__(self):
        return str(self)

    def __str__(self):
        return str(self.value)

    @property
    def v(self) -> Any:
        """Return the value."""
        return self.value

    @property
    def s(self) -> str:
        """Return the value as string."""
        return str(self)

    @property
    def fmt(self) -> Callable:
        """Return the formatting function given by the string value."""
        return self.s.format


class Lin(PrettyEnum):
    """Lineage aggregation type."""

    REST = "rest"
    OTHERS = "others"


class Direction(PrettyEnum):
    """Direction of the process."""

    FORWARD = "fwd"
    BACKWARD = "bwd"


class DirectionPlot(PrettyEnum):
    """Pretty direction names for plotting."""

    FORWARD = "forward"
    BACKWARD = "backward"


class DirPrefix(PrettyEnum):
    """Direction prefix."""

    FORWARD = "to"
    BACKWARD = "from"


class FinalStatesKey(PrettyEnum):
    """State key in `adata.obs`."""

    FORWARD = f"{_terminal}_states"
    BACKWARD = f"{_initial}_states"


class FinalStatesPlot(PrettyEnum):
    """Pretty state names for pl."""

    FORWARD = FinalStatesKey.FORWARD.s.replace("_", " ")
    BACKWARD = FinalStatesKey.BACKWARD.s.replace("_", " ")


class MetaKey(PrettyEnum):
    """Key for scaled, subsampled metastable state in `adata.obsm`."""

    FORWARD = f"metastable_states_{Direction.FORWARD}"
    BACKWARD = f"metastable_states_{Direction.BACKWARD}"


# FinalStatesKey and AbsProbKey must have the same suffix `_..._states` because of model.prepare
class AbsProbKey(PrettyEnum):
    """Lineage key in `adata.obsm`."""

    FORWARD = f"{DirPrefix.FORWARD}_{FinalStatesKey.FORWARD}"
    BACKWARD = f"{DirPrefix.BACKWARD}_{FinalStatesKey.BACKWARD}"


def _transition(d: Union[str, Direction]) -> str:
    return f"T_{d}"


def _lin_names(k: Union[str, AbsProbKey]) -> str:
    return f"{k}_names"


def _colors(k: Union[str, AbsProbKey, FinalStatesKey]) -> str:
    return f"{k}_colors"


def _probs(k: Union[str, FinalStatesKey]) -> str:
    return f"{k}_probs"


def _dp(k: Union[str, AbsProbKey]) -> str:
    return f"{k}_dp"


def _pretty_raise_enum(cls, fun):
    @wraps(fun)
    def wrapper(*args, **kwargs):
        try:
            return fun(*args, **kwargs)
        except ValueError as e:
            _cls, value, *_ = args
            e.args = (cls._format(value),)
            raise e

    if not issubclass(cls, ErrorFormatterABC):
        raise TypeError(f"Class `{cls}` must be subtype of `ErrorFormatterABC`.")
    elif not len(cls.__members__):
        # empty enum, for class hierarchy
        return fun

    return wrapper


class ABCEnumMeta(EnumMeta, ABCMeta):  # noqa
    def __call__(cls, *args, **kw):  # noqa
        if getattr(cls, "__error_format__", None) is None:
            raise TypeError(
                f"Can't instantiate class `{cls.__name__}` "
                f"without `__error_format__` class attribute."
            )
        return super().__call__(*args, **kw)

    def __new__(cls, clsname, superclasses, attributedict):  # noqa
        res = super().__new__(cls, clsname, superclasses, attributedict)
        res.__new__ = _pretty_raise_enum(res, res.__new__)
        return res


class ErrorFormatterABC(ABC):  # noqa
    __error_format__ = "Invalid option `{!r}` for `{}`. Valid options are: `{}`."

    @classmethod
    def _format(cls, value):
        return cls.__error_format__.format(
            value, cls.__name__, [m.value for m in cls.__members__.values()]
        )


class ModeEnum(ErrorFormatterABC, PrettyEnum, metaclass=ABCEnumMeta):  # noqa
    pass
