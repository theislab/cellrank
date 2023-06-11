import abc
import enum
import functools
from typing import Any, Callable, Dict, Literal, Tuple, Type

__all__ = ["ModeEnum", "DEFAULT_BACKEND"]

DEFAULT_BACKEND = "loky"
Backend_t = Literal["loky", "multiprocessing", "threading"]


class PrettyEnum(enum.Enum):
    """Enum with a pretty :meth:`__str__` and :meth:`__repr__`."""

    @property
    def v(self) -> Any:
        """Alias for :attr`value`."""
        return self.value

    def __repr__(self) -> str:
        return f"{self.value!r}"

    def __str__(self) -> str:
        return f"{self.value!s}"


def _pretty_raise_enum(cls: Type["ErrorFormatterABC"], func: Callable) -> Callable:
    @functools.wraps(func)
    def wrapper(*args: Any, **kwargs: Any) -> "ErrorFormatterABC":
        try:
            return func(*args, **kwargs)
        except ValueError as e:
            _cls, value, *_ = args
            e.args = (cls._format(value),)
            raise e

    if not issubclass(cls, ErrorFormatterABC):
        raise TypeError(f"Class `{cls}` must be subtype of `ErrorFormatterABC`.")
    if not len(cls.__members__):
        # empty enum, for class hierarchy
        return func

    return wrapper


class ABCEnumMeta(enum.EnumMeta, abc.ABCMeta):  # noqa: D101
    def __call__(cls, *args, **kwargs):  # noqa
        if getattr(cls, "__error_format__", None) is None:
            raise TypeError(f"Can't instantiate class `{cls.__name__}` " f"without `__error_format__` class attribute.")
        return super().__call__(*args, **kwargs)

    def __new__(cls, clsname: str, superclasses: Tuple[type], attributedict: Dict[str, Any]):  # noqa: D102
        res = super().__new__(cls, clsname, superclasses, attributedict)
        res.__new__ = _pretty_raise_enum(res, res.__new__)
        return res


class ErrorFormatterABC(abc.ABC):  # noqa: D101
    __error_format__ = "Invalid option `{!r}` for `{}`. Valid options are: `{}`."

    @classmethod
    def _format(cls, value) -> str:
        return cls.__error_format__.format(value, cls.__name__, [m.value for m in cls.__members__.values()])


class ModeEnum(str, ErrorFormatterABC, PrettyEnum, metaclass=ABCEnumMeta):  # noqa: D101
    def _generate_next_value_(self, start, count, last_values):
        return str(self).lower()
