from typing import Any, Optional


class ImportErrorMixin:
    """
    Mixin class which always raises :class:`ImportError`.

    Subclasses can modify the message by overriding `__import_error_message__`.

    Parameters
    ----------
    args
        Ignored.
    kwargs
        Ignored.

    Raises
    ------
    ImportError
        Always.
    """

    __import_error_message__ = "Unable to import the class."

    def __init__(self, *args: Any, **kwargs: Any):
        raise ImportError(self.__import_error_message__) from self.__error__

    def __init_subclass__(cls, error: Optional[Exception] = None, **kwargs: Any):
        super().__init_subclass__(**kwargs)
        cls.__error__ = error
