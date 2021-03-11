from typing import Any


class ImportErrorMixin:
    """
    Mixin class which always raises :class:`ImportError`.

    Subclasses can modify the message by overriding `__import_error_message__`.

    Parameters
    ----------
    args
        Ignored.
    kwargs
        Ignore.

    Raises
    ------
    ImportError
        Always.
    """

    __import_error_message__ = "Unable to import the kernel."

    def __init__(self, *args: Any, **kwargs: Any):
        raise ImportError(self.__import_error_message__) from None
