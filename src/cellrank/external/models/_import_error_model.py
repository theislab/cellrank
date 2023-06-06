from typing import Any, Optional

from cellrank.models import BaseModel
from cellrank.external._error_mixin import ImportErrorMixin

import numpy as np

__all__ = ["ErroredModel"]


class ErroredModel(ImportErrorMixin, BaseModel):
    """
    Utility model class which always throw :class:`ImportError` when instantiated.

    Subclasses can modify the message by overriding the class attribute ``__import_error_message__``.
    """

    __import_error_message__ = "Unable to import external model."

    def fit(
        self,
        x: Optional[np.ndarray] = None,
        y: Optional[np.ndarray] = None,
        w: Optional[np.ndarray] = None,
        **kwargs: Any,
    ) -> "ErroredModel":
        """Not implemented."""
        raise NotImplementedError

    def predict(
        self,
        x_test: Optional[np.ndarray] = None,
        key_added: Optional[str] = "_x_test",
        **kwargs: Any,
    ) -> np.ndarray:
        """Not implemented."""
        raise NotImplementedError

    def confidence_interval(
        self, x_test: Optional[np.ndarray] = None, **kwargs: Any
    ) -> np.ndarray:
        """Not implemented."""
        raise NotImplementedError

    def copy(self) -> "ErroredModel":
        """Not implemented."""
        raise NotImplementedError
