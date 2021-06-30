from typing import Optional

from cellrank.ul.models import BaseModel
from cellrank.external._error_mixin import ImportErrorMixin

import numpy as np


class ErroredModel(ImportErrorMixin, BaseModel):
    """
    Utility model class which always throw :class:`ImportError` when instantiated.

    Subclasses can modify the message by overriding `__import_error_message__`.
    """

    __import_error_message__ = "Unable to import the model."

    def fit(  # noqa: D102
        self,
        x: Optional[np.ndarray] = None,
        y: Optional[np.ndarray] = None,
        w: Optional[np.ndarray] = None,
        **kwargs,
    ) -> "BaseModel":
        raise NotImplementedError

    def predict(  # noqa: D102
        self,
        x_test: Optional[np.ndarray] = None,
        key_added: Optional[str] = "_x_test",
        **kwargs,
    ) -> np.ndarray:
        raise NotImplementedError

    def confidence_interval(  # noqa: D102
        self, x_test: Optional[np.ndarray] = None, **kwargs
    ) -> np.ndarray:
        raise NotImplementedError

    def copy(self) -> BaseModel:  # noqa: D102
        raise NotImplementedError
