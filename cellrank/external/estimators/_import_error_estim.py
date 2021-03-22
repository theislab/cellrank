from typing import Any

from cellrank.tl.estimators import BaseEstimator
from cellrank.external._error_mixin import ImportErrorMixin


class ErroredEstimator(ImportErrorMixin, BaseEstimator):
    """
    Utility estimator class which always throw :class:`ImportError` when instantiated.

    Subclasses can modify the message by overriding `__import_error_message__`.
    """

    __import_error_message__ = "Unable to import the estimator."

    def _fit_terminal_states(self, *args: Any, **kwargs: Any) -> None:
        raise NotImplementedError

    def compute_terminal_states(self, *args: Any, **kwargs: Any) -> None:  # noqa: D102
        raise NotImplementedError
