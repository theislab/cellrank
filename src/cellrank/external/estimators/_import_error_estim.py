from typing import Any

from anndata import AnnData
from cellrank.estimators import BaseEstimator
from cellrank.external._error_mixin import ImportErrorMixin

__all__ = ["ErroredEstimator"]


class ErroredEstimator(ImportErrorMixin, BaseEstimator):
    """
    Utility estimator class which always throw :class:`ImportError` when instantiated.

    Subclasses can modify the message by overriding the class attribute ``__import_error_message__``.
    """

    __import_error_message__ = "Unable to import external estimator."

    def _fit_terminal_states(self, *args: Any, **kwargs: Any) -> None:
        raise NotImplementedError

    def fit(self, *args: Any, **kwargs: Any) -> "ErroredEstimator":
        """Not implemented."""
        raise NotImplementedError

    def predict(self, *args: Any, **kwargs: Any) -> None:
        """Not implemented."""
        raise NotImplementedError

    def _read_from_adata(self, adata: AnnData, **kwargs: Any) -> bool:
        raise NotImplementedError
