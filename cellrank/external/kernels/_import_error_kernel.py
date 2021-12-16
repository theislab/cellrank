from typing import Any

from abc import ABC

from cellrank.external._error_mixin import ImportErrorMixin
from cellrank.tl.kernels._base_kernel import Kernel, KernelExpression


class ErroredKernel(ImportErrorMixin, Kernel, ABC):
    """
    Utility kernel class which always throw :class:`ImportError` when instantiated.

    Subclasses can modify the message by overriding `__import_error_message__`.
    """

    __import_error_message__ = "Unable to import the kernel."

    def compute_transition_matrix(  # noqa: D102
        self, *args: Any, **kwargs: Any
    ) -> KernelExpression:
        raise NotImplementedError

    def copy(self, deep: bool = False) -> KernelExpression:  # noqa: D102
        raise NotImplementedError
