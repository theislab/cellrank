from typing import Any, Optional

from cellrank.kernels._base_kernel import Kernel
from cellrank.external._error_mixin import ImportErrorMixin

__all__ = ["ErroredKernel"]


# can't subclass UnidirectionalKernel since StationaryOTKernel have to use it
class ErroredKernel(ImportErrorMixin, Kernel):
    """
    Utility kernel class which always throw :class:`ImportError` when instantiated.

    Subclasses can modify the message by overriding the class attribute ``__import_error_message__``.
    """

    __import_error_message__ = "Unable to import external kernel."

    def compute_transition_matrix(self, *args: Any, **kwargs: Any) -> "ErroredKernel":
        """Not implemented."""
        raise NotImplementedError

    def copy(self, deep: bool = False) -> "ErroredKernel":
        """Not implemented."""
        raise NotImplementedError

    @property
    def backward(self) -> Optional[bool]:
        """None."""
        return None
