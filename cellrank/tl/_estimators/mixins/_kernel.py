from typing import Union

from abc import ABC

from cellrank.tl.kernels._base_kernel import KernelExpression

import numpy as np
from scipy.sparse import csr_matrix


class KernelMixin(ABC):
    """TODO."""

    def __init__(self, kernel: KernelExpression):
        self._kernel = kernel

    @property
    def kernel(self) -> "KernelExpression":
        """TODO."""
        return self._kernel

    @property
    def transition_matrix(self) -> Union[np.ndarray, csr_matrix]:
        """TODO."""
        return self.kernel.transition_matrix

    @property
    def backward(self) -> bool:
        """TODO."""
        return self.kernel.backward
