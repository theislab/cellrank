from typing import Any, Optional, Tuple, TypeVar, Union

import numpy as np
import scipy.sparse as sp

from anndata import AnnData

__all__ = ["KernelMixin"]

KernelExpression = TypeVar("KernelExpression", bound="KernelMixin")


class KernelMixin:
    """Mixin that exposes various properties of :class:`cellrank.kernels.KernelExpression`."""

    def __init__(self, kernel: "KernelExpression", **kwargs: Any):
        super().__init__(**kwargs)
        self._kernel = kernel
        self._n_obs = self.kernel.adata.n_obs

    @property
    def kernel(self) -> "KernelExpression":
        """Underlying kernel expression."""
        return self._kernel

    @property
    def adata(self) -> AnnData:
        """Annotated data object."""
        return self.kernel.adata

    @property
    def shape(self) -> Tuple[int, int]:
        """Shape of the kernel."""
        return self.kernel.shape

    @adata.setter
    def adata(self, adata: Optional[AnnData]) -> None:
        self.kernel.adata = adata

    def __len__(self) -> int:
        return self._n_obs

    @property
    def transition_matrix(self) -> Union[np.ndarray, sp.spmatrix]:
        """Transition matrix of the :attr:`kernel`."""
        return self.kernel.transition_matrix

    @property
    def backward(self) -> bool:
        """Direction of the :attr:`kernel`."""
        return self.kernel.backward
