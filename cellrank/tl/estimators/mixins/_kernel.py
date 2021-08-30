from typing import Any, Union, Optional

from abc import ABC

from anndata import AnnData
from cellrank.tl.kernels._base_kernel import KernelExpression

import numpy as np
from scipy.sparse import spmatrix


class KernelMixin(ABC):
    """Mixin that exposes various properties of :class:`cellrank.kernels.KernelExpression`."""

    def __init__(self, kernel: KernelExpression, **kwargs: Any):
        super().__init__(**kwargs)
        self._kernel = kernel
        self._n_obs = self.kernel.adata.n_obs

    @property
    def kernel(self) -> "KernelExpression":
        """Underlying kernel."""
        return self._kernel

    @property
    def adata(self) -> AnnData:
        """Annotated data object."""
        return self.kernel.adata

    @adata.setter
    def adata(self, adata: Optional[AnnData]) -> None:
        self.kernel.adata = adata

    def __len__(self) -> int:
        return self._n_obs

    @property
    def transition_matrix(self) -> Union[np.ndarray, spmatrix]:
        """Transition matrix of :attr:`kernel`."""
        return self.kernel.transition_matrix

    @property
    def backward(self) -> bool:
        """Direction of :attr:`kernel`."""
        return self.kernel.backward
