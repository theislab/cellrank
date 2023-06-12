import abc
from typing import Any, Union

import numpy as np
import scipy.sparse as sp

from cellrank import logging as logg
from cellrank._utils._utils import _connected, _get_neighs, _symmetric

__all__ = ["ConnectivityMixin", "UnidirectionalMixin", "BidirectionalMixin"]


class ConnectivityMixin:
    """Mixin class that reads kNN connectivities and allows for density normalization."""

    def _read_from_adata(
        self,
        conn_key: str = "connectivities",
        check_connectivity: bool = False,
        check_symmetric: bool = True,
        **kwargs: Any,
    ) -> None:
        super()._read_from_adata(**kwargs)
        self._conn_key = conn_key
        conn = _get_neighs(self.adata, mode="connectivities", key=conn_key)
        self._conn = sp.csr_matrix(conn).astype(np.float64, copy=False)

        if check_connectivity and not _connected(self.connectivities):
            logg.warning("kNN graph is not connected")
        if check_symmetric and not _symmetric(self.connectivities):
            logg.warning("kNN graph is not symmetric")

    def _density_normalize(self, matrix: Union[np.ndarray, sp.spmatrix]) -> Union[np.ndarray, sp.spmatrix]:
        """
        Density normalization by the underlying kNN graph.

        Parameters
        ----------
        matrix
            Matrix to normalize.

        Returns
        -------
        Density normalized matrix.
        """
        logg.debug("Density normalizing the transition matrix")

        q = np.asarray(self.connectivities.sum(axis=0)).squeeze()
        Q = sp.spdiags(1.0 / q, 0, matrix.shape[0], matrix.shape[0])

        return Q @ matrix @ Q

    @property
    def connectivities(self) -> sp.csr_matrix:
        """Underlying connectivity matrix."""
        return self._conn


class UnidirectionalMixin:
    """Mixin specifying that its kernel doesn't is directionless."""

    @property
    def backward(self) -> None:
        """None."""


class BidirectionalMixin(abc.ABC):
    """Mixin specifying that its kernel has forward or backward directions."""

    def __init__(self, *args: Any, backward: bool = False, **kwargs: Any):
        super().__init__(*args, **kwargs)
        if not isinstance(backward, (bool, np.bool_)):
            raise TypeError(f"Expected `backward` to be `bool`, found `{type(backward).__name__}`.")
        self._backward = bool(backward)
        self._init_kwargs["backward"] = backward

    @abc.abstractmethod
    def __invert__(self) -> "BidirectionalMixin":
        pass

    @property
    def backward(self) -> bool:
        """Direction of the process."""
        return self._backward
