from typing import Any, Union, Optional

from abc import ABC, abstractmethod

from cellrank import logging as logg
from cellrank.tl._utils import _connected, _symmetric
from cellrank.ul._utils import _get_neighs

import numpy as np
from scipy.sparse import spdiags, spmatrix


class ConnectivityMixin:
    def _read_from_adata(
        self,
        conn_key: str = "connectivities",
        check_connectivity: bool = False,
        **kwargs: Any,
    ) -> None:
        super()._read_from_adata(**kwargs)

        self._conn = _get_neighs(
            self.adata, mode="connectivities", key=conn_key
        ).astype(np.float64)
        if check_connectivity:
            if not _connected(self._conn):
                logg.warning("KNN graph is not connected")
            else:
                logg.debug("KNN graph is connected")

        if not _symmetric(self._conn):
            logg.warning("KNN graph is not symmetric")
        else:
            logg.debug("KNN graph is symmetric")

    def _density_normalize(
        self, matrix: Union[np.ndarray, spmatrix]
    ) -> Union[np.ndarray, spmatrix]:
        """
        Density normalization by the underlying KNN graph.

        Parameters
        ----------
        matrix
            Matrix to normalize.

        Returns
        -------
        Density normalized transition matrix.
        """
        logg.debug("Density-normalizing the transition matrix")

        q = np.asarray(self._conn.sum(axis=0))
        Q = spdiags(1.0 / q, 0, matrix.shape[0], matrix.shape[0])

        return Q @ matrix @ Q


class UnidirectionalMixin:
    @property
    def backward(self) -> None:
        """None."""
        return None


class BidirectionalMixin(ABC):
    @abstractmethod
    def __invert__(self) -> "BidirectionalMixin":
        pass

    @property
    def backward(self) -> Optional[bool]:
        """Direction of the process."""
        return self._backward
