from abc import ABC, abstractmethod
from typing import Any, Callable

import numpy as np
from scipy.sparse import csr_matrix

from cellrank.ul._docs import d


class PseudotimeScheme(ABC):
    """TODO."""

    @d.get_sections(base="pt_scheme", sections=["Parameters", "Returns"])
    @abstractmethod
    def __call__(
        self,
        cell_pseudotime: float,
        neigh_pseudotime: np.ndarray,
        neigh_dist: np.ndarray,
        **kwargs: Any,
    ) -> np.ndarray:
        """
        TODO.

        Parameters
        ----------
        cell_pseudotime
            Pseudotime of the current cell.
        neigh_pseudotime
            Array of shape ``(n_neighbors,)`` containing pseudotimes of neighbors
        neigh_dist
            Array of shape ``(n_neighbors,)`` containing distances from the current cells to its neighbors.

        Returns
        -------
        Array of shape ``(n_neighbors,)`` containing the biased distances.
        """

    def bias_knn(
        self, conn: csr_matrix, pseudotime: np.ndarray, **kwargs: Any
    ) -> csr_matrix:
        """
        Palantir Kernel utility function.

        TODO: update me.

        This function takes in symmetric connectivities and a pseudotime and removes edges that point "against"
        pseudotime, in this way creating a directed graph. For each node, it always keeps the closest neighbors,
        making sure the graph remains connected.

        Parameters
        ----------
        conn
            The nearest neighbor connectivities.
        pseudotime
            Pseudotemporal ordering of cells.
        kwargs
            TODO.

        Returns
        -------
        TODO.
        """
        conn_biased = conn.copy()
        for i in range(conn.shape[0]):
            row, start, end = conn[i], conn.indptr[i], conn.indptr[i + 1]

            biased_row = self(
                pseudotime[i], pseudotime[row.indices], row.data, **kwargs
            )
            if np.shape(biased_row) != row.data.shape:
                raise ValueError(
                    f"Expected row of shape `{row.data.shape}`, found `{np.shape(biased_row)}`."
                )
            conn_biased.data[start:end] = biased_row

        conn_biased.eliminate_zeros()
        return conn_biased


class HardThresholdScheme(PseudotimeScheme):
    """TODO."""

    @d.dedent
    def __call__(
        self,
        cell_pseudotime: float,
        neigh_pseudotime: np.ndarray,
        neigh_dist: np.ndarray,
        n_neighs: int,
        k: int = 3,
    ) -> np.ndarray:
        """
        TODO.

        Parameters
        ----------
        %(pt_scheme.parameters)s
        n_neighs
            Number of neighbors to keep.
        k
            Number, alongside with ``n_neighbors`` which determined the threshold for candidate indices.

        Returns
        -------
        %(pt_scheme.returns)s
        """
        k_thresh = max(0, min(30, int(np.floor(n_neighs / k)) - 1))

        # below code does not work with argpartition
        ixs = np.flip(np.argsort(neigh_dist))
        close_ixs, far_ixs = ixs[:k_thresh], ixs[k_thresh:]

        mask_keep = cell_pseudotime < neigh_pseudotime[far_ixs]
        far_ixs_keep = far_ixs[mask_keep]

        biased_dist = np.zeros_like(neigh_dist)
        biased_dist[close_ixs] = neigh_dist[close_ixs]
        biased_dist[far_ixs_keep] = neigh_dist[far_ixs_keep]

        return biased_dist


class SoftThresholdScheme(PseudotimeScheme):
    """TODO."""

    @d.dedent
    def __call__(
        self,
        cell_pseudotime: float,
        neigh_pseudotime: np.ndarray,
        neigh_dist: np.ndarray,
        b: float = 20.0,
        nu: float = 0.5,
        perc: int = 95,
    ) -> np.ndarray:
        """
        TODO.

        Parameters
        ----------
        %(pt_scheme.parameters)s

        Returns
        -------
        %(pt_scheme.returns)s
        """
        past_ixs = np.where(neigh_pseudotime < cell_pseudotime)[0]
        if not len(past_ixs):
            return neigh_dist * 0.5

        # TODO: should all of them be clipped, or just the ones in the past?
        neigh_dist[past_ixs] = np.clip(
            neigh_dist[past_ixs],
            np.percentile(neigh_dist[past_ixs], 100 - perc),
            np.percentile(neigh_dist[past_ixs], perc),
        )

        weights = np.full_like(neigh_dist, fill_value=0.5)

        dt = cell_pseudotime - neigh_pseudotime[past_ixs]
        weights[past_ixs] = 1.0 / ((1.0 + np.exp(b * dt)) ** (1.0 / nu))

        return neigh_dist * weights


class CustomThresholdScheme(PseudotimeScheme):
    """TODO."""

    def __init__(
        self,
        callback: Callable[
            [float, np.ndarray, np.ndarray, np.ndarray, Any], np.ndarray
        ],
    ):
        super().__init__()
        self._callback = callback

    def __call__(
        self,
        cell_pseudotime: float,
        neigh_pseudotime: np.ndarray,
        neigh_dist: np.ndarray,
        **kwargs: Any,
    ) -> np.ndarray:
        """
        TODO.

        Parameters
        ----------
        %(pt_scheme.parameters)s
        kwargs
            Additional keyword arguments.

        Returns
        -------
        %(pt_scheme.returns)s
        """
        return self._callback(cell_pseudotime, neigh_pseudotime, neigh_dist, **kwargs)
