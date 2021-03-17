from abc import ABC, abstractmethod
from typing import Any, Callable, Optional

import numpy as np
from scipy.sparse import csr_matrix

from cellrank.ul._docs import d


class ThresholdSchemeABC(ABC):
    """Base class for all connectivity biasing schemes."""

    @d.get_summary(base="pt_scheme")
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
        Calculate biased connections for a given cell.

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
        Bias cell-cell connectivities of a KNN graph.

        Parameters
        ----------
        conn
            Sparse matrix of shape ``(n_cells, n_cells)`` containing the nearest neighbor connectivities.
        pseudotime
            Pseudotemporal order of cells.

        Returns
        -------
        The biased connectivities.
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

    def __repr__(self):
        return f"<{self.__class__.__name__}>"

    def __str__(self):
        return repr(self)


class HardThresholdScheme(ThresholdSchemeABC):
    """
    Thresholding scheme inspired by Palantir [Setty19]_.

    Note that this won't exactly reproduce the original Palantir results, for three reasons:

        - Palantir computes the KNN graph in a scaled space of diffusion components.
        - Palantir uses its own pseudotime to bias the KNN graph which is not implemented here.
        - Palantir uses a slightly different mechanism to ensure the graph remains connected when removing edges
          that point into the "pseudotime past".
    """

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
        Bias the connections by removing ones to past cells.

        It takes in symmetric connectivities and a pseudotime and removes edges that point "against" pseudotime,
        in this way creating a directed graph. For each node, it always keeps the closest neighbors,
        making sure the graph remains connected.

        Parameters
        ----------
        %(pt_scheme.parameters)s
        n_neighs
            Number of neighbors to keep.
        k
            Number, alongside with ``n_neighbors`` which determines the threshold for candidate indices.

        Returns
        -------
        %(pt_scheme.returns)s
        """
        k_thresh = max(0, min(30, int(np.floor(n_neighs / k)) - 1))

        # below code does not work with argpartition
        ixs = np.flip(np.argsort(neigh_dist))
        close_ixs, far_ixs = ixs[:k_thresh], ixs[k_thresh:]

        mask_keep = cell_pseudotime <= neigh_pseudotime[far_ixs]
        far_ixs_keep = far_ixs[mask_keep]

        biased_dist = np.zeros_like(neigh_dist)
        biased_dist[close_ixs] = neigh_dist[close_ixs]
        biased_dist[far_ixs_keep] = neigh_dist[far_ixs_keep]

        return biased_dist


class SoftThresholdScheme(ThresholdSchemeABC):
    """
    Thresholding scheme inspired by [VIA21]_.

    The idea is to downweight edges that points against the direction of increasing pseudotime. Essentially, the
    further "behind" a query cell is in pseudotime with respect to the current reference cell, the more penalized will
    be its graph-connectivity.
    """

    @d.dedent
    def __call__(
        self,
        cell_pseudotime: float,
        neigh_pseudotime: np.ndarray,
        neigh_dist: np.ndarray,
        b: float = 20.0,
        nu: float = 1.0,
        perc: Optional[int] = 95,
    ) -> np.ndarray:
        """
        Bias the connections by downweighting ones to past cells.

        This function uses `generalized logistic regression
        <https://en.wikipedia.org/wiki/Generalized_logistic_function>`_ to do the weighting.

        Parameters
        ----------
        %(pt_scheme.parameters)s
        %(soft_scheme)s

        Returns
        -------
        %(pt_scheme.returns)s
        """
        if perc is not None:
            neigh_dist = np.clip(
                neigh_dist,
                np.percentile(neigh_dist, 100 - perc),
                np.percentile(neigh_dist, perc),
            )

        past_ixs = np.where(neigh_pseudotime < cell_pseudotime)[0]
        if not len(past_ixs):
            return neigh_dist

        weights = np.ones_like(neigh_dist)

        dt = cell_pseudotime - neigh_pseudotime[past_ixs]
        weights[past_ixs] = 2.0 / ((1.0 + np.exp(b * dt)) ** (1.0 / nu))

        return neigh_dist * weights


class CustomThresholdScheme(ThresholdSchemeABC):
    """
    Class that wraps a user supplied scheme.

    Parameters
    ----------
    callback
        Function which returns the biased connectivities.
    """

    def __init__(
        self,
        callback: Callable[
            [float, np.ndarray, np.ndarray, np.ndarray, Any], np.ndarray
        ],
    ):
        super().__init__()
        self._callback = callback

    @d.dedent
    def __call__(
        self,
        cell_pseudotime: float,
        neigh_pseudotime: np.ndarray,
        neigh_dist: np.ndarray,
        **kwargs: Any,
    ) -> np.ndarray:
        """
        %(pt_scheme.summary)s

        Parameters
        ----------
        %(pt_scheme.parameters)s
        kwargs
            Additional keyword arguments.

        Returns
        -------
        %(pt_scheme.returns)s
        """  # noqa: D400
        return self._callback(cell_pseudotime, neigh_pseudotime, neigh_dist, **kwargs)
