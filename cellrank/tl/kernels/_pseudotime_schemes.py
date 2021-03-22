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
        neigh_conn: np.ndarray,
        **kwargs: Any,
    ) -> np.ndarray:
        """
        Calculate biased connections for a given cell.

        Parameters
        ----------
        cell_pseudotime
            Pseudotime of the current cell.
        neigh_pseudotime
            Array of shape ``(n_neighbors,)`` containing pseudotimes of neighbors.
        neigh_conn
            Array of shape ``(n_neighbors,)`` containing connectivities of the current cell and its neighbors.

        Returns
        -------
        Array of shape ``(n_neighbors,)`` containing the biased connectivities.
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
            Pseudotemporal ordering of cells.

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
        neigh_conn: np.ndarray,
        n_neighs: int,
        frac_to_keep: float = 0.3,
    ) -> np.ndarray:
        """
        Convert the undirected graph of cell-cell similarities into a directed one by removing "past" edges.

        This uses a pseudotemporal measure to remove graph-edges that point into the pseudotime-past. For each cell,
        it keeps the closest neighbors, even if they are in the pseudotime past, to make sure the graph remains
        connected.

        Parameters
        ----------
        %(pt_scheme.parameters)s
        n_neighs
            Number of neighbors to keep.
        frac_to_keep
            The `fract_to_keep` * n_neighbors closest neighbors (according to graph connectivities) are kept, no matter
            whether they lie in the pseudotemporal past or future.

        Returns
        -------
        %(pt_scheme.returns)s
        """
        k_thresh = max(0, min(30, int(np.floor(n_neighs * frac_to_keep)) - 1))

        # below code does not work with argpartition
        ixs = np.flip(np.argsort(neigh_conn))
        close_ixs, far_ixs = ixs[:k_thresh], ixs[k_thresh:]

        mask_keep = cell_pseudotime <= neigh_pseudotime[far_ixs]
        far_ixs_keep = far_ixs[mask_keep]

        biased_conn = np.zeros_like(neigh_conn)
        biased_conn[close_ixs] = neigh_conn[close_ixs]
        biased_conn[far_ixs_keep] = neigh_conn[far_ixs_keep]

        return biased_conn


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
        neigh_conn: np.ndarray,
        b: float = 20.0,
        nu: float = 1.0,
        perc: Optional[int] = 95,
    ) -> np.ndarray:
        """
        Bias the connectivities by downweighting ones to past cells.

        This function uses `generalized logistic regression
        <https://en.wikipedia.org/wiki/Generalized_logistic_function>`_ to weight the past connectivities.

        Parameters
        ----------
        %(pt_scheme.parameters)s
        %(soft_scheme)s

        Returns
        -------
        %(pt_scheme.returns)s
        """
        if perc is not None:
            neigh_conn = np.clip(
                neigh_conn,
                np.percentile(neigh_conn, 100 - perc),
                np.percentile(neigh_conn, perc),
            )

        past_ixs = np.where(neigh_pseudotime < cell_pseudotime)[0]
        if not len(past_ixs):
            return neigh_conn

        weights = np.ones_like(neigh_conn)

        dt = cell_pseudotime - neigh_pseudotime[past_ixs]
        weights[past_ixs] = 2.0 / ((1.0 + np.exp(b * dt)) ** (1.0 / nu))

        return neigh_conn * weights


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
        neigh_conn: np.ndarray,
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
        return self._callback(cell_pseudotime, neigh_pseudotime, neigh_conn, **kwargs)
