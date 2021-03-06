from abc import ABC, abstractmethod
from typing import Any, Callable

import numpy as np

from cellrank.ul._docs import d


class PseudotimeScheme(ABC):
    """TODO."""

    @d.get_sections("pt_scheme", sections=["Parameters", "Returns"])
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


class HardThresholdScheme(PseudotimeScheme):
    """TODO."""

    @d.dedent
    def __call__(
        self,
        cell_pseudotime: float,
        neigh_pseudotime: np.ndarray,
        neigh_dist: np.ndarray,
        n_neigh: int,
        k: int = 3,
    ) -> np.ndarray:
        """
        TODO.

        Parameters
        ----------
        %(pt_scheme.parameters)s
        n_neigh
            Number of neighbors to keep.
        k
            Number, alongside with ``n_neighbors`` which determined the threshold for candidate indices.

        Returns
        -------
        %(pt_scheme.returns)

        Notes
        -----
        It is up to the implementor of this function to ensure that the graph stays connected after removing some of
        the cell's neighbors.
        """


class SoftThresholdScheme(PseudotimeScheme):
    """TODO."""

    @d.dedent
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
        TODO: remove kwargs, determine which args we need.

        Returns
        -------
        %(pt_scheme.returns)
        """


class CustomScheme(PseudotimeScheme):
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
        %(pt_scheme.returns)
        """
        return self._callback(cell_pseudotime, neigh_pseudotime, neigh_dist, **kwargs)
