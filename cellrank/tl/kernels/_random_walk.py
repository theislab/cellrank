from typing import Any, List, Union, Optional, Sequence

import numpy as np
from scipy.sparse import issparse, spmatrix


class RandomWalk:
    """TODO."""

    def __init__(
        self,
        transition_matrix: Union[np.ndarray, spmatrix],
        max_iter: int = 1000,
        starting_ixs: Optional[Sequence[int]] = None,
        barrier: Optional[Sequence[int]] = None,
    ):
        self._tmat = transition_matrix
        self._ixs = np.arange(self._tmat.shape[0])
        self._is_sparse = issparse(self._tmat)
        self._max_iter = max_iter
        self._barrier = set([] if barrier is None else barrier)
        self._starting_dist = (
            np.ones_like(self._ixs)
            if starting_ixs is None
            else np.isin(self._ixs, starting_ixs)
        )
        self._starting_dist = self._starting_dist.astype(np.float64) / np.sum(
            self._starting_dist
        )

    def _should_stop(self, ix: int) -> bool:
        return ix in self._barrier

    def _sample(self, ix: int) -> int:
        return np.random.choice(
            self._ixs,
            p=self._tmat[ix].A.squeeze() if self._is_sparse else self._tmat[ix],
        )

    def simulate_one(
        self, seed: Optional[int] = None, threshold: int = 0
    ) -> np.ndarray:
        """TODO."""
        rs = np.random.RandomState(seed)

        ix = rs.choice(self._ixs, p=self._starting_dist)
        sim, cnt = [ix], -1

        for _ in range(self._max_iter):
            ix = self._sample(ix)
            sim.append(ix)
            cnt = (cnt + 1) if self._should_stop(ix) else -1
            if cnt >= threshold:
                break

        return np.array(sim)

    def simulate_many(
        self, n_sims: int, seed: Optional[int] = None, threshold: int = 0, **kwargs: Any
    ) -> List[np.ndarray]:
        """TODO."""
