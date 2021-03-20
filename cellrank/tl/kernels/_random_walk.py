from typing import Any, List, Union, Optional, Sequence
from itertools import chain

import numpy as np
from scipy.sparse import issparse, spmatrix

from cellrank.ul._parallelize import parallelize


class RandomWalk:
    """TODO."""

    def __init__(
        self,
        transition_matrix: Union[np.ndarray, spmatrix],
        starting_ixs: Optional[Sequence[int]] = None,
        barrier: Optional[Sequence[int]] = None,
    ):
        if not np.allclose(transition_matrix.sum(1), 1.0):
            raise ValueError()

        self._tmat = transition_matrix
        self._ixs = np.arange(self._tmat.shape[0])
        self._is_sparse = issparse(self._tmat)
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

    def _sample(self, ix: int, *, rs: np.random.RandomState) -> int:
        return rs.choice(
            self._ixs,
            p=self._tmat[ix].A.squeeze() if self._is_sparse else self._tmat[ix],
        )

    def simulate_one(
        self,
        max_iter: Union[int, float] = 0.25,
        seed: Optional[int] = None,
        successive_hits: int = 0,
    ) -> np.ndarray:
        """TODO."""
        if isinstance(max_iter, float):
            max_iter = int(np.ceil(max_iter * len(self._ixs)))
        if max_iter <= 0:
            raise ValueError()
        if successive_hits < 0:
            raise ValueError()

        rs = np.random.RandomState(seed)
        ix = rs.choice(self._ixs, p=self._starting_dist)
        sim, cnt = [ix], -1

        for _ in range(max_iter):
            ix = self._sample(ix, rs=rs)
            sim.append(ix)
            cnt = (cnt + 1) if self._should_stop(ix) else -1
            if cnt >= successive_hits:
                break

        return np.array(sim)

    def _simulate_many(
        self,
        sims: np.ndarray,
        max_iter: Union[int, float] = 0.25,
        seed: Optional[int] = None,
        successive_hits: int = 0,
        queue: Optional[Any] = None,
    ) -> List[np.ndarray]:
        res = []
        for s in sims:
            sim = self.simulate_one(
                max_iter=max_iter,
                seed=None if seed is None else seed + s,
                successive_hits=successive_hits,
            )
            res.append(sim)
            if queue is not None:
                queue.put(1)

        queue.put(None)

        return res

    def simulate_many(
        self,
        n_sims: int,
        max_iter: Union[int, float] = 0.25,
        seed: Optional[int] = None,
        successive_hits: int = 0,
        n_jobs: Optional[int] = None,
        backend: str = "loky",
        show_progress_bar: bool = True,
    ) -> List[np.ndarray]:
        """TODO."""
        # TODO: logging
        simss = parallelize(
            self._simulate_many,
            collection=np.arange(n_sims),
            n_jobs=n_jobs,
            backend=backend,
            show_progress_bar=show_progress_bar,
            as_array=False,
            unit="sim",
        )(max_iter=max_iter, seed=seed, successive_hits=successive_hits)

        return list(chain.from_iterable(simss))
