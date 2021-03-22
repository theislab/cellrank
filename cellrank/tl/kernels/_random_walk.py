from typing import Any, List, Union, Optional, Sequence
from itertools import chain

import numpy as np
from scipy.sparse import issparse, spmatrix

from cellrank import logging as logg
from cellrank.ul._docs import d
from cellrank.ul._parallelize import parallelize


class RandomWalk:
    """
    Class that simulates a random walk on a Markov chain.

    Parameters
    ----------
    transition_matrix
        Row-stochastic transition matrix.
    start_ixs
        Indices from which to uniformly sample the starting points. If `None`, use all points.
    stop_ixs
        Indices ...
        See ``successive_hits`` in :meth:`simulate_one`.
    """

    def __init__(
        self,
        transition_matrix: Union[np.ndarray, spmatrix],
        start_ixs: Optional[Sequence[int]] = None,
        stop_ixs: Optional[Sequence[int]] = None,
    ):
        if not np.allclose(transition_matrix.sum(1), 1.0):
            raise ValueError("Transition matrix is not row-stochastic.")

        self._tmat = transition_matrix
        self._ixs = np.arange(self._tmat.shape[0])
        self._is_sparse = issparse(self._tmat)
        self._stop_ixs = set([] if stop_ixs is None else stop_ixs)
        self._starting_dist = (
            np.ones_like(self._ixs)
            if start_ixs is None
            else np.isin(self._ixs, start_ixs)
        )
        self._starting_dist = self._starting_dist.astype(np.float64) / np.sum(
            self._starting_dist
        )

    def _should_stop(self, ix: int) -> bool:
        return ix in self._stop_ixs

    def _sample(self, ix: int, *, rs: np.random.RandomState) -> int:
        return rs.choice(
            self._ixs,
            p=self._tmat[ix].A.squeeze() if self._is_sparse else self._tmat[ix],
        )

    @d.get_sections(base="rw_sim", sections=["Parameters"])
    def simulate_one(
        self,
        max_iter: Union[int, float] = 0.25,
        seed: Optional[int] = None,
        successive_hits: int = 0,
    ) -> np.ndarray:
        """
        Simulate one random walk.

        Parameters
        ----------
        max_iter
            Maximum number of steps of a random walk.
        seed
            Random seed.
        successive_hits
            Number of successive hits in the ``stop_ixs`` required to stop prematurely.

        Returns
        -------
        Array of shape ``(max_iter,)`` of states that have been visited. If ``stop_ixs`` was specified, the array
        may have smaller shape.
        """
        if isinstance(max_iter, float):
            max_iter = int(np.ceil(max_iter * len(self._ixs)))
        if max_iter <= 1:
            raise ValueError(
                f"Expected number of iteration to be > 1, found `{max_iter}`."
            )
        if successive_hits < 0:
            raise ValueError(
                f"Expected number of successive hits to be positive, found `{successive_hits}`."
            )

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

    @d.dedent
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
        """
        Simulate many random walks.

        Parameters
        ----------
        n_sims
            Number of random walks to simulate.
        %(rw_sim.params)s
        %(parallel)s

        Returns
        -------
        List of arrays of shape ``(max_iter,)`` of states that have been visited. If ``stop_ixs`` was specified,
        the arrays may have smaller shape.
        """
        if n_sims <= 0:
            raise ValueError(f"Expected `n_sims` to be positive, found `{n_sims}`.")
        start = logg.info(
            f"Simulating `{n_sims}` random walks of maximum length `{max_iter}`"
        )

        simss = parallelize(
            self._simulate_many,
            collection=np.arange(n_sims),
            n_jobs=n_jobs,
            backend=backend,
            show_progress_bar=show_progress_bar,
            as_array=False,
            unit="sim",
        )(max_iter=max_iter, seed=seed, successive_hits=successive_hits)
        simss = list(chain.from_iterable(simss))

        logg.info("    Finish", time=start)

        return simss
