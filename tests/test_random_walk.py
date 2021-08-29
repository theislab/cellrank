import pytest

from cellrank.tl.kernels._random_walk import RandomWalk

import numpy as np


class TestRandomWalk:
    def test_matrix_not_row_stochastic(self, test_matrix_4: np.ndarray):
        with pytest.raises(
            ValueError, match=r"Transition matrix is not row-stochastic."
        ):
            _ = RandomWalk(test_matrix_4)

    def test_starting_dist_does_not_sum_to_1(self, test_matrix_1: np.ndarray):
        with pytest.raises(
            ValueError, match=r"No starting indices have been selected."
        ):
            _ = RandomWalk(test_matrix_1, start_ixs=[1000])

    @pytest.mark.parametrize("kind", ["simulations", "iterations", "hits"])
    def test_invalid_numbers(self, test_matrix_1: np.ndarray, kind: str):
        rw = RandomWalk(test_matrix_1)
        with pytest.raises(ValueError, match=kind):
            if kind == "simulations":
                rw.simulate_many(n_sims=0)
            elif kind == "iterations":
                rw.simulate_one(max_iter=1)
            elif kind == "hits":
                rw.simulate_one(successive_hits=-1)
            else:
                raise NotImplementedError(kind)

    @pytest.mark.parametrize("n_sims", [1, 10])
    def test_reproducibility(self, test_matrix_1: np.ndarray, n_sims: int):
        rw1 = RandomWalk(test_matrix_1)
        rw2 = RandomWalk(test_matrix_1)

        if n_sims == 1:
            r1, r2 = rw1.simulate_one(seed=42), rw2.simulate_one(seed=42)
        else:
            r1, r2 = rw1.simulate_many(10, seed=42), rw2.simulate_many(10, seed=42)

        np.testing.assert_array_equal(r1, r2)

    def test_simulate_one(self, test_matrix_1: np.ndarray):
        res = RandomWalk(test_matrix_1).simulate_one(max_iter=10)

        assert isinstance(res, np.ndarray)
        assert np.issubdtype(res.dtype, np.integer)
        np.testing.assert_array_equal(res.shape, (11,))

    def test_start_ixs(self, test_matrix_1: np.ndarray):
        res = RandomWalk(test_matrix_1, start_ixs=[1]).simulate_one(max_iter=10)

        assert res[0] == 1

    def test_stop_ixs(self, test_matrix_1: np.ndarray):
        res = RandomWalk(test_matrix_1, stop_ixs=[1]).simulate_one(
            max_iter=1000, seed=42
        )

        assert len(res) <= 1001
        assert res[-1] == 1
        assert res[-2] != 1

    def test_successive_hits(self, test_matrix_1: np.ndarray):
        test_matrix_1[0, 0] = 0.2
        test_matrix_1[0, 2] = 0.0
        res = RandomWalk(test_matrix_1, stop_ixs=[0]).simulate_one(
            max_iter=1000, seed=42, successive_hits=1
        )

        assert len(res) <= 1001
        assert res[-1] == 0
        assert res[-2] == 0
        assert res[-3] != 0

    @pytest.mark.parametrize("n_jobs", [1, 2])
    @pytest.mark.parametrize("backend", ["threading", "loky"])
    def test_simulate_many(self, test_matrix_1: np.ndarray, backend: str, n_jobs: int):
        res = RandomWalk(test_matrix_1).simulate_many(
            10, max_iter=100, seed=42, backend=backend, n_jobs=n_jobs
        )

        assert isinstance(res, list)
        assert len(res) == 10
        for i, r in enumerate(res):
            assert isinstance(r, np.ndarray), i
            assert np.issubdtype(r.dtype, np.integer), i
            np.testing.assert_array_equal(r.shape, (101,))
