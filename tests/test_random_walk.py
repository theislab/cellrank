import pytest


class TestRandomWalk:
    def test_matrix_not_row_stochastic(self):
        pass

    def test_starting_dist_does_not_sum_to_1(self):
        pass

    @pytest.mark("kind", ["sims", "iters", "hits"])
    def test_invalid_numbers(self):
        pass

    @pytest.mark.parametrize("n_sims", [1, 10])
    def test_reproducibility(self, n_sims: int):
        pass

    def test_simulate_one(self):
        pass

    def test_start_ixs(self):
        pass

    def test_stop_ixs(self):
        pass

    def successive_hits(self):
        pass

    @pytest.mark.parametrize("n_jobs", [1, 2])
    @pytest.mark.parametrize("backend", ["threading", "loky"])
    def test_simulate_many(self, backend: str, n_jobs: int):
        pass
