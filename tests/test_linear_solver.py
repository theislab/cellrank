import pytest

import numpy as np
import scipy.sparse as sp

from cellrank._utils._linear_solver import (
    _create_petsc_matrix,
    _is_petsc_slepc_available,
    _petsc_direct_solve,
    _solve_lin_system,
)

petsc_slepc_skip = pytest.mark.skipif(not _is_petsc_slepc_available(), reason="PETSc or SLEPc is not installed.")


def _create_a_b_matrices(seed: int, sparse: bool) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    if sparse:
        A = sp.random(20, 20, density=0.8, random_state=rng.integers(0, 100), format="csr")
        B = sp.random(20, 10, density=1, random_state=rng.integers(0, 100), format="csr")
    else:
        A = rng.normal(size=(20, 20))
        B = rng.normal(size=(20, 10))

    return A, B


class TestScipyLinearSolver:
    def test_invalid_solver(self):
        rng = np.random.default_rng(42)
        A = rng.normal(size=(20, 10))
        B = rng.normal(size=(20, 10))

        with pytest.raises(ValueError, match="Invalid solver"):
            _solve_lin_system(A, B, solver="foobar", use_petsc=False)

    @pytest.mark.parametrize(("seed", "sparse"), zip(range(10), [False] * 5 + [True] * 5))
    def test_gmres(self, seed: int, sparse: bool):
        A, B = _create_a_b_matrices(seed, sparse)

        sol = _solve_lin_system(
            A,
            B,
            solver="gmres",
            use_petsc=False,
            use_eye=False,
            show_progress_bar=False,
            tol=1e-6,
        )
        assert sol.ndim == 2

        if sparse:
            A = A.toarray()
            B = B.toarray()

        np.testing.assert_allclose(A @ sol, B, rtol=1e-6, atol=1e-10)

    @pytest.mark.parametrize(("seed", "sparse"), zip(range(10, 20), [False] * 5 + [True] * 5))
    def test_direct_solver_dense(self, seed: int, sparse: bool):
        A, B = _create_a_b_matrices(seed, sparse)

        sol = _solve_lin_system(
            A,
            B,
            solver="direct",
            use_petsc=False,
            use_eye=False,
            show_progress_bar=False,
            tol=1e-6,
        )
        assert sol.ndim == 2

        if sparse:
            A = A.toarray()
            B = B.toarray()

        np.testing.assert_allclose(A @ sol, B, rtol=1e-6, atol=1e-10)

    @pytest.mark.parametrize(("seed", "sparse"), zip(range(30, 40), [False] * 5 + [True] * 5))
    def test_eye(self, seed: int, sparse: bool):
        A, B = _create_a_b_matrices(seed, sparse)

        sol = _solve_lin_system(
            A,
            B,
            solver="gmres",
            use_petsc=False,
            use_eye=True,
            show_progress_bar=False,
            tol=1e-6,
        )
        assert sol.ndim == 2

        if sparse:
            A = A.toarray()
            B = B.toarray()
        A = np.eye(20, 20) - A

        np.testing.assert_allclose(A @ sol, B, rtol=1e-6, atol=1e-10)


@petsc_slepc_skip
class TestLinearSolverPETSc:
    def test_create_petsc_matrix_no_a_matrix(self):
        with pytest.raises(TypeError, match="integer is required"):
            _create_petsc_matrix(np.empty((100,)))

    def test_create_petsc_matrix_from_dense(self):
        rng = np.random.default_rng(42)
        x = rng.normal(size=(10, 2))
        res = _create_petsc_matrix(x)

        assert res.assembled
        assert res.getType() == "seqdense"
        np.testing.assert_array_equal(res.getDenseArray(), x)

    def test_create_petsc_matrix_from_sparse_as_dense(self):
        x = sp.random(100, 10, format="csr", density=0.1)
        res = _create_petsc_matrix(x, as_dense=True)

        assert res.assembled
        assert res.getType() == "seqdense"
        np.testing.assert_array_equal(res.getDenseArray(), x.toarray())

    def test_create_petsc_matrix_from_sparse_as_not_dense(self):
        x = sp.random(100, 10, format="csr", density=0.1)
        res = _create_petsc_matrix(x, as_dense=False)

        assert res.assembled
        assert res.getType() == "seqaij"
        np.testing.assert_array_equal(sp.csr_matrix(res.getValuesCSR()[::-1]).toarray(), x.toarray())

    def test_create_petsc_matrix_from_sparse_not_csr(self):
        x = sp.random(100, 10, format="coo", density=0.1)
        res = _create_petsc_matrix(x, as_dense=False)

        assert res.assembled
        assert res.getType() == "seqaij"

        np.testing.assert_array_equal(sp.csr_matrix(res.getValuesCSR()[::-1]).toarray(), x.toarray())

    def test_create_solver_invalid_solver(self):
        from petsc4py.PETSc import Error

        rng = np.random.default_rng(42)

        A = rng.normal(size=(20, 20))
        B = rng.normal(size=(20, 10))

        with pytest.raises(Error, match=r".*foobar"):
            _solve_lin_system(A, B, solver="foobar", use_petsc=True)

    def test_create_solver_invalid_preconditioner(self):
        from petsc4py.PETSc import Error

        rng = np.random.default_rng(42)

        A = rng.normal(size=(20, 20))
        B = rng.normal(size=(20, 10))

        with pytest.raises(Error, match=r".*foobar"):
            _solve_lin_system(A, B, preconditioner="foobar", use_petsc=True)

    def test_solve_invalid_dimension(self):
        from petsc4py.PETSc import Error

        rng = np.random.default_rng()

        A = rng.normal(size=(20, 10))
        B = rng.normal(size=(20, 10))

        with pytest.raises(Error, match=r".* square"):
            _solve_lin_system(A, B, use_petsc=True)

    @pytest.mark.parametrize(
        ("seed", "solver", "sparse"),
        zip(
            range(42, 62),
            ["direct"] * 5 + ["gmres"] * 5 + ["direct"] * 5 + ["gmres"] * 5,
            [False] * 10 + [True] * 10,
        ),
    )
    def test_petsc_scipy_matches(self, seed: int, solver: str, sparse: bool):
        A, B = _create_a_b_matrices(seed, sparse)

        sol_petsc = _solve_lin_system(
            A,
            B,
            solver=solver,
            use_petsc=True,
            use_eye=True,
            show_progress_bar=False,
            preconditioner="lu",
            tol=1e-6,
        )
        sol_scipy = _solve_lin_system(
            A,
            B,
            solver=solver,
            use_petsc=False,
            use_eye=True,
            show_progress_bar=False,
            tol=1e-6,
        )
        assert sol_petsc.ndim == 2
        assert sol_scipy.ndim == 2
        if sparse:
            A = A.toarray()
            B = B.toarray()
        A = np.eye(A.shape[0]) - A

        np.testing.assert_allclose(A @ sol_scipy, B, rtol=1e-6, atol=1e-8)
        np.testing.assert_allclose(A @ sol_petsc, B, rtol=1e-6, atol=1e-8)
        np.testing.assert_allclose(sol_petsc, sol_scipy, rtol=1e-6, atol=1e-8)

    @pytest.mark.parametrize(("seed", "sparse"), zip(range(10), [False] * 5 + [True] * 5))
    def test_gmres(self, seed: int, sparse: bool):
        A, B = _create_a_b_matrices(seed, sparse)

        sol = _solve_lin_system(
            A,
            B,
            solver="gmres",
            use_petsc=True,
            use_eye=sparse,
            show_progress_bar=False,
            preconditioner="lu",
            tol=1e-6,
        )
        assert sol.ndim == 2

        if sparse:
            A = A.toarray()
            A = np.eye(A.shape[0]) - A
            B = B.toarray()

        np.testing.assert_allclose(A @ sol, B, rtol=1e-6, atol=1e-8)

    @pytest.mark.parametrize(("seed", "sparse"), zip(range(10, 20), [False] * 5 + [True] * 5))
    def test_direct_solver(self, seed: int, sparse: bool):
        A, B = _create_a_b_matrices(seed, sparse)

        sol = _solve_lin_system(
            A,
            B,
            solver="direct",
            use_petsc=True,
            use_eye=True,
            show_progress_bar=False,
            tol=1e-6,
        )
        assert sol.ndim == 2

        if sparse:
            A = A.toarray()
            B = B.toarray()
        A = np.eye(A.shape[0]) - A

        np.testing.assert_allclose(A @ sol, B, rtol=1e-6, atol=1e-10)

    @pytest.mark.parametrize(("seed", "sparse"), zip(range(10, 20), [False] * 5 + [True] * 5))
    def test_mat_solve(self, seed: int, sparse: bool):
        A, B = _create_a_b_matrices(seed, sparse)
        A = (sp.eye(*A.shape) if sparse else np.eye(*A.shape)) - A

        X = _petsc_direct_solve(A, B, tol=1e-8)
        assert X.ndim == 2
        if sparse:
            A = A.toarray()
            B = B.toarray()

        np.testing.assert_allclose(A @ X, B, rtol=1e-6)

    @pytest.mark.parametrize(("seed", "sparse"), zip(range(10, 20), [False] * 5 + [True] * 5))
    def test_mat_solve_1_dim_b(self, seed: int, sparse: bool):
        A, B = _create_a_b_matrices(seed, sparse)
        A = (sp.eye(*A.shape) if sparse else np.eye(*A.shape)) - A
        B = B[:, 0]

        X = _petsc_direct_solve(A, B, tol=1e-12)
        assert X.ndim == 2
        assert X.shape[1] == 1
        if sparse:
            A = A.toarray()
            B = B.toarray()

        np.testing.assert_allclose(((A @ X).squeeze()), B.squeeze(), rtol=1e-6, atol=1e-6)

    @pytest.mark.parametrize(("seed", "sparse"), zip(range(10, 20), [False] * 5 + [True] * 5))
    def test_mat_solve_inversion(self, seed: int, sparse: bool):
        A, _ = _create_a_b_matrices(seed, sparse)
        A = (sp.eye(*A.shape) if sparse else np.eye(*A.shape)) - A

        X = _petsc_direct_solve(A, tol=1e-8)
        assert X.ndim == 2
        if sparse:
            A = A.toarray()

        np.testing.assert_allclose(A @ X, np.eye(*A.shape), rtol=1e-6, atol=1e-6)
