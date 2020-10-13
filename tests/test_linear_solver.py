# -*- coding: utf-8 -*-
import pytest

import numpy as np
from scipy.sparse import eye as speye
from scipy.sparse import random, csr_matrix

from cellrank.tl._linear_solver import _solve_lin_system, _create_petsc_matrix


def _petsc_not_installed() -> bool:
    try:
        import petsc4py
        import slepc4py

        return False
    except ImportError:
        return True


petsc_slepc_skip = pytest.mark.skipif(
    _petsc_not_installed(), reason="PETSc or SLEPc is not installed."
)


class TestLinearSolverScipy:
    def test_invalid_solver(self):
        np.random.seed(42)
        A = np.random.normal(size=(20, 10))
        B = np.random.normal(size=(20, 10))

        with pytest.raises(ValueError):
            _solve_lin_system(A, B, solver="foobar", use_petsc=False)

    @pytest.mark.parametrize("seed", range(5))
    def test_gmres_dense(self, seed: int):
        np.random.seed(seed)
        A = np.random.normal(size=(20, 20))
        B = np.random.normal(size=(20, 10))

        sol = _solve_lin_system(
            A,
            B,
            solver="gmres",
            use_petsc=False,
            use_eye=False,
            show_progress_bar=False,
            tol=1e-6,
        )
        for i in range(B.shape[1]):
            np.testing.assert_allclose(A @ sol[:, i], B[:, i], rtol=1e-6, atol=1e-10)

    @pytest.mark.parametrize("seed", range(5, 10))
    def test_eye_dense(self, seed: int):
        np.random.seed(seed)
        A = np.random.normal(size=(20, 20))
        B = np.random.normal(size=(20, 10))
        I = np.eye(20, 20)

        sol = _solve_lin_system(
            A,
            B,
            solver="gmres",
            use_petsc=False,
            use_eye=True,
            show_progress_bar=False,
            tol=1e-6,
        )
        for i in range(B.shape[1]):
            np.testing.assert_allclose(
                (I - A) @ sol[:, i], B[:, i], rtol=1e-6, atol=1e-10
            )

    @pytest.mark.parametrize("seed", range(10, 15))
    def test_direct_solver_dense(self, seed: int):
        np.random.seed(seed)
        A = np.random.normal(size=(20, 20))
        B = np.random.normal(size=(20, 10))

        sol = _solve_lin_system(
            A,
            B,
            solver="direct",
            use_petsc=False,
            use_eye=False,
            show_progress_bar=False,
            tol=1e-6,
        )
        for i in range(B.shape[1]):
            np.testing.assert_allclose(A @ sol[:, i], B[:, i], rtol=1e-6, atol=1e-10)

    @pytest.mark.parametrize("seed", range(5))
    def test_gmres_sparse(self, seed: int):
        np.random.seed(seed)
        A = random(
            20, 20, density=0.5, random_state=np.random.randint(0, 100), format="csr"
        )
        B = random(
            20, 10, density=0.3, random_state=np.random.randint(0, 100), format="csr"
        )

        sol = _solve_lin_system(
            A,
            B,
            solver="gmres",
            use_petsc=False,
            use_eye=False,
            show_progress_bar=False,
            tol=1e-6,
        )
        for i in range(B.shape[1]):
            np.testing.assert_allclose(
                A.A @ sol[:, i], B.A[:, i], rtol=1e-6, atol=1e-10
            )

    @pytest.mark.parametrize("seed", range(5, 10))
    def test_eye_sparse(self, seed: int):
        np.random.seed(seed)
        A = random(
            20, 20, density=0.5, random_state=np.random.randint(0, 100), format="coo"
        )
        B = random(
            20, 10, density=0.3, random_state=np.random.randint(0, 100), format="csc"
        )
        I = speye(20, 20)

        sol = _solve_lin_system(
            A,
            B,
            solver="gmres",
            use_petsc=False,
            use_eye=True,
            show_progress_bar=False,
            tol=1e-6,
        )
        for i in range(B.shape[1]):
            np.testing.assert_allclose(
                (I - A) @ sol[:, i], B.A[:, i], rtol=1e-6, atol=1e-10
            )

    @pytest.mark.parametrize("seed", range(10, 15))
    def test_direct_solver_sparse(self, seed: int):
        np.random.seed(seed)
        A = random(
            20, 20, density=0.5, random_state=np.random.randint(0, 100), format="csr"
        )
        B = random(
            20, 10, density=0.3, random_state=np.random.randint(0, 100), format="csr"
        )

        sol = _solve_lin_system(
            A,
            B,
            solver="direct",
            use_petsc=False,
            use_eye=False,
            show_progress_bar=False,
            tol=1e-6,
        )
        for i in range(B.shape[1]):
            np.testing.assert_allclose(
                A.A @ sol[:, i], B.A[:, i], rtol=1e-6, atol=1e-10
            )


@petsc_slepc_skip
class TestLinearSolverPETSc:
    def test_create_petsc_matrix_no_a_matrix(self):
        with pytest.raises(TypeError):
            _create_petsc_matrix(np.empty((100,)))

    def test_create_petsc_matrix_from_dense(self):
        x = np.random.normal(size=(10, 2))
        res = _create_petsc_matrix(x)

        assert res.assembled
        assert res.getType() == "seqdense"
        np.testing.assert_array_equal(res.getDenseArray(), x)

    def test_create_petsc_matrix_from_sparse_as_dense(self):
        x = random(100, 10, format="csr", density=0.1)
        res = _create_petsc_matrix(x, as_dense=True)

        assert res.assembled
        assert res.getType() == "seqdense"
        np.testing.assert_array_equal(res.getDenseArray(), x.A)

    def test_create_petsc_matrix_from_sparse_as_not_dense(self):
        x = random(100, 10, format="csr", density=0.1)
        res = _create_petsc_matrix(x, as_dense=False)

        assert res.assembled
        assert res.getType() == "seqaij"
        np.testing.assert_array_equal(csr_matrix(res.getValuesCSR()[::-1]).A, x.A)

    def test_create_petsc_matrix_from_sparse_not_csr(self):
        x = random(100, 10, format="coo", density=0.1)
        res = _create_petsc_matrix(x, as_dense=False)

        assert res.assembled
        assert res.getType() == "seqaij"

        np.testing.assert_array_equal(csr_matrix(res.getValuesCSR()[::-1]).A, x.A)

    def test_create_solver_invalid_solver(self):
        from petsc4py.PETSc import Error

        np.random.seed(42)
        A = np.random.normal(size=(20, 10))
        B = np.random.normal(size=(20, 10))

        with pytest.raises(Error):
            _solve_lin_system(A, B, solver="foobar", use_petsc=True)

    def test_create_solver_invalid_preconditioner(self):
        pass

    def test_many_solves_sparse_matrix(self):
        pass

    def test_many_solves_array(self):
        pass

    def test_mat_mat_solve_invalid_dimension(self):
        pass

    def test_mat_mat_solve_normal_run(self):
        pass

    def test_invert_matrix_not_square(self):
        pass

    def test_invert_matrix_normal_run(self):
        pass

    def test_petsc_scipy_direct_matches(self):
        pass

    def test_petsc_scipy_gmres_matches(self):
        pass
