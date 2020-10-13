# -*- coding: utf-8 -*-
import pytest

from cellrank.tl._linear_solver import _create_petsc_matrix


def _petsc_slepc_installed() -> True:
    try:
        import petsc4py
        import slepc4py

        return True
    except ImportError:
        return False


petsc_slepc_skip = pytest.mark.skipif(
    _petsc_slepc_installed(), reason="PETSc or SLEPc is not installed."
)


class TestLinearSolverScipy:
    def test_gmres(self):
        pass

    def test_direct_solver(self):
        pass


@petsc_slepc_skip
class TestLinearSolverPETSc:
    def test_create_petsc_matrix_from_sparse_not_as_dense(self):
        pass

    def test_create_petsc_matrix_from_sparse_as_dense(self):
        pass

    def test_create_petsc_matrix_from_sparse_not_csr(self):
        pass

    def test_create_solver_invalid_tolerance(self):
        pass

    def test_create_solver_invalid_solver(self):
        pass

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
