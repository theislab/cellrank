# -*- coding: utf-8 -*-
import unittest
import numpy as np

from cellrank.tools._constants import Direction, _transition
from cellrank.tools.kernels import VelocityKernel, ConnectivityKernel, PalantirKernel
from cellrank.tools.kernels._kernel import Constant, KernelAdd, KernelMul, _is_bin_mult
from cellrank.tools._utils import _normalize
from _test_helper import (
    transition_matrix,
    bias_knn,
    density_normalization,
    create_dummy_adata,
    create_kernels,
)

np.random.seed(42)
_adata = create_dummy_adata(50)
_rtol = 1e-14


class KernelInitializeTestCase(unittest.TestCase):
    def test_none_transition_matrix(self):
        adata = _adata.copy()
        vk = VelocityKernel(adata)
        ck = ConnectivityKernel(adata)
        pk = PalantirKernel(adata, time_key="latent_time")

        self.assertIsNone(vk._transition_matrix)
        self.assertIsNone(ck._transition_matrix)
        self.assertIsNone(pk._transition_matrix)

    def test_not_none_transition_matrix_compute(self):
        adata = _adata.copy()
        vk = VelocityKernel(adata).compute_transition_matrix()
        ck = ConnectivityKernel(adata).compute_transition_matrix()
        pk = PalantirKernel(adata, time_key="latent_time").compute_transition_matrix()

        self.assertIsNotNone(vk._transition_matrix)
        self.assertIsNotNone(ck._transition_matrix)
        self.assertIsNotNone(pk._transition_matrix)

    def test_not_none_transition_matrix_accessor(self):
        adata = _adata.copy()
        vk = VelocityKernel(adata)
        ck = ConnectivityKernel(adata)
        pk = PalantirKernel(adata, time_key="latent_time")

        self.assertIsNotNone(vk.transition_matrix)
        self.assertIsNotNone(ck.transition_matrix)
        self.assertIsNotNone(pk.transition_matrix)

    def test_adding_hidden_constants(self):
        adata = _adata.copy()
        k: KernelAdd = VelocityKernel(adata) + ConnectivityKernel(adata)

        self.assertTrue(_is_bin_mult(k[0]))
        self.assertIsInstance(k[0], KernelMul)
        self.assertIsInstance(k[0][0], Constant)
        self.assertIsInstance(k[0][1], VelocityKernel)
        self.assertEqual(k[0][0].transition_matrix, 1.0)

        self.assertTrue(_is_bin_mult(k[1]))
        self.assertIsInstance(k[1], KernelMul)
        self.assertIsInstance(k[1][0], Constant)
        self.assertIsInstance(k[1][1], ConnectivityKernel)
        self.assertEqual(k[1][0].transition_matrix, 1.0)

    def test_length(self):
        adata = _adata.copy()
        k: KernelAdd = VelocityKernel(adata) + ConnectivityKernel(adata)
        self.assertTrue(len(k) == 2)

    def test_accessor_out_of_range(self):
        adata = _adata.copy()
        k: KernelAdd = VelocityKernel(adata) + ConnectivityKernel(adata)

        with self.assertRaises(IndexError):
            _ = k[2]

    def test_parent(self):
        adata = _adata.copy()
        vk = VelocityKernel(adata)
        ck = ConnectivityKernel(adata)
        k = vk + ck

        self.assertIs(vk._parent._parent, k)  # invisible constants
        self.assertIs(ck._parent._parent, k)
        self.assertIs(k._parent, None)

    def test_uninitialized_both(self):
        adata = _adata.copy()
        k = VelocityKernel(adata) + ConnectivityKernel(adata)

        with self.assertRaises(RuntimeError):
            k.compute_transition_matrix()

    def test_uninitialized_one(self):
        adata = _adata.copy()
        k = (
            VelocityKernel(adata)
            + ConnectivityKernel(adata).compute_transition_matrix()
        )

        with self.assertRaises(RuntimeError):
            k.compute_transition_matrix()

    def test_initialized(self):
        adata = _adata.copy()
        k = (
            VelocityKernel(adata).compute_transition_matrix()
            + ConnectivityKernel(adata).compute_transition_matrix()
        )
        k.compute_transition_matrix()

        self.assertIsNotNone(k.transition_matrix)

    def test_invalida_type(self):
        adata = _adata.copy()

        with self.assertRaises(TypeError):
            _ = None * VelocityKernel(adata)

    def test_negative_constant(self):
        adata = _adata.copy()

        with self.assertRaises(ValueError):
            _ = -1 * VelocityKernel(adata)

    def test_invalid_constant(self):
        with self.assertRaises(TypeError):
            _ = Constant(None, None)

    def test_inversion(self):
        adata = _adata.copy()
        c = ConnectivityKernel(adata, backward=False)
        self.assertFalse(c.backward)

        nc = ~c
        self.assertTrue(nc.backward)

    def test_inversion_inplace(self):
        adata = _adata.copy()
        c = ConnectivityKernel(adata, backward=False)

        self.assertFalse(c.backward)
        _ = ~c
        self.assertTrue(c.backward)

    def test_inversion_propagation(self):
        adata = _adata.copy()
        c = ConnectivityKernel(adata, backward=False)
        v = VelocityKernel(adata, backward=False)
        k = ~(c + v)

        self.assertTrue(c.backward)
        self.assertTrue(v.backward)
        self.assertTrue(k.backward)

    def test_inversion_recalculation(self):
        adata = _adata.copy()
        c = ConnectivityKernel(adata).compute_transition_matrix()
        z = ~(c + c)
        with self.assertRaises(RuntimeError):
            z.compute_transition_matrix()

    def test_inversion_preservation_of_constants(self):
        adata = _adata.copy()
        c = ConnectivityKernel(adata).compute_transition_matrix()
        a = (3 * c + 1 * c).compute_transition_matrix()
        b = ~a
        c.compute_transition_matrix()

        self.assertEqual(a[0][0].transition_matrix, 3 / 4)
        self.assertEqual(b[0][0].transition_matrix, 3 / 4)
        self.assertEqual(a[1][0].transition_matrix, 1 / 4)
        self.assertEqual(b[1][0].transition_matrix, 1 / 4)

    def test_addition_simple(self):
        adata = _adata.copy()
        k = VelocityKernel(adata) + ConnectivityKernel(adata)

        self.assertIsInstance(k, KernelAdd)

    def test_multiplication_simple(self):
        adata = _adata.copy()
        k = 10 * VelocityKernel(adata)
        c = _is_bin_mult(k)

        self.assertIsInstance(c, Constant)
        self.assertEqual(c.transition_matrix, 10)

    def test_multiplication_simple_normalization(self):
        adata = _adata.copy()
        k = 10 * VelocityKernel(adata).compute_transition_matrix()
        c = _is_bin_mult(k)

        self.assertEqual(c.transition_matrix, 10)

    def test_constant(self):
        adata = _adata.copy()
        k = 9 * VelocityKernel(adata) + 1 * ConnectivityKernel(adata)
        c1, c2 = _is_bin_mult(k[0]), _is_bin_mult(k[1])

        self.assertEqual(c1.transition_matrix, 9)
        self.assertEqual(c2.transition_matrix, 1)

    def test_constant_normalize_2(self):
        adata = _adata.copy()
        k = (
            9 * VelocityKernel(adata).compute_transition_matrix()
            + 1 * ConnectivityKernel(adata).compute_transition_matrix()
        )
        k.compute_transition_matrix()
        c1, c2 = _is_bin_mult(k[0]), _is_bin_mult(k[1])

        self.assertEqual(c1.transition_matrix, 9 / 10)
        self.assertEqual(c2.transition_matrix, 1 / 10)

    def test_constant_normalize_3(self):
        adata = _adata.copy()
        k = (
            VelocityKernel(adata).compute_transition_matrix()
            + ConnectivityKernel(adata).compute_transition_matrix()
            + ConnectivityKernel(adata).compute_transition_matrix()
        )
        k.compute_transition_matrix()
        c1, c2, c3 = _is_bin_mult(k[0]), _is_bin_mult(k[1]), _is_bin_mult(k[2])

        self.assertEqual(c1.transition_matrix, 1 / 3)
        self.assertEqual(c2.transition_matrix, 1 / 3)
        self.assertEqual(c3.transition_matrix, 1 / 3)

    def test_constant_wrong_parentheses(self):
        adata = _adata.copy()
        k = VelocityKernel(adata).compute_transition_matrix() + (
            ConnectivityKernel(adata).compute_transition_matrix()
            + ConnectivityKernel(adata).compute_transition_matrix()
        )
        k.compute_transition_matrix()
        c1, c2, c3 = _is_bin_mult(k[0]), _is_bin_mult(k[1]), _is_bin_mult(k[2])

        self.assertEqual(c1.transition_matrix, 1 / 3)
        self.assertEqual(c2.transition_matrix, 1 / 3)
        self.assertEqual(c3.transition_matrix, 1 / 3)

    def test_constant_correct_parentheses(self):
        adata = _adata.copy()
        k = 1 * VelocityKernel(adata).compute_transition_matrix() + 1 * (
            ConnectivityKernel(adata).compute_transition_matrix()
            + ConnectivityKernel(adata).compute_transition_matrix()
        )
        k.compute_transition_matrix()
        c1, c2, c3 = (
            _is_bin_mult(k[0]),
            _is_bin_mult(k[1][1][0]),
            _is_bin_mult(k[1][1][1]),
        )

        self.assertEqual(c1.transition_matrix, 1 / 2)
        self.assertEqual(c2.transition_matrix, 1 / 2)
        self.assertEqual(c3.transition_matrix, 1 / 2)

    def test_adaptive_kernel_constants(self):
        adata = _adata.copy()
        k = (3 * ConnectivityKernel(adata).compute_transition_matrix()) ^ (
            1 * ConnectivityKernel(adata).compute_transition_matrix()
        )
        k.compute_transition_matrix()

        self.assertEqual(k[0][0]._value, 3 / 4)
        self.assertEqual(k[1][0]._value, 1 / 4)

    def test_adaptive_kernel_complex(self):
        adata = _adata.copy()
        k = (
            4
            * (
                (3 * ConnectivityKernel(adata).compute_transition_matrix())
                ^ (1 * ConnectivityKernel(adata).compute_transition_matrix())
            )
            + 2 * ConnectivityKernel(adata).compute_transition_matrix()
        )
        k.compute_transition_matrix()

        self.assertEqual(k[0][0].transition_matrix, 4 / 6)
        self.assertEqual(k[1][0].transition_matrix, 2 / 6)
        self.assertEqual(k[0][1][0][0]._value, 3 / 4)
        self.assertEqual(k[0][1][1][0]._value, 1 / 4)


class KernelTestCase(unittest.TestCase):
    def test_row_normalized(self):
        adata = _adata.copy()
        vk = VelocityKernel(adata)

        vk.compute_transition_matrix(density_normalize=False)
        T = vk.transition_matrix
        np.testing.assert_allclose(T.sum(1), 1, rtol=_rtol)

    def test_row_normalized_dense_norm(self):
        adata = _adata.copy()
        vk = VelocityKernel(adata)
        vk.compute_transition_matrix(density_normalize=True)
        T = vk.transition_matrix

        np.testing.assert_allclose(T.sum(1), 1, rtol=_rtol)

    def test_transition_forward(self):
        adata = _adata.copy()
        backward = False

        vk = VelocityKernel(adata, backward=backward).compute_transition_matrix(
            density_normalize=True
        )
        T_1 = vk.transition_matrix

        transition_matrix(adata, density_normalize=True, backward=backward)
        T_2 = adata.uns[_transition(Direction.FORWARD)]["T"]

        np.testing.assert_allclose(T_1.A, T_2.A, rtol=_rtol)

    def test_transition_forward_dense_norm(self):
        adata = _adata.copy()
        backward = False

        vk = VelocityKernel(adata, backward=backward).compute_transition_matrix(
            density_normalize=False
        )
        T_1 = vk.transition_matrix

        transition_matrix(adata, density_normalize=False, backward=backward)
        T_2 = adata.uns[_transition(Direction.FORWARD)]["T"]

        np.testing.assert_allclose(T_1.A, T_2.A, rtol=_rtol)

    def test_transition_forward_differ_dense_norm(self):
        adata = _adata.copy()
        backward = False

        vk = VelocityKernel(adata, backward=backward).compute_transition_matrix(
            density_normalize=True
        )
        T_1 = vk.transition_matrix

        transition_matrix(adata, density_normalize=False, backward=backward)
        T_2 = adata.uns[_transition(Direction.FORWARD)]["T"]

        self.assertFalse(np.allclose(T_1.A, T_2.A, rtol=_rtol))

    def test_transition_backward(self):
        adata = _adata.copy()
        backward = True

        vk = VelocityKernel(adata, backward=backward).compute_transition_matrix(
            density_normalize=True
        )
        T_1 = vk.transition_matrix

        transition_matrix(adata, density_normalize=True, backward=backward)
        T_2 = adata.uns[_transition(Direction.BACKWARD)]["T"]

        np.testing.assert_allclose(T_1.A, T_2.A, rtol=_rtol)

    def test_transition_backward_dense_norm(self):
        adata = _adata.copy()
        backward = True

        vk = VelocityKernel(adata, backward=backward).compute_transition_matrix(
            density_normalize=False
        )
        T_1 = vk.transition_matrix

        transition_matrix(adata, density_normalize=False, backward=backward)
        T_2 = adata.uns[_transition(Direction.BACKWARD)]["T"]

        np.testing.assert_allclose(T_1.A, T_2.A, rtol=_rtol)

    def test_transition_backward_differ_dense_norm(self):
        adata = _adata.copy()
        backward = True

        vk = VelocityKernel(adata, backward=backward).compute_transition_matrix(
            density_normalize=True
        )
        T_1 = vk.transition_matrix

        transition_matrix(adata, density_normalize=False, backward=backward)
        T_2 = adata.uns[_transition(Direction.BACKWARD)]["T"]

        self.assertFalse(np.allclose(T_1.A, T_2.A, rtol=_rtol))

    def test_backward_negate(self):
        adata = _adata.copy()
        backward = True
        dense_norm = False
        backward_mode = "negate"
        vk = VelocityKernel(adata, backward=backward)

        vk.compute_transition_matrix(
            density_normalize=dense_norm, backward_mode=backward_mode
        )
        T_1 = vk.transition_matrix
        transition_matrix(
            adata,
            density_normalize=dense_norm,
            backward=backward,
            backward_mode=backward_mode,
        )
        T_2 = adata.uns[_transition(Direction.BACKWARD)]["T"]

        np.testing.assert_allclose(T_1.A, T_2.A, rtol=_rtol)

    def test_backward_negate_dense_norm(self):
        adata = _adata.copy()
        backward = True
        dense_norm = True
        backward_mode = "negate"
        vk = VelocityKernel(adata, backward=backward)

        vk.compute_transition_matrix(
            density_normalize=dense_norm, backward_mode=backward_mode
        )
        T_1 = vk.transition_matrix
        transition_matrix(
            adata,
            density_normalize=dense_norm,
            backward=backward,
            backward_mode=backward_mode,
        )
        T_2 = adata.uns["T_bwd"]["T"]

        np.testing.assert_allclose(T_1.A, T_2.A, rtol=_rtol)

    def test_backward_negate_differ(self):
        adata = _adata.copy()
        backward = True
        backward_mode = "negate"
        vk = VelocityKernel(adata, backward=backward)

        vk.compute_transition_matrix(
            density_normalize=True, backward_mode=backward_mode
        )
        T_1 = vk.transition_matrix
        transition_matrix(
            adata,
            density_normalize=False,
            backward=backward,
            backward_mode=backward_mode,
        )
        T_2 = adata.uns["T_bwd"]["T"]

        self.assertFalse(np.allclose(T_1.A, T_2.A, rtol=_rtol))

    def test_palantir(self):
        adata = _adata.copy()
        conn = adata.uns["neighbors"]["connectivities"]
        n_neighbors = adata.uns["neighbors"]["params"]["n_neighbors"]
        pseudotime = adata.obs["latent_time"]

        conn_biased = bias_knn(conn, pseudotime, n_neighbors)
        T_1 = _normalize(conn_biased)

        pk = PalantirKernel(adata, time_key="latent_time").compute_transition_matrix(
            density_normalize=False
        )
        T_2 = pk.transition_matrix

        np.testing.assert_allclose(T_1.A, T_2.A, rtol=_rtol)

    def test_palantir_dense_norm(self):
        adata = _adata.copy()
        conn = adata.uns["neighbors"]["connectivities"]
        n_neighbors = adata.uns["neighbors"]["params"]["n_neighbors"]
        pseudotime = adata.obs["latent_time"]

        conn_biased = bias_knn(conn, pseudotime, n_neighbors)
        T_1 = density_normalization(conn_biased, conn)
        T_1 = _normalize(T_1)

        pk = PalantirKernel(adata, time_key="latent_time").compute_transition_matrix(
            density_normalize=True
        )
        T_2 = pk.transition_matrix

        np.testing.assert_allclose(T_1.A, T_2.A, rtol=_rtol)

    def test_palantir_differ_dense_norm(self):
        adata = _adata.copy()
        conn = adata.uns["neighbors"]["connectivities"]
        n_neighbors = adata.uns["neighbors"]["params"]["n_neighbors"]
        pseudotime = adata.obs["latent_time"]

        conn_biased = bias_knn(conn, pseudotime, n_neighbors)
        T_1 = density_normalization(conn_biased, conn)
        T_1 = _normalize(T_1)

        pk = PalantirKernel(adata, time_key="latent_time").compute_transition_matrix(
            density_normalize=False
        )
        T_2 = pk.transition_matrix

        self.assertFalse(np.allclose(T_1.A, T_2.A, rtol=_rtol))

    def test_manual_combination(self):
        adata = _adata.copy()
        density_normalize = False

        vk = VelocityKernel(adata).compute_transition_matrix(
            density_normalize=density_normalize
        )
        ck = ConnectivityKernel(adata).compute_transition_matrix(
            density_normalize=density_normalize
        )

        T_vk = vk.transition_matrix
        T_ck = ck.transition_matrix
        T_comb_manual = 0.8 * T_vk + 0.2 * T_ck

        comb_kernel = 0.8 * vk + 0.2 * ck
        T_comb_kernel = comb_kernel.transition_matrix

        np.testing.assert_allclose(T_comb_kernel.A, T_comb_manual.A, rtol=_rtol)

    def test_manual_combination_dense_norm(self):
        adata = _adata.copy()
        density_normalize = True

        vk = VelocityKernel(adata).compute_transition_matrix(
            density_normalize=density_normalize
        )
        ck = ConnectivityKernel(adata).compute_transition_matrix(
            density_normalize=density_normalize
        )

        T_vk = vk.transition_matrix
        T_ck = ck.transition_matrix
        T_comb_manual = 0.8 * T_vk + 0.2 * T_ck

        comb_kernel = 0.8 * vk + 0.2 * ck
        T_comb_kernel = comb_kernel.transition_matrix

        np.testing.assert_allclose(T_comb_kernel.A, T_comb_manual.A, rtol=_rtol)

    def test_manual_combination_no_precomputed(self):
        adata = _adata.copy()
        density_normalize = False

        vk = VelocityKernel(adata).compute_transition_matrix(
            density_normalize=density_normalize
        )
        ck = ConnectivityKernel(adata).compute_transition_matrix(
            density_normalize=density_normalize
        )

        T_vk = vk.transition_matrix
        T_ck = ck.transition_matrix
        T_comb_manual = 0.8 * T_vk + 0.2 * T_ck

        vk = VelocityKernel(adata).compute_transition_matrix(
            density_normalize=density_normalize
        )
        ck = ConnectivityKernel(adata).compute_transition_matrix(
            density_normalize=density_normalize
        )
        comb_kernel = 0.8 * vk + 0.2 * ck
        comb_kernel.compute_transition_matrix()
        T_comb_kernel = comb_kernel.transition_matrix

        np.testing.assert_allclose(T_comb_manual.A, T_comb_kernel.A, rtol=_rtol)

    def test_manual_combination_backward(self):
        adata = _adata.copy()
        backward, density_normalize = True, False

        vk = VelocityKernel(adata, backward=backward).compute_transition_matrix(
            density_normalize=density_normalize
        )
        ck = ConnectivityKernel(adata, backward=backward).compute_transition_matrix(
            density_normalize=density_normalize
        )

        T_vk = vk.transition_matrix
        T_ck = ck.transition_matrix
        T_comb_manual = 0.8 * T_vk + 0.2 * T_ck

        comb_kernel = 0.8 * vk + 0.2 * ck
        T_comb_kernel = comb_kernel.transition_matrix

        np.testing.assert_allclose(T_comb_manual.A, T_comb_kernel.A, rtol=_rtol)

    def test_manual_combination_backward_dense_norm(self):
        adata = _adata.copy()
        backward, density_normalize = True, True

        vk = VelocityKernel(adata, backward=backward).compute_transition_matrix(
            density_normalize=density_normalize
        )
        ck = ConnectivityKernel(adata, backward=backward).compute_transition_matrix(
            density_normalize=density_normalize
        )

        T_vk = vk.transition_matrix
        T_ck = ck.transition_matrix
        T_comb_manual = 0.8 * T_vk + 0.2 * T_ck

        comb_kernel = 0.8 * vk + 0.2 * ck
        T_comb_kernel = comb_kernel.transition_matrix

        np.testing.assert_allclose(T_comb_manual.A, T_comb_kernel.A, rtol=_rtol)


class PreviousImplTestCase(unittest.TestCase):
    def test_foward(self):
        adata = _adata.copy()
        density_normalize = False
        vk = VelocityKernel(adata).compute_transition_matrix(
            density_normalize=density_normalize
        )
        ck = ConnectivityKernel(adata).compute_transition_matrix(
            density_normalize=density_normalize
        )

        comb = 0.8 * vk + 0.2 * ck
        T_1 = comb.transition_matrix

        transition_matrix(
            adata, diff_kernel="sum", weight_diffusion=0.2, density_normalize=False
        )
        T_2 = adata.uns[_transition(Direction.FORWARD)]["T"]

        np.testing.assert_allclose(T_1.A, T_2.A, rtol=_rtol)

    def test_forward_manual_dense_norm(self):
        adata = _adata.copy()
        vk = VelocityKernel(adata).compute_transition_matrix(density_normalize=False)
        ck = ConnectivityKernel(adata).compute_transition_matrix(
            density_normalize=False
        )

        comb = 0.8 * vk + 0.2 * ck
        T_1 = comb.transition_matrix
        conn = adata.uns["neighbors"]["connectivities"]
        T_1 = density_normalization(T_1, conn)
        T_1 = _normalize(T_1)

        transition_matrix(
            adata, diff_kernel="sum", weight_diffusion=0.2, density_normalize=True
        )
        T_2 = adata.uns[_transition(Direction.FORWARD)]["T"]

        np.testing.assert_allclose(T_1.A, T_2.A, rtol=_rtol)

    def test_backward(self):
        adata = _adata.copy()
        density_norm, backward = False, True
        vk = VelocityKernel(adata, backward=backward).compute_transition_matrix(
            density_normalize=density_norm
        )
        ck = ConnectivityKernel(adata, backward=backward).compute_transition_matrix(
            density_normalize=density_norm
        )

        # combine the kernels
        comb = 0.8 * vk + 0.2 * ck
        T_1 = comb.transition_matrix

        transition_matrix(
            adata,
            diff_kernel="sum",
            weight_diffusion=0.2,
            density_normalize=False,
            backward=backward,
        )
        T_2 = adata.uns[_transition(Direction.BACKWARD)]["T"]

        np.testing.assert_allclose(T_1.A, T_2.A, rtol=_rtol)

    def test_backward_manual_dense_norm(self):
        adata = _adata.copy()
        backward = True
        vk = VelocityKernel(adata, backward=backward).compute_transition_matrix(
            density_normalize=False
        )
        ck = ConnectivityKernel(adata, backward=backward).compute_transition_matrix(
            density_normalize=False
        )

        # combine the kernels
        comb = 0.8 * vk + 0.2 * ck
        T_1 = comb.transition_matrix
        conn = adata.uns["neighbors"]["connectivities"]
        T_1 = density_normalization(T_1, conn)
        T_1 = _normalize(T_1)

        transition_matrix(
            adata,
            diff_kernel="sum",
            weight_diffusion=0.2,
            density_normalize=True,
            backward=backward,
        )
        T_2 = adata.uns[_transition(Direction.BACKWARD)]["T"]

        np.testing.assert_allclose(T_1.A, T_2.A, rtol=_rtol)


class AdditionTestCase(unittest.TestCase):
    def test_simple_addition(self):
        vk, ck = create_kernels(_adata)  # diagonal + upper diag

        k = (vk + ck).compute_transition_matrix()
        expected = np.eye(_adata.n_obs) * 0.75 + np.eye(_adata.n_obs, k=1) * 0.25
        expected[-1, -1] = 1

        np.testing.assert_allclose(k.transition_matrix.A, expected)

    def test_addtion_with_constant(self):
        vk, ck = create_kernels(_adata)  # diagonal + upper diag

        k = (2 * vk + 3 * ck).compute_transition_matrix()
        expected = (
            np.eye(_adata.n_obs) * (2 / 5)
            + np.eye(_adata.n_obs) * (3 / 5) * 0.5
            + np.eye(_adata.n_obs, k=1) * (3 / 5) * 0.5
        )
        expected[-1, -1] = 1

        np.testing.assert_allclose(k.transition_matrix.A, expected)

    def test_addition_3_kernels(self):
        adata = _adata.copy()
        vk, ck = create_kernels(adata)  # diagonal + upper diag
        vk1 = VelocityKernel(adata)
        vk1._transition_matrix = np.eye(adata.n_obs, k=-1) / 2 + np.eye(adata.n_obs) / 2
        vk1._transition_matrix[0, 0] = 1
        np.testing.assert_allclose(
            np.sum(ck._transition_matrix, axis=1), 1
        )  # sanity check

        k = (vk + ck + vk1).compute_transition_matrix()
        expected = (
            np.eye(adata.n_obs) * (1 / 3 + 1 / 6 + 1 / 6)
            + np.eye(adata._n_obs, k=1) * 1 / 6
            + np.eye(adata.n_obs, k=-1) * 1 / 6
        )
        expected[0, 0] = expected[-1, -1] = 2 / 3 + 1 / 3 * 0.5
        expected[0, 1] = expected[-1, -2] = 1 - expected[0, 0]

        np.testing.assert_allclose(k.transition_matrix.A, expected)

    def test_addition_adaptive(self):
        adata = _adata.copy()
        adata.uns["velocity_variances"] = vv = np.random.random(
            size=(adata.n_obs, adata.n_obs)
        )
        adata.uns["connectivity_variances"] = cv = np.random.random(
            size=(adata.n_obs, adata.n_obs)
        )
        vk, ck = create_kernels(adata)

        k = vk ^ ck
        expected = _normalize(
            0.5 * vv * vk.transition_matrix + 0.5 * cv * ck.transition_matrix
        )

        np.testing.assert_allclose(k.transition_matrix.A, expected)

    def test_addition_adataptive_constants(self):
        adata = _adata.copy()
        a, b = np.random.uniform(0, 10, 2)
        s = a + b
        adata.uns["velocity_variances"] = vv = np.random.random(
            size=(adata.n_obs, adata.n_obs)
        )
        adata.uns["connectivity_variances"] = cv = np.random.random(
            size=(adata.n_obs, adata.n_obs)
        )
        vk, ck = create_kernels(adata)  # diagonal + upper diag

        k = a * vk ^ b * ck
        expected = _normalize(
            a / s * vv * vk.transition_matrix + b / s * cv * ck.transition_matrix
        )

        np.testing.assert_allclose(k.transition_matrix.A, expected)

    def test_addition_adaptive_wrong_variances(self):
        adata = _adata.copy()
        a, b = np.random.uniform(0, 10, 2)
        s = a + b
        adata.uns["velocity_variances"] = vv = np.random.random(
            size=(adata.n_obs, adata.n_obs)
        )
        adata.uns["connectivity_variances"] = cv = np.random.random(
            size=(adata.n_obs, adata.n_obs)
        )
        vk, ck = create_kernels(adata)

        k = a * vk ^ b * ck
        expected = _normalize(
            a / s * vk.transition_matrix + b / s * ck.transition_matrix
        )

        self.assertFalse(np.allclose(k.transition_matrix.A, expected))

    def test_addition_adaptive_4_kernels(self):
        adata = _adata.copy()
        a, b, c, d = np.random.uniform(0, 10, 4)
        s = a + b + c + d
        adata.uns["velocity_variances"] = vv = np.random.random(
            size=(adata.n_obs, adata.n_obs)
        )
        adata.uns["connectivity_variances"] = cv = np.random.random(
            size=(adata.n_obs, adata.n_obs)
        )
        vk, ck = create_kernels(adata)
        vk1, ck1 = create_kernels(adata)

        k = a * vk ^ b * ck ^ c * vk1 ^ d * ck1
        expected = _normalize(
            a / s * vv * vk.transition_matrix
            + b / s * cv * ck.transition_matrix
            + c / s * vv * vk1.transition_matrix
            + d / s * cv * ck1.transition_matrix
        )

        np.testing.assert_allclose(k.transition_matrix.A, expected)


if __name__ == "__main__":
    unittest.main()
