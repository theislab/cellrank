# -*- coding: utf-8 -*-
import numpy as np
import pytest

import anndata
from scanpy import Neighbors

import cellrank as cr
from _helpers import (
    bias_knn,
    create_kernels,
    transition_matrix,
    density_normalization,
    jax_not_installed_skip,
    random_transition_matrix,
)
from cellrank.tools._utils import _normalize
from cellrank.utils._utils import _get_neighs, _get_neighs_params
from cellrank.tools.kernels import (
    Constant,
    PalantirKernel,
    VelocityKernel,
    PrecomputedKernel,
    ConnectivityKernel,
)
from cellrank.tools._constants import Direction, _transition
from cellrank.tools.kernels._base_kernel import KernelAdd, KernelMul, _is_bin_mult

_rtol = 1e-6


class TestInitializeKernel:
    def test_none_transition_matrix(self, adata):
        vk = VelocityKernel(adata)
        ck = ConnectivityKernel(adata)
        pk = PalantirKernel(adata, time_key="latent_time")

        assert vk._transition_matrix is None
        assert ck._transition_matrix is None
        assert pk._transition_matrix is None

    def test_not_none_transition_matrix_compute(self, adata):
        vk = VelocityKernel(adata).compute_transition_matrix()
        ck = ConnectivityKernel(adata).compute_transition_matrix()
        pk = PalantirKernel(adata, time_key="latent_time").compute_transition_matrix()

        assert vk.transition_matrix is not None
        assert ck.transition_matrix is not None
        assert pk.transition_matrix is not None

    def test_not_none_transition_matrix_accessor(self, adata):
        vk = VelocityKernel(adata)
        ck = ConnectivityKernel(adata)
        pk = PalantirKernel(adata, time_key="latent_time")

        assert vk.transition_matrix is not None
        assert ck.transition_matrix is not None
        assert pk.transition_matrix is not None

    def test_adding_hidden_constants(self, adata):
        k: KernelAdd = VelocityKernel(adata) + ConnectivityKernel(adata)

        assert _is_bin_mult(k[0])
        assert isinstance(k[0], KernelMul)
        assert isinstance(k[0][0], Constant)
        assert isinstance(k[0][1], VelocityKernel)
        assert k[0][0].transition_matrix == 1.0

        assert _is_bin_mult(k[1])
        assert isinstance(k[1], KernelMul)
        assert isinstance(k[1][0], Constant)
        assert isinstance(k[1][1], ConnectivityKernel)
        assert k[1][0].transition_matrix == 1.0

    def test_length(self, adata):
        k: KernelAdd = VelocityKernel(adata) + ConnectivityKernel(adata)
        assert len(k) == 2

    def test_accessor_out_of_range(self, adata):
        k: KernelAdd = VelocityKernel(adata) + ConnectivityKernel(adata)

        with pytest.raises(IndexError):
            _ = k[2]

    def test_parent(self, adata):
        vk = VelocityKernel(adata)
        ck = ConnectivityKernel(adata)
        k = vk + ck

        assert vk._parent._parent is k  # invisible constants
        assert ck._parent._parent is k
        assert k._parent is None

    def test_uninitialized_both(self, adata):
        k = VelocityKernel(adata) + ConnectivityKernel(adata)

        with pytest.raises(RuntimeError):
            k.compute_transition_matrix()

    def test_uninitialized_one(self, adata):
        k = (
            VelocityKernel(adata)
            + ConnectivityKernel(adata).compute_transition_matrix()
        )

        with pytest.raises(RuntimeError):
            k.compute_transition_matrix()

    def test_initialized(self, adata):
        k = (
            VelocityKernel(adata).compute_transition_matrix()
            + ConnectivityKernel(adata).compute_transition_matrix()
        )
        k.compute_transition_matrix()

        assert k.transition_matrix is not None

    def test_invalida_type(self, adata):
        with pytest.raises(TypeError):
            _ = None * VelocityKernel(adata)

    def test_negative_constant(self, adata):
        with pytest.raises(ValueError):
            _ = -1 * VelocityKernel(adata)

    def test_invalid_constant(self, adata):
        with pytest.raises(TypeError):
            _ = Constant(None, None)

    def test_inversion(self, adata):
        c = ConnectivityKernel(adata, backward=False)
        assert not c.backward

        nc = ~c
        assert nc.backward

    def test_inversion_inplace(self, adata):
        c = ConnectivityKernel(adata, backward=False)

        assert not c.backward
        _ = ~c
        assert c.backward

    def test_inversion_propagation(self, adata):
        c = ConnectivityKernel(adata, backward=False)
        v = VelocityKernel(adata, backward=False)
        k = ~(c + v)

        assert c.backward
        assert v.backward
        assert k.backward

    def test_inversion_recalculation(self, adata):
        c = ConnectivityKernel(adata).compute_transition_matrix()
        z = ~(c + c)
        with pytest.raises(RuntimeError):
            z.compute_transition_matrix()

    def test_inversion_preservation_of_constants(self, adata):
        c = ConnectivityKernel(adata).compute_transition_matrix()
        a = (3 * c + 1 * c).compute_transition_matrix()
        b = ~a
        c.compute_transition_matrix()

        assert a[0][0].transition_matrix == 3 / 4
        assert b[0][0].transition_matrix == 3 / 4
        assert a[1][0].transition_matrix == 1 / 4
        assert b[1][0].transition_matrix == 1 / 4

    def test_addition_simple(self, adata):
        k = VelocityKernel(adata) + ConnectivityKernel(adata)

        assert isinstance(k, KernelAdd)

    def test_multiplication_simple(self, adata):
        k = 10 * VelocityKernel(adata)
        c = _is_bin_mult(k)

        assert isinstance(c, Constant)
        assert c.transition_matrix == 10

    def test_multiplication_simple_normalization(self, adata):
        k = 10 * VelocityKernel(adata).compute_transition_matrix()
        c = _is_bin_mult(k)

        assert c.transition_matrix == 10

    def test_constant(self, adata):
        k = 9 * VelocityKernel(adata) + 1 * ConnectivityKernel(adata)
        c1, c2 = _is_bin_mult(k[0]), _is_bin_mult(k[1])

        assert c1.transition_matrix == 9
        assert c2.transition_matrix == 1

    def test_constant_normalize_2(self, adata):
        k = (
            9 * VelocityKernel(adata).compute_transition_matrix()
            + 1 * ConnectivityKernel(adata).compute_transition_matrix()
        )
        k.compute_transition_matrix()
        c1, c2 = _is_bin_mult(k[0]), _is_bin_mult(k[1])

        assert c1.transition_matrix == 9 / 10
        assert c2.transition_matrix == 1 / 10

    def test_constant_normalize_3(self, adata):
        k = (
            VelocityKernel(adata).compute_transition_matrix()
            + ConnectivityKernel(adata).compute_transition_matrix()
            + ConnectivityKernel(adata).compute_transition_matrix()
        )
        k.compute_transition_matrix()
        c1, c2, c3 = _is_bin_mult(k[0]), _is_bin_mult(k[1]), _is_bin_mult(k[2])

        assert c1.transition_matrix == 1 / 3
        assert c2.transition_matrix == 1 / 3
        assert c3.transition_matrix == 1 / 3

    def test_constant_wrong_parentheses(self, adata):
        k = VelocityKernel(adata).compute_transition_matrix() + (
            ConnectivityKernel(adata).compute_transition_matrix()
            + ConnectivityKernel(adata).compute_transition_matrix()
        )
        k.compute_transition_matrix()
        c1, c2, c3 = _is_bin_mult(k[0]), _is_bin_mult(k[1]), _is_bin_mult(k[2])

        assert c1.transition_matrix == 1 / 3
        assert c2.transition_matrix == 1 / 3
        assert c3.transition_matrix == 1 / 3

    def test_constant_correct_parentheses(self, adata):
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

        assert c1.transition_matrix == 1 / 2
        assert c2.transition_matrix == 1 / 2
        assert c3.transition_matrix == 1 / 2

    def test_adaptive_kernel_constants(self, adata):
        ck1 = ConnectivityKernel(adata).compute_transition_matrix()
        ck1._mat_scaler = np.random.normal(size=(adata.n_obs, adata.n_obs))

        ck2 = ConnectivityKernel(adata).compute_transition_matrix()
        ck2._mat_scaler = np.random.normal(size=(adata.n_obs, adata.n_obs))

        k = (3 * ck1) ^ (1 * ck2)
        k.compute_transition_matrix()

        assert k[0][0]._value == 3 / 4
        assert k[1][0]._value == 1 / 4

    def test_adaptive_kernel_complex(self, adata):
        ck1 = ConnectivityKernel(adata).compute_transition_matrix()
        ck1._mat_scaler = np.random.normal(size=(adata.n_obs, adata.n_obs))

        ck2 = ConnectivityKernel(adata).compute_transition_matrix()
        ck2._mat_scaler = np.random.normal(size=(adata.n_obs, adata.n_obs))

        ck3 = ConnectivityKernel(adata).compute_transition_matrix()
        ck3._mat_scaler = np.random.normal(size=(adata.n_obs, adata.n_obs))

        k = 4 * ((3 * ck1) ^ (1 * ck2)) + 2 * ck3
        k.compute_transition_matrix()

        assert k[0][0].transition_matrix == 4 / 6
        assert k[1][0].transition_matrix == 2 / 6
        assert k[0][1][0][0]._value == 3 / 4
        assert k[0][1][1][0]._value == 1 / 4

    def test_repr(self, adata):
        rpr = repr(VelocityKernel(adata))

        assert rpr == "<Velo>"

    def test_repr_inv(self, adata):
        rpr = repr(~VelocityKernel(adata))

        assert rpr == "~<Velo>"

    def test_repr_inv_comb(self, adata):
        rpr = repr(~(VelocityKernel(adata) + ConnectivityKernel(adata)))

        assert rpr == "~((1 * <Velo>) + (1 * <Conn>))"

    def test_str_repr_equiv_no_transition_matrix(self, adata):
        vk = VelocityKernel(adata)
        string = str(vk)
        rpr = repr(vk)

        assert string == rpr
        assert string == "<Velo>"

    def test_str(self, adata):
        string = str(ConnectivityKernel(adata).compute_transition_matrix())

        assert string == "<Conn[dnorm=True]>"

    def test_str_inv(self, adata):
        string = str(
            ConnectivityKernel(adata, backward=True).compute_transition_matrix()
        )

        assert string == "~<Conn[dnorm=True]>"


class TestKernel:
    def test_precomputed_not_array(self):
        with pytest.raises(TypeError):
            _ = PrecomputedKernel([[1, 0], [0, 1]])

    def test_precomputed_not_square(self):
        with pytest.raises(ValueError):
            _ = PrecomputedKernel(np.random.normal(size=(10, 9)))

    def test_precomputed_not_a_transition_matrix(self):
        mat = random_transition_matrix(100)
        mat[0, 0] = 0xDEADBEEF
        with pytest.raises(ValueError):
            _ = PrecomputedKernel(mat)

    def test_precomputed_no_adata(self):
        pk = PrecomputedKernel(random_transition_matrix(50))
        pk.write_to_adata()

        assert isinstance(pk.adata, anndata.AnnData)
        assert pk.adata.shape == (50, 1)
        assert pk.adata.obs.shape == (50, 0)
        assert pk.adata.var.shape == (1, 0)
        assert "T_fwd_params" in pk.adata.uns.keys()
        np.testing.assert_array_equal(
            pk.adata.obsp["T_fwd"].toarray(), pk.transition_matrix.toarray()
        )

    def test_precomputed_adata(self, adata):
        pk = PrecomputedKernel(random_transition_matrix(adata.n_obs), adata=adata)

        assert pk.adata is adata

    def test_precomputed_transition_matrix(self, adata):
        mat = random_transition_matrix(adata.n_obs)
        pk = PrecomputedKernel(mat)

        np.testing.assert_array_equal(mat, pk.transition_matrix.toarray())

    def test_precomputed_sum(self, adata):
        mat = random_transition_matrix(adata.n_obs)
        pk = PrecomputedKernel(mat)
        vk = VelocityKernel(adata).compute_transition_matrix()

        expected = (0.5 * vk.transition_matrix) + (0.5 * pk.transition_matrix)
        actual = (pk + vk).compute_transition_matrix()

        np.testing.assert_array_almost_equal(
            expected.toarray(), actual.transition_matrix.toarray()
        )

    def test_write_adata(self, adata):
        vk = VelocityKernel(adata).compute_transition_matrix()
        vk.write_to_adata()

        assert adata is vk.adata
        assert "T_fwd_params" in adata.uns.keys()
        np.testing.assert_array_equal(
            adata.obsp["T_fwd"].toarray(), vk.transition_matrix.toarray()
        )

    def test_row_normalized(self, adata):
        vk = VelocityKernel(adata)

        vk.compute_transition_matrix(density_normalize=False)
        T = vk.transition_matrix
        np.testing.assert_allclose(T.sum(1), 1, rtol=_rtol)

    def test_row_normalized_dense_norm(self, adata):
        vk = VelocityKernel(adata)
        vk.compute_transition_matrix(mode="deterministic")
        T_det = vk.transition_matrix
        vk.compute_transition_matrix(mode="sampling")
        T_sam = vk.transition_matrix

        np.testing.assert_allclose(T_det.sum(1), 1, rtol=_rtol)
        np.testing.assert_allclose(T_sam.sum(1), 1, rtol=_rtol)

    @jax_not_installed_skip
    def test_row_normalized_dense_norm_stoch(self, adata):
        vk = VelocityKernel(adata)
        vk.compute_transition_matrix(mode="stochastic")
        T_stoch = vk.transition_matrix

        np.testing.assert_allclose(T_stoch.sum(1), 1, rtol=_rtol)

    @jax_not_installed_skip
    def test_transition_forward_stoch(self, adata):
        # TODO: self-referential test (i.e. useless)
        backward = False

        vk = VelocityKernel(adata, backward=backward).compute_transition_matrix(
            mode="stochastic"
        )
        T_1 = vk.transition_matrix

        T_2 = cr.tl.transition_matrix(
            adata, mode="stochastic", backward=backward
        ).transition_matrix

        np.testing.assert_allclose(T_1.A, T_2.A, rtol=_rtol)

    @pytest.mark.xfail(reason="implementation difference")
    def test_transition_forward_det(self, adata):
        backward = False

        vk = VelocityKernel(adata, backward=backward).compute_transition_matrix(
            mode="deterministic"
        )
        T_1 = vk.transition_matrix

        transition_matrix(adata, backward=backward)
        T_2 = adata.uns[_transition(Direction.FORWARD)]["T"]

        np.testing.assert_allclose(T_1.A, T_2.A, rtol=_rtol)

    @pytest.mark.xfail(reason="implementation difference")
    def test_transition_backward_det(self, adata):
        backward = True

        vk = VelocityKernel(adata, backward=backward).compute_transition_matrix(
            mode="deterministic"
        )
        T_1 = vk.transition_matrix

        transition_matrix(adata, backward=backward)
        T_2 = adata.uns[_transition(Direction.BACKWARD)]["T"]

        np.testing.assert_allclose(T_1.A, T_2.A, rtol=_rtol)

    @jax_not_installed_skip
    def test_transition_forward_differ_mode(self, adata):
        # TODO: doesn't make sense - different implementation
        backward = False

        vk = VelocityKernel(adata, backward=backward).compute_transition_matrix(
            mode="stochastic"
        )
        T_1 = vk.transition_matrix

        # this is the very old deterministic approach
        transition_matrix(adata, backward=backward)
        T_2 = adata.uns[_transition(Direction.FORWARD)]["T"]

        assert not np.allclose(T_1.A, T_2.A, rtol=_rtol)

    @pytest.mark.xfail(reason="implementation difference")
    def test_backward_negate(self, adata):
        backward = True
        backward_mode = "negate"
        vk = VelocityKernel(adata, backward=backward)

        vk.compute_transition_matrix(backward_mode=backward_mode, mode="deterministic")
        T_1 = vk.transition_matrix
        transition_matrix(
            adata, backward=backward, backward_mode=backward_mode,
        )
        T_2 = adata.uns[_transition(Direction.BACKWARD)]["T"]

        np.testing.assert_allclose(T_1.A, T_2.A, rtol=_rtol)

    def test_palantir(self, adata):
        conn = _get_neighs(adata, "connectivities")
        n_neighbors = _get_neighs_params(adata)["n_neighbors"]
        pseudotime = adata.obs["latent_time"]

        conn_biased = bias_knn(conn, pseudotime, n_neighbors)
        T_1 = _normalize(conn_biased)

        pk = PalantirKernel(adata, time_key="latent_time").compute_transition_matrix(
            density_normalize=False
        )
        T_2 = pk.transition_matrix

        np.testing.assert_allclose(T_1.A, T_2.A, rtol=_rtol)

    def test_palantir_dense_norm(self, adata):
        conn = _get_neighs(adata, "connectivities")
        n_neighbors = _get_neighs_params(adata)["n_neighbors"]
        pseudotime = adata.obs["latent_time"]

        conn_biased = bias_knn(conn, pseudotime, n_neighbors)
        T_1 = density_normalization(conn_biased, conn)
        T_1 = _normalize(T_1)

        pk = PalantirKernel(adata, time_key="latent_time").compute_transition_matrix(
            density_normalize=True
        )
        T_2 = pk.transition_matrix

        np.testing.assert_allclose(T_1.A, T_2.A, rtol=_rtol)

    def test_palantir_differ_dense_norm(self, adata):
        conn = _get_neighs(adata, "connectivities")
        n_neighbors = _get_neighs_params(adata)["n_neighbors"]
        pseudotime = adata.obs["latent_time"]

        conn_biased = bias_knn(conn, pseudotime, n_neighbors)
        T_1 = density_normalization(conn_biased, conn)
        T_1 = _normalize(T_1)

        pk = PalantirKernel(adata, time_key="latent_time").compute_transition_matrix(
            density_normalize=False
        )
        T_2 = pk.transition_matrix

        assert not np.allclose(T_1.A, T_2.A, rtol=_rtol)

    def test_manual_combination_det(self, adata):
        mode = "deterministic"

        vk = VelocityKernel(adata).compute_transition_matrix(mode=mode)
        ck = ConnectivityKernel(adata).compute_transition_matrix()

        T_vk = vk.transition_matrix
        T_ck = ck.transition_matrix
        T_comb_manual = 0.8 * T_vk + 0.2 * T_ck

        comb_kernel = 0.8 * vk + 0.2 * ck
        T_comb_kernel = comb_kernel.transition_matrix

        np.testing.assert_allclose(T_comb_kernel.A, T_comb_manual.A, rtol=_rtol)

    @jax_not_installed_skip
    def test_manual_combination_stoch(self, adata):
        mode = "stochastic"

        vk = VelocityKernel(adata).compute_transition_matrix(mode=mode)
        ck = ConnectivityKernel(adata).compute_transition_matrix()

        T_vk = vk.transition_matrix
        T_ck = ck.transition_matrix
        T_comb_manual = 0.8 * T_vk + 0.2 * T_ck

        comb_kernel = 0.8 * vk + 0.2 * ck
        T_comb_kernel = comb_kernel.transition_matrix

        np.testing.assert_allclose(T_comb_kernel.A, T_comb_manual.A, rtol=_rtol)

    def test_manual_combination_sam(self, adata):
        mode = "sampling"

        vk = VelocityKernel(adata).compute_transition_matrix(mode=mode)
        ck = ConnectivityKernel(adata).compute_transition_matrix()

        T_vk = vk.transition_matrix
        T_ck = ck.transition_matrix
        T_comb_manual = 0.8 * T_vk + 0.2 * T_ck

        comb_kernel = 0.8 * vk + 0.2 * ck
        T_comb_kernel = comb_kernel.transition_matrix

        np.testing.assert_allclose(T_comb_kernel.A, T_comb_manual.A, rtol=_rtol)

    def test_manual_combination_no_precomputed(self, adata):
        density_normalize = False

        vk = VelocityKernel(adata).compute_transition_matrix()
        ck = ConnectivityKernel(adata).compute_transition_matrix(
            density_normalize=density_normalize
        )

        T_vk = vk.transition_matrix
        T_ck = ck.transition_matrix
        T_comb_manual = 0.8 * T_vk + 0.2 * T_ck

        vk = VelocityKernel(adata).compute_transition_matrix()
        ck = ConnectivityKernel(adata).compute_transition_matrix(
            density_normalize=density_normalize
        )
        comb_kernel = 0.8 * vk + 0.2 * ck
        comb_kernel.compute_transition_matrix()
        T_comb_kernel = comb_kernel.transition_matrix

        np.testing.assert_allclose(T_comb_manual.A, T_comb_kernel.A, rtol=_rtol)

    def test_manual_combination_backward(self, adata):
        backward, density_normalize = True, False

        vk = VelocityKernel(adata, backward=backward).compute_transition_matrix()

        ck = ConnectivityKernel(adata, backward=backward).compute_transition_matrix(
            density_normalize=density_normalize
        )

        T_vk = vk.transition_matrix
        T_ck = ck.transition_matrix
        T_comb_manual = 0.8 * T_vk + 0.2 * T_ck

        comb_kernel = 0.8 * vk + 0.2 * ck
        T_comb_kernel = comb_kernel.transition_matrix

        np.testing.assert_allclose(T_comb_manual.A, T_comb_kernel.A, rtol=_rtol)

    def test_manual_combination_backward_dense_norm(self, adata):
        backward, density_normalize = True, True

        vk = VelocityKernel(adata, backward=backward).compute_transition_matrix()

        ck = ConnectivityKernel(adata, backward=backward).compute_transition_matrix(
            density_normalize=density_normalize
        )

        T_vk = vk.transition_matrix
        T_ck = ck.transition_matrix
        T_comb_manual = 0.8 * T_vk + 0.2 * T_ck

        comb_kernel = 0.8 * vk + 0.2 * ck
        T_comb_kernel = comb_kernel.transition_matrix

        np.testing.assert_allclose(T_comb_manual.A, T_comb_kernel.A, rtol=_rtol)

    def compare_with_scanpy_density_normalize(self, adata):
        # check whether cellrank's transition matrix matches scanpy's
        density_normalize = True
        ck = ConnectivityKernel(adata).compute_transition_matrix(
            density_normalize=density_normalize
        )
        T_cr = ck.transition_matrix

        neigh = Neighbors(adata)
        neigh.compute_transitions(density_normalize=density_normalize)
        T_sc = neigh.transitions

        # check whether these are the same while leaving them sparse
        assert T_sc.shape == T_cr.shape
        assert len(T_sc.indices) == len(T_cr.indices)
        assert np.allclose((T_cr - T_sc).data, 0)

    def compare_with_scanpy(self, adata):
        # check whether cellrank's transition matrix matches scanpy's
        density_normalize = False
        ck = ConnectivityKernel(adata).compute_transition_matrix(
            density_normalize=density_normalize
        )
        T_cr = ck.transition_matrix

        neigh = Neighbors(adata)
        neigh.compute_transitions(density_normalize=density_normalize)
        T_sc = neigh.transitions

        # check whether these are the same while leaving them sparse
        assert T_sc.shape == T_cr.shape
        assert len(T_sc.indices) == len(T_cr.indices)
        assert np.allclose((T_cr - T_sc).data, 0)


class TestPreviousImplementation:
    @pytest.mark.xfail(reason="implementation difference")
    def test_forward(self, adata):
        density_normalize = False
        vk = VelocityKernel(adata).compute_transition_matrix(mode="deterministic")
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

    @pytest.mark.xfail(reason="velocity kernel does not have dense norm")
    def test_forward_manual_dense_norm(self, adata):
        vk = VelocityKernel(adata).compute_transition_matrix(mode="deterministic")
        ck = ConnectivityKernel(adata).compute_transition_matrix(
            density_normalize=False
        )

        comb = 0.8 * vk + 0.2 * ck
        T_1 = comb.transition_matrix
        conn = _get_neighs(adata, "connectivities")
        T_1 = density_normalization(T_1, conn)
        T_1 = _normalize(T_1)

        transition_matrix(
            adata, diff_kernel="sum", weight_diffusion=0.2, density_normalize=True
        )
        T_2 = adata.uns[_transition(Direction.FORWARD)]["T"]

        np.testing.assert_allclose(T_1.A, T_2.A, rtol=_rtol)

    @pytest.mark.xfail(reason="implementation difference")
    def test_backward(self, adata):
        density_norm, backward = False, True
        vk = VelocityKernel(adata, backward=backward).compute_transition_matrix()

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

    @pytest.mark.xfail(reason="velocity kernel does not have dense norm")
    def test_backward_manual_dense_norm(self, adata):
        backward = True
        vk = VelocityKernel(adata, backward=backward).compute_transition_matrix()

        ck = ConnectivityKernel(adata, backward=backward).compute_transition_matrix(
            density_normalize=False
        )

        # combine the kernels
        comb = 0.8 * vk + 0.2 * ck
        T_1 = comb.transition_matrix
        conn = _get_neighs(adata, "connectivities")
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


class TestKernelAddition:
    def test_simple_addition(self, adata):
        vk, ck = create_kernels(adata)  # diagonal + upper diag

        k = (vk + ck).compute_transition_matrix()
        expected = np.eye(adata.n_obs) * 0.75 + np.eye(adata.n_obs, k=1) * 0.25
        expected[-1, -1] = 1

        np.testing.assert_allclose(k.transition_matrix.A, expected)

    def test_addtion_with_constant(self, adata):
        vk, ck = create_kernels(adata)  # diagonal + upper diag

        k = (2 * vk + 3 * ck).compute_transition_matrix()
        expected = (
            np.eye(adata.n_obs) * (2 / 5)
            + np.eye(adata.n_obs) * (3 / 5) * 0.5
            + np.eye(adata.n_obs, k=1) * (3 / 5) * 0.5
        )
        expected[-1, -1] = 1

        np.testing.assert_allclose(k.transition_matrix.A, expected)

    def test_addition_3_kernels(self, adata):
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

    def test_addition_adaptive(self, adata):
        adata.uns["velocity_variances"] = vv = np.random.random(
            size=(adata.n_obs, adata.n_obs)
        )
        adata.uns["connectivity_variances"] = cv = np.random.random(
            size=(adata.n_obs, adata.n_obs)
        )
        vk, ck = create_kernels(
            adata,
            velocity_variances="velocity_variances",
            connectivity_variances="connectivity_variances",
        )

        k = vk ^ ck
        expected = _normalize(
            0.5 * vv * vk.transition_matrix + 0.5 * cv * ck.transition_matrix
        )

        np.testing.assert_allclose(k.transition_matrix.A, expected)

    def test_addition_adaptive_constants(self, adata):
        a, b = np.random.uniform(0, 10, 2)
        s = a + b
        adata.uns["velocity_variances"] = vv = np.random.random(
            size=(adata.n_obs, adata.n_obs)
        )
        adata.uns["connectivity_variances"] = cv = np.random.random(
            size=(adata.n_obs, adata.n_obs)
        )
        vk, ck = create_kernels(
            adata,
            velocity_variances="velocity_variances",
            connectivity_variances="connectivity_variances",
        )

        k = a * vk ^ b * ck
        expected = _normalize(
            a / s * vv * vk.transition_matrix + b / s * cv * ck.transition_matrix
        )

        np.testing.assert_allclose(k.transition_matrix.A, expected)

    def test_addition_adaptive_wrong_variances(self, adata):
        a, b = np.random.uniform(0, 10, 2)
        s = a + b
        adata.uns["velocity_variances"] = np.random.random(
            size=(adata.n_obs, adata.n_obs)
        )
        adata.uns["connectivity_variances"] = np.random.random(
            size=(adata.n_obs, adata.n_obs)
        )
        vk, ck = create_kernels(
            adata,
            velocity_variances="velocity_variances",
            connectivity_variances="connectivity_variances",
        )

        k = a * vk ^ b * ck
        expected = _normalize(
            a / s * vk.transition_matrix + b / s * ck.transition_matrix
        )

        assert not np.allclose(k.transition_matrix.A, expected.A)

    def test_addition_adaptive_4_kernels(self, adata):
        a, b, c, d = np.random.uniform(0, 10, 4)
        s = a + b + c + d
        adata.uns["velocity_variances"] = vv = np.random.random(
            size=(adata.n_obs, adata.n_obs)
        )
        adata.uns["connectivity_variances"] = cv = np.random.random(
            size=(adata.n_obs, adata.n_obs)
        )
        vk, ck = create_kernels(
            adata,
            velocity_variances="velocity_variances",
            connectivity_variances="connectivity_variances",
        )
        vk1, ck1 = create_kernels(
            adata,
            velocity_variances="velocity_variances",
            connectivity_variances="connectivity_variances",
        )

        k = a * vk ^ b * ck ^ c * vk1 ^ d * ck1
        expected = _normalize(
            a / s * vv * vk.transition_matrix
            + b / s * cv * ck.transition_matrix
            + c / s * vv * vk1.transition_matrix
            + d / s * cv * ck1.transition_matrix
        )

        np.testing.assert_allclose(k.transition_matrix.A, expected)


class TestKernelCopy:
    def test_copy_simple(self, adata):
        vk1 = VelocityKernel(adata)
        vk2 = vk1.copy()

        assert vk1 is not vk2

    def test_copy_no_adata_copy(self, adata):
        vk1 = VelocityKernel(adata)
        vk2 = vk1.copy()

        assert vk1.adata is adata
        assert vk2.adata is adata

    def test_copy_transition_matrix(self, adata):
        vk1 = VelocityKernel(adata).compute_transition_matrix()
        vk2 = vk1.copy()

        np.testing.assert_array_equal(vk1.transition_matrix.A, vk2.transition_matrix.A)

    def test_copy_params(self, adata):
        vk1 = VelocityKernel(adata).compute_transition_matrix()
        vk2 = vk1.copy()

        assert vk1.params == vk2.params

    def test_copy_cond_num(self, adata):
        for KernelClass in [
            VelocityKernel,
            ConnectivityKernel,
            PalantirKernel,
            PrecomputedKernel,
        ]:
            if KernelClass is PrecomputedKernel:
                k1 = KernelClass(
                    random_transition_matrix(adata.n_obs), compute_cond_num=True
                )
            else:
                k1 = KernelClass(
                    adata, compute_cond_num=True
                ).compute_transition_matrix()
            k2 = k1.copy()

            assert k1.condition_number == k2.condition_number

    def test_copy_velocity_kernel(self, adata):
        vk1 = VelocityKernel(adata).compute_transition_matrix()
        vk2 = vk1.copy()

        np.testing.assert_array_equal(vk1.transition_matrix.A, vk2.transition_matrix.A)
        np.testing.assert_array_equal(
            vk1.pearson_correlations.A, vk2.pearson_correlations.A
        )
        if vk1._tmats is not None:
            for m1, m2 in zip(vk1._tmats, vk2._tmats):
                np.testing.assert_array_equal(m1.A, m2.A)
        if vk1._pcors is not None:
            for m1, m2 in zip(vk1._pcors, vk2._pcors):
                np.testing.assert_array_equal(m1.A, m2.A)

        assert vk1.params == vk2.params
        assert vk1.backward == vk2.backward

    def test_copy_connectivity_kernel(self, adata):
        ck1 = ConnectivityKernel(adata).compute_transition_matrix()
        ck2 = ck1.copy()

        np.testing.assert_array_equal(ck1.transition_matrix.A, ck2.transition_matrix.A)
        assert ck1.params == ck2.params
        assert ck1.backward == ck2.backward

    def test_copy_palantir_kernel(self, adata):
        pk1 = PalantirKernel(adata).compute_transition_matrix()
        pk2 = pk1.copy()

        np.testing.assert_array_equal(pk1.transition_matrix.A, pk2.transition_matrix.A)
        assert pk1.params == pk2.params
        assert pk1.backward == pk2.backward

    def test_copy_works(self, adata):
        ck1 = ConnectivityKernel(adata)
        ck2 = ck1.copy()
        ck1.compute_transition_matrix()

        assert (
            ck1._transition_matrix is not None
        )  # calling the property would trigger the calculation
        assert ck2._transition_matrix is None


class TestGeneral:
    def test_kernels(self, adata):
        vk = VelocityKernel(adata)

        assert len(vk.kernels) == 1
        assert vk.kernels[0] is vk

    def test_kernels_multiple(self, adata):
        vk = VelocityKernel(adata)
        ck = ConnectivityKernel(adata)
        v = vk + ck

        assert len(v.kernels) == 2
        assert v.kernels == [vk, ck]

    def test_kernels_multiple_constant(self, adata):
        vk = VelocityKernel(adata)
        ck = ConnectivityKernel(adata)
        v = 100 * vk + 42 * ck

        assert len(v.kernels) == 2
        assert v.kernels == [vk, ck]

    def test_no_comp_cond_num(self, adata):
        vk = VelocityKernel(adata).compute_transition_matrix()

        assert vk.condition_number is None

    def test_comp_cond_num(self, adata):
        vk = VelocityKernel(adata, compute_cond_num=True).compute_transition_matrix()

        assert isinstance(vk.condition_number, float)

    def test_comp_cond_num_or_policy(self, adata):
        vk = VelocityKernel(adata, compute_cond_num=True).compute_transition_matrix()
        ck = ConnectivityKernel(
            adata, compute_cond_num=False
        ).compute_transition_matrix()
        v = (vk + ck).compute_transition_matrix()

        assert isinstance(vk.condition_number, float)
        assert ck.condition_number is None
        assert isinstance(v.condition_number, float)


@pytest.mark.xfail(reason="unfinished")
class TestTransitionProbabilities:
    def test_pearson_correlations_fwd(self, adata):
        # test whether pearson correlations in cellrank match those from scvelo, forward case
        backward = False

        # compute pearson correlations using scvelo
        velo_graph = adata.uns["velocity_graph"] + adata.uns["velocity_graph_neg"]

        # compute pearson correlations using cellrank
        vk = VelocityKernel(adata, backward=backward)
        vk.compute_transition_matrix(mode="deterministic")
        pearson_correlations_cr = vk.pearson_correlations

        assert np.allclose((velo_graph - pearson_correlations_cr).data, 0, atol=1e-5)

    def test_pearson_correlations_bwd(self, adata):
        # test whether pearson correlations in cellrank match those from scvelo, backward case
        backward = True

        # compute pearson correlations using scvelo
        velo_graph = (adata.uns["velocity_graph"] + adata.uns["velocity_graph_neg"]).T

        # compute pearson correlations using cellrak
        vk = VelocityKernel(adata, backward=backward)
        vk.compute_transition_matrix(mode="deterministic")
        pearson_correlations_cr = vk.pearson_correlations

        assert np.allclose((velo_graph - pearson_correlations_cr).data, 0, atol=1e-5)

    def test_transition_probabilities_fwd(self, adata):
        # test whether transition probabilities in cellrank match those from scvelo, forward case
        sigma_test = 3

        # compute transition probabilities using the scvelo graph
        velo_graph = adata.uns["velocity_graph"] + adata.uns["velocity_graph_neg"]
        T_scv = np.expm1(velo_graph * sigma_test)
        T_scv.data += 1
        T_scv = _normalize(T_scv)

        # compute transition probabilities using cellrank
        vk = VelocityKernel(adata)
        vk.compute_transition_matrix(softmax_scale=sigma_test, mode="deterministic")
        T_cr = vk.transition_matrix

        # check them for point-wise equality
        assert np.allclose((T_scv - T_cr).data, 0, atol=1e-5)

    def test_transition_probabilities_bwd(self, adata):
        # test whether transition probabilities in cellrank match those from scvelo, backward case
        sigma_test = 3

        # compute transition probabilities using the scvelo graph
        velo_graph = adata.uns["velocity_graph"] + adata.uns["velocity_graph_neg"]
        T_scv = np.expm1(velo_graph.T * sigma_test)
        T_scv.data += 1
        T_scv = _normalize(T_scv)

        # compute transition probabilities using cellrank
        vk = VelocityKernel(adata, backward=True)
        vk.compute_transition_matrix(softmax_scale=sigma_test, mode="deterministic")
        T_cr = vk.transition_matrix

        # check them for point-wise equality
        assert np.allclose((T_scv - T_cr).data, 0, atol=1e-5)
