from typing import Tuple, Callable, Optional

import pickle
import pytest
from copy import copy
from pathlib import Path
from _helpers import (
    bias_knn,
    create_kernels,
    density_normalization,
    jax_not_installed_skip,
    random_transition_matrix,
)

import scanpy as sc
import cellrank as cr
from scanpy import Neighbors
from anndata import AnnData
from cellrank._key import Key
from cellrank.tl._utils import _normalize
from cellrank.ul._utils import _get_neighs, _get_neighs_params
from cellrank.tl.kernels import (
    VelocityKernel,
    CytoTRACEKernel,
    PseudotimeKernel,
    PrecomputedKernel,
    ConnectivityKernel,
)
from cellrank.tl.kernels._base_kernel import (
    Kernel,
    Constant,
    KernelAdd,
    KernelMul,
    _dtype,
    _is_bin_mult,
)
from cellrank.tl.kernels._cytotrace_kernel import CytoTRACEAggregation

import numpy as np
from scipy.sparse import eye as speye
from scipy.sparse import isspmatrix_csr
from pandas.core.dtypes.common import is_bool_dtype, is_integer_dtype

_rtol = 1e-6


class CustomFunc(cr.tl.kernels.SimilaritySchemeABC):
    def __call__(
        self, v: np.ndarray, D: np.ndarray, softmax_scale: float = 1.0
    ) -> Tuple[np.ndarray, np.ndarray]:
        probs, logits = np.zeros((D.shape[0],), dtype=np.float64), np.zeros(
            (D.shape[0],), dtype=np.float64
        )
        probs[0] = 1.0

        return probs, logits


class CustomFuncHessian(CustomFunc):
    def hessian(
        self, v: np.ndarray, D: np.ndarray, _softmax_scale: float = 1.0
    ) -> np.ndarray:
        # should be either (n, g, g) or (n, g), will be (g, g)
        return np.zeros((D.shape[0], v.shape[0], v.shape[0]))


class CustomKernel(Kernel):
    def compute_transition_matrix(
        self, sparse: bool = False, dnorm: bool = False
    ) -> "KernelExpression":
        if sparse:
            tmat = speye(self.adata.n_obs, dtype=np.float32)
        else:
            tmat = np.eye(self.adata.n_obs, dtype=np.float32)

        self._compute_transition_matrix(tmat, density_normalize=dnorm)
        return self

    def copy(self) -> "KernelExpression":
        return copy(self)


class InvalidFuncProbs(cr.tl.kernels.SimilaritySchemeABC):
    def __call__(
        self, v: np.ndarray, D: np.ndarray, _softmax_scale: float = 1.0
    ) -> Tuple[np.ndarray, np.ndarray]:
        return np.ones((D.shape[0],), dtype=np.float64), np.zeros(
            (D.shape[0],), dtype=np.float64
        )


class InvalidFuncHessianShape(CustomFunc):
    def __call__(
        self, v: np.ndarray, D: np.ndarray, _softmax_scale: float = 1.0
    ) -> Tuple[np.ndarray, np.ndarray]:
        probs, logits = np.zeros((D.shape[0],), dtype=np.float64), np.zeros(
            (D.shape[0],), dtype=np.float64
        )
        probs[-1] = 1.0

        return probs, logits

    def hessian(
        self, v: np.ndarray, _D: np.ndarray, _softmax_scale: float = 1.0
    ) -> np.ndarray:
        # should be either (n, g, g) or (n, g), will be (g, g)
        return np.zeros((v.shape[0], v.shape[0]))


class TestInitializeKernel:
    def test_none_transition_matrix(self, adata: AnnData):
        vk = VelocityKernel(adata)
        ck = ConnectivityKernel(adata)
        pk = PseudotimeKernel(adata, time_key="latent_time")

        assert vk._transition_matrix is None
        assert ck._transition_matrix is None
        assert pk._transition_matrix is None

    def test_not_none_transition_matrix_compute(self, adata: AnnData):
        vk = VelocityKernel(adata).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata).compute_transition_matrix()
        pk = PseudotimeKernel(adata, time_key="latent_time").compute_transition_matrix()

        assert vk.transition_matrix is not None
        assert ck.transition_matrix is not None
        assert pk.transition_matrix is not None

    def test_not_none_transition_matrix_accessor(self, adata: AnnData):
        vk = VelocityKernel(adata)
        ck = ConnectivityKernel(adata)
        pk = PseudotimeKernel(adata, time_key="latent_time")

        assert vk.transition_matrix is not None
        assert ck.transition_matrix is not None
        assert pk.transition_matrix is not None

    def test_adding_hidden_constants(self, adata: AnnData):
        k = VelocityKernel(adata) + ConnectivityKernel(adata)

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

    def test_length(self, adata: AnnData):
        k = VelocityKernel(adata) + ConnectivityKernel(adata)
        assert len(k) == 2

    def test_accessor_out_of_range(self, adata: AnnData):
        k = VelocityKernel(adata) + ConnectivityKernel(adata)

        with pytest.raises(IndexError):
            _ = k[2]

    def test_parent(self, adata: AnnData):
        vk = VelocityKernel(adata)
        ck = ConnectivityKernel(adata)
        k = vk + ck

        assert vk._parent._parent is k  # invisible constants
        assert ck._parent._parent is k
        assert k._parent is None

    def test_uninitialized_both(self, adata: AnnData):
        k = VelocityKernel(adata) + ConnectivityKernel(adata)

        with pytest.raises(RuntimeError):
            k.compute_transition_matrix()

    def test_uninitialized_one(self, adata: AnnData):
        k = (
            VelocityKernel(adata)
            + ConnectivityKernel(adata).compute_transition_matrix()
        )

        with pytest.raises(RuntimeError):
            k.compute_transition_matrix()

    def test_initialized(self, adata: AnnData):
        k = (
            VelocityKernel(adata).compute_transition_matrix(softmax_scale=4)
            + ConnectivityKernel(adata).compute_transition_matrix()
        )
        k.compute_transition_matrix()

        assert k.transition_matrix is not None

    def test_invalida_type(self, adata: AnnData):
        with pytest.raises(TypeError):
            _ = None * VelocityKernel(adata)

    def test_negative_constant(self, adata: AnnData):
        with pytest.raises(ValueError):
            _ = -1 * VelocityKernel(adata)

    def test_invalid_constant(self, adata: AnnData):
        with pytest.raises(TypeError):
            _ = Constant(adata, None)

    def test_inversion(self, adata: AnnData):
        c = ConnectivityKernel(adata, backward=False)
        assert not c.backward

        nc = ~c
        assert nc.backward

    def test_inversion_inplace(self, adata: AnnData):
        c = ConnectivityKernel(adata, backward=False)

        assert not c.backward
        _ = ~c
        assert c.backward

    def test_inversion_propagation(self, adata: AnnData):
        c = ConnectivityKernel(adata, backward=False)
        v = VelocityKernel(adata, backward=False)
        k = ~(c + v)

        assert c.backward
        assert v.backward
        assert k.backward

    def test_inversion_recalculation(self, adata: AnnData):
        c = ConnectivityKernel(adata).compute_transition_matrix()
        z = ~(c + c)
        with pytest.raises(RuntimeError):
            z.compute_transition_matrix()

    def test_inversion_preservation_of_constants(self, adata: AnnData):
        c = ConnectivityKernel(adata).compute_transition_matrix()
        a = (3 * c + 1 * c).compute_transition_matrix()
        b = ~a
        c.compute_transition_matrix()

        assert a[0][0].transition_matrix == 3 / 4
        assert b[0][0].transition_matrix == 3 / 4
        assert a[1][0].transition_matrix == 1 / 4
        assert b[1][0].transition_matrix == 1 / 4

    def test_addition_simple(self, adata: AnnData):
        k = VelocityKernel(adata) + ConnectivityKernel(adata)

        assert isinstance(k, KernelAdd)

    def test_multiplication_simple(self, adata: AnnData):
        k = 10 * VelocityKernel(adata)
        c = _is_bin_mult(k)

        assert isinstance(c, Constant)
        assert c.transition_matrix == 10

    def test_multiplication_simple_normalization(self, adata: AnnData):
        k = 10 * VelocityKernel(adata).compute_transition_matrix(softmax_scale=4)
        c = _is_bin_mult(k)

        assert c.transition_matrix == 10

    def test_constant(self, adata: AnnData):
        k = 9 * VelocityKernel(adata) + 1 * ConnectivityKernel(adata)
        c1, c2 = _is_bin_mult(k[0]), _is_bin_mult(k[1])

        assert c1.transition_matrix == 9
        assert c2.transition_matrix == 1

    def test_constant_normalize_2(self, adata: AnnData):
        k = (
            9 * VelocityKernel(adata).compute_transition_matrix(softmax_scale=4)
            + 1 * ConnectivityKernel(adata).compute_transition_matrix()
        )
        k.compute_transition_matrix()
        c1, c2 = _is_bin_mult(k[0]), _is_bin_mult(k[1])

        assert c1.transition_matrix == 9 / 10
        assert c2.transition_matrix == 1 / 10

    def test_constant_normalize_3(self, adata: AnnData):
        k = (
            VelocityKernel(adata).compute_transition_matrix(softmax_scale=4)
            + ConnectivityKernel(adata).compute_transition_matrix()
            + ConnectivityKernel(adata).compute_transition_matrix()
        )
        k.compute_transition_matrix()
        c1, c2, c3 = _is_bin_mult(k[0]), _is_bin_mult(k[1]), _is_bin_mult(k[2])

        assert c1.transition_matrix == 1 / 3
        assert c2.transition_matrix == 1 / 3
        assert c3.transition_matrix == 1 / 3

    def test_constant_wrong_parentheses(self, adata: AnnData):
        k = VelocityKernel(adata).compute_transition_matrix(softmax_scale=4) + (
            ConnectivityKernel(adata).compute_transition_matrix()
            + ConnectivityKernel(adata).compute_transition_matrix()
        )
        k.compute_transition_matrix()
        c1, c2, c3 = _is_bin_mult(k[0]), _is_bin_mult(k[1]), _is_bin_mult(k[2])

        assert c1.transition_matrix == 1 / 3
        assert c2.transition_matrix == 1 / 3
        assert c3.transition_matrix == 1 / 3

    def test_constant_correct_parentheses(self, adata: AnnData):
        k = 1 * VelocityKernel(adata).compute_transition_matrix(softmax_scale=4) + 1 * (
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

    def test_adaptive_kernel_constants(self, adata: AnnData):
        ck1 = ConnectivityKernel(adata).compute_transition_matrix()
        ck1._mat_scaler = np.random.normal(size=(adata.n_obs, adata.n_obs))

        ck2 = ConnectivityKernel(adata).compute_transition_matrix()
        ck2._mat_scaler = np.random.normal(size=(adata.n_obs, adata.n_obs))

        k = (3 * ck1) ^ (1 * ck2)
        k.compute_transition_matrix()

        assert k[0][0]._value == 3 / 4
        assert k[1][0]._value == 1 / 4

    def test_adaptive_kernel_complex(self, adata: AnnData):
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

    def test_repr(self, adata: AnnData):
        rpr = repr(VelocityKernel(adata))

        assert rpr == f"<{VelocityKernel.__name__}>"

    def test_repr_inv(self, adata: AnnData):
        rpr = repr(~VelocityKernel(adata))

        assert rpr == f"~<{VelocityKernel.__name__}>"

    def test_repr_inv_comb(self, adata: AnnData):
        rpr = repr(~(VelocityKernel(adata) + ConnectivityKernel(adata)))

        assert (
            rpr
            == f"~((1 * <{VelocityKernel.__name__}>) + (1 * <{ConnectivityKernel.__name__}>))"
        )

    def test_str_repr_equiv_no_transition_matrix(self, adata: AnnData):
        vk = VelocityKernel(adata)
        string = str(vk)
        rpr = repr(vk)

        assert string == rpr
        assert string == f"<{VelocityKernel.__name__}>"

    def test_str(self, adata: AnnData):
        string = str(ConnectivityKernel(adata).compute_transition_matrix())

        assert (
            string == f"<{ConnectivityKernel.__name__}[dnorm=True, key=connectivities]>"
        )

    def test_str_inv(self, adata: AnnData):
        string = str(
            ConnectivityKernel(adata, backward=True).compute_transition_matrix()
        )

        assert (
            string
            == f"~<{ConnectivityKernel.__name__}[dnorm=True, key=connectivities]>"
        )

    def test_combination_correct_parameters(self, adata: AnnData):
        from cellrank.tl.kernels import CosineScheme

        k = VelocityKernel(adata).compute_transition_matrix(
            softmax_scale=4,
            seed=42,
            scheme="cosine",
        ) + (
            ConnectivityKernel(adata).compute_transition_matrix(density_normalize=False)
            + ConnectivityKernel(adata).compute_transition_matrix(
                density_normalize=True
            )
        )
        k.compute_transition_matrix()

        assert isinstance(k.params, dict)
        assert len(k.params) == 3
        assert {"dnorm": True, "key": "connectivities"} in k.params.values()
        assert {"dnorm": False, "key": "connectivities"} in k.params.values()
        assert {
            "softmax_scale": 4,
            "mode": "deterministic",
            "seed": 42,
            "scheme": str(CosineScheme()),
        } in k.params.values()


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

    def test_precomputed_from_kernel_no_transition(self, adata: AnnData):
        vk = VelocityKernel(adata)

        with pytest.raises(ValueError):
            PrecomputedKernel(vk)

    @pytest.mark.parametrize(
        "clazz",
        [
            ConnectivityKernel,
            VelocityKernel,
            PseudotimeKernel,
            CytoTRACEKernel,
            PrecomputedKernel,
        ],
    )
    @pytest.mark.parametrize("key_added", [None, "foo"])
    def test_kernel_reads_correct_connectivities(
        self, adata: AnnData, key_added: Optional[str], clazz: type
    ):
        del adata.uns["neighbors"]
        del adata.obsp["connectivities"]
        del adata.obsp["distances"]

        sc.pp.neighbors(adata, key_added=key_added)
        kwargs = {"adata": adata, "conn_key": key_added}

        if clazz == PseudotimeKernel:
            kwargs["time_key"] = "latent_time"
        elif clazz == PrecomputedKernel:
            adata.obsp["foo"] = np.eye(adata.n_obs)
            kwargs["transition_matrix"] = "foo"
        conn = (
            adata.obsp["connectivities"]
            if key_added is None
            else adata.obsp[f"{key_added}_connectivities"]
        )

        k = clazz(**kwargs)

        if isinstance(k, PrecomputedKernel):
            assert k._conn is None
        else:
            np.testing.assert_array_equal(k._conn.A, conn.A)

    def test_precomputed_from_kernel(self, adata: AnnData):
        vk = VelocityKernel(adata).compute_transition_matrix(
            mode="deterministic",
            softmax_scale=4,
        )

        pk = PrecomputedKernel(vk)
        pk.write_to_adata()

        assert pk.adata is vk.adata
        assert pk._origin == str(vk).strip("~<>")
        assert pk.params is not vk.params
        assert pk.params == vk.params
        assert pk.transition_matrix is not vk.transition_matrix
        np.testing.assert_array_equal(pk.transition_matrix.A, vk.transition_matrix.A)

    def test_precomputed_no_adata(self):
        pk = PrecomputedKernel(random_transition_matrix(50))
        pk.write_to_adata()

        assert isinstance(pk.adata, AnnData)
        assert pk._origin == "'array'"
        assert pk.adata.shape == (50, 1)
        assert pk.adata.obs.shape == (50, 0)
        assert pk.adata.var.shape == (1, 0)
        assert "T_fwd_params" in pk.adata.uns.keys()
        assert pk.adata.uns["T_fwd_params"] == {"params": pk.params}
        np.testing.assert_array_equal(
            pk.adata.obsp["T_fwd"].toarray(), pk.transition_matrix.toarray()
        )

    def test_precomputed_different_adata(self, adata: AnnData):
        vk = VelocityKernel(adata).compute_transition_matrix(
            mode="deterministic", softmax_scale=4
        )
        bdata = adata.copy()

        pk = PrecomputedKernel(vk, adata=bdata)

        assert pk.adata is adata
        assert pk.adata is vk.adata
        assert pk.adata is not bdata

    def test_precomputed_adata_origin(self, adata: AnnData):
        vk = VelocityKernel(adata).compute_transition_matrix(
            mode="deterministic", softmax_scale=4
        )
        vk.write_to_adata("foo")

        pk = PrecomputedKernel("foo", adata=adata)

        assert pk._origin == "adata.obsp['foo']"

    def test_precomputed_adata(self, adata: AnnData):
        pk = PrecomputedKernel(random_transition_matrix(adata.n_obs), adata=adata)

        assert pk.adata is adata

    def test_precomputed_transition_matrix(self, adata: AnnData):
        mat = random_transition_matrix(adata.n_obs)
        pk = PrecomputedKernel(mat)

        np.testing.assert_array_equal(mat, pk.transition_matrix.toarray())

    def test_precomputed_sum(self, adata: AnnData):
        mat = random_transition_matrix(adata.n_obs)
        pk = PrecomputedKernel(mat)
        vk = VelocityKernel(adata).compute_transition_matrix(softmax_scale=4)

        expected = (0.5 * vk.transition_matrix) + (0.5 * pk.transition_matrix)
        actual = (pk + vk).compute_transition_matrix()

        np.testing.assert_array_almost_equal(
            expected.toarray(), actual.transition_matrix.toarray()
        )

    @pytest.mark.parametrize("dnorm", [False, True])
    @pytest.mark.parametrize("sparse", [False, True])
    def test_custom_preserves_type(self, adata: AnnData, sparse: bool, dnorm: bool):
        c = CustomKernel(adata).compute_transition_matrix(sparse=sparse, dnorm=dnorm)

        if sparse:
            assert isspmatrix_csr(c.transition_matrix)
        else:
            assert isinstance(c.transition_matrix, np.ndarray)

        assert c.transition_matrix.dtype == _dtype

    def test_write_adata(self, adata: AnnData):
        vk = VelocityKernel(adata).compute_transition_matrix(softmax_scale=4)
        vk.write_to_adata()

        assert adata is vk.adata
        assert "T_fwd_params" in adata.uns.keys()
        np.testing.assert_array_equal(
            adata.obsp["T_fwd"].toarray(), vk.transition_matrix.toarray()
        )

    def test_write_adata_key(self, adata: AnnData):
        vk = VelocityKernel(adata).compute_transition_matrix(softmax_scale=4)
        vk.write_to_adata(key="foo")

        assert adata is vk.adata
        assert "foo_params" in adata.uns.keys()
        np.testing.assert_array_equal(
            adata.obsp["foo"].toarray(), vk.transition_matrix.toarray()
        )

    @pytest.mark.parametrize("mode", ["deterministic", "stochastic"])
    def test_vk_row_normalized(self, adata: AnnData, mode: str):
        if mode == "stochastic":
            pytest.importorskip("jax")
            pytest.importorskip("jaxlib")
        vk = VelocityKernel(adata)
        vk.compute_transition_matrix(mode="stochastic", softmax_scale=4)

        np.testing.assert_allclose(vk.transition_matrix.sum(1), 1, rtol=_rtol)

    # only to 15 because in kernel, if a row sums to 0, abs. states are created
    # this happens because `k_thresh = frac_to_keep = 0`
    @pytest.mark.parametrize("k", range(1, 15))
    def test_pseudotime_frac_to_keep(self, adata: AnnData, k: int):
        conn = _get_neighs(adata, "connectivities")
        n_neighbors = _get_neighs_params(adata)["n_neighbors"]
        pseudotime = adata.obs["latent_time"]
        k_thresh = max(0, min(int(np.floor(n_neighbors / k)) - 1, 30))
        frac_to_keep = k_thresh / float(n_neighbors)

        conn_biased = bias_knn(
            conn.copy(), pseudotime, n_neighbors, k=k, frac_to_keep=frac_to_keep
        )
        T_1 = _normalize(conn_biased)

        pk = PseudotimeKernel(adata, time_key="latent_time").compute_transition_matrix(
            frac_to_keep=frac_to_keep,
            threshold_scheme="hard",
        )
        T_2 = pk.transition_matrix

        np.testing.assert_allclose(T_1.A, T_2.A, rtol=_rtol)

    def test_pseudotime_parallelize(self, adata: AnnData):
        pk1 = PseudotimeKernel(adata, time_key="latent_time").compute_transition_matrix(
            n_jobs=None
        )
        pk2 = PseudotimeKernel(adata, time_key="latent_time").compute_transition_matrix(
            n_jobs=2
        )

        np.testing.assert_allclose(
            pk1.transition_matrix.A, pk2.transition_matrix.A, rtol=_rtol
        )

    def test_pseudotime_inverse(self, adata: AnnData):
        pk = PseudotimeKernel(adata, time_key="latent_time")
        pt = pk.pseudotime.copy()

        pk_inv = ~pk

        assert pk_inv is pk
        assert pk_inv.backward
        np.testing.assert_allclose(pt, 1 - pk_inv.pseudotime)

    @pytest.mark.parametrize("mode", ["deterministic", "stochastic", "sampling"])
    def test_manual_combination(self, adata: AnnData, mode: str):
        if mode == "stochastic":
            pytest.importorskip("jax")
            pytest.importorskip("jaxlib")
        vk = VelocityKernel(adata).compute_transition_matrix(mode=mode, softmax_scale=4)
        ck = ConnectivityKernel(adata).compute_transition_matrix()

        T_vk = vk.transition_matrix
        T_ck = ck.transition_matrix
        T_comb_manual = 0.8 * T_vk + 0.2 * T_ck

        comb_kernel = 0.8 * vk + 0.2 * ck
        T_comb_kernel = comb_kernel.transition_matrix

        np.testing.assert_allclose(T_comb_kernel.A, T_comb_manual.A, rtol=_rtol)

    def test_manual_combination_no_precomputed(self, adata: AnnData):
        density_normalize = False

        vk = VelocityKernel(adata).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata).compute_transition_matrix(
            density_normalize=density_normalize
        )

        T_vk = vk.transition_matrix
        T_ck = ck.transition_matrix
        T_comb_manual = 0.8 * T_vk + 0.2 * T_ck

        vk = VelocityKernel(adata).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata).compute_transition_matrix(
            density_normalize=density_normalize
        )
        comb_kernel = 0.8 * vk + 0.2 * ck
        comb_kernel.compute_transition_matrix()
        T_comb_kernel = comb_kernel.transition_matrix

        np.testing.assert_allclose(T_comb_manual.A, T_comb_kernel.A, rtol=_rtol)

    def test_manual_combination_backward(self, adata: AnnData):
        backward, density_normalize = True, False

        vk = VelocityKernel(adata, backward=backward).compute_transition_matrix(
            softmax_scale=4
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

    def test_manual_combination_backward_dense_norm(self, adata: AnnData):
        backward, density_normalize = True, True

        vk = VelocityKernel(adata, backward=backward).compute_transition_matrix(
            softmax_scale=4
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

    def compare_with_scanpy_density_normalize(self, adata: AnnData):
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

    def compare_with_scanpy(self, adata: AnnData):
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

    def test_connectivities_key_kernel(self, adata: AnnData):
        key = "foobar"
        assert key not in adata.obsp
        adata.obsp[key] = np.eye(adata.n_obs)

        ck = ConnectivityKernel(adata, conn_key=key).compute_transition_matrix()
        T_cr = ck.transition_matrix

        assert key == ck.params["key"]
        np.testing.assert_array_equal(T_cr, adata.obsp[key])

        del adata.obsp[key]


class TestKernelAddition:
    def test_simple_addition(self, adata: AnnData):
        vk, ck = create_kernels(adata)  # diagonal + upper diag

        k = (vk + ck).compute_transition_matrix()
        expected = np.eye(adata.n_obs) * 0.75 + np.eye(adata.n_obs, k=1) * 0.25
        expected[-1, -1] = 1

        np.testing.assert_allclose(k.transition_matrix.A, expected)

    def test_addtion_with_constant(self, adata: AnnData):
        vk, ck = create_kernels(adata)  # diagonal + upper diag

        k = (2 * vk + 3 * ck).compute_transition_matrix()
        expected = (
            np.eye(adata.n_obs) * (2 / 5)
            + np.eye(adata.n_obs) * (3 / 5) * 0.5
            + np.eye(adata.n_obs, k=1) * (3 / 5) * 0.5
        )
        expected[-1, -1] = 1

        np.testing.assert_allclose(k.transition_matrix.A, expected)

    def test_addition_3_kernels(self, adata: AnnData):
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

    def test_addition_adaptive(self, adata: AnnData):
        adata.obsp["velocity_variances"] = vv = np.random.random(
            size=(adata.n_obs, adata.n_obs)
        )
        adata.obsp["connectivity_variances"] = cv = np.random.random(
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

    def test_addition_adaptive_constants(self, adata: AnnData):
        a, b = np.random.uniform(0, 10, 2)
        s = a + b
        adata.obsp["velocity_variances"] = vv = np.random.random(
            size=(adata.n_obs, adata.n_obs)
        )
        adata.obsp["connectivity_variances"] = cv = np.random.random(
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

    def test_addition_adaptive_wrong_variances(self, adata: AnnData):
        a, b = np.random.uniform(0, 10, 2)
        s = a + b
        adata.obsp["velocity_variances"] = np.random.random(
            size=(adata.n_obs, adata.n_obs)
        )
        adata.obsp["connectivity_variances"] = np.random.random(
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

    def test_addition_adaptive_4_kernels(self, adata: AnnData):
        a, b, c, d = np.random.uniform(0, 10, 4)
        s = a + b + c + d
        adata.obsp["velocity_variances"] = vv = np.random.random(
            size=(adata.n_obs, adata.n_obs)
        )
        adata.obsp["connectivity_variances"] = cv = np.random.random(
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
    def test_copy_simple(self, adata: AnnData):
        vk1 = VelocityKernel(adata)
        vk2 = vk1.copy()

        assert vk1 is not vk2

    def test_copy_no_adata_copy(self, adata: AnnData):
        vk1 = VelocityKernel(adata)
        vk2 = vk1.copy()

        assert vk1.adata is adata
        assert vk2.adata is adata

    def test_copy_transition_matrix(self, adata: AnnData):
        vk1 = VelocityKernel(adata).compute_transition_matrix(softmax_scale=4)
        vk2 = vk1.copy()

        np.testing.assert_array_equal(vk1.transition_matrix.A, vk2.transition_matrix.A)

    def test_copy_params(self, adata: AnnData):
        vk1 = VelocityKernel(adata).compute_transition_matrix(softmax_scale=4)
        vk2 = vk1.copy()

        assert vk1.params == vk2.params

    def test_copy_cond_num(self, adata: AnnData):
        for KernelClass in [
            VelocityKernel,
            ConnectivityKernel,
            PseudotimeKernel,
            PrecomputedKernel,
        ]:
            if KernelClass is PrecomputedKernel:
                k1 = KernelClass(
                    random_transition_matrix(adata.n_obs), compute_cond_num=True
                )
            elif KernelClass is VelocityKernel:
                k1 = KernelClass(
                    adata, compute_cond_num=True
                ).compute_transition_matrix(softmax_scale=4)
            else:
                k1 = KernelClass(
                    adata, compute_cond_num=True
                ).compute_transition_matrix()
            k2 = k1.copy()

            assert k1.condition_number == k2.condition_number

    def test_copy_velocity_kernel(self, adata: AnnData):
        vk1 = VelocityKernel(adata).compute_transition_matrix(softmax_scale=4)
        vk2 = vk1.copy()

        np.testing.assert_array_equal(vk1.transition_matrix.A, vk2.transition_matrix.A)
        np.testing.assert_array_equal(vk1.logits.A, vk2.logits.A)

        assert vk1.params == vk2.params
        assert vk1.backward == vk2.backward

    def test_copy_connectivity_kernel(self, adata: AnnData):
        ck1 = ConnectivityKernel(adata).compute_transition_matrix()
        ck2 = ck1.copy()

        np.testing.assert_array_equal(ck1.transition_matrix.A, ck2.transition_matrix.A)
        assert ck1.params == ck2.params
        assert ck1.backward == ck2.backward

    def test_copy_palantir_kernel(self, adata: AnnData):
        pk1 = PseudotimeKernel(adata).compute_transition_matrix()
        pk2 = pk1.copy()

        np.testing.assert_array_equal(pk1.transition_matrix.A, pk2.transition_matrix.A)
        assert pk1.params == pk2.params
        assert pk1.backward == pk2.backward

    def test_copy_works(self, adata: AnnData):
        ck1 = ConnectivityKernel(adata)
        ck2 = ck1.copy()
        ck1.compute_transition_matrix()

        assert (
            ck1._transition_matrix is not None
        )  # calling the property would trigger the calculation
        assert ck2._transition_matrix is None


class TestGeneral:
    def test_kernels(self, adata: AnnData):
        vk = VelocityKernel(adata)

        assert len(vk.kernels) == 1
        assert vk.kernels[0] is vk

    def test_kernels_multiple(self, adata: AnnData):
        vk = VelocityKernel(adata)
        ck = ConnectivityKernel(adata)
        v = vk + ck

        assert len(v.kernels) == 2
        assert vk in v.kernels
        assert ck in v.kernels

    def test_kernels_multiple_constant(self, adata: AnnData):
        vk = VelocityKernel(adata)
        ck = ConnectivityKernel(adata)
        v = 100 * vk + 42 * ck

        assert len(v.kernels) == 2
        assert vk in v.kernels
        assert ck in v.kernels

    def test_kernels_unique(self, adata: AnnData):
        vk = VelocityKernel(adata)
        v = vk + vk + vk + vk

        assert len(v.kernels) == 1
        assert v.kernels[0] is vk

    def test_no_comp_cond_num(self, adata: AnnData):
        vk = VelocityKernel(adata).compute_transition_matrix(softmax_scale=4)

        assert vk.condition_number is None

    def test_comp_cond_num(self, adata: AnnData):
        vk = VelocityKernel(adata, compute_cond_num=True).compute_transition_matrix(
            softmax_scale=4
        )

        assert isinstance(vk.condition_number, float)

    def test_comp_cond_num_or_policy(self, adata: AnnData):
        vk = VelocityKernel(adata, compute_cond_num=True).compute_transition_matrix(
            softmax_scale=4
        )
        ck = ConnectivityKernel(
            adata, compute_cond_num=False
        ).compute_transition_matrix()
        v = (vk + ck).compute_transition_matrix()

        assert isinstance(vk.condition_number, float)
        assert ck.condition_number is None
        assert isinstance(v.condition_number, float)


class TestTransitionProbabilities:
    def test_pearson_correlations_fwd(self, adata: AnnData):
        # test whether pearson correlations in cellrank match those from scvelo, forward case
        backward = False

        # compute pearson correlations using scvelo
        velo_graph = adata.obsp["velocity_graph"] + adata.obsp["velocity_graph_neg"]

        # compute pearson correlations using cellrank
        vk = VelocityKernel(adata, backward=backward)
        vk.compute_transition_matrix(mode="deterministic", softmax_scale=4)
        pearson_correlations_cr = vk.logits

        pc_r = velo_graph.copy()
        pc_r.data = np.array(pearson_correlations_cr[(velo_graph != 0)]).squeeze()

        assert np.max(np.abs((pc_r - velo_graph).data)) < _rtol

    def test_pearson_correlations_bwd(self, adata: AnnData):
        # test whether pearson correlations in cellrank match those from scvelo, backward case
        backward = True

        # compute pearson correlations using scvelo
        velo_graph = (adata.obsp["velocity_graph"] + adata.obsp["velocity_graph_neg"]).T

        # compute pearson correlations using cellrak
        vk = VelocityKernel(adata, backward=backward)
        vk.compute_transition_matrix(
            mode="deterministic", backward_mode="transpose", softmax_scale=4
        )
        pearson_correlations_cr = vk.logits

        pc_r = velo_graph.copy()
        pc_r.data = np.array(pearson_correlations_cr[(velo_graph != 0)]).squeeze()

        assert np.max(np.abs((pc_r - velo_graph.T).data)) < _rtol

    def test_transition_probabilities_fwd(self, adata: AnnData):
        # test whether transition probabilities in cellrank match those from scvelo, forward case
        sigma_test = 3

        # compute transition probabilities using cellrank
        vk = VelocityKernel(adata)
        vk.compute_transition_matrix(softmax_scale=sigma_test, mode="deterministic")
        T_cr = vk.transition_matrix

        pearson_correlation = vk.logits
        T_exp = np.expm1(pearson_correlation * sigma_test)
        T_exp.data += 1
        T_exp = _normalize(T_exp)

        np.testing.assert_allclose(T_exp.A, T_cr.A)  # don't use data, can be reordered

    def test_transition_probabilities_bwd(self, adata: AnnData):
        # test whether transition probabilities in cellrank match those from scvelo, backward case
        sigma_test = 3

        # compute transition probabilities using cellrank
        vk = VelocityKernel(adata, backward=True)
        vk.compute_transition_matrix(softmax_scale=sigma_test, mode="deterministic")
        T_cr = vk.transition_matrix

        pearson_correlation = vk.logits
        T_exp = np.expm1(pearson_correlation * sigma_test)
        T_exp.data += 1
        T_exp = _normalize(T_exp)

        np.testing.assert_allclose(T_exp.A, T_cr.A)  # don't use data, can be reordered

    def test_estimate_softmax_scale(self, adata: AnnData):
        vk = VelocityKernel(adata)
        vk.compute_transition_matrix(
            mode="deterministic", show_progress_bar=False, softmax_scale=None
        )

        assert isinstance(vk.params["softmax_scale"], float)


class TestMonteCarlo:
    def test_mc_and_mc_fwd_1k(self, adata: AnnData):
        vk_mc = VelocityKernel(adata, backward=False)
        vk_mc.compute_transition_matrix(
            mode="monte_carlo",
            show_progress_bar=False,
            n_samples=1000,
            n_jobs=4,
            softmax_scale=4,
        )

        vk_s = VelocityKernel(adata, backward=False)
        vk_s.compute_transition_matrix(
            mode="monte_carlo",
            show_progress_bar=False,
            n_samples=1000,
            n_jobs=4,
            softmax_scale=4,
        )

        val = np.mean(
            np.abs(vk_mc.transition_matrix.data - vk_s.transition_matrix.data)
        )
        assert val < 1e-5, val

    def test_monte_carlo_5k(self, adata: AnnData):
        vk_mc = VelocityKernel(adata, backward=False)
        vk_mc.compute_transition_matrix(
            mode="monte_carlo",
            show_progress_bar=False,
            n_samples=5000,
            n_jobs=4,
            softmax_scale=4,
            seed=42,
        )

        vk_s = VelocityKernel(adata, backward=False)
        vk_s.compute_transition_matrix(
            mode="monte_carlo",
            show_progress_bar=False,
            n_samples=5000,
            n_jobs=4,
            softmax_scale=4,
            seed=43,
        )

        val = np.mean(
            np.abs(vk_mc.transition_matrix.data - vk_s.transition_matrix.data)
        )
        assert val < 1e-5, val

    @jax_not_installed_skip
    def test_monte_carlo_and_stochastic(self, adata: AnnData):
        vk_mc = VelocityKernel(adata, backward=False)
        vk_mc.compute_transition_matrix(
            mode="monte_carlo",
            show_progress_bar=False,
            n_samples=1000,
            n_jobs=4,
            softmax_scale=4,
        )

        vk_s = VelocityKernel(adata, backward=False)
        vk_s.compute_transition_matrix(
            mode="stochastic", show_progress_bar=False, n_jobs=4, softmax_scale=4
        )

        val = np.mean(
            np.abs(vk_mc.transition_matrix.data - vk_s.transition_matrix.data)
        )
        assert val < 1e-3, val


class TestVelocityScheme:
    def test_invalid_string_key(self, adata: AnnData):
        vk = VelocityKernel(adata)
        with pytest.raises(ValueError):
            vk.compute_transition_matrix(scheme="foobar")

    def test_not_callable(self, adata: AnnData):
        vk = VelocityKernel(adata)
        with pytest.raises(
            TypeError, match="Expected `scheme` to be a function, found"
        ):
            vk.compute_transition_matrix(scheme=1311)

    def test_custom_function_not_sum_to_1(self, adata: AnnData):
        vk = VelocityKernel(adata)
        with pytest.raises(ValueError, match=r"Matrix is not row-stochastic."):
            vk.compute_transition_matrix(scheme=InvalidFuncProbs())

    def test_custom_function_invalid_hessian(self, adata: AnnData):
        vk = VelocityKernel(adata)
        with pytest.raises(ValueError, match=r"Expected full Hessian matrix"):
            vk.compute_transition_matrix(
                mode="stochastic", scheme=InvalidFuncHessianShape(), softmax_scale=4
            )

    @pytest.mark.parametrize("backward", [True, False])
    def test_implementations_differ(self, adata: AnnData, backward: bool):
        vk_dot = VelocityKernel(adata, backward=backward)
        vk_dot.compute_transition_matrix(
            mode="deterministic", softmax_scale=4, scheme="dot_product"
        )
        vk_cos = VelocityKernel(adata, backward=backward)
        vk_cos.compute_transition_matrix(
            mode="deterministic", softmax_scale=4, scheme="cosine"
        )
        vk_cor = VelocityKernel(adata, backward=backward)
        vk_cor.compute_transition_matrix(
            mode="deterministic", softmax_scale=4, scheme="correlation"
        )

        np.testing.assert_allclose(vk_dot.transition_matrix.sum(1), 1.0)
        np.testing.assert_allclose(vk_cor.transition_matrix.sum(1), 1.0)
        np.testing.assert_allclose(vk_cor.transition_matrix.sum(1), 1.0)

        assert not np.allclose(vk_dot.transition_matrix.A, vk_cos.transition_matrix.A)
        assert not np.allclose(vk_cos.transition_matrix.A, vk_cor.transition_matrix.A)
        assert not np.allclose(vk_cor.transition_matrix.A, vk_dot.transition_matrix.A)

    @pytest.mark.parametrize(
        "key,fn",
        zip(
            ["dot_product", "cosine", "correlation"],
            [
                cr.tl.kernels.DotProductScheme(),
                cr.tl.kernels.CosineScheme(),
                cr.tl.kernels.CorrelationScheme(),
            ],
        ),
    )
    def test_function_and_string_key(self, adata: AnnData, key: str, fn: Callable):
        vk_k = VelocityKernel(adata)
        vk_fn = VelocityKernel(adata)

        vk_k.compute_transition_matrix(
            mode="deterministic", softmax_scale=4, scheme=key
        )
        vk_fn.compute_transition_matrix(
            mode="deterministic", softmax_scale=4, scheme=fn
        )

        np.testing.assert_allclose(vk_k.transition_matrix.A, vk_fn.transition_matrix.A)

    @pytest.mark.parametrize("backward", [True, False])
    def test_custom_function(self, adata: AnnData, backward: bool):
        vk = VelocityKernel(adata, backward=backward)
        vk.compute_transition_matrix(
            mode="deterministic", softmax_scale=4, scheme=CustomFuncHessian()
        )

        assert vk.params["scheme"] == str(CustomFuncHessian())

    def test_custom_function_stochastic_no_hessian(self, adata: AnnData):
        vk = VelocityKernel(adata)
        vk.compute_transition_matrix(
            mode="stochastic", scheme=CustomFunc(), softmax_scale=4, n_samples=10
        )

        assert vk.params["mode"] == "monte_carlo"
        assert vk.params["scheme"] == str(CustomFunc())


class TestComputeProjection:
    def test_no_transition_matrix(self, adata: AnnData):
        with pytest.raises(RuntimeError, match=r"Compute transition matrix first as"):
            cr.tl.kernels.ConnectivityKernel(adata).compute_projection()

    def test_no_basis(self, adata: AnnData):
        ck = cr.tl.kernels.ConnectivityKernel(adata).compute_transition_matrix()
        with pytest.raises(KeyError, match=r"Unable to find a basis in"):
            ck.compute_projection(basis="foo")

    def test_basis_prefix(self, adata: AnnData):
        ck = cr.tl.kernels.ConnectivityKernel(adata).compute_transition_matrix()
        ck.compute_projection(basis="X_umap")

    @pytest.mark.parametrize("write_first", [True, False])
    def test_write_to_adata(self, adata: AnnData, write_first: bool):
        ck = cr.tl.kernels.ConnectivityKernel(adata).compute_transition_matrix()
        if write_first:
            ck.write_to_adata()
            ck.compute_projection(basis="umap")
        else:
            ck.compute_projection(basis="umap")
            ck.write_to_adata()

        assert adata.uns[Key.uns.kernel(ck.backward) + "_params"] == {
            "params": ck.params,
            "embeddings": ["umap"],
        }

    @pytest.mark.parametrize("key_added", [None, "foo"])
    def test_key_added(self, adata: AnnData, key_added: Optional[str]):
        ck = cr.tl.kernels.ConnectivityKernel(adata).compute_transition_matrix()
        ck.compute_projection(basis="umap", copy=False, key_added=key_added)

        key = Key.uns.kernel(ck.backward, key=key_added)
        ukey = f"{key}_params"
        key = f"{key}_umap"

        assert adata.uns[ukey] == {"embeddings": ["umap"]}
        np.testing.assert_array_equal(adata.obsm[key].shape, adata.obsm["X_umap"].shape)

    @pytest.mark.parametrize("copy", [True, False])
    def test_copy(self, adata: AnnData, copy: bool):
        ck = cr.tl.kernels.ConnectivityKernel(adata).compute_transition_matrix()
        res = ck.compute_projection(basis="umap", copy=copy)

        if copy:
            assert isinstance(res, np.ndarray)
            np.testing.assert_array_equal(res.shape, adata.obsm["X_umap"].shape)
        else:
            assert res is None
            key = Key.uns.kernel(ck.backward) + "_umap"
            np.testing.assert_array_equal(
                adata.obsm[key].shape, adata.obsm["X_umap"].shape
            )

    def test_nan_in_embedding(self, adata: AnnData):
        adata.obsm["X_umap"][-1] = np.nan

        ck = cr.tl.kernels.ConnectivityKernel(adata).compute_transition_matrix()
        res = ck.compute_projection(basis="umap", copy=True)

        assert not np.all(np.isnan(res))
        assert np.all(np.isnan(res[-1, :]))


class TestPseudotimeKernelScheme:
    def test_invalid_scheme(self, adata: AnnData):
        pk = PseudotimeKernel(adata)
        with pytest.raises(ValueError, match="foo"):
            pk.compute_transition_matrix(threshold_scheme="foo")

    def test_invalid_custom_scheme(self, adata: AnnData):
        pk = PseudotimeKernel(adata)
        with pytest.raises(ValueError, match="Expected row of shape"):
            pk.compute_transition_matrix(
                threshold_scheme=lambda cpt, npt, ndist: np.ones(
                    (len(ndist) - 1), dtype=np.float64
                ),
            )

    def test_custom_scheme(self, adata: AnnData):
        pk = PseudotimeKernel(adata)
        pk.compute_transition_matrix(
            threshold_scheme=lambda cpt, npt, ndist: np.ones(
                (len(ndist)), dtype=np.float64
            ),
        )

        np.testing.assert_allclose(pk.transition_matrix.sum(1), 1.0)
        for row in pk.transition_matrix:
            np.testing.assert_allclose(row.data, 1 / len(row.data))

    @pytest.mark.parametrize("scheme", ["hard", "soft"])
    def test_scheme(self, adata: AnnData, scheme: str):
        pk = PseudotimeKernel(adata)
        pk.compute_transition_matrix(
            threshold_scheme=scheme, frac_to_keep=0.3, b=10, nu=0.5
        )

        np.testing.assert_allclose(pk.transition_matrix.sum(1), 1.0)
        assert pk.params["scheme"] == scheme
        if scheme == "hard":
            assert pk.params["frac_to_keep"] == 0.3
            assert "b" not in pk.params
            assert "nu" not in pk.params
        elif scheme == "soft":
            assert pk.params["b"] == 10
            assert pk.params["nu"] == 0.5
            assert "k" not in pk.params


class TestCytoTRACEKernel:
    @pytest.mark.parametrize("layer", ["X", "Ms", "foo"])
    def test_layer(self, adata: AnnData, layer: str):
        if layer == "foo":
            with pytest.raises(KeyError, match=layer):
                _ = CytoTRACEKernel(adata, layer=layer)
        else:
            _ = CytoTRACEKernel(adata, layer=layer)
            assert adata.uns[Key.cytotrace("params")]["layer"] == layer

    @pytest.mark.parametrize("agg", list(CytoTRACEAggregation))
    def test_aggregation(self, adata: AnnData, agg: CytoTRACEAggregation):
        _ = CytoTRACEKernel(adata, aggregation=agg)
        assert adata.uns[Key.cytotrace("params")]["aggregation"] == agg

    @pytest.mark.parametrize("use_raw", [False, True])
    def test_raw(self, adata: AnnData, use_raw: bool):
        _ = CytoTRACEKernel(adata, use_raw=use_raw)
        assert adata.uns[Key.cytotrace("params")]["use_raw"] == (
            adata.raw.n_vars == adata.n_vars if use_raw else False
        )

    def test_correct_class(self, adata: AnnData):
        k = CytoTRACEKernel(adata)
        assert isinstance(k, PseudotimeKernel)
        assert k._time_key == Key.cytotrace("pseudotime")

    def test_writes_params(self, adata: AnnData):
        k = CytoTRACEKernel(adata, use_raw=False, layer="X", aggregation="mean")

        assert adata.uns[Key.cytotrace("params")] == {
            "layer": "X",
            "aggregation": "mean",
            "use_raw": False,
        }

        assert np.all(adata.var[Key.cytotrace("gene_corr")] <= 1.0)
        assert np.all(-1 <= adata.var[Key.cytotrace("gene_corr")])
        assert is_bool_dtype(adata.var[Key.cytotrace("correlates")])
        assert adata.var[Key.cytotrace("correlates")].sum() == min(200, adata.n_vars)

        assert Key.cytotrace("score") in adata.obs
        assert Key.cytotrace("pseudotime") in adata.obs
        assert Key.cytotrace("num_exp_genes") in adata.obs
        assert is_integer_dtype(adata.obs[Key.cytotrace("num_exp_genes")])
        np.testing.assert_array_equal(
            k.pseudotime, adata.obs[Key.cytotrace("pseudotime")].values
        )
        np.testing.assert_array_equal(k.pseudotime.min(), 0.0)
        np.testing.assert_array_equal(k.pseudotime.max(), 1.0)

    def test_compute_transition_matrix(self, adata: AnnData):
        k = CytoTRACEKernel(adata, use_raw=False, layer="X", aggregation="mean")
        k.compute_transition_matrix()

        np.testing.assert_allclose(k.transition_matrix.sum(1), 1.0)

    def test_inversion(self, adata: AnnData):
        k = ~CytoTRACEKernel(adata, use_raw=False, layer="X", aggregation="mean")

        pt = adata.obs[Key.cytotrace("pseudotime")].values
        np.testing.assert_array_equal(np.max(pt) - pt, k.pseudotime)


class TestSingleFlow:
    def test_no_transition_matrix(self, kernel: Kernel):
        kernel._transition_matrix = None
        with pytest.raises(RuntimeError, match=r"Compute transition matrix first as"):
            kernel.plot_single_flow("Astrocytes", "clusters", "age(days)")

    def test_invalid_cluster_key(self, kernel: Kernel):
        with pytest.raises(KeyError, match=r"Unable to find clusters in"):
            kernel.plot_single_flow("Astrocytes", "foo", "age(days)")

    def test_invalid_source_cluster(self, kernel: Kernel):
        with pytest.raises(ValueError, match=r"Invalid source cluster"):
            kernel.plot_single_flow("foo", "clusters", "age(days)")

    def test_too_few_invalid_clusters(self, kernel: Kernel):
        with pytest.raises(ValueError, match=r"Expected at least `2` clusters"):
            kernel.plot_single_flow(
                "Astrocytes", "clusters", "age(days)", clusters=["foo", "bar", "baz"]
            )

    def test_all_invalid_clusters(self, kernel: Kernel):
        with pytest.raises(ValueError, match=r"No valid clusters have been selected."):
            kernel.plot_single_flow(
                "quux", "clusters", "age(days)", clusters=["foo", "bar", "baz"]
            )

    def test_invalid_time_key(self, kernel: Kernel):
        with pytest.raises(
            KeyError, match=r"Unable to find data in `adata.obs\['foo'\]`."
        ):
            kernel.plot_single_flow("Astrocytes", "clusters", "foo")

    def test_too_few_valid_timepoints(self, kernel: Kernel):
        with pytest.raises(ValueError, match=r"Expected at least `2` time points"):
            kernel.plot_single_flow(
                "Astrocytes", "clusters", "age(days)", time_points=["35"]
            )

    def test_all_invalid_times(self, kernel: Kernel):
        with pytest.raises(ValueError, match=r"No valid time points"):
            kernel.plot_single_flow(
                "Astrocytes", "clusters", "age(days)", time_points=[0, 1, 2]
            )

    def test_time_key_cannot_be_coerced_to_numeric(self, kernel: Kernel):
        with pytest.raises(TypeError, match=r"Unable to convert .* to `float`."):
            kernel.plot_single_flow("Astrocytes", "clusters", "clusters")

    def test_remove_empty_clusters_none_remain(self, kernel: Kernel):
        with pytest.raises(ValueError, match=r"After removing clusters with no"):
            kernel.plot_single_flow(
                "Astrocytes",
                "clusters",
                "age(days)",
                min_flow=np.inf,
                remove_empty_clusters=True,
            )


class TestKernelIO:
    @pytest.mark.parametrize("copy", [False, True])
    @pytest.mark.parametrize("write_adata", [False, True])
    def test_read_write(self, kernel: Kernel, tmpdir, write_adata: bool, copy: bool):
        path = Path(tmpdir) / "kernel.pickle"
        kernel.write(path, write_adata=write_adata)

        if write_adata:
            k: Kernel = type(kernel).read(path)
            assert k.adata is not None
        else:
            with open(path, "rb") as fin:
                k: Kernel = pickle.load(fin)
                assert k.adata is None
                assert k.shape == (kernel.adata.n_obs, kernel.adata.n_obs)
            k: Kernel = type(kernel).read(path, adata=kernel.adata, copy=copy)
            if copy:
                assert k.adata is not kernel.adata
            else:
                assert k.adata is kernel.adata

        np.testing.assert_array_equal(k.transition_matrix.A, kernel.transition_matrix.A)
