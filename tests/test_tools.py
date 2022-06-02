import pytest

import cellrank as cr
from anndata import AnnData
from cellrank.tl import Lineage
from cellrank._key import Key
from cellrank.tl.estimators import GPCCA
from cellrank.tl.kernels._base_kernel import KernelAdd

import numpy as np
import pandas as pd
from scipy.sparse import isspmatrix_csr


class TestLineageDrivers:
    @pytest.mark.parametrize("use_raw", [False, True])
    def test_normal_run(self, g: GPCCA, use_raw: bool):
        key = Key.varm.lineage_drivers(False)
        names = g.absorption_probabilities.names
        if use_raw:
            g.adata.raw = g.adata.copy()

        g.compute_lineage_drivers(use_raw=use_raw)

        adata = g.adata.raw if use_raw else g.adata

        assert isinstance(adata.varm[key], pd.DataFrame)
        for name in names:
            assert np.all(adata.varm[key][f"{name}_corr"] >= -1.0)
            assert np.all(adata.varm[key][f"{name}_corr"] <= 1.0)

            assert np.all(adata.varm[key][f"{name}_qval"] >= 0)
            assert np.all(adata.varm[key][f"{name}_qval"] <= 1.0)

    def test_invalid_method(self, g: GPCCA):
        with pytest.raises(ValueError, match=r".*foobar.*"):
            g.compute_lineage_drivers(method="foobar")

    def test_invalid_n_perms_value(self, g: GPCCA):
        with pytest.raises(ValueError, match=r".*n_perms.*"):
            g.compute_lineage_drivers(n_perms=0, method="perm_test")

    def test_seed_reproducible(self, g: GPCCA):
        res_a = g.compute_lineage_drivers(
            n_perms=10,
            n_jobs=1,
            seed=0,
            method="perm_test",
        )
        res_b = g.compute_lineage_drivers(
            n_perms=10,
            n_jobs=1,
            seed=0,
            method="perm_test",
        )
        res_diff_seed = g.compute_lineage_drivers(
            n_perms=10,
            n_jobs=1,
            seed=1,
            method="perm_test",
        )

        assert res_a is not res_b
        np.testing.assert_array_equal(res_a.index, res_b.index)
        np.testing.assert_array_equal(res_a.columns, res_b.columns)
        np.testing.assert_allclose(res_a.values, res_b.values)

        assert not np.allclose(res_a.values, res_diff_seed.values)

    def test_seed_reproducible_parallel(self, g: GPCCA):
        res_a = g.compute_lineage_drivers(
            n_perms=10,
            n_jobs=2,
            backend="threading",
            seed=42,
            method="perm_test",
        )
        res_b = g.compute_lineage_drivers(
            n_perms=10,
            n_jobs=2,
            backend="threading",
            seed=42,
            method="perm_test",
        )

        assert res_a is not res_b
        np.testing.assert_array_equal(res_a.index, res_b.index)
        np.testing.assert_array_equal(res_a.columns, res_b.columns)
        np.testing.assert_allclose(res_a.values, res_b.values)

    def test_confidence_level(self, g: GPCCA):
        res_narrow = g.compute_lineage_drivers(confidence_level=0.95)
        res_wide = g.compute_lineage_drivers(confidence_level=0.99)

        for name in ["0", "1"]:
            assert np.all(res_narrow[f"{name}_ci_low"] >= res_wide[f"{name}_ci_low"])
            assert np.all(res_narrow[f"{name}_ci_high"] <= res_wide[f"{name}_ci_high"])


class TestTransitionMatrix:
    def test_invalid_velocity_key(self, adata: AnnData):
        with pytest.raises(KeyError):
            cr.tl.transition_matrix(adata, vkey="foo")

    def test_invalid_weight(self, adata: AnnData):
        with pytest.raises(ValueError):
            cr.tl.transition_matrix(adata, weight_connectivities=-1)

    def test_invalid_conn_key(self, adata: AnnData):
        with pytest.raises(KeyError):
            cr.tl.transition_matrix(adata, conn_key="foo")

    def test_forward(self, adata: AnnData):
        kernel_add = cr.tl.transition_matrix(adata, backward=False, softmax_scale=None)

        assert isinstance(kernel_add, KernelAdd)
        assert not kernel_add.backward

    def test_only_connectivities(self, adata: AnnData):
        ck = cr.tl.transition_matrix(adata, weight_connectivities=1)

        assert isinstance(ck, cr.tl.kernels.ConnectivityKernel)

    def test_connectivities_conn_key(self, adata: AnnData):
        key = "foobar"
        assert key not in adata.obsp
        adata.obsp[key] = np.eye(adata.n_obs)

        ck = cr.tl.transition_matrix(adata, weight_connectivities=1, conn_key=key)

        assert isinstance(ck, cr.tl.kernels.ConnectivityKernel)
        assert isspmatrix_csr(ck.transition_matrix)

        np.testing.assert_array_equal(ck.transition_matrix.A, adata.obsp[key])

    def test_only_velocity(self, adata: AnnData):
        vk = cr.tl.transition_matrix(adata, weight_connectivities=0)

        assert isinstance(vk, cr.tl.kernels.VelocityKernel)

    def test_backward(self, adata: AnnData):
        kernel_add = cr.tl.transition_matrix(adata, backward=True, softmax_scale=None)

        assert isinstance(kernel_add, KernelAdd)
        assert kernel_add.backward
