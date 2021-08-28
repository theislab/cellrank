import pytest

import cellrank as cr
from anndata import AnnData
from cellrank.tl import Lineage
from cellrank._key import Key
from cellrank.tl.kernels._base_kernel import KernelAdd

import numpy as np
import pandas as pd


class TestLineages:
    def test_no_root_cells(self, adata: AnnData):
        with pytest.raises(RuntimeError):
            cr.tl.lineages(adata)

    def test_normal_run(self, adata_cflare):
        cr.tl.lineages(adata_cflare)

        key = Key.obsm.abs_probs(False)
        assert isinstance(adata_cflare.obsm[key], Lineage)

    def test_normal_run_copy(self, adata_cflare):
        adata_cr2 = cr.tl.lineages(adata_cflare, copy=True)

        assert isinstance(adata_cr2, AnnData)
        assert adata_cflare is not adata_cr2


class TestLineageDrivers:
    def test_no_abs_probs(self, adata: AnnData):
        with pytest.raises(KeyError):
            cr.tl.lineage_drivers(adata)

    @pytest.mark.parametrize("use_raw", [False, True])
    def test_normal_run(self, adata_cflare: AnnData, use_raw: bool):
        if use_raw:
            adata_cflare.raw = adata_cflare.copy()
        cr.tl.lineages(adata_cflare)
        cr.tl.lineage_drivers(adata_cflare, use_raw=use_raw)

        bwd = False
        key = Key.varm.lineage_drivers(bwd)
        names = adata_cflare.obsm[Key.obsm.abs_probs(bwd)].names

        if use_raw:
            adata_cflare = adata_cflare.raw
        assert isinstance(adata_cflare.varm[key], pd.DataFrame)
        for name in names:
            assert np.all(adata_cflare.varm[key][f"{name}_corr"] >= -1.0)
            assert np.all(adata_cflare.varm[key][f"{name}_corr"] <= 1.0)

            assert np.all(adata_cflare.varm[key][f"{name}_qval"] >= 0)
            assert np.all(adata_cflare.varm[key][f"{name}_qval"] <= 1.0)

    def test_invalid_mode(self, adata_cflare: AnnData):
        cr.tl.lineages(adata_cflare)
        with pytest.raises(ValueError):
            cr.tl.lineage_drivers(adata_cflare, use_raw=False, method="foobar")

    def test_invalid_n_perms_value(self, adata_cflare: AnnData):
        cr.tl.lineages(adata_cflare)
        with pytest.raises(ValueError):
            cr.tl.lineage_drivers(
                adata_cflare, use_raw=False, n_perms=0, method="perm_test"
            )

    def test_seed_reproducible(self, adata_cflare: AnnData):
        cr.tl.lineages(adata_cflare)
        res_a = cr.tl.lineage_drivers(
            adata_cflare,
            use_raw=False,
            n_perms=10,
            n_jobs=1,
            seed=0,
            method="perm_test",
        )
        res_b = cr.tl.lineage_drivers(
            adata_cflare,
            use_raw=False,
            n_perms=10,
            n_jobs=1,
            seed=0,
            method="perm_test",
        )
        res_diff_seed = cr.tl.lineage_drivers(
            adata_cflare,
            use_raw=False,
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

    def test_seed_reproducible_parallel(self, adata_cflare: AnnData):
        cr.tl.lineages(adata_cflare)
        res_a = cr.tl.lineage_drivers(
            adata_cflare,
            use_raw=False,
            n_perms=10,
            n_jobs=2,
            backend="threading",
            seed=42,
            method="perm_test",
        )
        res_b = cr.tl.lineage_drivers(
            adata_cflare,
            use_raw=False,
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

    def test_confidence_level(self, adata_cflare: AnnData):
        cr.tl.lineages(adata_cflare)
        res_narrow = cr.tl.lineage_drivers(
            adata_cflare,
            use_raw=False,
            confidence_level=0.95,
        )
        res_wide = cr.tl.lineage_drivers(
            adata_cflare,
            use_raw=False,
            confidence_level=0.99,
        )

        for name in ["0", "1"]:
            assert np.all(res_narrow[f"{name}_ci_low"] >= res_wide[f"{name}_ci_low"])
            assert np.all(res_narrow[f"{name}_ci_high"] <= res_wide[f"{name}_ci_high"])


class TestRootFinal:
    def test_find_root(self, adata: AnnData):
        cr.tl.initial_states(adata)

        key = Key.obs.term_states(True)
        assert key in adata.obs
        assert Key.obs.probs(key) in adata.obs

    def test_find_final(self, adata: AnnData):
        cr.tl.terminal_states(adata, n_states=5, fit_kwargs=dict(n_cells=5))

        key = Key.obs.term_states(False)
        assert key in adata.obs
        assert Key.obs.probs(key) in adata.obs

    def test_invalid_cluster_key(self, adata: AnnData):
        with pytest.raises(KeyError):
            cr.tl.initial_states(adata, cluster_key="foo")

    def test_invalid_weight(self, adata: AnnData):
        with pytest.raises(ValueError):
            cr.tl.initial_states(adata, weight_connectivities=10)

    def test_invalid_invalid_mode(self, adata: AnnData):
        with pytest.raises(ValueError):
            cr.tl.initial_states(adata, mode="foobar")

    @pytest.mark.parametrize("force_recompute", [False, True])
    def test_force_recompute(self, adata: AnnData, force_recompute: bool):
        cr.tl.terminal_states(adata, n_states=5, fit_kwargs=dict(n_cells=5))
        tmat1 = adata.obsp["T_fwd"]

        cr.tl.terminal_states(
            adata,
            n_states=5,
            fit_kwargs=dict(n_cells=5),
            force_recompute=force_recompute,
        )
        tmat2 = adata.obsp["T_fwd"]
        if force_recompute:
            assert tmat1 is not tmat2
            np.testing.assert_allclose(tmat1.A, tmat2.A)
            np.testing.assert_allclose(tmat1.sum(1), 1.0)
        else:
            assert tmat1 is tmat2


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

        np.testing.assert_array_equal(ck.transition_matrix, adata.obsp[key])

    def test_only_velocity(self, adata: AnnData):
        vk = cr.tl.transition_matrix(adata, weight_connectivities=0)

        assert isinstance(vk, cr.tl.kernels.VelocityKernel)

    def test_backward(self, adata: AnnData):
        kernel_add = cr.tl.transition_matrix(adata, backward=True, softmax_scale=None)

        assert isinstance(kernel_add, KernelAdd)
        assert kernel_add.backward
