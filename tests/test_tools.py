import pytest

from anndata import AnnData

import numpy as np
import pandas as pd

import cellrank as cr
from cellrank.tl._constants import AbsProbKey, TermStatesKey, _probs
from cellrank.tl.kernels._base_kernel import KernelAdd


class TestLineages:
    def test_no_root_cells(self, adata: AnnData):
        with pytest.raises(RuntimeError):
            cr.tl.lineages(adata)

    def test_normal_run(self, adata_cflare):
        cr.tl.lineages(adata_cflare)

        assert str(AbsProbKey.FORWARD) in adata_cflare.obsm.keys()

    def test_normal_run_copy(self, adata_cflare):
        adata_cr2 = cr.tl.lineages(adata_cflare, copy=True)

        assert isinstance(adata_cr2, AnnData)
        assert adata_cflare is not adata_cr2


class TestLineageDrivers:
    def test_no_abs_probs(self, adata: AnnData):
        with pytest.raises(RuntimeError):
            cr.tl.lineage_drivers(adata)

    def test_normal_run(self, adata_cflare):
        cr.tl.lineages(adata_cflare)
        cr.tl.lineage_drivers(adata_cflare, use_raw=False)

        for name in adata_cflare.obsm[AbsProbKey.FORWARD.s].names:
            assert np.all(adata_cflare.var[f"to {name} corr"] >= -1.0)
            assert np.all(adata_cflare.var[f"to {name} corr"] <= 1.0)

            assert np.all(adata_cflare.var[f"to {name} qval"] >= 0)
            assert np.all(adata_cflare.var[f"to {name} qval"] <= 1.0)

    def test_normal_run_raw(self, adata_cflare):
        adata_cflare.raw = adata_cflare.copy()
        cr.tl.lineages(adata_cflare)
        cr.tl.lineage_drivers(adata_cflare, use_raw=True)

        for name in adata_cflare.obsm[AbsProbKey.FORWARD.s].names:
            assert np.all(adata_cflare.raw.var[f"to {name} corr"] >= -1.0)
            assert np.all(adata_cflare.raw.var[f"to {name} corr"] <= 1.0)

            assert np.all(adata_cflare.raw.var[f"to {name} qval"] >= 0)
            assert np.all(adata_cflare.raw.var[f"to {name} qval"] <= 1.0)

    def test_return_drivers(self, adata_cflare):
        cr.tl.lineages(adata_cflare)
        res = cr.tl.lineage_drivers(adata_cflare, return_drivers=True, use_raw=False)

        assert isinstance(res, pd.DataFrame)
        assert res.shape[0] == adata_cflare.n_vars
        assert {
            "0 corr",
            "0 pval",
            "0 qval",
            "0 ci low",
            "0 ci high",
            "1 corr",
            "1 pval",
            "1 qval",
            "1 ci low",
            "1 ci high",
        } == set(res.columns)

        for name in ["0", "1"]:
            assert np.all(adata_cflare.var[f"to {name} corr"] >= -1.0)
            assert np.all(adata_cflare.var[f"to {name} corr"] <= 1.0)
            # res is sorted
            np.testing.assert_allclose(
                adata_cflare.var[f"to {name} corr"],
                res[f"{name} corr"].loc[adata_cflare.var_names],
            )

            assert np.all(res[f"{name} pval"] >= 0)
            assert np.all(res[f"{name} pval"] <= 1.0)

            assert np.all(adata_cflare.var[f"to {name} qval"] >= 0)
            assert np.all(adata_cflare.var[f"to {name} qval"] <= 1.0)
            # res is sorted
            np.testing.assert_allclose(
                adata_cflare.var[f"to {name} qval"],
                res[f"{name} qval"].loc[adata_cflare.var_names],
            )

            assert np.all(res[f"{name} corr"] >= res[f"{name} ci low"])
            assert np.all(res[f"{name} corr"] <= res[f"{name} ci high"])

    def test_invalid_mode(self, adata_cflare: AnnData):
        cr.tl.lineages(adata_cflare)
        with pytest.raises(ValueError):
            cr.tl.lineage_drivers(adata_cflare, use_raw=False, method="foobar")

    def test_invalid_n_perms_type(self, adata_cflare: AnnData):
        cr.tl.lineages(adata_cflare)
        with pytest.raises(TypeError):
            cr.tl.lineage_drivers(
                adata_cflare, use_raw=False, n_perms=None, method="perm_test"
            )

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
            return_drivers=True,
            n_perms=10,
            n_jobs=1,
            seed=0,
            method="perm_test",
        )
        res_b = cr.tl.lineage_drivers(
            adata_cflare,
            use_raw=False,
            return_drivers=True,
            n_perms=10,
            n_jobs=1,
            seed=0,
            method="perm_test",
        )
        res_diff_seed = cr.tl.lineage_drivers(
            adata_cflare,
            use_raw=False,
            return_drivers=True,
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
            return_drivers=True,
            n_perms=10,
            n_jobs=2,
            backend="threading",
            seed=42,
            method="perm_test",
        )
        res_b = cr.tl.lineage_drivers(
            adata_cflare,
            use_raw=False,
            return_drivers=True,
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
            return_drivers=True,
            confidence_level=0.95,
        )
        res_wide = cr.tl.lineage_drivers(
            adata_cflare,
            use_raw=False,
            return_drivers=True,
            confidence_level=0.99,
        )

        for name in ["0", "1"]:
            assert np.all(res_narrow[f"{name} ci low"] >= res_wide[f"{name} ci low"])
            assert np.all(res_narrow[f"{name} ci high"] <= res_wide[f"{name} ci high"])


class TestRootFinal:
    def test_find_root(self, adata: AnnData):
        cr.tl.initial_states(adata)

        assert str(TermStatesKey.BACKWARD) in adata.obs.keys()
        assert _probs(TermStatesKey.BACKWARD) in adata.obs.keys()

    def test_find_final(self, adata: AnnData):
        cr.tl.terminal_states(adata, n_states=5, fit_kwargs=dict(n_cells=5))

        assert str(TermStatesKey.FORWARD) in adata.obs.keys()
        assert _probs(TermStatesKey.FORWARD) in adata.obs.keys()

    def test_invalid_cluster_key(self, adata: AnnData):
        with pytest.raises(KeyError):
            cr.tl.initial_states(adata, cluster_key="foo")

    def test_invalid_weight(self, adata: AnnData):
        with pytest.raises(ValueError):
            cr.tl.initial_states(adata, weight_connectivities=10)

    def test_invalid_invalid_mode(self, adata: AnnData):
        with pytest.raises(ValueError):
            cr.tl.initial_states(adata, mode="foobar")


class TestTransitionMatrix:
    def test_invalid_velocity_key(self, adata: AnnData):
        with pytest.raises(KeyError):
            cr.tl.transition_matrix(adata, vkey="foo")

    def test_invalid_weight(self, adata: AnnData):
        with pytest.raises(ValueError):
            cr.tl.transition_matrix(adata, weight_connectivities=-1)

    def test_forward(self, adata: AnnData):
        kernel_add = cr.tl.transition_matrix(adata, backward=False, softmax_scale=None)

        assert isinstance(kernel_add, KernelAdd)
        assert not kernel_add.backward

    def test_only_connectivities(self, adata: AnnData):
        ck = cr.tl.transition_matrix(adata, weight_connectivities=1)

        assert isinstance(ck, cr.tl.kernels.ConnectivityKernel)

    def test_only_velocity(self, adata: AnnData):
        vk = cr.tl.transition_matrix(adata, weight_connectivities=0)

        assert isinstance(vk, cr.tl.kernels.VelocityKernel)

    def test_backward(self, adata: AnnData):
        kernel_add = cr.tl.transition_matrix(adata, backward=True, softmax_scale=None)

        assert isinstance(kernel_add, KernelAdd)
        assert kernel_add.backward
