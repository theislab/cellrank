import pytest

import cellrank as cr
from anndata import AnnData
from cellrank._key import Key
from cellrank.tl.kernels import VelocityKernel, ConnectivityKernel

import numpy as np
import pandas as pd
from pandas.api.types import is_categorical_dtype


def _assert_has_all_keys(adata: AnnData, bwd: bool = False) -> None:
    # fmt: off
    # term states
    key = Key.obs.term_states(bwd)
    assert is_categorical_dtype(adata.obs[key])
    assert Key.obs.probs(key) in adata.obs
    assert Key.uns.colors(key) in adata.uns

    # lineages
    abs_probs = adata.obsm[Key.obsm.abs_probs(bwd)]
    assert isinstance(abs_probs, cr.tl.Lineage)
    np.testing.assert_array_equal(abs_probs.names, adata.obs[key].cat.categories)
    np.testing.assert_array_equal(abs_probs.colors, adata.uns[Key.uns.colors(key)])
    np.testing.assert_allclose(abs_probs.X.sum(1), 1.0, rtol=1e-3)

    # drivers
    assert isinstance(adata.varm[Key.varm.lineage_drivers(bwd)], pd.DataFrame)
    # fmt: on


class TestHighLevelPipeline:
    def test_plot_states_not_computed(self, adata: AnnData):
        with pytest.raises(KeyError):
            cr.pl.initial_states(adata)
        with pytest.raises(KeyError):
            cr.pl.terminal_states(adata)

    def test_write_transition_matrix(self, adata: AnnData):
        cr.tl.transition_matrix(adata, key="foo")

        assert "foo" in adata.obsp
        np.testing.assert_allclose(adata.obsp["foo"].A.sum(1), 1.0)
        assert "foo_params" in adata.uns

    def test_states_no_precomputed_transition_matrix(self, adata: AnnData):
        cr.tl.terminal_states(adata, key="foo")

        np.testing.assert_allclose(adata.obsp[Key.uns.kernel(False)].A.sum(1), 1.0)

    def test_states_use_precomputed_transition_matrix(self, adata: AnnData):
        cr.tl.transition_matrix(adata, key="foo")
        obsp_keys = set(adata.obsp.keys())

        cr.tl.terminal_states(adata, key="foo")

        assert obsp_keys == set(adata.obsp.keys())

    def test_fwd_pipeline_cflare(self, adata: AnnData):
        cr.tl.terminal_states(
            adata,
            estimator=cr.tl.estimators.CFLARE,
            cluster_key="clusters",
            method="kmeans",
            show_plots=True,
        )
        cr.pl.terminal_states(adata)
        cr.tl.lineages(adata)
        cr.pl.lineages(adata)
        cr.tl.lineage_drivers(adata, use_raw=False)

        ln = adata.obsm[Key.obsm.abs_probs(False)].names[0]
        cr.pl.lineage_drivers(adata, ln, use_raw=False, backward=False)

        _assert_has_all_keys(adata)

    def test_fwd_pipeline_invalid_raw_requested(self, adata: AnnData):
        cr.tl.terminal_states(
            adata,
            estimator=cr.tl.estimators.CFLARE,
            cluster_key="clusters",
            method="kmeans",
            show_plots=True,
        )
        cr.pl.terminal_states(adata)
        cr.tl.lineages(adata)
        cr.pl.lineages(adata)
        cr.tl.lineage_drivers(adata, use_raw=False)

        ln = adata.obsm[Key.obsm.abs_probs(False)].names[0]
        with pytest.raises(RuntimeError):
            cr.pl.lineage_drivers(adata, ln, use_raw=True, backward=False)

    def test_bwd_pipeline_cflare(self, adata: AnnData):
        cr.tl.initial_states(
            adata,
            estimator=cr.tl.estimators.CFLARE,
            cluster_key="clusters",
            method="leiden",
            show_plots=True,
        )
        cr.pl.initial_states(adata)
        cr.tl.lineages(adata, backward=True)
        cr.pl.lineages(adata, backward=True)
        cr.tl.lineage_drivers(adata, use_raw=False, backward=True)

        ln = adata.obsm[Key.obsm.abs_probs(True)].names[0]
        cr.pl.lineage_drivers(adata, ln, use_raw=False, backward=True)

        _assert_has_all_keys(adata, bwd=True)

    def test_fwd_pipeline_gpcca(self, adata: AnnData):
        cr.tl.terminal_states(
            adata,
            estimator=cr.tl.estimators.GPCCA,
            cluster_key="clusters",
            method="brandts",
            show_plots=True,
        )
        cr.pl.terminal_states(adata)
        cr.tl.lineages(adata)
        cr.pl.lineages(adata)
        cr.tl.lineage_drivers(adata, use_raw=False)
        ln = adata.obsm[Key.obsm.abs_probs(False)].names[0]
        cr.pl.lineage_drivers(adata, ln, use_raw=False, backward=False)

        _assert_has_all_keys(adata)

    def test_fwd_pipeline_gpcca_invalid_raw_requested(self, adata: AnnData):
        cr.tl.terminal_states(
            adata,
            estimator=cr.tl.estimators.GPCCA,
            cluster_key="clusters",
            method="brandts",
            show_plots=True,
        )
        cr.pl.terminal_states(adata)
        cr.tl.lineages(adata)
        cr.pl.lineages(adata)
        cr.tl.lineage_drivers(adata, use_raw=False)
        ln = adata.obsm[Key.obsm.abs_probs(False)].names[0]
        with pytest.raises(RuntimeError):
            cr.pl.lineage_drivers(adata, ln, use_raw=True, backward=False)

    def test_bwd_pipeline_gpcca(self, adata: AnnData):
        cr.tl.initial_states(
            adata,
            estimator=cr.tl.estimators.GPCCA,
            cluster_key="clusters",
            method="brandts",
            show_plots=True,
        )
        cr.pl.initial_states(adata)
        cr.tl.lineages(adata, backward=True)
        cr.pl.lineages(adata, backward=True)
        cr.tl.lineage_drivers(adata, use_raw=False, backward=True)
        ln = adata.obsm[Key.obsm.abs_probs(True)].names[0]
        cr.pl.lineage_drivers(adata, ln, use_raw=False, backward=True)

        _assert_has_all_keys(adata, bwd=True)

    def test_multiple_read_write_diff(self, adata: AnnData):
        for n_states in range(2, 5):
            e = cr.tl.terminal_states(
                adata,
                estimator=cr.tl.estimators.GPCCA,
                cluster_key="clusters",
                method="brandts",
                show_plots=True,
                n_states=n_states,
                fit_kwargs={"n_cells": 5},
                return_estimator=True,
            )
            assert e.absorption_probabilities is None

            cr.tl.lineages(adata, backward=False)
            cr.pl.lineages(adata)
            cr.tl.lineage_drivers(adata)

            _assert_has_all_keys(adata)
            probs = adata.obsm[Key.obsm.abs_probs(False)]

            assert probs.shape[1] == n_states
            for name in probs.names:
                assert "Lineage" not in name, probs.names

    def test_multiple_read_write_same(self, adata: AnnData):
        e = cr.tl.terminal_states(
            adata,
            estimator=cr.tl.estimators.GPCCA,
            cluster_key="clusters",
            method="brandts",
            show_plots=True,
            n_states=2,
            fit_kwargs={"n_cells": 5},
            return_estimator=True,
        )
        cr.tl.lineages(adata, backward=False)
        cr.pl.lineages(adata)

        e.set_terminal_states({"foo": adata.obs_names[:20]})

        e = cr.tl.estimators.GPCCA(adata, obsp_key="T_fwd")
        assert e.macrostates is None
        assert e.absorption_probabilities is None


class TestLowLevelPipeline:
    def test_fwd_pipeline_cflare(self, adata: AnnData):
        vk = VelocityKernel(adata).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        estimator_fwd = cr.tl.estimators.CFLARE(final_kernel)

        estimator_fwd.compute_eigendecomposition()
        estimator_fwd.plot_spectrum()
        estimator_fwd.plot_spectrum(real_only=True)

        estimator_fwd.compute_terminal_states(use=1, method="leiden")
        estimator_fwd.plot_terminal_states()

        estimator_fwd.compute_absorption_probabilities()
        estimator_fwd.plot_absorption_probabilities()

        estimator_fwd.compute_lineage_drivers(cluster_key="clusters", use_raw=False)

        _assert_has_all_keys(adata)

    def test_bwd_pipeline_cflare(self, adata: AnnData):
        vk = VelocityKernel(adata, backward=True).compute_transition_matrix(
            softmax_scale=4
        )
        ck = ConnectivityKernel(adata, backward=True).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        estimator_bwd = cr.tl.estimators.CFLARE(final_kernel)

        estimator_bwd.compute_eigendecomposition()
        estimator_bwd.plot_spectrum()
        estimator_bwd.plot_spectrum(real_only=True)

        estimator_bwd.compute_terminal_states(use=1, method="kmeans")
        estimator_bwd.plot_terminal_states()

        estimator_bwd.compute_absorption_probabilities()
        estimator_bwd.plot_absorption_probabilities()

        estimator_bwd.compute_lineage_drivers(cluster_key="clusters", use_raw=False)

        _assert_has_all_keys(adata, bwd=True)

    def test_fwd_pipeline_gpcca(self, adata: AnnData):
        vk = VelocityKernel(adata).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        estimator_fwd = cr.tl.estimators.GPCCA(final_kernel)

        estimator_fwd.compute_eigendecomposition()
        estimator_fwd.plot_spectrum()
        estimator_fwd.plot_spectrum(real_only=True)

        estimator_fwd.compute_schur(5, method="brandts")

        estimator_fwd.compute_macrostates(3, n_cells=10)
        estimator_fwd.plot_macrostates()
        estimator_fwd.plot_coarse_T(show_initial_dist=True, show_stationary_dist=True)
        estimator_fwd.plot_schur_matrix()

        # select all states
        estimator_fwd.set_terminal_states_from_macrostates(n_cells=10)
        estimator_fwd.plot_terminal_states()

        estimator_fwd.compute_absorption_probabilities()
        estimator_fwd.compute_lineage_drivers(cluster_key="clusters", use_raw=False)

        _assert_has_all_keys(adata)

        # select a subset of states
        estimator_fwd.set_terminal_states_from_macrostates(
            n_cells=16,
            names=estimator_fwd.macrostates.cat.categories[:2],
        )
        estimator_fwd.plot_terminal_states()

        estimator_fwd.compute_absorption_probabilities()
        estimator_fwd.compute_lineage_drivers(cluster_key="clusters", use_raw=False)

        _assert_has_all_keys(adata)

    def test_bwd_pipeline_gpcca(self, adata: AnnData):
        vk = VelocityKernel(adata, backward=True).compute_transition_matrix(
            softmax_scale=4
        )
        ck = ConnectivityKernel(adata, backward=True).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        estimator_bwd = cr.tl.estimators.GPCCA(final_kernel)

        estimator_bwd.compute_eigendecomposition()
        estimator_bwd.plot_spectrum()
        estimator_bwd.plot_spectrum(real_only=True)

        estimator_bwd.compute_schur(5, method="brandts")

        estimator_bwd.compute_macrostates(3, n_cells=16)
        estimator_bwd.plot_macrostates()
        estimator_bwd.plot_coarse_T(show_initial_dist=True, show_stationary_dist=True)
        estimator_bwd.plot_schur_matrix()

        # select all cells
        estimator_bwd.set_terminal_states_from_macrostates(n_cells=16)
        estimator_bwd.plot_terminal_states()

        estimator_bwd.compute_absorption_probabilities()
        estimator_bwd.compute_lineage_drivers(cluster_key="clusters", use_raw=False)

        _assert_has_all_keys(adata, bwd=True)

        # select a subset of states
        estimator_bwd.set_terminal_states_from_macrostates(
            n_cells=16,
            names=estimator_bwd.macrostates.cat.categories[:2],
        )
        estimator_bwd.plot_terminal_states()

        estimator_bwd.compute_absorption_probabilities()
        estimator_bwd.compute_lineage_drivers(cluster_key="clusters", use_raw=False)

        _assert_has_all_keys(adata, bwd=True)
