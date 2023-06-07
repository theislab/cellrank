import numpy as np
import pandas as pd
from pandas.api.types import is_categorical_dtype

from anndata import AnnData

import cellrank as cr
from cellrank._utils._key import Key
from cellrank.kernels import ConnectivityKernel, VelocityKernel


def _assert_has_all_keys(adata: AnnData, bwd: bool = False) -> None:
    # fmt: off
    # term states
    key = Key.obs.term_states(bwd)
    assert is_categorical_dtype(adata.obs[key])
    assert Key.obs.probs(key) in adata.obs
    assert Key.uns.colors(key) in adata.uns

    # fate probabilities
    fate_probs = adata.obsm[Key.obsm.fate_probs(bwd)]
    assert isinstance(fate_probs, cr.Lineage)
    np.testing.assert_array_equal(fate_probs.names, adata.obs[key].cat.categories)
    np.testing.assert_array_equal(fate_probs.colors, adata.uns[Key.uns.colors(key)])
    np.testing.assert_allclose(fate_probs.X.sum(1), 1.0, rtol=1e-3)

    # drivers
    assert isinstance(adata.varm[Key.varm.lineage_drivers(bwd)], pd.DataFrame)
    # fmt: on


class TestLowLevelPipeline:
    def test_fwd_pipeline_cflare(self, adata: AnnData):
        vk = VelocityKernel(adata).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        estimator_fwd = cr.estimators.CFLARE(final_kernel)

        estimator_fwd.compute_eigendecomposition()
        estimator_fwd.plot_spectrum()
        estimator_fwd.plot_spectrum(real_only=True)

        estimator_fwd.predict(use=1, method="leiden")
        estimator_fwd.plot_macrostates(which="terminal")

        estimator_fwd.compute_fate_probabilities()
        estimator_fwd.plot_fate_probabilities()

        estimator_fwd.compute_lineage_drivers(cluster_key="clusters", use_raw=False)

        _assert_has_all_keys(adata)

    def test_bwd_pipeline_cflare(self, adata: AnnData):
        vk = VelocityKernel(adata, backward=True).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        estimator_bwd = cr.estimators.CFLARE(final_kernel)

        estimator_bwd.compute_eigendecomposition()
        estimator_bwd.plot_spectrum()
        estimator_bwd.plot_spectrum(real_only=True)

        estimator_bwd.predict(use=1, method="kmeans")
        estimator_bwd.plot_macrostates(which="terminal")

        estimator_bwd.compute_fate_probabilities()
        estimator_bwd.plot_fate_probabilities()

        estimator_bwd.compute_lineage_drivers(cluster_key="clusters", use_raw=False)

        _assert_has_all_keys(adata, bwd=True)

    def test_fwd_pipeline_gpcca(self, adata: AnnData):
        vk = VelocityKernel(adata).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        estimator_fwd = cr.estimators.GPCCA(final_kernel)

        estimator_fwd.compute_eigendecomposition()
        estimator_fwd.plot_spectrum()
        estimator_fwd.plot_spectrum(real_only=True)

        estimator_fwd.compute_schur(5, method="brandts")

        estimator_fwd.compute_macrostates(3, n_cells=10)
        estimator_fwd.plot_macrostates(which="all")
        estimator_fwd.plot_coarse_T(show_initial_dist=True, show_stationary_dist=True)
        estimator_fwd.plot_schur_matrix()

        # select all states
        estimator_fwd.set_terminal_states(n_cells=10)
        estimator_fwd.plot_macrostates(which="terminal")

        estimator_fwd.compute_fate_probabilities()
        estimator_fwd.compute_lineage_drivers(cluster_key="clusters", use_raw=False)

        _assert_has_all_keys(adata)

        # select a subset of states
        estimator_fwd.set_terminal_states(
            n_cells=16,
            states=estimator_fwd.macrostates.cat.categories[:2],
        )
        estimator_fwd.plot_macrostates(which="terminal")

        estimator_fwd.compute_fate_probabilities()
        estimator_fwd.compute_lineage_drivers(cluster_key="clusters", use_raw=False)

        _assert_has_all_keys(adata)

    def test_bwd_pipeline_gpcca(self, adata: AnnData):
        vk = VelocityKernel(adata, backward=True).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        estimator_bwd = cr.estimators.GPCCA(final_kernel)

        estimator_bwd.compute_eigendecomposition()
        estimator_bwd.plot_spectrum()
        estimator_bwd.plot_spectrum(real_only=True)

        estimator_bwd.compute_schur(5, method="brandts")

        estimator_bwd.compute_macrostates(3, n_cells=16)
        estimator_bwd.plot_macrostates(which="all")
        estimator_bwd.plot_coarse_T(show_initial_dist=True, show_stationary_dist=True)
        estimator_bwd.plot_schur_matrix()

        # select all cells
        estimator_bwd.set_terminal_states(n_cells=16)
        estimator_bwd.plot_macrostates(which="terminal")

        estimator_bwd.compute_fate_probabilities()
        estimator_bwd.compute_lineage_drivers(cluster_key="clusters", use_raw=False)

        _assert_has_all_keys(adata, bwd=True)

        # select a subset of states
        estimator_bwd.set_terminal_states(
            n_cells=16,
            states=estimator_bwd.macrostates.cat.categories[:2],
        )
        estimator_bwd.plot_macrostates(which="terminal")

        estimator_bwd.compute_fate_probabilities()
        estimator_bwd.compute_lineage_drivers(cluster_key="clusters", use_raw=False)

        _assert_has_all_keys(adata, bwd=True)
