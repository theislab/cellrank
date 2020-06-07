# -*- coding: utf-8 -*-

from pandas.api.types import is_categorical_dtype

from anndata import AnnData

import cellrank as cr
from cellrank.tools.kernels import VelocityKernel, ConnectivityKernel
from cellrank.tools._constants import (
    LinKey,
    StateKey,
    Direction,
    _probs,
    _colors,
    _lin_names,
    _transition,
)


def _assert_has_all_keys(adata: AnnData, direction: Direction):
    assert _transition(direction) in adata.uns.keys()

    if direction == Direction.FORWARD:
        assert str(LinKey.FORWARD) in adata.obsm
        assert isinstance(adata.obsm[str(LinKey.FORWARD)], cr.tl.Lineage)

        assert _colors(LinKey.FORWARD) in adata.uns.keys()
        assert _lin_names(LinKey.FORWARD) in adata.uns.keys()

        assert str(StateKey.FORWARD) in adata.obs
        assert is_categorical_dtype(adata.obs[str(StateKey.FORWARD)])

        assert _probs(StateKey.FORWARD) in adata.obs
    else:
        assert str(LinKey.BACKWARD) in adata.obsm
        assert isinstance(adata.obsm[str(LinKey.BACKWARD)], cr.tl.Lineage)

        assert _colors(LinKey.BACKWARD) in adata.uns.keys()
        assert _lin_names(LinKey.BACKWARD) in adata.uns.keys()

        assert str(StateKey.BACKWARD) in adata.obs
        assert is_categorical_dtype(adata.obs[str(StateKey.BACKWARD)])

        assert _probs(StateKey.BACKWARD) in adata.obs


class TestHighLevelPipeline:
    def test_fwd_pipeline_cflare(self, adata):
        cr.tl.final_states(adata, estimator=cr.tl.CFLARE, cluster_key="clusters")
        cr.tl.lineages(adata)
        cr.pl.lineages(adata)

        _assert_has_all_keys(adata, Direction.FORWARD)

    def test_bwd_pipeline_cflare(self, adata):
        cr.tl.root_states(adata, estimator=cr.tl.CFLARE, cluster_key="clusters")
        cr.tl.lineages(adata, final=False)
        cr.pl.lineages(adata, final=False)

        _assert_has_all_keys(adata, Direction.BACKWARD)

    def test_fwd_pipeline_gpcca(self, adata):
        cr.tl.final_states(adata, estimator=cr.tl.GPCCA, cluster_key="clusters")
        cr.tl.lineages(adata)
        cr.pl.lineages(adata)

        _assert_has_all_keys(adata, Direction.FORWARD)

    def test_bwd_pipeline_gpcca(self, adata):
        cr.tl.root_states(adata, estimator=cr.tl.GPCCA, cluster_key="clusters")
        cr.tl.lineages(adata, final=False)
        cr.pl.lineages(adata, final=False)

        _assert_has_all_keys(adata, Direction.BACKWARD)


class TestLowLevelPipeline:
    def test_fwd_pipelne_cflare(self, adata):
        vk = VelocityKernel(adata).compute_transition_matrix()
        ck = ConnectivityKernel(adata).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc_fwd = cr.tl.CFLARE(final_kernel)

        mc_fwd.compute_partition()

        mc_fwd.compute_eig()
        mc_fwd.plot_spectrum()
        mc_fwd.plot_spectrum(real_only=True)
        mc_fwd.plot_eig_embedding()
        mc_fwd.plot_eig_embedding(left=False)

        mc_fwd.compute_metastable_states(use=1)
        mc_fwd.plot_metastable_states()

        mc_fwd.compute_lin_probs()
        mc_fwd.plot_lin_probs()

        mc_fwd.compute_lineage_drivers(cluster_key="clusters", use_raw=False)

        _assert_has_all_keys(adata, Direction.FORWARD)

    def test_bwd_pipelne_cflare(self, adata):
        vk = VelocityKernel(adata, backward=True).compute_transition_matrix()
        ck = ConnectivityKernel(adata, backward=True).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc_bwd = cr.tl.CFLARE(final_kernel)

        mc_bwd.compute_partition()

        mc_bwd.compute_eig()
        mc_bwd.plot_spectrum()
        mc_bwd.plot_spectrum(real_only=True)
        mc_bwd.plot_eig_embedding()
        mc_bwd.plot_eig_embedding(left=False)

        mc_bwd.compute_metastable_states(use=1)
        mc_bwd.plot_metastable_states()

        mc_bwd.compute_lin_probs()
        mc_bwd.plot_lin_probs()

        mc_bwd.compute_lineage_drivers(cluster_key="clusters", use_raw=False)

        _assert_has_all_keys(adata, Direction.BACKWARD)

    def test_fwd_pipelne_gpcca(self, adata):
        vk = VelocityKernel(adata).compute_transition_matrix()
        ck = ConnectivityKernel(adata).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc_fwd = cr.tl.GPCCA(final_kernel)

        mc_fwd.compute_partition()

        mc_fwd.compute_eig()
        mc_fwd.plot_spectrum()
        mc_fwd.plot_spectrum(real_only=True)

        mc_fwd.compute_schur(5)
        mc_fwd.plot_schur_embedding()

        mc_fwd.compute_metastable_states(2)
        mc_fwd.plot_metastable_states()
        mc_fwd.plot_coarse_T(show_initial_dist=True, show_stationary_dist=True)
        mc_fwd.plot_schur_matrix()

        mc_fwd.set_main_states()
        mc_fwd.plot_main_states()

        mc_fwd.compute_lineage_drivers(cluster_key="clusters", use_raw=False)

        _assert_has_all_keys(adata, Direction.FORWARD)

    def test_bwd_pipelne_gpcca(self, adata):
        vk = VelocityKernel(adata, backward=True).compute_transition_matrix()
        ck = ConnectivityKernel(adata, backward=True).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc_bwd = cr.tl.GPCCA(final_kernel)

        mc_bwd.compute_eig()
        mc_bwd.plot_spectrum()
        mc_bwd.plot_spectrum(real_only=True)

        mc_bwd.compute_schur(5)
        mc_bwd.plot_schur_embedding()

        mc_bwd.compute_metastable_states(2)
        mc_bwd.plot_metastable_states()
        mc_bwd.plot_coarse_T(show_initial_dist=True, show_stationary_dist=True)
        mc_bwd.plot_schur_matrix()

        mc_bwd.set_main_states()
        mc_bwd.plot_main_states()

        mc_bwd.compute_lineage_drivers(cluster_key="clusters", use_raw=False)

        _assert_has_all_keys(adata, Direction.BACKWARD)
