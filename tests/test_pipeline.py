# -*- coding: utf-8 -*-

import cellrank as cr

from cellrank.tools._constants import (
    Direction,
    RcKey,
    LinKey,
    _transition,
    _colors,
    _lin_names,
    _probs,
)
from cellrank.tools.kernels import VelocityKernel, ConnectivityKernel
from anndata import AnnData

from pandas.api.types import is_categorical_dtype


def _assert_has_all_keys(adata: AnnData, direction: Direction):
    assert _transition(direction) in adata.uns.keys()

    if direction == Direction.FORWARD:
        assert str(LinKey.FORWARD) in adata.obsm
        assert isinstance(adata.obsm[str(LinKey.FORWARD)], cr.tl.Lineage)

        assert _colors(LinKey.FORWARD) in adata.uns.keys()
        assert _lin_names(LinKey.FORWARD) in adata.uns.keys()

        assert str(RcKey.FORWARD) in adata.obs
        assert is_categorical_dtype(adata.obs[str(RcKey.FORWARD)])

        assert _probs(RcKey.FORWARD) in adata.obs
    else:
        assert str(LinKey.BACKWARD) in adata.obsm
        assert isinstance(adata.obsm[str(LinKey.BACKWARD)], cr.tl.Lineage)

        assert _colors(LinKey.BACKWARD) in adata.uns.keys()
        assert _lin_names(LinKey.BACKWARD) in adata.uns.keys()

        assert str(RcKey.BACKWARD) in adata.obs
        assert is_categorical_dtype(adata.obs[str(RcKey.BACKWARD)])

        assert _probs(RcKey.BACKWARD) in adata.obs


class TestHighLevelPipeline:
    def test_fwd_pipeline(self, adata):
        cr.tl.find_final(adata, cluster_key="clusters")
        cr.tl.lineages(adata)
        cr.pl.lineages(adata)

        _assert_has_all_keys(adata, Direction.FORWARD)

    def test_bwd_pipeline(self, adata):
        cr.tl.find_root(adata, cluster_key="clusters")
        cr.tl.lineages(adata, final=False)
        cr.pl.lineages(adata, final=False)

        _assert_has_all_keys(adata, Direction.BACKWARD)


class TestLowLevelPipeline:
    def test_fwd_pipelne(self, adata):
        vk = VelocityKernel(adata).compute_transition_matrix()
        ck = ConnectivityKernel(adata).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc_fwd = cr.tl.CFLARE(final_kernel)

        mc_fwd.compute_partition()

        mc_fwd.compute_eig()
        mc_fwd.plot_spectrum()
        mc_fwd.plot_real_spectrum()
        mc_fwd.plot_eig_embedding()
        mc_fwd.plot_eig_embedding(left=False)

        mc_fwd.compute_approx_rcs(use=1)
        mc_fwd.plot_approx_rcs()

        mc_fwd.compute_lin_probs()
        mc_fwd.plot_lin_probs()

        mc_fwd.compute_lineage_drivers(cluster_key="clusters", use_raw=False)

        _assert_has_all_keys(adata, Direction.FORWARD)

    def test_bwd_pipelne(self, adata):
        vk = VelocityKernel(adata, backward=True).compute_transition_matrix()
        ck = ConnectivityKernel(adata, backward=True).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc_bwd = cr.tl.CFLARE(final_kernel)

        mc_bwd.compute_partition()

        mc_bwd.compute_eig()
        mc_bwd.plot_spectrum()
        mc_bwd.plot_real_spectrum()
        mc_bwd.plot_eig_embedding()
        mc_bwd.plot_eig_embedding(left=False)

        mc_bwd.compute_approx_rcs(use=1)
        mc_bwd.plot_approx_rcs()

        mc_bwd.compute_lin_probs()
        mc_bwd.plot_lin_probs()

        mc_bwd.compute_lineage_drivers(cluster_key="clusters", use_raw=False)

        _assert_has_all_keys(adata, Direction.BACKWARD)
