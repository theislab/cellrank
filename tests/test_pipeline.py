# -*- coding: utf-8 -*-
import pytest

from anndata import AnnData

import numpy as np
from pandas.api.types import is_categorical_dtype

import cellrank as cr
from cellrank.tl.kernels import VelocityKernel, ConnectivityKernel
from cellrank.tl._constants import (
    Direction,
    DirPrefix,
    AbsProbKey,
    FinalStatesKey,
    _probs,
    _colors,
    _lin_names,
    _transition,
)


def _assert_has_all_keys(adata: AnnData, direction: Direction):
    assert _transition(direction) in adata.obsp.keys()
    assert f"{_transition(direction)}_params" in adata.uns.keys()

    if direction == Direction.FORWARD:
        assert str(AbsProbKey.FORWARD) in adata.obsm
        assert isinstance(adata.obsm[str(AbsProbKey.FORWARD)], cr.tl.Lineage)

        assert _colors(AbsProbKey.FORWARD) in adata.uns.keys()
        assert _lin_names(AbsProbKey.FORWARD) in adata.uns.keys()

        assert str(FinalStatesKey.FORWARD) in adata.obs
        assert is_categorical_dtype(adata.obs[str(FinalStatesKey.FORWARD)])

        assert _probs(FinalStatesKey.FORWARD) in adata.obs

        # check the correlations with all lineages have been computed
        lin_probs = adata.obsm[str(AbsProbKey.FORWARD)]
        np.in1d(
            [f"{str(DirPrefix.FORWARD)} {key}" for key in lin_probs.names],
            adata.var.keys(),
        ).all()

    else:
        assert str(AbsProbKey.BACKWARD) in adata.obsm
        assert isinstance(adata.obsm[str(AbsProbKey.BACKWARD)], cr.tl.Lineage)

        assert _colors(AbsProbKey.BACKWARD) in adata.uns.keys()
        assert _lin_names(AbsProbKey.BACKWARD) in adata.uns.keys()

        assert str(FinalStatesKey.BACKWARD) in adata.obs
        assert is_categorical_dtype(adata.obs[str(FinalStatesKey.BACKWARD)])

        assert _probs(FinalStatesKey.BACKWARD) in adata.obs

        # check the correlations with all lineages have been computed
        lin_probs = adata.obsm[str(AbsProbKey.BACKWARD)]
        np.in1d(
            [f"{str(DirPrefix.BACKWARD)} {key}" for key in lin_probs.names],
            adata.var.keys(),
        ).all()


class TestHighLevelPipeline:
    def test_fwd_pipeline_cflare(self, adata: AnnData):
        cr.tl.final_states(
            adata,
            estimator=cr.tl.estimators.CFLARE,
            cluster_key="clusters",
            method="kmeans",
            show_plots=True,
        )
        cr.tl.lineages(adata)
        cr.pl.lineages(adata)
        cr.tl.lineage_drivers(adata, use_raw=False)

        ln = adata.obsm[str(AbsProbKey.FORWARD)].names[0]
        cr.pl.lineage_drivers(adata, ln, use_raw=False, backward=False)

        _assert_has_all_keys(adata, Direction.FORWARD)

    def test_fwd_pipeline_invalid_raw_requested(self, adata: AnnData):
        cr.tl.final_states(
            adata,
            estimator=cr.tl.estimators.CFLARE,
            cluster_key="clusters",
            method="kmeans",
            show_plots=True,
        )
        cr.tl.lineages(adata)
        cr.pl.lineages(adata)
        cr.tl.lineage_drivers(adata, use_raw=False)

        ln = adata.obsm[str(AbsProbKey.FORWARD)].names[0]
        with pytest.raises(RuntimeError):
            cr.pl.lineage_drivers(adata, ln, use_raw=True, backward=False)

    def test_bwd_pipeline_cflare(self, adata: AnnData):
        cr.tl.root_states(
            adata,
            estimator=cr.tl.estimators.CFLARE,
            cluster_key="clusters",
            method="louvain",
            show_plots=True,
        )
        cr.tl.lineages(adata, backward=True)
        cr.pl.lineages(adata, backward=True)
        cr.tl.lineage_drivers(adata, use_raw=False, backward=True)

        ln = adata.obsm[str(AbsProbKey.BACKWARD)].names[0]
        cr.pl.lineage_drivers(adata, ln, use_raw=False, backward=True)

        _assert_has_all_keys(adata, Direction.BACKWARD)

    def test_fwd_pipeline_gpcca(self, adata: AnnData):
        cr.tl.final_states(
            adata,
            estimator=cr.tl.estimators.GPCCA,
            cluster_key="clusters",
            method="brandts",
            show_plots=True,
        )
        cr.tl.lineages(adata)
        cr.pl.lineages(adata)
        cr.tl.lineage_drivers(adata, use_raw=False)
        ln = adata.obsm[str(AbsProbKey.FORWARD)].names[0]
        cr.pl.lineage_drivers(adata, ln, use_raw=False, backward=False)

        _assert_has_all_keys(adata, Direction.FORWARD)

    def test_fwd_pipeline_gpcca_invalid_raw_requested(self, adata: AnnData):
        cr.tl.final_states(
            adata,
            estimator=cr.tl.estimators.GPCCA,
            cluster_key="clusters",
            method="brandts",
            show_plots=True,
        )
        cr.tl.lineages(adata)
        cr.pl.lineages(adata)
        cr.tl.lineage_drivers(adata, use_raw=False)
        ln = adata.obsm[str(AbsProbKey.FORWARD)].names[0]
        with pytest.raises(RuntimeError):
            cr.pl.lineage_drivers(adata, ln, use_raw=True, backward=False)

    def test_bwd_pipeline_gpcca(self, adata: AnnData):
        cr.tl.root_states(
            adata,
            estimator=cr.tl.estimators.GPCCA,
            cluster_key="clusters",
            method="brandts",
            show_plots=True,
        )
        cr.tl.lineages(adata, backward=True)
        cr.pl.lineages(adata, backward=True)
        cr.tl.lineage_drivers(adata, use_raw=False, backward=True)
        ln = adata.obsm[str(AbsProbKey.BACKWARD)].names[0]
        cr.pl.lineage_drivers(adata, ln, use_raw=False, backward=True)

        _assert_has_all_keys(adata, Direction.BACKWARD)


class TestLowLevelPipeline:
    def test_fwd_pipelne_cflare(self, adata: AnnData):
        vk = VelocityKernel(adata).compute_transition_matrix()
        ck = ConnectivityKernel(adata).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        estimator_fwd = cr.tl.estimators.CFLARE(final_kernel)

        estimator_fwd.compute_partition()

        estimator_fwd.compute_eigendecomposition()
        estimator_fwd.plot_spectrum()
        estimator_fwd.plot_spectrum(real_only=True)
        estimator_fwd.plot_eigendecomposition()
        estimator_fwd.plot_eigendecomposition(left=False)

        estimator_fwd.compute_final_states(use=1)
        estimator_fwd.plot_final_states()

        estimator_fwd.compute_absorption_probabilities()
        estimator_fwd.plot_absorption_probabilities()

        estimator_fwd.compute_lineage_drivers(cluster_key="clusters", use_raw=False)

        _assert_has_all_keys(adata, Direction.FORWARD)

    def test_bwd_pipelne_cflare(self, adata: AnnData):
        vk = VelocityKernel(adata, backward=True).compute_transition_matrix()
        ck = ConnectivityKernel(adata, backward=True).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        estimator_bwd = cr.tl.estimators.CFLARE(final_kernel)

        estimator_bwd.compute_partition()

        estimator_bwd.compute_eigendecomposition()
        estimator_bwd.plot_spectrum()
        estimator_bwd.plot_spectrum(real_only=True)
        estimator_bwd.plot_eigendecomposition()
        estimator_bwd.plot_eigendecomposition(left=False)

        estimator_bwd.compute_final_states(use=1)
        estimator_bwd.plot_final_states()

        estimator_bwd.compute_absorption_probabilities()
        estimator_bwd.plot_absorption_probabilities()

        estimator_bwd.compute_lineage_drivers(cluster_key="clusters", use_raw=False)

        _assert_has_all_keys(adata, Direction.BACKWARD)

    def test_fwd_pipelne_gpcca(self, adata: AnnData):
        vk = VelocityKernel(adata).compute_transition_matrix()
        ck = ConnectivityKernel(adata).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        estimator_fwd = cr.tl.estimators.GPCCA(final_kernel)

        estimator_fwd.compute_partition()

        estimator_fwd.compute_eigendecomposition()
        estimator_fwd.plot_spectrum()
        estimator_fwd.plot_spectrum(real_only=True)

        estimator_fwd.compute_schur(5, method="brandts")
        estimator_fwd.plot_schur()

        estimator_fwd.compute_metastable_states(3, n_cells=10)
        estimator_fwd.plot_metastable_states()
        estimator_fwd.plot_coarse_T(show_initial_dist=True, show_stationary_dist=True)
        estimator_fwd.plot_schur_matrix()

        # select all states
        estimator_fwd.set_final_states_from_metastable_states(n_cells=10)
        estimator_fwd.plot_final_states()

        estimator_fwd.compute_absorption_probabilities()
        estimator_fwd.compute_lineage_drivers(cluster_key="clusters", use_raw=False)

        _assert_has_all_keys(adata, Direction.FORWARD)

        # select a subset of states
        estimator_fwd.set_final_states_from_metastable_states(
            n_cells=16, names=estimator_fwd.metastable_states.cat.categories[:2],
        )
        estimator_fwd.plot_final_states()

        estimator_fwd.compute_absorption_probabilities()
        estimator_fwd.compute_lineage_drivers(cluster_key="clusters", use_raw=False)

        _assert_has_all_keys(adata, Direction.FORWARD)

    def test_bwd_pipelne_gpcca(self, adata: AnnData):
        vk = VelocityKernel(adata, backward=True).compute_transition_matrix()
        ck = ConnectivityKernel(adata, backward=True).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        estimator_bwd = cr.tl.estimators.GPCCA(final_kernel)

        estimator_bwd.compute_eigendecomposition()
        estimator_bwd.plot_spectrum()
        estimator_bwd.plot_spectrum(real_only=True)

        estimator_bwd.compute_schur(5, method="brandts")
        estimator_bwd.plot_schur()

        estimator_bwd.compute_metastable_states(3, n_cells=16)
        estimator_bwd.plot_metastable_states()
        estimator_bwd.plot_coarse_T(show_initial_dist=True, show_stationary_dist=True)
        estimator_bwd.plot_schur_matrix()

        # select all cells
        estimator_bwd.set_final_states_from_metastable_states(n_cells=16)
        estimator_bwd.plot_final_states()

        estimator_bwd.compute_absorption_probabilities()
        estimator_bwd.compute_lineage_drivers(cluster_key="clusters", use_raw=False)

        _assert_has_all_keys(adata, Direction.BACKWARD)

        # select a subset of states
        estimator_bwd.set_final_states_from_metastable_states(
            n_cells=16, names=estimator_bwd.metastable_states.cat.categories[:2],
        )
        estimator_bwd.plot_final_states()

        estimator_bwd.compute_absorption_probabilities()
        estimator_bwd.compute_lineage_drivers(cluster_key="clusters", use_raw=False)

        _assert_has_all_keys(adata, Direction.BACKWARD)
