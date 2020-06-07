# -*- coding: utf-8 -*-

from copy import deepcopy

import numpy as np
import pandas as pd
import pytest

from anndata import AnnData

import cellrank as cr
from cellrank.tools.kernels import VelocityKernel, ConnectivityKernel
from cellrank.tools._constants import StateKey, Direction


def _check_eigdecomposition(mc: cr.tl.GPCCA) -> None:
    assert isinstance(mc.eigendecomposition, dict)
    assert set(mc.eigendecomposition.keys()) == {
        "D",
        "eigengap",
        "params",
    }
    assert "V_l" not in mc.eigendecomposition
    assert "V_r" not in mc.eigendecomposition
    assert "stationary_dist" not in mc.eigendecomposition

    assert f"eig_{Direction.FORWARD}" in mc.adata.uns.keys()


def _check_compute_schur(mc: cr.tl.GPCCA) -> None:
    assert isinstance(mc.metastable_states, pd.Series)
    if "stationary_dist" in mc.eigendecomposition:
        assert mc._coarse_init_dist is None
        assert mc._schur_matrix is None
        assert mc.coarse_stationary_distribution is None
        assert mc.coarse_T is None
        assert mc.schur_vectors is None
    else:
        assert isinstance(mc._coarse_init_dist, pd.Series)
        assert isinstance(mc._schur_matrix, np.ndarray)
        assert mc.coarse_stationary_distribution is None or isinstance(
            mc.coarse_stationary_distribution, pd.Series
        )
        assert isinstance(mc.coarse_T, pd.DataFrame)
        assert isinstance(mc.schur_vectors, np.ndarray)


class TestCGPCCA:
    def test_compute_partition(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.GPCCA(final_kernel)
        mc.compute_partition()

        assert isinstance(mc.irreducible, bool)
        if not mc.irreducible:
            assert isinstance(mc.recurrent_classes, list)
            assert isinstance(mc.transient_classes, list)
            assert f"{StateKey.FORWARD}_rec_classes" in mc.adata.obs
            assert f"{StateKey.FORWARD}_trans_classes" in mc.adata.obs
        else:
            assert mc.recurrent_classes is None
            assert mc.transient_classes is None

    def test_compute_eig(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.GPCCA(final_kernel)
        mc.compute_eig(k=2)

        _check_eigdecomposition(mc)

        assert f"eig_{Direction.FORWARD}" in mc.adata.uns.keys()

    def test_compute_schur_invalid_n_comps(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.GPCCA(final_kernel)
        with pytest.raises(ValueError):
            mc.compute_schur(n_components=1)

    def test_compute_schur_invalid_method(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.GPCCA(final_kernel)
        with pytest.raises(ValueError):
            mc.compute_schur(method="foobar")

    def test_compute_schur_invalid_eig_sort(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.GPCCA(final_kernel)
        with pytest.raises(ValueError):
            mc.compute_schur(which="foobar")

    def test_compute_schur_write_eigvals(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="brandts")

        _check_eigdecomposition(mc)

    def test_compute_schur_write_eigvals_similar_to_orig_eigdecomp(
        self, adata_large: AnnData
    ):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.GPCCA(final_kernel)
        mc.compute_eig(k=11)

        _check_eigdecomposition(mc)
        orig_ed = deepcopy(mc.eigendecomposition)

        mc._eigendecomposition = None
        mc.compute_schur(n_components=10, method="brandts")

        _check_eigdecomposition(mc)
        schur_ed = mc.eigendecomposition

        assert orig_ed.keys() == schur_ed.keys()
        assert orig_ed["eigengap"] == schur_ed["eigengap"]
        assert orig_ed["params"] == schur_ed["params"]
        np.testing.assert_array_almost_equal(orig_ed["D"].real, schur_ed["D"].real)
        np.testing.assert_array_almost_equal(
            np.abs(orig_ed["D"].imag), np.abs(schur_ed["D"].imag)
        )  # complex conj.

    def test_compute_metastable_states_no_eig(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.GPCCA(final_kernel)
        with pytest.raises(RuntimeError):
            mc.compute_metastable_states(n_states=None)

    def test_compute_metastable_states_1_state_no_eig(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.GPCCA(final_kernel)
        mc.compute_metastable_states(n_states=1)

    def test_compute_metastable_none_states(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.GPCCA(final_kernel)
        mc.compute_eig()
        mc.compute_metastable_states(n_states=None)

        _check_compute_schur(mc)

    def test_compute_metastable_states_schur(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.GPCCA(final_kernel)
        mc.compute_schur(n_components=10)
        mc.compute_metastable_states(n_states=2)

        _check_compute_schur(mc)

    def test_compute_metastable_invalid_cluster_key(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.GPCCA(final_kernel)
        mc.compute_schur(n_components=10)
        with pytest.raises(KeyError):
            mc.compute_metastable_states(n_states=2, cluster_key="foobar")
