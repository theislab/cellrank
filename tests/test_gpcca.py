# -*- coding: utf-8 -*-

from copy import deepcopy
from typing import Tuple

import numpy as np
import pandas as pd
import pytest

from anndata import AnnData

import cellrank as cr
from _helpers import assert_array_nan_equal
from cellrank.tools.kernels import VelocityKernel, ConnectivityKernel
from cellrank.tools._constants import (
    LinKey,
    Prefix,
    MetaKey,
    StateKey,
    Direction,
    _dp,
    _probs,
    _colors,
    _lin_names,
)


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


def _check_compute_meta(mc: cr.tl.GPCCA) -> None:
    assert mc.lineage_probabilities is None
    assert isinstance(mc._meta_lin_probs, cr.tl.Lineage)

    assert isinstance(mc.metastable_states, pd.Series)
    assert_array_nan_equal(mc.metastable_states, mc.adata.obs[str(MetaKey.FORWARD)])

    np.testing.assert_array_equal(mc._meta_states_colors, mc._meta_lin_probs.colors)
    np.testing.assert_array_equal(
        mc._meta_states_colors, mc.adata.uns[_colors(str(MetaKey.FORWARD))]
    )

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


def _check_main_states(mc: cr.tl.GPCCA, has_main_states: bool = True):
    if has_main_states:
        assert isinstance(mc.main_states, pd.Series)
        assert_array_nan_equal(mc.adata.obs[str(StateKey.FORWARD)], mc.main_states)
        np.testing.assert_array_equal(
            mc.adata.uns[_colors(StateKey.FORWARD)],
            mc.lineage_probabilities[list(mc.main_states.cat.categories)].colors,
        )

    assert isinstance(mc.diff_potential, np.ndarray)
    assert isinstance(mc.lineage_probabilities, cr.tl.Lineage)

    np.testing.assert_array_equal(
        mc.adata.obsm[str(LinKey.FORWARD)], mc.lineage_probabilities.X
    )
    np.testing.assert_array_equal(
        mc.adata.uns[_lin_names(LinKey.FORWARD)], mc.lineage_probabilities.names
    )
    np.testing.assert_array_equal(
        mc.adata.uns[_colors(LinKey.FORWARD)], mc.lineage_probabilities.colors
    )

    np.testing.assert_array_equal(mc.adata.obs[_dp(LinKey.FORWARD)], mc.diff_potential)
    np.testing.assert_array_equal(
        mc.adata.obs[_probs(StateKey.FORWARD)], mc.main_states_probabilities
    )


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
            mc.compute_schur(n_components=1, method="brandts")

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
            mc.compute_schur(which="foobar", method="brandts")

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
        mc.compute_eig(k=10)

        _check_eigdecomposition(mc)
        orig_ed = deepcopy(mc.eigendecomposition)

        mc._eigendecomposition = None
        mc.compute_schur(n_components=10, method="brandts")

        _check_eigdecomposition(mc)
        schur_ed = mc.eigendecomposition

        assert orig_ed.keys() == schur_ed.keys()
        assert orig_ed["eigengap"] == schur_ed["eigengap"]
        n = min(orig_ed["params"]["k"], schur_ed["params"]["k"])
        np.testing.assert_array_almost_equal(
            orig_ed["D"].real[:n], schur_ed["D"].real[:n]
        )
        np.testing.assert_array_almost_equal(
            np.abs(orig_ed["D"].imag[:n]), np.abs(schur_ed["D"].imag[:n])
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

        _check_compute_meta(mc)

    def test_compute_metastable_states_schur(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="brandts")
        mc.compute_metastable_states(n_states=2)

        _check_compute_meta(mc)

    def test_compute_metastable_invalid_cluster_key(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="brandts")
        with pytest.raises(KeyError):
            mc.compute_metastable_states(n_states=2, cluster_key="foobar")

    def test_compute_metastable_cache(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="brandts")

        mc.compute_metastable_states(n_states=2)

        assert mc.schur_vectors.shape[1] == 10
        assert mc._schur_matrix.shape == (10, 10)

    def test_set_main_states(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="brandts")

        mc.compute_metastable_states(n_states=2, n_cells=5)
        mc.set_main_states()

        _check_main_states(mc)

    def test_set_main_states_no_cells(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="brandts")

        mc.compute_metastable_states(n_states=2, n_cells=None)
        with pytest.raises(TypeError):
            mc.set_main_states(n_cells=None)

    def test_set_main_states_non_positive_cells(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="brandts")

        mc.compute_metastable_states(n_states=2, n_cells=None)
        with pytest.raises(ValueError):
            mc.set_main_states(n_cells=0)

    def test_set_main_states_invalid_name(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="brandts")

        mc.compute_metastable_states(n_states=2)
        with pytest.raises(KeyError):
            mc.set_main_states(names=["foobar"])

    def test_compute_main_states_invalid_method(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="brandts")

        mc.compute_metastable_states(n_states=2)
        with pytest.raises(ValueError):
            mc.compute_main_states(method="foobar")

    def test_compute_main_states_no_cells(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="brandts")

        mc.compute_metastable_states(n_states=2)
        with pytest.raises(TypeError):
            mc.compute_main_states(n_cells=None)

    def test_compute_main_states_non_positive_cells(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="brandts")

        mc.compute_metastable_states(n_states=2)
        with pytest.raises(ValueError):
            mc.compute_main_states(n_cells=0)

    def test_compute_main_states_eigengap(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="brandts")

        mc.compute_metastable_states(n_states=2)
        mc.compute_main_states(n_cells=5, method="eigengap")

        _check_main_states(mc)

    def test_compute_main_states_n_main_states(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="brandts")

        mc.compute_metastable_states(n_states=2)
        mc.compute_main_states(n_cells=5, method="top_n", n_main_states=1)

        _check_main_states(mc)

    def test_compute_main_states_min_self_prob(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="brandts")

        mc.compute_metastable_states(n_states=2)
        mc.compute_main_states(n_cells=5, method="min_self_prob", min_self_prob=0.5)

        _check_main_states(mc)

    def test_compute_main_states(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="brandts")

        mc.compute_metastable_states(n_states=2)
        mc.compute_main_states(n_cells=5)

        _check_main_states(mc)

    def test_compute_gdpt_no_schur(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.GPCCA(final_kernel)

        mc.compute_gdpt(method="brandts")

        assert "gdpt_pseudotime" in mc.adata.obs

    def test_compute_gdpt_no_iroot(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.GPCCA(final_kernel)
        mc.adata.uns.pop("iroot", None)

        with pytest.raises(KeyError):
            mc.compute_gdpt()

    def test_compute_gdpt_invalid_n_comps(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.GPCCA(final_kernel)

        with pytest.raises(ValueError):
            mc.compute_gdpt(n_components=1)

    def test_compute_gdpt_cellname_key_added(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="brandts")

        mc.compute_gdpt(key_added="foobar")

        assert "foobar" in mc.adata.obs

    def test_compute_gdpt_cellname_iroot(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.GPCCA(final_kernel)
        mc.adata.uns["iroot"] = mc.adata.obs_names[0]

        mc.compute_gdpt()

        assert "gdpt_pseudotime" in mc.adata.obs

    def test_compute_lineage_drivers_no_lineages(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="brandts")
        mc.compute_metastable_states(n_states=2)

        with pytest.raises(RuntimeError):
            mc.compute_lineage_drivers()

    def test_compute_lineage_drivers_invalid_lineages(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="brandts")
        mc.compute_metastable_states(n_states=2)
        mc.set_main_states()

        with pytest.raises(KeyError):
            mc.compute_lineage_drivers(use_raw=False, lin_names=["foo"])

    def test_compute_lineage_drivers_invalid_clusters(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="brandts")
        mc.compute_metastable_states(n_states=2)
        mc.set_main_states()

        with pytest.raises(KeyError):
            mc.compute_lineage_drivers(
                use_raw=False, cluster_key="clusters", clusters=["foo"]
            )

    def test_compute_lineage_drivers_normal_run(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="brandts")
        mc.compute_metastable_states(n_states=2)
        mc.set_main_states()
        mc.compute_lineage_drivers(use_raw=False, cluster_key="clusters")

        for lineage in ["0", "1"]:
            assert f"{Prefix.FORWARD} {lineage} corr" in mc.adata.var.keys()


class TestGPCCACopy:
    def test_copy_simple(self, adata_gpcca_fwd: Tuple[AnnData, cr.tl.GPCCA]):
        _, mc1 = adata_gpcca_fwd
        mc2 = mc1.copy()

        assert mc1 is not mc2
        assert mc1.adata is not mc2.adata
        assert mc1.kernel is not mc2.kernel

    def test_copy_deep(self, adata_gpcca_fwd: Tuple[AnnData, cr.tl.GPCCA]):
        _, mc1 = adata_gpcca_fwd
        mc2 = mc1.copy()

        assert mc1.irreducible == mc2.irreducible
        np.testing.assert_array_equal(mc1.recurrent_classes, mc2.recurrent_classes)
        np.testing.assert_array_equal(mc1.transient_classes, mc2.transient_classes)

        for k, v in mc1.eigendecomposition.items():
            if isinstance(v, np.ndarray):
                assert_array_nan_equal(v, mc2.eigendecomposition[k])
            else:
                assert v == mc2.eigendecomposition[k]

        assert mc1._gpcca is not mc2._gpcca

        np.testing.assert_array_equal(mc1.schur_vectors, mc2.schur_vectors)
        assert mc1.schur_vectors is not mc2.schur_vectors

        np.testing.assert_array_equal(mc1._schur_matrix, mc2._schur_matrix)
        assert mc1._schur_matrix is not mc2._schur_matrix

        np.testing.assert_array_equal(mc1.coarse_T.values, mc2.coarse_T.values)
        assert mc1.coarse_T is not mc2.coarse_T

        np.testing.assert_array_equal(
            mc1._coarse_init_dist.values, mc2._coarse_init_dist.values
        )
        assert mc1._coarse_init_dist is not mc2._coarse_init_dist

        np.testing.assert_array_equal(
            mc1.coarse_stationary_distribution.values,
            mc2.coarse_stationary_distribution.values,
        )
        assert (
            mc1.coarse_stationary_distribution is not mc2.coarse_stationary_distribution
        )

        assert_array_nan_equal(mc1.metastable_states, mc2.metastable_states)
        assert mc1.metastable_states is not mc2.metastable_states

        np.testing.assert_array_equal(mc1._meta_states_colors, mc2._meta_states_colors)
        assert mc1._meta_states_colors is not mc2._meta_states_colors

        np.testing.assert_array_equal(mc1._meta_lin_probs, mc2._meta_lin_probs)
        assert mc1._meta_lin_probs is not mc2._meta_lin_probs

        assert_array_nan_equal(mc1.main_states, mc2.main_states)
        assert mc1.main_states is not mc2.main_states

        np.testing.assert_array_equal(
            mc1.main_states_probabilities, mc2.main_states_probabilities
        )
        assert mc1.main_states_probabilities is not mc2.main_states_probabilities

        np.testing.assert_array_equal(
            mc1.lineage_probabilities, mc2.lineage_probabilities
        )
        assert mc1.lineage_probabilities is not mc2.lineage_probabilities

        np.testing.assert_array_equal(mc1.diff_potential, mc2.diff_potential)
        assert mc1.diff_potential is not mc2.diff_potential

        assert mc1._G2M_score == mc2._G2M_score
        assert mc1._S_score == mc2._S_score
        assert mc1._g2m_key == mc2._g2m_key
        assert mc1._s_key == mc2._s_key
        assert mc1._key_added == mc2._key_added
        assert mc1._is_sparse == mc2._is_sparse
