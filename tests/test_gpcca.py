# -*- coding: utf-8 -*-

from copy import deepcopy
from typing import Tuple

import pytest
from _helpers import assert_array_nan_equal

from anndata import AnnData

import numpy as np
import pandas as pd

import cellrank as cr
from cellrank.tl.kernels import VelocityKernel, ConnectivityKernel
from cellrank.tl._constants import (
    Direction,
    DirPrefix,
    AbsProbKey,
    FinalStatesKey,
    _dp,
    _probs,
    _colors,
    _lin_names,
)
from cellrank.tl.estimators._constants import A, P


def _check_eigdecomposition(mc: cr.tl.estimators.GPCCA) -> None:
    assert isinstance(mc._get(P.EIG), dict)
    assert set(mc._get(P.EIG).keys()) == {
        "D",
        "eigengap",
        "params",
    }
    assert "V_l" not in mc._get(P.EIG)
    assert "V_r" not in mc._get(P.EIG)
    assert "stationary_dist" not in mc._get(P.EIG)

    assert f"eig_{Direction.FORWARD}" in mc.adata.uns.keys()


def _check_compute_meta(mc: cr.tl.estimators.GPCCA) -> None:
    assert isinstance(mc._get(P.META), pd.Series)
    assert len(mc._get(A.META_COLORS)) == len(mc._get(P.META).cat.categories)

    if "stationary_dist" in mc._get(P.EIG):  # one state
        assert isinstance(mc._get(P.META_PROBS), cr.tl.Lineage)
        assert mc._get(P.META_PROBS).shape[1] == 1
        np.testing.assert_array_almost_equal(mc._get(P.META_PROBS).X.sum(), 1.0)

        assert mc._get(P.COARSE_INIT_D) is None
        assert mc._get(P.SCHUR_MAT) is None
        assert mc._get(P.COARSE_STAT_D) is None
        assert mc._get(P.COARSE_T) is None
        assert mc._get(P.SCHUR) is None
    else:
        assert isinstance(mc._get(P.META_PROBS), cr.tl.Lineage)
        if mc._get(P.META_PROBS).shape[1] > 1:
            np.testing.assert_array_almost_equal(mc._get(P.META_PROBS).sum(1), 1.0)

        assert isinstance(mc._get(P.COARSE_INIT_D), pd.Series)
        assert isinstance(mc._get(P.SCHUR_MAT), np.ndarray)
        assert mc._get(P.COARSE_STAT_D) is None or isinstance(
            mc._get(P.COARSE_STAT_D), pd.Series
        )
        assert isinstance(mc._get(P.COARSE_T), pd.DataFrame)
        assert isinstance(mc._get(P.SCHUR), np.ndarray)


def _check_abs_probs(mc: cr.tl.estimators.GPCCA, has_main_states: bool = True):
    if has_main_states:
        assert isinstance(mc._get(P.FIN), pd.Series)
        assert_array_nan_equal(
            mc.adata.obs[str(FinalStatesKey.FORWARD)], mc._get(P.FIN)
        )
        np.testing.assert_array_equal(
            mc.adata.uns[_colors(FinalStatesKey.FORWARD)],
            mc._get(A.FIN_ABS_PROBS)[list(mc._get(P.FIN).cat.categories)].colors,
        )

    assert isinstance(mc._get(P.DIFF_POT), pd.Series)
    assert isinstance(mc._get(P.ABS_PROBS), cr.tl.Lineage)
    np.testing.assert_array_almost_equal(mc._get(P.ABS_PROBS).sum(1), 1.0)

    np.testing.assert_array_equal(
        mc.adata.obsm[str(AbsProbKey.FORWARD)], mc._get(P.ABS_PROBS).X
    )
    np.testing.assert_array_equal(
        mc.adata.uns[_lin_names(AbsProbKey.FORWARD)], mc._get(P.ABS_PROBS).names
    )
    np.testing.assert_array_equal(
        mc.adata.uns[_colors(AbsProbKey.FORWARD)], mc._get(P.ABS_PROBS).colors
    )

    np.testing.assert_array_equal(
        mc.adata.obs[_dp(AbsProbKey.FORWARD)], mc._get(P.DIFF_POT)
    )

    assert_array_nan_equal(mc.adata.obs[FinalStatesKey.FORWARD.s], mc._get(P.FIN))
    np.testing.assert_array_equal(
        mc.adata.obs[_probs(FinalStatesKey.FORWARD)], mc._get(P.FIN_PROBS)
    )


class TestGPCCA:
    def test_compute_partition(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(final_kernel)
        mc.compute_partition()

        assert isinstance(mc.is_irreducible, bool)
        if not mc.is_irreducible:
            assert isinstance(mc.recurrent_classes, list)
            assert isinstance(mc.transient_classes, list)
            assert f"{FinalStatesKey.FORWARD}_rec_classes" in mc.adata.obs
            assert f"{FinalStatesKey.FORWARD}_trans_classes" in mc.adata.obs
        else:
            assert mc.recurrent_classes is None
            assert mc.transient_classes is None

    def test_compute_eigendecomposition(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(final_kernel)
        mc.compute_eigendecomposition(k=2, only_evals=True)

        _check_eigdecomposition(mc)

    def test_compute_schur_invalid_n_comps(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(final_kernel)
        with pytest.raises(ValueError):
            mc.compute_schur(n_components=1, method="krylov")

    def test_compute_schur_invalid_method(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(final_kernel)
        with pytest.raises(ValueError):
            mc.compute_schur(method="foobar")

    def test_compute_schur_invalid_eig_sort(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(final_kernel)
        with pytest.raises(ValueError):
            mc.compute_schur(which="foobar", method="krylov")

    def test_compute_schur_write_eigvals(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        _check_eigdecomposition(mc)

    def test_compute_schur_write_eigvals_similar_to_orig_eigdecomp(
        self, adata_large: AnnData
    ):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(final_kernel)
        mc.compute_eigendecomposition(k=10, only_evals=True)

        _check_eigdecomposition(mc)
        orig_ed = deepcopy(mc._get(P.EIG))

        mc._set(A.EIG, None)
        mc.compute_schur(n_components=10, method="krylov")

        _check_eigdecomposition(mc)
        schur_ed = mc._get(P.EIG)

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

        mc = cr.tl.estimators.GPCCA(final_kernel)
        with pytest.raises(RuntimeError):
            mc.compute_metastable_states(n_states=None)

    def test_compute_metastable_states_1_state_no_eig(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(final_kernel)
        mc.compute_metastable_states(n_states=1)

    def test_compute_metastable_none_states(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(final_kernel)
        mc.compute_eigendecomposition(only_evals=True)
        mc.compute_metastable_states(n_states=None)

        _check_compute_meta(mc)

    def test_compute_metastable_states_schur(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_metastable_states(n_states=2)

        _check_compute_meta(mc)

    def test_compute_metastable_invalid_cluster_key(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        with pytest.raises(KeyError):
            mc.compute_metastable_states(n_states=2, cluster_key="foobar")

    def test_compute_metastable_cache(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(final_kernel)
        mc.compute_schur(n_components=11, method="krylov")

        mc.compute_metastable_states(n_states=2)

        assert mc._get(P.SCHUR).shape[1] == 11
        assert mc._get(P.SCHUR_MAT).shape == (11, 11)

    def test_set_final_states_from_metastable_states(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_metastable_states(n_states=2, n_cells=5)
        mc.set_final_states_from_metastable_states()
        mc.compute_absorption_probabilities()

        _check_abs_probs(mc)

    def test_set_final_states_from_metastable_states_no_cells(
        self, adata_large: AnnData
    ):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_metastable_states(n_states=2, n_cells=None)
        with pytest.raises(TypeError):
            mc.set_final_states_from_metastable_states(n_cells=None)

    def test_set_final_states_from_metastable_states_non_positive_cells(
        self, adata_large: AnnData
    ):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_metastable_states(n_states=2, n_cells=None)
        with pytest.raises(ValueError):
            mc.set_final_states_from_metastable_states(n_cells=0)

    def test_set_final_states_from_metastable_states_invalid_name(
        self, adata_large: AnnData
    ):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_metastable_states(n_states=2)
        with pytest.raises(KeyError):
            mc.set_final_states_from_metastable_states(names=["foobar"])

    def test_compute_final_states_invalid_method(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_metastable_states(n_states=2)
        with pytest.raises(ValueError):
            mc.compute_final_states(method="foobar")

    def test_compute_final_states_no_cells(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_metastable_states(n_states=2)
        with pytest.raises(TypeError):
            mc.compute_final_states(n_cells=None)

    def test_compute_final_states_non_positive_cells(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_metastable_states(n_states=2)
        with pytest.raises(ValueError):
            mc.compute_final_states(n_cells=0)

    def test_compute_final_states_eigengap(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_metastable_states(n_states=2)
        mc.compute_final_states(n_cells=5, method="eigengap")
        mc.compute_absorption_probabilities()

        _check_abs_probs(mc)

    def test_compute_final_states_n_main_states(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_metastable_states(n_states=2)
        mc.compute_final_states(n_cells=5, method="top_n", n_main_states=1)
        mc.compute_absorption_probabilities()

        _check_abs_probs(mc)

    def test_compute_final_states_min_self_prob(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_metastable_states(n_states=2)
        mc.compute_final_states(n_cells=5, method="min_self_prob", min_self_prob=0.5)
        mc.compute_absorption_probabilities()

        _check_abs_probs(mc)

    def test_compute_final_states(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_metastable_states(n_states=2)
        mc.compute_final_states(n_cells=5)
        mc.compute_absorption_probabilities()

        _check_abs_probs(mc)

    def test_compute_gdpt_no_schur(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(final_kernel)

        mc.compute_gdpt(method="krylov")

        assert "gdpt_pseudotime" in mc.adata.obs

    def test_compute_gdpt_no_iroot(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(final_kernel)
        mc.adata.uns.pop("iroot", None)

        with pytest.raises(KeyError):
            mc.compute_gdpt()

    def test_compute_gdpt_invalid_n_comps(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(final_kernel)

        with pytest.raises(ValueError):
            mc.compute_gdpt(n_components=1)

    def test_compute_gdpt_cellname_key_added(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_gdpt(key_added="foobar")

        assert "foobar" in mc.adata.obs

    def test_compute_gdpt_cellname_iroot(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(final_kernel)
        mc.adata.uns["iroot"] = mc.adata.obs_names[0]

        mc.compute_gdpt()

        assert "gdpt_pseudotime" in mc.adata.obs

    def test_compute_lineage_drivers_invalid_lineages(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_metastable_states(n_states=2)
        mc.set_final_states_from_metastable_states()
        mc.compute_absorption_probabilities()

        with pytest.raises(KeyError):
            mc.compute_lineage_drivers(use_raw=False, lineages=["foo"])

    def test_compute_lineage_drivers_invalid_clusters(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_metastable_states(n_states=2)
        mc.set_final_states_from_metastable_states()
        mc.compute_absorption_probabilities()

        with pytest.raises(KeyError):
            mc.compute_lineage_drivers(
                use_raw=False, cluster_key="clusters", clusters=["foo"]
            )

    def test_compute_lineage_drivers_normal_run(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_metastable_states(n_states=2)
        mc.set_final_states_from_metastable_states()
        mc.compute_absorption_probabilities()
        mc.compute_lineage_drivers(use_raw=False, cluster_key="clusters")

        for lineage in ["0", "1"]:
            assert f"{DirPrefix.FORWARD} {lineage}" in mc.adata.var.keys()

    def test_plot_lineage_drivers_not_computed(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_metastable_states(n_states=2)
        mc.set_final_states_from_metastable_states()
        mc.compute_absorption_probabilities()

        with pytest.raises(RuntimeError):
            mc.plot_lineage_drivers("0")

    def test_plot_lineage_drivers_invalid_name(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_metastable_states(n_states=2)
        mc.set_final_states_from_metastable_states()
        mc.compute_absorption_probabilities()
        mc.compute_lineage_drivers(use_raw=False, cluster_key="clusters")

        with pytest.raises(KeyError):
            mc.plot_lineage_drivers("foo", use_raw=False)

    def test_plot_lineage_drivers_invalid_n_genes(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_metastable_states(n_states=2)
        mc.set_final_states_from_metastable_states()
        mc.compute_absorption_probabilities()
        mc.compute_lineage_drivers(use_raw=False, cluster_key="clusters")

        with pytest.raises(ValueError):
            mc.plot_lineage_drivers("0", use_raw=False, n_genes=0)

    def test_plot_lineage_drivers_normal_run(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(final_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_metastable_states(n_states=2)
        mc.set_final_states_from_metastable_states()
        mc.compute_absorption_probabilities()
        mc.compute_lineage_drivers(use_raw=False, cluster_key="clusters")

        mc.plot_lineage_drivers("0", use_raw=False)


class TestGPCCACopy:
    def test_copy_simple(self, adata_gpcca_fwd: Tuple[AnnData, cr.tl.estimators.GPCCA]):
        _, mc1 = adata_gpcca_fwd
        mc2 = mc1.copy()

        assert mc1 is not mc2
        assert mc1.adata is not mc2.adata
        assert mc1.kernel is not mc2.kernel
