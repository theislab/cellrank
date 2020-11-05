# -*- coding: utf-8 -*-

import os
from copy import deepcopy
from typing import Tuple
from tempfile import TemporaryDirectory

import pytest
from _helpers import assert_array_nan_equal, assert_estimators_equal

from anndata import AnnData

import numpy as np
import pandas as pd

import cellrank as cr
from cellrank.tl.kernels import VelocityKernel, ConnectivityKernel
from cellrank.tl._constants import (
    Direction,
    DirPrefix,
    AbsProbKey,
    TermStatesKey,
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


def _check_compute_macro(mc: cr.tl.estimators.GPCCA) -> None:
    assert isinstance(mc._get(P.MACRO), pd.Series)
    assert len(mc._get(A.MACRO_COLORS)) == len(mc._get(P.MACRO).cat.categories)

    if "stationary_dist" in mc._get(P.EIG):  # one state
        assert isinstance(mc._get(P.MACRO_MEMBER), cr.tl.Lineage)
        assert mc._get(P.MACRO_MEMBER).shape[1] == 1
        np.testing.assert_array_almost_equal(mc._get(P.MACRO_MEMBER).X.sum(), 1.0)

        assert mc._get(P.COARSE_INIT_D) is None
        assert mc._get(P.SCHUR_MAT) is None
        assert mc._get(P.COARSE_STAT_D) is None
        assert mc._get(P.COARSE_T) is None
        assert mc._get(P.SCHUR) is None
    else:
        assert isinstance(mc._get(P.MACRO_MEMBER), cr.tl.Lineage)
        if mc._get(P.MACRO_MEMBER).shape[1] > 1:
            np.testing.assert_array_almost_equal(mc._get(P.MACRO_MEMBER).sum(1), 1.0)

        assert isinstance(mc._get(P.COARSE_INIT_D), pd.Series)
        assert isinstance(mc._get(P.SCHUR_MAT), np.ndarray)
        assert mc._get(P.COARSE_STAT_D) is None or isinstance(
            mc._get(P.COARSE_STAT_D), pd.Series
        )
        assert isinstance(mc._get(P.COARSE_T), pd.DataFrame)
        assert isinstance(mc._get(P.SCHUR), np.ndarray)


def _check_renaming_no_write_terminal(mc: cr.tl.estimators.GPCCA) -> None:
    assert mc._get(P.TERM) is None
    assert mc._get(P.TERM_PROBS) is None
    assert mc._get(A.TERM_ABS_PROBS) is None

    assert TermStatesKey.FORWARD.s not in mc.adata.obs
    assert _probs(TermStatesKey.FORWARD.s) not in mc.adata.obs
    assert _colors(TermStatesKey.FORWARD.s) not in mc.adata.uns
    assert _lin_names(TermStatesKey.FORWARD.s) not in mc.adata.uns


def _check_abs_probs(mc: cr.tl.estimators.GPCCA, has_main_states: bool = True):
    if has_main_states:
        assert isinstance(mc._get(P.TERM), pd.Series)
        assert_array_nan_equal(
            mc.adata.obs[str(TermStatesKey.FORWARD)], mc._get(P.TERM)
        )
        np.testing.assert_array_equal(
            mc.adata.uns[_colors(TermStatesKey.FORWARD)],
            mc._get(A.TERM_ABS_PROBS)[list(mc._get(P.TERM).cat.categories)].colors,
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

    assert_array_nan_equal(mc.adata.obs[TermStatesKey.FORWARD.s], mc._get(P.TERM))
    np.testing.assert_array_equal(
        mc.adata.obs[_probs(TermStatesKey.FORWARD)], mc._get(P.TERM_PROBS)
    )


class TestGPCCA:
    def test_compute_partition(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_partition()

        assert isinstance(mc.is_irreducible, bool)
        if not mc.is_irreducible:
            assert isinstance(mc.recurrent_classes, list)
            assert isinstance(mc.transient_classes, list)
            assert f"{TermStatesKey.FORWARD}_rec_classes" in mc.adata.obs
            assert f"{TermStatesKey.FORWARD}_trans_classes" in mc.adata.obs
        else:
            assert mc.recurrent_classes is None
            assert mc.transient_classes is None

    def test_compute_eigendecomposition(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_eigendecomposition(k=2, only_evals=True)

        _check_eigdecomposition(mc)

    def test_compute_schur_invalid_n_comps(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        with pytest.raises(ValueError):
            mc.compute_schur(n_components=1, method="krylov")

    def test_compute_schur_invalid_method(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        with pytest.raises(ValueError):
            mc.compute_schur(method="foobar")

    def test_compute_schur_invalid_eig_sort(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        with pytest.raises(ValueError):
            mc.compute_schur(which="foobar", method="krylov")

    def test_compute_schur_write_eigvals(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        _check_eigdecomposition(mc)

    def test_compute_schur_write_eigvals_similar_to_orig_eigdecomp(
        self, adata_large: AnnData
    ):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
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

    def test_compute_macrostates_no_eig(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        with pytest.raises(RuntimeError):
            mc.compute_macrostates(n_states=None)

    def test_compute_macrostates_1_state_no_eig(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_macrostates(n_states=1)

    def test_compute_macro_none_states(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_eigendecomposition(only_evals=True)
        mc.compute_macrostates(n_states=None)

        _check_compute_macro(mc)

    def test_compute_macrostates_schur(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_macrostates(n_states=2)

        _check_compute_macro(mc)

    def test_compute_macrostates_min_chi_too_low_min(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        with pytest.raises(ValueError):
            mc.compute_macrostates(n_states=[1, 4], use_min_chi=True)

    def test_compute_macrostates_min_chi_inverted_range(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_macrostates(n_states=[4, 2], use_min_chi=True)

        _check_compute_macro(mc)

    def test_compute_macrostates_min_chi_range_same_values(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_macrostates(n_states=[2, 2], use_min_chi=True)

        _check_compute_macro(mc)

    def test_compute_macrostates_min_chi_dict_wrong_keys(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        with pytest.raises(KeyError):
            mc.compute_macrostates(n_states={"foo": 2, "max": 3}, use_min_chi=True)

    def test_compute_macrostates_min_chi_normal_run(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_macrostates(n_states=[2, 4], use_min_chi=True)

        _check_compute_macro(mc)

    def test_compute_macro_invalid_cluster_key(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        with pytest.raises(KeyError):
            mc.compute_macrostates(n_states=2, cluster_key="foobar")

    def test_compute_macro_cache(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=11, method="krylov")

        mc.compute_macrostates(n_states=2)

        assert mc._get(P.SCHUR).shape[1] == 11
        assert mc._get(P.SCHUR_MAT).shape == (11, 11)

    def test_set_initial_states_from_forward(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large, backward=False).compute_transition_matrix(
            softmax_scale=4
        )
        ck = ConnectivityKernel(adata_large, backward=False).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=2, n_cells=5)
        obsm_keys = set(mc.adata.obsm.keys())
        mc._set_initial_states_from_macrostates("0")

        key = TermStatesKey.BACKWARD.s

        assert key in mc.adata.obs
        np.testing.assert_array_equal(mc.adata.obs[key].cat.categories, ["0"])
        assert _probs(key) in mc.adata.obs
        assert _colors(key) in mc.adata.uns
        assert _lin_names(key) in mc.adata.uns

        # make sure that we don't write anything there - it's useless
        assert set(mc.adata.obsm.keys()) == obsm_keys

    def test_compute_initial_states_from_forward_no_macro(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large, backward=False).compute_transition_matrix(
            softmax_scale=4
        )
        ck = ConnectivityKernel(adata_large, backward=False).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        with pytest.raises(RuntimeError):
            mc._compute_initial_states(1)

    def test_compute_initial_states_from_forward_too_many_states(
        self, adata_large: AnnData
    ):
        vk = VelocityKernel(adata_large, backward=False).compute_transition_matrix(
            softmax_scale=4
        )
        ck = ConnectivityKernel(adata_large, backward=False).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=2, n_cells=5)
        with pytest.raises(ValueError):
            mc._compute_initial_states(42)

    def test_compute_initial_states_from_forward_too_few_states(
        self, adata_large: AnnData
    ):
        vk = VelocityKernel(adata_large, backward=False).compute_transition_matrix(
            softmax_scale=4
        )
        ck = ConnectivityKernel(adata_large, backward=False).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=2, n_cells=5)
        with pytest.raises(ValueError):
            mc._compute_initial_states(0)

    def test_compute_initial_states_from_forward_no_stat_dist(
        self, adata_large: AnnData
    ):
        vk = VelocityKernel(adata_large, backward=False).compute_transition_matrix(
            softmax_scale=4
        )
        ck = ConnectivityKernel(adata_large, backward=False).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=2, n_cells=5)
        mc._coarse_stat_dist = None

        with pytest.raises(ValueError):
            mc._compute_initial_states(0)

    def test_compute_initial_states_from_forward_normal_run(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large, backward=False).compute_transition_matrix(
            softmax_scale=4
        )
        ck = ConnectivityKernel(adata_large, backward=False).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=2, n_cells=5)
        obsm_keys = set(mc.adata.obsm.keys())
        expected = mc._get(P.COARSE_STAT_D).index[np.argmin(mc._get(P.COARSE_STAT_D))]

        mc._compute_initial_states(1)

        key = TermStatesKey.BACKWARD.s

        assert key in mc.adata.obs
        np.testing.assert_array_equal(mc.adata.obs[key].cat.categories, [expected])
        assert _probs(key) in mc.adata.obs
        assert _colors(key) in mc.adata.uns
        assert _lin_names(key) in mc.adata.uns

        # make sure that we don't write anything there - it's useless
        assert set(mc.adata.obsm.keys()) == obsm_keys

    def test_set_terminal_states_from_macrostates(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=2, n_cells=5)
        mc.set_terminal_states_from_macrostates()
        mc.compute_absorption_probabilities()

        _check_abs_probs(mc)

    def test_set_terminal_states_from_macrostates_rename_states_invalid(
        self, adata_large: AnnData
    ):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=2, n_cells=5)
        with pytest.raises(ValueError):
            mc.set_terminal_states_from_macrostates({"0": "1"})
        _check_renaming_no_write_terminal(mc)

    def test_set_terminal_states_from_macrostates_rename_not_unique_new_names(
        self, adata_large: AnnData
    ):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=3, n_cells=5)
        with pytest.raises(ValueError):
            mc.set_terminal_states_from_macrostates({"0": "1", "1": "1"})
        _check_renaming_no_write_terminal(mc)

    def test_set_terminal_states_from_macrostates_rename_overlapping_old_keys(
        self, adata_large: AnnData
    ):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=3, n_cells=5)
        with pytest.raises(ValueError):
            mc.set_terminal_states_from_macrostates({"0, 1": "foo", "1, 2": "bar"})
        _check_renaming_no_write_terminal(mc)

    def test_set_terminal_states_from_macrostates_rename_states(
        self, adata_large: AnnData
    ):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=2, n_cells=5)
        mc.set_terminal_states_from_macrostates({"0": "foo"})
        mc.compute_absorption_probabilities()

        _check_abs_probs(mc)

        np.testing.assert_array_equal(mc._get(P.TERM).cat.categories, ["foo"])
        np.testing.assert_array_equal(mc._get(A.TERM_ABS_PROBS).names, ["foo"])

    def test_set_terminal_states_from_macrostates_join_and_rename_states(
        self, adata_large: AnnData
    ):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=3, n_cells=5)
        mc.set_terminal_states_from_macrostates({"0, 1": "foo", "2": "bar"})
        mc.compute_absorption_probabilities()

        _check_abs_probs(mc)

        np.testing.assert_array_equal(mc._get(P.TERM).cat.categories, ["foo", "bar"])
        np.testing.assert_array_equal(mc._get(A.TERM_ABS_PROBS).names, ["foo", "bar"])

    def test_set_terminal_states_from_macrostates_no_cells(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=2, n_cells=None)
        with pytest.raises(TypeError):
            mc.set_terminal_states_from_macrostates(n_cells=None)

    def test_set_terminal_states_from_macrostates_non_positive_cells(
        self, adata_large: AnnData
    ):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=2, n_cells=None)
        with pytest.raises(ValueError):
            mc.set_terminal_states_from_macrostates(n_cells=0)

    def test_set_terminal_states_from_macrostates_invalid_name(
        self, adata_large: AnnData
    ):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=2)
        with pytest.raises(KeyError):
            mc.set_terminal_states_from_macrostates(names=["foobar"])

    def test_compute_terminal_states_invalid_method(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=2)
        with pytest.raises(ValueError):
            mc.compute_terminal_states(method="foobar")

    def test_compute_terminal_states_no_cells(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=2)
        with pytest.raises(TypeError):
            mc.compute_terminal_states(n_cells=None)

    def test_compute_terminal_states_non_positive_cells(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=2)
        with pytest.raises(ValueError):
            mc.compute_terminal_states(n_cells=0)

    def test_compute_terminal_states_eigengap(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=2)
        mc.compute_terminal_states(n_cells=5, method="eigengap")
        mc.compute_absorption_probabilities()

        _check_abs_probs(mc)

    def test_compute_terminal_states_n_main_states(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=2)
        mc.compute_terminal_states(n_cells=5, method="top_n", n_states=1)
        mc.compute_absorption_probabilities()

        _check_abs_probs(mc)

    def test_compute_terminal_states_stability(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck
        thresh = 0.5

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=5)
        mc.compute_terminal_states(
            n_cells=5, method="stability", stability_threshold=thresh
        )
        mc.compute_absorption_probabilities()

        coarse_T = mc._get(P.COARSE_T)
        self_probs = pd.Series(np.diag(coarse_T), index=coarse_T.columns)
        names = self_probs[self_probs.values >= thresh].index

        np.testing.assert_array_equal(set(names), set(mc._get(P.TERM).cat.categories))
        _check_abs_probs(mc)

    def test_compute_terminal_states_no_selected(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=2)
        with pytest.raises(ValueError):
            mc.compute_terminal_states(
                n_cells=5, method="stability", stability_threshold=42
            )

    def test_compute_terminal_states_too_many_cells(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=2)
        with pytest.raises(ValueError):
            mc.compute_terminal_states(n_cells=4200)

    def test_compute_terminal_states_default(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=2)
        mc.compute_terminal_states(n_cells=5)
        mc.compute_absorption_probabilities()

        _check_abs_probs(mc)

    def test_compute_gdpt_no_schur(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)

        mc.compute_gdpt(method="krylov")

        assert "gdpt_pseudotime" in mc.adata.obs

    def test_compute_gdpt_no_iroot(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.adata.uns.pop("iroot", None)

        with pytest.raises(KeyError):
            mc.compute_gdpt()

    def test_compute_gdpt_invalid_n_comps(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)

        with pytest.raises(ValueError):
            mc.compute_gdpt(n_components=1)

    def test_compute_gdpt_cellname_key_added(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_gdpt(key_added="foobar")

        assert "foobar" in mc.adata.obs

    def test_compute_gdpt_cellname_iroot(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.adata.uns["iroot"] = mc.adata.obs_names[0]

        mc.compute_gdpt()

        assert "gdpt_pseudotime" in mc.adata.obs

    def test_compute_lineage_drivers_invalid_lineages(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_macrostates(n_states=2)
        mc.set_terminal_states_from_macrostates()
        mc.compute_absorption_probabilities()

        with pytest.raises(KeyError):
            mc.compute_lineage_drivers(use_raw=False, lineages=["foo"])

    def test_compute_lineage_drivers_invalid_clusters(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_macrostates(n_states=2)
        mc.set_terminal_states_from_macrostates()
        mc.compute_absorption_probabilities()

        with pytest.raises(KeyError):
            mc.compute_lineage_drivers(
                use_raw=False, cluster_key="clusters", clusters=["foo"]
            )

    def test_compute_lineage_drivers_normal_run(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_macrostates(n_states=2)
        mc.set_terminal_states_from_macrostates()
        mc.compute_absorption_probabilities()
        mc.compute_lineage_drivers(use_raw=False, cluster_key="clusters")

        for lineage in ["0", "1"]:
            assert np.all(mc.adata.var[f"{DirPrefix.FORWARD} {lineage} corr"] >= -1.0)
            assert np.all(mc.adata.var[f"{DirPrefix.FORWARD} {lineage} corr"] <= 1.0)

            assert np.all(mc.adata.var[f"{DirPrefix.FORWARD} {lineage} qval"] >= 0)
            assert np.all(mc.adata.var[f"{DirPrefix.FORWARD} {lineage} qval"] <= 1.0)

    def test_plot_lineage_drivers_not_computed(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_macrostates(n_states=2)
        mc.set_terminal_states_from_macrostates()
        mc.compute_absorption_probabilities()

        with pytest.raises(RuntimeError):
            mc.plot_lineage_drivers("0")

    def test_plot_lineage_drivers_invalid_name(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_macrostates(n_states=2)
        mc.set_terminal_states_from_macrostates()
        mc.compute_absorption_probabilities()
        mc.compute_lineage_drivers(use_raw=False, cluster_key="clusters")

        with pytest.raises(KeyError):
            mc.plot_lineage_drivers("foo", use_raw=False)

    def test_plot_lineage_drivers_invalid_n_genes(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_macrostates(n_states=2)
        mc.set_terminal_states_from_macrostates()
        mc.compute_absorption_probabilities()
        mc.compute_lineage_drivers(use_raw=False, cluster_key="clusters")

        with pytest.raises(ValueError):
            mc.plot_lineage_drivers("0", use_raw=False, n_genes=0)

    def test_plot_lineage_drivers_normal_run(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_macrostates(n_states=2)
        mc.set_terminal_states_from_macrostates()
        mc.compute_absorption_probabilities()
        mc.compute_lineage_drivers(use_raw=False, cluster_key="clusters")

        mc.plot_lineage_drivers("0", use_raw=False)


class TestGPCCAIO:
    def test_copy(self, adata_gpcca_fwd: Tuple[AnnData, cr.tl.estimators.GPCCA]):
        _, mc1 = adata_gpcca_fwd
        mc2 = mc1.copy()

        assert_estimators_equal(mc1, mc2, copy=True)

    def test_write_ext(self, adata_gpcca_fwd: Tuple[AnnData, cr.tl.estimators.GPCCA]):
        _, mc = adata_gpcca_fwd

        with TemporaryDirectory() as tmpdir:
            fname = "foo"
            mc.write(os.path.join(tmpdir, fname), ext="bar")

            assert os.path.isfile(os.path.join(tmpdir, f"foo.bar"))

    def test_write_no_ext(
        self, adata_gpcca_fwd: Tuple[AnnData, cr.tl.estimators.GPCCA]
    ):
        _, mc = adata_gpcca_fwd

        with TemporaryDirectory() as tmpdir:
            fname = "foo"
            mc.write(os.path.join(tmpdir, fname), ext=None)

            assert os.path.isfile(os.path.join(tmpdir, f"foo"))

    def test_write_ext_with_dot(
        self, adata_gpcca_fwd: Tuple[AnnData, cr.tl.estimators.GPCCA]
    ):
        _, mc = adata_gpcca_fwd

        with TemporaryDirectory() as tmpdir:
            fname = "foo"
            mc.write(os.path.join(tmpdir, fname), ext=".bar")

            assert os.path.isfile(os.path.join(tmpdir, f"foo.bar"))

    def test_read(self, adata_gpcca_fwd: Tuple[AnnData, cr.tl.estimators.GPCCA]):
        _, mc1 = adata_gpcca_fwd

        with TemporaryDirectory() as tmpdir:
            mc1.write(os.path.join(tmpdir, "foo"))
            mc2 = cr.tl.estimators.GPCCA.read(os.path.join(tmpdir, "foo.pickle"))

        assert_estimators_equal(mc1, mc2)
