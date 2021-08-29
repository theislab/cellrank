from typing import List, Tuple, Union

import os
import pytest
from copy import deepcopy
from _helpers import assert_array_nan_equal, assert_estimators_equal
from tempfile import TemporaryDirectory

import cellrank as cr
from anndata import AnnData
from cellrank._key import Key
from cellrank.tl.kernels import VelocityKernel, ConnectivityKernel

import numpy as np
import pandas as pd
from pandas.testing import assert_series_equal


def _check_eigdecomposition(mc: cr.tl.estimators.GPCCA) -> None:
    assert isinstance(mc.eigendecomposition, dict)
    assert set(mc.eigendecomposition.keys()) == {
        "D",
        "eigengap",
        "params",
    }
    assert "V_l" not in mc.eigendecomposition
    assert "V_r" not in mc.eigendecomposition
    assert "stationary_dist" not in mc.eigendecomposition

    assert Key.uns.eigen(mc.backward) in mc.adata.uns


def _check_compute_macro(mc: cr.tl.estimators.GPCCA) -> None:
    assert isinstance(mc.macrostates, pd.Series)
    assert len(mc._macrostates_colors) == len(mc.macrostates.cat.categories)

    if "stationary_dist" in mc.eigendecomposition:  # one state
        assert isinstance(mc.macrostates_memberships, cr.tl.Lineage)
        assert mc.macrostates_memberships.shape[1] == 1
        np.testing.assert_allclose(mc.macrostates_memberships.X.sum(), 1.0)

        assert mc.schur_matrix is None
        assert mc.schur_vectors is None
        assert mc.coarse_initial_distribution is None
        assert mc.coarse_stationary_distribution is None
        assert mc.coarse_T is None
    else:
        assert isinstance(mc.macrostates_memberships, cr.tl.Lineage)
        if mc.macrostates_memberships.shape[1] > 1:
            np.testing.assert_allclose(mc.macrostates_memberships.sum(1), 1.0)

        assert isinstance(mc.schur_matrix, np.ndarray)
        assert isinstance(mc.schur_vectors, np.ndarray)
        assert isinstance(mc.coarse_initial_distribution, pd.Series)
        assert isinstance(mc.coarse_T, pd.DataFrame)
        np.testing.assert_array_equal(mc.coarse_T.index, mc.coarse_T.columns)
        np.testing.assert_array_equal(
            mc.coarse_T.index, mc.coarse_stationary_distribution.index
        )
        if mc.coarse_stationary_distribution is not None:
            assert isinstance(mc.coarse_stationary_distribution, pd.Series)
            np.testing.assert_array_equal(
                mc.coarse_T.index, mc.coarse_stationary_distribution.index
            )


def _check_renaming_no_write_terminal(mc: cr.tl.estimators.GPCCA) -> None:
    assert mc.terminal_states is None
    assert mc.terminal_states_probabilities is None
    assert mc.terminal_states_memberships is None

    key = Key.obs.term_states(mc.backward)
    assert key not in mc.adata.obs
    assert Key.obs.probs(key) not in mc.adata.obs
    assert Key.uns.colors(key) not in mc.adata.uns


def _check_abs_probs(mc: cr.tl.estimators.GPCCA) -> None:
    # fmt: off
    # macrostates
    assert isinstance(mc.macrostates, pd.Series)
    assert isinstance(mc.macrostates_memberships, cr.tl.Lineage)
    np.testing.assert_array_equal(mc._macrostates_colors, mc.macrostates_memberships.colors)

    # term states
    key = Key.obs.term_states(mc.backward)
    assert isinstance(mc.terminal_states, pd.Series)
    # TODO: assert series equal
    assert_array_nan_equal(mc.adata.obs[key], mc.terminal_states)
    np.testing.assert_array_equal(mc.adata.uns[Key.uns.colors(key)], mc.absorption_probabilities.colors)
    np.testing.assert_array_equal(mc.adata.uns[Key.uns.colors(key)], mc._term_states_colors)
    assert isinstance(mc.terminal_states_probabilities, pd.Series)
    np.testing.assert_array_equal(mc.adata.obs[Key.obs.probs(key)], mc.terminal_states_probabilities)

    # abs probs
    key = Key.obsm.abs_probs(mc.backward)
    assert isinstance(mc.absorption_probabilities, cr.tl.Lineage)
    np.testing.assert_array_almost_equal(mc.absorption_probabilities.sum(1), 1.0)
    assert isinstance(mc.adata.obsm[key], cr.tl.Lineage)
    np.testing.assert_array_equal(mc.adata.obsm[key], mc.absorption_probabilities.X)

    # priming
    key = Key.obs.priming_degree(mc.backward)
    assert isinstance(mc.priming_degree, pd.Series)
    np.testing.assert_array_equal(mc.adata.obs[key], mc.priming_degree)
    # fmt: on


class TestGPCCA:
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
        orig_ed = deepcopy(mc.eigendecomposition)

        mc._eigendecomposition = None
        mc.compute_schur(n_components=10, method="krylov")

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
        mc.compute_macrostates(n_states=[1, 1])

        _check_compute_macro(mc)

    def test_compute_macrostates_min_chi_inverted_range(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_macrostates(n_states=[4, 2])

        _check_compute_macro(mc)

    def test_compute_macrostates_min_chi_range_same_values(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_macrostates(n_states=[2, 2])

        _check_compute_macro(mc)

    def test_compute_macrostates_min_chi_normal_run(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_macrostates(n_states=[2, 4])

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

        assert mc.schur_vectors.shape[1] == 11
        np.testing.assert_array_equal(mc.schur_matrix.shape, [11, 11])

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
        # TODO
        mc._set_initial_states_from_macrostates("0")

        key = Key.obs.term_states(not mc.backward)

        assert key in mc.adata.obs
        np.testing.assert_array_equal(mc.adata.obs[key].cat.categories, ["0"])
        assert Key.obs.probs(key) in mc.adata.obs
        assert Key.uns.colors(key) in mc.adata.uns

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
            # TODO
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
        csd = mc.coarse_stationary_distribution
        expected = csd.index[np.argmin(csd)]

        mc._compute_initial_states(1)

        key = Key.obs.term_states(not mc.backward)

        assert key in mc.adata.obs
        np.testing.assert_array_equal(mc.adata.obs[key].cat.categories, [expected])
        assert Key.obs.probs(key) in mc.adata.obs
        assert Key.uns.colors(key) in mc.adata.uns

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
        mc.compute_lineage_priming()

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
        mc.compute_lineage_priming()

        _check_abs_probs(mc)

        np.testing.assert_array_equal(mc.terminal_states.cat.categories, ["foo"])
        np.testing.assert_array_equal(mc.terminal_states_memberships.names, ["foo"])

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
        mc.compute_lineage_priming()

        _check_abs_probs(mc)

        np.testing.assert_array_equal(mc.terminal_states.cat.categories, ["foo", "bar"])
        np.testing.assert_array_equal(
            mc.terminal_states_memberships.names, ["foo", "bar"]
        )

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

    @pytest.mark.parametrize("values", ["Astrocytes", ["Astrocytes", "OPC"]])
    def test_set_terminal_states_clusters(
        self, adata_large: AnnData, values: Union[str, List[str]]
    ):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        to_remove = list(
            set(adata_large.obs["clusters"].cat.categories)
            - ({values} if isinstance(values, str) else set(values))
        )
        expected = adata_large.obs["clusters"].cat.remove_categories(to_remove)

        mc = cr.tl.estimators.GPCCA(terminal_kernel)

        mc.set_terminal_states({"clusters": values})
        pd.testing.assert_series_equal(
            expected, mc.terminal_states, check_category_order=False, check_names=False
        )

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
        mc.compute_lineage_priming()

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
        mc.compute_lineage_priming()

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
        mc.compute_lineage_priming()

        coarse_T = mc.coarse_T
        self_probs = pd.Series(np.diag(coarse_T), index=coarse_T.columns)
        names = self_probs[self_probs.values >= thresh].index

        np.testing.assert_array_equal(
            set(names), set(mc.terminal_states.cat.categories)
        )
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
        mc.compute_lineage_priming()

        _check_abs_probs(mc)

    def test_compute_gdpt_no_schur(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)

        res = mc.compute_gdpt(method="krylov")
        assert isinstance(res, pd.Series)
        np.testing.assert_array_equal(res.index, adata_large.obs_names)

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

        # auto-bumped to 2
        res = mc.compute_gdpt(n_components=1)
        assert isinstance(res, pd.Series)
        np.testing.assert_array_equal(res.index, adata_large.obs_names)

    def test_compute_gdpt_cellname_iroot(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.adata.uns["iroot"] = mc.adata.obs_names[0]

        res = mc.compute_gdpt()
        assert isinstance(res, pd.Series)
        np.testing.assert_array_equal(res.index, adata_large.obs_names)

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

        key = Key.varm.lineage_drivers(False)
        for lineage in ["0", "1"]:
            assert np.all(mc.adata.varm[key][f"{lineage}_corr"] >= -1.0)
            assert np.all(mc.adata.varm[key][f"{lineage}_corr"] <= 1.0)

            assert np.all(mc.adata.varm[key][f"{lineage}_qval"] >= 0)
            assert np.all(mc.adata.varm[key][f"{lineage}_qval"] <= 1.0)

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

    def test_compute_priming_clusters(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_macrostates(n_states=2)
        mc.set_terminal_states_from_macrostates()
        mc.compute_absorption_probabilities()

        cat = adata_large.obs["clusters"].cat.categories[0]
        deg1 = mc.compute_lineage_priming(
            method="kl_divergence", early_cells={"clusters": [cat]}
        )
        deg2 = mc.compute_lineage_priming(
            method="kl_divergence",
            early_cells=(adata_large.obs["clusters"] == cat).values,
        )

        assert_series_equal(deg1, deg2)
        # because passing it to a dataframe changes its name
        key = Key.obs.priming_degree(mc.backward)
        assert_series_equal(adata_large.obs[key], deg1, check_names=False)
        assert_series_equal(mc.priming_degree, deg1)

    def test_recompute_terminal_states_different_n_states(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_macrostates(n_states=2)
        mc.set_terminal_states_from_macrostates()

        key = Key.obsm.memberships(Key.obs.term_states(mc.backward))
        assert len(mc.terminal_states.cat.categories) == 2
        assert adata_large.obsm[key].shape == (adata_large.n_obs, 2)

        mc.set_terminal_states({"foo": adata_large.obs_names[:20]})

        assert len(mc.terminal_states.cat.categories) == 1
        assert mc.terminal_states_memberships is None
        assert key not in adata_large.obsm


class TestGPCCAIO:
    @pytest.mark.parametrize("deep", [False, True])
    def test_copy(
        self, adata_gpcca_fwd: Tuple[AnnData, cr.tl.estimators.GPCCA], deep: bool
    ):
        _, mc1 = adata_gpcca_fwd
        mc2 = mc1.copy(deep=deep)

        assert_estimators_equal(mc1, mc2, copy=True, deep=deep)

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
