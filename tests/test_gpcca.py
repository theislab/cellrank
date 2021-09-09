from typing import List, Tuple, Union, Optional, Sequence

import os
import pytest
from copy import deepcopy
from enum import Enum
from _helpers import assert_array_nan_equal, assert_estimators_equal

import cellrank as cr
from anndata import AnnData
from cellrank.tl import Lineage
from cellrank._key import Key
from cellrank.tl.kernels import VelocityKernel, ConnectivityKernel

import numpy as np
import pandas as pd
from scipy.sparse import issparse
from pandas.testing import assert_frame_equal, assert_series_equal


# fmt: off
class State(str, Enum):
    SCHUR = "schur"
    MACRO = "macro"
    TERM = "term"
    ABS = "abs"
    LIN = "lin"

    @property
    def prev(self) -> Optional["State"]:
        if self.value == "schur":
            return None
        if self.value == "macro":
            return State.SCHUR
        if self.value == "term":
            return State.MACRO
        if self.value == "abs":
            return State.TERM
        if self.value == "lin":
            return State.ABS
        return None

    @property
    def next(self) -> Optional["State"]:
        if self.value == "schur":
            return State.MACRO
        if self.value == "macro":
            return State.TERM
        if self.value == "term":
            return State.ABS
        if self.value == "abs":
            return State.LIN
        if self.value == "lin":
            return None
        return None

    @property
    def key(self) -> Optional[str]:
        bwd = False
        if self.value == "schur":
            return f"schur_decomposition_{Key.backward(bwd)}"
        if self.value == "macro":
            return Key.obs.macrostates(bwd)
        if self.value == "term":
            return Key.obs.term_states(bwd)
        if self.value == "abs":
            return Key.obsm.abs_probs(bwd)
        if self.value == "lin":
            return Key.varm.lineage_drivers(bwd)
        return None

    @property
    def attr_keys(self) -> Optional[Sequence[Tuple[str, str]]]:
        bwd = False
        if self.value == "schur":
            key1 = Key.uns.eigen(bwd)
            key2 = Key.obsm.schur_vectors(bwd)
            key3 = Key.uns.schur_matrix(bwd)
            return ("uns", key1), ("obsm", key2), ("uns", key3)
        if self.value == "macro":
            key1 = Key.obs.macrostates(bwd)
            key2 = Key.uns.colors(key1)
            key3 = Key.obsm.memberships(key1)
            return ("obs", key1), ("uns", key2), ("obsm", key3)
        if self.value == "term":
            key1 = Key.obs.term_states(bwd)
            key2 = Key.obs.probs(key1)
            key3 = Key.uns.colors(key1)
            return ("obs", key1), ("obs", key2), ("uns", key3)
        if self.value == "abs":
            key1 = Key.obsm.abs_probs(bwd)
            key2 = Key.obsm.abs_times(bwd)
            key3 = Key.obs.priming_degree(bwd)
            return ("obsm", key1), ("obsm", key2), ("obs", key3)
        if self.value == "lin":
            key1 = Key.varm.lineage_drivers(bwd)
            return ("varm", key1),
        return None

    @property
    def attrs(self) -> Optional[Sequence[Tuple[str, type]]]:
        if self.value == "schur":
            return ("_eigendecomposition", dict), ("_schur_vectors", np.ndarray), ("_schur_matrix", np.ndarray)
        if self.value == "macro":
            return (
                ("_macrostates", pd.Series), ("_macrostates_colors", np.ndarray),
                ("_macrostates_memberships", Lineage), ("_coarse_tmat", pd.DataFrame),
                ("_coarse_init_dist", pd.Series), ("_coarse_stat_dist", pd.Series)
            )
        if self.value == "term":
            return ("_term_states", pd.Series), ("_term_states_probs", pd.Series), ("_term_states_colors", np.ndarray)
        if self.value == "abs":
            return (
                ("_absorption_probabilities", Lineage), ("_absorption_times", pd.DataFrame),
                ("_priming_degree", pd.Series)
            )
        if self.value == "lin":
            return ("_lineage_drivers", pd.DataFrame),
        return None

# fmt: on


def shares_mem(x, y) -> bool:
    if isinstance(x, Lineage):
        return np.shares_memory(x.X, y.X)
    if issparse(x):
        return np.shares_memory(x.data, y.data)
    return np.shares_memory(x, y)


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
    # TODO(michalk8): assert series equal from pandas
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


def _fit_gpcca(adata, state: str, backward: bool = False) -> cr.estimators.GPCCA:
    vk = VelocityKernel(adata, backward=backward).compute_transition_matrix(
        softmax_scale=4
    )
    ck = ConnectivityKernel(adata, backward=backward).compute_transition_matrix()
    terminal_kernel = 0.8 * vk + 0.2 * ck

    mc = cr.tl.estimators.GPCCA(terminal_kernel)
    mc.compute_eigendecomposition()
    mc.compute_schur(n_components=10, method="krylov")
    if state == State.SCHUR:
        return mc
    mc.compute_macrostates(n_states=2)
    if state == State.MACRO:
        return mc
    mc.set_terminal_states_from_macrostates()
    if state == State.TERM:
        return mc
    mc.compute_absorption_probabilities(time_to_absorption="all")
    mc.compute_lineage_priming()
    if state == State.ABS:
        return mc
    mc.compute_lineage_drivers(use_raw=False, cluster_key="clusters")
    if state == State.LIN:
        return mc

    raise NotImplementedError(state)


def _assert_params(
    g: cr.tl.estimators.GPCCA,
    state: Optional[State],
    fwd: bool = False,
    init: bool = True,
) -> None:
    if state is None:
        return

    key = state.key
    if fwd:
        if init:
            assert isinstance(g.params[key], dict), state
        else:
            assert g.params.get(key, {}) == {}, state
        _assert_params(g, state.next, fwd=True, init=False)
    else:
        assert isinstance(g.params[key], dict), state
        _assert_params(g, state.prev, fwd=False)


def _assert_adata(
    adata: AnnData, state: Optional[State], fwd: bool = False, init: bool = True
) -> None:
    if state is None:
        return

    if fwd:
        if init:
            for attr, key in state.attr_keys:
                obj = getattr(adata, attr)
                assert key in obj, sorted(obj.keys())
        else:
            for attr, key in state.attr_keys:
                obj = getattr(adata, attr)
                assert key not in obj, sorted(obj.keys())
        _assert_adata(adata, state.next, fwd=True, init=False)
    else:
        for attr, key in state.attr_keys:
            obj = getattr(adata, attr)
            assert key in obj, sorted(obj.keys())
        _assert_adata(adata, state.prev, fwd=False)


def _assert_gpcca_attrs(
    g: cr.tl.estimators.GPCCA,
    state: Optional[State] = None,
    fwd: bool = False,
    init: bool = True,
) -> None:
    if state is None:
        return
    if fwd:
        if init:
            for attr, dtype in state.attrs:
                obj = getattr(g, attr)
                assert isinstance(obj, dtype)
        else:
            for attr, dtype in state.attrs:
                obj = getattr(g, attr)
                assert obj is None, attr
        _assert_gpcca_attrs(g, state.next, fwd=True, init=False)
    else:
        for attr, dtype in state.attrs:
            obj = getattr(g, attr)
            assert isinstance(obj, dtype), attr
        _assert_gpcca_attrs(g, state.prev, fwd=False)


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
        mc.compute_schur(n_components=1, method="krylov")

        assert mc.schur_vectors.shape[1] == 2

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

    def test_fit_predict(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        g = cr.tl.estimators.GPCCA(terminal_kernel)
        tmp = g.fit(n_states=2)

        assert tmp is g
        _assert_params(g, State.MACRO, fwd=False)
        _assert_params(g, State.MACRO, fwd=True)

        res = g.predict()
        assert res is None
        _assert_params(g, State.TERM, fwd=False)
        _assert_params(g, State.TERM, fwd=True)

    @pytest.mark.parametrize("state", list(State))
    def test_params(self, adata_large: AnnData, state: State):
        g = _fit_gpcca(adata_large, state)
        _assert_params(g, state, fwd=False)
        _assert_params(g, state, fwd=True)


class TestGPCCASerialization:
    @pytest.mark.parametrize("state", list(State))
    def test_to_adata(self, adata_large: AnnData, state: State):
        g = _fit_gpcca(adata_large, state)
        adata = g.to_adata()
        _assert_adata(adata, state, fwd=False)
        _assert_adata(adata, state, fwd=True)

    @pytest.mark.parametrize("copy", [False, True])
    @pytest.mark.parametrize("keep", ["X", ("obs", "obsm"), ("layers")])
    def test_to_adata_keep(
        self, adata_large: AnnData, keep: Union[str, Sequence[str]], copy: bool
    ):
        g = _fit_gpcca(adata_large, State.MACRO)
        adata = g.to_adata(keep=keep, copy=copy)
        if "X" in keep:
            res = shares_mem(adata.X, g.adata.X)
            assert not res if copy else res
        if "obs" in keep:
            # we don't save macrostates in the underlying object, only during serialization
            key = Key.obs.macrostates(g.backward)
            cols = [c for c in adata.obs.columns if c != key]
            obs = adata.obs[cols]
            np.testing.assert_array_equal(obs.shape, g.adata.obs.shape)
            assert_frame_equal(obs, g.adata.obs[cols])
        if "obsm" in keep:
            for key in g.adata.obsm.keys():
                res = shares_mem(adata.obsm[key], g.adata.obsm[key])
                assert not res if copy else res, key
        if "layers" in keep:
            for key in g.adata.layers.keys():
                res = shares_mem(adata.layers[key], g.adata.layers[key])
                assert not res if copy else res, key

    @pytest.mark.parametrize("copy", [False, True])
    def test_to_adata_copy_finegrained(self, adata_large: AnnData, copy):
        g = _fit_gpcca(adata_large, State.MACRO)
        adata = g.to_adata(keep=["obsp", "varm", "var"], copy=["obsp"])

        assert_frame_equal(adata.var, g.adata.var)
        for key in g.adata.obsp.keys():
            assert not shares_mem(adata.obsp[key], g.adata.obsp[key]), key
        for key in g.adata.varm.keys():
            assert shares_mem(adata.varm[key], g.adata.varm[key]), key

    @pytest.mark.parametrize("state", list(State))
    def test_from_adata(self, adata_large: AnnData, state: State):
        g1 = _fit_gpcca(adata_large, state)
        adata = g1.to_adata()
        g2 = cr.tl.estimators.GPCCA.from_adata(adata, obsp_key="T_fwd")

        _assert_adata(g2.adata, state, fwd=False)
        _assert_adata(g2.adata, state, fwd=True)
        assert_estimators_equal(g1, g2, from_adata=True)

    @pytest.mark.parametrize("state", list(State))
    def test_from_adata_incomplete(self, adata_large: AnnData, state: State):
        if state in (State.SCHUR, State.MACRO):
            pytest.skip("See the reintroduce TODO in GPCCA")
        g_orig = _fit_gpcca(adata_large, State.LIN)
        adata = g_orig.to_adata()

        for attr, key in state.attr_keys:
            obj = getattr(adata, attr)
            del obj[key]

        g = cr.tl.estimators.GPCCA.from_adata(adata, obsp_key="T_fwd")
        _assert_gpcca_attrs(g, state, fwd=True, init=False)
        _assert_gpcca_attrs(g, state.prev, fwd=False)
        _assert_adata(g._shadow_adata, state, fwd=True, init=False)
        _assert_adata(g._shadow_adata, state.prev, fwd=False)

    @pytest.mark.parametrize("bwd", [False, True])
    def test_reconstruct_lineage(self, adata_large: AnnData, bwd: bool):
        g_orig = _fit_gpcca(adata_large, State.LIN, backward=bwd)
        adata = g_orig.to_adata()
        bdata = adata.copy()
        keys = [
            Key.obsm.memberships(Key.obs.macrostates(bwd)),
            Key.obsm.memberships(Key.obs.term_states(bwd)),
            Key.obsm.abs_probs(bwd),
        ]
        attrs = [
            "macrostates_memberships",
            "terminal_states_memberships",
            "absorption_probabilities",
        ]

        for key in keys:
            assert isinstance(adata.obsm[key], Lineage)
            adata.obsm[key] = np.array(adata.obsm[key])
            assert isinstance(adata.obsm[key], np.ndarray)

        g = cr.tl.estimators.GPCCA.from_adata(adata, obsp_key=Key.uns.kernel(bwd))
        for attr, key in zip(attrs, keys):
            obj = getattr(g, attr)
            assert isinstance(obj, Lineage)
            np.testing.assert_allclose(obj.X, bdata.obsm[key].X)
            np.testing.assert_array_equal(obj.names, bdata.obsm[key].names)
            np.testing.assert_array_equal(obj.colors, bdata.obsm[key].colors)

        assert_estimators_equal(g_orig, g, from_adata=True)


class TestGPCCAIO:
    @pytest.mark.parametrize("deep", [False, True])
    def test_copy(
        self, adata_gpcca_fwd: Tuple[AnnData, cr.tl.estimators.GPCCA], deep: bool
    ):
        _, mc1 = adata_gpcca_fwd
        mc2 = mc1.copy(deep=deep)

        assert_estimators_equal(mc1, mc2, copy=True, deep=deep)

    def test_write_ext(
        self, adata_gpcca_fwd: Tuple[AnnData, cr.tl.estimators.GPCCA], tmpdir
    ):
        _, mc = adata_gpcca_fwd

        fname = "foo"
        mc.write(os.path.join(tmpdir, fname), ext="bar")

        assert os.path.isfile(os.path.join(tmpdir, f"foo.bar"))

    def test_write_no_ext(
        self,
        adata_gpcca_fwd: Tuple[AnnData, cr.tl.estimators.GPCCA],
        tmpdir,
    ):
        _, mc = adata_gpcca_fwd
        fname = "foo"
        mc.write(os.path.join(tmpdir, fname), ext=None)

        assert os.path.isfile(os.path.join(tmpdir, f"foo"))

    def test_write_ext_with_dot(
        self,
        adata_gpcca_fwd: Tuple[AnnData, cr.tl.estimators.GPCCA],
        tmpdir,
    ):
        _, mc = adata_gpcca_fwd

        fname = "foo"
        mc.write(os.path.join(tmpdir, fname), ext=".bar")

        assert os.path.isfile(os.path.join(tmpdir, f"foo.bar"))

    def test_read(
        self, adata_gpcca_fwd: Tuple[AnnData, cr.tl.estimators.GPCCA], tmpdir
    ):
        _, mc1 = adata_gpcca_fwd

        mc1.write(os.path.join(tmpdir, "foo"))
        mc2 = cr.tl.estimators.GPCCA.read(os.path.join(tmpdir, "foo.pickle"))

        assert_estimators_equal(mc1, mc2)

    @pytest.mark.parametrize("copy", [False, True])
    def test_write_no_adata(
        self,
        adata_gpcca_fwd: Tuple[AnnData, cr.tl.estimators.GPCCA],
        copy: bool,
        tmpdir,
    ):
        adata, mc1 = adata_gpcca_fwd

        mc1.write(os.path.join(tmpdir, "foo"), write_adata=False)
        mc2 = cr.tl.estimators.GPCCA.read(
            os.path.join(tmpdir, "foo.pickle"), adata=adata, copy=copy
        )

        if copy:
            assert adata is not mc2.adata
        else:
            assert adata is mc2.adata
        assert_estimators_equal(mc1, mc2)

    def test_write_no_adata_read_none_supplied(
        self, adata_gpcca_fwd: Tuple[AnnData, cr.tl.estimators.GPCCA], tmpdir
    ):
        _, mc1 = adata_gpcca_fwd

        mc1.write(os.path.join(tmpdir, "foo"), write_adata=False)
        with pytest.raises(TypeError, match="This object was saved without"):
            _ = cr.tl.estimators.GPCCA.read(
                os.path.join(tmpdir, "foo.pickle"), adata=None
            )

    def test_write_no_adata_read_wrong_length(
        self, adata_gpcca_fwd: Tuple[AnnData, cr.tl.estimators.GPCCA], tmpdir
    ):
        _, mc1 = adata_gpcca_fwd
        adata = AnnData(np.random.normal(size=(len(mc1) + 1, 1)))

        mc1.write(os.path.join(tmpdir, "foo"), write_adata=False)
        with pytest.raises(ValueError, match="Expected `adata` to be of length"):
            _ = cr.tl.estimators.GPCCA.read(
                os.path.join(tmpdir, "foo.pickle"), adata=adata
            )
