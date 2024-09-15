import copy
import enum
import os
from collections.abc import Sequence
from typing import Optional, Union

import pytest
from _helpers import assert_array_nan_equal, assert_estimators_equal

import numpy as np
import pandas as pd
import scipy.sparse as sp
from pandas.testing import assert_frame_equal, assert_series_equal

from anndata import AnnData

import cellrank as cr
from cellrank._utils import Lineage
from cellrank._utils._key import Key
from cellrank.kernels import ConnectivityKernel, VelocityKernel


# fmt: off
class State(str, enum.Enum):
    SCHUR = "schur"
    MACRO = "macro"
    TERM = "term"
    FATE_PROBS = "fate"
    DRIVERS = "lin"

    @property
    def prev(self) -> Optional["State"]:
        if self.value == "schur":
            return None
        if self.value == "macro":
            return State.SCHUR
        if self.value == "term":
            return State.MACRO
        if self.value == "fate":
            return State.TERM
        if self.value == "lin":
            return State.FATE_PROBS
        return None

    @property
    def next(self) -> Optional["State"]:
        if self.value == "schur":
            return State.MACRO
        if self.value == "macro":
            return State.TERM
        if self.value == "term":
            return State.FATE_PROBS
        if self.value == "fate":
            return State.DRIVERS
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
            return Key.obs.term_states(False, bwd=bwd)
        if self.value == "fate":
            return Key.obsm.fate_probs(bwd)
        if self.value == "lin":
            return Key.varm.lineage_drivers(bwd)
        return None

    @property
    def attr_keys(self) -> Optional[Sequence[tuple[str, str]]]:
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
            key1 = Key.obs.term_states(False, bwd=bwd)
            key2 = Key.obs.probs(key1)
            key3 = Key.uns.colors(key1)
            return ("obs", key1), ("obs", key2), ("uns", key3)
        if self.value == "fate":
            key1 = Key.obsm.fate_probs(bwd)
            key2 = Key.obsm.abs_times(bwd)
            key3 = Key.obs.priming_degree(bwd)
            return ("obsm", key1), ("obsm", key2), ("obs", key3)
        if self.value == "lin":
            key1 = Key.varm.lineage_drivers(bwd)
            return ("varm", key1),
        return None

    @property
    def attrs(self) -> Optional[Sequence[tuple[str, type]]]:
        if self.value == "schur":
            return ("_eigendecomposition", dict), ("_schur_vectors", np.ndarray), ("_schur_matrix", np.ndarray)
        if self.value == "macro":
            return (
                ("macrostates", pd.Series), ("macrostates_memberships", Lineage), ("_coarse_tmat", pd.DataFrame),
                ("_coarse_init_dist", pd.Series), ("_coarse_stat_dist", pd.Series)
            )
        if self.value == "term":
            return ("terminal_states", pd.Series), ("terminal_states_probabilities", pd.Series)
        if self.value == "fate":
            return (
                ("_fate_probabilities", Lineage), ("_absorption_times", pd.DataFrame),
                ("_priming_degree", pd.Series)
            )
        if self.value == "lin":
            return ("_lineage_drivers", pd.DataFrame),
        return None

# fmt: on


def shares_mem(x, y) -> bool:
    if isinstance(x, Lineage):
        return np.shares_memory(x.X, y.X)
    if sp.issparse(x):
        return np.shares_memory(x.data, y.data)
    return np.shares_memory(x, y)


def _check_eigdecomposition(mc: cr.estimators.GPCCA) -> None:
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


def _check_compute_macro(mc: cr.estimators.GPCCA) -> None:
    assert isinstance(mc.macrostates, pd.Series)
    assert len(mc._macrostates.colors) == len(mc.macrostates.cat.categories)

    if "stationary_dist" in mc.eigendecomposition:  # one state
        assert isinstance(mc.macrostates_memberships, cr.Lineage)
        assert mc.macrostates_memberships.shape[1] == 1
        np.testing.assert_allclose(mc.macrostates_memberships.X.sum(), 1.0)

        assert mc.schur_matrix is None
        assert mc.schur_vectors is None
        assert mc.coarse_initial_distribution is None
        assert mc.coarse_stationary_distribution is None
        assert mc.coarse_T is None
    else:
        assert isinstance(mc.macrostates_memberships, cr.Lineage)
        if mc.macrostates_memberships.shape[1] > 1:
            np.testing.assert_allclose(mc.macrostates_memberships.X.sum(1), 1.0)

        assert isinstance(mc.schur_matrix, np.ndarray)
        assert isinstance(mc.schur_vectors, np.ndarray)
        assert isinstance(mc.coarse_initial_distribution, pd.Series)
        assert isinstance(mc.coarse_T, pd.DataFrame)
        np.testing.assert_array_equal(mc.coarse_T.index, mc.coarse_T.columns)
        np.testing.assert_array_equal(mc.coarse_T.index, mc.coarse_stationary_distribution.index)
        if mc.coarse_stationary_distribution is not None:
            assert isinstance(mc.coarse_stationary_distribution, pd.Series)
            np.testing.assert_array_equal(mc.coarse_T.index, mc.coarse_stationary_distribution.index)


def _check_renaming_no_write_terminal(mc: cr.estimators.GPCCA) -> None:
    assert mc.terminal_states is None
    assert mc.terminal_states_probabilities is None
    assert mc.terminal_states_memberships is None

    key = Key.obs.term_states(mc.backward)
    assert key not in mc.adata.obs
    assert Key.obs.probs(key) not in mc.adata.obs
    assert Key.uns.colors(key) not in mc.adata.uns


def _check_fate_probs(mc: cr.estimators.GPCCA) -> None:
    # fmt: off
    # macrostates
    assert isinstance(mc.macrostates, pd.Series)
    assert isinstance(mc.macrostates_memberships, cr.Lineage)
    np.testing.assert_array_equal(mc._macrostates.colors, mc.macrostates_memberships.colors)

    # term states
    key = Key.obs.term_states(mc.backward)
    assert isinstance(mc.terminal_states, pd.Series)
    # TODO(michalk8): assert series equal from pandas
    assert_array_nan_equal(mc.adata.obs[key], mc.terminal_states)
    np.testing.assert_array_equal(mc.adata.uns[Key.uns.colors(key)], mc.fate_probabilities.colors)
    np.testing.assert_array_equal(mc.adata.uns[Key.uns.colors(key)], mc._term_states.colors)
    assert isinstance(mc.terminal_states_probabilities, pd.Series)
    np.testing.assert_array_equal(mc.adata.obs[Key.obs.probs(key)], mc.terminal_states_probabilities)

    # fate probs
    key = Key.obsm.fate_probs(mc.backward)
    assert isinstance(mc.fate_probabilities, cr.Lineage)
    np.testing.assert_array_almost_equal(mc.fate_probabilities.sum(1), 1.0)
    assert isinstance(mc.adata.obsm[key], cr.Lineage)
    np.testing.assert_array_equal(mc.adata.obsm[key], mc.fate_probabilities.X)

    # priming
    key = Key.obs.priming_degree(mc.backward)
    assert isinstance(mc.priming_degree, pd.Series)
    np.testing.assert_array_equal(mc.adata.obs[key], mc.priming_degree)
    # fmt: on


def _fit_gpcca(adata, state: str, backward: bool = False) -> cr.estimators.GPCCA:
    vk = VelocityKernel(adata, backward=backward).compute_transition_matrix(softmax_scale=4)
    ck = ConnectivityKernel(adata).compute_transition_matrix()
    terminal_kernel = 0.8 * vk + 0.2 * ck

    mc = cr.estimators.GPCCA(terminal_kernel)
    mc.compute_eigendecomposition()
    mc.compute_schur(n_components=10, method="krylov")
    if state == State.SCHUR:
        return mc
    mc.compute_macrostates(n_states=2)
    if state == State.MACRO:
        return mc
    mc.set_terminal_states()
    if state == State.TERM:
        return mc
    mc.compute_fate_probabilities()
    mc.compute_absorption_times()
    mc.compute_lineage_priming()
    if state == State.FATE_PROBS:
        return mc
    mc.compute_lineage_drivers(use_raw=False, cluster_key="clusters")
    if state == State.DRIVERS:
        return mc

    raise NotImplementedError(state)


def _assert_params(
    g: cr.estimators.GPCCA,
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


def _assert_adata(adata: AnnData, state: Optional[State], fwd: bool = False, init: bool = True) -> None:
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
                assert key not in obj, (state, attr, key)
        _assert_adata(adata, state.next, fwd=True, init=False)
    else:
        for attr, key in state.attr_keys:
            obj = getattr(adata, attr)
            assert key in obj, (state, attr, key)
        _assert_adata(adata, state.prev, fwd=False)


def _assert_gpcca_attrs(
    g: cr.estimators.GPCCA,
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
            for attr, _ in state.attrs:
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

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_eigendecomposition(k=2, only_evals=True)

        _check_eigdecomposition(mc)

    def test_compute_schur_invalid_n_comps(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=1, method="krylov")

        assert mc.schur_vectors.shape[1] == 2

    def test_compute_schur_invalid_method(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        with pytest.raises(ValueError, match=r"Invalid method"):
            mc.compute_schur(method="foobar")

    def test_compute_schur_invalid_eig_sort(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        with pytest.raises(ValueError, match=r".* sorting criterion"):
            mc.compute_schur(which="foobar", method="krylov")

    def test_compute_schur_write_eigvals(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        _check_eigdecomposition(mc)

    def test_compute_schur_write_eigvals_similar_to_orig_eigdecomp(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_eigendecomposition(k=10, only_evals=True)

        _check_eigdecomposition(mc)
        orig_ed = copy.deepcopy(mc.eigendecomposition)

        mc._eigendecomposition = None
        mc.compute_schur(n_components=10, method="krylov")

        _check_eigdecomposition(mc)
        schur_ed = mc.eigendecomposition

        assert orig_ed.keys() == schur_ed.keys()
        assert orig_ed["eigengap"] == schur_ed["eigengap"]
        n = min(orig_ed["params"]["k"], schur_ed["params"]["k"])
        np.testing.assert_array_almost_equal(orig_ed["D"].real[:n], schur_ed["D"].real[:n])
        np.testing.assert_array_almost_equal(
            np.abs(orig_ed["D"].imag[:n]), np.abs(schur_ed["D"].imag[:n])
        )  # complex conj.

    def test_compute_macrostates_no_eig(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        with pytest.raises(RuntimeError, match=r"Compute eigendecomposition"):
            mc.compute_macrostates(n_states=None)

    def test_compute_macrostates_1_state_no_eig(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_macrostates(n_states=1)

    def test_compute_macro_none_states(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_eigendecomposition(only_evals=True)
        mc.compute_macrostates(n_states=None)

        _check_compute_macro(mc)

    def test_compute_macrostates_schur(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_macrostates(n_states=2)

        _check_compute_macro(mc)

    def test_compute_macrostates_min_chi_too_low_min(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_macrostates(n_states=[1, 1])

        _check_compute_macro(mc)

    def test_compute_macrostates_min_chi_inverted_range(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_macrostates(n_states=[4, 2])

        _check_compute_macro(mc)

    def test_compute_macrostates_min_chi_range_same_values(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_macrostates(n_states=[2, 2])

        _check_compute_macro(mc)

    def test_compute_macrostates_min_chi_normal_run(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_macrostates(n_states=[2, 4])

        _check_compute_macro(mc)

    def test_compute_macro_invalid_cluster_key(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        with pytest.raises(KeyError, match=r"Unable to find clusters"):
            mc.compute_macrostates(n_states=2, cluster_key="foobar")

    def test_compute_macro_cache(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=11, method="krylov")

        mc.compute_macrostates(n_states=2)

        assert mc.schur_vectors.shape[1] == 11
        np.testing.assert_array_equal(mc.schur_matrix.shape, [11, 11])

    def test_set_initial_states_from_forward(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large, backward=False).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=2, n_cells=5)
        mc.set_initial_states("0")

        key = Key.obs.term_states(mc.backward, bwd=True)

        assert key in mc.adata.obs
        np.testing.assert_array_equal(mc.adata.obs[key].cat.categories, ["0"])
        assert Key.obs.probs(key) in mc.adata.obs
        assert Key.uns.colors(key) in mc.adata.uns

    def test_compute_initial_states_from_forward_no_macro(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large, backward=False).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        with pytest.raises(RuntimeError, match=r"Compute macro"):
            mc.predict_initial_states(n_states=1)

    def test_compute_initial_states_from_forward(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=4, n_cells=5)
        mc.predict_initial_states(n_states=3)

        assert mc.terminal_states is None
        assert len(mc.initial_states.cat.categories) == 3
        assert mc.initial_states_memberships.shape == (adata_large.n_obs, 3)

    def test_compute_initial_states_from_forward_normal_run(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large, backward=False).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=2, n_cells=5)

        mc.predict_initial_states(n_states=1)

        key = Key.obs.term_states(mc.backward, bwd=True)

        assert key in mc.adata.obs
        np.testing.assert_array_equal(mc.adata.obs[key].cat.categories, ["0"])
        assert Key.obs.probs(key) in mc.adata.obs
        assert Key.uns.colors(key) in mc.adata.uns

    def test_set_terminal_states_from_macrostates(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=2, n_cells=5)
        mc.set_terminal_states()
        mc.compute_fate_probabilities()
        mc.compute_lineage_priming()

        _check_fate_probs(mc)

    def test_set_terminal_states_from_macrostates_non_positive_cells(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=2, n_cells=None)
        with pytest.raises(ValueError, match=r"Expected .* to be positive"):
            mc.set_terminal_states(n_cells=0)

    def test_set_terminal_states_from_macrostates_invalid_name(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=2)
        with pytest.raises(KeyError, match=r"Invalid lineage name"):
            mc.set_terminal_states(states=["foobar"])

    @pytest.mark.parametrize("values", ["Astrocytes", ["Astrocytes", "OPC"]])
    def test_set_terminal_states_clusters(self, adata_large: AnnData, values: Union[str, list[str]]):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        to_remove = list(
            set(adata_large.obs["clusters"].cat.categories) - ({values} if isinstance(values, str) else set(values))
        )
        expected = adata_large.obs["clusters"].cat.remove_categories(to_remove)

        mc = cr.estimators.GPCCA(terminal_kernel)

        mc.set_terminal_states({"clusters": values})
        pd.testing.assert_series_equal(expected, mc.terminal_states, check_category_order=False, check_names=False)

    def test_compute_terminal_states_invalid_method(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=2)
        with pytest.raises(ValueError, match=r"Invalid option"):
            mc.predict(method="foobar")

    def test_compute_terminal_states_non_positive_cells(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=2)
        with pytest.raises(ValueError, match=r".* to be positive"):
            mc.predict(n_cells=0)

    def test_compute_terminal_states_eigengap(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=2)
        mc.predict(n_cells=5, method="eigengap")
        mc.compute_fate_probabilities()
        mc.compute_lineage_priming()

        _check_fate_probs(mc)

    def test_compute_terminal_states_n_main_states(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=2)
        mc.predict(n_cells=5, method="top_n", n_states=1)
        mc.compute_fate_probabilities()
        mc.compute_lineage_priming()

        _check_fate_probs(mc)

    def test_compute_terminal_states_stability(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck
        thresh = 0.5

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=5)
        mc.predict(n_cells=5, method="stability", stability_threshold=thresh)
        mc.compute_fate_probabilities()
        mc.compute_lineage_priming()

        coarse_T = mc.coarse_T
        self_probs = pd.Series(np.diag(coarse_T), index=coarse_T.columns)
        names = self_probs[self_probs.values >= thresh].index

        np.testing.assert_array_equal(set(names), set(mc.terminal_states.cat.categories))
        _check_fate_probs(mc)

    def test_compute_terminal_states_no_selected(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=2)
        with pytest.raises(ValueError, match=r"No macrostates"):
            mc.predict(n_cells=5, method="stability", stability_threshold=42)

    def test_compute_terminal_states_too_many_cells(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=2)
        with pytest.raises(ValueError, match=r".* decrease this to at most"):
            mc.predict(n_cells=4200)

    def test_compute_terminal_states_default(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")

        mc.compute_macrostates(n_states=2)
        mc.predict(n_cells=5)
        mc.compute_fate_probabilities()
        mc.compute_lineage_priming()

        _check_fate_probs(mc)

    def test_compute_lineage_drivers_invalid_lineages(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_macrostates(n_states=2)
        mc.set_terminal_states()
        mc.compute_fate_probabilities()

        with pytest.raises(KeyError, match=r"Invalid lineage name"):
            mc.compute_lineage_drivers(use_raw=False, lineages=["foo"])

    def test_compute_lineage_drivers_invalid_clusters(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_macrostates(n_states=2)
        mc.set_terminal_states()
        mc.compute_fate_probabilities()

        with pytest.raises(KeyError, match=r"Clusters .* not found"):
            mc.compute_lineage_drivers(use_raw=False, cluster_key="clusters", clusters=["foo"])

    def test_compute_lineage_drivers_normal_run(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_macrostates(n_states=2)
        mc.set_terminal_states()
        mc.compute_fate_probabilities()
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

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_macrostates(n_states=2)
        mc.set_terminal_states()
        mc.compute_fate_probabilities()

        with pytest.raises(RuntimeError, match=r".*lineage_drivers"):
            mc.plot_lineage_drivers("0")

    def test_plot_lineage_drivers_invalid_name(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_macrostates(n_states=2)
        mc.set_terminal_states()
        mc.compute_fate_probabilities()
        mc.compute_lineage_drivers(use_raw=False, cluster_key="clusters")

        with pytest.raises(KeyError, match=r"Lineage .* not found"):
            mc.plot_lineage_drivers("foo", use_raw=False)

    def test_plot_lineage_drivers_invalid_n_genes(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_macrostates(n_states=2)
        mc.set_terminal_states()
        mc.compute_fate_probabilities()
        mc.compute_lineage_drivers(use_raw=False, cluster_key="clusters")

        with pytest.raises(ValueError, match=r".* to be positive"):
            mc.plot_lineage_drivers("0", use_raw=False, n_genes=0)

    def test_plot_lineage_drivers_normal_run(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_macrostates(n_states=2)
        mc.set_terminal_states()
        mc.compute_fate_probabilities()
        mc.compute_lineage_drivers(use_raw=False, cluster_key="clusters")

        mc.plot_lineage_drivers("0", use_raw=False)

    def test_tsi(self, adata_large: AnnData):
        groundtruth_adata = adata_large.uns["tsi"].copy()

        vk = VelocityKernel(adata_large).compute_transition_matrix()
        estimator = cr.estimators.GPCCA(vk)
        estimator.compute_schur(n_components=5)

        terminal_states = ["Neuroblast", "Astrocyte", "Granule mature"]
        cluster_key = "clusters"
        tsi_score = estimator.tsi(n_macrostates=3, terminal_states=terminal_states, cluster_key=cluster_key, n_cells=10)

        np.testing.assert_almost_equal(tsi_score, groundtruth_adata.uns["score"])
        assert isinstance(estimator._tsi.uns["terminal_states"], list)
        assert len(estimator._tsi.uns["terminal_states"]) == len(groundtruth_adata.uns["terminal_states"])
        assert (estimator._tsi.uns["terminal_states"] == groundtruth_adata.uns["terminal_states"]).all()
        assert estimator._tsi.uns["cluster_key"] == groundtruth_adata.uns["cluster_key"]
        pd.testing.assert_frame_equal(estimator._tsi.to_df(), groundtruth_adata.to_df())

    def test_compute_priming_clusters(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_macrostates(n_states=2)
        mc.set_terminal_states()
        mc.compute_fate_probabilities()

        cat = adata_large.obs["clusters"].cat.categories[0]
        deg1 = mc.compute_lineage_priming(method="kl_divergence", early_cells={"clusters": [cat]})
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

        mc = cr.estimators.GPCCA(terminal_kernel)
        mc.compute_schur(n_components=10, method="krylov")
        mc.compute_macrostates(n_states=2)
        mc.set_terminal_states()

        key = Key.obsm.memberships(Key.obs.term_states(estim_bwd=mc.backward))
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

        g = cr.estimators.GPCCA(terminal_kernel)
        tmp = g.fit(n_states=2)

        assert tmp is g
        _assert_params(g, State.MACRO, fwd=False)
        _assert_params(g, State.MACRO, fwd=True)

        res = g.predict()
        assert res is g
        _assert_params(g, State.TERM, fwd=False)
        _assert_params(g, State.TERM, fwd=True)

    @pytest.mark.parametrize("state", list(State))
    def test_params(self, adata_large: AnnData, state: State):
        g = _fit_gpcca(adata_large, state)
        _assert_params(g, state, fwd=False)
        _assert_params(g, state, fwd=True)

    def test_drivers_constant_gene_qvalues(self, adata_large: AnnData):
        ix = 0
        gene = adata_large.var_names[ix]
        adata_large.X[:, ix] = 0
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        g = cr.estimators.GPCCA(ck)
        g.fit(n_states=2)
        g.predict()
        g.compute_fate_probabilities()

        drivers = g.compute_lineage_drivers(use_raw=False)

        np.testing.assert_array_equal(drivers.loc[gene], np.nan)
        assert np.asarray(drivers.iloc[drivers.index != gene].isnull()).sum() == 0

    @pytest.mark.parametrize("keys", ["0", ["0", "1"], ["0", "1, 2"]])
    def test_compute_time_to_absorption(self, adata_large: AnnData, keys: Union[str, Sequence[str]]):
        n_states = 3
        n_expected = 1 if isinstance(keys, str) else len(keys)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        g = cr.estimators.GPCCA(ck)
        g.fit(n_states=n_states)
        g.set_terminal_states()
        g.compute_fate_probabilities()
        g.compute_absorption_times(keys=keys)

        np.testing.assert_array_equal(g.absorption_times.shape, (len(g), n_expected))


class TestGPCCASerialization:
    @pytest.mark.parametrize("state", list(State))
    def test_to_adata(self, adata_large: AnnData, state: State):
        g = _fit_gpcca(adata_large, state)
        adata = g.to_adata()
        _assert_adata(adata, state, fwd=False)
        _assert_adata(adata, state, fwd=True)

    @pytest.mark.parametrize("copy", [False, True])
    @pytest.mark.parametrize("keep", ["X", ("obs", "obsm"), ("layers",)])
    def test_to_adata_keep(self, adata_large: AnnData, keep: Union[str, Sequence[str]], copy: bool):
        g = _fit_gpcca(adata_large, State.MACRO)
        adata = g.to_adata(keep=keep, copy=copy)
        if "X" in keep:
            res = shares_mem(adata.X, g.adata.X)
            assert not res if copy else res
        if "obs" in keep:
            columns = adata.obs.columns
            np.testing.assert_array_equal(adata.obs.shape, g.adata.obs.shape)
            # reorder columns
            assert_frame_equal(adata.obs, g.adata.obs[columns])
        if "obsm" in keep:
            for key in g.adata.obsm:
                res = shares_mem(adata.obsm[key], g.adata.obsm[key])
                assert not res if copy else res, key
        if "layers" in keep:
            for key in g.adata.layers:
                res = shares_mem(adata.layers[key], g.adata.layers[key])
                assert not res if copy else res, key

    @pytest.mark.parametrize("copy", [False, True])
    def test_to_adata_copy_finegrained(self, adata_large: AnnData, copy):
        g = _fit_gpcca(adata_large, State.MACRO)
        adata = g.to_adata(keep=["obsp", "varm", "var"], copy=["obsp"])

        assert_frame_equal(adata.var, g.adata.var)
        for key in g.adata.obsp:
            assert not shares_mem(adata.obsp[key], g.adata.obsp[key]), key
        for key in g.adata.varm:
            assert shares_mem(adata.varm[key], g.adata.varm[key]), key

    @pytest.mark.parametrize("state", list(State))
    def test_from_adata(self, adata_large: AnnData, state: State):
        g1 = _fit_gpcca(adata_large, state)
        adata = g1.to_adata()
        g2 = cr.estimators.GPCCA.from_adata(adata, obsp_key="T_fwd")

        _assert_adata(g2.adata, state, fwd=False)
        _assert_adata(g2.adata, state, fwd=True)
        assert_estimators_equal(g1, g2, from_adata=True)

    # TODO(michalk8): parametrize by bwd, needs attr_keys modification
    @pytest.mark.parametrize("state", list(State))
    def test_from_adata_incomplete(self, adata_large: AnnData, state: State):
        if state == State.SCHUR:
            pytest.xfail("Schur decomposition is not needed.")
        g_orig = _fit_gpcca(adata_large, State.DRIVERS)
        adata = g_orig.to_adata()

        for attr, key in state.attr_keys:
            obj = getattr(adata, attr)
            del obj[key]

        g = cr.estimators.GPCCA.from_adata(adata, obsp_key="T_fwd")
        _assert_gpcca_attrs(g, state, fwd=True, init=False)
        _assert_gpcca_attrs(g, state.prev, fwd=False)
        _assert_adata(g._shadow_adata, state, fwd=True, init=False)
        _assert_adata(g._shadow_adata, state.prev, fwd=False)

    @pytest.mark.parametrize("bwd", [False, True])
    def test_reconstruct_lineage(self, adata_large: AnnData, bwd: bool):
        g_orig = _fit_gpcca(adata_large, State.DRIVERS, backward=bwd)
        adata = g_orig.to_adata()
        bdata = adata.copy()
        keys = [
            Key.obsm.memberships(Key.obs.macrostates(bwd)),
            Key.obsm.memberships(Key.obs.term_states(bwd)),
            Key.obsm.fate_probs(bwd),
        ]
        attrs = [
            "macrostates_memberships",
            "terminal_states_memberships",
            "fate_probabilities",
        ]

        for key in keys:
            assert isinstance(adata.obsm[key], Lineage)
            adata.obsm[key] = np.array(adata.obsm[key], copy=True)
            assert isinstance(adata.obsm[key], np.ndarray)

        g = cr.estimators.GPCCA.from_adata(adata, obsp_key=Key.uns.kernel(bwd))
        for attr, key in zip(attrs, keys):
            obj = getattr(g, attr)
            assert isinstance(obj, Lineage), attr
            np.testing.assert_allclose(obj.X, bdata.obsm[key].X)
            np.testing.assert_array_equal(obj.names, bdata.obsm[key].names)
            np.testing.assert_array_equal(obj.colors, bdata.obsm[key].colors)

        assert_estimators_equal(g_orig, g, from_adata=True)

    def test_overlapping_states_predict(self, adata_large: AnnData):
        g = _fit_gpcca(adata_large, State.MACRO, backward=False)
        n_states = g.macrostates_memberships.shape[1]

        g.predict_initial_states(n_states=n_states)
        with pytest.raises(ValueError, match=r"Found \`\d+\` overlapping"):
            g.predict_terminal_states(method="top_n", n_states=n_states)

        g.predict_terminal_states(method="top_n", n_states=n_states, allow_overlap=True)

        pd.testing.assert_series_equal(g.initial_states, g.terminal_states)

    def test_overlapping_states_from_macro(self, adata_large: AnnData):
        g = _fit_gpcca(adata_large, State.TERM, backward=False)
        initial_states = g.terminal_states.cat.categories

        with pytest.raises(ValueError, match=r"Found \`\d+\` overlapping"):
            g.set_initial_states(initial_states)
        assert g.initial_states is None

        g.set_initial_states(initial_states, allow_overlap=True)

        pd.testing.assert_series_equal(g.initial_states, g.terminal_states)

    def test_overlapping_states_explicit(self, adata_large: AnnData):
        g = _fit_gpcca(adata_large, State.TERM, backward=False)
        initial_states = {"overlapping_initial_state": g.adata.obs_names[:20]}

        with pytest.raises(ValueError, match=r"Found \`\d+\` overlapping"):
            g.set_initial_states(initial_states)
        assert g.initial_states is None

        g.set_initial_states(initial_states, allow_overlap=True)
        assert g.initial_states is not None

    @pytest.mark.parametrize("initial", [False, True])
    def test_rename_states_normal_run(self, adata_large: AnnData, initial: bool):
        g = _fit_gpcca(adata_large, State.MACRO)
        expected = ["foo", "bar"]

        if initial:
            g = g.set_initial_states(states=None)
            g = g.rename_initial_states({"0": expected[0], "1": expected[1]})
            np.testing.assert_array_equal(g.initial_states.cat.categories, expected)
            np.testing.assert_array_equal(g.initial_states_memberships.names, expected)
        else:
            g = g.set_terminal_states(states=None)
            g = g.rename_terminal_states({"0": expected[0], "1": expected[1]})
            np.testing.assert_array_equal(g.terminal_states.cat.categories, expected)
            np.testing.assert_array_equal(g.terminal_states_memberships.names, expected)

        np.testing.assert_array_equal(
            g.adata.obs[Key.obs.term_states(g.backward, bwd=initial)].cat.categories,
            expected,
        )


class TestGPCCAIO:
    @pytest.mark.parametrize("deep", [False, True])
    def test_copy(self, adata_gpcca_fwd: tuple[AnnData, cr.estimators.GPCCA], deep: bool):
        _, mc1 = adata_gpcca_fwd
        mc2 = mc1.copy(deep=deep)

        assert_estimators_equal(mc1, mc2, copy=True, deep=deep)

    def test_read(self, adata_gpcca_fwd: tuple[AnnData, cr.estimators.GPCCA], tmpdir):
        _, mc1 = adata_gpcca_fwd

        mc1.write(os.path.join(tmpdir, "foo.pkl"))
        mc2 = cr.estimators.GPCCA.read(os.path.join(tmpdir, "foo.pkl"))

        assert_estimators_equal(mc1, mc2)

    @pytest.mark.parametrize("copy", [False, True])
    def test_write_no_adata(
        self,
        adata_gpcca_fwd: tuple[AnnData, cr.estimators.GPCCA],
        copy: bool,
        tmpdir,
    ):
        adata, mc1 = adata_gpcca_fwd

        mc1.write(os.path.join(tmpdir, "foo.bar"), write_adata=False)
        mc2 = cr.estimators.GPCCA.read(os.path.join(tmpdir, "foo.bar"), adata=adata, copy=copy)

        if copy:
            assert adata is not mc2.adata
        else:
            assert adata is mc2.adata
        assert_estimators_equal(mc1, mc2)

    def test_write_no_adata_read_none_supplied(self, adata_gpcca_fwd: tuple[AnnData, cr.estimators.GPCCA], tmpdir):
        _, mc1 = adata_gpcca_fwd

        mc1.write(os.path.join(tmpdir, "foo.pkl"), write_adata=False)
        with pytest.raises(TypeError, match="This object was saved without"):
            _ = cr.estimators.GPCCA.read(os.path.join(tmpdir, "foo.pkl"), adata=None)

    def test_write_no_adata_read_wrong_length(self, adata_gpcca_fwd: tuple[AnnData, cr.estimators.GPCCA], tmpdir):
        rng = np.random.default_rng()
        _, mc1 = adata_gpcca_fwd
        adata = AnnData(rng.normal(size=(len(mc1) + 1, 1)))

        mc1.write(os.path.join(tmpdir, "foo.pkl"), write_adata=False)
        with pytest.raises(ValueError, match="Expected `adata` to be of length"):
            _ = cr.estimators.GPCCA.read(os.path.join(tmpdir, "foo.pkl"), adata=adata)

    @pytest.mark.parametrize("verbose", [None, False])
    def test_compute_schur_verbosity(self, adata_large: AnnData, verbose: Optional[bool], capsys):
        _ = pytest.importorskip("petsc4py")
        _ = pytest.importorskip("slepc4py")

        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4.0)
        g = cr.estimators.GPCCA(vk)

        g.compute_schur(n_components=10, method="krylov", verbose=verbose)
        out, _ = capsys.readouterr()

        assert not len(out)
