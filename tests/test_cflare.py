from typing import Tuple

import os
import pytest
from _helpers import assert_estimators_equal

import cellrank as cr
from anndata import AnnData
from cellrank._key import Key
from cellrank.tl.kernels import VelocityKernel, ConnectivityKernel

import numpy as np
import pandas as pd
from pandas.api.types import is_categorical_dtype

EPS = np.finfo(np.float64).eps


class TestCFLARE:
    def test_compute_eigendecomposition(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(terminal_kernel)
        mc.compute_eigendecomposition(k=2)

        assert isinstance(mc.eigendecomposition, dict)
        assert set(mc.eigendecomposition.keys()) == {
            "D",
            "V_l",
            "V_r",
            "eigengap",
            "params",
            "stationary_dist",
        }
        assert Key.uns.eigen(mc.backward) in mc.adata.uns

    def test_compute_terminal_states_no_eig(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(terminal_kernel)
        with pytest.raises(RuntimeError):
            mc.compute_terminal_states(use=2)

    def test_compute_terminal_states_too_large_use(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(terminal_kernel)
        mc.compute_eigendecomposition(k=2)
        with pytest.raises(ValueError):
            mc.compute_terminal_states(use=1000)

    def test_compute_approx_normal_run(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(terminal_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.compute_terminal_states(use=2)

        assert is_categorical_dtype(mc.terminal_states)
        assert mc.terminal_states_probabilities is not None

        key = Key.obs.term_states(mc.backward)
        assert key in mc.adata.obs
        assert Key.obs.probs(key) in mc.adata.obs
        assert Key.uns.colors(key) in mc.adata.uns

    def test_rename_terminal_states_no_terminal_states(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(terminal_kernel)
        mc.compute_eigendecomposition(k=5)
        with pytest.raises(RuntimeError):
            mc.rename_terminal_states({"foo": "bar"})

    def test_rename_terminal_states_invalid_old_name(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(terminal_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.compute_terminal_states(use=2)
        with pytest.raises(ValueError):
            mc.rename_terminal_states({"foo": "bar"})

    def test_rename_terminal_states_invalid_new_name(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(terminal_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.compute_terminal_states(use=2, method="kmeans")
        with pytest.raises(ValueError):
            mc.rename_terminal_states({"0": "1"})

    def test_rename_terminal_states_try_joining_states(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(terminal_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.compute_terminal_states(use=2)
        with pytest.raises(ValueError):
            mc.rename_terminal_states({"0, 1": "foo"})

    def test_rename_terminal_states_empty_mapping(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(terminal_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.compute_terminal_states(use=2)

        orig_cats = list(mc.terminal_states.cat.categories)
        mc.rename_terminal_states({})

        np.testing.assert_array_equal(mc.terminal_states.cat.categories, orig_cats)

    def test_rename_terminal_states_normal_run(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(terminal_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.compute_terminal_states(use=2, method="kmeans")
        mc.rename_terminal_states({"0": "foo", "1": "bar"})

        np.testing.assert_array_equal(mc.terminal_states.cat.categories, ["foo", "bar"])
        np.testing.assert_array_equal(
            mc.adata.obs[Key.obs.term_states(mc.backward)].cat.categories,
            ["foo", "bar"],
        )

    def test_compute_absorption_probabilities_no_args(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(terminal_kernel)
        with pytest.raises(RuntimeError):
            mc.compute_absorption_probabilities()

    def test_compute_absorption_probabilities_normal_run(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(terminal_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.compute_terminal_states(use=2, method="kmeans")
        mc.compute_absorption_probabilities()
        mc.compute_lineage_priming()

        key = Key.obs.priming_degree(mc.backward)
        assert isinstance(mc.priming_degree, pd.Series)
        assert key in mc.adata.obs
        np.testing.assert_array_equal(mc.priming_degree, mc.adata.obs[key])

        key = Key.obsm.abs_probs(mc.backward)
        assert isinstance(mc.absorption_probabilities, cr.tl.Lineage)
        assert mc.absorption_probabilities.shape == (mc.adata.n_obs, 2)
        assert key in mc.adata.obsm
        assert isinstance(mc.adata.obsm[key], cr.tl.Lineage)
        np.testing.assert_array_equal(mc.absorption_probabilities.X, mc.adata.obsm[key])

        np.testing.assert_array_equal(
            mc.absorption_probabilities.names,
            mc.adata.obs[Key.obs.term_states(mc.backward)].cat.categories,
        )

        key = Key.uns.colors(Key.obs.term_states(mc.backward))
        assert key in mc.adata.uns
        np.testing.assert_array_equal(
            mc.absorption_probabilities.colors, mc.adata.uns[key]
        )
        np.testing.assert_array_equal(mc._term_states_colors, mc.adata.uns[key])
        np.testing.assert_allclose(mc.absorption_probabilities.X.sum(1), 1, rtol=1e-6)

    def test_compute_absorption_probabilities_solver(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck
        tol = 1e-6

        mc = cr.tl.estimators.CFLARE(terminal_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.compute_terminal_states(use=2)

        # compute lin probs using direct solver
        mc.compute_absorption_probabilities(solver="direct")
        l_direct = mc.absorption_probabilities.copy()

        # compute lin probs using iterative solver
        mc.compute_absorption_probabilities(solver="gmres", tol=tol)
        l_iterative = mc.absorption_probabilities.copy()

        assert not np.shares_memory(l_direct.X, l_iterative.X)  # sanity check
        np.testing.assert_allclose(l_direct.X, l_iterative.X, rtol=0, atol=tol)

    def test_compute_absorption_probabilities_solver_petsc(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck
        tol = 1e-6

        mc = cr.tl.estimators.CFLARE(terminal_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.compute_terminal_states(use=2, method="kmeans")

        # compute lin probs using direct solver
        mc.compute_absorption_probabilities(solver="gmres", use_petsc=False, tol=tol)
        l_iter = mc.absorption_probabilities.copy()

        # compute lin probs using petsc iterative solver
        mc.compute_absorption_probabilities(solver="gmres", use_petsc=True, tol=tol)
        l_iter_petsc = mc.absorption_probabilities.copy()

        assert not np.shares_memory(l_iter.X, l_iter_petsc.X)  # sanity check
        np.testing.assert_allclose(l_iter.X, l_iter_petsc.X, rtol=0, atol=tol)

    def test_compute_absorption_probabilities_lineage_absorption_mean(
        self, adata_large: AnnData
    ):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck
        tol = 1e-6

        mc = cr.tl.estimators.CFLARE(terminal_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.compute_terminal_states(use=2)

        # compute lin probs using direct solver
        mc.compute_absorption_probabilities(
            solver="gmres",
            use_petsc=False,
            tol=tol,
            time_to_absorption="0",
        )
        at = mc.absorption_times

        assert isinstance(at, pd.DataFrame)
        np.testing.assert_array_equal(at.index, adata_large.obs_names)
        np.testing.assert_array_equal(at.columns, ["0 mean"])

    def test_compute_absorption_probabilities_lineage_absorption_var(
        self, adata_large: AnnData
    ):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck
        tol = 1e-6

        mc = cr.tl.estimators.CFLARE(terminal_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.compute_terminal_states(use=2)

        # compute lin probs using direct solver
        mc.compute_absorption_probabilities(
            solver="gmres", use_petsc=False, tol=tol, time_to_absorption={"0": "var"}
        )
        at = mc.absorption_times

        assert isinstance(at, pd.DataFrame)
        np.testing.assert_array_equal(at.index, adata_large.obs_names)
        np.testing.assert_array_equal(at.columns, ["0 mean", "0 var"])

    def test_compute_absorption_probabilities_lineage_absorption_all(
        self, adata_large: AnnData
    ):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck
        tol = 1e-6

        mc = cr.tl.estimators.CFLARE(terminal_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.compute_terminal_states(use=2)

        # compute lin probs using direct solver
        mc.compute_absorption_probabilities(
            solver="gmres", use_petsc=False, tol=tol, time_to_absorption={"all": "var"}
        )
        name = ", ".join(mc.absorption_probabilities.names)
        at = mc.absorption_times

        assert isinstance(at, pd.DataFrame)
        np.testing.assert_array_equal(at.index, adata_large.obs_names)
        np.testing.assert_array_equal(at.columns, [f"{name} mean", f"{name} var"])

    def test_compute_lineage_drivers_no_lineages(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(terminal_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.compute_terminal_states(use=2)
        with pytest.raises(RuntimeError):
            mc.compute_lineage_drivers()

    def test_compute_lineage_drivers_invalid_lineages(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(terminal_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.compute_terminal_states(use=2)
        mc.compute_absorption_probabilities()
        with pytest.raises(KeyError):
            mc.compute_lineage_drivers(use_raw=False, lineages=["foo"])

    def test_compute_lineage_drivers_invalid_clusters(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(terminal_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.compute_terminal_states(use=2)
        mc.compute_absorption_probabilities()
        with pytest.raises(KeyError):
            mc.compute_lineage_drivers(
                use_raw=False, cluster_key="clusters", clusters=["foo"]
            )

    def test_compute_lineage_drivers_normal_run(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(terminal_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.compute_terminal_states(use=2, method="kmeans")
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

        mc = cr.tl.estimators.CFLARE(terminal_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.compute_terminal_states(use=2)
        mc.compute_absorption_probabilities()

        with pytest.raises(RuntimeError):
            mc.plot_lineage_drivers("0")

    def test_plot_lineage_drivers_invalid_name(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(terminal_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.compute_terminal_states(use=2)
        mc.compute_absorption_probabilities()
        mc.compute_lineage_drivers(use_raw=False, cluster_key="clusters")

        with pytest.raises(KeyError):
            mc.plot_lineage_drivers("foo", use_raw=False)

    def test_plot_lineage_drivers_invalid_n_genes(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(terminal_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.compute_terminal_states(use=2)
        mc.compute_absorption_probabilities()
        mc.compute_lineage_drivers(use_raw=False, cluster_key="clusters")

        with pytest.raises(ValueError):
            mc.plot_lineage_drivers("0", use_raw=False, n_genes=0)

    def test_plot_lineage_drivers_normal_run(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(terminal_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.compute_terminal_states(use=2)
        mc.compute_absorption_probabilities()
        mc.compute_lineage_drivers(use_raw=False, cluster_key="clusters")
        mc.plot_lineage_drivers("0", use_raw=False)

    def test_compute_absorption_probabilities_keys_colors(self, adata_large: AnnData):
        adata = adata_large
        vk = VelocityKernel(adata).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc_fwd = cr.tl.estimators.CFLARE(terminal_kernel)
        mc_fwd.compute_eigendecomposition()
        mc_fwd.compute_terminal_states(use=3, method="kmeans")

        arcs = ["0", "2"]
        arc_colors = [
            c
            for arc, c in zip(
                mc_fwd.terminal_states.cat.categories, mc_fwd._term_states_colors
            )
            if arc in arcs
        ]

        mc_fwd.compute_absorption_probabilities(keys=arcs)
        lin_colors = mc_fwd.absorption_probabilities[arcs].colors

        np.testing.assert_array_equal(arc_colors, lin_colors)

    def test_compare_absorption_probabilities_with_reference(self):
        # define a reference transition matrix. This is an absorbing MC with 2 absorbing states
        transition_matrix = np.array(
            [
                # 0.   1.   2.   3.   4.   5.   6.   7.   8.   9.   10.  11.
                [0.0, 0.8, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # 0.
                [0.2, 0.0, 0.6, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # 1.
                [0.6, 0.2, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # 2.
                [0.0, 0.05, 0.05, 0.0, 0.45, 0.45, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # 3.
                [0.0, 0.0, 0.0, 0.25, 0.0, 0.25, 0.4, 0.0, 0.0, 0.1, 0.0, 0.0],  # 4.
                [0.0, 0.0, 0.0, 0.25, 0.25, 0.0, 0.1, 0.0, 0.0, 0.4, 0.0, 0.0],  # 5.
                [0.0, 0.0, 0.0, 0.0, 0.05, 0.05, 0.0, 0.7, 0.2, 0.0, 0.0, 0.0],  # 6.
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],  # 7.
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8, 0.2, 0.0, 0.0, 0.0, 0.0],  # 8.
                [0.0, 0.0, 0.0, 0.0, 0.05, 0.05, 0.0, 0.0, 0.0, 0.0, 0.7, 0.2],  # 9.
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],  # 10.
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8, 0.2, 0.0],  # 11.
            ]
        )

        absorption_probabilities_reference = np.array(
            [
                [0.5, 0.5],
                [0.5, 0.5],
                [0.5, 0.5],
                [0.5, 0.5],
                [0.60571429, 0.39428571],
                [0.39428571, 0.60571429],
                [0.94047619, 0.05952381],
                [0.95238095, 0.04761905],
                [0.05952381, 0.94047619],
                [0.04761905, 0.95238095],
            ]
        )

        c = cr.tl.estimators.CFLARE(cr.tl.kernels.PrecomputedKernel(transition_matrix))

        state_annotation = pd.Series(index=range(len(c)))
        state_annotation[7] = "terminal_1"
        state_annotation[10] = "terminal_2"
        state_annotation = state_annotation.astype("category")
        c._term_states = state_annotation
        c._term_states_colors = np.array(["#000000", "#ffffff"])

        # compute absorption probabilities
        c.compute_absorption_probabilities()
        absorption_probabilities_query = c.absorption_probabilities[
            state_annotation.isna()
        ]

        np.allclose(absorption_probabilities_query, absorption_probabilities_reference)

    def test_manual_approx_rc_set(self, adata_large):
        adata = adata_large
        vk = VelocityKernel(adata).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc_fwd = cr.tl.estimators.CFLARE(terminal_kernel)
        mc_fwd.compute_eigendecomposition()
        key = Key.obs.term_states(mc_fwd.backward)

        mc_fwd.compute_terminal_states(use=3)
        original = np.array(adata.obs[key].copy())
        zero_mask = original == "0"

        cells = list(adata[zero_mask].obs_names)
        mc_fwd.set_terminal_states({"foo": cells})

        assert (adata.obs[key][zero_mask] == "foo").all()
        assert pd.isna(adata.obs[key][~zero_mask]).all()

    def test_abs_probs_do_not_sum_to_1(self, adata_large: AnnData, mocker):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(terminal_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.set_terminal_states(
            {"x": adata_large.obs_names[:3], "y": adata_large.obs_names[3:6]}
        )

        n_term = np.sum(~pd.isnull(mc.terminal_states))
        abs_prob = np.zeros((adata_large.n_obs - n_term, n_term))
        abs_prob[:, 0] = 1.0
        abs_prob[0, 0] = 1.01
        mocker.patch(
            "cellrank.tl.estimators.mixins._absorption_probabilities._solve_lin_system",
            return_value=abs_prob,
        )

        with pytest.raises(
            ValueError, match=r"`1` value\(s\) do not sum to 1 \(rtol=1e-3\)."
        ):
            mc.compute_absorption_probabilities()

    def test_abs_probs_negative(self, adata_large: AnnData, mocker):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        terminal_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(terminal_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.set_terminal_states(
            {"x": adata_large.obs_names[:3], "y": adata_large.obs_names[3:6]}
        )

        n_term = np.sum(~pd.isnull(mc.terminal_states))
        abs_prob = np.zeros((adata_large.n_obs - n_term, n_term))
        abs_prob[:, 0] = 1.0
        abs_prob[0, 0] = -0.5
        abs_prob[0, 1] = -1.5
        mocker.patch(
            "cellrank.tl.estimators.mixins._absorption_probabilities._solve_lin_system",
            return_value=abs_prob,
        )

        with pytest.raises(ValueError, match=r"`2` value\(s\) are negative."):
            mc.compute_absorption_probabilities()


class TestCFLAREIO:
    @pytest.mark.parametrize("deep", [False, True])
    def test_copy(
        self, adata_cflare_fwd: Tuple[AnnData, cr.tl.estimators.CFLARE], deep: bool
    ):
        _, mc1 = adata_cflare_fwd
        mc2 = mc1.copy(deep=deep)

        assert_estimators_equal(mc1, mc2, copy=True, deep=deep)

    def test_read(
        self, adata_cflare_fwd: Tuple[AnnData, cr.tl.estimators.CFLARE], tmpdir
    ):
        _, mc1 = adata_cflare_fwd

        mc1.write(os.path.join(tmpdir, "foo.pickle"))
        mc2 = cr.tl.estimators.CFLARE.read(os.path.join(tmpdir, "foo.pickle"))

        assert_estimators_equal(mc1, mc2)
