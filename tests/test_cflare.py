# -*- coding: utf-8 -*-
from typing import Tuple

import pytest

from anndata import AnnData

import numpy as np
import pandas as pd
from pandas.api.types import is_categorical_dtype

import cellrank as cr
import cellrank.tl.kernels._precomputed_kernel
from cellrank.tl.kernels import VelocityKernel, ConnectivityKernel
from cellrank.tl._constants import (
    Direction,
    DirPrefix,
    AbsProbKey,
    FinalStatesKey,
    _probs,
    _colors,
    _lin_names,
)
from cellrank.tl.estimators._constants import A, P

EPS = np.finfo(np.float64).eps


class TestCFLARE:
    def test_compute_partition(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(final_kernel)
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
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(final_kernel)
        mc.compute_eigendecomposition(k=2)

        assert isinstance(mc.eigendecomposition, dict)
        assert set(mc._get(P.EIG).keys()) == {
            "D",
            "V_l",
            "V_r",
            "eigengap",
            "params",
            "stationary_dist",
        }
        assert f"eig_{Direction.FORWARD}" in mc.adata.uns.keys()

    def test_compute_final_states_no_eig(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(final_kernel)
        with pytest.raises(RuntimeError):
            mc.compute_final_states(use=2)

    def test_compute_final_states_too_large_use(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(final_kernel)
        mc.compute_eigendecomposition(k=2)
        with pytest.raises(ValueError):
            mc.compute_final_states(use=1000)

    def test_compute_approx_normal_run(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(final_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.compute_final_states(use=2)

        assert is_categorical_dtype(mc._get(P.FIN))
        assert mc._get(P.FIN_PROBS) is not None

        assert FinalStatesKey.FORWARD.s in mc.adata.obs.keys()
        assert _probs(FinalStatesKey.FORWARD) in mc.adata.obs.keys()
        assert _colors(FinalStatesKey.FORWARD) in mc.adata.uns.keys()

    def test_compute_absorption_probabilities_no_args(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(final_kernel)
        with pytest.raises(RuntimeError):
            mc.compute_absorption_probabilities()

    def test_compute_absorption_probabilities_normal_run(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(final_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.compute_final_states(use=2)
        mc.compute_absorption_probabilities()

        assert isinstance(mc._get(P.DIFF_POT), pd.Series)
        assert f"{AbsProbKey.FORWARD}_dp" in mc.adata.obs.keys()
        np.testing.assert_array_equal(
            mc._get(P.DIFF_POT), mc.adata.obs[f"{AbsProbKey.FORWARD}_dp"]
        )

        assert isinstance(mc._get(P.ABS_PROBS), cr.tl.Lineage)
        assert mc._get(P.ABS_PROBS).shape == (mc.adata.n_obs, 2)
        assert f"{AbsProbKey.FORWARD}" in mc.adata.obsm.keys()
        np.testing.assert_array_equal(
            mc._get(P.ABS_PROBS).X, mc.adata.obsm[f"{AbsProbKey.FORWARD}"]
        )

        assert _lin_names(AbsProbKey.FORWARD) in mc.adata.uns.keys()
        np.testing.assert_array_equal(
            mc._get(P.ABS_PROBS).names,
            mc.adata.uns[_lin_names(AbsProbKey.FORWARD)],
        )

        assert _colors(AbsProbKey.FORWARD) in mc.adata.uns.keys()
        np.testing.assert_array_equal(
            mc._get(P.ABS_PROBS).colors,
            mc.adata.uns[_colors(AbsProbKey.FORWARD)],
        )
        np.testing.assert_allclose(mc._get(P.ABS_PROBS).X.sum(1), 1)

    def test_compute_absorption_probabilities_solver(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck
        tol = 1e-6

        mc = cr.tl.estimators.CFLARE(final_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.compute_final_states(use=2)

        # compute lin probs using direct solver
        mc.compute_absorption_probabilities(solver="direct")
        l_direct = mc._get(P.ABS_PROBS).copy()

        # compute lin probs using iterative solver
        mc.compute_absorption_probabilities(solver="gmres", tol=tol)
        l_iterative = mc._get(P.ABS_PROBS).copy()

        assert not np.shares_memory(l_direct.X, l_iterative.X)  # sanity check
        np.testing.assert_allclose(l_direct.X, l_iterative.X, rtol=0, atol=tol)

    def test_compute_absorption_probabilities_solver_petsc(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck
        tol = 1e-6

        mc = cr.tl.estimators.CFLARE(final_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.compute_final_states(use=2)

        # compute lin probs using direct solver
        mc.compute_absorption_probabilities(solver="gmres", use_petsc=False, tol=tol)
        l_iter = mc._get(P.ABS_PROBS).copy()

        # compute lin probs using petsc iterative solver
        mc.compute_absorption_probabilities(solver="gmres", use_petsc=True, tol=tol)
        l_iter_petsc = mc._get(P.ABS_PROBS).copy()

        assert not np.shares_memory(l_iter.X, l_iter_petsc.X)  # sanity check
        np.testing.assert_allclose(l_iter.X, l_iter_petsc.X, rtol=0, atol=tol)

    def test_compute_absorption_probabilities_lineage_absorption_mean(
        self, adata_large: AnnData
    ):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck
        tol = 1e-6

        mc = cr.tl.estimators.CFLARE(final_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.compute_final_states(use=2)

        # compute lin probs using direct solver
        mc.compute_absorption_probabilities(
            solver="gmres",
            use_petsc=False,
            tol=tol,
            time_to_absorption="0",
        )
        at = mc._get(P.LIN_ABS_TIMES)

        assert isinstance(at, pd.DataFrame)
        np.testing.assert_array_equal(at.index, adata_large.obs_names)
        np.testing.assert_array_equal(at.columns, ["0 mean"])

    def test_compute_absorption_probabilities_lineage_absorption_var(
        self, adata_large: AnnData
    ):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck
        tol = 1e-6

        mc = cr.tl.estimators.CFLARE(final_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.compute_final_states(use=2)

        # compute lin probs using direct solver
        mc.compute_absorption_probabilities(
            solver="gmres", use_petsc=False, tol=tol, time_to_absorption={"0": "var"}
        )
        at = mc._get(P.LIN_ABS_TIMES)

        assert isinstance(at, pd.DataFrame)
        np.testing.assert_array_equal(at.index, adata_large.obs_names)
        np.testing.assert_array_equal(at.columns, ["0 mean", "0 var"])

    def test_compute_absorption_probabilities_lineage_absorption_all(
        self, adata_large: AnnData
    ):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck
        tol = 1e-6

        mc = cr.tl.estimators.CFLARE(final_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.compute_final_states(use=2)

        # compute lin probs using direct solver
        mc.compute_absorption_probabilities(
            solver="gmres", use_petsc=False, tol=tol, time_to_absorption={"all": "var"}
        )
        name = ", ".join(mc._get(P.ABS_PROBS).names)
        at = mc._get(P.LIN_ABS_TIMES)

        assert isinstance(at, pd.DataFrame)
        np.testing.assert_array_equal(at.index, adata_large.obs_names)
        np.testing.assert_array_equal(at.columns, [f"{name} mean", f"{name} var"])

    def test_compute_lineage_drivers_no_lineages(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(final_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.compute_final_states(use=2)
        with pytest.raises(RuntimeError):
            mc.compute_lineage_drivers()

    def test_compute_lineage_drivers_invalid_lineages(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(final_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.compute_final_states(use=2)
        mc.compute_absorption_probabilities()
        with pytest.raises(KeyError):
            mc.compute_lineage_drivers(use_raw=False, lineages=["foo"])

    def test_compute_lineage_drivers_invalid_clusters(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(final_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.compute_final_states(use=2)
        mc.compute_absorption_probabilities()
        with pytest.raises(KeyError):
            mc.compute_lineage_drivers(
                use_raw=False, cluster_key="clusters", clusters=["foo"]
            )

    def test_compute_lineage_drivers_normal_run(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(final_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.compute_final_states(use=2)
        mc.compute_absorption_probabilities()
        mc.compute_lineage_drivers(use_raw=False, cluster_key="clusters")

        for lineage in ["0", "1"]:
            assert f"{DirPrefix.FORWARD} {lineage}" in mc.adata.var.keys()

    def test_plot_lineage_drivers_not_computed(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(final_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.compute_final_states(use=2)
        mc.compute_absorption_probabilities()

        with pytest.raises(RuntimeError):
            mc.plot_lineage_drivers("0")

    def test_plot_lineage_drivers_invalid_name(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(final_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.compute_final_states(use=2)
        mc.compute_absorption_probabilities()
        mc.compute_lineage_drivers(use_raw=False, cluster_key="clusters")

        with pytest.raises(KeyError):
            mc.plot_lineage_drivers("foo", use_raw=False)

    def test_plot_lineage_drivers_invalid_n_genes(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(final_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.compute_final_states(use=2)
        mc.compute_absorption_probabilities()
        mc.compute_lineage_drivers(use_raw=False, cluster_key="clusters")

        with pytest.raises(ValueError):
            mc.plot_lineage_drivers("0", use_raw=False, n_genes=0)

    def test_plot_lineage_drivers_normal_run(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.estimators.CFLARE(final_kernel)
        mc.compute_eigendecomposition(k=5)
        mc.compute_final_states(use=2)
        mc.compute_absorption_probabilities()
        mc.compute_lineage_drivers(use_raw=False, cluster_key="clusters")
        mc.plot_lineage_drivers("0", use_raw=False)

    def test_compute_absorption_probabilities_keys_colors(self, adata_large: AnnData):
        adata = adata_large
        vk = VelocityKernel(adata).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc_fwd = cr.tl.estimators.CFLARE(final_kernel)
        mc_fwd.compute_partition()
        mc_fwd.compute_eigendecomposition()
        mc_fwd.compute_final_states(use=3)

        arcs = ["0", "2"]
        arc_colors = [
            c
            for arc, c in zip(
                mc_fwd._get(P.FIN).cat.categories, mc_fwd._get(A.FIN_COLORS)
            )
            if arc in arcs
        ]

        mc_fwd.compute_absorption_probabilities(keys=arcs)
        lin_colors = mc_fwd._get(P.ABS_PROBS)[arcs].colors

        np.testing.assert_array_equal(arc_colors, lin_colors)

    def compare_absorption_probabilites_with_reference(self):
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

        # initialise a pre-computed kernel and CFLARE estimator object
        c = cr.tl.estimators.CFLARE(
            cellrank.tl.kernels._precomputed_kernel.PrecomputedKernel(transition_matrix)
        )

        # define the set of metastable states
        state_annotation = pd.Series(index=range(p.shape[0]))
        state_annotation[7] = "terminal_1"
        state_annotation[10] = "terminal_2"
        state_annotation = state_annotation.astype("category")
        c._set(A.FIN, state_annotation)

        # compute absorption probabilities
        c.compute_absorption_probabilities()
        absorption_probabilities_query = c.get(P.ABS_PROBS)[state_annotation.isna()]

        # check whether these two agree
        np.allclose(absorption_probabilities_query, absorption_probabilities_reference)

    def test_manual_approx_rc_set(self, adata_large):
        adata = adata_large
        vk = VelocityKernel(adata).compute_transition_matrix(softmax_scale=4)
        ck = ConnectivityKernel(adata).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc_fwd = cr.tl.estimators.CFLARE(final_kernel)
        mc_fwd.compute_partition()
        mc_fwd.compute_eigendecomposition()

        mc_fwd.compute_final_states(use=3)
        original = np.array(adata.obs[f"{FinalStatesKey.FORWARD}"].copy())
        zero_mask = original == "0"

        cells = list(adata[zero_mask].obs_names)
        mc_fwd.set_final_states({"foo": cells})

        assert (adata.obs[f"{FinalStatesKey.FORWARD}"][zero_mask] == "foo").all()
        assert pd.isna(adata.obs[f"{FinalStatesKey.FORWARD}"][~zero_mask]).all()


class TestCFLARECopy:
    def test_copy_simple(
        self, adata_cflare_fwd: Tuple[AnnData, cr.tl.estimators.CFLARE]
    ):
        _, mc1 = adata_cflare_fwd
        mc2 = mc1.copy()

        assert mc1 is not mc2
        assert mc1.adata is not mc2.adata
        assert mc1.kernel is not mc2.kernel

    # TODO - copy tests, saving tests, fit tests
    # careful in 3.6 - can't pickle enums
