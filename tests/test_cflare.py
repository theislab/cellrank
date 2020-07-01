# -*- coding: utf-8 -*-
from typing import Tuple

import numpy as np
import pandas as pd
import pytest
from pandas.api.types import is_categorical_dtype

from anndata import AnnData

import cellrank as cr
from _helpers import assert_array_nan_equal
from cellrank.tools._colors import _create_categorical_colors
from cellrank.tools.kernels import VelocityKernel, ConnectivityKernel
from cellrank.tools._constants import (
    LinKey,
    Prefix,
    StateKey,
    Direction,
    _probs,
    _colors,
    _lin_names,
)

EPS = np.finfo(np.float64).eps


class TestCFLARE:
    def test_compute_partition(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.CFLARE(final_kernel)
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

        mc = cr.tl.CFLARE(final_kernel)
        mc.compute_eig(k=2)

        assert isinstance(mc.eigendecomposition, dict)
        assert set(mc.eigendecomposition.keys()) == {
            "D",
            "V_l",
            "V_r",
            "eigengap",
            "params",
            "stationary_dist",
        }
        assert f"eig_{Direction.FORWARD}" in mc.adata.uns.keys()

    def test_compute_metastable_states_no_eig(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.CFLARE(final_kernel)
        with pytest.raises(RuntimeError):
            mc.compute_metastable_states(use=2)

    def test_compute_metastable_states_too_large_use(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.CFLARE(final_kernel)
        mc.compute_eig(k=2)
        with pytest.raises(ValueError):
            mc.compute_metastable_states(use=1000)

    def test_compute_approx_normal_run(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.CFLARE(final_kernel)
        mc.compute_eig(k=5)
        mc.compute_metastable_states(use=2)

        assert is_categorical_dtype(mc.metastable_states)
        assert mc.metastable_states_probabilities is not None
        assert _colors(StateKey.FORWARD) in mc.adata.uns.keys()
        assert _probs(StateKey.FORWARD) in mc.adata.obs.keys()

    def test_compute_lin_probs_no_arcs(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.CFLARE(final_kernel)
        with pytest.raises(RuntimeError):
            mc.compute_lin_probs()

    def test_compute_lin_probs_normal_run(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.CFLARE(final_kernel)
        mc.compute_eig(k=5)
        mc.compute_metastable_states(use=2)
        mc.compute_lin_probs()

        assert isinstance(mc.diff_potential, np.ndarray)
        assert f"{LinKey.FORWARD}_dp" in mc.adata.obs.keys()
        np.testing.assert_array_equal(
            mc.diff_potential, mc.adata.obs[f"{LinKey.FORWARD}_dp"]
        )

        assert isinstance(mc.lineage_probabilities, cr.tl.Lineage)
        assert mc.lineage_probabilities.shape == (mc.adata.n_obs, 2)
        assert f"{LinKey.FORWARD}" in mc.adata.obsm.keys()
        np.testing.assert_array_equal(
            mc.lineage_probabilities.X, mc.adata.obsm[f"{LinKey.FORWARD}"]
        )

        assert _lin_names(LinKey.FORWARD) in mc.adata.uns.keys()
        np.testing.assert_array_equal(
            mc.lineage_probabilities.names, mc.adata.uns[_lin_names(LinKey.FORWARD)]
        )

        assert _colors(LinKey.FORWARD) in mc.adata.uns.keys()
        np.testing.assert_array_equal(
            mc.lineage_probabilities.colors, mc.adata.uns[_colors(LinKey.FORWARD)]
        )
        np.testing.assert_allclose(mc.lineage_probabilities.X.sum(1), 1)

    def test_compute_lin_probs_solver(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck
        tol = 1e-6

        mc = cr.tl.CFLARE(final_kernel)
        mc.compute_eig(k=5)
        mc.compute_metastable_states(use=2)

        # compute lin probs using direct solver
        mc.compute_lin_probs(use_iterative_solver=False)
        l_direct = mc.lineage_probabilities

        # comptue lin probs using iterative solver
        mc.compute_lin_probs(use_iterative_solver=True, tol=tol)
        l_iterative = mc.lineage_probabilities

        np.testing.assert_allclose(l_direct.X, l_iterative.X, rtol=0, atol=tol)

    def test_compute_lin_probs_initialization(self, adata_large: AnnData):
        # check whether changing the initial x0 vector for the iterative solver makes a difference
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck
        tol = 1e-6

        mc = cr.tl.CFLARE(final_kernel)
        mc.compute_eig(k=5)
        mc.compute_metastable_states(use=2)

        # compute lin probs using direct solver
        mc.compute_lin_probs(
            use_iterative_solver=True, use_initialization=True, tol=tol
        )
        l_0 = mc.lineage_probabilities

        # comptue lin probs using iterative solver
        mc.compute_lin_probs(
            use_iterative_solver=True, use_initialization=False, tol=tol
        )
        l_1 = mc.lineage_probabilities

        np.testing.assert_allclose(l_0.X, l_1.X, rtol=0, atol=tol)

    def test_compute_lineage_drivers_no_lineages(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.CFLARE(final_kernel)
        mc.compute_eig(k=5)
        mc.compute_metastable_states(use=2)
        with pytest.raises(RuntimeError):
            mc.compute_lineage_drivers()

    def test_compute_lineage_drivers_invalid_lineages(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.CFLARE(final_kernel)
        mc.compute_eig(k=5)
        mc.compute_metastable_states(use=2)
        mc.compute_lin_probs()
        with pytest.raises(KeyError):
            mc.compute_lineage_drivers(use_raw=False, lin_names=["foo"])

    def test_compute_lineage_drivers_invalid_clusters(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.CFLARE(final_kernel)
        mc.compute_eig(k=5)
        mc.compute_metastable_states(use=2)
        mc.compute_lin_probs()
        with pytest.raises(KeyError):
            mc.compute_lineage_drivers(
                use_raw=False, cluster_key="clusters", clusters=["foo"]
            )

    def test_compute_lineage_drivers_normal_run(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.CFLARE(final_kernel)
        mc.compute_eig(k=5)
        mc.compute_metastable_states(use=2)
        mc.compute_lin_probs()
        mc.compute_lineage_drivers(use_raw=False, cluster_key="clusters")

        for lineage in ["0", "1"]:
            assert f"{Prefix.FORWARD} {lineage} corr" in mc.adata.var.keys()

    def test_compute_lin_probs_keys_colors(self, adata_large: AnnData):
        adata = adata_large
        vk = VelocityKernel(adata).compute_transition_matrix()
        ck = ConnectivityKernel(adata).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc_fwd = cr.tl.CFLARE(final_kernel)
        mc_fwd.compute_partition()
        mc_fwd.compute_eig()
        mc_fwd.compute_metastable_states(use=3)

        arcs = ["0", "2"]
        arc_colors = [
            c
            for arc, c in zip(
                mc_fwd.metastable_states.cat.categories, mc_fwd._meta_states_colors
            )
            if arc in arcs
        ]

        mc_fwd.compute_lin_probs(keys=arcs)
        lin_colors = mc_fwd.lineage_probabilities[arcs].colors

        np.testing.assert_array_equal(arc_colors, lin_colors)

    def test_manual_approx_rc_set(self, adata_large):
        adata = adata_large
        vk = VelocityKernel(adata).compute_transition_matrix()
        ck = ConnectivityKernel(adata).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc_fwd = cr.tl.CFLARE(final_kernel)
        mc_fwd.compute_partition()
        mc_fwd.compute_eig()

        mc_fwd.compute_metastable_states(use=3)
        original = np.array(adata.obs[f"{StateKey.FORWARD}"].copy())
        zero_mask = original == "0"

        cells = list(adata[zero_mask].obs_names)
        mc_fwd.set_metastable_states({"foo": cells})

        assert (adata.obs[f"{StateKey.FORWARD}"][zero_mask] == "foo").all()
        assert pd.isna(adata.obs[f"{StateKey.FORWARD}"][~zero_mask]).all()

    def test_check_and_create_colors(self, adata_large):
        adata = adata_large
        vk = VelocityKernel(adata).compute_transition_matrix()
        ck = ConnectivityKernel(adata).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc_fwd = cr.tl.CFLARE(final_kernel)
        mc_fwd.compute_partition()
        mc_fwd.compute_eig()

        mc_fwd.compute_metastable_states(use=3)

        mc_fwd._meta_states_colors = None
        del mc_fwd.adata.uns[_colors(StateKey.FORWARD)]

        mc_fwd._check_and_create_colors()

        assert _colors(StateKey.FORWARD) in mc_fwd.adata.uns
        np.testing.assert_array_equal(
            mc_fwd.adata.uns[_colors(StateKey.FORWARD)], _create_categorical_colors(3)
        )
        np.testing.assert_array_equal(
            mc_fwd.adata.uns[_colors(StateKey.FORWARD)], mc_fwd._meta_states_colors
        )


class TestCFLARECopy:
    def test_copy_simple(self, adata_cflare_fwd: Tuple[AnnData, cr.tl.CFLARE]):
        _, mc1 = adata_cflare_fwd
        mc2 = mc1.copy()

        assert mc1 is not mc2
        assert mc1.adata is not mc2.adata
        assert mc1.kernel is not mc2.kernel

    def test_copy_deep(self, adata_cflare_fwd: Tuple[AnnData, cr.tl.CFLARE]):
        _, mc1 = adata_cflare_fwd
        mc2 = mc1.copy()

        assert mc1.irreducible == mc2.irreducible
        np.testing.assert_array_equal(mc1.recurrent_classes, mc2.recurrent_classes)
        np.testing.assert_array_equal(mc1.transient_classes, mc2.transient_classes)
        for k, v in mc1.eigendecomposition.items():
            if isinstance(v, np.ndarray):
                assert_array_nan_equal(v, mc2.eigendecomposition[k])
            else:
                assert v == mc2.eigendecomposition[k]
        np.testing.assert_array_equal(
            mc1.lineage_probabilities, mc2.lineage_probabilities
        )
        np.testing.assert_array_equal(mc1.diff_potential, mc2.diff_potential)
        assert_array_nan_equal(mc1.metastable_states, mc2.metastable_states)
        np.testing.assert_array_equal(
            mc1.metastable_states_probabilities, mc2.metastable_states_probabilities
        )
        np.testing.assert_array_equal(mc1._meta_states_colors, mc2._meta_states_colors)
        assert mc1._G2M_score == mc2._G2M_score
        assert mc1._S_score == mc2._S_score
        assert mc1._g2m_key == mc2._g2m_key
        assert mc1._s_key == mc2._s_key
        assert mc1._is_sparse == mc2._is_sparse

    def test_copy_works(self, adata_cflare_fwd: Tuple[AnnData, cr.tl.CFLARE]):
        _, mc1 = adata_cflare_fwd
        mc2 = mc1.copy()

        mc1._is_irreducible = not mc2.irreducible
        mc1.eigendecomposition["foo"] = "bar"
        mc1._rec_classes = "foo"
        mc1._trans_classes = "bar"
        mc1._G2M_score = "foo"
        mc1._S_score = "bar"

        assert mc1.irreducible != mc2.irreducible
        assert mc1.recurrent_classes is not mc2.recurrent_classes
        assert mc1.transient_classes is not mc2.transient_classes
        assert mc1.eigendecomposition != mc2.eigendecomposition
        assert mc1.lineage_probabilities is not mc2.lineage_probabilities
        assert mc1.diff_potential is not mc2.diff_potential
        assert mc1.metastable_states is not mc2.metastable_states
        assert (
            mc1.metastable_states_probabilities
            is not mc2.metastable_states_probabilities
        )
        assert mc1._meta_states_colors is not mc2._meta_states_colors
        assert mc1._G2M_score != mc2._G2M_score
        assert mc1._S_score != mc2._S_score
