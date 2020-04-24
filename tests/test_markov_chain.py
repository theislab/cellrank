# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import pytest
import cellrank as cr

from anndata import AnnData
from pandas.api.types import is_categorical_dtype
from cellrank.tools._constants import (
    Prefix,
    RcKey,
    Direction,
    _colors,
    _probs,
    LinKey,
    _lin_names,
)
from cellrank.tools.kernels import VelocityKernel, ConnectivityKernel
from _helpers import assert_array_nan_equal


class TestMarkovChain:
    def foo(self, adata):
        vk = VelocityKernel(adata).compute_transition_matrix()
        ck = ConnectivityKernel(adata).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc_fwd = cr.tl.MarkovChain(final_kernel)

        mc_fwd.compute_lin_probs()

        mc_fwd.compute_lineage_drivers(cluster_key="clusters", use_raw=False)

    def test_compute_partition(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.MarkovChain(final_kernel)
        mc.compute_partition()

        assert isinstance(mc.irreducible, bool)
        if not mc.irreducible:
            assert isinstance(mc.recurrent_classes, list)
            assert isinstance(mc.transient_classes, list)
            assert f"{RcKey.FORWARD}_rec_classes" in mc.adata.obs
            assert f"{RcKey.FORWARD}_trans_classes" in mc.adata.obs
        else:
            assert mc.recurrent_classes is None
            assert mc.transient_classes is None

    def test_compute_eig(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.MarkovChain(final_kernel)
        mc.compute_eig(k=2)

        assert isinstance(mc.eigendecomposition, dict)
        assert set(mc.eigendecomposition.keys()) == {
            "D",
            "V_l",
            "V_r",
            "eigengap",
            "params",
        }
        assert f"eig_{Direction.FORWARD}" in mc.adata.uns.keys()

    def test_compute_approx_rcs_no_eig(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.MarkovChain(final_kernel)
        with pytest.raises(RuntimeError):
            mc.compute_approx_rcs(use=2)

    def test_compute_approx_rcs_too_large_use(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.MarkovChain(final_kernel)
        mc.compute_eig(k=2)
        with pytest.raises(ValueError):
            mc.compute_approx_rcs(use=1000)

    def test_compute_approx_normal_run(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.MarkovChain(final_kernel)
        mc.compute_eig(k=5)
        mc.compute_approx_rcs(use=2)

        assert is_categorical_dtype(mc.approx_recurrent_classes)
        assert mc.approx_recurrent_classes_probabilities is not None
        assert _colors(RcKey.FORWARD) in mc.adata.uns.keys()
        assert _probs(RcKey.FORWARD) in mc.adata.obs.keys()

    def test_compute_lin_probs_no_arcs(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.MarkovChain(final_kernel)
        with pytest.raises(RuntimeError):
            mc.compute_lin_probs()

    def test_compute_lin_probs_normal_run(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.MarkovChain(final_kernel)
        mc.compute_eig(k=5)
        mc.compute_approx_rcs(use=2)
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

    def test_compute_lineage_drivers_no_lineages(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.MarkovChain(final_kernel)
        mc.compute_eig(k=5)
        mc.compute_approx_rcs(use=2)
        with pytest.raises(RuntimeError):
            mc.compute_lineage_drivers()

    def test_compute_lineage_drivers_invalid_lineages(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.MarkovChain(final_kernel)
        mc.compute_eig(k=5)
        mc.compute_approx_rcs(use=2)
        mc.compute_lin_probs()
        with pytest.raises(KeyError):
            mc.compute_lineage_drivers(use_raw=False, lin_names=["foo"])

    def test_compute_lineage_drivers_invalid_clusters(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.MarkovChain(final_kernel)
        mc.compute_eig(k=5)
        mc.compute_approx_rcs(use=2)
        mc.compute_lin_probs()
        with pytest.raises(KeyError):
            mc.compute_lineage_drivers(
                use_raw=False, cluster_key="clusters", clusters=["foo"]
            )

    def test_compute_lineage_drivers_normal_run(self, adata_large: AnnData):
        vk = VelocityKernel(adata_large).compute_transition_matrix()
        ck = ConnectivityKernel(adata_large).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.MarkovChain(final_kernel)
        mc.compute_eig(k=5)
        mc.compute_approx_rcs(use=2)
        mc.compute_lin_probs()
        mc.compute_lineage_drivers(use_raw=False, cluster_key="clusters")

        for lineage in ["0", "1"]:
            assert f"{Prefix.FORWARD} {lineage} corr" in mc.adata.var.keys()

    def test_compute_lin_probs_keys_colors(self, adata_large: AnnData):
        adata = adata_large
        vk = VelocityKernel(adata).compute_transition_matrix()
        ck = ConnectivityKernel(adata).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc_fwd = cr.tl.MarkovChain(final_kernel)
        mc_fwd.compute_partition()
        mc_fwd.compute_eig()
        mc_fwd.compute_approx_rcs(use=3)

        arcs = ["0", "2"]
        arc_colors = [
            c
            for arc, c in zip(
                mc_fwd._approx_rcs.cat.categories, mc_fwd._approx_rcs_colors
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

        mc_fwd = cr.tl.MarkovChain(final_kernel)
        mc_fwd.compute_partition()
        mc_fwd.compute_eig()

        mc_fwd.compute_approx_rcs(use=3)
        original = np.array(adata.obs["final_cells"].copy())
        zero_mask = original == "0"

        cells = list(adata[zero_mask].obs_names)
        mc_fwd.set_approx_rcs({"foo": cells})

        assert (adata.obs["final_cells"][zero_mask] == "foo").all()
        assert pd.isna(adata.obs["final_cells"][~zero_mask]).all()

    def test_check_and_create_colors(self, adata_large):
        adata = adata_large
        vk = VelocityKernel(adata).compute_transition_matrix()
        ck = ConnectivityKernel(adata).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc_fwd = cr.tl.MarkovChain(final_kernel)
        mc_fwd.compute_partition()
        mc_fwd.compute_eig()

        mc_fwd.compute_approx_rcs(use=3)

        del mc_fwd._approx_rcs_colors
        del mc_fwd.uns["to_final_cells_colors"]

        mc_fwd._check_and_create_colors()

        assert "to_final_cells_colors" in mc_fwd.adata.uns
        np.testing.assert_array_equal(
            mc_fwd.adata.uns["fo_final_cells_colors"], _create_categorical_colors(3)
        )
        np.testing.assert_array_equal(
            mc_fwd.adata.uns["fo_final_cells_colors"], mc_fwd._approx_rcs_colors
        )


class TestMarkovChainCopy:
    def test_copy_simple(self, adata_mc_fwd):
        _, mc1 = adata_mc_fwd
        mc2 = mc1.copy()

        assert mc1 is not mc2
        assert mc1.adata is not mc2.adata
        assert mc1.kernel is not mc2.kernel

    def test_copy_deep(self, adata_mc_fwd):
        _, mc1 = adata_mc_fwd
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
        assert_array_nan_equal(
            mc1.approx_recurrent_classes, mc2.approx_recurrent_classes
        )
        np.testing.assert_array_equal(
            mc1.approx_recurrent_classes_probabilities,
            mc2.approx_recurrent_classes_probabilities,
        )
        np.testing.assert_array_equal(mc1._approx_rcs_colors, mc2._approx_rcs_colors)
        assert mc1._G2M_score == mc2._G2M_score
        assert mc1._S_score == mc2._S_score

    def test_copy_works(self, adata_mc_fwd):
        _, mc1 = adata_mc_fwd
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
        assert mc1.approx_recurrent_classes is not mc2.approx_recurrent_classes
        assert (
            mc1.approx_recurrent_classes_probabilities
            is not mc2.approx_recurrent_classes_probabilities
        )
        assert mc1._approx_rcs_colors is not mc2._approx_rcs_colors
        assert mc1._G2M_score != mc2._G2M_score
        assert mc1._S_score != mc2._S_score
