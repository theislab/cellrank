# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import cellrank as cr

from cellrank.tools.kernels import VelocityKernel, ConnectivityKernel


class TestMarkovChain:
    def test_compute_lin_probs_keys_colors(self, adata_large):
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
