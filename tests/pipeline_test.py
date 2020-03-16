# -*- coding: utf-8 -*-
import unittest
import cellrank as cr
import matplotlib

matplotlib.use("Agg")

from cellrank.tools.kernels import VelocityKernel, ConnectivityKernel
from _helpers import create_dummy_adata


_adata = create_dummy_adata(50)


class PipeLineTestCase(unittest.TestCase):
    def test_default_fwd_pipelne(self):
        adata = _adata.copy()
        vk = VelocityKernel(adata).compute_transition_matrix()
        ck = ConnectivityKernel(adata).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc_fwd = cr.tl.MarkovChain(final_kernel)

        mc_fwd.compute_partition()

        mc_fwd.compute_eig()
        mc_fwd.plot_real_spectrum()
        mc_fwd.plot_eig_embedding()
        mc_fwd.plot_eig_embedding(left=False)

        mc_fwd.compute_approx_rcs(use=1)
        mc_fwd.plot_approx_rcs()

        mc_fwd.compute_lin_probs()
        mc_fwd.plot_lin_probs()

        mc_fwd.compute_lineage_drivers(cluster_key="clusters", use_raw=False)

    def test_default_bwd_pipelne(self):
        adata = _adata.copy()
        vk = VelocityKernel(adata, backward=True).compute_transition_matrix()
        ck = ConnectivityKernel(adata, backward=True).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc_bwd = cr.tl.MarkovChain(final_kernel)

        mc_bwd.compute_partition()

        mc_bwd.compute_eig()
        mc_bwd.plot_eig_embedding()
        mc_bwd.plot_eig_embedding(left=False)

        mc_bwd.compute_approx_rcs(use=1)
        mc_bwd.plot_approx_rcs()

        mc_bwd.compute_lin_probs()
        mc_bwd.plot_lin_probs()

        mc_bwd.compute_lineage_drivers(cluster_key="clusters", use_raw=False)


if __name__ == "__main__":
    unittest.main()
