import unittest
import cellrank as cr

from cellrank.tools.kernels import VelocityKernel, ConnectivityKernel
from _test_helper import create_dummy_adata


_adata = create_dummy_adata(50)


class MarkovChainTestCase(unittest.TestCase):
    def test_compute_lin_probs_keys_colors(self):
        adata = _adata.copy()
        vk = VelocityKernel(adata).compute_transition_matrix()
        ck = ConnectivityKernel(adata).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc_fwd = cr.tl.MarkovChain(final_kernel)

        mc_fwd.compute_partition()

        mc_fwd.compute_eig()
        mc_fwd.compute_approx_rcs(use=1)
        print(adata)

        mc_fwd.compute_lin_probs()


if __name__ == '__main__':
    unittest.main()
