# -*- coding: utf-8 -*-
from cellrank.tools import MarkovChain
from cellrank.tools.kernels import VelocityKernel, ConnectivityKernel
from typing import Tuple
from anndata import AnnData

import pytest
import scanpy as sc
import scvelo as scv
import cellrank as cr
import numpy as np
import matplotlib


matplotlib.use("Agg")
np.random.seed(42)


def _create_dummy_adata(n_obs: int) -> AnnData:
    adata = scv.datasets.toy_data(n_obs=n_obs)
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=1000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    scv.tl.recover_dynamics(adata)
    scv.tl.velocity(adata, mode="dynamical")
    scv.tl.velocity_graph(adata)
    scv.tl.latent_time(adata)

    adata.uns["iroot"] = 0
    sc.tl.dpt(adata)

    adata.uns["connectivity_variances"] = np.ones((n_obs, n_obs), dtype=np.float64)
    adata.uns["velocity_variances"] = np.ones((n_obs, n_obs), dtype=np.float64)

    return adata


def _create_cellrank_adata(
    n_obs: int, *, backward: bool = False
) -> Tuple[AnnData, MarkovChain]:
    adata = _create_dummy_adata(n_obs)
    sc.tl.paga(adata, groups="clusters")
    try:
        vk = VelocityKernel(adata, backward=backward).compute_transition_matrix()
        ck = ConnectivityKernel(adata, backward=backward).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.MarkovChain(final_kernel)

        mc.compute_partition()
        mc.compute_eig()
        mc.compute_approx_rcs(use=2)  # can fail for small #cells
        mc.compute_lin_probs()
        mc.compute_lineage_drivers(cluster_key="clusters", use_raw=False)

        assert adata is mc.adata
        if backward:
            assert "from_root_cells" in adata.obsm
        else:
            assert "to_final_cells" in adata.obsm
    except:
        mc = None

    return adata, mc


@pytest.fixture
def adata(adata=_create_dummy_adata(50)) -> AnnData:
    return adata.copy()


@pytest.fixture
def adata_large(adata=_create_dummy_adata(200)) -> AnnData:
    return adata.copy()


@pytest.fixture
def adata_mc_fwd(
    adata_mc=_create_cellrank_adata(100, backward=False)
) -> Tuple[AnnData, MarkovChain]:
    adata, mc = adata_mc
    return adata.copy(), mc


@pytest.fixture
def adata_cr(adata_mc=_create_cellrank_adata(100, backward=False)) -> AnnData:
    return adata_mc[0].copy()
