# -*- coding: utf-8 -*-
import pytest
import scvelo as scv
import numpy as np
import matplotlib


matplotlib.use("Agg")
np.random.seed(42)


def _create_dummy_adata(n_obs: int):
    adata = scv.datasets.toy_data(n_obs=n_obs)
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=1000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    scv.tl.recover_dynamics(adata)
    scv.tl.velocity(adata, mode="dynamical")
    scv.tl.velocity_graph(adata)
    scv.tl.latent_time(adata)
    adata.uns["connectivity_variances"] = np.ones((50, 50), dtype=np.float64)
    adata.uns["velocity_variances"] = np.ones((50, 50), dtype=np.float64)

    return adata


@pytest.fixture
def adata(adata=_create_dummy_adata(50)):
    return adata.copy()


@pytest.fixture
def adata_large(adata=_create_dummy_adata(200)):
    return adata.copy()
