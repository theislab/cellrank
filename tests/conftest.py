# -*- coding: utf-8 -*-
from typing import Tuple

import numpy as np
import pytest

import matplotlib

import scanpy as sc
import scvelo as scv
from anndata import AnnData

import cellrank as cr
from cellrank.tools import CFLARE
from cellrank.tools.kernels import VelocityKernel, ConnectivityKernel
from cellrank.tools._constants import LinKey

matplotlib.use("Agg")
np.random.seed(42)


def _create_dummy_adata(n_obs: int) -> AnnData:
    np.random.seed(42)
    adata = scv.datasets.toy_data(n_obs=n_obs)
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=1000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    scv.tl.recover_dynamics(adata)
    scv.tl.velocity(adata, mode="dynamical")
    scv.tl.velocity_graph(adata)
    scv.tl.latent_time(adata)

    adata.uns["iroot"] = 0
    sc.pp.neighbors(adata, n_pcs=15)
    sc.tl.dpt(adata)

    adata.uns["connectivity_variances"] = np.ones((n_obs, n_obs), dtype=np.float64)
    adata.uns["velocity_variances"] = np.ones((n_obs, n_obs), dtype=np.float64)

    return adata


def _create_cflare(n_obs: int, *, backward: bool = False) -> Tuple[AnnData, CFLARE]:
    adata = _create_dummy_adata(n_obs)
    sc.tl.paga(adata, groups="clusters")
    try:
        vk = VelocityKernel(adata, backward=backward).compute_transition_matrix()
        ck = ConnectivityKernel(adata, backward=backward).compute_transition_matrix()
        final_kernel = 0.8 * vk + 0.2 * ck

        mc = cr.tl.CFLARE(final_kernel)

        mc.compute_partition()
        mc.compute_eig()
        mc.compute_metastable_states(use=2)  # can fail for small #cells
        mc.compute_lin_probs()
        mc.compute_lineage_drivers(cluster_key="clusters", use_raw=False)

        assert adata is mc.adata
        if backward:
            assert str(LinKey.BACKWARD) in adata.obsm
        else:
            assert str(LinKey.FORWARD) in adata.obsm
    except:
        mc = None

    return adata, mc


def _create_gpcca(n_obs: int, *, backward: bool = False) -> Tuple[AnnData, CFLARE]:
    adata = _create_dummy_adata(n_obs)
    sc.tl.paga(adata, groups="clusters")
    vk = VelocityKernel(adata, backward=backward).compute_transition_matrix()
    ck = ConnectivityKernel(adata, backward=backward).compute_transition_matrix()
    final_kernel = 0.8 * vk + 0.2 * ck

    mc = cr.tl.GPCCA(final_kernel)

    mc.compute_partition()
    mc.compute_eig()
    mc.compute_schur(method="brandts")
    mc.compute_metastable_states(n_states=2)
    mc.set_main_states()
    mc.compute_lineage_drivers(cluster_key="clusters", use_raw=False)

    assert adata is mc.adata
    if backward:
        assert str(LinKey.BACKWARD) in adata.obsm
    else:
        assert str(LinKey.FORWARD) in adata.obsm

    return adata, mc


@pytest.fixture
def adata(adata=_create_dummy_adata(50)) -> AnnData:
    return adata.copy()


@pytest.fixture
def adata_large(adata=_create_dummy_adata(200)) -> AnnData:
    return adata.copy()


@pytest.fixture
def adata_cflare_fwd(
    adata_cflare=_create_cflare(100, backward=False)
) -> Tuple[AnnData, CFLARE]:
    adata, cflare = adata_cflare
    return adata.copy(), cflare


@pytest.fixture
def adata_gpcca_fwd(
    adata_gpcca=_create_gpcca(100, backward=False)
) -> Tuple[AnnData, CFLARE]:
    adata, gpcca = adata_gpcca
    return adata.copy(), gpcca


@pytest.fixture
def adata_cflare(adata_cflare=_create_cflare(100, backward=False)) -> AnnData:
    return adata_cflare[0].copy()


@pytest.fixture
def lineage():
    return cr.tl.Lineage(
        np.array(
            [
                [1.23459664e-01, 1.29965675e-01, -1.92828002e-01, 9.39402664e-01],
                [1.05635239e00, -4.45833459e-01, 2.29080759e00, -1.90132652e00],
                [6.77880737e-02, -4.97556864e-02, 1.18428661e00, -2.02318999e-01],
                [-4.87500398e-01, 1.00657498e00, -2.20834882e-02, 5.03008905e-01],
                [-6.27190917e00, -7.27864781e00, -1.03978903e00, 1.55903460e01],
                [3.85149269e-01, 3.54765380e-01, 1.77871487e-01, 8.22138648e-02],
                [7.06618729e00, -1.33133671e01, 1.44904591e00, 5.79813391e00],
                [-8.18005744e-02, 5.36844933e-01, -1.86646162e00, 2.41141727e00],
                [1.44892035e-01, 2.34036215e-01, 6.32392890e-01, -1.13211403e-02],
                [2.44926466e-01, 2.50293183e-01, 1.77540208e-01, 3.27240144e-01],
            ]
        ),
        names=["foo", "bar", "baz", "quux"],
    )
