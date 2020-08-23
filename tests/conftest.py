# -*- coding: utf-8 -*-
import os
from typing import Tuple

os.environ["NUMBA_NUM_THREADS"] = "4"

import pytest

import scanpy as sc
from anndata import AnnData

import numpy as np

import matplotlib

import cellrank as cr
from cellrank.tl.kernels import VelocityKernel, ConnectivityKernel
from cellrank.tl._constants import AbsProbKey
from cellrank.tl.estimators import GPCCA, CFLARE

matplotlib.use("Agg")
np.random.seed(42)

_adata_small = sc.read("tests/_ground_truth_adatas/adata_50.h5ad")
_adata_medium = sc.read("tests/_ground_truth_adatas/adata_100.h5ad")
_adata_large = sc.read("tests/_ground_truth_adatas/adata_200.h5ad")


def _create_cflare(*, backward: bool = False) -> Tuple[AnnData, CFLARE]:
    adata = _adata_medium.copy()

    sc.tl.paga(adata, groups="clusters")

    vk = VelocityKernel(adata, backward=backward).compute_transition_matrix()
    ck = ConnectivityKernel(adata, backward=backward).compute_transition_matrix()
    final_kernel = 0.8 * vk + 0.2 * ck

    mc = CFLARE(final_kernel)

    mc.compute_partition()
    mc.compute_eigendecomposition()
    mc.compute_final_states(use=2)
    mc.compute_absorption_probabilities(use_petsc=False)
    mc.compute_lineage_drivers(cluster_key="clusters", use_raw=False)

    assert adata is mc.adata
    if backward:
        assert str(AbsProbKey.BACKWARD) in adata.obsm
    else:
        assert str(AbsProbKey.FORWARD) in adata.obsm
    np.testing.assert_array_almost_equal(mc.absorption_probabilities.sum(1), 1)

    return adata, mc


def _create_gpcca(*, backward: bool = False) -> Tuple[AnnData, GPCCA]:
    adata = _adata_medium.copy()

    sc.tl.paga(adata, groups="clusters")

    vk = VelocityKernel(adata, backward=backward).compute_transition_matrix()
    ck = ConnectivityKernel(adata, backward=backward).compute_transition_matrix()
    final_kernel = 0.8 * vk + 0.2 * ck

    mc = GPCCA(final_kernel)

    mc.compute_partition()
    mc.compute_eigendecomposition()
    mc.compute_schur(method="krylov")
    mc.compute_metastable_states(n_states=2)
    mc.set_final_states_from_metastable_states()
    mc.compute_absorption_probabilities()
    mc.compute_lineage_drivers(cluster_key="clusters", use_raw=False)

    assert adata is mc.adata
    if backward:
        assert str(AbsProbKey.BACKWARD) in adata.obsm
    else:
        assert str(AbsProbKey.FORWARD) in adata.obsm
    np.testing.assert_array_almost_equal(mc.absorption_probabilities.sum(1), 1)

    return adata, mc


@pytest.fixture
def adata() -> AnnData:
    return _adata_small.copy()


@pytest.fixture
def adata_large() -> AnnData:
    return _adata_large.copy()


@pytest.fixture
def adata_cflare_fwd(
    adata_cflare=_create_cflare(backward=False),
) -> Tuple[AnnData, CFLARE]:
    adata, cflare = adata_cflare
    return adata.copy(), cflare


@pytest.fixture
def adata_gpcca_fwd(adata_gpcca=_create_gpcca(backward=False)) -> Tuple[AnnData, GPCCA]:
    adata, gpcca = adata_gpcca
    return adata.copy(), gpcca


@pytest.fixture
def adata_cflare(adata_cflare=_create_cflare(backward=False)) -> AnnData:
    return adata_cflare[0].copy()


@pytest.fixture
def lineage():
    x = cr.tl.Lineage(
        np.array(
            [
                [1.23459664e-01, 1.29965675e-01, 1.92828002e-01, 9.39402664e-01],
                [1.05635239e00, 4.45833459e-01, 2.29080759e00, 1.90132652e00],
                [6.77880737e-02, 4.97556864e-02, 1.18428661e00, 2.02318999e-01],
                [4.87500398e-01, 1.00657498e00, 2.20834882e-02, 5.03008905e-01],
                [6.27190917e00, 7.27864781e00, 1.03978903e00, 1.55903460e01],
                [3.85149269e-01, 3.54765380e-01, 1.77871487e-01, 8.22138648e-02],
                [7.06618729e00, 1.33133671e01, 1.44904591e00, 5.79813391e00],
                [8.18005744e-02, 5.36844933e-01, 1.86646162e00, 2.41141727e00],
                [1.44892035e-01, 2.34036215e-01, 6.32392890e-01, 1.13211403e-02],
                [2.44926466e-01, 2.50293183e-01, 1.77540208e-01, 3.27240144e-01],
            ]
        ),
        names=["foo", "bar", "baz", "quux"],
    )
    return x / x.sum(1)
