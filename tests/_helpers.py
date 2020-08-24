# -*- coding: utf-8 -*-
import os
from typing import Tuple, Union, Optional
from pathlib import Path

import pytest
from PIL import Image

import scanpy as sc
import scvelo as scv
from scanpy import logging as logg
from anndata import AnnData

import numpy as np
import pandas as pd
from sklearn.svm import SVR
from scipy.sparse import spdiags, issparse, csr_matrix
from scipy.sparse.linalg import norm

import cellrank as cr
from cellrank.tl._utils import _normalize
from cellrank.ul._utils import _get_neighs, _get_neighs_params
from cellrank.tl.kernels import VelocityKernel, ConnectivityKernel
from cellrank.tl._constants import Direction, _transition


def _jax_not_installed() -> bool:
    try:
        import jax
        import jaxlib

        return False
    except ImportError:
        return True


def bias_knn(conn, pseudotime, n_neighbors, k=3):
    k_thresh = np.min([int(np.floor(n_neighbors / k)) - 1, 30])
    conn_biased = conn.copy()

    # check whether the original graph was connected
    assert _is_connected(conn), "The underlying KNN graph is disconnected."

    for i in range(conn.shape[0]):
        # get indices, values and current pseudo t
        row_data = conn[i, :].data
        row_ixs = conn[i, :].indices
        current_t = pseudotime[i]

        # get the 'candidates' - ixs of nodes not in the k_thresh closest neighbors
        p = np.flip(np.argsort(row_data))
        sorted_ixs = row_ixs[p]
        cand_ixs = sorted_ixs[k_thresh:]

        # compare pseudotimes and set indices to zero
        cand_t = pseudotime[cand_ixs]
        rem_ixs = cand_ixs[cand_t < current_t]
        conn_biased[i, rem_ixs] = 0

    conn_biased.eliminate_zeros()

    # check whether the biased graph is still connected
    assert _is_connected(conn_biased), "The biased KNN graph has become disconnected."

    return conn_biased


def density_normalization(velo_graph, trans_graph):
    # function copied from scanpy
    q = np.asarray(trans_graph.sum(axis=0))
    if not issparse(trans_graph):
        Q = np.diag(1.0 / q)
    else:
        Q = spdiags(1.0 / q, 0, trans_graph.shape[0], trans_graph.shape[0])
    velo_graph = Q @ velo_graph @ Q

    return velo_graph


def _is_connected(c) -> bool:
    import networkx as nx
    from scipy.sparse import issparse

    G = nx.from_scipy_sparse_matrix(c) if issparse(c) else nx.from_numpy_array(c)

    return nx.is_connected(G)


def create_kernels(
    adata: AnnData,
    velocity_variances: Optional[str] = None,
    connectivity_variances: Optional[str] = None,
) -> Tuple[VelocityKernel, ConnectivityKernel]:
    vk = VelocityKernel(adata)
    vk._mat_scaler = adata.uns.get(
        velocity_variances, np.random.normal(size=(adata.n_obs, adata.n_obs))
    )

    ck = ConnectivityKernel(adata)
    ck._mat_scaler = adata.uns.get(
        connectivity_variances, np.random.normal(size=(adata.n_obs, adata.n_obs))
    )

    vk._transition_matrix = csr_matrix(np.eye(adata.n_obs))
    ck._transition_matrix = np.eye(adata.n_obs, k=1) / 2 + np.eye(adata.n_obs) / 2
    ck._transition_matrix[-1, -1] = 1
    ck._transition_matrix = csr_matrix(ck._transition_matrix)

    np.testing.assert_allclose(
        np.sum(ck._transition_matrix.A, axis=1), 1
    )  # sanity check

    return vk, ck


def create_model(adata: AnnData) -> cr.ul.models.SKLearnModel:
    return cr.ul.models.SKLearnModel(adata, SVR(kernel="rbf"))


def resize_images_to_same_sizes(
    expected_image_path: Union[str, Path],
    actual_image_path: Union[str, Path],
    kind: str = "actual_to_expected",
) -> None:
    if not os.path.isfile(actual_image_path):
        raise OSError(f"Actual image path `{actual_image_path!r}` does not exist.")
    if not os.path.isfile(expected_image_path):
        raise OSError(f"Expected image path `{expected_image_path!r}` does not exist.")
    expected_image = Image.open(expected_image_path)
    actual_image = Image.open(actual_image_path)
    if expected_image.size != actual_image.size:
        if kind == "actual_to_expected":
            actual_image.resize(expected_image.size).save(actual_image_path)
        elif kind == "expected_to_actual":
            expected_image.resize(actual_image.size).save(expected_image)
        else:
            raise ValueError(
                f"Invalid kind of conversion `{kind!r}`. Valid options are `'actual_to_expected'`, `'expected_to_actual'`."
            )


def assert_array_nan_equal(
    actual: Union[np.ndarray, pd.Series], expected: Union[np.ndarray, pd.Series]
) -> None:
    """
    Test is 2 arrays or :class:`pandas.Series` are equal.

    Params
    ------
    actual
        The actual data.
    expected
        The expected result.

    Returns
    -------
    None
        Nothing, but raises an exception if arrays are not equal, including the locations of NaN values.
    """

    mask1 = ~(pd.isnull(actual) if isinstance(actual, pd.Series) else np.isnan(actual))
    mask2 = ~(
        pd.isnull(expected) if isinstance(expected, pd.Series) else np.isnan(expected)
    )
    np.testing.assert_array_equal(np.where(mask1), np.where(mask2))
    np.testing.assert_array_equal(actual[mask1], expected[mask2])


def random_transition_matrix(n: int) -> np.ndarray:
    """
    Parameters
    ----------
    n
        Number of states.

    Returns
    -------
    :class:`numpy.ndarray`
        Row-normalized transition matrix.
    """

    x = np.abs(np.random.normal(size=(n, n)))
    rsum = x.sum(axis=1)
    return x / rsum[:, np.newaxis]


def _create_dummy_adata(n_obs: int) -> AnnData:
    """
    Create a testing :class:`anndata.AnnData` object.

    Call this function to regenerate the above objects.

    Params
    ------
    n_obs
        Number of cells.

    Returns
    -------
    :class:`anndata.AnnData`
        The created adata object.
    """

    np.random.seed(42)
    adata = scv.datasets.toy_data(n_obs=n_obs)
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=1000)
    adata.raw = adata[:, 42 : 42 + 50].copy()
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    scv.tl.recover_dynamics(adata)
    scv.tl.velocity(adata, mode="dynamical")
    scv.tl.velocity_graph(adata, mode_neighbors="connectivities")
    scv.tl.latent_time(adata)

    adata.uns["iroot"] = 0
    sc.tl.dpt(adata)

    adata.uns["connectivity_variances"] = np.ones((n_obs, n_obs), dtype=np.float64)
    adata.uns["velocity_variances"] = np.ones((n_obs, n_obs), dtype=np.float64)

    sc.write(f"tests/_ground_truth_adatas/adata_{n_obs}.h5ad", adata)

    return adata


jax_not_installed_skip = pytest.mark.skipif(
    _jax_not_installed(), reason="JAX is not installed."
)


if __name__ == "__main__":
    for size in [50, 100, 200]:
        _ = _create_dummy_adata(size)
