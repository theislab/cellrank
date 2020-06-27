# -*- coding: utf-8 -*-
import os
from typing import Tuple, Union, Optional
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image
from sklearn.svm import SVR
from scipy.sparse import spdiags, issparse, csr_matrix
from scipy.sparse.linalg import norm

from scanpy import logging as logg
from anndata import AnnData

import cellrank as cr
from cellrank.tools._utils import _normalize
from cellrank.utils._utils import _get_neighs, _get_neighs_params
from cellrank.tools.kernels import VelocityKernel, ConnectivityKernel
from cellrank.tools._constants import Direction, _transition


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


def transition_matrix(
    adata: AnnData,
    vkey: str = "velocity",
    backward: bool = False,
    self_transitions: Optional[str] = None,
    sigma_corr: Optional[float] = None,
    diff_kernel: Optional[str] = None,
    weight_diffusion: float = 0.2,
    density_normalize: bool = True,
    backward_mode: str = "transpose",
    inplace: bool = True,
) -> csr_matrix:
    """
    Computes transition probabilities from velocity graph.

    THIS FUNCTION HAS BEEN DEPRECATED.
    Interact with kernels via the Kernel class or via cellrank.tools_transition_matrix.transition_matrix

    Employs ideas of both scvelo as well as velocyto.

    Parameters
    --------
    adata : :class:`anndata.AnnData`
        Annotated Data Matrix
    vkey
        Name of the velocity estimates to be used
    backward
        Whether to use the transition matrix to push forward (`False`) or to pull backward (`True`)
    self_transitions
        How to fill the diagonal. Can be either 'velocyto' or 'scvelo'. Two diffent
        heuristics are used. Can prevent dividing by zero in unlucky sitatuations for the
        reverse process
    sigma_corr
        Kernel width for exp kernel to be used to compute transition probabilities
        from the velocity graph. If None, the median cosine correlation of all
        potisive cosine correlations will be used.
    diff_kernel
        Whether to multiply the velocity connectivities with transcriptomic distances to make them more robust.
        Options are ('sum', 'mult', 'both')
    weight_diffusion
        Relative weight given to the diffusion kernel. Must be in [0, 1]. Only matters when using 'sum' or 'both'
        for the diffusion kernel.
    density_normalize
        Whether to use the transcriptomic KNN graph for density normalization as performed in scanpy when
        computing diffusion maps
    backward_mode
        Options are ['transpose', 'negate'].
    inplace
        If True, adds to adata. Otherwise returns.

    Returns
    --------
    T: :class:`scipy.sparse.csr_matrix`
        Transition matrix
    """
    logg.info("Computing transition probability from velocity graph")

    from datetime import datetime

    print(datetime.now())

    # get the direction of the process
    direction = Direction.BACKWARD if backward else Direction.FORWARD

    # get the velocity correlations
    if (vkey + "_graph" not in adata.uns.keys()) or (
        vkey + "_graph_neg" not in adata.uns.keys()
    ):
        raise ValueError(
            "You need to run `tl.velocity_graph` first to compute cosine correlations"
        )
    velo_corr, velo_corr_neg = (
        csr_matrix(adata.uns[vkey + "_graph"]).copy(),
        csr_matrix(adata.uns[vkey + "_graph_neg"]).copy(),
    )
    velo_corr_comb_ = (velo_corr + velo_corr_neg).astype(np.float64)
    if backward:
        if backward_mode == "negate":
            velo_corr_comb = velo_corr_comb_.multiply(-1)
        elif backward_mode == "transpose":
            velo_corr_comb = velo_corr_comb_.T
        else:
            raise ValueError(f"Unknown backward_mode `{backward_mode}`.")
    else:
        velo_corr_comb = velo_corr_comb_
    med_corr = np.median(np.abs(velo_corr_comb.data))

    # compute the raw transition matrix. At the moment, this is just an exponential kernel
    logg.debug("DEBUG: Computing the raw transition matrix")
    if sigma_corr is None:
        sigma_corr = 1 / med_corr
    velo_graph = velo_corr_comb.copy()
    velo_graph.data = np.exp(velo_graph.data * sigma_corr)

    # should I row-_normalize the transcriptomic connectivities?
    if diff_kernel is not None or density_normalize:
        params = _get_neighs_params(adata)
        logg.debug(
            f'DEBUG: Using KNN graph computed in basis {params.get("use_rep", "Unknown")!r} '
            'with {params["n_neighbors"]} neighbors'
        )
        trans_graph = _get_neighs(adata, "connectivities")
        dev = norm((trans_graph - trans_graph.T), ord="fro")
        if dev > 1e-4:
            logg.warning("KNN base graph not symmetric, `dev={dev}`")

    # KNN smoothing
    if diff_kernel is not None:
        logg.debug("DEBUG: Smoothing KNN graph with diffusion kernel")
        velo_graph = _knn_smooth(diff_kernel, velo_graph, trans_graph, weight_diffusion)
    # return velo_graph

    # set the diagonal elements. This is important especially for the backwards direction
    logg.debug("DEBUG: Setting diagonal elements")
    velo_graph = _self_loops(self_transitions, velo_graph)

    # density normalisation - taken from scanpy
    if density_normalize:
        logg.debug("DEBUG: Density correcting the velocity graph")
        velo_graph = density_normalization(velo_graph, trans_graph)

    # normalize
    T = _normalize(velo_graph)

    if not inplace:
        logg.info("Computed transition matrix")
        return T

    if _transition(direction) in adata.uns.keys():
        logg.warning(
            f"`.uns` already contains a field `{_transition(direction)!r}`. Overwriting"
        )

    params = {
        "backward": backward,
        "self_transitions": self_transitions,
        "sigma_corr": np.round(sigma_corr, 3),
        "diff_kernel": diff_kernel,
        "weight_diffusion": weight_diffusion,
        "density_normalize": density_normalize,
    }

    adata.uns[_transition(direction)] = {"T": T, "params": params}
    logg.info(
        f"Computed transition matrix and added the key `{_transition(direction)!r}` to `adata.uns`"
    )


def _knn_smooth(diff_kernel, velo_graph, trans_graph, weight_diffusion):
    # utility function for combining KNN kernel and velocity kernel
    assert weight_diffusion >= 0, "Weight diffusion must be non-negative."
    assert weight_diffusion <= 1, "Weight diffusion must be <= 1."

    # this is necessary because I don't want to normalize this graph (density correction)
    G_sim = trans_graph.copy()

    if diff_kernel == "mult":
        logg.debug("DEBUG: Using a multiplicative diffusion kernel")
        # element wise multiplication
        velo_graph = velo_graph.multiply(G_sim)
    elif diff_kernel == "sum":
        logg.debug("DEBUG: Using an additive diffusion kernel")
        # G_sim  = G_sim.multiply(velo_graph>0)
        velo_graph, trans_graph = _normalize(velo_graph), _normalize(G_sim)
        velo_graph = (
            1 - weight_diffusion
        ) * velo_graph + weight_diffusion * trans_graph
    elif diff_kernel == "both":
        logg.debug(
            "DEBUG: Using first a multiplicative and then an additive diffusion kernel"
        )
        G_sim = G_sim.multiply(velo_graph > 0)
        velo_graph = velo_graph.multiply(G_sim)
        velo_graph, trans_grap = _normalize(velo_graph), _normalize(G_sim)
        velo_graph = (1 - weight_diffusion) * velo_graph + weight_diffusion * G_sim
    else:
        raise ValueError(
            f"Invalid kernel type `{diff_kernel}`. Valid options are: `'mult', 'sum', 'both'`."
        )

    return velo_graph


def _self_loops(self_transitions, velo_graph):
    # set the diagonal elements.
    if self_transitions is not None:
        logg.info(f"Self transitions using {self_transitions!r}")
    if self_transitions == "scvelo":
        confidence = velo_graph.max(1).A.flatten()
        ub = np.percentile(confidence, 98)
        self_prob = np.clip(ub - confidence, 0, 1)
        velo_graph.setdiag(self_prob)
    if self_transitions == "velocyto":
        self_prob = velo_graph.max(1).A.flatten()
        velo_graph.setdiag(self_prob)

    return velo_graph


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
    from scipy.sparse import issparse
    import networkx as nx

    G = nx.from_scipy_sparse_matrix(c) if issparse(c) else nx.from_numpy_array(c)

    return nx.is_connected(G)


def create_kernels(
    adata: AnnData,
    var_key_connectivities: str = "connectivity_variances",
    var_key_velocities: str = "velocity_variances",
) -> Tuple[VelocityKernel, ConnectivityKernel]:
    vk = VelocityKernel(adata, var_key=var_key_velocities)
    ck = ConnectivityKernel(adata, var_key=var_key_connectivities)
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
