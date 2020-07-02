# -*- coding: utf-8 -*-
"""Utility functions for CellRank tools."""

import os
import warnings
from typing import (
    Any,
    Dict,
    List,
    Tuple,
    Union,
    TypeVar,
    Hashable,
    Iterable,
    Optional,
    Sequence,
)
from itertools import tee, product, combinations

import matplotlib.colors as mcolors

import scanpy as sc
from scanpy import logging as logg
from anndata import AnnData

import numpy as np
import pandas as pd
import networkx as nx
from pandas import Series
from numpy.linalg import norm as d_norm
from scipy.sparse import issparse, spmatrix, csr_matrix
from sklearn.cluster import KMeans
from pandas.api.types import infer_dtype, is_categorical_dtype
from sklearn.neighbors import NearestNeighbors
from scipy.sparse.linalg import norm as s_norm
from cellrank.utils._utils import _get_neighs, _has_neighs, _get_neighs_params
from cellrank.tools._colors import _convert_to_hex_colors, _insert_categorical_colors

ColorLike = TypeVar("ColorLike")
GPCCA = TypeVar("GPCCA")
CFLARE = TypeVar("CFLARE")
EPS = np.finfo(np.float64).eps


def _create_root_final_annotations(
    adata: AnnData,
    fwd: Union[GPCCA, CFLARE],
    bwd: Union[GPCCA, CFLARE],
    final_pref: Optional[str] = "final",
    root_pref: Optional[str] = "root",
    key_added: Optional[str] = "root_final",
):
    """
    Create categorical annotations of both root and final states.

    Params
    ------
    adata
        AnnData object to write to (`.obs[key_added]`)
    fwd
        Estimator object modelling forward process.
    bwd
        Estimator object modelling backward process.
    final_pref, root_pref
        Prefix used in the annotations.
    key_added
        Key added to `adata.obs`.

    Returns
    -------
    None
        Nothing, just writes to AnnData.
    """

    assert not fwd.kernel.backward, "Forward estimator object is in fact backward."
    assert bwd.kernel.backward, "Backward estimator object is in fact forward."

    # get restricted categories and colors
    try:
        # this will work for GPCCA
        cats_final, colors_final = fwd.main_states, fwd.lineage_probabilities.colors
        cats_root, colors_root = bwd.main_states, bwd.lineage_probabilities.colors
    except AttributeError:
        # this works for CFLARE
        cats_final, colors_final = fwd._get_restriction_to_main()
        cats_root, colors_root = bwd._get_restriction_to_main()

    # merge
    cats_merged, colors_merged = _merge_categorical_series(
        cats_final, cats_root, list(colors_final), list(colors_root)
    )

    # adjust the names
    final_names = cats_final.cat.categories
    final_labels = [
        f"{final_pref if key in final_names else root_pref}: {key}"
        for key in cats_merged.cat.categories
    ]
    cats_merged.cat.rename_categories(final_labels, inplace=True)

    # write to AnnData
    adata.obs[key_added] = cats_merged
    adata.uns[f"{key_added}_colors"] = colors_merged


def _process_series(
    series: pd.Series, keys: Optional[List[str]], colors: Optional[np.array] = None
) -> Union[pd.Series, Tuple[pd.Series, List[str]]]:
    """
    Process :class:`pandas.Series` categorical objects.

    Categories in :paramref:`series` are combined/removed according to :paramref:`keys`,
    the same transformation is applied to the corresponding colors.

    Params
    ------
    series
        Input data, must be a pd.series of categorical type.
    keys
        Keys could be e.g. `['cat_1, cat_2', 'cat_4']`. If originally,
        there were 4 categories in `series`, then this would combine the first
        and the second and remove the third. The same would be done to `colors`,
        i.e. the first and second color would be merged (average color), while
        the third would be removed.
    colors
        List of colors which aligns with the order of the categories.

    Returns
    -------
    :class:`pandas.Series`
        Categorical updated annotation. Each cell is assigned to either `NaN`
        or one of updated approximate recurrent classes.
    list
        Color list processed according to keys.
    """

    # determine whether we want to process colors as well
    process_colors = colors is not None

    # if keys is None, just return
    if keys is None:
        if process_colors:
            return series, colors
        return series

    # assert dtype of the series
    if not is_categorical_dtype(series):
        raise TypeError(f"Series must be `categorical`, found `{infer_dtype(series)}`.")

    # initialize a copy of the series object
    series_in = series.copy()
    if process_colors:
        colors_in = np.array(colors.copy())
        if len(colors_in) != len(series_in.cat.categories):
            raise ValueError(
                f"Length of colors ({len(colors_in)}) does not match length of "
                f"categories ({len(series_in.cat.categories)})."
            )
        if not all(mcolors.is_color_like(c) for c in colors_in):
            raise ValueError("Not all colors are color-like.")

    # define a set of keys
    keys_ = {
        tuple(sorted({key.strip(" ") for key in rc.strip(" ,").split(",")}))
        for rc in keys
    }

    # check the `keys` are unique
    overlap = [set(ks) for ks in keys_]
    for c1, c2 in combinations(overlap, 2):
        overlap = c1 & c2
        if overlap:
            raise ValueError(f"Found overlapping keys: `{list(overlap)}`.")

    # check the `keys` are all proper categories
    remaining_cat = [b for a in keys_ for b in a]
    if not np.all(np.in1d(remaining_cat, series_in.cat.categories)):
        raise ValueError(
            "Not all keys are proper categories. Check for spelling mistakes in `keys`."
        )

    # remove cats and colors according to `keys`
    n_remaining = len(remaining_cat)
    removed_cat = list(set(series_in.cat.categories) - set(remaining_cat))
    if process_colors:
        mask = np.in1d(series_in.cat.categories, remaining_cat)
        colors_temp = colors_in[mask].copy()
    series_temp = series_in.cat.remove_categories(removed_cat)

    # loop over all indiv. or combined rc's
    colors_mod = {}
    for cat in keys_:
        # if there are more than two keys in this category, combine them
        if len(cat) > 1:
            new_cat_name = " or ".join(cat)
            mask = np.repeat(False, len(series_temp))
            for key in cat:
                mask = np.logical_or(mask, series_temp == key)
                remaining_cat.remove(key)
            series_temp.cat.add_categories(new_cat_name, inplace=True)
            remaining_cat.append(new_cat_name)
            series_temp[mask] = new_cat_name

            if process_colors:
                # apply the same to the colors array. We just append new colors at the end
                color_mask = np.in1d(series_temp.cat.categories[:n_remaining], cat)
                colors_merge = np.array(colors_temp)[:n_remaining][color_mask]
                colors_mod[new_cat_name] = _compute_mean_color(colors_merge)
        elif process_colors:
            color_mask = np.in1d(series_temp.cat.categories[:n_remaining], cat[0])
            colors_mod[cat[0]] = np.array(colors_temp)[:n_remaining][color_mask][0]

    # Since we have just appended colors at the end, we must now delete the unused ones
    series_temp.cat.remove_unused_categories(inplace=True)
    series_temp.cat.reorder_categories(remaining_cat, inplace=True)

    if process_colors:
        # original colors can still be present, convert to hex
        colors_temp = _convert_to_hex_colors(
            [colors_mod[c] for c in series_temp.cat.categories]
        )
        return series_temp, colors_temp

    return series_temp


def _complex_warning(
    X: np.array, use: Union[list, int, tuple, range], use_imag: bool = False
) -> np.ndarray:
    """
    Check for imaginary components in columns of X specified by `use`.

    Params
    ------
    X
        Matrix containing the eigenvectors
    use
        Selection of columns of `X`
    use_imag
        For eigenvectors that are complex, use real or imaginary part

    Returns
    -------
    class:`numpy.ndarray`
        X_
    """

    complex_mask = np.sum(X.imag != 0, axis=0) > 0
    complex_ixs = np.array(use)[np.where(complex_mask)[0]]
    complex_key = "imaginary" if use_imag else "real"
    if len(complex_ixs) > 0:
        logg.warning(
            f"The eigenvectors with indices {complex_ixs} have an imaginary part. Showing their {complex_key} part."
        )
    X_ = X.real
    if use_imag:
        X_[:, complex_mask] = X.imag[:, complex_mask]

    return X_


def bias_knn(conn, pseudotime, n_neighbors, k=3):
    """
    Palantir Kernel utility function.

    This function takes in symmetric connectivities and a pseudotime and removes edges that point "against" pseudotime,
    in this way creating a directed graph. For each node, it always keeps the closest neighbors, making sure the graph
    remains connected.
    """

    # set a threshold for the neighbors to keep
    k_thresh = np.min([int(np.floor(n_neighbors / k)) - 1, 30])
    conn_biased = conn.copy()

    # loop over rows in the adjacency matrix
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

    return conn_biased


def _vec_mat_corr(X: Union[np.ndarray, spmatrix], y: np.ndarray) -> np.ndarray:
    """
    Compute the correlation between columns in matrix X and a vector y.

    Return NaN for genes which don't vary across cells.

    Params
    ------
    X
        Matrix of `NxM` elements.
    y:
        Vector of `M` elements.

    Returns
    -------
    :class:`numpy.ndarray`
        The computed correlation.
    """

    X_bar, y_std, n = np.array(X.mean(axis=0)).reshape(-1), np.std(y), X.shape[0]
    denom = X.T.dot(y) - n * X_bar * np.mean(y)
    nom = (
        (n - 1) * np.std(X.A, axis=0) * y_std
        if issparse(X)
        else (X.shape[0] - 1) * np.std(X, axis=0) * y_std
    )

    if np.sum(nom == 0) > 0:
        logg.warning(
            f"No variation found in `{np.sum(nom==0)}` genes. Setting correlation for these to `NaN`"
        )

    return denom / nom


def cyto_trace(
    adata: AnnData, layer: str = "Ms", copy: bool = False, use_median: bool = False
) -> Optional[AnnData]:
    """
    Re-implementation of the CytoTrace algorithm by *Gulati et al.* to infer cell plasticity.

    Finds the top 200 genes correlated with #genes/cell and computes their (imputed) mean or median expression.
    For more references, see [Cyto20]_.

    Workflow
    In *scanpy*, take your raw :paramref:`adata` object and run :func:`scvelo.pp.moments` on it. Then run this function.

    Params
    ------
    adata : :class:`anndata.AnnData`
        Annotated data object.
    copy
        Whether to write directly to :paramref:`adata` or to a copy.
    use_median
        If `True`, use *median*, otherwise *mean*.

    Returns
    -------
    :class:`anndata.AnnData` or :class:`NoneType`
        Depending on :paramref:`copy`, either updates :paramref:`adata` or returns a copy.
    """

    # check use_raw and copy
    adata_comp = adata.copy() if copy else adata
    if layer not in adata_comp.layers:
        raise KeyError(f"Compute layer `{layer!r}` first")

    start = logg.info(f"Computing CytoTrace score with `{adata.n_vars}` genes")
    if adata_comp.n_vars < 10000:
        logg.warning("Consider using more genes")

    # compute number of expressed genes per cell
    logg.debug("Computing number of genes expressed per cell")
    num_exp_genes = np.array((adata_comp.X > 0).sum(axis=1)).reshape(-1)
    adata_comp.obs["num_exp_genes"] = num_exp_genes

    # compute correlation with all genes
    logg.debug("Correlating all genes with number of genes expressed per cell")
    gene_corr = _vec_mat_corr(adata_comp.X, num_exp_genes)

    # annotate the top 200 genes in terms of correlation
    logg.debug("Finding the top `200` most correlated genes")
    adata_comp.var["gene_corr"] = gene_corr
    top_200 = adata_comp.var.sort_values(by="gene_corr", ascending=False).index[:200]
    adata_comp.var["correlates"] = False
    adata_comp.var.loc[top_200, "correlates"] = True

    # compute mean/median over top 200 genes, aggregate over genes and shift to [0, 1] range
    logg.debug("Aggregating imputed gene expression")
    corr_mask = adata_comp.var["correlates"]
    imputed_exp = (
        adata_comp[:, corr_mask].X
        if layer == "X"
        else adata_comp[:, corr_mask].layers[layer]
    )
    gcs = np.median(imputed_exp, axis=1) if use_median else np.mean(imputed_exp, axis=1)
    gcs /= np.max(gcs)
    adata_comp.obs["gcs"] = gcs

    logg.info("    Finish", time=start)

    if copy:
        return adata_comp


def _make_cat(
    labels: List[List[Any]], n_states: int, state_names: Sequence[str]
) -> Series:
    """Get categorical from list of lists."""

    labels_new = np.repeat(np.nan, n_states)
    for i, c in enumerate(labels):
        labels_new[c] = i
    labels_new = Series(labels_new, index=state_names, dtype="category")
    labels_new.cat.categories = labels_new.cat.categories.astype("int")

    return labels_new


def _compute_comm_classes(
    A: Union[np.ndarray, spmatrix]
) -> Tuple[List[List[Any]], bool]:
    """Compute communication classes for a graph given by A."""

    di_graph = (
        nx.from_scipy_sparse_matrix(A, create_using=nx.DiGraph)
        if issparse(A)
        else nx.from_numpy_array(A, create_using=nx.DiGraph)
    )

    nx.strongly_connected_components(di_graph)
    comm_classes = sorted(
        nx.strongly_connected_components(di_graph), key=len, reverse=True
    )
    is_irreducible = len(comm_classes) == 1

    return comm_classes, is_irreducible


def _filter_cells(distances: np.ndarray, rc_labels: Series, n_matches_min: int):
    """Filter out some cells that look like transient states based on their neighbors."""

    if not is_categorical_dtype(rc_labels):
        raise TypeError("`categories` must be a categorical variable.")

    # retrieve knn graph
    rows, cols = distances.nonzero()
    cls = rc_labels.cat.categories
    freqs_orig = np.array([np.sum(rc_labels == cl) for cl in cls])

    # loop over cells and check whether they have neighbors from the same class
    for cl in cls:
        cells = np.where(rc_labels == cl)[0]
        for cell in cells:
            own_cl = rc_labels[cell]
            neighbors = cols[rows == cell]
            n_cls = rc_labels[neighbors]
            n_matches = np.sum(np.in1d(n_cls, own_cl))
            if n_matches < n_matches_min:
                rc_labels[cell] = None

    freqs_new = np.array([np.sum(rc_labels == cl) for cl in cls])

    if any(freqs_new / freqs_orig < 0.5):
        print(
            "Warning: consider lowering  'n_matches_min' or "
            "increasing 'n_neighbors_filtering'. This filters out too many cells."
        )

    return rc_labels


def _cluster_X(
    X: Union[np.ndarray, spmatrix],
    n_clusters: int,
    method: str = "kmeans",
    n_neighbors: int = 15,
    resolution: float = 1.0,
) -> List[Any]:
    """Cluster the rows of the matrix X.

    Parameters
    --------
    X
        Data matrix of shape `n_samples x n_features`
    n_clusters
        Number of clusters to use
    method
        Method to use for clustering. Options are ['kmeans', 'louvain', 'leiden']
    n_neighbors
        If using a community-detection based clustering algorithm, number of neighbors for KNN construction
    resolution
        Resolution parameter for ['louvain', 'leiden']

    Returns
    --------
    labels
        List of cluster labels of length `n_samples`

    """
    # make sure data is at least 2D
    if X.ndim == 1:
        X = X[:, None]

    if method == "kmeans":
        kmeans = KMeans(n_clusters=n_clusters).fit(X)
        labels = kmeans.labels_
    elif method in ["louvain", "leiden"]:
        adata_dummy = sc.AnnData(X=X)
        sc.pp.neighbors(adata_dummy, use_rep="X", n_neighbors=n_neighbors)
        if method == "louvain":
            sc.tl.louvain(adata_dummy, resolution=resolution)
        elif method == "leiden":
            sc.tl.leiden(adata_dummy, resolution=resolution)
        labels = adata_dummy.obs[method]
    else:
        raise NotImplementedError(
            f"Invalid method `{method!r}`. Valid options are: `'kmeans', 'louvain'`."
        )

    return list(labels)


def _compute_mean_color(color_list: List[str]) -> str:
    """Compute the mean color."""

    if not all(map(lambda c: mcolors.is_color_like(c), color_list)):
        raise ValueError(f"Not all values are valid colors `{color_list}`.")

    color_list = np.array([mcolors.to_rgb(c) for c in color_list])

    return mcolors.to_hex(np.mean(color_list, axis=0))


def _eigengap(evals: np.ndarray, alpha: float) -> int:
    """
    Compute the eigengap among the top eigenvalues of a matrix.

    Params
    ------
    evals
        Must be real numbers.
    alpha
        Determines how much weight is given to the deviation of an eigenvalue from one.

    Returns
    -------
    int
        Number of eigenvectors to be used.
    """

    gap, eps = evals[:-1] - evals[1:], (1 - evals)[:-1]
    J = gap - alpha * eps

    return int(np.argmax(J))


def partition(
    conn: Union[nx.DiGraph, np.ndarray, spmatrix], sort: bool = True
) -> Tuple[List[List[Any]], List[List[Any]]]:
    """
    Partition a directed graph into its transient and recurrent classes.

    In a directed graph *G*, node *j* is accessible from node *i* if there exists a path from *i* to *j*.
    If *i* is accessible from *j* and the converse holds as well, then *i* and *j* communicate.
    Communication forms and equivalence relation on directed graphs, so every directed graph can be uniquely partitioned
    into its communication classes (also called strongly connected components).

    If *G* describes the state space of a Markov chain, then communication classes are often
    characterized as either recurrent or transient. Intuitively, once the process enters a recurrent class, it will
    never leave it again. See [Tolver16]_ for more formal definition.

    Params
    ------
    conn
        Directed graph to partition.

    Returns
    -------
    (:class:`list`, :class:`list`)
        Recurrent and transient classes respectively.
    """

    start = logg.debug("Partitioning the graph into current and transient classes")

    def partition(g):
        yield from (
            (
                (sorted(scc) if sort else scc),
                all((not nx.has_path(g, s, t) for s, t in product(scc, g.nodes - scc))),
            )
            for scc in nx.strongly_connected_components(g)
        )

    def maybe_sort(iterable):
        return (
            sorted(iterable, key=lambda x: (-len(x), x[0]))
            if sort
            else list(map(list, iterable))
        )

    rec_classes, trans_classes = tee(
        partition(nx.DiGraph(conn) if not isinstance(conn, nx.DiGraph) else conn), 2
    )

    rec_classes = (node for node, is_rec in rec_classes if is_rec)
    trans_classes = (node for node, is_rec in trans_classes if not is_rec)

    logg.debug("    Finish", time=start)

    return maybe_sort(rec_classes), maybe_sort(trans_classes)


def is_connected(c):
    """Check whether the undirected graph encoded by c is connected."""

    G = nx.from_scipy_sparse_matrix(c) if issparse(c) else nx.from_numpy_array(c)

    return nx.is_connected(G)


def is_symmetric(
    matrix: Union[spmatrix, np.ndarray],
    ord: str = "fro",
    eps: float = 1e-4,
    only_check_sparsity_pattern: bool = False,
):
    """Check whether the graph encoded by `matrix` is symmetric."""
    if only_check_sparsity_pattern:
        if issparse(matrix):
            is_sym = len(((matrix != 0) - (matrix != 0).T).data) == 0
        else:
            is_sym = ((matrix != 0) == (matrix != 0).T).all()
    else:
        if issparse(matrix):
            is_sym = s_norm((matrix - matrix.T), ord=ord) < eps
        else:
            is_sym = d_norm((matrix - matrix.T), ord=ord) < eps

    return is_sym


def _subsample_embedding(
    data: Union[np.ndarray, AnnData],
    basis: str = "umap",
    n_dim: int = 2,
    n_grid_points_total: Optional[int] = None,
    n_grid_points_dim: Optional[int] = None,
) -> Tuple[np.ndarray, float]:
    """
    Subsample cells to uniformly cover an embedding.

    If using default parameter settings, this will get very slow for more than 4 embedding dimensions.

    Params
    ------
    data
        Either the embedding or an annotated data object containing an embedding.
    basis
        Key to use to get the embedding from `adata.obsm`.
        Ignored when data is an :class:`np.ndarray`.
    n_dim:
        Number of dimensions in the embedding to use for subsampling.
    n_grid_points_total
        Determines how many gridpoints to use in total.
    n_grid_points_dim
        Determines how many gridpoints to use in each dimension.
        Only one of :paramref:`n_grid_points_total` and :paramref:`n_grid_points_dim` can be specified.

    Returns
    -------
    :class:`np.ndarray`
        Contains the indices of the subsampled cells.
    float
        Euclidian distance between neighboring grid points. Can be used downstream to select a kernel width.
    """

    # check whether we are given an AnnData object
    if isinstance(data, AnnData):
        if f"X_{basis}" not in data.obsm.keys():
            raise ValueError(f"Basis {basis} not found.")
        X_em = data.obsm[f"X_{basis}"][:, :n_dim]
    else:
        X_em = data[:, n_dim]

    # handle grid specification
    if (n_grid_points_total is not None) and (n_grid_points_dim is not None):
        raise ValueError(
            "Can only specify one of `n_grid_points_total` and `n_grid_points_dim`"
        )
    if n_grid_points_total is not None:
        n_grid_points_dim = int((n_grid_points_total) ** (1 / n_dim))
    if (n_grid_points_total is None) and (n_grid_points_dim is None):
        n_grid_points_total = 10 ** (2 + n_dim)
        n_grid_points_dim = int((n_grid_points_total) ** (1 / n_dim))

    # set up the grid
    steps = np.repeat(n_grid_points_dim, n_dim)
    grs = []
    for dim_i in range(n_dim):
        m, M = np.min(X_em[:, dim_i]), np.max(X_em[:, dim_i])
        m = m - 0.025 * np.abs(M - m)
        M = M + 0.025 * np.abs(M - m)
        gr = np.linspace(m, M, steps[dim_i])
        grs.append(gr)
    meshes_tuple = np.meshgrid(*grs)
    gridpoints_coordinates = np.vstack([i.flat for i in meshes_tuple]).T

    # fit nearest neighbors classifier to embedding and determine nearest neighbors of each grid point
    nn = NearestNeighbors()
    nn.fit(X_em)
    dist, ixs = nn.kneighbors(gridpoints_coordinates, 1)

    # reduce this to indices coming from gridpoints covered by the data
    diag_step_dist = np.linalg.norm(
        [grs[dim_i][1] - grs[dim_i][0] for dim_i in range(n_dim)]
    )
    min_dist = diag_step_dist / 2
    ixs = ixs[dist < min_dist]

    cells_ixs = np.unique(ixs)

    return cells_ixs, diag_step_dist


def _gaussian_kernel(
    X: Union[np.ndarray, spmatrix], mu: float = 0, sigma: float = 1
) -> np.ndarray:
    """Compute a Gaussian kernel."""

    if issparse(X):
        G = X.copy()
        G.data = np.exp(-((G.data - mu) ** 2) / (2 * sigma ** 2)) / np.sqrt(
            2 * np.pi * sigma ** 2
        )
    else:
        G = np.exp(-((X - mu) ** 2) / (2 * sigma ** 2)) / np.sqrt(
            2 * np.pi * sigma ** 2
        )

    return G


def _normalize(X: Union[np.ndarray, spmatrix]) -> Union[np.ndarray, spmatrix]:
    """
    Row-normalizes an array to sum to 1.

    Params
    ------
    X
        Array to be normalized.

    Returns
    -------
    :class:`numpy.ndarray` or :class:`scipy.sparse.spmatrx`
        The normalized array.
    """

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if issparse(X):
            X = X.multiply(csr_matrix(1.0 / np.abs(X).sum(1)))
        else:
            X = np.array(X)
            X /= X.sum(1)[:, None]
    return X


def _get_connectivities(
    adata: AnnData, mode: str = "connectivities", n_neighbors: Optional[int] = None
) -> Optional[spmatrix]:
    # utility function, copied from scvelo
    if _has_neighs(adata):
        C = _get_neighs(adata, mode)
        if (
            n_neighbors is not None
            and n_neighbors <= _get_neighs_params(adata)["n_neighbors"]
        ):
            C = (
                _select_connectivities(C, n_neighbors)
                if mode == "connectivities"
                else _select_distances(C, n_neighbors)
            )

        return C.tocsr().astype(np.float32)


def _select_connectivities(
    connectivities: spmatrix, n_neighbors: Optional[int] = None
) -> spmatrix:
    # utility function, copied from scvelo
    C = connectivities.copy()
    n_counts = (C > 0).sum(1).A1 if issparse(C) else (C > 0).sum(1)
    n_neighbors = (
        n_counts.min() if n_neighbors is None else min(n_counts.min(), n_neighbors)
    )
    rows = np.where(n_counts > n_neighbors)[0]
    cumsum_neighs = np.insert(n_counts.cumsum(), 0, 0)
    dat = C.data

    for row in rows:
        n0, n1 = cumsum_neighs[row], cumsum_neighs[row + 1]
        rm_idx = n0 + dat[n0:n1].argsort()[::-1][n_neighbors:]
        dat[rm_idx] = 0
    C.eliminate_zeros()

    return C


def _select_distances(dist, n_neighbors: Optional[int] = None) -> spmatrix:
    # utility funtion, copied from scvelo
    D = dist.copy()
    n_counts = (D > 0).sum(1).A1 if issparse(D) else (D > 0).sum(1)
    n_neighbors = (
        n_counts.min() if n_neighbors is None else min(n_counts.min(), n_neighbors)
    )
    rows = np.where(n_counts > n_neighbors)[0]
    cumsum_neighs = np.insert(n_counts.cumsum(), 0, 0)
    dat = D.data

    for row in rows:
        n0, n1 = cumsum_neighs[row], cumsum_neighs[row + 1]
        rm_idx = n0 + dat[n0:n1].argsort()[n_neighbors:]
        dat[rm_idx] = 0
    D.eliminate_zeros()

    return D


def _maybe_create_dir(dirpath: Union[str, os.PathLike]) -> None:
    """
    Try creating a directory if it does not already exist.

    Params
    ------
    dirpath
        Path of the directory to create.

    Returns
    -------
    None
        Nothing, just creates a directory if it doesn't exist.
    """

    if not os.path.exists(dirpath) or not os.path.isdir(dirpath):
        try:
            os.makedirs(dirpath, exist_ok=True)
        except OSError:
            pass


def save_fig(
    fig, path: Union[str, os.PathLike], make_dir: bool = True, ext: str = "png"
) -> None:
    """
    Save a plot.

    Params
    ------
    fig: :class:`matplotlib.figure.Figure`
        Figure to save.
    path:
        Path where to save the figure.
        If path is relative, save it under `sc.settings.figdir`.
    make_dir:
        Whether to try making the directory if it does not exist.
    ext:
        Extension to use.

    Returns
    -------
    None
        Just saves the plot.
    """

    if os.path.splitext(path)[1] == "":
        path = f"{path}.{ext}"

    if not os.path.isabs(path):
        path = os.path.join(sc.settings.figdir, path)

    if make_dir:
        _maybe_create_dir(os.path.split(path)[0])

    logg.debug(f"Saving figure to `{path!r}`")

    fig.savefig(path, bbox_inches="tight", transparent=True)


def _convert_to_categorical_series(
    rc_classes: Dict[Union[int, str], Iterable[Union[int, str]]], cell_names: List[str]
) -> Series:
    """
    Convert a mapping of recurrent classes to cells to a :class:`pandas.Series`.

    Params
    ------
    rc_classes
        Recurrent classes in the following format: `{'rc_0': ['cell_0', 'cell_1', ...], ...}`.
    cell_names
        List of valid cell names, usually taken from `adata.obs_names`.

    Returns
    -------
    :class:`pandas.Series`
        Categorical series where `NaN` mark cells which do not belong to any recurrent class.
    """

    cnames = set(cell_names)
    mapper, expected_size = {}, 0
    for rc, cells in rc_classes.items():
        if not cells:
            continue
        cells = [c if isinstance(c, str) else cell_names[c] for c in cells]
        rest = set(cells) - cnames
        if rest:
            raise ValueError(f"Invalid cell names: `{list(rest)}`.")
        mapper[str(rc)] = cells
        expected_size += 1

    if len(mapper) != expected_size:
        raise ValueError(
            "All recurrent class labels are being converted to strings, ensure "
            "that there are no conflicting keys, such as `0` and `'0'`."
        )

    rc_labels = Series([np.nan] * len(cell_names), index=cell_names)
    for rc, cells in mapper.items():
        rc_labels[cells] = rc

    return rc_labels.astype("category")


def _merge_categorical_series(
    old: pd.Series,
    new: pd.Series,
    colors_old: Union[List[ColorLike], np.ndarray, Dict[Any, ColorLike]] = None,
    colors_new: Union[List[ColorLike], np.ndarray, Dict[Any, ColorLike]] = None,
    inplace: bool = False,
    color_overwrite: bool = False,
) -> Optional[Union[pd.Series, np.ndarray, Tuple[pd.Series, np.ndarray]]]:
    """
    Update categorical :class:`pandas.Series.` with new information.

    It **can never remove** old categories, only add to the existing ones.
    Optionally, new colors can be created or merged.

    Params
    ------
    old
        Old categories to be updated.
    new
        New categories used to update the old ones.
    colors_old
        Colors associated with old categories.
    colors_new
        Colors associated with new categories.
    color_overwrite
        If `True`, overwrite the old colors with new ones for overlapping categories.
    inplace
        Whether to update :paramref:`old` or create a copy.

    Returns
    -------
    :class:`pandas.Series`
        If paramref:`inplace` is `False`, returns the modified approximate recurrent classes and if
        :paramref:`colors_old` and :paramref:`colors_new` are both `None`.
    :class:`numpy.ndarray`
        If :paramref:`inplace` is `True` and any of :paramref:`colors_old`, :paramref:`colors_new`
        containing the new colors.
    :class:`pandas.Series`, :class:`numpy.ndarray`
        The same as above, but with :paremref:`inplace` is `False`.
    """

    def get_color_mapper(
        series: pd.Series,
        colors: Union[List[ColorLike], np.ndarray, Dict[Any, ColorLike]],
    ):
        if len(series.cat.categories) != len(colors):
            raise ValueError(
                f"Series ({len(series.cat.categories)}) and colors ({len(colors_new)}) differ in length."
            )

        if isinstance(colors, dict):
            if set(series.cat.categories) != set(colors.keys()):
                raise ValueError(
                    "Color mapper and series' categories don't share the keys."
                )
        else:
            colors = dict(zip(series.cat.categories, colors))

        for color in colors.values():
            if not mcolors.is_color_like(color):
                raise ValueError(f"Color `{color}` is not color-like.")

        return colors

    if not is_categorical_dtype(old):
        raise TypeError(
            f"Expected old approx. recurrent classes to be categorical, found "
            f"`{infer_dtype(old)}`."
        )

    if not is_categorical_dtype(new):
        raise TypeError(
            f"Expected new approx. recurrent classes to be categorical, found "
            f"`{infer_dtype(new)}`."
        )

    if (old.index != new.index).any():
        raise ValueError("Index for old and new approx. recurrent classes differ.")

    if not inplace:
        old = old.copy()

    mask = ~new.isna()

    if np.sum(mask) == 0:
        return old if not inplace else None

    old_cats = old.cat.categories
    new_cats = new.cat.categories
    cats_to_add = (
        pd.CategoricalIndex(new.loc[mask]).remove_unused_categories().categories
    )

    if not colors_old and colors_new:
        colors_old = _insert_categorical_colors(
            list(colors_new.values()) if isinstance(colors_new, dict) else colors_new,
            len(old_cats),
        )
    if not colors_new and colors_old:
        colors_new = _insert_categorical_colors(
            list(colors_old.values()) if isinstance(colors_old, dict) else colors_old,
            len(new_cats),
        )

    if colors_old:
        colors_old = get_color_mapper(old, colors_old)
    if colors_new:
        colors_new = get_color_mapper(new, colors_new)

    old.cat.set_categories(old_cats | cats_to_add, inplace=True)
    new.cat.set_categories(old_cats | cats_to_add, inplace=True)

    old.loc[mask] = new.loc[mask]
    old.cat.remove_unused_categories(inplace=True)

    new.cat.set_categories(new_cats, inplace=True)  # return to previous state

    if not colors_old and not colors_new:
        return old if not inplace else None

    colors_merged = (
        {**colors_old, **colors_new}
        if color_overwrite
        else {**colors_new, **colors_old}
    )
    colors_merged = np.array([colors_merged[c] for c in old.cat.categories])

    return (old, colors_merged) if not inplace else colors_merged


def _unique_order_preserving(iterable: Iterable[Hashable]) -> List[Hashable]:
    """Remove items from an iterable while preserving the order."""
    seen = set()
    return [i for i in iterable if i not in seen and not seen.add(i)]


def _generate_random_keys(adata: AnnData, where: str, n: Optional[int] = None):
    def generator():
        return f"CELLRANK_RANDOM_COL_{np.random.randint(2**16)}"

    if n is None:
        n = 1

    where = getattr(adata, where)
    names, seen = [], set(where.keys())

    while len(names) != n:
        name = generator()
        if name not in seen:
            seen.add(name)
            names.append(name)

    return names


def _convert_lineage_name(names: str) -> Tuple[str, ...]:
    sep = "or" if "or" in names else ","
    return tuple(
        sorted({name.strip(" ") for name in names.strip(f" {sep}").split(sep)})
    )


def _long_form_frequencies(
    adata: AnnData,
    query_var: str = "clusters",
    query_var_groups: Optional[Union[Iterable, str]] = None,
    groupby: str = "identifier",
    x_label: Optional[str] = None,
):
    """
    Compute frequencies of a `query_var` over groups defined by `groupby`.

    Params
    --------
    adata
        Annotated Data Matrix
    query_var
        Key from `adata.obs` to a categorical variable whose frequencies with respect to groups defined by `groupby`
        we want to compute
    query_var_groups
        Subset of the categories from `query_var`. These are the categories whose frequencies we are intersted in
    groupby
        Key from `adata.obs`. This defined the categorical variable with respect to which we are computing frequencies
    x_label
        Optional annotation from `adata.obs` that's mapped to `groupby`. Mapping must be unique.

    Returns
    --------
    sub_frequs
        Long-form pandas DataFrame that's convenient for plotting with seaborn
    """

    # input checks
    if query_var not in adata.obs.keys():
        raise ValueError(f"`{query_var}` not found in `adata.obs`")
    if groupby not in adata.obs.keys():
        raise ValueError(f"`{groupby}` not found in `adata.obs`")
    if x_label is None:
        x_label = groupby
    else:
        if x_label not in adata.obs.keys():
            raise ValueError(f"{x_label} not in `adata.obs.keys`")
    if isinstance(query_var_groups, str):
        query_var_groups = [query_var_groups]

    # compute frequencies, and get unique annotations
    frequs = adata.obs.groupby([groupby, query_var]).size()
    samples = adata.obs[groupby].cat.categories.to_numpy()
    ind = adata.obs[query_var].cat.categories.to_numpy()

    # compute relative frequencies (for all categories in query_var)
    rel_frequs = [frequs[ident] / np.sum(frequs[ident]) for ident in samples]
    rel_frequs = pd.DataFrame(rel_frequs, columns=ind, index=samples).fillna(0)
    rel_frequs["x_label"] = np.array(
        [
            np.unique(adata.obs[adata.obs[groupby] == sample][x_label].to_numpy())[0]
            for sample in rel_frequs.index
        ],
        dtype="int",
    )

    # subset to groups of interest and bring into long-form
    sub_frequs = rel_frequs.loc[:, query_var_groups + ["x_label"]]

    return pd.melt(sub_frequs, id_vars="x_label")


def _info_if_obs_key_categorical_present(
    adata: AnnData, key: str, msg_fmt: str
) -> bool:
    if key in adata.obs.keys() and is_categorical_dtype(adata.obs[key]):
        logg.info(msg_fmt.format(key))
        return True


def _info_if_obs_keys_categorical_present(
    adata: AnnData, keys: Iterable[str], msg_fmt: str, warn_once: bool = True
) -> None:
    for key in keys:
        if (
            _info_if_obs_key_categorical_present(adata, key, msg_fmt=msg_fmt)
            and warn_once
        ):
            break


def _one_hot(n, cat: Optional[int] = None) -> np.ndarray:
    """
    One-hot encode cat to a vector of length n.

    If cat is None, return a vector of zeros.
    """

    out = np.zeros(n, dtype=np.bool)
    if cat is not None:
        out[cat] = True

    return out


def _fuzzy_to_discrete(
    a_fuzzy: np.array,
    n_most_likely: int = 10,
    remove_overlap: bool = True,
    raise_threshold: Optional[float] = 0.2,
    check_row_sums: bool = True,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Map fuzzy clustering to discrete clustering.

    Given a fuzzy clustering of `n_samples` samples represented by a matrix `a_fuzzy` of shape
    `(n_samples x n_clusters)` where rows sum to one and indicate cluster membership to each of
    the `c_clusters` clusters, we compute an assignment of a subset of samples to clusters such
    that each cluster is represented by its `n_most_likely` most likely samples. In case a sample
    is assigned more than once, it can either be removed (`remove_overlap = True`) or it can be
    assigned to the cluster it most likely belongs to (`remove_overlap = False`). In case this
    leaves a cluster with less than `raise_threshold x n_most_likely` samples, we raise an exception.
    In case this leaves clusters c_1, ..., c_m with less than `n_most_likely` samples, but more than
    `raise_threshold x n_most_likely` samples, we append c_1, ..., c_m to a list `critical_clusters`,
    which we return.

    We return a boolean matrix `a_discrete` of the same shape as `a_fuzzy`, where 1 in position
    `i, j` indicates that sample `i` is assigned to cluster `j`. Note that we don't assign all samples
    to clusters (most entries in `a_discrete` will be 0) - this is meant to only assign a small
    subset of the samples, which we are most confident in.

    Params
    ------
    a_fuzzy
        Numpy array of shape (`n_samples x n_clusters`) representing a fuzzy clustering. Rows must sum to one.
    n_most_likely
        Number of samples we want to assign to each cluster.
    remove_overlap
        If `True`, remove ambigious samples. Otherwise, assign them to the most likely cluster.
    raise_threshold
        If a cluster is assigned less than `raise_threshold x n_most_likely` samples, raise an
        exception. Set to `None` if you only want to raise if there is an empty cluster.
    check_row_sums
        Check whether rows in `a_fuzzy` sum to one. The one situation where we don't do this is when
        we have selected a couple of main states and we don't want to re-distribute probability mass

    Returns
    -------
    a_discrete
        Boolean matrix of the same shape as `a_fuzzy`, assigning a subset of the samples to clusters.
    critical_clusters
        Array of clusters with less than `n_most_likely` samples assigned.
    """
    # check the inputs
    n_samples, n_clusters = a_fuzzy.shape
    if not isinstance(a_fuzzy, np.ndarray):
        raise TypeError(
            f"Expected `a_fuzzy` to be of type `numpy.ndarray`, got `{type(a_fuzzy).__name__!r}`."
        )
    a_fuzzy = np.asarray(a_fuzzy)  # convert to array from lineage classs, don't copy
    if check_row_sums:
        if n_clusters != 1 and not np.allclose(
            a_fuzzy.sum(1), 1, rtol=1e6 * EPS, atol=1e6 * EPS
        ):
            raise ValueError("Rows in `a_fuzzy` do not sum to one.")
    if n_most_likely > int(n_samples / n_clusters):
        raise ValueError(
            f"You've selected {n_most_likely} cells, please decrease this to at most "
            f"{int(n_samples / n_clusters)} cells for your dataset."
        )

    # initialise
    n_raise = (
        1
        if raise_threshold is None
        else np.max([int(raise_threshold * n_most_likely), 1])
    )
    logg.debug(
        f"DEBUG: Raising an exception if if there are less than `{n_raise}` cells."
    )

    # initially select `n_most_likely` samples per cluster
    sample_assignment = {
        cl: fuzzy_assignment.argpartition(-n_most_likely)[-n_most_likely:]
        for cl, fuzzy_assignment in enumerate(a_fuzzy.T)
    }

    # create the one-hot encoded discrete clustering
    a_discrete = np.zeros(
        a_fuzzy.shape, dtype=np.bool
    )  # don't use `zeros_like` - it also copies the dtype
    for ix in range(n_clusters):
        a_discrete[sample_assignment[ix], ix] = True

    # handle samples assigned to more than one cluster
    critical_samples = np.where(a_discrete.sum(1) > 1)[0]
    for sample_ix in critical_samples:
        if remove_overlap:
            a_discrete[sample_ix, :] = _one_hot(n_clusters)
        else:
            candidate_ixs = np.where(a_discrete[sample_ix, :])[0]
            most_likely_ix = candidate_ixs[
                np.argmax(a_fuzzy[sample_ix, list(a_discrete[sample_ix, :])])
            ]
            a_discrete[sample_ix, :] = _one_hot(n_clusters, most_likely_ix)

    # check how many samples this left for each cluster
    n_samples_per_cluster = a_discrete.sum(0)
    if raise_threshold is not None:
        if (n_samples_per_cluster < n_raise).any():
            min_samples = np.min(n_samples_per_cluster)
            raise ValueError(
                f"Discretizing leads to a cluster with `{min_samples}` samples, less than the threshold which is "
                f"`{n_raise}` samples. Consider recomputing the fuzzy clustering."
            )
    if (n_samples_per_cluster > n_most_likely).any():
        raise ValueError("Assigned more samples than requested.")
    critical_clusters = np.where(n_samples_per_cluster < n_most_likely)[0]

    return a_discrete, critical_clusters


def _series_from_one_hot_matrix(
    a: np.array, index: Optional[Iterable] = None, names: Optional[Iterable] = None
) -> pd.Series:
    """
    Create a pandas Series based on a one-hot encoded matrix.

    Params
    ------
    a
        One-hot encoded membership matrix, of shape (`n_samples x n_clusters`) i.e. a `1` in position `i, j`
        signifies that sample `i` belongs to cluster `j`.
    index
        Index for the Series. Careful, if this is not given, categories are removed when writing to AnnData.

    Returns
    -------
    cluster_series
        Pandas Series, indicating cluster membership for each sample. The dtype of the categories is `str`
        and samples that belong to no cluster are assigned `NaN`.
    """
    n_samples, n_clusters = a.shape
    if not isinstance(a, np.ndarray):
        raise TypeError(
            f"Expected `a` to be of type `numpy.ndarray`, found `{type(a).__name__!r}`."
        )
    a = np.asarray(a)  # change the type in case a lineage object was passed.
    if a.dtype != np.bool:
        raise TypeError(
            f"Expected `a`'s elements to be boolean, found `{a.dtype.name}`."
        )

    if not np.all(a.sum(axis=1) <= 1):
        raise ValueError("Not all items are one-hot encoded or empty.")
    if (a.sum(0) == 0).any():
        logg.warning(f"Detected {np.sum((a.sum(0) == 0))} empty categorie(s) ")

    if index is None:
        index = range(n_samples)
    if names is not None:
        if len(names) != n_clusters:
            raise ValueError(
                f"Shape mismatch, length of `names` is `{len(names)}`, but `n_clusters` = {n_clusters}."
            )
    else:
        names = np.arange(n_clusters).astype("str")

    target_series = pd.Series(index=index, dtype="category")
    for (vec, name) in zip(a.T, names):
        target_series.cat.add_categories(name, inplace=True)
        target_series[np.where(vec)[0]] = name

    return target_series


def _colors_in_order(
    adata: AnnData,
    clusters: Optional[Iterable[str]] = None,
    cluster_key: str = "clusters",
):
    """Get list of colors from AnnData in defined order.

    This will extract a list of colors from `adata.uns[cluster_key]` corresponding to the `clusters`, in the
    order defined by the `clusters`

    Parameters
    --------
    adata
        Annotated data matrix
    clusters
        Subset of the clusters we want the color for. Must be a subset of `adata.obs[cluster_key].cat.categories`
    cluster_key
        Key from `adata.obs`

    Returns
    --------
    color_list
        List of colors in order defined by `clusters`
    """
    assert (
        cluster_key in adata.obs.keys()
    ), f"Could not find {cluster_key} in `adata.obs`"
    assert np.in1d(
        clusters, adata.obs[cluster_key].cat.categories
    ).all(), "Not all `clusters` found"
    assert (
        f"{cluster_key}_colors" in adata.uns.keys()
    ), f"No colors associated to {cluster_key} in `adata.uns`"

    if clusters is None:
        clusters = adata.obs[cluster_key].cat.categories

    color_list = []
    all_clusters = adata.obs[cluster_key].cat.categories
    for cl in clusters:
        mask = np.in1d(all_clusters, cl)
        color_list.append(adata.uns[f"{cluster_key}_colors"][mask][0])

    return color_list
