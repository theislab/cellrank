# -*- coding: utf-8 -*-
"""
Utility functions for the cellrank tools
"""

from scipy.sparse.linalg import norm as s_norm
from numpy.linalg import norm as d_norm
from itertools import product, tee
from typing import Optional, Any, Union, Tuple, List, Sequence, Dict, Iterable

import os
import warnings
import pandas as pd
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import networkx as nx
import numpy as np
import scanpy as sc

from anndata import AnnData
from pandas import Series, DataFrame, to_numeric
from pandas.api.types import is_categorical_dtype, infer_dtype
from scanpy import logging as logg
from scipy.sparse import csr_matrix, spmatrix
from scipy.stats import entropy
from scipy.sparse import issparse
from sklearn.cluster import KMeans
from sklearn.neighbors import NearestNeighbors

from cellrank.utils._utils import has_neighs, get_neighs, get_neighs_params

from typing import List, Optional
from itertools import combinations


def _map_names_and_colors(
    series_reference: Series,
    series_query: Series,
    colors_reference: Optional[np.array] = None,
    en_cutoff: Optional[float] = None,
) -> Tuple[Series, List[Any]]:
    """
    Utility function to map annotations and colors from one series to another

    Params
    ------
    series_reference
        Series object with categorical annotations
    series_query
        Series for which we would like to query the category names
    colors_reference
        If given, colors for the query categories are pulled from this color array.
    en_cutoff
        In case of a non-perfect overlap between categories of the two series, this decides when to label a categoy in
        the query as 'Unknown'

    Returns
    -------
    :class:`pandas.Series`, :class:`list`
        Series with updated category names and a corresponding array of colors
    """

    # create dataframe to store approx_rc associtation with ts_clusters
    approx_rc_clusters = series_query.cat.categories
    ts_clusters = series_reference.cat.categories
    rc_df = DataFrame(None, index=approx_rc_clusters, columns=ts_clusters)

    # populate the df - compute the overlap
    for cl in approx_rc_clusters:
        row = [
            np.sum(series_reference.loc[np.array(series_query == cl)] == key)
            for key in ts_clusters
        ]
        rc_df.loc[cl] = row
    rc_df = rc_df.apply(to_numeric)

    # label endpoints and add uncertainty through entropy
    rc_df["entropy"] = entropy(rc_df.T)
    rc_df["name"] = rc_df.T.idxmax()

    # add cluster colors
    # clusters_colors = self.lin_probs.colors  # _convert_to_hex_colors(self._adata.uns[_colors(cluster_key)])
    clusters_colors = _convert_to_hex_colors(colors_reference)
    names = rc_df["name"]
    lin_colors = []
    for name in names:
        mask = ts_clusters == name
        color = np.array(clusters_colors)[mask][0]
        lin_colors.append(color)
    rc_df["color"] = lin_colors

    # make series_query unique
    names = Series(rc_df["name"], dtype="category")
    lin_colors = Series(lin_colors, dtype="category")  # colors must be hex
    frequ = {key: np.sum(names == key) for key in names.cat.categories}

    # Create unique names by adding suffixes "..._1, ..._2" etc and unique colors by shifting the original color
    names_new = np.array(names.copy())
    lin_colors_new = np.array(lin_colors.copy())

    for key, value in frequ.items():
        if value == 1:
            continue  # already unique, skip

        # deal with non-unique names
        suffix = list(np.arange(1, value + 1).astype("str"))
        unique_names = [f"{key}_{rep}" for rep in suffix]
        names_new[names == key] = unique_names

        color = rc_df[rc_df["name"] == key]["color"].values[0]
        shifted_colors = _create_colors(color, value, saturation_range=None)
        lin_colors_new[lin_colors == color] = shifted_colors

    rc_df["name"] = names_new
    rc_df["color"] = lin_colors_new

    # issue a warning for rcs with high entropy
    if en_cutoff is not None:
        critical_rcs = list(rc_df.loc[rc_df["entropy"] > en_cutoff, "name"].values)
        if len(critical_rcs) > 0:
            logg.warning(
                f"The following groups could not be mapped uniquely: `{critical_rcs}`"
            )

    return rc_df["name"], list(rc_df["color"])


def _process_series(series: pd.Series, keys: List, colors: Optional[np.array] = None):
    """
    Utility function to process pd.Series categorical objects

    Categories in `series` are combined/removed according to `keys`,
    the same transformation is applied to the corresponding colors.

    Parameters
    --------
    series
        Input data, must be a pd.series of categorical type
    keys
        Keys could be e.g. ['cat_1, cat_2', 'cat_4']. If originally,
        there were 4 categories in `series`, then this would combine the first
        and the second and remove the third. The same would be done to `colors`,
        i.e. the first and second color would be merged (average color), while
        the third would be removed.
    colors
        List of colors which aligns with the order of the categories

    Returns
    --------
    :class:`pandas.Series`
        Categorical updated annotation. Each cell is assigned to either `NaN`
        or one of updated approximate recurrent classes.
    list
        Color list processed according to keys
    """

    if colors is not None:
        process_colors = True
    else:
        process_colors = False

    # if keys is None, just return
    if keys is None:
        if process_colors:
            return series, colors
        else:
            return series

    # assert dtype of the series
    assert series.dtype.name == "category", "Type of the series is not categorical"

    # initialize a copy of the series object
    series_in = series.copy()
    if process_colors:
        colors_in = np.array(colors.copy())
        assert len(colors_in) == len(
            series_in.cat.categories
        ), "Length of colors does not match length of categories"

    # define a set of keys
    keys_ = {tuple((key.strip() for key in rc.strip(" ,").split(","))) for rc in keys}

    # check the `keys` are unique
    overlap = [set(ks) for ks in keys_]
    for c1, c2 in combinations(overlap, 2):
        overlap = c1 & c2
        if overlap:
            raise ValueError(f"Found overlapping keys: `{list(overlap)}`.")

    # check the `keys` are all proper categories
    remaining_cat = [b for a in keys_ for b in a]
    assert all(np.in1d(remaining_cat, series_in.cat.categories)), (
        "Not all keys are proper categories. " "Check for spelling mistakes in `keys`."
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
        colors_temp = [colors_mod[c] for c in series_temp.cat.categories]
        return series_temp, colors_temp
    else:
        return series_temp


def _complex_warning(
    X: np.array, use: Union[list, int, tuple, range], use_imag: bool = False
):
    """
    Check for imaginary components in columns of X specified by `use`

    Params
    --------
    X
        Matrix containing the eigenvectors
    use
        Selection of columns of `X`
    use_imag
        For eigenvectors that are complex, use real or imaginary part

    Returns
    --------
    X_
        np.array
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
    Utility function for the Palantir Kernel.

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
    Computes the correlation between columns in matrix X and a vector y

    Returns NaN for genes which don't vary across cells

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
    """
    Get categorical from list of lists.
    """

    labels_new = np.repeat(np.nan, n_states)
    for i, c in enumerate(labels):
        labels_new[c] = i
    labels_new = Series(labels_new, index=state_names, dtype="category")
    labels_new.cat.categories = labels_new.cat.categories.astype("int")

    return labels_new


def _compute_comm_classes(
    A: Union[np.ndarray, spmatrix]
) -> Tuple[List[List[Any]], bool]:
    """
    Utility function to compute communication classes for a graph given by A.
    """

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
    """
    Utility function which filters out some cells that look like transient states based on their neighbors.
    """

    if not is_categorical_dtype(rc_labels):
        raise TypeError("`rc_labels` must be a categorical variable.")

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
    method: str,
    n_clusters_kmeans: int,
    percentile: Optional[float],
    use: Union[Tuple[int], List[int]],
    n_neighbors_louvain: int,
    resolution_louvain: float,
) -> List[Any]:
    """
    Utility function which clusters the rows of the matrix X.
    """

    if method == "kmeans":
        if n_clusters_kmeans is None:
            if percentile is not None:
                n_clusters_kmeans = len(use)
            else:
                n_clusters_kmeans = len(use) + 1
        kmeans = KMeans(n_clusters=n_clusters_kmeans).fit(X)
        labels = kmeans.labels_
    elif method == "louvain":
        if len(use) <= 1:
            raise ValueError(
                f"Number of eigenvector must be larger than `1` for method `{method!r}`, found `{len(use)}`."
            )
        adata_dummy = sc.AnnData(X=X)
        sc.pp.neighbors(adata_dummy, use_rep="X", n_neighbors=n_neighbors_louvain)
        sc.tl.louvain(adata_dummy, resolution=resolution_louvain)
        labels = adata_dummy.obs["louvain"]
    else:
        raise ValueError(
            f"Invalid method `{method!r}`. Valid options are: `'kmeans', 'louvain'`."
        )

    return list(labels)


def _compute_mean_color(color_list: List[str]) -> str:
    """
    Utility function to compute the mean color.
    """

    if not all(map(lambda c: mcolors.is_color_like(c), color_list)):
        raise ValueError("Not all values are valid colors.")

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
    """
    Utility function to check whether the undirected graph encoded by c is connected.
    """

    G = nx.from_scipy_sparse_matrix(c) if issparse(c) else nx.from_numpy_array(c)

    return nx.is_connected(G)


def is_symmetric(c: Union[spmatrix, np.ndarray], ord: str = "fro", eps: float = 1e-4):
    """
    Utility function to check whether the graph encoded by c is symmetric.
    """

    dev = s_norm((c - c.T), ord=ord) if issparse(c) else d_norm((c - c.T), ord=ord)
    return dev < eps


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
    """
    Computes a gaussian kernel.
    """

    if issparse(X):
        G = X.copy()
        G.data = np.exp(-(G.data - mu) ** 2 / (2 * sigma ** 2)) / np.sqrt(
            2 * np.pi * sigma ** 2
        )
    else:
        G = np.exp(-(X - mu) ** 2 / (2 * sigma ** 2)) / np.sqrt(2 * np.pi * sigma ** 2)

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
    if has_neighs(adata):
        C = get_neighs(adata, mode)
        if (
            n_neighbors is not None
            and n_neighbors <= get_neighs_params(adata)["n_neighbors"]
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


def _create_colors(
    base_color: Union[str, Tuple[float, float, float]],
    n: int,
    hue_range: Optional[Tuple[float, float]] = (-0.1, 0.1),
    saturation_range: Optional[Tuple[float, float]] = (-0.3, 0.3),
    value_range: Optional[Tuple[float, float]] = (-0.3, 0.3),
    convert_to_rgb: bool = True,
    as_hex: bool = True,
) -> List[Any]:
    """
    Create variations of colors from base color.

    Params
    ------
    base_color
        Base color which serves as a starting point.
    n
        Number of colors to create.
    hue_range
        Minimum and maximum value to add to the base color's hue.
        If `None`, don't adjust the hue.
    saturation_range
        Minimum and maximum value to add to the base color's saturation.
        If `None`, don't adjust the saturation.
    value_range
        Minimum and maximum value to add to the base color's value.
        If `None`, don't adjust the value.
    convert_to_rgb
        Whether to convert colors from HSV to RGB.
    as_hex:
        Whether to return colors as hex string.

    Returns
    -------
        List of colors, either as a hex string or an array.
    """

    if not mcolors.is_color_like(base_color):
        raise ValueError(f"Base color is not color-like.")
    if n <= 0:
        raise ValueError(f"Number of colors must be > 0, found `{n}`.")

    base_color = mcolors.rgb_to_hsv(mcolors.to_rgb(base_color))

    if n == 1:
        colors = [base_color]
    else:
        n *= 2  # sometimes the colors are too similar, we take every 2nd one
        colors = np.repeat(base_color[..., np.newaxis], n, axis=1).T

        for i, r in enumerate((hue_range, saturation_range, value_range)):
            if r is None:
                continue
            r_low, r_high = sorted(r)
            c = base_color[i]

            colors[:, i] = np.linspace(max(c + r_low, 0), min(c + r_high, 1), n)

    if convert_to_rgb:
        colors = map(mcolors.hsv_to_rgb, colors)
    if as_hex:
        colors = map(mcolors.to_hex, colors)

    return list(colors)[::2]  # we've created twice as much colors, select every other


def _convert_to_hex_colors(colors: Sequence[Any]) -> List[str]:
    if not all(mcolors.is_color_like(c) for c in colors):
        raise ValueError("Not all colors are color-like.")

    return [mcolors.to_hex(c) for c in colors]


def _create_categorical_colors(n_categories: int):
    if n_categories > 51:
        raise ValueError(f"Maximum number of colors (51) exceeded: `{n_categories}`.")
    colors = [cm.Set1(i) for i in range(cm.Set1.N)][:n_categories]
    colors += [cm.Set2(i) for i in range(cm.Set2.N)][: n_categories - len(colors)]
    colors += [cm.Set3(i) for i in range(cm.Set3.N)][: n_categories - len(colors)]
    colors += [cm.tab10(i) for i in range(cm.tab10.N)][: n_categories - len(colors)]
    colors += [cm.Paired(i) for i in range(cm.Paired.N)][: n_categories - len(colors)]

    return _convert_to_hex_colors(colors)


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


def _merge_approx_rcs(
    rc_old: pd.Series, rc_new: pd.Series, inplace: bool = False
) -> Optional[pd.Series]:
    """
    Update approximate recurrent classes with new information. It **can never remove** old categories, only
    add to the existing ones.

    Params
    ------
    rc_old
        Old approximate recurrent classes.
    rc_new
        New approximate recurrent classes.
    inplace
        Whether to update :paramref:`rc_old` or create a copy.

    Returns
    -------
    :class:`pd.Series`
        If paramref:`inplace` is `False`, returns the modified approximate recurrent classes.
    """

    if not is_categorical_dtype(rc_old):
        raise TypeError(
            f"Expected old approx. recurrent classes to be categorical, found "
            f"`{infer_dtype(rc_old)}`."
        )

    if not is_categorical_dtype(rc_new):
        raise TypeError(
            f"Expected new approx. recurrent classes to be categorical, found "
            f"`{infer_dtype(rc_new)}`."
        )

    if (rc_old.index != rc_new.index).any():
        raise ValueError(f"Index for old and new approx. recurrent classes differ.")

    if not inplace:
        rc_old = rc_old.copy()

    mask = ~rc_new.isna()

    if np.sum(mask) == 0:
        return rc_old if not inplace else None

    old_cats = rc_old.cat.categories
    new_cats = rc_new.cat.categories
    cats_to_add = (
        pd.CategoricalIndex(rc_new.loc[mask]).remove_unused_categories().categories
    )

    rc_old.cat.set_categories(old_cats | cats_to_add, inplace=True)
    rc_new.cat.set_categories(old_cats | cats_to_add, inplace=True)

    rc_old.loc[mask] = rc_new.loc[mask]
    rc_old.cat.remove_unused_categories(inplace=True)

    rc_new.cat.set_categories(new_cats, inplace=True)  # return to previous state

    return rc_old if not inplace else None
