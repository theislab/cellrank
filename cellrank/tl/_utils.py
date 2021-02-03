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

from statsmodels.stats.multitest import multipletests

import numpy as np
import pandas as pd
from pandas import Series
from scipy.stats import norm
from numpy.linalg import norm as d_norm
from scipy.sparse import eye as speye
from scipy.sparse import diags, issparse, spmatrix, csr_matrix, isspmatrix_csr
from sklearn.cluster import KMeans
from pandas.api.types import infer_dtype, is_categorical_dtype
from scipy.sparse.linalg import norm as sparse_norm

import matplotlib.colors as mcolors

from cellrank import logging as logg
from cellrank.ul._docs import d
from cellrank.ul._utils import _get_neighs, _has_neighs, _get_neighs_params
from cellrank.tl._colors import (
    _compute_mean_color,
    _convert_to_hex_colors,
    _insert_categorical_colors,
)
from cellrank.tl._constants import ModeEnum
from cellrank.ul._parallelize import parallelize
from cellrank.tl._linear_solver import _solve_lin_system
from cellrank.tl.kernels._utils import np_std, np_mean, _filter_kwargs

AnnData = TypeVar("AnnData")
ColorLike = TypeVar("ColorLike")
GPCCA = TypeVar("GPCCA")
CFLARE = TypeVar("CFLARE")
DiGraph = TypeVar("DiGraph")
KernelDirection = TypeVar("KernelDirection")

EPS = np.finfo(np.float64).eps


class TestMethod(ModeEnum):  # noqa
    FISCHER = "fischer"
    PERM_TEST = "perm_test"


class RandomKeys:
    """
    Create random keys inside an :class:`anndata.AnnData` object.

    Parameters
    ----------
    adata
        Annotated data object.
    n
        Number of keys, If `None`, create just 1 keys.
    where
        Attribute of ``adata``. If `'obs'`, also clean up `'{key}_colors'` for each generated key.
    """

    def __init__(self, adata: AnnData, n: Optional[int] = None, where: str = "obs"):
        self._adata = adata
        self._where = where
        self._n = n
        self._keys = []

    def __enter__(self):
        self._keys = _generate_random_keys(self._adata, self._where, self._n)
        return self._keys

    def __exit__(self, exc_type, exc_val, exc_tb):
        for key in self._keys:
            try:
                del self._adata.obs[key]
            except KeyError:
                pass
            if self._where == "obs":
                try:
                    del self._adata.uns[f"{key}_colors"]
                except KeyError:
                    pass


def _pairwise(iterable: Iterable) -> zip:
    """Return pairs of elements from an iterable."""
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def _min_max_scale(x: np.ndarray) -> np.ndarray:
    """
    Scale a 1D array to 0-1 range.

    Parameters
    ----------
    x
        Array to be scaled.

    Returns
    -------
        The scaled array.
    """
    minn, maxx = np.nanmin(x), np.nanmax(x)
    return (x - minn) / (maxx - minn)


def _process_series(
    series: pd.Series, keys: Optional[List[str]], colors: Optional[np.array] = None
) -> Union[pd.Series, Tuple[pd.Series, List[str]]]:
    """
    Process :class:`pandas.Series` categorical objects.

    Categories in ``series`` are combined/removed according to ``keys``,
    the same transformation is applied to the corresponding colors.

    Parameters
    ----------
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
        Categorical updated annotation. Each cell is assigned to either
        `NaN` or one of updated approximate recurrent classes.
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

    # check that the keys are unique
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

    # remove cats and colors according to keys
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
    Check for imaginary components in columns of X specified by ``use``.

    Parameters
    ----------
    X
        Matrix containing the eigenvectors.
    use
        Selection of columns of `X`.
    use_imag
        For eigenvectors that are complex, use real or imaginary part.

    Returns
    -------
    class:`numpy.ndarray`
        An array containing either only real eigenvectors or also complex ones.
    """

    complex_mask = np.sum(X.imag != 0, axis=0) > 0
    complex_ixs = np.array(use)[np.where(complex_mask)[0]]
    complex_key = "imaginary" if use_imag else "real"
    if len(complex_ixs) > 0:
        logg.warning(
            f"The eigenvectors with indices `{list(complex_ixs)}` have an imaginary part. "
            f"Showing their {complex_key} part"
        )
    X_ = X.real
    if use_imag:
        X_[:, complex_mask] = X.imag[:, complex_mask]

    return X_


def _bias_knn(
    conn: csr_matrix, pseudotime: np.ndarray, n_neighbors: int, k: int = 3
) -> csr_matrix:
    """
    Palantir Kernel utility function.

    This function takes in symmetric connectivities and a pseudotime and removes edges that point "against" pseudotime,
    in this way creating a directed graph. For each node, it always keeps the closest neighbors, making sure the graph
    remains connected.

    Parameters
    ----------
    conn
        The nearest neighbor connectivities.
    pseudotime
        Pseudotemporal ordering of cells.
    n_neighbors
        Number of neighbors to keep
    k
        Number, alongside with ``n_neighbors`` which determined the threshold for candidate indices.

    Returns
    -------
    :class:`scipy.sparse.csr_matrix`.
        Biased connectivities according to the ``pseudotime``.
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


def _mat_mat_corr_sparse(
    X: csr_matrix,
    Y: np.ndarray,
) -> np.ndarray:
    n = X.shape[1]

    X_bar = np.reshape(np.array(X.mean(axis=1)), (-1, 1))
    X_std = np.reshape(
        np.sqrt(np.array(X.power(2).mean(axis=1)) - (X_bar ** 2)), (-1, 1)
    )

    y_bar = np.reshape(np.mean(Y, axis=0), (1, -1))
    y_std = np.reshape(np.std(Y, axis=0), (1, -1))

    with np.warnings.catch_warnings():
        np.warnings.filterwarnings(
            "ignore", r"invalid value encountered in true_divide"
        )
        return (X @ Y - (n * X_bar * y_bar)) / ((n - 1) * X_std * y_std)


def _mat_mat_corr_dense(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    n = X.shape[1]

    X_bar = np.reshape(np_mean(X, axis=1), (-1, 1))
    X_std = np.reshape(np_std(X, axis=1), (-1, 1))

    y_bar = np.reshape(np_mean(Y, axis=0), (1, -1))
    y_std = np.reshape(np_std(Y, axis=0), (1, -1))

    with np.warnings.catch_warnings():
        np.warnings.filterwarnings(
            "ignore", r"invalid value encountered in true_divide"
        )
        return (X @ Y - (n * X_bar * y_bar)) / ((n - 1) * X_std * y_std)


def _perm_test(
    ixs: np.ndarray,
    corr: np.ndarray,
    X: Union[np.ndarray, spmatrix],
    Y: np.ndarray,
    seed: Optional[int] = None,
    queue=None,
) -> Tuple[np.ndarray, np.ndarray]:
    rs = np.random.RandomState(None if seed is None else seed + ixs[0])
    cell_ixs = np.arange(X.shape[1])
    pvals = np.zeros_like(corr, dtype=np.float64)
    corr_bs = np.zeros((len(ixs), X.shape[0], Y.shape[1]))  # perms x genes x lineages

    mmc = _mat_mat_corr_sparse if issparse(X) else _mat_mat_corr_dense

    for i, _ in enumerate(ixs):
        rs.shuffle(cell_ixs)
        corr_i = mmc(X, Y[cell_ixs, :])
        pvals += np.abs(corr_i) >= np.abs(corr)

        bootstrap_ixs = rs.choice(cell_ixs, replace=True, size=len(cell_ixs))
        corr_bs[i, :, :] = mmc(X[:, bootstrap_ixs], Y[bootstrap_ixs, :])

        if queue is not None:
            queue.put(1)

    if queue is not None:
        queue.put(None)

    return pvals, corr_bs


@d.get_sections(base="correlation_test", sections=["Returns"])
def _correlation_test(
    X: Union[np.ndarray, spmatrix],
    Y: "Lineage",  # noqa: F821
    gene_names: Sequence[str],
    method: TestMethod = TestMethod.FISCHER,
    confidence_level: float = 0.95,
    n_perms: Optional[int] = None,
    seed: Optional[int] = None,
    **kwargs,
) -> pd.DataFrame:
    """
    Perform a stastical test.

    Return NaN for genes which don't vary across cells.

    Parameters
    ----------
    X
        Array or sparse matrix of shape ``(n_cells, n_genes)`` containing the expression.
    Y
        Array of shape ``(n_cells, n_lineages)`` containing the absorption probabilities.
    gene_names
        Sequence of shape ``(n_genes,)`` containing the gene names.
    method
        Method for p-value calculation.
    confidence_level
        Confidence level for the confidence interval calculation. Must be in `[0, 1]`.
    n_perms
        Number of permutations if ``method='perm_test'``.
    seed
        Random seed if ``method='perm_test'``.
    kwargs
        Keyword arguments for :func:`cellrank.ul._parallelize.parallelize`.

    Returns
    -------
    :class:`pandas.DataFrame`
        Dataframe of shape ``(n_genes, n_lineages * 5)`` containing the following columns, 1 for each lineage:

            - ``{lineage} corr`` - correlation between the gene expression and absorption probabilities.
            - ``{lineage} pval`` - calulated p-values for double-sided test.
            - ``{lineage} qval`` - corrected p-values using Benjamini-Hochberg method at level `0.05`.
            - ``{lineage} ci low`` - lower bound of the ``confidence_level`` correlation confidence interval.
            - ``{lineage} ci high`` - upper bound of the ``confidence_level`` correlation confidence interval.
    """

    corr, pvals, ci_low, ci_high = _correlation_test_helper(
        X.T,
        Y.X,
        method=method,
        n_perms=n_perms,
        seed=seed,
        confidence_level=confidence_level,
        **kwargs,
    )

    res = pd.DataFrame(corr, index=gene_names, columns=[f"{c} corr" for c in Y.names])
    for idx, c in enumerate(Y.names):
        res[f"{c} pval"] = pvals[:, idx]
        res[f"{c} qval"] = multipletests(pvals[:, idx], alpha=0.05, method="fdr_bh")[1]
        res[f"{c} ci low"] = ci_low[:, idx]
        res[f"{c} ci high"] = ci_high[:, idx]

    res = res[
        [
            f"{c} {stat}"
            for c in Y.names
            for stat in ("corr", "pval", "qval", "ci low", "ci high")
        ]
    ]
    res.sort_values(by=[f"{c} corr" for c in Y.names], ascending=False, inplace=True)

    return res


def _correlation_test_helper(
    X: Union[np.ndarray, spmatrix],
    Y: np.ndarray,
    method: TestMethod = TestMethod.FISCHER,
    n_perms: Optional[int] = None,
    seed: Optional[int] = None,
    confidence_level: float = 0.95,
    **kwargs,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute the correlation between rows in matrix ``X`` columns of matrix ``Y``.

    Parameters
    ----------
    X
        Array or matrix of `(M, N)` elements.
    Y
        Array of `(N, K)` elements.
    method
        Method for p-value calculation.
    n_perms
        Number of permutations if ``method='perm_test'``.
    seed
        Random seed if ``method='perm_test'``.
    confidence_level
        Confidence level for the confidence interval calculation. Must be in `[0, 1]`.
    kwargs
        Keyword arguments for :func:`cellrank.ul._parallelize.parallelize`.

    Returns
    -------
        Correlations, p-values, corrected p-values, lower and upper bound of 95% confidence interval.
        Each array if of shape ``(n_genes, n_lineages)``.
    """

    def perm_test_extractor(
        res: Sequence[Tuple[np.ndarray, np.ndarray]]
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        pvals, corr_bs = zip(*res)
        pvals = np.sum(pvals, axis=0) / float(n_perms)

        corr_bs = np.concatenate(corr_bs, axis=0)
        corr_ci_low, corr_ci_high = np.quantile(corr_bs, q=ql, axis=0), np.quantile(
            corr_bs, q=qh, axis=0
        )

        return pvals, corr_ci_low, corr_ci_high

    if not (0 <= confidence_level <= 1):
        raise ValueError(
            f"Expected `confidence_level` to be in interval `[0, 1]`, found `{confidence_level}`."
        )

    n = X.shape[1]  # genes x cells
    ql = 1 - confidence_level - (1 - confidence_level) / 2.0
    qh = confidence_level + (1 - confidence_level) / 2.0

    if issparse(X) and not isspmatrix_csr(X):
        X = csr_matrix(X)

    corr = _mat_mat_corr_sparse(X, Y) if issparse(X) else _mat_mat_corr_dense(X, Y)

    if method == TestMethod.FISCHER:
        # see: https://en.wikipedia.org/wiki/Pearson_correlation_coefficient#Using_the_Fisher_transformation
        mean, se = np.arctanh(corr), 1.0 / np.sqrt(n - 3)
        z_score = (np.arctanh(corr) - np.arctanh(0)) * np.sqrt(n - 3)

        z = norm.ppf(qh)
        corr_ci_low = np.tanh(mean - z * se)
        corr_ci_high = np.tanh(mean + z * se)
        pvals = 2 * norm.cdf(-np.abs(z_score))

    elif method == TestMethod.PERM_TEST:
        if not isinstance(n_perms, int):
            raise TypeError(
                f"Expected `n_perms` to be an integer, found `{type(n_perms).__name__!r}`."
            )
        if n_perms <= 0:
            raise ValueError(f"Expcted `n_perms` to be positive, found `{n_perms}`.")

        pvals, corr_ci_low, corr_ci_high = parallelize(
            _perm_test,
            np.arange(n_perms),
            as_array=False,
            unit="permutation",
            extractor=perm_test_extractor,
            **kwargs,
        )(corr, X, Y, seed=seed)

    else:
        raise NotImplementedError(method)

    return corr, pvals, corr_ci_low, corr_ci_high


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


def _filter_cells(distances: np.ndarray, rc_labels: Series, n_matches_min: int):
    """Filter out some cells that look like transient states based on their neighbors."""

    if not is_categorical_dtype(rc_labels):
        raise TypeError(
            f"Argument `categories` must be a categorical variable, found `{infer_dtype(rc_labels)}`."
        )

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
        logg.warning(
            "Consider lowering  'n_matches_min' or "
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
    """
    Cluster the rows of the matrix X.

    Parameters
    ----------
    X
        Data matrix of shape `n_samples x n_features`.
    n_clusters
        Number of clusters to use.
    method
        Method to use for clustering. Options are `'kmeans', 'louvain', 'leiden'`.
    n_neighbors
        If using a community-detection based clustering algorithm, number of neighbors for KNN construction.
    resolution
        Resolution parameter for `'louvain', 'leiden'`.

    Returns
    -------
    :class:`list`
        List of cluster labels of length `n_samples`.
    """

    import scanpy as sc

    # make sure data is at least 2D
    if X.ndim == 1:
        X = X.reshape((1, -1))

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
            f"Invalid method `{method!r}`. Valid options are: `'kmeans'`, `'louvain'` or `'leiden'`."
        )

    return list(labels)


def _eigengap(evals: np.ndarray, alpha: float) -> int:
    """
    Compute the eigengap among the top eigenvalues of a matrix.

    Parameters
    ----------
    evals
        Sorted array of real numbers. If complex, take their real part.
    alpha
        Determines how much weight is given to the deviation of an eigenvalue from one.

    Returns
    -------
    int
        Number of eigenvectors to be used.
    """

    if np.iscomplexobj(evals):
        evals = evals.real
    evals = np.sort(evals)[::-1]  # they could be ordered by LM, not LR

    gap, eps = evals[:-1] - evals[1:], (1 - evals)[:-1]
    J = gap - alpha * eps

    return int(np.argmax(J))


def _partition(
    conn: Union[DiGraph, np.ndarray, spmatrix], sort: bool = True
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

    Parameters
    ----------
    conn
        Directed graph to _partition.

    Returns
    -------
    :class:`list`, :class:`list`
        Recurrent and transient classes, respectively.
    """

    import networkx as nx

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


def _connected(c: Union[spmatrix, np.ndarray]) -> bool:
    """Check whether the undirected graph encoded by c is connected."""

    import networkx as nx

    G = nx.from_scipy_sparse_matrix(c) if issparse(c) else nx.from_numpy_array(c)

    return nx.is_connected(G)


def _symmetric(
    matrix: Union[spmatrix, np.ndarray],
    ord: str = "fro",
    eps: float = 1e-4,
    only_check_sparsity_pattern: bool = False,
) -> bool:
    """Check whether the graph encoded by `matrix` is symmetric."""
    if only_check_sparsity_pattern:
        if issparse(matrix):
            return len(((matrix != 0) - (matrix != 0).T).data) == 0
        return ((matrix != 0) == (matrix != 0).T).all()

    if issparse(matrix):
        return sparse_norm((matrix - matrix.T), ord=ord) < eps
    return d_norm((matrix - matrix.T), ord=ord) < eps


def _normalize(X: Union[np.ndarray, spmatrix]) -> Union[np.ndarray, spmatrix]:
    """
    Row-normalizes an array to sum to 1.

    Parameters
    ----------
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
            return X.multiply(csr_matrix(1.0 / np.abs(X).sum(1)))
        X = np.array(X)
        return X / X.sum(1)[:, None]


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

    Parameters
    ----------
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

    Parameters
    ----------
    fig: :class:`matplotlib.figure.Figure`
        Figure to save.
    path:
        Path where to save the figure. If path is relative, save it under ``cellrank.settings.figdir``.
    make_dir:
        Whether to try making the directory if it does not exist.
    ext:
        Extension to use.

    Returns
    -------
    None
        Just saves the plot.
    """

    from cellrank import settings

    if os.path.splitext(path)[1] == "":
        path = f"{path}.{ext}"

    if not os.path.isabs(path):
        path = os.path.join(settings.figdir, path)

    if make_dir:
        _maybe_create_dir(os.path.split(path)[0])

    logg.debug(f"Saving figure to `{path!r}`")

    fig.savefig(path, bbox_inches="tight", transparent=True)


def _convert_to_categorical_series(
    rc_classes: Dict[Union[int, str], Sequence[Union[int, str]]], cell_names: List[str]
) -> Series:
    """
    Convert a mapping of recurrent classes to cells to a :class:`pandas.Series`.

    Parameters
    ----------
    rc_classes
        Recurrent classes in the following format: `{'rc_0': ['cell_0', 'cell_1', ...], ...}`.
    cell_names
        List of valid cell names, usually taken from ``adata.obs_names``.

    Returns
    -------
    :class:`pandas.Series`
        Categorical series where `NaN` mark cells which do not belong to any recurrent class.
    """

    cnames = set(cell_names)
    mapper, expected_size = {}, 0
    for rc, cells in rc_classes.items():
        if not len(cells):
            logg.warning(f"No cells selected for category `{rc!r}`.")
            continue
        cells = [c if isinstance(c, str) else cell_names[c] for c in cells]
        rest = set(cells) - cnames
        if rest:
            raise ValueError(f"Invalid cell names `{list(rest)}`.")
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

    Parameters
    ----------
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
        Whether to update ``old`` or create a copy.

    Returns
    -------
    :class:`pandas.Series`
        If ``inplace`` is `False`, returns the modified approximate recurrent classes and if
        ``colors_old`` and ``colors_new`` are both `None`.
    :class:`numpy.ndarray`
        If ``inplace`` is `True` and any of ``colors_old``, ``colors_new`` containing the new colors.
    :class:`pandas.Series`, :class:`numpy.ndarray`
        The same as above, but with ``inplaces` is `False`.
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


def _info_if_obs_keys_categorical_present(
    adata: AnnData, keys: Iterable[str], msg_fmt: str, warn_once: bool = True
) -> None:
    for key in keys:
        if key in adata.obs.keys() and is_categorical_dtype(adata.obs[key]):
            logg.info(msg_fmt.format(key))
            if warn_once:
                break


def _one_hot(n, cat: Optional[int] = None) -> np.ndarray:
    """
    One-hot encode cat to a vector of length n.

    If cat is `None`, return a vector of zeros.
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

    Given a fuzzy clustering of `n_samples` samples represented by a matrix ``a_fuzzy`` of shape
    `(n_samples x n_clusters)` where rows sum to one and indicate cluster membership to each of
    the `c_clusters` clusters, we compute an assignment of a subset of samples to clusters such
    that each cluster is represented by its ``n_most_likely`` most likely samples. In case a sample
    is assigned more than once, it can either be removed (``remove_overlap=True``) or it can be
    assigned to the cluster it most likely belongs to (``remove_overlap=False``). In case this
    leaves a cluster with less than ``raise_threshold x n_most_likely`` samples, we raise an exception.
    In case this leaves clusters c_1, ..., c_m with less than ``n_most_likely`` samples, but more than
    ``raise_threshold x n_most_likely`` samples, we append c_1, ..., c_m to a list `critical_clusters`,
    which we return.

    We return a boolean matrix `a_discrete` of the same shape as ``a_fuzzy`;`, where `1` in position
    `i, j` indicates that sample `i` is assigned to cluster `j`. Note that we don't assign all samples
    to clusters (most entries in `a_discrete` will be `0`) - this is meant to only assign a small
    subset of the samples, which we are most confident in.

    Parameters
    ----------
    a_fuzzy
        Numpy array of shape `(n_samples x n_clusters)` representing a fuzzy clustering.
        Rows must sum to one.
    n_most_likely
        Number of samples we want to assign to each cluster.
    remove_overlap
        If `True`, remove ambiguous samples. Otherwise, assign them to the most likely cluster.
    raise_threshold
        If a cluster is assigned less than ``raise_threshold x n_most_likely`` samples, raise an
        exception. Set to `None` if you only want to raise if there is an empty cluster.
    check_row_sums
        Check whether rows in `a_fuzzy` sum to one. The one situation where we don't do this is when
        we have selected a couple of main states and we don't want to re-distribute probability mass.

    Returns
    -------
    :class:`numpy.ndarray`m :class:`numpy.ndarray`
        Boolean matrix of the same shape as `a_fuzzy`, assigning a subset of the samples to clusters and
        an rray of clusters with less than `n_most_likely` samples assigned, respectively.
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
            raise ValueError("Rows in `a_fuzzy` do not sum to `1`.")
    if n_most_likely > int(n_samples / n_clusters):
        raise ValueError(
            f"You've selected `{n_most_likely}` cells, please decrease this to at most "
            f"`{int(n_samples / n_clusters)}` cells for your dataset."
        )

    # initialise
    n_raise = (
        1
        if raise_threshold is None
        else np.max([int(raise_threshold * n_most_likely), 1])
    )
    logg.debug(f"Raising an exception if there are less than `{n_raise}` cells.")

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
    membership: np.array,
    index: Optional[Iterable] = None,
    names: Optional[Iterable] = None,
) -> pd.Series:
    """
    Create a pandas Series based on a one-hot encoded matrix.

    Parameters
    ----------
    membership
        One-hot encoded membership matrix, of shape `(n_samples x n_clusters)` i.e. a `1` in position `i, j`
        signifies that sample `i` belongs to cluster `j`.
    index
        Index for the Series. Careful, if this is not given, categories are removed when writing to AnnData.

    Returns
    -------
    :class:`pandas.Series`
        Series, indicating cluster membership for each sample. The data type of the categories is :class:`str`
        and samples that belong to no cluster are assigned `NaN`.
    """
    n_samples, n_clusters = membership.shape
    if not isinstance(membership, np.ndarray):
        raise TypeError(
            f"Expected `membership` to be of type `numpy.ndarray`, found `{type(membership).__name__!r}`."
        )
    membership = np.asarray(
        membership
    )  # change the type in case a lineage object was passed.
    if membership.dtype != np.bool:
        raise TypeError(
            f"Expected `membership`'s elements to be boolean, found `{membership.dtype.name!r}`."
        )

    if not np.all(membership.sum(axis=1) <= 1):
        raise ValueError("Not all items are one-hot encoded or empty.")
    if (membership.sum(0) == 0).any():
        logg.warning(f"Detected {np.sum((membership.sum(0) == 0))} empty categories")

    if index is None:
        index = range(n_samples)
    if names is not None:
        if len(names) != n_clusters:
            raise ValueError(
                f"Shape mismatch, length of `names` is `{len(names)}`, but `n_clusters={n_clusters}`."
            )
    else:
        names = np.arange(n_clusters).astype("str")

    target_series = pd.Series(index=index, dtype="category")
    for vec, name in zip(membership.T, names):
        target_series.cat.add_categories(name, inplace=True)
        target_series[np.where(vec)[0]] = name

    return target_series


def _get_cat_and_null_indices(
    cat_series: Series,
) -> Tuple[np.ndarray, np.ndarray, Dict[Any, np.ndarray]]:
    """
    Given a categorical :class:`pandas.Series`, get the indices corresponding to categories and `NaNs`.

    Parameters
    ----------
    cat_series
        Series that contains categorical annotations.

    Returns
    -------
    :class: `numpy.ndarray`
        Array containing the indices of elements corresponding to categories in ``cat_series``.
    :class: `numpy.ndarray`
        Array containing the indices of elements corresponding to NaNs in ``cat_series``.
    :class:`dict`
        Dict containing categories of ``cat_series`` as keys and an array of corresponding indices as values.
    """

    # check the dtype
    if cat_series.dtype != "category":
        raise TypeError(
            f"Expected `cat_series` to be categorical, found `{cat_series.dtype.name!r}`."
        )

    # define a dict that has category names as keys and arrays of indices as values
    lookup_dict = {
        cat: np.where(cat_series == cat)[0] for cat in cat_series.cat.categories
    }
    all_indices = np.arange(len(cat_series))

    # collect all category indices
    cat_indices = np.concatenate(list(lookup_dict.values()))

    # collect all null indices (the ones where we have NaN in `cat_series`)
    null_indices = np.array(list(set(all_indices) - set(cat_indices)))

    # check that null indices and cat indices are unique
    assert (
        np.unique(cat_indices, return_counts=True)[1] == 1
    ).all(), "Cat indices are not unique."
    assert (
        np.unique(null_indices, return_counts=True)[1] == 1
    ).all(), "Null indices are not unique."

    # check that there is no overlap
    assert (
        len(set(cat_indices).intersection(set(null_indices))) == 0
    ), "Cat and null indices overlap."

    # check that their untion is the set of all indices
    assert set(cat_indices).union(set(null_indices)) == set(
        all_indices
    ), "Some indices got lost on the way."

    return cat_indices, null_indices, lookup_dict


def _check_estimator_type(estimator: Any) -> None:
    # prevents cyclic import
    from cellrank.tl.estimators._base_estimator import BaseEstimator

    if not isinstance(estimator, type):
        raise TypeError(
            f"Expected estimator to be a class, found `{type(estimator).__name__!r}`."
        )

    if not issubclass(estimator, BaseEstimator):
        raise TypeError(
            f"Expected estimator to be a subclass of `cellrank.tl.estimators.BaseEstimator`, "
            f"found `{type(estimator).__name__!r}`."
        )


def _calculate_absorption_time_moments(
    Q: Union[np.ndarray, spmatrix],
    trans_indices: np.ndarray,
    n: int,
    calculate_variance: bool = False,
    **kwargs,
) -> Tuple[np.ndarray, Optional[np.ndarray]]:
    """
    Calculate the mean time until absorption and optionally its variance.

    Parameters
    ----------
    Q
        Transient-transient submatrix of the transition matrix.
    trans_indices
        Transient indices.
    n
        Number of states of the full transition matrix.
    calculate_variance
        Whether to calculate also the variance of time to absorption, not only mean.
    kwargs
        Keyword arguments for :func:`cellrank.tl._lin_solver._solver_lin_system`.

    Returns
    -------
        Mean time until absorption and optionally its variance, based on ``calculate_variance``.
    """

    n_jobs = kwargs.pop("n_jobs", None)
    solve_kwargs = _filter_kwargs(_solve_lin_system, **kwargs)

    logg.debug("Calculating mean time to absorption to any absorbing state")
    m = _solve_lin_system(
        Q,
        np.ones((Q.shape[0],), dtype=np.float32),
        n_jobs=1,
        use_eye=True,
        **solve_kwargs,
    ).squeeze()
    mean = np.zeros(n, dtype=np.float32)
    var = None
    mean[trans_indices] = m

    if calculate_variance:
        logg.debug(
            "Calculating variance of mean time to absorption to any absorbing state"
        )

        I = speye(Q.shape[0]) if issparse(Q) else np.eye(Q.shape[0])  # noqa
        A_t = (I + Q).T
        B_t = (I - Q).T

        logg.debug("Solving equation (1/2)")
        X = _solve_lin_system(A_t, B_t, n_jobs=n_jobs, **kwargs).T
        y = m - X @ (m ** 2)

        logg.debug("Solving equation (2/2)")
        v = _solve_lin_system(X, y, use_eye=False, n_jobs=1, **solve_kwargs).squeeze()
        assert np.all(v >= 0), f"Encountered negative variance: `{v[v < 0]}`."

        var = np.zeros(n, dtype=np.float32)
        var[trans_indices] = v

    return mean, var


def _calculate_lineage_absorption_time_means(
    Q: csr_matrix,
    R: csr_matrix,
    trans_indices: np.ndarray,
    n: int,
    ixs: Dict[str, np.ndarray],
    lineages: Dict[Sequence[str], str],
    **kwargs,
) -> pd.DataFrame:
    """
    Calculate the mean time until absorption and optionally its variance for specific lineages or their combinations.

    Parameters
    ----------
    Q
        Transient-transient submatrix of the transition matrix.
    R
        Transient-recurrent submatrix of the transition matrix.
    trans_indices
        Transient indices.
    n
        Number of states of the full transition matrix.
    ixs
        Mapping of names of absorbing states and their indices in the full transition matrix.
    lineages
        Lineages for which to calculate the mean time until absorption moments.
    kwargs
        Keyword arguments for :func:`cellrank.tl._lin_solver._solver_lin_system`.

    Returns
    -------
    :class:`pandas.DataFrame`
        A :class:`pandas.DataFrame. with means and optionally variances of
        mean time to absorption for each lineage in ``lineages``.

        Uses more efficient implementation if compute the time for all lineages.
    """

    res = pd.DataFrame()
    if len(lineages) == 1 and set(next(iter(lineages.keys()))) == set(ixs.keys()):
        # use faster implementation in this case
        name = ", ".join(ixs.keys())
        res[f"{name} mean"], var = _calculate_absorption_time_moments(
            Q,
            trans_indices,
            n,
            calculate_variance=next(iter(lineages.values())) == "var",
            **kwargs,
        )
        if var is not None:
            res[f"{name} var"] = var

        return res

    res = pd.DataFrame()
    tmp_ixs, cnt = {}, 0
    for k, ix in ixs.items():
        # get the indices to B matrix
        tmp_ixs[k] = np.arange(cnt, cnt + len(ix), dtype=np.int32)
        cnt += len(ix)

    I = speye(Q.shape[0]) if issparse(Q) else np.eye(Q.shape)  # noqa
    N_inv = I - Q

    logg.debug("Solving equation for `B`")
    B = _solve_lin_system(Q, R, use_eye=True, **kwargs)

    no_jobs_kwargs = kwargs.copy()
    _ = no_jobs_kwargs.pop("n_jobs", None)

    for lns, moment in lineages.items():
        name = ", ".join(lns)
        ix = np.concatenate([ixs[ln] for ln in lns])

        D_j = diags(np.sum(B[:, np.concatenate([tmp_ixs[ln] for ln in lns])], axis=1))
        D_j_inv = D_j.copy()
        D_j_inv.data = 1.0 / D_j.data

        logg.debug(f"Calculating mean time to absorption to `{name!r}`")
        m = _solve_lin_system(
            D_j_inv @ N_inv @ D_j, np.ones(Q.shape[0]), **kwargs
        ).squeeze()

        mean = np.empty(n, dtype=np.float64)
        mean[:] = np.inf
        mean[ix] = 0
        mean[trans_indices] = m

        res[f"{name} mean"] = mean

        if moment == "var":
            logg.debug(f"Calculating variance of mean time to absorption to `{name!r}`")

            logg.debug("Solving equation (1/2)")
            X = _solve_lin_system(D_j + Q @ D_j, N_inv @ D_j, use_eye=False, **kwargs)
            y = m - X @ (m ** 2)

            logg.debug("Solving equation (2/2)")
            v = _solve_lin_system(
                X, y, use_eye=False, n_jobs=1, **no_jobs_kwargs
            ).squeeze()
            assert np.all(v >= 0), f"Encountered negative variance: `{v[v < 0]}`."

            var = np.empty(n, dtype=np.float64)
            var[:] = np.inf
            var[ix] = 0
            var[trans_indices] = v

            res[f"{name} var"] = var

    return res


def _create_initial_terminal_annotations(
    adata: AnnData,
    terminal_key: str = "terminal_states",
    initial_key: str = "initial_states",
    terminal_prefix: Optional[str] = "terminal",
    initial_prefix: Optional[str] = "initial",
    key_added: Optional[str] = "initial_terminal",
) -> None:
    """
    Create categorical annotations of both initial and terminal states.

    This is a utility function for creating a categorical :class:`pandas.Series` object which combines
    the information about initial and terminal states. The :class:`pandas.Series` is written directly
    to the :class:`anndata.AnnData`object. This can for example be used to create a scatter plot in :mod:`scvelo`.

    Parameters
    ----------
    adata
        AnnData object to write to ``.obs[key_added]``.
    terminal_key
        Key from ``adata.obs`` where terminal states have been saved.
    initial_key
        Key from ``adata.obs`` where initial states have been saved.
    terminal_prefix
        Forward direction prefix used in the annotations.
    initial_prefix
        Backward direction prefix used in the annotations.
    key_added
        Key added to ``adata.obs``.

    Returns
    -------
    None
        Nothing, just writes to ``adata``.
    """

    # get both Series objects
    cats_final, colors_final = (
        adata.obs[terminal_key],
        adata.uns[f"{terminal_key}_colors"],
    )
    cats_root, colors_root = adata.obs[initial_key], adata.uns[f"{initial_key}_colors"]

    # merge
    cats_merged, colors_merged = _merge_categorical_series(
        cats_final, cats_root, list(colors_final), list(colors_root)
    )

    # adjust the names
    final_names = cats_final.cat.categories
    final_labels = [
        f"{terminal_prefix if key in final_names else initial_prefix}: {key}"
        for key in cats_merged.cat.categories
    ]
    cats_merged.cat.rename_categories(final_labels, inplace=True)

    # write to AnnData
    adata.obs[key_added] = cats_merged
    adata.uns[f"{key_added}_colors"] = colors_merged
