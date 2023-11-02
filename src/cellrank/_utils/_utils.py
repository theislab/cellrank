import contextlib
import functools
import inspect
import itertools
import os
import types
import warnings
from typing import (
    Any,
    Callable,
    Dict,
    Hashable,
    Iterable,
    List,
    Literal,
    Optional,
    Sequence,
    Tuple,
    TypeVar,
    Union,
)

import wrapt

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.stats as st
from pandas.api.types import infer_dtype
from sklearn.cluster import KMeans
from statsmodels.stats.multitest import multipletests

from matplotlib import colors

import scanpy as sc
from anndata import AnnData
from anndata.utils import make_index_unique

from cellrank import logging as logg
from cellrank._utils._colors import (
    _compute_mean_color,
    _convert_to_hex_colors,
    _insert_categorical_colors,
)
from cellrank._utils._docs import d
from cellrank._utils._enum import ModeEnum
from cellrank._utils._linear_solver import _solve_lin_system

ColorLike = TypeVar("ColorLike")
GPCCA = TypeVar("GPCCA")
CFLARE = TypeVar("CFLARE")
DiGraph = TypeVar("DiGraph")

EPS = np.finfo(np.float64).eps


class TestMethod(ModeEnum):
    FISHER = "fisher"
    PERM_TEST = "perm_test"


class RandomKeys:
    """Create random keys inside an :class:`~anndata.AnnData` object.

    Parameters
    ----------
    adata
        Annotated data object.
    n
        Number of keys, If `None`, create just 1 keys.
    where
        Attribute of ``adata``. If `'obs'`, also clean up `'{key}_colors'` for each generated key.
    seed
        Random seed.
    """

    def __init__(self, adata: AnnData, n: Optional[int] = None, where: str = "obs", seed: int = 0):
        self._adata = adata
        self._where = where
        self._n = n or 1
        self._seed = seed
        self._keys = []

    def _generate_random_keys(self):
        rng = np.random.default_rng(self._seed)
        where = getattr(self._adata, self._where)
        names, seen = [], set(where.keys())

        while len(names) != self._n:
            name = f"RNG_COL_{rng.integers(0, 2 ** 32 - 1)}"
            if name not in seen:
                seen.add(name)
                names.append(name)

        return names

    def __enter__(self):
        self._keys = self._generate_random_keys()
        return self._keys

    def __exit__(self, exc_type, exc_val, exc_tb):
        for key in self._keys:
            with contextlib.suppress(KeyError):
                getattr(self._adata, self._where).drop(key, axis="columns", inplace=True)

            if self._where == "obs":
                with contextlib.suppress(KeyError):
                    del self._adata.uns[f"{key}_colors"]


def _filter_kwargs(_fn: Callable, **kwargs: Any) -> dict:
    """Filter keyword arguments.

    Parameters
    ----------
    _fn
        Function for which to filter keyword arguments.
    kwargs
        Keyword arguments to filter

    Returns
    -------
    dict
        Filtered keyword arguments for the given function.
    """
    sig = inspect.signature(_fn).parameters
    return {k: v for k, v in kwargs.items() if k in sig}


def _pairwise(iterable: Iterable) -> zip:
    """Return pairs of elements from an iterable."""
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)


def _min_max_scale(x: np.ndarray) -> np.ndarray:
    """Scale a 1D array to 0-1 range.

    Parameters
    ----------
    x
        Array to be scaled.

    Returns
    -------
    The scaled array.
    """
    minn, maxx = np.nanmin(x), np.nanmax(x)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return (x - minn) / (maxx - minn)


def _process_series(
    series: pd.Series, keys: Optional[List[str]], cols: Optional[np.array] = None
) -> Union[pd.Series, Tuple[pd.Series, List[str]]]:
    """Process :class:`~pandas.Series` of categorical objects.

    Categories in ``series`` are combined/removed according to ``keys``,
    the same transformation is applied to the corresponding colors.

    Parameters
    ----------
    series
        Input categorical data.
    keys
        Keys could be e.g. `['cat_1, cat_2', 'cat_4']`. If originally,
        there were 4 categories in `series`, then this would combine the first
        and the second and remove the third. The same would be done to `colors`,
        i.e. the first and second color would be merged (average color), while
        the third would be removed.
    cols
        List of colors which aligns with the order of the categories.

    Returns
    -------
    - Categorical updated annotation. Each cell is assigned to either
      `NaN` or one of updated approximate recurrent classes.
    - Color list processed according to keys.
    """
    # determine whether we want to process colors as well
    process_colors = cols is not None

    # assert dtype of the series
    if not isinstance(series.dtype, pd.CategoricalDtype):
        raise TypeError(f"Series must be `categorical`, found `{infer_dtype(series)}`.")

    # if keys is None, just return
    if keys is None:
        if process_colors:
            return series, cols
        return series

    # initialize a copy of the series object
    series_in = series.copy()
    if process_colors:
        colors_in = np.array(cols.copy())
        if len(colors_in) != len(series_in.cat.categories):
            raise ValueError(
                f"Length of colors ({len(colors_in)}) does not match length of "
                f"categories ({len(series_in.cat.categories)})."
            )
        if not all(colors.is_color_like(c) for c in colors_in):
            raise ValueError("Not all colors are color-like.")

    # define a set of keys
    keys_ = {tuple(sorted({key.strip(" ") for key in rc.strip(" ,").split(",")})) for rc in keys}

    # check that the keys are unique
    overlap = [set(ks) for ks in keys_]
    for c1, c2 in itertools.combinations(overlap, 2):
        overlap = c1 & c2
        if overlap:
            raise ValueError(f"Found overlapping keys: `{list(overlap)}`.")

    # check the `keys` are all proper categories
    remaining_cat = [b for a in keys_ for b in a]
    if not np.all(np.in1d(remaining_cat, series_in.cat.categories)):
        raise ValueError("Not all keys are proper categories. Check for spelling mistakes in `keys`.")

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
            new_cat_name = ", ".join(cat)
            mask = np.repeat(False, len(series_temp))
            for key in cat:
                mask = np.logical_or(mask, series_temp == key)
                remaining_cat.remove(key)
            series_temp = series_temp.cat.add_categories(new_cat_name)
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
    series_temp = series_temp.cat.remove_unused_categories()
    series_temp = series_temp.cat.reorder_categories(remaining_cat)

    if process_colors:
        # original colors can still be present, convert to hex
        colors_temp = _convert_to_hex_colors([colors_mod[c] for c in series_temp.cat.categories])
        return series_temp, colors_temp

    return series_temp


def _complex_warning(X: np.array, use: Union[list, int, tuple, range], use_imag: bool = False) -> np.ndarray:
    """Check for imaginary components in columns of X specified by ``use``.

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


def _mat_mat_corr_sparse(
    X: sp.csr_matrix,
    Y: np.ndarray,
) -> np.ndarray:
    n = X.shape[1]

    X_bar = np.reshape(np.array(X.mean(axis=1)), (-1, 1))
    X_std = np.reshape(np.sqrt(np.array(X.power(2).mean(axis=1)) - (X_bar**2)), (-1, 1))

    y_bar = np.reshape(np.mean(Y, axis=0), (1, -1))
    y_std = np.reshape(np.std(Y, axis=0), (1, -1))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return (X @ Y - (n * X_bar * y_bar)) / ((n - 1) * X_std * y_std)


def _mat_mat_corr_dense(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    from cellrank.kernels._utils import np_mean, np_std

    n = X.shape[1]

    X_bar = np.reshape(np_mean(X, axis=1), (-1, 1))
    X_std = np.reshape(np_std(X, axis=1), (-1, 1))

    y_bar = np.reshape(np_mean(Y, axis=0), (1, -1))
    y_std = np.reshape(np_std(Y, axis=0), (1, -1))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return (X @ Y - (n * X_bar * y_bar)) / ((n - 1) * X_std * y_std)


def _perm_test(
    ixs: np.ndarray,
    corr: np.ndarray,
    X: Union[np.ndarray, sp.spmatrix],
    Y: np.ndarray,
    seed: Optional[int] = None,
    queue=None,
) -> Tuple[np.ndarray, np.ndarray]:
    rs = np.random.RandomState(None if seed is None else seed + ixs[0])
    cell_ixs = np.arange(X.shape[1])
    pvals = np.zeros_like(corr, dtype=np.float64)
    corr_bs = np.zeros((len(ixs), X.shape[0], Y.shape[1]))  # perms x genes x lineages

    mmc = _mat_mat_corr_sparse if sp.issparse(X) else _mat_mat_corr_dense

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
@d.dedent
def _correlation_test(
    X: Union[np.ndarray, sp.spmatrix],
    Y: "Lineage",  # noqa: F821
    gene_names: Sequence[str],
    method: TestMethod = TestMethod.FISHER,
    confidence_level: float = 0.95,
    n_perms: Optional[int] = None,
    seed: Optional[int] = None,
    **kwargs: Any,
) -> pd.DataFrame:
    """Perform a statistical test.

    Return NaN for genes which don't vary across cells.

    Parameters
    ----------
    X
        Array or sparse matrix of shape ``(n_cells, n_genes)`` containing the expression.
    Y
        Array of shape ``(n_cells, n_lineages)`` containing the fate probabilities.
    gene_names
        Sequence of shape ``(n_genes,)`` containing the gene names.
    method
        Method for p-value calculation.
    confidence_level
        Confidence level for the confidence interval calculation. Must be in :math:`[0, 1]`.
    n_perms
        Number of permutations if ``method = 'perm_test'``.
    seed
        Random seed if ``method = 'perm_test'``.
    %(parallel)s

    Returns
    -------
    Dataframe of shape ``(n_genes, n_lineages * 5)`` containing the following columns, one for each lineage:

    - ``'{lineage}_corr'`` - correlation between the gene expression and fate probabilities.
    - ``'{lineage}_pval'`` - calculated p-values for double-sided test.
    - ``'{lineage}_qval'`` - corrected p-values using Benjamini-Hochberg method at level `0.05`.
    - ``'{lineage}_ci_low'`` - lower bound of the ``confidence_level`` correlation confidence interval.
    - ``'{lineage}_ci_high'`` - upper bound of the ``confidence_level`` correlation confidence interval.
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
    invalid = (corr < -1) | (corr > 1)
    if np.any(invalid):
        logg.warning(
            f"Found `{np.sum(invalid)}` correlation(s) that are not in `[0, 1]`. "
            f"This usually happens when gene expression is constant across all cells. "
            f"Setting to `NaN`"
        )
        corr[invalid] = np.nan
        pvals[invalid] = np.nan
        ci_low[invalid] = np.nan
        ci_high[invalid] = np.nan

    res = pd.DataFrame(corr, index=gene_names, columns=[f"{c}_corr" for c in Y.names])
    for idx, c in enumerate(Y.names):
        p = pvals[:, idx]
        valid_mask = ~np.isnan(p)

        res[f"{c}_pval"] = p
        res[f"{c}_qval"] = np.nan
        if np.any(valid_mask):
            res.loc[gene_names[valid_mask], f"{c}_qval"] = multipletests(p[valid_mask], alpha=0.05, method="fdr_bh")[1]
        res[f"{c}_ci_low"] = ci_low[:, idx]
        res[f"{c}_ci_high"] = ci_high[:, idx]

    # fmt: off
    res = res[[f"{c}_{stat}" for c in Y.names for stat in ("corr", "pval", "qval", "ci_low", "ci_high")]]
    return res.sort_values(by=[f"{c}_corr" for c in Y.names], ascending=False)
    # fmt: on


def _correlation_test_helper(
    X: Union[np.ndarray, sp.spmatrix],
    Y: np.ndarray,
    method: TestMethod = TestMethod.FISHER,
    n_perms: Optional[int] = None,
    seed: Optional[int] = None,
    confidence_level: float = 0.95,
    **kwargs: Any,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute the correlation between rows in matrix ``X`` columns of matrix ``Y``.

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
        Random seed if ``method = 'perm_test'``.
    confidence_level
        Confidence level for the confidence interval calculation. Must be in `[0, 1]`.
    kwargs
        Keyword arguments for :func:`cellrank._utils._parallelize.parallelize`.

    Returns
    -------
        Correlations, p-values, corrected p-values, lower and upper bound of 95% confidence interval.
        Each array if of shape ``(n_genes, n_lineages)``.
    """
    from cellrank._utils._parallelize import parallelize

    def perm_test_extractor(res: Sequence[Tuple[np.ndarray, np.ndarray]]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        pvals, corr_bs = zip(*res)
        pvals = np.sum(pvals, axis=0) / float(n_perms)

        corr_bs = np.concatenate(corr_bs, axis=0)
        corr_ci_low, corr_ci_high = np.quantile(corr_bs, q=ql, axis=0), np.quantile(corr_bs, q=qh, axis=0)

        return pvals, corr_ci_low, corr_ci_high

    if not (0 <= confidence_level <= 1):
        raise ValueError(f"Expected `confidence_level` to be in interval `[0, 1]`, found `{confidence_level}`.")

    n = X.shape[1]  # genes x cells
    ql = 1 - confidence_level - (1 - confidence_level) / 2.0
    qh = confidence_level + (1 - confidence_level) / 2.0

    if sp.issparse(X) and not sp.isspmatrix_csr(X):
        X = sp.csr_matrix(X)

    corr = _mat_mat_corr_sparse(X, Y) if sp.issparse(X) else _mat_mat_corr_dense(X, Y)

    if method == TestMethod.FISHER:
        # see: https://en.wikipedia.org/wiki/Pearson_correlation_coefficient#Using_the_Fisher_transformation
        mean, se = np.arctanh(corr), 1.0 / np.sqrt(n - 3)
        z_score = (np.arctanh(corr) - np.arctanh(0)) * np.sqrt(n - 3)

        z = st.norm.ppf(qh)
        corr_ci_low = np.tanh(mean - z * se)
        corr_ci_high = np.tanh(mean + z * se)
        pvals = 2 * st.norm.cdf(-np.abs(z_score))

    elif method == TestMethod.PERM_TEST:
        if not isinstance(n_perms, int):
            raise TypeError(f"Expected `n_perms` to be an integer, found `{type(n_perms).__name__}`.")
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


def _filter_cells(distances: sp.spmatrix, rc_labels: pd.Series, n_matches_min: int) -> pd.Series:
    """Filter out some cells that look like transient states based on their neighbors."""
    if not isinstance(rc_labels.dtype, pd.CategoricalDtype):
        raise TypeError(f"Expected `categories` be `categorical`, found `{infer_dtype(rc_labels)}`.")

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

    if np.any((freqs_new / freqs_orig) < 0.5):
        logg.warning(
            "Consider lowering  'n_matches_min' or "
            "increasing 'n_neighbors_filtering'. This filters out too many cells"
        )

    return rc_labels


def _cluster_X(
    X: Union[np.ndarray, sp.spmatrix],
    n_clusters: int,
    method: Literal["leiden", "kmeans"] = "leiden",
    n_neighbors: int = 20,
    resolution: float = 1.0,
) -> List[Any]:
    """Cluster the rows of ``X``.

    Parameters
    ----------
    X
        Matrix of shape ``n_samples x n_features``.
    n_clusters
        Number of clusters to use.
    method
        Method to use for clustering. Options are `'kmeans'`, `'leiden'`.
    n_neighbors
        If using a community-detection based clustering algorithm, number of neighbors for KNN construction.
    resolution
        Resolution parameter for `'leiden'` clustering.

    Returns
    -------
    List of cluster labels of length `n_samples`.
    """
    if X.shape[0] == 1:
        # leiden issue
        return [0]

    if method == "kmeans":
        kmeans = KMeans(n_clusters=n_clusters).fit(X)
        labels = kmeans.labels_
    elif method == "leiden":
        adata_dummy = sc.AnnData(X=X)
        sc.pp.neighbors(adata_dummy, use_rep="X", n_neighbors=n_neighbors)
        sc.tl.leiden(adata_dummy, resolution=resolution)
        labels = adata_dummy.obs[method]
    else:
        raise NotImplementedError(f"Invalid method `{method}`. Valid options are `kmeans` or `leiden`.")

    return list(labels)


def _eigengap(evals: np.ndarray, alpha: float) -> int:
    """Compute the eigengap among the top eigenvalues of a matrix.

    Parameters
    ----------
    evals
        Sorted array of real numbers. If complex, take their real part.
    alpha
        Determines how much weight is given to the deviation of an eigenvalue from one.

    Returns
    -------
    Number of eigenvectors to be used.
    """
    if np.iscomplexobj(evals):
        evals = evals.real
    evals = np.sort(evals)[::-1]  # they could be ordered by LM, not LR

    gap, eps = evals[:-1] - evals[1:], (1 - evals)[:-1]
    J = gap - alpha * eps

    return int(np.argmax(J))


def _partition(
    conn: Union[DiGraph, np.ndarray, sp.spmatrix], sort: bool = True
) -> Tuple[List[List[Any]], List[List[Any]]]:
    """Partition a directed graph into its transient and recurrent classes.

    In a directed graph *G*, node *j* is accessible from node *i* if there exists a path from *i* to *j*.
    If *i* is accessible from *j* and the converse holds as well, then *i* and *j* communicate.
    Communication forms and equivalence relation on directed graphs, so every directed graph can be uniquely partitioned
    into its communication classes (also called strongly connected components).

    If *G* describes the state space of a Markov chain, then communication classes are often
    characterized as either recurrent or transient. Intuitively, once the process enters a recurrent class, it will
    never leave it again.

    Parameters
    ----------
    conn
        Directed graph to partition.

    Returns
    -------
    Recurrent and transient classes, respectively.
    """
    import networkx as nx

    start = logg.debug("Partitioning the graph into current and transient classes")

    def partition(g):
        yield from (
            (
                (sorted(scc) if sort else scc),
                all((not nx.has_path(g, s, t) for s, t in itertools.product(scc, g.nodes - scc))),
            )
            for scc in nx.strongly_connected_components(g)
        )

    def maybe_sort(iterable):
        return sorted(iterable, key=lambda x: (-len(x), x[0])) if sort else list(map(list, iterable))

    rec_classes, trans_classes = itertools.tee(
        partition(nx.DiGraph(conn) if not isinstance(conn, nx.DiGraph) else conn), 2
    )

    rec_classes = (node for node, is_rec in rec_classes if is_rec)
    trans_classes = (node for node, is_rec in trans_classes if not is_rec)

    logg.debug("    Finish", time=start)

    return maybe_sort(rec_classes), maybe_sort(trans_classes)


def _connected(c: Union[sp.spmatrix, np.ndarray]) -> bool:
    """Check whether the undirected graph is connected."""
    import networkx as nx

    if sp.issparse(c):
        try:
            G = nx.from_scipy_sparse_array(c)
        except AttributeError:
            # bwd compatibility for `networkx <2.7`
            G = nx.from_scipy_sparse_matrix(c)
    else:
        G = nx.from_numpy_array(c)

    return nx.is_connected(G)


def _irreducible(d: Union[sp.spmatrix, np.ndarray]) -> bool:
    """Check whether the undirected graph encoded by d is irreducible."""
    import networkx as nx

    G = nx.DiGraph(d) if not isinstance(d, nx.DiGraph) else d
    try:
        it = iter(nx.strongly_connected_components(G))
        _ = next(it)
        _ = next(it)
        return False
    except StopIteration:
        return True


def _symmetric(
    matrix: Union[sp.spmatrix, np.ndarray],
    ord: str = "fro",
    eps: float = 1e-4,
    only_check_sparsity_pattern: bool = False,
) -> bool:
    """Check whether the graph encoded by `matrix` is symmetric."""
    if only_check_sparsity_pattern:
        if sp.issparse(matrix):
            return len(((matrix != 0) - (matrix != 0).T).data) == 0
        return ((matrix != 0) == (matrix != 0).T).all()

    if sp.issparse(matrix):
        return sp.linalg.norm((matrix - matrix.T), ord=ord) < eps
    return np.linalg.norm((matrix - matrix.T), ord=ord) < eps


def _normalize(
    X: Union[np.ndarray, sp.spmatrix],
) -> Union[np.ndarray, sp.spmatrix]:
    """Row-normalizes an array to sum to :math:`1`.

    Parameters
    ----------
    X
        Array to be normalized.

    Returns
    -------
    :class:`numpy.ndarray` or :class:`scipy.sparse.spmatrix`
        The normalized array.
    """
    with np.errstate(divide="ignore"):
        if sp.issparse(X):
            return X.multiply(sp.csr_matrix(1.0 / np.abs(X).sum(1)))
        X = np.array(X)
        return X / (X.sum(1)[:, None])


def _get_connectivities(
    adata: AnnData, mode: str = "connectivities", n_neighbors: Optional[int] = None
) -> Optional[sp.spmatrix]:
    # utility function, copied from scvelo
    if _has_neighs(adata):
        C = _get_neighs(adata, mode)
        if n_neighbors is not None and n_neighbors <= _get_neighs_params(adata)["n_neighbors"]:
            C = (
                _select_connectivities(C, n_neighbors)
                if mode == "connectivities"
                else _select_distances(C, n_neighbors)
            )

        return C.tocsr().astype(np.float32)


def _select_connectivities(connectivities: sp.spmatrix, n_neighbors: Optional[int] = None) -> sp.spmatrix:
    # utility function, copied from scvelo
    C = connectivities.copy()
    n_counts = (C > 0).sum(1).A1 if sp.issparse(C) else (C > 0).sum(1)
    n_neighbors = n_counts.min() if n_neighbors is None else min(n_counts.min(), n_neighbors)
    rows = np.where(n_counts > n_neighbors)[0]
    cumsum_neighs = np.insert(n_counts.cumsum(), 0, 0)
    dat = C.data

    for row in rows:
        n0, n1 = cumsum_neighs[row], cumsum_neighs[row + 1]
        rm_idx = n0 + dat[n0:n1].argsort()[::-1][n_neighbors:]
        dat[rm_idx] = 0
    C.eliminate_zeros()

    return C


def _select_distances(dist, n_neighbors: Optional[int] = None) -> sp.spmatrix:
    # utility function, copied from scvelo
    D = dist.copy()
    n_counts = (D > 0).sum(1).A1 if sp.issparse(D) else (D > 0).sum(1)
    n_neighbors = n_counts.min() if n_neighbors is None else min(n_counts.min(), n_neighbors)
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
    """Try creating a directory if it does not already exist.

    Parameters
    ----------
    dirpath
        Path of the directory to create.

    Returns
    -------
    Nothing, just creates a directory if it doesn't exist.
    """
    if not os.path.exists(dirpath) or not os.path.isdir(dirpath):
        with contextlib.suppress(OSError):
            os.makedirs(dirpath, exist_ok=True)


def save_fig(fig, path: Union[str, os.PathLike], make_dir: bool = True, ext: str = "png") -> None:
    """Save a plot.

    Parameters
    ----------
    fig
        Figure to save.
    path
        Path where to save the figure. If path is relative, save it under ``cellrank.settings.figdir``.
    make_dir
        Whether to try making the directory if it does not exist.
    ext
        Extension to use.

    Returns
    -------
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
    term_states: Dict[Union[int, str], Sequence[Union[int, str]]], cell_names: List[str]
) -> pd.Series:
    """Convert a mapping of terminal states to cells to a :class:`~pandas.Series`.

    Parameters
    ----------
    term_states
        Terminal states in the following format: `{'state_0': ['cell_0', 'cell_1', ...], ...}`.
    cell_names
        List of valid cell names, usually taken from :attr:`~anndata.AnnData.obs_names`.

    Returns
    -------
    Categorical series where `NaN` mark cells which do not belong to any recurrent class.
    """
    cnames = set(cell_names)
    mapper, expected_size = {}, 0
    for ts, cells in term_states.items():
        if not len(cells):
            logg.warning(f"No cells selected for category `{ts!r}`. Skipping")
            continue
        cells = [c if isinstance(c, str) else cell_names[c] for c in cells]
        rest = set(cells) - cnames
        if rest:
            raise ValueError(f"Invalid cell names `{list(rest)}`.")
        mapper[str(ts)] = cells
        expected_size += 1

    if len(mapper) != expected_size:
        raise ValueError(
            "All terminal states are being converted to strings, ensure "
            "that there are no conflicting keys, such as `0` and `'0'`."
        )

    term_states = pd.Series([None] * len(cell_names), index=cell_names, dtype=str)
    for ts, cells in mapper.items():
        term_states[cells] = ts

    term_states = term_states.astype("category")
    if not len(term_states.cat.categories):
        raise ValueError("No categories have been selected.")

    return term_states


def _merge_categorical_series(
    old: pd.Series,
    new: pd.Series,
    colors_old: Union[List[ColorLike], np.ndarray, Dict[Any, ColorLike]] = None,
    colors_new: Union[List[ColorLike], np.ndarray, Dict[Any, ColorLike]] = None,
    color_overwrite: bool = False,
) -> Optional[Union[pd.Series, Tuple[pd.Series, np.ndarray]]]:
    """Update categorical :class:`~pandas.Series` with new information.

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

    Returns
    -------
    - Returns the modified approximate recurrent classes and if ``colors_old`` and ``colors_new`` are both `None`.
    - If any of ``colors_old``, ``colors_new`` contain the new colors.
    """

    def get_color_mapper(
        series: pd.Series,
        cols: Union[List[ColorLike], np.ndarray, Dict[Any, ColorLike]],
    ):
        if len(series.cat.categories) != len(cols):
            raise ValueError(f"Series ({len(series.cat.categories)}) and colors ({len(colors_new)}) differ in length.")

        if isinstance(cols, dict):
            if set(series.cat.categories) != set(cols.keys()):
                raise ValueError("Color mapper and series' categories don't share the keys.")
        else:
            cols = dict(zip(series.cat.categories, cols))

        for color in cols.values():
            if not colors.is_color_like(color):
                raise ValueError(f"Color `{color}` is not color-like.")

        return cols

    if not isinstance(old.dtype, pd.CategoricalDtype):
        raise TypeError(f"Expected old approx. recurrent classes to be categorical, found " f"`{infer_dtype(old)}`.")

    if not isinstance(new.dtype, pd.CategoricalDtype):
        raise TypeError(f"Expected new approx. recurrent classes to be categorical, found " f"`{infer_dtype(new)}`.")

    if (old.index != new.index).any():
        raise ValueError("Index for old and new approx. recurrent classes differ.")

    old, new = old.copy(), new.copy()
    mask = ~new.isna()
    if np.sum(mask) == 0:
        return old

    old_cats = old.cat.categories
    new_cats = new.cat.categories
    cats_to_add = pd.CategoricalIndex(new.loc[mask]).remove_unused_categories().categories

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

    tmp = pd.CategoricalIndex(old_cats).union(pd.CategoricalIndex(cats_to_add))
    old = old.cat.set_categories(tmp)
    new = new.cat.set_categories(tmp)

    old.loc[mask] = new.loc[mask]
    old = old.cat.remove_unused_categories()

    if not colors_old and not colors_new:
        return old

    colors_merged = {**colors_old, **colors_new} if color_overwrite else {**colors_new, **colors_old}
    colors_merged = np.array([colors_merged[c] for c in old.cat.categories])

    return old, colors_merged


def _unique_order_preserving(iterable: Iterable[Hashable]) -> List[Hashable]:
    """Remove items from an iterable while preserving the order."""
    seen = set()
    return [i for i in iterable if i not in seen and not seen.add(i)]


def _one_hot(n, cat: Optional[int] = None) -> np.ndarray:
    """
    One-hot encode cat to a vector of length n.

    If cat is `None`, return a vector of zeros.
    """
    out = np.zeros(n, dtype=bool)
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
    """Map fuzzy clustering to discrete clustering.

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
    Boolean matrix of the same shape as `a_fuzzy`, assigning a subset of the samples to clusters and
    an array of clusters with less than `n_most_likely` samples assigned, respectively.
    """
    # check the inputs
    n_samples, n_clusters = a_fuzzy.shape
    if not isinstance(a_fuzzy, np.ndarray):
        raise TypeError(f"Expected `a_fuzzy` to be of type `numpy.ndarray`, got `{type(a_fuzzy).__name__}`.")
    a_fuzzy = np.asarray(a_fuzzy)  # convert to array from lineage classs, don't copy
    if check_row_sums and n_clusters != 1 and not np.allclose(a_fuzzy.sum(1), 1, rtol=1e6 * EPS, atol=1e6 * EPS):
        raise ValueError("Rows in `a_fuzzy` do not sum to `1`.")
    if n_most_likely > int(n_samples / n_clusters):
        raise ValueError(
            f"You've selected `{n_most_likely}` cells, please decrease this to at most "
            f"`{int(n_samples / n_clusters)}` cells for your dataset."
        )

    # initialise
    n_raise = 1 if raise_threshold is None else np.max([int(raise_threshold * n_most_likely), 1])
    logg.debug(f"Raising an exception if there are less than `{n_raise}` cells.")

    # initially select `n_most_likely` samples per cluster
    sample_assignment = {
        cl: fuzzy_assignment.argpartition(-n_most_likely)[-n_most_likely:]
        for cl, fuzzy_assignment in enumerate(a_fuzzy.T)
    }

    # create the one-hot encoded discrete clustering
    a_discrete = np.zeros(a_fuzzy.shape, dtype=bool)  # don't use `zeros_like` - it also copies the dtype
    for ix in range(n_clusters):
        a_discrete[sample_assignment[ix], ix] = True

    # handle samples assigned to more than one cluster
    critical_samples = np.where(a_discrete.sum(1) > 1)[0]
    for sample_ix in critical_samples:
        if remove_overlap:
            a_discrete[sample_ix, :] = _one_hot(n_clusters)
        else:
            candidate_ixs = np.where(a_discrete[sample_ix, :])[0]
            most_likely_ix = candidate_ixs[np.argmax(a_fuzzy[sample_ix, list(a_discrete[sample_ix, :])])]
            a_discrete[sample_ix, :] = _one_hot(n_clusters, most_likely_ix)

    # check how many samples this left for each cluster
    n_samples_per_cluster = a_discrete.sum(0)
    if raise_threshold is not None and (n_samples_per_cluster < n_raise).any():
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
    """Create a pandas Series based on a one-hot encoded matrix.

    Parameters
    ----------
    membership
        One-hot encoded membership matrix, of shape `(n_samples x n_clusters)` i.e. a `1` in position `i, j`
        signifies that sample `i` belongs to cluster `j`.
    index
        Index for the Series. Careful, if this is not given, categories are removed when writing to AnnData.

    Returns
    -------
    Series, indicating cluster membership for each sample. The data type of the categories is :class:`str`
    and samples that belong to no cluster are assigned `NaN`.
    """
    n_samples, n_clusters = membership.shape
    if not isinstance(membership, np.ndarray):
        raise TypeError(f"Expected `membership` to be of type `numpy.ndarray`, found `{type(membership).__name__}`.")
    membership = np.asarray(membership)  # change the type in case a lineage object was passed.
    if membership.dtype != bool:
        raise TypeError(f"Expected `membership`'s elements to be boolean, found `{membership.dtype.name!r}`.")

    if not np.all(membership.sum(axis=1) <= 1):
        raise ValueError("Not all items are one-hot encoded or empty.")
    if (membership.sum(0) == 0).any():
        logg.warning(f"Detected {np.sum((membership.sum(0) == 0))} empty categories")

    if index is None:
        index = range(n_samples)
    if names is not None:
        if len(names) != n_clusters:
            raise ValueError(f"Shape mismatch, length of `names` is `{len(names)}`, but `n_clusters={n_clusters}`.")
    else:
        names = np.arange(n_clusters).astype("str")

    target_series = pd.Series(index=index, dtype="category")
    for vec, name in zip(membership.T, names):
        target_series = target_series.cat.add_categories(name)
        target_series.iloc[np.where(vec)[0]] = name

    return target_series


def _get_cat_and_null_indices(
    cat_series: pd.Series,
) -> Tuple[np.ndarray, np.ndarray, Dict[Any, np.ndarray]]:
    """Given a categorical :class:`~pandas.Series`, get the indices corresponding to categories and `NaNs`.

    Parameters
    ----------
    cat_series
        Series that contains categorical annotations.

    Returns
    -------
    - Array containing the indices of elements corresponding to categories in ``cat_series``.
    - Array containing the indices of elements corresponding to NaNs in ``cat_series``.
    - Dict containing categories of ``cat_series`` as keys and an array of corresponding indices as values.
    """
    # check the dtype
    if cat_series.dtype != "category":
        raise TypeError(f"Expected `cat_series` to be categorical, found `{cat_series.dtype.name!r}`.")

    # define a dict that has category names as keys and arrays of indices as values
    lookup_dict = {cat: np.where(cat_series == cat)[0] for cat in cat_series.cat.categories}
    all_indices = np.arange(len(cat_series))

    # collect all category indices
    cat_indices = np.concatenate(list(lookup_dict.values()))

    # collect all null indices (the ones where we have NaN in `cat_series`)
    null_indices = np.array(list(set(all_indices) - set(cat_indices)))

    # check that null indices and cat indices are unique
    assert (np.unique(cat_indices, return_counts=True)[1] == 1).all(), "Cat indices are not unique."
    assert (np.unique(null_indices, return_counts=True)[1] == 1).all(), "Null indices are not unique."

    # check that there is no overlap
    assert len(set(cat_indices).intersection(set(null_indices))) == 0, "Cat and null indices overlap."

    # check that their untion is the set of all indices
    assert set(cat_indices).union(set(null_indices)) == set(all_indices), "Some indices got lost on the way."

    return cat_indices, null_indices, lookup_dict


def _calculate_absorption_time_moments(
    Q: Union[np.ndarray, sp.spmatrix],
    trans_indices: np.ndarray,
    n: int,
    calculate_variance: bool = False,
    **kwargs: Any,
) -> Tuple[np.ndarray, Optional[np.ndarray]]:
    """Calculate the mean time until absorption and optionally its variance.

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
        Keyword arguments for :func:`cellrank._utils._linear_solver._solver_lin_system`.

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
        logg.debug("Calculating variance of mean time to absorption to any absorbing state")

        I = sp.eye(Q.shape[0]) if sp.issparse(Q) else np.eye(Q.shape[0])
        A_t = (I + Q).T
        B_t = (I - Q).T

        logg.debug("Solving equation (1/2)")
        X = _solve_lin_system(A_t, B_t, n_jobs=n_jobs, **kwargs).T
        y = m - X @ (m**2)

        logg.debug("Solving equation (2/2)")
        v = _solve_lin_system(X, y, use_eye=False, n_jobs=1, **solve_kwargs).squeeze()
        assert np.all(v >= 0), f"Encountered negative variance: `{v[v < 0]}`."

        var = np.zeros(n, dtype=np.float32)
        var[trans_indices] = v

    return mean, var


def _calculate_lineage_absorption_time_means(
    Q: sp.csr_matrix,
    R: sp.csr_matrix,
    trans_indices: np.ndarray,
    ixs: Dict[str, np.ndarray],
    index: pd.Index,
    calculate_variance: bool = False,
    **kwargs: Any,
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
    calculate_variance
        Whether to calculate variance of the mean time to absorption.
    kwargs
        Keyword arguments for :func:`cellrank._utils._linear_solver._solver_lin_system`.

    Returns
    -------
    A dataframe with means and optionally variances of the time to absorption for each lineage in ``ixs``.
    Uses more efficient implementation if compute the time for all lineages.
    """
    n = len(index)
    res = pd.DataFrame(index=index)

    if len(ixs) == 1:
        # use faster implementation in this case
        name = ", ".join(ixs.keys())
        res[f"{name}_mean"], var = _calculate_absorption_time_moments(
            Q,
            trans_indices,
            n,
            calculate_variance=calculate_variance,
            **kwargs,
        )
        if var is not None:
            res[f"{name}_var"] = var

        return res

    I = sp.eye(Q.shape[0]) if sp.issparse(Q) else np.eye(Q.shape)
    N_inv = I - Q

    logg.debug("Solving equation for `B`")
    B = _solve_lin_system(Q, R, use_eye=True, **kwargs)

    no_jobs_kwargs = kwargs.copy()
    _ = no_jobs_kwargs.pop("n_jobs", None)

    for i, (lineage, indices) in enumerate(ixs.items()):
        D_j = sp.diags(B[:, i])  # use `i`, since `B` is already aggregated
        D_j_inv = D_j.copy()
        D_j_inv.data = 1.0 / D_j.data

        # fmt: off
        logg.debug(f"Calculating mean time to absorption to `{lineage!r}`")
        m = _solve_lin_system(D_j_inv @ N_inv @ D_j, np.ones(Q.shape[0]), **kwargs).squeeze()
        # fmt: on

        mean = np.full(n, fill_value=np.inf, dtype=np.float64)
        mean[indices] = 0
        mean[trans_indices] = m
        res[f"{lineage}_mean"] = mean

        if calculate_variance:
            # fmt: off
            logg.debug(f"Calculating variance of the mean time to absorption to `{lineage!r}`")

            X = _solve_lin_system(D_j + Q @ D_j, N_inv @ D_j, use_eye=False, **kwargs)
            y = m - X @ (m**2)

            v = _solve_lin_system(X, y, use_eye=False, n_jobs=1, **no_jobs_kwargs).squeeze()
            assert np.all(v >= 0), f"Encountered negative variance: `{v[v < 0]}`."

            var = np.full(n, fill_value=np.nan, dtype=np.float64)
            var[indices] = 0
            var[trans_indices] = v
            res[f"{lineage}_var"] = var
            # fmt: on
    return res


def _check_collection(
    adata: AnnData,
    needles: Iterable[str],
    attr_name: str,
    key_name: str = "Gene",
    use_raw: bool = False,
    raise_exc: bool = True,
) -> List[str]:
    """Check if given collection contains all the keys.

    Parameters
    ----------
    adata
        Annotated data object.
    needles
        Keys to check.
    attr_name
        Attribute of ``adata`` where the needles are stored.
    key_name
        Pretty name of the key which will be displayed when error is found.
    use_raw
        Whether to access :attr:`anndata.AnnData.raw` or just ``adata``.

    Returns
    -------
    Nothing, but raises and :class:`KeyError` if one of the needles is not found.
    """
    adata_name = "adata"

    if use_raw and adata.raw is None:
        logg.warning("Argument `use_raw` was set to `True`, but no `raw` attribute is found. Ignoring")
        use_raw = False
    if use_raw:
        adata_name = "adata.raw"
        adata = adata.raw

    haystack, res = getattr(adata, attr_name), []
    for needle in needles:
        if needle not in haystack:
            if raise_exc:
                raise KeyError(f"{key_name} `{needle}` not found in `{adata_name}.{attr_name}`.")
        else:
            res.append(needle)

    return res


def _minmax(data: np.ndarray, perc: Optional[Tuple[float, float]] = None) -> Tuple[float, float]:
    """Return minimum and maximum value of the data.

    Parameters
    ----------
    data
        Values for which to return the minimum and maximum.
    perc
        If not `None`, clip the values by the percentiles before getting the result.

    Returns
    -------
    Minimum and maximum values, respectively.
    """
    if perc is not None:
        data = np.clip(data, *np.percentile(data, sorted(perc)))

    return float(np.nanmin(data)), float(np.nanmax(data))


def _modify_neigh_key(key: Optional[str]) -> str:
    if key in (None, "connectivities", "distances"):
        return "neighbors"

    if key.endswith("_connectivities"):
        key = key[:-15]
    elif key.endswith("_distances"):
        key = key[:-10]

    return key


def _get_neighs(adata: AnnData, mode: str = "distances", key: Optional[str] = None) -> Union[np.ndarray, sp.spmatrix]:
    if key is None:
        res = _read_graph_data(adata, mode)  # legacy behavior
    else:
        try:
            res = _read_graph_data(adata, key)
            assert isinstance(res, (np.ndarray, sp.spmatrix))
        except (KeyError, AssertionError):
            res = _read_graph_data(adata, f"{_modify_neigh_key(key)}_{mode}")

    if not isinstance(res, (np.ndarray, sp.spmatrix)):
        raise TypeError(f"Expected to find `numpy.ndarray` or `scipy.sparse.spmatrix`, found `{type(res)}`.")

    return res


def _has_neighs(adata: AnnData, key: Optional[str] = None) -> bool:
    return _modify_neigh_key(key) in adata.uns


def _get_neighs_params(adata: AnnData, key: str = "neighbors") -> Dict[str, Any]:
    return adata.uns.get(key, {}).get("params", {})


def _read_graph_data(adata: AnnData, key: str) -> Union[np.ndarray, sp.spmatrix]:
    """Read graph data from :class:`~anndata.AnnData`.

    Parameters
    ----------
    adata
        Annotated data object.
    key
        Key in :attr:`~anndata.AnnData.obsp`.

    Returns
    -------
    The graph data.
    """
    if key in adata.obsp:
        return adata.obsp[key]

    raise KeyError(f"Unable to find data in `adata.obsp[{key!r}]`.")


def valuedispatch(func):
    """Dispatch a function based on the first value."""
    registry = {}

    def dispatch(value):
        return registry.get(value, func)

    def register(value, func=None):
        if func is None:
            return lambda f: register(value, f)
        registry[value] = func
        return func

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return dispatch(args[0])(*args[1:], **kwargs)

    wrapper.register = register
    wrapper.dispatch = dispatch
    wrapper.registry = types.MappingProxyType(registry)

    functools.update_wrapper(wrapper, func)

    return wrapper


def _densify_squeeze(x: Union[sp.spmatrix, np.ndarray], dtype=np.float32) -> np.ndarray:
    if sp.issparse(x):
        x = x.toarray()
    # use np.array instead of asarray to create a copy
    x = np.array(x, dtype=dtype)
    if x.ndim == 2 and x.shape[1] == 1:
        x = np.squeeze(x, axis=1)

    return x


@contextlib.contextmanager
@d.dedent
def _gene_symbols_ctx(
    adata: AnnData,
    *,
    key: Optional[str] = None,
    use_raw: bool = False,
    make_unique: bool = False,
) -> AnnData:
    """Set gene names from a column in :attr:`~anndata.AnnData.var`.

    Parameters
    ----------
    %(adata)s
    key
        Key in :attr:`~anndata.AnnData.var` where the gene symbols are stored. If `None`, this operation is a no-op.
    use_raw
        Whether to change the gene names in :attr:`~anndata.AnnData.raw`.
    make_unique
        Whether to make the newly assigned gene names unique.

    Yields
    ------
    The same ``adata`` with modified :attr:`~anndata.AnnData.var_names`, depending on the ``use_raw``.
    """

    def key_present() -> bool:
        if use_raw:
            if adata.raw is None:
                raise AttributeError("No `.raw` attribute found. Try specifying `use_raw=False`.")
            return key in adata.raw.var
        return key in adata.var

    if key is None:
        yield adata
    elif not key_present():
        raise KeyError(f"Unable to find gene symbols in `adata.{'raw.' if use_raw else ''}var[{key!r}]`.")
    else:
        adata_orig = adata
        if use_raw:
            adata = adata.raw

        var_names = adata.var_names.copy()
        try:
            # TODO(michalk8): doesn't update varm (niche)
            adata.var.index = make_index_unique(adata.var[key]) if make_unique else adata.var[key]
            yield adata_orig
        finally:
            # in principle we assume the callee doesn't change the index
            # otherwise, would need to check whether it has been changed and add an option to determine what to do
            adata.var.index = var_names


@wrapt.decorator
def _genesymbols(wrapped: Callable[..., Any], instance: Any, args: Any, kwargs: Any) -> Any:
    with _gene_symbols_ctx(
        args[0],
        key=kwargs.pop("gene_symbols", None),
        make_unique=True,
        use_raw=kwargs.get("use_raw", False),
    ):
        return wrapped(*args, **kwargs)
