# -*- coding: utf-8 -*-
"""Exact permutation test for Markov chains."""

from typing import Dict, List, Tuple, Union, TypeVar, Callable, Optional, Sequence

import numpy as np
from scipy.stats import chi2_contingency
from pandas.api.types import infer_dtype, is_categorical_dtype
from scipy.spatial.distance import euclidean

import matplotlib.pyplot as plt

from cellrank import logging as logg
from cellrank.ul._docs import d
from cellrank.tl._constants import AbsProbKey

AnnData = TypeVar("AnnData")


def _get_counts(pd: Union[np.ndarray, List[float]], n: int) -> List[float]:
    """
    Generate a list of counts that follows a given probability distribution.

    Parameters
    ----------
    pd
        Probability distribution used to generate the samples.
    n
        Total number of samples.

    Returns
    -------
    :class:`list`
        The counts.
    """

    lst = list(np.random.choice(np.arange(len(pd)), n, p=pd))
    freq = [lst.count(i) for i in np.arange(len(pd))]  # list of counts
    freq = [
        i if i > 0 else i + 1e-5 for i in freq
    ]  # replace end points with zero counts by a small number

    return freq


@d.dedent
def _permutation_test(
    adata: AnnData,
    cluster_key: str,
    cluster1: str,
    cluster2: str,
    backward: bool = False,
    dist_measure: Callable[[Sequence[float], Sequence[float]], float] = euclidean,
    n_perms: int = 1000,
    use_counts: bool = False,
    n_counts: int = 1000,
    n_bins: int = 200,
    seed: Optional[int] = None,
    **kwargs,
) -> Tuple[List[float], float, float]:
    """
    Permutation test implemented for both probability distributions and count distributions.

    Get as input two clusters, then calculate its average probability distribution and calculate the distance
    between both averages. Then permute the elements on the clusters and repeat the process.

    Parameters
    ----------
    %(adata)s
    cluster_key
        Key in ``adata.obs`` that contains the clusters.
    cluster1
        Name of the first cluster to compare.
    cluster2
        Name of the second cluster to compare.
    %(backward)s
    dist_measure : :class:`Callable`
        Distance measure to use. Default is `'euclidean'`.
    n_perms
        Number of permutations to perform.
    use_counts
        Whether to use counts distribution or probability distribution.
    n_counts
        Number of total counts used to calculate the counts distribution.
    n_bins
        Number of bins for the histogram.
    seed
        Random seed.
    **kwargs
        Keyword arguments for ``dist_measure``.

    Returns
    -------
    :class:`list`, :class:`float`, :class:`float`
        List containing the distance value between average of the probability distributions of both clusters,
        the observed distance and the corresponding p-value.
    """

    if cluster_key not in adata.obs:
        raise ValueError(f"Cluster key `{cluster_key!r}` not found in `adata.obs`.")

    if not is_categorical_dtype(adata.obs[cluster_key]):
        raise TypeError(
            f"Expected `adata.obs[{cluster_key!r}]` to be categorical, "
            f"found `{infer_dtype(adata.obs[cluster_key])}`."
        )

    if cluster1 not in adata.obs[cluster_key].cat.categories:
        raise ValueError(f"Cluster `{cluster1!r}` is not a valid cluster.")

    if cluster2 not in adata.obs[cluster_key].cat.categories:
        raise ValueError(f"Cluster `{cluster2!r}` is not a valid cluster.")

    if n_perms <= 0:
        raise ValueError(f"Argument `n_perms` must be positive, found `{n_perms}`.")

    if n_counts <= 0:
        raise ValueError(f"Argument `n_counts` must be positive, found `{n_counts}`.")

    start = logg.info("Starting exact permutation test")
    np.random.seed(seed)

    # consider the two possible directions
    lin_key = str(AbsProbKey.BACKWARD if backward else AbsProbKey.FORWARD)

    xs = adata[adata.obs[cluster_key] == cluster1].obsm[lin_key]
    ys = adata[adata.obs[cluster_key] == cluster2].obsm[lin_key]

    diff = []  # list of distances
    n = len(xs)
    xs_av = np.nanmean(xs, axis=0)  # average of both distributions ignoring nan values
    ys_av = np.nanmean(ys, axis=0)

    if use_counts:
        logg.debug("Using counts distribution")
        freq_x = _get_counts(xs_av, n_counts)
        freq_y = _get_counts(ys_av, n_counts)
        dist = dist_measure(freq_x, freq_y, **kwargs)
    else:
        logg.debug("Using probability distribution")
        dist = dist_measure(xs_av, ys_av, **kwargs)  # Distance between both averages
    diff.append(dist)

    zs = np.concatenate(
        (xs, ys), axis=0
    )  # create a extended list with all the distributions

    for _ in range(n_perms):
        np.random.shuffle(zs)  # randomly permute the clusters
        xs_sh = np.nanmean(zs[:n], axis=0)  # cluster 1
        ys_sh = np.nanmean(zs[n:], axis=0)  # cluster 2

        if use_counts:
            freq_x = _get_counts(xs_sh, n_counts)
            freq_y = _get_counts(ys_sh, n_counts)
            dist = dist_measure(freq_x, freq_y, **kwargs)
        else:
            dist = dist_measure(xs_sh, ys_sh, **kwargs)

        diff.append(dist)  # calculate distance

    # plot a histogram of the distances and the CDF
    gs = plt.GridSpec(nrows=1, ncols=2, figure=plt.figure(None, (10, 5)))

    plt.subplot(gs[0])
    plt.hist(diff, density=True, bins=n_bins)
    plt.title("Histogram")

    plt.subplot(gs[1])
    cum = plt.hist(diff, density=True, bins=n_bins, cumulative=True)
    plt.title("CDF")

    # calculate the pvalue
    thr = [n for n, i in enumerate(cum[1]) if i >= diff[0]][0]
    p_value = 1 - cum[0][thr - 1]
    logg.info("    Finish", time=start)

    return diff[1:], diff[0], p_value


@d.dedent
def _counts(
    adata: AnnData,
    cluster_key: str,
    clusters: Optional[List[str]] = None,
    n_samples: int = 1000,
    backward: bool = False,
) -> Dict[str, List[float]]:
    """
    Calculate the counts for each endpoint per cluster.

    Randomly choose with replacement *n* cells in each cluster and for each choice samples one point using its
    probability distribution. Then count the occurrences per each point.

    Parameters
    ----------
    %(adata)s
    cluster_key
        Key in ``adata.obs`` to access the cluster names.
    clusters
        List of clusters to consider. If `None`, all clusters are considered.
    n_samples
        Number of cells to sample from per cluster.
    %(backward)s

    Returns
    -------
    :class:`dict`
        Dictionary with keys as cluster names and values as lists with counts per endpoint.
    """

    if cluster_key not in adata.obs:
        raise KeyError(f"Unable to find cluster key `{cluster_key!r}` in `adata.obs`.")

    if n_samples <= 0:
        raise ValueError(f"Number of samples must be positive, found `{n_samples}`.")

    if clusters is not None:
        for cname in clusters:
            if cname not in adata.obs[cluster_key].cat.categories:
                raise ValueError(f"Key `{cname!r}` is not a categorical observation.")
        cluster_names = clusters
    else:
        cluster_names = adata.obs[cluster_key].cat.categories

    # Consider the two possible directions
    lin_key = str(AbsProbKey.BACKWARD if backward else AbsProbKey.FORWARD)
    if lin_key not in adata.obsm:
        raise KeyError(f"Lineages key `{lin_key!r}` not found in `adata.obsm`.")

    logg.debug("Calculating counts distribution of endpoint per cluster")
    res = {}

    for name in cluster_names:
        data = adata[adata.obs[cluster_key] == name].obsm[lin_key].X
        dim = data.shape[1]
        index = np.random.randint(data.shape[0], size=n_samples)
        lst = [
            np.random.choice(np.arange(dim), 1, p=np.atleast_1d(data[ind]))[0]
            if ~np.isnan(data[ind]).any()
            else 0
            for ind in index
        ]
        freq = (lst.count(i) for i in range(dim))
        res[name] = [i if i > 0 else 1e-5 for i in freq]

    return res


def _cramers_v(x: List[float], y: List[float]) -> float:
    """
    Calculate Cramer's V statistic for categorical-categorical association.

    Uses correction from **Bergsma and Wicher, Journal of the Korean Statistical Society 42 (2013): 323-328**.
    This is a symmetric coefficient: :math:`V(x,y) = V(y,x)`.

    Parameters
    ----------
    x
        A sequence of categorical measurements.
    y
        A sequence of categorical measurements.

    Returns
    -------
    float
        Cramer's v statistic for the ``x`, ``y` association.
    """

    #  Taken from: https://github.com/shakedzy/dython
    #  Original function taken from: https://stackoverflow.com/a/46498792/5863503
    #  Wikipedia: https://en.wikipedia.org/wiki/Cram%C3%A9r%27s_V

    confusion_matrix = np.array([x, y])
    chi2 = chi2_contingency(confusion_matrix)[0]
    n = confusion_matrix.sum().sum()
    phi2 = chi2 / n
    r, k = confusion_matrix.shape
    phi2corr = max(0, phi2 - ((k - 1) * (r - 1)) / (n - 1))
    rcorr = r - ((r - 1) ** 2) / (n - 1)
    kcorr = k - ((k - 1) ** 2) / (n - 1)

    return np.sqrt(phi2corr / min((kcorr - 1), (rcorr - 1)))
