# -*- coding: utf-8 -*-
"""Cluster lineages module."""

from types import MappingProxyType
from typing import Dict, Tuple, Union, TypeVar, Optional, Sequence
from pathlib import Path
from collections import Iterable

import numpy as np
from sklearn.preprocessing import StandardScaler

import matplotlib.pyplot as plt

from cellrank import logging as logg
from cellrank.ul._docs import d
from cellrank.pl._utils import (
    _model_type,
    _get_backend,
    _callback_type,
    _create_models,
    _time_range_type,
    _create_callbacks,
)
from cellrank.tl._utils import save_fig, _unique_order_preserving
from cellrank.ul._utils import _get_n_cores, _check_collection
from cellrank.tl._constants import _DEFAULT_BACKEND, AbsProbKey
from cellrank.ul._parallelize import parallelize

AnnData = TypeVar("AnnData")
Queue = TypeVar("Queue")


def _cluster_lineages_helper(
    genes: Sequence[str],
    models: _model_type,
    callbacks: _callback_type,
    lineage: str,
    time_range: Optional[Union[float, Tuple[float, float]]],
    queue: Optional[Queue],
    **kwargs,
) -> np.ndarray:
    """
    Fit models to genes in given lineages. Used by :func:`cellrank.pl.cluster_lineage`.

    Parameters
    ----------
    genes
        Genes or observations for which to fit the models.
    models
        Gene and lineage specific models.
    callbacks
        Gene and lineage specific prepare callbacks.
    lineage
        Name of the lineage for which to fit the models.
    time_range
        Minimum and maximum pseudotime.
    queue
        Signalling queue in the parent process/thread used to update the progress bar.
    kwargs
        Keyword arguments for :meth:`cellrank.ul.models.BaseModel.prepare`.

    Returns
    -------
    :class:`numpy.ndarray`
        The model predictions, optionally normed.
    """

    res = []
    for gene in genes:
        cb = callbacks[gene][lineage]

        model = cb(
            models[gene][lineage],
            gene=gene,
            lineage=lineage,
            time_range=time_range,
            **kwargs,
        ).fit()
        res.append(model.predict())

        if queue is not None:
            queue.put(1)

    if queue is not None:
        queue.put(None)

    return np.array(res)


@d.dedent
def cluster_lineage(
    adata: AnnData,
    model: _model_type,
    genes: Sequence[str],
    lineage: str,
    backward: bool = False,
    time_range: _time_range_type = None,
    clusters: Optional[Sequence[str]] = None,
    n_points: int = 200,
    time_key: str = "latent_time",
    cluster_key: str = "clusters",
    norm: bool = True,
    recompute: bool = False,
    callback: _callback_type = None,
    ncols: int = 3,
    sharey: bool = False,
    key_added: Optional[str] = None,
    show_progress_bar: bool = True,
    n_jobs: Optional[int] = 1,
    backend: str = _DEFAULT_BACKEND,
    figsize: Optional[Tuple[float, float]] = None,
    dpi: Optional[int] = None,
    save: Optional[Union[str, Path]] = None,
    pca_kwargs: Dict = MappingProxyType({"svd_solver": "arpack"}),
    neighbors_kwargs: Dict = MappingProxyType({"use_rep": "X"}),
    louvain_kwargs: Dict = MappingProxyType({}),
    **kwargs,
) -> None:
    """
    Cluster gene expression trends within a lineage and plot the clusters.

    This function is based on Palantir, see [Setty19]_. It can be used to discover modules of genes that drive
    development along a given lineage. Consider running this function on a subset of genes which are potential lineage
    drivers, identified e.g. by running :func:`cellrank.tl.lineage_drivers`.

    .. image:: https://raw.githubusercontent.com/theislab/cellrank/master/resources/images/cluster_lineage.png
       :width: 400px
       :align: center

    Parameters
    ----------
    %(adata)s
    %(model)s
    %(genes)s
    lineage
        Name of the lineage for which to cluster the genes.
    %(backward)s
    %(time_ranges)s
    clusters
        Cluster identifiers to plot. If `None`, all clusters will be considered.
        Useful when plotting previously computed clusters.
    n_points
        Number of points used for prediction.
    time_key
        Key in ``adata.obs`` where the pseudotime is stored.
    cluster_key
        Key in ``adata.obs`` where the clustering is stored.
    norm
        Whether to z-normalize each trend to have zero mean, unit variance.
    recompute
        If `True`, recompute the clustering, otherwise try to find already existing one.
    %(model_callback)s
    ncols
        Number of columns for the plot.
    sharey
        Whether to share y-axis across multiple plots.
    key_added
        Postfix to add when saving the results to ``adata.uns``.
    %(parallel)s
    %(plotting)s
    pca_kwargs
        Keyword arguments for :func:`scanpy.pp.pca`.
    neighbors_kwargs
        Keyword arguments for :func:`scanpy.pp.neighbors`.
    louvain_kwargs
        Keyword arguments for :func:`scanpy.tl.louvain`.
    **kwargs:
        Keyword arguments for :meth:`cellrank.ul.models.BaseModel.prepare`.

    Returns
    -------
    %(just_plots)s

        Updates ``adata.uns`` with the following key:

            - ``lineage_{lineage}_trend_{key_added}`` - an :class:`anndata.AnnData` object
              of shape ``(n_genes, n_points)`` containing the clustered genes.
    """

    import scanpy as sc
    from anndata import AnnData as _AnnData

    lineage_key = str(AbsProbKey.BACKWARD if backward else AbsProbKey.FORWARD)
    if lineage_key not in adata.obsm:
        raise KeyError(f"Lineages key `{lineage_key!r}` not found in `adata.obsm`.")

    _ = adata.obsm[lineage_key][lineage]

    genes = _unique_order_preserving(genes)
    _check_collection(adata, genes, "var_names", kwargs.get("use_raw", False))

    key_to_add = f"lineage_{lineage}_trend"
    if key_added is not None:
        logg.debug(f"Adding key `{key_added!r}`")
        key_to_add += f"_{key_added}"

    if recompute or key_to_add not in adata.uns:
        kwargs["time_key"] = time_key  # kwargs for the model.prepare
        kwargs["n_test_points"] = n_points
        kwargs["backward"] = backward

        models = _create_models(model, genes, [lineage])
        callbacks = _create_callbacks(adata, callback, genes, [lineage])

        backend = _get_backend(model, backend)
        n_jobs = _get_n_cores(n_jobs, len(genes))

        start = logg.info(f"Computing gene trends using `{n_jobs}` core(s)")
        trends = parallelize(
            _cluster_lineages_helper,
            genes,
            as_array=True,
            unit="gene",
            n_jobs=n_jobs,
            backend=backend,
            extractor=np.vstack,
            show_progress_bar=show_progress_bar,
        )(models, callbacks, lineage, time_range, **kwargs)
        logg.info("    Finish", time=start)

        trends = trends.T
        if norm:
            logg.debug("Normalizing using `StandardScaler`")
            _ = StandardScaler(copy=False).fit_transform(trends)

        trends = _AnnData(trends.T)
        trends.obs_names = genes

        # sanity check
        if trends.n_obs != len(genes):
            raise RuntimeError(
                f"Expected to find `{len(genes)}` genes, found `{trends.n_obs}`."
            )
        if n_points is not None and trends.n_vars != n_points:
            raise RuntimeError(
                f"Expected to find `{n_points}` points, found `{trends.n_vars}`."
            )

        pca_kwargs = dict(pca_kwargs)
        n_comps = pca_kwargs.pop(
            "n_comps", min(50, kwargs.get("n_test_points"), len(genes)) - 1
        )  # default value

        sc.pp.pca(trends, n_comps=n_comps, **pca_kwargs)
        sc.pp.neighbors(trends, **neighbors_kwargs)

        louvain_kwargs = dict(louvain_kwargs)
        louvain_kwargs["key_added"] = cluster_key
        sc.tl.louvain(trends, **louvain_kwargs)

        adata.uns[key_to_add] = trends
    else:
        logg.info(f"Loading data from `adata.uns[{key_to_add!r}]`")
        trends = adata.uns[key_to_add]

    if clusters is None:
        if cluster_key not in trends.obs:
            raise KeyError(f"Invalid cluster key `{cluster_key!r}`.")
        clusters = trends.obs[cluster_key].cat.categories

    nrows = int(np.ceil(len(clusters) / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(ncols * 10, nrows * 10) if figsize is None else figsize,
        sharey=sharey,
        dpi=dpi,
    )

    if not isinstance(axes, Iterable):
        axes = [axes]
    axes = np.ravel(axes)

    j = 0
    for j, (ax, c) in enumerate(zip(axes, clusters)):  # noqa
        data = trends[trends.obs[cluster_key] == c].X
        mean, sd = np.mean(data, axis=0), np.var(data, axis=0)
        sd = np.sqrt(sd)

        for i in range(data.shape[0]):
            ax.plot(data[i], color="gray", lw=0.5)

        ax.plot(mean, lw=2, color="black")
        ax.plot(mean - sd, lw=1.5, color="black", linestyle="--")
        ax.plot(mean + sd, lw=1.5, color="black", linestyle="--")
        ax.fill_between(
            range(len(mean)), mean - sd, mean + sd, color="black", alpha=0.1
        )

        ax.set_title(f"Cluster {c}")
        ax.set_xticks([])

        if not sharey:
            ax.set_yticks([])

    for j in range(j + 1, len(axes)):
        axes[j].remove()

    if save is not None:
        save_fig(fig, save)
