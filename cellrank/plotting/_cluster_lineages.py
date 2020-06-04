# -*- coding: utf-8 -*-
"""Cluster lineages module."""

from types import MappingProxyType
from typing import Dict, Tuple, Union, Optional, Sequence
from pathlib import Path
from collections import Iterable

import numpy as np

import matplotlib.pyplot as plt

import scanpy as sc
from scanpy import logging as logg
from anndata import AnnData

from cellrank.tools._utils import save_fig
from cellrank.utils._utils import _get_n_cores, check_collection
from cellrank.plotting._utils import _model_type, _create_models, _is_any_gam_mgcv
from cellrank.tools._constants import LinKey
from cellrank.utils._parallelize import parallelize
from cellrank.utils.models._models import Model


def _cl_process(
    genes: Sequence[str],
    models: Dict[str, Dict[str, Model]],
    lineage_name: str,
    norm: str,
    queue,
    **kwargs,
) -> np.ndarray:
    """
    Fit models to genes in given lineages. Used by :func:`cellrank.pl.cluster_lineage`.

    Params
    ------
    genes
        Genes or observations for which to fit the models.
    models
        Gene and lineage specific models.
    lineage_name
        Name of the lineage for which to fit the models.
    norm
        Whether to z-normalize the fitted values.
    queue
        Signalling queue in the parent process/thread used to update the progress bar.
    kwargs
        Keyword arguments for :meth:`cellrank.ul.models.Model.prepare`.

    Returns
    -------
    :class:`numpy.ndarray`
        The model predictions, optionally normed.
    """

    res = []
    for gene in genes:
        model = models[gene][lineage_name].prepare(gene, lineage_name, **kwargs).fit()
        res.append(model.predict())
        queue.put(1)
    queue.put(None)

    res = np.squeeze(np.array(res))

    if not norm:
        return res

    mean = np.expand_dims(np.mean(res, 1), -1)
    sd = np.expand_dims(np.sqrt(np.var(res, 1)), -1)

    return (res - mean) / sd


def cluster_lineage(
    adata: AnnData,
    model: _model_type,
    genes: Sequence[str],
    lineage: str,
    final: bool = True,
    clusters: Optional[Sequence[str]] = None,
    n_points: int = 200,
    time_key: str = "latent_time",
    cluster_key: str = "louvain",
    norm: bool = True,
    recompute: bool = False,
    ncols: int = 3,
    sharey: bool = False,
    n_jobs: Optional[int] = 1,
    backend: str = "multiprocessing",
    pca_kwargs: Dict = MappingProxyType({"svd_solver": "arpack"}),
    neighbors_kwargs: Dict = MappingProxyType({"use_rep": "X"}),
    louvain_kwargs: Dict = MappingProxyType({}),
    key_added: Optional[str] = None,
    save: Optional[Union[str, Path]] = None,
    figsize: Optional[Tuple[float, float]] = None,
    dpi: Optional[int] = None,
    show_progress_bar: bool = True,
    **kwargs,
) -> None:
    """
    Cluster gene expression trends within a lineage and plot the clusters.

    This function is based on Palantir, see [Setty19]_. It can be used to discover modules of genes that drive
    development along a given lineage. Consider running this function on a subset of genes which are potential lineage
    drivers, identified e.g. by running :func:`cellrank.tl.gene_importance`.

    .. image:: https://raw.githubusercontent.com/theislab/cellrank/master/resources/images/cluster_lineage.png
       :width: 400px
       :align: center

    Params
    ------
    adata : :class:`anndata.AnnData`
        Annotated data object.
    model
        Model to fit.

        - If a :class:`dict`, gene and lineage specific models can be specified. Use `'*'` to indicate
        all genes or lineages, for example `{'Map2': {'*': ...}, 'Dcx': {'Alpha': ..., '*': ...}}`.
    genes
        Genes in :paramref:`adata`.var_names to cluster.
    lineage_name
        Name of the lineage along which to cluster the genes.
    final
        Whether to consider cells going to final states or vice versa.
    clusters
        Cluster identifiers to plot. If `None`, all clusters will be considered.
        Useful when plotting previously computed clusters.
    n_points
        Number of points used for prediction.
    time_key
        Key in :paramref:`adata` `.obs` where the pseudotime is stored.
    cluster_key
        Key in :paramref:`adata` `.obs` where the clustering is stored.
    norm
        Whether to z-normalize each trend to have `0` mean, `1` variance.
    recompute
        If `True`, recompute the clustering, otherwise try to find already existing one.
    ncols
        Number of columns for the plot.
    sharey
        Whether to share y-axis across multiple plots.
    n_jobs
        Number of parallel jobs. If `-1`, use all available cores. If `None` or `1`, the execution is sequential.
    backend
        Which backend to use for multiprocessing.
        See :class:`joblib.Parallel` for valid options.
    pca_kwargs
        Keyword arguments for :func:`scanpy.pp.pca`.
    neighbors_kwargs
        Keyword arguments for :func:`scanpy.pp.neighbors`.
    louvain_kwargs
        Keyword arguments for :func:`scanpy.tl.louvain`.
    save
        Filename where to save the plot.
        If `None`, just shows the plot.
    figsize
        Size of the figure. If `None`, it will be set automatically.
    dpi
        Dots per inch.
    show_progress_bar
        Whether to show a progress bar tracking models fitted.
    kwargs:
        Keyword arguments for :meth:`cellrank.ul.models.Model.prepare`.

    Returns
    -------
    None
        Plots the clusters of :paramref:`genes` for the given :paramref:`lineage_name`.
        Optionally saves the figure based on :paramref:`save`.

        Updates :paramref:`adata` `.uns` with the following key:

        lineage_{:paramref:`lineage_name`}_trend_{:paramref:`key_added`}:
            - :class:`anndata.AnnData` object of shape `len` (:paramref:`genes`) x :paramref:`n_points`
              containing the clustered genes.
    """

    lineage_key = str(LinKey.FORWARD if final else LinKey.BACKWARD)
    if lineage_key not in adata.obsm:
        raise KeyError(f"Lineages key `{lineage_key!r}` not found in `adata.obsm`.")

    _ = adata.obsm[lineage_key][lineage]

    check_collection(adata, genes, "var_names")

    key_to_add = f"lineage_{lineage}_trend"
    if key_added is not None:
        logg.debug(f"DEBUG: Adding key `{key_added!r}`")
        key_to_add += f"_{key_added}"

    if recompute or key_to_add not in adata.uns:
        kwargs["time_key"] = time_key  # kwargs for the model.prepare
        kwargs["n_test_points"] = n_points
        kwargs["final"] = final

        models = _create_models(model, genes, [lineage])
        if _is_any_gam_mgcv(models):
            backend = "multiprocessing"

        n_jobs = _get_n_cores(n_jobs, len(genes))

        start = logg.info(f"Computing trends using `{n_jobs}` core(s)")
        trends = parallelize(
            _cl_process,
            genes,
            as_array=True,
            unit="gene",
            n_jobs=n_jobs,
            backend=backend,
            show_progress_bar=show_progress_bar,
        )(models, lineage, norm, **kwargs)
        logg.info("    Finish", time=start)

        trends = AnnData(np.vstack(trends))
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
        n_comps = pca_kwargs.pop("n_comps", 50)  # default value
        if n_comps > len(genes):
            n_comps = len(genes) - 1

        sc.pp.pca(trends, n_comps=n_comps, **pca_kwargs)
        sc.pp.neighbors(trends, **neighbors_kwargs)

        louvain_kwargs = dict(louvain_kwargs)
        louvain_kwargs["key_added"] = cluster_key
        sc.tl.louvain(trends, **louvain_kwargs)

        adata.uns[key_to_add] = trends
    else:
        logg.info(f"Loading data from `adata.uns[{key_to_add}!r]`")
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
