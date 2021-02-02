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
    _fit_bulk,
    _get_backend,
    _callback_type,
    _create_models,
    _time_range_type,
    _create_callbacks,
    _input_model_type,
    _return_model_type,
)
from cellrank.tl._utils import save_fig, _unique_order_preserving
from cellrank.ul._utils import _get_n_cores, _check_collection
from cellrank.tl._constants import _DEFAULT_BACKEND, AbsProbKey

AnnData = TypeVar("AnnData")
Queue = TypeVar("Queue")


@d.dedent
def cluster_lineage(
    adata: AnnData,
    model: _input_model_type,
    genes: Sequence[str],
    lineage: str,
    backward: bool = False,
    time_range: _time_range_type = None,
    clusters: Optional[Sequence[str]] = None,
    n_points: int = 200,
    time_key: str = "latent_time",
    norm: bool = True,
    recompute: bool = False,
    callback: _callback_type = None,
    ncols: int = 3,
    sharey: Union[str, bool] = False,
    key: Optional[str] = None,
    random_state: Optional[int] = None,
    use_leiden: bool = False,
    show_progress_bar: bool = True,
    n_jobs: Optional[int] = 1,
    backend: str = _DEFAULT_BACKEND,
    figsize: Optional[Tuple[float, float]] = None,
    dpi: Optional[int] = None,
    save: Optional[Union[str, Path]] = None,
    pca_kwargs: Dict = MappingProxyType({"svd_solver": "arpack"}),
    neighbors_kwargs: Dict = MappingProxyType({"use_rep": "X"}),
    clustering_kwargs: Dict = MappingProxyType({}),
    return_models: bool = False,
    **kwargs,
) -> Optional[_return_model_type]:
    """
    Cluster gene expression trends within a lineage and plot the clusters.

    This function is based on Palantir, see [Setty19]_. It can be used to discover modules of genes that drive
    development along a given lineage. Consider running this function on a subset of genes which are potential
    lineage drivers, identified e.g. by running :func:`cellrank.tl.lineage_drivers`.

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
        Cluster identifiers to plot. If `None`, all clusters will be considered. Useful when
        plotting previously computed clusters.
    n_points
        Number of points used for prediction.
    time_key
        Key in ``adata.obs`` where the pseudotime is stored.
    norm
        Whether to z-normalize each trend to have zero mean, unit variance.
    recompute
        If `True`, recompute the clustering, otherwise try to find already existing one.
    %(model_callback)s
    ncols
        Number of columns for the plot.
    sharey
        Whether to share y-axis across multiple plots.
    key
        Key in ``adata.uns`` where to save the results. If `None`, it will be saved as ``lineage_{lineage}_trend`` .
    random_state
        Random seed for reproducibility.
    use_leiden
        Whether to use :func:`scanpy.tl.leiden` for clustering or :func:`scanpy.tl.louvain`.
    %(parallel)s
    %(plotting)s
    pca_kwargs
        Keyword arguments for :func:`scanpy.pp.pca`.
    neighbors_kwargs
        Keyword arguments for :func:`scanpy.pp.neighbors`.
    clustering_kwargs
        Keyword arguments for :func:`scanpy.tl.louvain` or :func:`scanpy.tl.leiden`.
    %(return_models)s
    kwargs:
        Keyword arguments for :meth:`cellrank.ul.models.BaseModel.prepare`.

    Returns
    -------
    %(plots_or_returns_models)s

        Also updates ``adata.uns`` with the following:

            - ``key`` or ``lineage_{lineage}_trend`` - an :class:`anndata.AnnData` object of
              shape `(n_genes, n_points)` containing the clustered genes.
    """

    import scanpy as sc
    from anndata import AnnData as _AnnData

    lineage_key = str(AbsProbKey.BACKWARD if backward else AbsProbKey.FORWARD)
    if lineage_key not in adata.obsm:
        raise KeyError(f"Lineages key `{lineage_key!r}` not found in `adata.obsm`.")

    _ = adata.obsm[lineage_key][lineage]

    genes = _unique_order_preserving(genes)
    _check_collection(adata, genes, "var_names", kwargs.get("use_raw", False))

    if key is None:
        key = f"lineage_{lineage}_trend"

    if recompute or key not in adata.uns:
        kwargs["backward"] = backward
        kwargs["time_key"] = time_key
        kwargs["n_test_points"] = n_points
        models = _create_models(model, genes, [lineage])
        all_models, models, genes, _ = _fit_bulk(
            models,
            _create_callbacks(adata, callback, genes, [lineage], **kwargs),
            genes,
            lineage,
            time_range,
            return_models=True,  # always return (better error messages)
            filter_all_failed=True,
            parallel_kwargs={
                "show_progress_bar": show_progress_bar,
                "n_jobs": _get_n_cores(n_jobs, len(genes)),
                "backend": _get_backend(models, backend),
            },
            **kwargs,
        )

        # `n_genes, n_test_points`
        trends = np.vstack([model[lineage].y_test for model in models.values()]).T

        if norm:
            logg.debug("Normalizing trends")
            _ = StandardScaler(copy=False).fit_transform(trends)

        trends = _AnnData(trends.T)
        trends.obs_names = genes

        # sanity check
        if trends.n_obs != len(genes):
            raise RuntimeError(
                f"Expected to find `{len(genes)}` genes, found `{trends.n_obs}`."
            )
        if trends.n_vars != n_points:
            raise RuntimeError(
                f"Expected to find `{n_points}` points, found `{trends.n_vars}`."
            )

        random_state = np.random.mtrand.RandomState(random_state).randint(2 ** 16)

        pca_kwargs = dict(pca_kwargs)
        pca_kwargs.setdefault("n_comps", min(50, n_points, len(genes)) - 1)
        pca_kwargs.setdefault("random_state", random_state)
        sc.pp.pca(trends, **pca_kwargs)

        neighbors_kwargs = dict(neighbors_kwargs)
        neighbors_kwargs.setdefault("random_state", random_state)
        sc.pp.neighbors(trends, **neighbors_kwargs)

        clustering_kwargs = dict(clustering_kwargs)
        clustering_kwargs["key_added"] = "clusters"
        clustering_kwargs.setdefault("random_state", random_state)
        try:
            if use_leiden:
                sc.tl.leiden(trends, **clustering_kwargs)
            else:
                sc.tl.louvain(trends, **clustering_kwargs)
        except ImportError as e:
            logg.warning(str(e))
            if use_leiden:
                sc.tl.louvain(trends, **clustering_kwargs)
            else:
                sc.tl.leiden(trends, **clustering_kwargs)

        logg.info(f"Saving data to `adata.uns[{key!r}]`")
        adata.uns[key] = trends
    else:
        all_models = None
        logg.info(f"Loading data from `adata.uns[{key!r}]`")
        trends = adata.uns[key]

    if "clusters" not in trends.obs:
        raise KeyError("Unable to find the clustering in `trends.obs['clusters']`.")

    if clusters is None:
        clusters = trends.obs["clusters"].cat.categories
    for c in clusters:
        if c not in trends.obs["clusters"].cat.categories:
            raise ValueError(
                f"Invalid cluster name `{c!r}`. "
                f"Valid options are `{list(trends.obs['clusters'].cat.categories)}`."
            )

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
        data = trends[trends.obs["clusters"] == c].X
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

    if return_models:
        return all_models
