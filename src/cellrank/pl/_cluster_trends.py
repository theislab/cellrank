import pathlib
import types
from typing import Any, Dict, Optional, Sequence, Tuple, Union

import numpy as np
from sklearn.preprocessing import StandardScaler

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, is_color_like
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

import scanpy as sc
from anndata import AnnData

from cellrank import logging as logg
from cellrank._utils import Lineage
from cellrank._utils._docs import d
from cellrank._utils._enum import DEFAULT_BACKEND, Backend_t
from cellrank._utils._parallelize import _get_n_cores
from cellrank._utils._utils import (
    _check_collection,
    _genesymbols,
    _unique_order_preserving,
    save_fig,
)
from cellrank.pl._utils import (
    _callback_type,
    _create_callbacks,
    _create_models,
    _fit_bulk,
    _get_backend,
    _get_sorted_colors,
    _input_model_type,
    _return_model_type,
    _time_range_type,
)

__all__ = ["cluster_trends"]


@d.dedent
@_genesymbols
def cluster_trends(
    adata: AnnData,
    model: _input_model_type,
    genes: Sequence[str],
    lineage: str,
    time_key: str,
    backward: bool = False,
    time_range: _time_range_type = None,
    clusters: Optional[Sequence[str]] = None,
    n_points: int = 200,
    covariate_key: Optional[Union[str, Sequence[str]]] = None,
    ratio: float = 0.05,
    cmap: Optional[str] = "viridis",
    norm: bool = True,
    recompute: bool = False,
    callback: _callback_type = None,
    ncols: int = 3,
    sharey: Union[str, bool] = False,
    key: Optional[str] = None,
    random_state: Optional[int] = None,
    show_progress_bar: bool = True,
    n_jobs: Optional[int] = 1,
    backend: Backend_t = DEFAULT_BACKEND,
    figsize: Optional[Tuple[float, float]] = None,
    dpi: Optional[int] = None,
    save: Optional[Union[str, pathlib.Path]] = None,
    pca_kwargs: Dict = types.MappingProxyType({"svd_solver": "arpack"}),
    neighbors_kwargs: Dict = types.MappingProxyType({"use_rep": "X"}),
    clustering_kwargs: Dict = types.MappingProxyType({}),
    return_models: bool = False,
    **kwargs: Any,
) -> Optional[_return_model_type]:
    """Cluster and plot gene expression trends within a lineage.

    .. seealso::
        - See :doc:`../../../notebooks/tutorials/estimators/800_gene_trends` on how to
          visualize the gene trends.

    This function is based on *Palantir* :cite:`setty:19`. It can be used to discover modules of genes that drive
    development along a given lineage. Consider running this function on a subset of genes which are potential
    lineage drivers.

    Parameters
    ----------
    %(adata)s
    %(model)s
    %(genes)s
    lineage
        Name of the lineage for which to cluster the genes.
    time_key
        Key in :attr:`~anndata.AnnData.obs` where the pseudotime is stored.
    %(backward)s
    %(time_range)s
    clusters
        Cluster identifiers to plot. If :obj:`None`, all clusters will be considered. Useful when
        plotting previously computed clusters.
    n_points
        Number of points used for prediction.
    covariate_key
        Keys in :attr:`~anndata.AnnData.obs` containing observations to be plotted at the bottom of each plot.
    %(gene_symbols)s
    ratio
        Height ratio of each covariate in ``covariate_key``.
    cmap
        Colormap to use for continuous covariates in ``covariate_key``.
    norm
        Whether to z-normalize each trend to have zero mean, unit variance.
    recompute
        If :obj:`True`, recompute the clustering, otherwise try to find already existing one.
    %(model_callback)s
    ncols
        Number of columns for the plot.
    sharey
        Whether to share y-axis across multiple plots.
    key
        Key in :attr:`~anndata.AnnData.uns` where to save the results.
        If :obj:`None`, it will be saved as ``'lineage_{lineage}_trend'`` .
    random_state
        Random seed for reproducibility.
    %(parallel)s
    %(plotting)s
    pca_kwargs
        Keyword arguments for :func:`~scanpy.pp.pca`.
    neighbors_kwargs
        Keyword arguments for :func:`~scanpy.pp.neighbors`.
    clustering_kwargs
        Keyword arguments for :func:`~scanpy.tl.leiden`.
    %(return_models)s
    kwargs
        Keyword arguments for :meth:`~cellrank.models.BaseModel.prepare`.

    Returns
    -------
    %(plots_or_returns_models)s Also updates :attr:`adata.uns <anndata.AnnData.uns>` with the following:

    - ``key`` or ``'lineage_{lineage}_trend'`` - :class:`~anndata.AnnData` object of
      shape ``(n_genes, n_points)`` containing the clustered genes.
    """

    def plot_cluster(row: int, col: int, cluster: str, sharey_ax: Optional[str] = None) -> Optional[plt.Axes]:
        gss = GridSpecFromSubplotSpec(
            row_delta,
            1,
            subplot_spec=gs[row : row + row_delta, col],
            hspace=0,
            height_ratios=[1.0] + [ratio] * (row_delta - 1),
        )
        ax = fig.add_subplot(gss[0, 0], sharey=sharey_ax)

        data = trends[trends.obs["clusters"] == c].X
        mean, sd = np.mean(data, axis=0), np.std(data, axis=0)

        for i in range(data.shape[0]):
            ax.plot(data[i], color="gray", lw=0.5)

        ax.plot(mean, lw=2, color="black")
        ax.plot(mean - sd, lw=1.5, color="black", linestyle="--")
        ax.plot(mean + sd, lw=1.5, color="black", linestyle="--")
        ax.fill_between(range(len(mean)), mean - sd, mean + sd, color="black", alpha=0.1)

        ax.set_title(f"cluster {cluster}")
        ax.set_xticks([])
        if sharey:
            ax.set_yticks([])
        ax.margins(0)

        if covariate_colors is not None:
            for i, colors in enumerate(covariate_colors):
                ax_clusters = fig.add_subplot(gss[i + 1, 0])
                if is_color_like(colors[0]):  # e.g. categorical
                    cm = ListedColormap(colors, N=len(colors))
                    ax_clusters.imshow(np.arange(cm.N)[None, :], cmap=cm, aspect="auto")
                else:
                    cm = plt.get_cmap(cmap)
                    ax_clusters.imshow(colors[None, :], cmap=cm, aspect="auto")
                ax_clusters.set_xticks([])
                ax_clusters.set_yticks([])

            ax_clusters.set_xticks(np.linspace(0, len(colors), 5))
            ax_clusters.set_xticklabels([f"{v:.3f}" for v in np.linspace(tmin, tmax, 5)])

        return ax if sharey else None

    use_raw = kwargs.get("use_raw", False)
    _ = Lineage.from_adata(adata, backward=backward)[lineage]  # sanity check

    genes = _unique_order_preserving(genes)
    _check_collection(adata, genes, "var_names", use_raw=use_raw)

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

        mod = next(mod for tmp in all_models.values() for mod in tmp.values())
        trends = AnnData(trends.T)
        trends.obs_names = genes
        trends.var["x_test"] = x_test = mod.x_test

        # sanity check
        if trends.n_obs != len(genes):
            raise RuntimeError(f"Expected to find `{len(genes)}` genes, found `{trends.n_obs}`.")
        if trends.n_vars != n_points:
            raise RuntimeError(f"Expected to find `{n_points}` points, found `{trends.n_vars}`.")

        random_state = np.random.RandomState(random_state).randint(2**16)

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
        sc.tl.leiden(trends, **clustering_kwargs)

        logg.info(f"Saving data to `adata.uns[{key!r}]`")
        adata.uns[key] = trends
    else:
        all_models = None
        logg.info(f"Loading data from `adata.uns[{key!r}]`")
        trends = adata.uns[key]
        x_test = trends.var["x_test"]

    if "clusters" not in trends.obs:
        raise KeyError("Unable to find the clustering in `trends.obs['clusters']`.")

    if clusters is None:
        clusters = trends.obs["clusters"].cat.categories
    for c in clusters:
        if c not in trends.obs["clusters"].cat.categories:
            raise ValueError(
                f"Invalid cluster name `{c!r}`. " f"Valid options are `{list(trends.obs['clusters'].cat.categories)}`."
            )

    nrows = int(np.ceil(len(clusters) / ncols))
    fig = plt.figure(
        dpi=dpi,
        figsize=(ncols * 10, nrows * 10) if figsize is None else figsize,
        tight_layout=True,
    )

    if covariate_key is None:
        covariate_colors, row_delta = None, 1
    else:
        tmin, tmax = np.min(x_test), np.max(x_test)
        covariate_colors = _get_sorted_colors(
            adata,
            covariate_key,
            time_key,
            tmin=tmin,
            tmax=tmax,
        )
        row_delta = len(covariate_colors) + 1
    gs = GridSpec(nrows=nrows * row_delta, ncols=ncols, figure=fig)

    row, sharey_ax = -row_delta, None
    for i, c in enumerate(clusters):
        if i % ncols == 0:
            row += row_delta
        sharey_ax = plot_cluster(row, i % ncols, c, sharey_ax=sharey_ax)

    if save is not None:
        save_fig(fig, save)

    if return_models:
        return all_models
