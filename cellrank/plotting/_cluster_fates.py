# -*- coding: utf-8 -*-
from collections import OrderedDict as odict
from math import ceil
from types import MappingProxyType
from typing import Optional, Sequence, Tuple, List, Mapping, Any, Union

import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from anndata import AnnData
from scanpy import logging as logg

from cellrank.tools._cluster_fates import _cramers_v, _counts
from cellrank.tools._constants import LinKey
from cellrank.tools._utils import save_fig
from cellrank.utils._utils import _make_unique
from cellrank.tools._lineage import Lineage

_cluster_fates_modes = ("bar", "paga", "paga_pie", "violin")


def cluster_fates(
    adata: AnnData,
    cluster_key: Optional[str] = "louvain",
    clusters: Optional[Union[str, Sequence[str]]] = None,
    lineages: Optional[Union[str, Sequence[str]]] = None,
    mode: str = "bar",
    final: bool = True,
    show_cbar: bool = True,
    ncols: Optional[int] = None,
    sharey: bool = False,
    save: Optional[str] = None,
    legend_kwargs: Mapping[str, Any] = MappingProxyType({"loc": "best"}),
    figsize: Optional[Tuple[float, float]] = None,
    dpi: Optional[int] = None,
    **kwargs,
) -> None:
    """
    Produces plots that aggregate lineage probabilities to a cluster level.

    This can be used to investigate how likely a certain cluster is to go to one of the endpoints, or in turn to have
    descended from one of the starting points. For mode `paga` and `paga_pie`, we use PAGA, see [Wolf19]_.

    .. image:: https://raw.githubusercontent.com/theislab/cellrank/master/resources/images/cluster_fates.png
       :width: 400px
       :align: center

    Params
    ------
    adata : :class:`anndata.AnnData`
        Annotated data object.
    cluster_key
        Key in :paramref:`adata` `.obs` containing the clusters.
    clusters
        Clusters to visualize.
        If `None`, all clusters will be plotted.
    lineages
        Lineages for which to visualize absorption probabilities.
        If `None`, use all available lineages.
    mode
        Type of plot to show.

        - If `'bar'`, plot barplots for specified :paramref:`clusters` and :paramref:`endpoints`.
        - If `'paga'`, plot `N` :func:`scanpy.pl.paga` plots, one for each endpoint in :paramref:`endpoints`.
        - If `'paga_pie'`, visualize absorption probabilities as a pie chart for each cluster
          for the given :paramref:`endpoints`.
        - If `'violin'`, group the data by lineages and plot the fate distribution per cluster.

        Best for looking at the distribution of fates within one cluster.
    dpi
        Dots per inch.
    final
        Whether to consider cells going to final states or vice versa.
    show_cbar
        Whether to show colorbar when :paramref:`mode` is `'paga_pie'`.
    ncols
        Number of columns when :paramref:`mode` is `'bar'` or `'paga'`.
    sharey
        Whether to share y-axis when :paramref:`mode` is `'bar'`.
    figsize
        Size of the figure.
    save
        Filename where to save the plots.
        If `None`, just shows the plot.
    legend_kwargs
        Keyword arguments for :func:`matplotlib.axes.Axes.legend`, such as `'loc'` for legend position.
    figsize
        Size of the figure. If `None`, it will be set automatically.
    dpi
        Dots per inch.
    kwargs
        Keyword arguments for :func:`scanpy.pl.paga`, :func:`scanpy.pl.violin` or :func:`matplotlib.pyplot.bar`,
        depending on :paramref:`mode`.

    Returns
    -------
    None
        Nothing, just plots the fates for specified :paramref:`clusters` and :paramref:`endpoints`.
        Optionally saves the figure based on :paramref:`save`.
    """

    def plot_bar():
        cols = 4 if ncols is None else ncols
        n_rows = ceil(len(clusters) / cols)
        fig = plt.figure(
            None, (3.5 * cols, 5 * n_rows) if figsize is None else figsize, dpi=dpi
        )
        fig.tight_layout()

        gs = plt.GridSpec(n_rows, cols, figure=fig, wspace=0.7, hspace=0.9)

        ax = None
        colors = list(adata.obsm[lk][:, lin_names].colors)

        for g, k in zip(gs, d.keys()):
            current_ax = fig.add_subplot(g, sharey=ax)
            current_ax.bar(
                x=np.arange(len(lin_names)),
                height=d[k][0],
                color=colors,
                yerr=d[k][1],
                ecolor="black",
                capsize=10,
                **kwargs,
            )
            if sharey:
                ax = current_ax

            current_ax.set_xticks(np.arange(len(lin_names)))
            current_ax.set_xticklabels(lin_names, rotation="vertical")
            current_ax.set_title(k)
            current_ax.set_xlabel("Endpoints" if final else "Startpoints")
            current_ax.set_ylabel("Absorption probability")

        return fig

    def plot_paga():
        kwargs["save"] = None
        kwargs["show"] = False
        if "cmap" not in kwargs:
            kwargs["cmap"] = cm.viridis

        cols = len(lin_names) if ncols is None else ncols
        nrows = ceil(len(lin_names) / cols)
        fig, axes = plt.subplots(
            nrows,
            cols,
            figsize=(6 * cols, 4 * nrows) if figsize is None else figsize,
            constrained_layout=True,
            dpi=dpi,
        )

        i = 0
        axes = [axes] if not isinstance(axes, np.ndarray) else np.ravel(axes)
        vmin, vmax = np.inf, -np.inf

        for i, (ax, lineage_name) in enumerate(zip(axes, lin_names)):
            colors = [v[0][i] for v in d.values()]
            vmin, vmax = np.nanmin(colors + [vmin]), np.nanmax(colors + [vmax])
            kwargs["ax"] = ax
            kwargs["colors"] = tuple(colors)
            kwargs["title"] = lineage_name

            sc.pl.paga(adata, **kwargs)

        if show_cbar:
            norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
            cax, _ = mpl.colorbar.make_axes(ax, aspect=100)  # new matplotlib feature
            _ = mpl.colorbar.ColorbarBase(
                cax, norm=norm, cmap=kwargs["cmap"], label="Absorption probability"
            )

        for ax in axes[i + 1 :]:
            ax.remove()

        return fig

    def plot_paga_pie():
        colors = list(adata.obsm[lk][:, lin_names].colors)
        colors = {i: odict(zip(colors, mean)) for i, (mean, _) in enumerate(d.values())}

        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

        kwargs["ax"] = ax
        kwargs["show"] = False
        kwargs["colorbar"] = False  # has to be disabled

        kwargs["colors"] = colors
        kwargs.pop("save", None)  # we will handle saving

        sc.pl.paga(adata, **kwargs)
        dummy_pos = adata.uns["paga"]["pos"][0]

        if legend_kwargs.get("loc", None) is not None:
            for lineage_name, color in zip(lin_names, colors[0].keys()):
                ax.plot(*dummy_pos, label=lineage_name, color=color)
            if len(colors[0].keys()) != len(adata.obsm[lk].names):
                ax.plot(*dummy_pos, label="Rest", color="grey")

            ax.legend(**legend_kwargs)

        return fig

    def plot_violin():
        kwargs["show"] = False
        kwargs.pop("ax", None)
        kwargs.pop("keys", None)
        kwargs.pop("save", None)  # we will handle saving
        kwargs["groupby"] = cluster_key
        if kwargs.get("rotation", None) is None:
            kwargs["rotation"] = 90

        data = adata.obsm[lk]
        to_clean = []

        for i, name in enumerate(lin_names):
            if name not in adata.obs_keys():
                to_clean.append(name)
                adata.obs[name] = np.array(
                    data[:, name]
                )  # TODO: better approach - dummy adata

        cols = len(lin_names) if ncols is None else ncols
        nrows = ceil(len(lin_names) / cols)
        fig, axes = plt.subplots(
            nrows,
            cols,
            figsize=(6 * cols, 4 * nrows) if figsize is None else figsize,
            dpi=dpi,
        )
        if not isinstance(axes, np.ndarray):
            axes = [axes]
        axes = np.ravel(axes)

        i = 0
        for i, (name, ax) in enumerate(zip(lin_names, axes)):
            ax.set_title(name)  # ylabel not yet supported
            sc.pl.violin(adata, keys=[name], ax=ax, **kwargs)
        for ax in axes[i + 1 :]:
            ax.remove()
        for name in to_clean:
            del adata.obs[name]

        return fig

    def plot_violin_no_cluster_key():
        key = "Endpoints" if final else "Startpoints"
        kwargs["show"] = False
        kwargs.pop("ax", None)
        kwargs.pop("keys", None)  # don't care
        kwargs.pop("save", None)
        kwargs["groupby"] = key

        data = np.ravel(np.array(adata.obsm[lk]).T)[..., np.newaxis]
        dadata = AnnData(np.zeros_like(data))
        dadata.obs["Absorption probability"] = data
        dadata.obs[key] = (
            pd.Series(np.ravel([[n] * adata.n_obs for n in adata.obsm[lk].names]))
            .astype("category")
            .values
        )
        dadata.uns[f"{key}_colors"] = adata.obsm[lk].colors

        fig, ax = plt.subplots(
            figsize=figsize if figsize is not None else (8, 6), dpi=dpi
        )
        ax.set_title("All Clusters")
        sc.pl.violin(dadata, keys=["Absorption probability"], ax=ax, **kwargs)

        return fig

    if mode not in _cluster_fates_modes:
        raise ValueError(
            f"Invalid mode: `{mode!r}`. Valid options are: `{_cluster_fates_modes}`."
        )
    if cluster_key is not None:
        if cluster_key not in adata.obs:
            raise KeyError(f"Key `{cluster_key!r}` not found in `adata.obs`.")
    elif mode not in ("bar", "violin"):
        raise ValueError(
            f"Not specifying cluster key is only available for modes `'bar'` and `'violin'`."
        )

    if cluster_key is not None:
        if clusters is not None:
            if isinstance(clusters, str):
                clusters = [clusters]
            clusters = _make_unique(clusters)
            if mode in ("paga", "paga_pie"):
                logg.debug(
                    f"DEBUG: Setting `clusters` to all available ones because of `mode={mode!r}`"
                )
                clusters = list(adata.obs[cluster_key].cat.categories)
            else:
                for cname in clusters:
                    if cname not in adata.obs[cluster_key].cat.categories:
                        raise KeyError(
                            f"Cluster `{cname!r}` not found in `adata.obs[{cluster_key!r}]`"
                        )
        else:
            clusters = list(adata.obs[cluster_key].cat.categories)
    else:
        clusters = ["All"]

    lk = str(LinKey.FORWARD if final else LinKey.BACKWARD)
    if lk not in adata.obsm:
        raise KeyError(f"Lineages key `{lk!r}` not found in `adata.obsm`.")

    if lineages is not None:
        if isinstance(lineages, str):
            lineages = [lineages]
        lineages = _make_unique(lineages)
        for ep in lineages:
            if ep not in adata.obsm[lk].names:
                raise ValueError(
                    f"Endpoint `{ep!r}` not found in `adata.obsm[{lk!r}].names`."
                )
        lin_names = list(lineages)
    else:
        # must be list for sc.pl.violin, else cats str
        lin_names = list(adata.obsm[lk].names)

    if mode == "violin" and clusters != ["All"]:
        # TODO: temporary fix, until subclassing is made ready
        names, colors = adata.obsm[lk].names, adata.obsm[lk].colors
        adata = adata[np.isin(adata.obs[cluster_key], clusters)].copy()
        adata.obsm[lk] = Lineage(adata.obsm[lk], names=names, colors=colors)

    d = odict()
    for name in clusters:
        mask = (
            np.ones((adata.n_obs,), dtype=np.bool)
            if name == "All"
            else (adata.obs[cluster_key] == name).values
        )
        mask = list(np.array(mask, dtype=np.bool))
        data = adata.obsm[lk][mask, lin_names].X
        mean = np.nanmean(data, axis=0)
        std = np.nanstd(data, axis=0) / np.sqrt(data.shape[0])
        d[name] = [mean, std]

    logg.debug(f"DEBUG: Using mode: `{mode!r}`")
    if mode == "bar":
        fig = plot_bar()
    elif mode == "paga":
        if "paga" not in adata.uns:
            raise KeyError("Compute PAGA first as `scanpy.tl.paga()`.")
        fig = plot_paga()
    elif mode == "paga_pie":
        if "paga" not in adata.uns:
            raise KeyError("Compute PAGA first as `scanpy.tl.paga()`.")
        fig = plot_paga_pie()
    elif mode == "violin":
        fig = plot_violin_no_cluster_key() if cluster_key is None else plot_violin()
    else:
        raise ValueError(
            f"Invalid mode `{mode!r}`. Valid options are: `{_cluster_fates_modes}`."
        )

    if save is not None:
        save_fig(fig, save)

    fig.show()


def similarity_plot(
    adata: AnnData,
    cluster_key: str,
    clusters: Optional[List[str]] = None,
    n_samples: int = 1000,
    cmap: mpl.colors.ListedColormap = cm.viridis,
    fontsize: float = 14,
    rotation: float = 45,
    figsize: Tuple[float, float] = (12, 10),
    final: bool = True,
    save: Optional[str] = None,
) -> None:
    """
    Compare clusters with respect to their absorption probabilities in a heatmap.

    For each cluster, we compute how likely an 'average cell' is to go to each endpoint or to come from each
    starting point. We then compare these averaged probabilities using Cramér's V statistic, see
    `here <https://en.wikipedia.org/wiki/Cram%C3%A9r%27s_V>`_. The similarity is defined as :math:`1 - Cramér's V`.

    .. image:: https://raw.githubusercontent.com/theislab/cellrank/master/resources/images/similarity_plot.png
       :width: 400px
       :align: center

    Params
    ------
    adata: :class:`anndata.AnnData`
        Annotated data object.
    cluster_key
        Key in :paramref:`adata` `.obs` corresponding the the clustering.
    clusters
        Clusters in :paramref:`adata` `.obs` to consider.
        If `None`, all cluster will be considered.
    n_samples
        Number of samples per cluster.
    cmap
        Colormap to use.
    fontsize
        Font size of the labels.
    rotation
        Rotation of labels on x-axis.
    figsize
        Size of the figure.
    final
        Whether to consider cells going to final states or vice versa.
    save
        Filename where to save the plot.
        If `None`, just shows the plot.

    Returns
    -------
    None
        Nothing, just plots the similarity matrix.
        Optionally saves the figure based on :paramref:`save`.
    """

    logg.debug("DEBUG: Getting the counts")
    data = _counts(
        adata,
        cluster_key=cluster_key,
        clusters=clusters,
        n_samples=n_samples,
        final=final,
    )

    cluster_names = list(data.keys())
    logg.debug("DEBUG: Calculating Cramer`s V statistic")
    sim = [
        [1 - _cramers_v(data[name2], data[name]) for name in cluster_names]
        for name2 in cluster_names
    ]

    # Plotting function
    fig, ax = plt.subplots(figsize=figsize)

    im = ax.imshow(sim, cmap=cmap)

    ax.set_xticks(range(len(cluster_names)))
    ax.set_yticks(range(len(cluster_names)))

    ax.set_xticklabels(cluster_names, fontsize=fontsize, rotation=rotation)
    ax.set_yticklabels(cluster_names, fontsize=fontsize)

    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)

    cbar = ax.figure.colorbar(im, ax=ax, norm=mpl.colors.Normalize(vmin=0, vmax=1))
    cbar.set_ticks(np.linspace(0, 1, 10))
    cbar.set_label("Similarity", rotation=90, va="top")

    if save is not None:
        save_fig(fig, save)

    fig.show()
