# -*- coding: utf-8 -*-
"""Cluster fatess and similarity module."""

from math import ceil
from types import MappingProxyType
from typing import Any, List, Tuple, Union, Mapping, Optional, Sequence
from pathlib import Path
from collections import OrderedDict as odict

import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors
import matplotlib.pyplot as plt
from seaborn import heatmap, clustermap

import scanpy as sc
import scvelo as scv
from scanpy import logging as logg
from anndata import AnnData

from cellrank.tools._utils import save_fig
from cellrank.utils._utils import _make_unique
from cellrank.plotting._utils import _position_legend
from cellrank.tools._constants import LinKey
from cellrank.tools._exact_mc_test import _counts, _cramers_v

_cluster_fates_modes = ("bar", "paga", "paga_pie", "violin", "heatmap", "clustermap")


def cluster_fates(
    adata: AnnData,
    cluster_key: Optional[str] = "louvain",
    lineage_key: Optional[str] = None,
    clusters: Optional[Union[str, Sequence[str]]] = None,
    lineages: Optional[Union[str, Sequence[str]]] = None,
    mode: str = "bar",
    final: bool = True,
    basis: Optional[str] = None,
    show_cbar: bool = True,
    ncols: Optional[int] = None,
    sharey: bool = False,
    save: Optional[Union[str, Path]] = None,
    legend_kwargs: Mapping[str, Any] = MappingProxyType({"loc": "best"}),
    figsize: Optional[Tuple[float, float]] = None,
    dpi: Optional[int] = None,
    **kwargs,
) -> None:
    """
    Plot aggregate lineage probabilities at a cluster level.

    This can be used to investigate how likely a certain cluster is to go to the final states, or in turn to have
    descended from the root states. For mode `'paga'` and `'paga_pie'`, we use *PAGA*, see [Wolf19]_.

    .. image:: https://raw.githubusercontent.com/theislab/cellrank/master/resources/images/cluster_fates.png
       :width: 400px
       :align: center

    Params
    ------
    adata : :class:`anndata.AnnData`
        Annotated data object.
    cluster_key
        Key in :paramref:`adata` `.obs` containing the clusters.
    lineage_key
        Key in :paramref:`adata` `.obsm` containing fate probabilities.
    clusters
        Clusters to visualize.
        If `None`, all clusters will be plotted.
    lineages
        Lineages for which to visualize absorption probabilities.
        If `None`, use all available lineages.
    mode
        Type of plot to show.

        - `'bar'`: barplot, one panel per cluster
        - `'paga'`: scanpy's PAGA, one per root/final state, colored in by fate
        - `'paga_pie'`: scanpy's PAGA with pie charts indicating aggregated fates
        - `'violin'`: violin plots, one per root/final state
        - `'heatmap'`: seaborn heatmap, showing average fates per cluster
        - `'clustermap'`: same as heatmap, but with dendrogram
    dpi
        Dots per inch.
    final
        Whether to consider cells going to final states or vice versa.
    basis
        Basis for scatterplot to use when :paramref:`mode` `='paga_pie'`. If `None`, don't show the scatterplot.
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
        For `mode='paga_pie'` and `basis='...'`, this controls the placement of the absorption probabilities legend.
    figsize
        Size of the figure. If `None`, it will be set automatically.
    dpi
        Dots per inch.
    kwargs
        Keyword arguments for :func:`scvelo.pl.paga`, :func:`scanpy.pl.violin` or :func:`matplotlib.pyplot.bar`,
        depending on :paramref:`mode`.

    Returns
    -------
    None
        Nothing, just plots the fates for specified :paramref:`clusters` and :paramref:`lineages`.
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
            current_ax.set_xticklabels(
                lin_names, rotation=xrot if has_xrot else "vertical"
            )
            if not is_all:
                current_ax.set_xlabel(points)
            current_ax.set_ylabel("probability")
            current_ax.set_title(k)

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
            figsize=(7 * cols, 4 * nrows) if figsize is None else figsize,
            constrained_layout=True,
            dpi=dpi,
        )

        i = 0
        axes = [axes] if not isinstance(axes, np.ndarray) else np.ravel(axes)
        vmin, vmax = np.inf, -np.inf

        if basis is not None:
            kwargs["basis"] = basis
            kwargs["scatter_flag"] = True
            kwargs["color"] = cluster_key

        for i, (ax, lineage_name) in enumerate(zip(axes, lin_names)):
            colors = [v[0][i] for v in d.values()]
            kwargs["ax"] = ax
            kwargs["colors"] = tuple(colors)
            kwargs["title"] = f"{dir_prefix} {lineage_name}"

            scv.pl.paga(adata, **kwargs)

        if show_cbar:
            norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
            cax, _ = mpl.colorbar.make_axes(ax, aspect=100)  # new matplotlib feature
            _ = mpl.colorbar.ColorbarBase(
                cax, norm=norm, cmap=kwargs["cmap"], label="probability"
            )

        for ax in axes[i + 1 :]:  # noqa
            ax.remove()

        return fig

    def plot_paga_pie():
        colors = list(adata.obsm[lk][:, lin_names].colors)
        colors = {i: odict(zip(colors, mean)) for i, (mean, _) in enumerate(d.values())}

        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

        kwargs["ax"] = ax
        kwargs["show"] = False
        kwargs["colorbar"] = False  # has to be disabled
        kwargs["show"] = False

        kwargs["node_colors"] = colors
        kwargs.pop("save", None)  # we will handle saving

        kwargs["transitions"] = kwargs.get("transitions", "transitions_confidence")
        if "legend_loc" in kwargs:
            orig_ll = kwargs["legend_loc"]
            if orig_ll != "on data":
                kwargs["legend_loc"] = "none"  # we will handle legend
        else:
            orig_ll = None
            kwargs["legend_loc"] = "on data"

        if basis is not None:
            kwargs["basis"] = basis
            kwargs["scatter_flag"] = True
            kwargs["color"] = cluster_key

        ax = scv.pl.paga(adata, **kwargs)
        ax.set_title(kwargs.get("title", cluster_key))

        if basis is not None and orig_ll not in ("none", "on data", None):
            handles = []
            for cluster_name, color in zip(
                adata.obs[f"{cluster_key}"].cat.categories,
                adata.uns[f"{cluster_key}_colors"],
            ):
                handles += [ax.scatter([], [], label=cluster_name, c=color)]
            first_legend = _position_legend(
                ax,
                legend_loc=orig_ll,
                handles=handles,
                **{k: v for k, v in legend_kwargs.items() if k != "loc"},
                title=cluster_key,
            )
            fig.add_artist(first_legend)

        if legend_kwargs.get("loc", None) not in ("none", "on data", None):
            # we need to use these, because scvelo can have its own handles and
            # they would be plotted here
            handles = []
            for lineage_name, color in zip(lin_names, colors[0].keys()):
                handles += [ax.scatter([], [], label=lineage_name, c=color)]
            if len(colors[0].keys()) != len(adata.obsm[lk].names):
                handles += [ax.scatter([], [], label="Rest", c="grey")]

            second_legend = _position_legend(
                ax,
                legend_loc=legend_kwargs["loc"],
                handles=handles,
                **{k: v for k, v in legend_kwargs.items() if k != "loc"},
                title=points,
            )
            fig.add_artist(second_legend)

        return fig

    def plot_violin():
        kwargs.pop("ax", None)
        kwargs.pop("keys", None)
        kwargs.pop("save", None)  # we will handle saving

        kwargs["show"] = False
        kwargs["groupby"] = cluster_key
        kwargs["rotation"] = xrot

        data = adata.obsm[lk]
        to_clean = []

        for name in lin_names:
            # TODO: once ylabel is implemented, the prefix isn't necessary
            key = f"{dir_prefix} {name}"
            if key not in adata.obs_keys():
                to_clean.append(key)
                adata.obs[key] = np.array(
                    data[:, name]
                )  # TODO: better approach - dummy adata

        cols = len(lin_names) if ncols is None else ncols
        nrows = ceil(len(lin_names) / cols)
        fig, axes = plt.subplots(
            nrows,
            cols,
            figsize=(6 * cols, 4 * nrows) if figsize is None else figsize,
            sharey=sharey,
            dpi=dpi,
        )
        if not isinstance(axes, np.ndarray):
            axes = [axes]
        axes = np.ravel(axes)

        i = 0
        for i, (name, ax) in enumerate(zip(lin_names, axes)):
            key = f"{dir_prefix} {name}"
            ax.set_title(key)
            sc.pl.violin(
                adata, ylabel="" if i else "probability", keys=key, ax=ax, **kwargs
            )
        for ax in axes[i + 1 :]:  # noqa
            ax.remove()
        for name in to_clean:
            del adata.obs[name]

        return fig

    def plot_violin_no_cluster_key():
        kwargs.pop("ax", None)
        kwargs.pop("keys", None)  # don't care
        kwargs.pop("save", None)

        kwargs["show"] = False
        kwargs["groupby"] = points
        kwargs["xlabel"] = None
        kwargs["rotation"] = xrot

        data = np.ravel(np.array(adata.obsm[lk]).T)[..., np.newaxis]
        dadata = AnnData(np.zeros_like(data))
        dadata.obs["probability"] = data
        dadata.obs[points] = (
            pd.Series(
                np.ravel(
                    [
                        [f"{dir_prefix.lower()} {n}"] * adata.n_obs
                        for n in adata.obsm[lk].names
                    ]
                )
            )
            .astype("category")
            .values
        )
        dadata.uns[f"{points}_colors"] = adata.obsm[lk].colors

        fig, ax = plt.subplots(
            figsize=figsize if figsize is not None else (8, 6), dpi=dpi
        )
        ax.set_title(points.capitalize())
        sc.pl.violin(dadata, keys=["probability"], ax=ax, **kwargs)

        return fig

    def plot_heatmap():
        title = kwargs.pop("title", None)
        if not title:
            title = "average fate per cluster"
        data = pd.DataFrame(
            [mean for mean, _ in d.values()], columns=lin_names, index=clusters
        ).T

        if "cmap" not in kwargs:
            kwargs["cmap"] = "viridis"

        if use_clustermap:
            kwargs["cbar_pos"] = (0, 0.9, 0.025, 0.15) if show_cbar else None
            max_size = float(max(data.shape))

            g = clustermap(
                data,
                robust=True,
                annot=True,
                fmt=".2f",
                row_colors=adata.obsm[lk][lin_names].colors,
                dendrogram_ratio=(
                    0.15 * data.shape[0] / max_size,
                    0.15 * data.shape[1] / max_size,
                ),
                figsize=figsize,
                **kwargs,
            )
            g.ax_heatmap.set_xlabel(cluster_key)
            g.ax_heatmap.set_ylabel("lineage")
            g.ax_col_dendrogram.set_title(title)

            fig = g.fig
            g = g.ax_heatmap
        else:
            fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
            g = heatmap(
                data,
                robust=True,
                annot=True,
                fmt=".2f",
                cbar=show_cbar,
                ax=ax,
                **kwargs,
            )
            ax.set_title(title)
            ax.set_xlabel(cluster_key)
            ax.set_ylabel("lineage")

        g.set_xticklabels(g.get_xticklabels(), rotation=xrot)
        g.set_yticklabels(g.get_yticklabels(), rotation=0)

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
            f"Not specifying cluster key is only available for modes `'bar'` and `'violin'`, found `mode={mode!r}`."
        )
    if lineage_key is None:
        lk = str(LinKey.FORWARD if final else LinKey.BACKWARD)
    else:
        lk = lineage_key
    points = "final states" if final else "root states"
    dir_prefix = "To" if final else "From"

    if cluster_key is not None:
        is_all = False
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
        is_all = True
        clusters = [points]

    if lk not in adata.obsm:
        raise KeyError(f"Lineage key `{lk!r}` not found in `adata.obsm`.")

    if lineages is not None:
        if isinstance(lineages, str):
            lineages = [lineages]
        lineages = _make_unique(lineages)
        for ep in lineages:
            if ep not in adata.obsm[lk].names:
                raise ValueError(
                    f"State `{ep!r}` not found in `adata.obsm[{lk!r}].names`."
                )
        lin_names = list(lineages)
    else:
        # must be list for sc.pl.violin, else cats str
        lin_names = list(adata.obsm[lk].names)

    if mode == "violin" and not is_all:
        adata = adata[np.isin(adata.obs[cluster_key], clusters)].copy()

    d = odict()
    for name in clusters:
        mask = (
            np.ones((adata.n_obs,), dtype=np.bool)
            if is_all
            else (adata.obs[cluster_key] == name).values
        )
        mask = list(np.array(mask, dtype=np.bool))
        data = adata.obsm[lk][mask, lin_names].X
        mean = np.nanmean(data, axis=0)
        std = np.nanstd(data, axis=0) / np.sqrt(data.shape[0])
        d[name] = [mean, std]

    has_xrot = "xticks_rotation" in kwargs
    xrot = kwargs.pop("xticks_rotation", 45)

    logg.debug(f"DEBUG: Using mode: `{mode!r}`")

    use_clustermap = mode == "clustermap"
    if use_clustermap:
        mode = "heatmap"

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
    elif mode == "heatmap":
        fig = plot_heatmap()
    else:
        raise ValueError(
            f"Invalid mode `{mode!r}`. Valid options are: `{_cluster_fates_modes}`."
        )

    if save is not None:
        save_fig(fig, save)

    fig.show()


def similarity_plot(
    adata: AnnData,
    cluster_key: str = "clusters",
    clusters: Optional[List[str]] = None,
    n_samples: int = 1000,
    cmap: mpl.colors.ListedColormap = cm.viridis,
    fontsize: float = 14,
    rotation: float = 45,
    title: Optional[str] = "similarity",
    figsize: Tuple[float, float] = (12, 10),
    dpi: Optional[int] = None,
    final: bool = True,
    save: Optional[Union[str, Path]] = None,
) -> None:
    """
    Compare clusters with respect to their root/final probabilities.

    For each cluster, we compute how likely an 'average cell' is to go towards to final states/come from the root
    states. We then compare these averaged probabilities using Cramér's V statistic, see
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
    title
        Title of the figure.
    dpi
        Dots per inch.
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
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    im = ax.imshow(sim, cmap=cmap)

    ax.set_xticks(range(len(cluster_names)))
    ax.set_yticks(range(len(cluster_names)))

    ax.set_xticklabels(cluster_names, fontsize=fontsize, rotation=rotation)
    ax.set_yticklabels(cluster_names, fontsize=fontsize)

    ax.set_title(title)
    ax.tick_params(top=False, bottom=True, labeltop=False, labelbottom=True)

    cbar = ax.figure.colorbar(im, ax=ax, norm=mpl.colors.Normalize(vmin=0, vmax=1))
    cbar.set_ticks(np.linspace(0, 1, 10))

    if save is not None:
        save_fig(fig, save)

    fig.show()
