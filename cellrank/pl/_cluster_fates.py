# -*- coding: utf-8 -*-
"""Cluster fatess and similarity module."""

from math import ceil
from types import MappingProxyType
from typing import Any, Tuple, Union, Mapping, TypeVar, Optional, Sequence
from pathlib import Path
from collections import OrderedDict as odict

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors
import matplotlib.pyplot as plt

from cellrank import logging as logg
from cellrank.ul._docs import d, inject_docs
from cellrank.pl._utils import _position_legend
from cellrank.tl._utils import RandomKeys, save_fig, _unique_order_preserving
from cellrank.ul._utils import valuedispatch
from cellrank.tl._constants import ModeEnum, DirPrefix, AbsProbKey, TerminalStatesPlot

AnnData = TypeVar("AnnData")


class ClusterFatesMode(ModeEnum):  # noqa
    BAR = "bar"
    PAGA = "paga"
    PAGA_PIE = "paga_pie"
    VIOLIN = "violin"
    HEATMAP = "heatmap"
    CLUSTERMAP = "clustermap"


@d.dedent
@inject_docs(m=ClusterFatesMode)
def cluster_fates(
    adata: AnnData,
    mode: str = ClusterFatesMode.PAGA_PIE.s,
    backward: bool = False,
    lineages: Optional[Union[str, Sequence[str]]] = None,
    cluster_key: Optional[str] = "clusters",
    clusters: Optional[Union[str, Sequence[str]]] = None,
    basis: Optional[str] = None,
    cbar: bool = True,
    ncols: Optional[int] = None,
    sharey: bool = False,
    fmt: str = "0.2f",
    xrot: float = 90,
    legend_kwargs: Mapping[str, Any] = MappingProxyType({"loc": "best"}),
    figsize: Optional[Tuple[float, float]] = None,
    dpi: Optional[int] = None,
    save: Optional[Union[str, Path]] = None,
    **kwargs,
) -> None:
    """
    Plot aggregate lineage probabilities at a cluster level.

    This can be used to investigate how likely a certain cluster is to go to the %(terminal)s states,or in turn to have
    descended from the %(initial)s states. For mode `{m.PAGA.s!r}` and `{m.PAGA_PIE.s!r}`, we use *PAGA*, see [Wolf19]_.

    Parameters
    ----------
    %(adata)s
    mode
        Type of plot to show. Valid options are:

            - `{m.BAR.s!r}` - barplot, one panel per cluster.
            - `{m.PAGA.s!r}` - scanpy's PAGA, one per %(initial_or_terminal)s state, colored in by fate.
            - `{m.PAGA_PIE.s!r}` - scanpy's PAGA with pie charts indicating aggregated fates.
            - `{m.VIOLIN.s!r}` - violin plots, one per %(initial_or_terminal)s state.
            - `{m.HEATMAP.s!r}` - a heatmap, showing average fates per cluster.
            - `{m.CLUSTERMAP.s!r}` - same as a heatmap, but with a dendrogram.
    %(backward)s
    lineages
        Lineages for which to visualize absorption probabilities. If `None`, use all lineages.
    cluster_key
        Key in ``adata.obs`` containing the clusters.
    clusters
        Clusters to visualize. If `None`, all clusters will be plotted.
    basis
        Basis for scatterplot to use when ``mode={m.PAGA_PIE.s!r}``. If `None`, don't show the scatterplot.
    cbar
        Whether to show colorbar when ``mode={m.PAGA_PIE.s!r}``.
    ncols
        Number of columns when ``mode={m.BAR.s!r}`` or ``mode={m.PAGA.s!r}``.
    sharey
        Whether to share y-axis when ``mode={m.BAR.s!r}``.
    fmt
        Format when using ``mode={m.HEATMAP.s!r}`` or ``mode={m.CLUSTERMAP.s!r}``.
    xrot
        Rotation of the labels on the x-axis.
    figsize
        Size of the figure.
    legend_kwargs
        Keyword arguments for :func:`matplotlib.axes.Axes.legend`, such as `'loc'` for legend position.
        For ``mode={m.PAGA_PIE.s!r}`` and ``basis='...'``, this controls the placement of the
        absorption probabilities legend.
    %(plotting)s
    **kwargs
        Keyword arguments for :func:`scvelo.pl.paga`, :func:`scanpy.pl.violin` or :func:`matplotlib.pyplot.bar`,
        depending on the value of ``mode``.

    Returns
    -------
    %(just_plots)s
    """

    from scanpy.plotting import violin
    from scvelo.plotting import paga

    from seaborn import heatmap, clustermap

    @valuedispatch
    def plot(mode: ClusterFatesMode, *_args, **_kwargs):
        raise NotImplementedError(mode.value)

    @plot.register(ClusterFatesMode.BAR)
    def _():
        cols = 4 if ncols is None else ncols
        n_rows = ceil(len(clusters) / cols)
        fig = plt.figure(
            None, (3.5 * cols, 5 * n_rows) if figsize is None else figsize, dpi=dpi
        )
        fig.tight_layout()

        gs = plt.GridSpec(n_rows, cols, figure=fig, wspace=0.5, hspace=0.5)

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
            current_ax.set_xticklabels(lin_names, rotation=xrot)
            if not is_all:
                current_ax.set_xlabel(points)
            current_ax.set_ylabel("absorption probability")
            current_ax.set_title(k)

        return fig

    @plot.register(ClusterFatesMode.PAGA)
    def _():
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
        # fig.tight_layout()  can't use this because colorbar.make_axes fails

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

            vmin = np.min(colors + [vmin])
            vmax = np.max(colors + [vmax])

            paga(adata, **kwargs)

        if cbar:
            norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
            cax, _ = mpl.colorbar.make_axes(ax, aspect=60)
            _ = mpl.colorbar.ColorbarBase(
                cax,
                ticks=np.linspace(norm.vmin, norm.vmax, 5),
                norm=norm,
                cmap=kwargs["cmap"],
                label="average absorption probability",
            )

        for ax in axes[i + 1 :]:  # noqa
            ax.remove()

        return fig

    @plot.register(ClusterFatesMode.PAGA_PIE)
    def _():
        colors = list(adata.obsm[lk][:, lin_names].colors)
        colors = {i: odict(zip(colors, mean)) for i, (mean, _) in enumerate(d.values())}

        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        fig.tight_layout()

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

        ax = paga(adata, **kwargs)
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

    @plot.register(ClusterFatesMode.VIOLIN)
    def _():
        kwargs.pop("ax", None)
        kwargs.pop("keys", None)
        kwargs.pop("save", None)  # we will handle saving

        kwargs["show"] = False
        kwargs["groupby"] = cluster_key
        kwargs["rotation"] = xrot

        cols = len(lin_names) if ncols is None else ncols
        nrows = ceil(len(lin_names) / cols)

        fig, axes = plt.subplots(
            nrows,
            cols,
            figsize=(6 * cols, 4 * nrows) if figsize is None else figsize,
            sharey=sharey,
            dpi=dpi,
        )
        fig.tight_layout()
        fig.subplots_adjust(wspace=0.2, hspace=0.3)

        if not isinstance(axes, np.ndarray):
            axes = [axes]
        axes = np.ravel(axes)

        with RandomKeys(adata, len(lin_names), where="obs") as keys:
            _i = 0
            for _i, (name, key, ax) in enumerate(zip(lin_names, keys, axes)):
                adata.obs[key] = adata.obsm[lk][name].X
                ax.set_title(f"{dir_prefix} {name}")
                violin(
                    adata, ylabel="absorption probability", keys=key, ax=ax, **kwargs
                )
            for ax in axes[_i + 1 :]:  # noqa
                ax.remove()

        return fig

    def plot_violin_no_cluster_key():
        from anndata import AnnData as _AnnData

        kwargs.pop("ax", None)
        kwargs.pop("keys", None)  # don't care
        kwargs.pop("save", None)

        kwargs["show"] = False
        kwargs["groupby"] = points
        kwargs["xlabel"] = None
        kwargs["rotation"] = xrot

        data = np.ravel(adata.obsm[lk].X.T)[..., np.newaxis]
        tmp = _AnnData(csr_matrix(data.shape, dtype=np.float32))
        tmp.obs["absorption probability"] = data
        tmp.obs[points] = (
            pd.Series(
                np.concatenate(
                    [
                        [f"{dir_prefix.lower()} {n}"] * adata.n_obs
                        for n in adata.obsm[lk].names
                    ]
                )
            )
            .astype("category")
            .values
        )
        tmp.obs[points].cat.reorder_categories(
            [f"{dir_prefix.lower()} {n}" for n in adata.obsm[lk].names], inplace=True
        )
        tmp.uns[f"{points}_colors"] = adata.obsm[lk].colors

        fig, ax = plt.subplots(
            figsize=figsize if figsize is not None else (8, 6), dpi=dpi
        )
        ax.set_title(points.capitalize())

        violin(tmp, keys=["absorption probability"], ax=ax, **kwargs)

        return fig

    @plot.register(ClusterFatesMode.HEATMAP)
    def _():
        data = pd.DataFrame(
            [mean for mean, _ in d.values()], columns=lin_names, index=clusters
        ).T

        title = kwargs.pop("title", "average fate per cluster")
        vmin, vmax = data.values.min(), data.values.max()
        cbar_kws = {
            "label": "probability",
            "ticks": np.linspace(vmin, vmax, 5),
            "format": "%.3f",
        }
        kwargs.setdefault("cmap", "viridis")

        if use_clustermap:
            kwargs["cbar_pos"] = (0, 0.9, 0.025, 0.15) if cbar else None
            max_size = float(max(data.shape))

            g = clustermap(
                data,
                annot=True,
                vmin=vmin,
                vmax=vmax,
                fmt=fmt,
                row_colors=adata.obsm[lk][lin_names].colors,
                dendrogram_ratio=(
                    0.15 * data.shape[0] / max_size,
                    0.15 * data.shape[1] / max_size,
                ),
                cbar_kws=cbar_kws,
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
                vmin=vmin,
                vmax=vmax,
                annot=True,
                fmt=fmt,
                cbar=cbar,
                cbar_kws=cbar_kws,
                ax=ax,
                **kwargs,
            )
            ax.set_title(title)
            ax.set_xlabel(cluster_key)
            ax.set_ylabel("lineage")

        g.set_xticklabels(g.get_xticklabels(), rotation=xrot)
        g.set_yticklabels(g.get_yticklabels(), rotation=0)

        return fig

    mode = ClusterFatesMode(mode)

    if cluster_key is not None:
        if cluster_key not in adata.obs:
            raise KeyError(f"Key `{cluster_key!r}` not found in `adata.obs`.")
    elif mode not in (mode.BAR, mode.VIOLIN):
        raise ValueError(
            f"Not specifying cluster key is only available for modes "
            f"`{ClusterFatesMode.BAR!r}` and `{ClusterFatesMode.VIOLIN!r}`, found `mode={mode!r}`."
        )

    if backward:
        lk = AbsProbKey.BACKWARD.s
        points = TerminalStatesPlot.BACKWARD.s
        dir_prefix = DirPrefix.BACKWARD.s
    else:
        lk = AbsProbKey.FORWARD.s
        points = TerminalStatesPlot.FORWARD.s
        dir_prefix = DirPrefix.FORWARD.s

    if cluster_key is not None:
        is_all = False
        if clusters is not None:
            if isinstance(clusters, str):
                clusters = [clusters]
            clusters = _unique_order_preserving(clusters)
            if mode in (mode.PAGA, mode.PAGA_PIE):
                logg.debug(
                    f"Setting `clusters` to all available ones because of `mode={mode!r}`"
                )
                clusters = list(adata.obs[cluster_key].cat.categories)
            else:
                for cname in clusters:
                    if cname not in adata.obs[cluster_key].cat.categories:
                        raise KeyError(
                            f"Cluster `{cname!r}` not found in `adata.obs[{cluster_key!r}]`."
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
        lin_names = _unique_order_preserving(lineages)
    else:
        # must be list for `sc.pl.violin`, else cats str
        lin_names = list(adata.obsm[lk].names)
    _ = adata.obsm[lk][lin_names]

    if mode == mode.VIOLIN and not is_all:
        adata = adata[np.isin(adata.obs[cluster_key], clusters)].copy()

    d = odict()
    for name in clusters:
        mask = (
            np.ones((adata.n_obs,), dtype=np.bool)
            if is_all
            else (adata.obs[cluster_key] == name).values
        )
        mask = np.array(mask, dtype=np.bool)
        data = adata.obsm[lk][mask, lin_names].X
        mean = np.nanmean(data, axis=0)
        std = np.nanstd(data, axis=0) / np.sqrt(data.shape[0])
        d[name] = [mean, std]

    logg.debug(f"Plotting in mode `{mode!r}`")
    use_clustermap = False
    if mode == mode.CLUSTERMAP:
        use_clustermap = True
        mode = mode.HEATMAP
    elif (
        mode in (ClusterFatesMode.PAGA, ClusterFatesMode.PAGA_PIE)
        and "paga" not in adata.uns
    ):
        raise KeyError("Compute PAGA first as `scvelo.tl.paga()`.")

    fig = (
        plot_violin_no_cluster_key()
        if mode == ClusterFatesMode.VIOLIN and cluster_key is None
        else plot(mode)
    )

    if save is not None:
        save_fig(fig, save)

    fig.show()
