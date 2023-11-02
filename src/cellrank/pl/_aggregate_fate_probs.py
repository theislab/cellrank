import collections
import enum
import math
import pathlib
import types
from typing import Any, Literal, Mapping, Optional, Sequence, Tuple, Union

from scvelo.plotting import paga

import numpy as np
import pandas as pd
import scipy.sparse as sp

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import cm, colors

from anndata import AnnData
from scanpy.plotting import violin

from cellrank import logging as logg
from cellrank._utils import Lineage
from cellrank._utils._docs import d, inject_docs
from cellrank._utils._enum import ModeEnum
from cellrank._utils._key import Key
from cellrank._utils._utils import (
    RandomKeys,
    _unique_order_preserving,
    save_fig,
    valuedispatch,
)
from cellrank.pl._utils import _position_legend

__all__ = ["aggregate_fate_probabilities"]


class AggregationMode(ModeEnum):
    BAR = enum.auto()
    PAGA = enum.auto()
    PAGA_PIE = enum.auto()
    VIOLIN = enum.auto()
    HEATMAP = enum.auto()
    CLUSTERMAP = enum.auto()


@d.dedent
@inject_docs(m=AggregationMode)
def aggregate_fate_probabilities(
    adata: AnnData,
    mode: Literal["bar", "paga", "paga_pie", "violin", "heatmap", "clustermap"] = AggregationMode.PAGA_PIE,
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
    legend_kwargs: Mapping[str, Any] = types.MappingProxyType({"loc": "best"}),
    figsize: Optional[Tuple[float, float]] = None,
    dpi: Optional[int] = None,
    save: Optional[Union[str, pathlib.Path]] = None,
    **kwargs: Any,
) -> None:
    """Plot aggregate lineage probabilities at a cluster level.

    This can be used to investigate how likely a certain cluster is to go to the %(terminal)s states,
    or in turn to have descended from the %(initial)s states.
    For mode `{m.PAGA!r}` and `{m.PAGA_PIE!r}`, we use *PAGA* :cite:`wolf:19`.

    Parameters
    ----------
    %(adata)s
    mode
        Type of plot to show. Valid options are:

        - ``{m.BAR!r}`` - barplot, one panel per cluster. The whiskers correspond to the standard error of the mean.
        - ``{m.PAGA!r}`` - :func:`~scvelo.pl.paga`, one per %(initial_or_terminal)s state, colored in by fate.
        - ``{m.PAGA_PIE!r}`` - :func:`~scvelo.pl.paga` with pie charts indicating aggregated fates.
        - ``{m.VIOLIN!r}`` - violin plots, one per %(initial_or_terminal)s state.
        - ``{m.HEATMAP!r}`` - a heatmap, showing average fates per cluster.
        - ``{m.CLUSTERMAP!r}`` - same as a heatmap, but with a dendrogram.
    %(backward)s
    lineages
        Lineages for which to visualize the fate probabilities. If :obj:`None`, use all lineages.
    cluster_key
        Key in :attr:`~anndata.AnnData.obs` containing the clusters.
    clusters
        Clusters to visualize. If :obj:`None`, all clusters will be plotted.
    basis
        Basis for scatterplot to use when ``mode = {m.PAGA_PIE!r}``. If :obj:`None`, don't show the scatterplot.
    cbar
        Whether to show colorbar when ``mode = {m.PAGA_PIE!r}``.
    ncols
        Number of columns when ``mode = {m.BAR!r}`` or ``mode = {m.PAGA!r}``.
    sharey
        Whether to share y-axis when ``mode = {m.BAR!r}``.
    fmt
        Format when using ``mode = {m.HEATMAP!r}`` or ``mode = {m.CLUSTERMAP!r}``.
    xrot
        Rotation of the labels on the x-axis.
    figsize
        Size of the figure.
    legend_kwargs
        Keyword arguments for :meth:`~matplotlib.axes.Axes.legend`.
        For ``mode = {m.PAGA_PIE!r}`` and ``basis = '...'``, this controls the placement of the
        fate probabilities legend.
    %(plotting)s
    kwargs
        Keyword arguments for :func:`~scvelo.pl.paga`, :func:`~scanpy.pl.violin` or
        :func:`~matplotlib.pyplot.bar`, depending on the ``mode``.

    Returns
    -------
    %(just_plots)s
    """

    @valuedispatch
    def plot(mode: AggregationMode, *_args, **_kwargs):
        raise NotImplementedError(mode.value)

    @plot.register(AggregationMode.BAR)
    def _():
        cols = 4 if ncols is None else ncols
        n_rows = math.ceil(len(clusters) / cols)
        fig = plt.figure(None, (3.5 * cols, 5 * n_rows) if figsize is None else figsize, dpi=dpi)
        fig.tight_layout()

        gs = plt.GridSpec(n_rows, cols, figure=fig, wspace=0.5, hspace=0.5)

        ax = None
        colors = list(probs.colors)

        for g, k in zip(gs, d.keys()):
            current_ax = fig.add_subplot(g, sharey=ax)
            current_ax.bar(
                x=np.arange(probs.nlin),
                height=d[k][0],
                color=colors,
                yerr=d[k][1],
                ecolor="black",
                capsize=10,
                **kwargs,
            )
            if sharey:
                ax = current_ax

            current_ax.set_xticks(np.arange(probs.nlin))
            current_ax.set_xticklabels(probs.names, rotation=xrot)
            if not is_all:
                current_ax.set_xlabel(term_states)
            current_ax.set_ylabel("fate probability")
            current_ax.set_title(k)

        return fig

    @plot.register(AggregationMode.PAGA)
    def _():
        kwargs["save"] = None
        kwargs["show"] = False
        if "cmap" not in kwargs:
            kwargs["cmap"] = cm.viridis

        cols = probs.nlin if ncols is None else ncols
        nrows = math.ceil(probs.nlin / cols)
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

        for i, (ax, lineage_name) in enumerate(zip(axes, probs.names)):
            cols = [v[0][i] for v in d.values()]
            kwargs["ax"] = ax
            kwargs["colors"] = tuple(cols)
            kwargs["title"] = f"{direction} {lineage_name}"

            vmin = np.min(cols + [vmin])
            vmax = np.max(cols + [vmax])

            paga(adata, **kwargs)

        if cbar:
            norm = colors.Normalize(vmin=vmin, vmax=vmax)
            cax, _ = mpl.colorbar.make_axes(ax, aspect=60)
            _ = mpl.colorbar.ColorbarBase(
                cax,
                ticks=np.linspace(norm.vmin, norm.vmax, 5),
                norm=norm,
                cmap=kwargs["cmap"],
                label="average fate probability",
            )

        for ax in axes[i + 1 :]:
            ax.remove()

        return fig

    @plot.register(AggregationMode.PAGA_PIE)
    def _():
        colors = {i: collections.OrderedDict(zip(probs.colors, mean)) for i, (mean, _) in enumerate(d.values())}

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
            for lineage_name, color in zip(probs.names, colors[0].keys()):
                handles += [ax.scatter([], [], label=lineage_name, c=color)]
            if len(colors[0].keys()) != Lineage.from_adata(adata, backward=backward).nlin:
                handles += [ax.scatter([], [], label="Rest", c="grey")]

            second_legend = _position_legend(
                ax,
                legend_loc=legend_kwargs["loc"],
                handles=handles,
                **{k: v for k, v in legend_kwargs.items() if k != "loc"},
                title=term_states,
            )
            fig.add_artist(second_legend)

        return fig

    @plot.register(AggregationMode.VIOLIN)
    def _():
        kwargs.pop("ax", None)
        kwargs.pop("keys", None)
        kwargs.pop("save", None)  # we will handle saving

        kwargs["show"] = False
        kwargs["groupby"] = cluster_key
        kwargs["rotation"] = xrot

        cols = probs.nlin if ncols is None else ncols
        nrows = math.ceil(probs.nlin / cols)

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

        with RandomKeys(adata, probs.nlin, where="obs") as keys:
            _i = 0
            for _i, (name, key, ax) in enumerate(zip(probs.names, keys, axes)):
                adata.obs[key] = probs[name].X[:, 0]
                ax.set_title(f"{direction} {name}")
                violin(adata, ylabel="fate probability", keys=key, ax=ax, **kwargs)
            for ax in axes[_i + 1 :]:
                ax.remove()

        return fig

    def plot_violin_no_cluster_key():
        kwargs.pop("ax", None)
        kwargs.pop("keys", None)  # don't care
        kwargs.pop("save", None)

        kwargs["show"] = False
        kwargs["groupby"] = term_states
        kwargs["xlabel"] = None
        kwargs["rotation"] = xrot

        data = np.ravel(probs.X.T)[..., None]
        tmp = AnnData(sp.csr_matrix(data.shape, dtype=data.dtype))
        tmp.obs["fate probability"] = data
        tmp.obs[term_states] = (
            pd.Series(np.concatenate([[f"{direction.lower()} {n}"] * adata.n_obs for n in probs.names]))
            .astype("category")
            .values
        )
        tmp.obs[term_states] = tmp.obs[term_states].cat.reorder_categories(
            [f"{direction.lower()} {n}" for n in probs.names]
        )
        tmp.uns[f"{term_states}_colors"] = probs.colors

        fig, ax = plt.subplots(figsize=figsize if figsize is not None else (8, 6), dpi=dpi)
        ax.set_title(term_states.capitalize())

        violin(tmp, keys=["fate probability"], ax=ax, **kwargs)

        return fig

    @plot.register(AggregationMode.HEATMAP)
    def _():
        data = pd.DataFrame([mean for mean, _ in d.values()], columns=probs.names, index=clusters).T

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

            g = sns.clustermap(
                data=data,
                annot=True,
                vmin=vmin,
                vmax=vmax,
                fmt=fmt,
                row_colors=probs.colors,
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
            g = sns.heatmap(
                data=data,
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

    mode = AggregationMode(mode)

    if cluster_key is not None:
        if cluster_key not in adata.obs:
            raise KeyError(f"Key `{cluster_key!r}` not found in `adata.obs`.")
    elif mode not in (mode.BAR, mode.VIOLIN):
        raise ValueError(
            f"Not specifying cluster key is only available for modes "
            f"`{AggregationMode.BAR!r}` and `{AggregationMode.VIOLIN!r}`, found `mode={mode!r}`."
        )

    term_states = Key.obs.term_states(backward)
    direction = Key.where(backward)
    if cluster_key is not None:
        is_all = False
        if clusters is not None:
            if isinstance(clusters, str):
                clusters = [clusters]
            clusters = _unique_order_preserving(clusters)
            if mode in (mode.PAGA, mode.PAGA_PIE):
                logg.debug(f"Setting `clusters` to all available ones because of `mode={mode!r}`")
                clusters = list(adata.obs[cluster_key].cat.categories)
            else:
                for cname in clusters:
                    if cname not in adata.obs[cluster_key].cat.categories:
                        raise KeyError(f"Cluster `{cname!r}` not found in `adata.obs[{cluster_key!r}]`.")
        else:
            clusters = list(adata.obs[cluster_key].cat.categories)
    else:
        is_all = True
        clusters = [term_states]

    probs = Lineage.from_adata(adata, backward=backward)
    if lineages is None:
        # must be list for `sc.pl.violin`, else cats str
        lineages = list(probs.names)
    elif isinstance(lineages, str):
        lineages = [lineages]
    lineages = _unique_order_preserving(lineages)
    probs = probs[:, lineages]

    if mode == mode.VIOLIN and not is_all:
        mask = np.isin(adata.obs[cluster_key], clusters)
        adata = adata[mask].copy()
        probs = probs[mask]

    d = collections.OrderedDict()
    for name in clusters:
        data = probs.X if is_all else probs.X[(adata.obs[cluster_key] == name).to_numpy()]
        mean = np.nanmean(data, axis=0)
        std = np.nanstd(data, axis=0) / np.sqrt(data.shape[0])
        d[name] = [mean, std]

    logg.debug(f"Plotting in mode `{mode!r}`")
    use_clustermap = False
    if mode == mode.CLUSTERMAP:
        use_clustermap = True
        mode = mode.HEATMAP
    elif mode in (AggregationMode.PAGA, AggregationMode.PAGA_PIE) and "paga" not in adata.uns:
        raise KeyError("Compute PAGA first as `scvelo.tl.paga()` or `scanpy.tl.paga()`.")

    fig = plot_violin_no_cluster_key() if mode == AggregationMode.VIOLIN and cluster_key is None else plot(mode)

    if save is not None:
        save_fig(fig, save)
