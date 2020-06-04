# -*- coding: utf-8 -*-
"""Heatmap module."""

from typing import Tuple, Union, Optional, Sequence
from pathlib import Path
from collections import Iterable, defaultdict

import numpy as np
import pandas as pd

import seaborn as sns
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

import scanpy.logging as logg
from anndata import AnnData

from cellrank.tools._utils import save_fig
from cellrank.utils._utils import _get_n_cores, check_collection
from cellrank.plotting._utils import _fit, _model_type, _create_models, _is_any_gam_mgcv
from cellrank.tools._constants import LinKey
from cellrank.utils._parallelize import parallelize


def heatmap(
    adata: AnnData,
    model: _model_type,
    genes: Sequence[str],
    final: bool = True,
    kind: str = "lineages",
    lineages: Optional[Union[str, Sequence[str]]] = None,
    start_lineage: Optional[Union[str, Sequence[str]]] = None,
    end_lineage: Optional[Union[str, Sequence[str]]] = None,
    lineage_height: float = 0.1,
    cluster_genes: bool = False,
    xlabel: Optional[str] = None,
    cmap: colors.ListedColormap = cm.Spectral_r,
    n_jobs: Optional[int] = 1,
    backend: str = "multiprocessing",
    hspace: float = 0.25,
    figsize: Optional[Tuple[float, float]] = None,
    dpi: Optional[int] = None,
    save: Optional[Union[str, Path]] = None,
    show_progress_bar: bool = True,
    **kwargs,
) -> None:
    """
    Plot a heatmap of smoothed gene expression along specified lineages.

    .. image:: https://raw.githubusercontent.com/theislab/cellrank/master/resources/images/heatmap.png
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
        Genes in :paramref:`adata` `.var_names` to plot.
    final
        Whether to consider cells going to final states or vice versa.
    kind
        Variant of the heatmap.

        - If `'genes'`, group by :paramref:`genes` for each lineage in :paramref:`lineage_names`.
        - If `'lineages'`, group by :paramref:`lineage_names` for each gene in :paramref:`genes`.
    lineage_names
        Names of the lineages for which to plot.
    start_lineage
        Lineage from which to select cells with lowest pseudotime as starting points.
        If specified, the trends start at the earliest pseudotime point within that lineage,
        otherwise they start from time `0`.
    end_lineage
        Lineage from which to select cells with highest pseudotime as endpoints.
        If specified, the trends end at the latest pseudotime point within that lineage,
        otherwise, it is determined automatically.
    lineage_height
        Height of a bar when :paramref:`kind` ='lineages'.
    xlabel
        Label on the x-axis. If `None`, it is determined based on :paramref:`time_key`.
    cluster_genes
        Whether to use :func:`seaborn.clustermap` when :paramref:`kind` `='lineages'`.
    cmap
        Colormap to use when visualizing the smoothed expression.
    n_jobs
        Number of parallel jobs. If `-1`, use all available cores. If `None` or `1`, the execution is sequential.
    backend
        Which backend to use for multiprocessing.
        See :class:`joblib.Parallel` for valid options.
    figsize
        Size of the figure.
        If `None`, it will be set to (15, len(:paramref:`genes`) + len(:paramref:`lineage_names`)).
    dpi
        Dots per inch.
    save
        Filename where to save the plot.
        If `None`, just shows the plot.
    show_progress_bar
        Whether to show a progress bar tracking models fitted.
    **kwargs
        Keyword arguments for :meth:`cellrank.ul.models.Model.prepare`.

    Returns
    -------
    None
        Nothing, just plots the heatmap variant depending on :paramref:`kind`.
        Optionally saves the figure based on :paramref:`save`.
    """

    def gene_per_lineage():
        def color_fill_rec(ax, x, y1, y2, colors=None, cmap=cmap, **kwargs):
            dx = x[1] - x[0]

            for (color, x, y1, y2) in zip(cmap(colors), x, y1, y2):
                ax.add_patch(
                    plt.Rectangle((x, y1), dx, y2 - y1, color=color, ec=color, **kwargs)
                )

            ax.plot(x, y2, lw=0)

        fig, axes = plt.subplots(
            nrows=len(genes),
            figsize=(15, len(genes) + len(lineages)) if figsize is None else figsize,
            dpi=dpi,
        )

        if not isinstance(axes, Iterable):
            axes = [axes]
        axes = np.ravel(axes)

        for ax, (gene, models) in zip(axes, data.items()):
            c = np.array([m.y_test for m in models.values()])
            c_min, c_max = np.nanmin(c), np.nanmax(c)
            norm = colors.Normalize(vmin=c_min, vmax=c_max)

            ix = 0
            ys = [ix]

            for x, c in ((m.x_test, m.y_test) for m in models.values()):
                y = np.ones_like(x)
                color_fill_rec(ax, x, y * ix, y * (ix + lineage_height), colors=norm(c))
                ix += lineage_height
                ys.append(ix)

            xs = np.array([m.x_test for m in models.values()])
            x_min, x_max = np.nanmin(xs), np.nanmax(xs)
            ax.set_xticks(np.linspace(x_min, x_max, 11))

            ax.set_yticks(np.array(ys) + lineage_height / 2)
            ax.set_yticklabels(lineages)
            ax.set_title(gene)
            ax.set_ylabel("lineage")

            for pos in ["top", "bottom", "left", "right"]:
                ax.spines[pos].set_visible(False)

            cax, _ = mpl.colorbar.make_axes(ax)
            _ = mpl.colorbar.ColorbarBase(cax, norm=norm, cmap=cmap, label="Expression")

            ax.tick_params(
                top=False,
                bottom=False,
                left=False,
                right=False,
                labelleft=True,
                labelbottom=False,
            )

        ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
        ax.tick_params(
            top=False,
            bottom=True,
            left=False,
            right=False,
            labelleft=True,
            labelbottom=True,
        )
        ax.set_xlabel(xlabel)

        return fig

    def lineage_per_gene():
        data_t = defaultdict(dict)  # transpose
        for gene, lns in data.items():
            for ln, d in lns.items():
                data_t[ln][gene] = d

        fig, ax = None, None
        if not cluster_genes:
            fig, axes = plt.subplots(
                nrows=len(lineages),
                figsize=(15, 5 + len(genes)) if figsize is None else figsize,
                dpi=dpi,
            )
            fig.subplots_adjust(hspace=hspace, bottom=0)

            if not isinstance(axes, Iterable):
                axes = [axes]
            axes = np.ravel(axes)
        else:
            axes = [None] * len(data)

        for ax, (lname, models) in zip(axes, data_t.items()):
            df = pd.DataFrame([m.y_test for m in models.values()], index=genes)
            df.index.name = "Genes"
            if cluster_genes:
                g = sns.clustermap(
                    df,
                    cmap=cmap,
                    xticklabels=False,
                    cbar_kws={"label": "Expression"},
                    row_cluster=True,
                    col_cluster=False,
                )
                g.ax_heatmap.set_title(lname)
                fig = g.fig
            else:
                xs = np.array([m.x_test for m in models.values()])
                x_min, x_max = np.nanmin(xs), np.nanmax(xs)

                sns.heatmap(
                    df,
                    ax=ax,
                    cmap=cmap,
                    xticklabels=False,
                    cbar_kws={"label": "Expression"},
                )

                ax.set_title(lname)
                ax.set_xticks(np.linspace(0, len(df.columns), 10))
                ax.set_xticklabels(
                    list(map(lambda n: round(n, 3), np.linspace(x_min, x_max, 10))),
                    rotation=90,
                )

        if not cluster_genes:
            ax.set_xlabel(xlabel)

        return fig

    lineage_key = str(LinKey.FORWARD if final else LinKey.BACKWARD)
    if lineage_key not in adata.obsm:
        raise KeyError(f"Lineages key `{lineage_key!r}` not found in `adata.obsm`.")

    if lineages is None:
        lineages = adata.obsm[lineage_key].names

    for lineage_name in lineages:
        _ = adata.obsm[lineage_key][lineage_name]

    if isinstance(genes, str):
        genes = [genes]
    check_collection(adata, genes, "var_names")

    if isinstance(start_lineage, (str, type(None))):
        start_lineage = [start_lineage] * len(lineages)
    if isinstance(end_lineage, (str, type(None))):
        end_lineage = [end_lineage] * len(lineages)

    xlabel = kwargs.get("time_key", None) if xlabel is None else xlabel

    _ = kwargs.pop("start_lineage", None)
    _ = kwargs.pop("end_lineage", None)

    for typp, clusters in zip(["Start", "End"], [start_lineage, end_lineage]):
        for cl in filter(lambda c: c is not None, clusters):
            if cl not in lineages:
                raise ValueError(f"{typp} lineage `{cl!r}` not found in lineage names.")

    kwargs["models"] = _create_models(model, genes, lineages)
    if _is_any_gam_mgcv(kwargs["models"]):
        logg.debug(
            "DEBUG: Setting backend to multiprocessing because model is `GamMGCV`"
        )
        backend = "multiprocessing"

    n_jobs = _get_n_cores(n_jobs, len(genes))
    start = logg.info(f"Computing trends using `{n_jobs}` core(s)")
    data = parallelize(
        _fit,
        genes,
        unit="gene",
        backend=backend,
        n_jobs=n_jobs,
        extractor=lambda data: {k: v for d in data for k, v in d.items()},
        show_progress_bar=show_progress_bar,
    )(lineages, start_lineage, end_lineage, **kwargs)
    logg.info("    Finish", time=start)
    logg.debug(f"DEBUG: Plotting {kind} heatmap")

    if kind == "genes":
        fig = gene_per_lineage()
    elif kind == "lineages":
        fig = lineage_per_gene()
    else:
        raise ValueError(
            f"Unknown heatmap kind `{kind!r}`. Valid options are: `'lineages'`, `'genes'`."
        )

    if save is not None and fig is not None:
        save_fig(fig, save)
