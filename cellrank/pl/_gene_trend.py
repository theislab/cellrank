# -*- coding: utf-8 -*-
"""Gene trend module."""

from types import MappingProxyType
from typing import List, Tuple, Union, Mapping, TypeVar, Optional, Sequence
from pathlib import Path

import numpy as np

import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from cellrank import logging as logg
from cellrank.ul._docs import d
from cellrank.pl._utils import (
    _model_type,
    _get_backend,
    _callback_type,
    _create_models,
    _trends_helper,
    _fit_gene_trends,
    _time_range_type,
    _create_callbacks,
)
from cellrank.tl._utils import save_fig, _unique_order_preserving
from cellrank.ul._utils import _get_n_cores, _check_collection
from cellrank.tl._constants import _DEFAULT_BACKEND, AbsProbKey
from cellrank.ul._parallelize import parallelize

AnnData = TypeVar("AnnData")


@d.dedent
def gene_trends(
    adata: AnnData,
    model: _model_type,
    genes: Union[str, Sequence[str]],
    lineages: Optional[Union[str, Sequence[str]]] = None,
    backward: bool = False,
    data_key: str = "X",
    time_range: Optional[Union[_time_range_type, List[_time_range_type]]] = None,
    callback: _callback_type = None,
    conf_int: bool = True,
    same_plot: bool = False,
    hide_cells: bool = False,
    perc: Optional[Union[Tuple[float, float], Sequence[Tuple[float, float]]]] = None,
    lineage_cmap: Optional[matplotlib.colors.ListedColormap] = None,
    abs_prob_cmap: matplotlib.colors.ListedColormap = cm.viridis,
    cell_color: str = "black",
    cell_alpha: float = 0.6,
    lineage_alpha: float = 0.2,
    size: float = 15,
    lw: float = 2,
    show_cbar: bool = True,
    margins: float = 0.015,
    sharex: Optional[Union[str, bool]] = None,
    sharey: Optional[Union[str, bool]] = None,
    gene_as_title: Optional[bool] = None,
    legend_loc: Optional[str] = "best",
    ncols: int = 2,
    suptitle: Optional[str] = None,
    n_jobs: Optional[int] = 1,
    backend: str = _DEFAULT_BACKEND,
    show_progres_bar: bool = True,
    figsize: Optional[Tuple[float, float]] = None,
    dpi: Optional[int] = None,
    save: Optional[Union[str, Path]] = None,
    plot_kwargs: Mapping = MappingProxyType({}),
    **kwargs,
) -> None:
    """
    Plot gene expression trends along lineages.

    Each lineage is defined via it's lineage weights which we compute using :func:`cellrank.tl.lineages`. This
    function accepts any model based off :class:`cellrank.ul.models.BaseModel` to fit gene expression,
    where we take the lineage weights into account in the loss function.

    .. image:: https://raw.githubusercontent.com/theislab/cellrank/master/resources/images/gene_trends.png
       :width: 400px
       :align: center

    Parameters
    ----------
    %(adata)s
    %(model)s
    %(genes)s
    lineages
        Names of the lineages to plot. If `None`, plot all lineages.
    %(backward)s
    data_key
        Key in ``adata.layers`` or `'X'` for ``adata.X`` where the data is stored.
    %(time_ranges)s
    %(model_callback)s
    conf_int
        Whether to compute and show confidence intervals.
    same_plot
        Whether to plot all lineages for each gene in the same plot.
    hide_cells
        If `True`, hide all the cells.
    perc
        Percentile for colors. Valid values are in `[0, 100]`.
        This can improve visualization. Can be specified indivudually for each lineage.
    lineage_cmap
        Colormap to use when coloring in the lineages. If `None` and ``same_plot``, use the corresponding colors
        in ``adata.uns``, otherwise use `'black'`.
    abs_prob_cmap
        Colormap to use when visualizing the absorption probabilities for each lineage.
        Only used when ``same_plot=False``.
    cell_color
        Color of the cells when not visualizing absorption probabilities. Only used when ``same_plot=True``.
    cell_alpha
        Alpha channel for cells.
    lineage_alpha
        Alpha channel for lineage confidence intervals.
    size
        Size of the points.
    lw
        Line width of the smoothed values.
    show_cbar
        Whether to show colorbar. Always shown when percentiles for lineages differ.
        Only used when ``same_plot=False``.
    margins
        Margins around the plot.
    sharex
        Whether to share x-axis. Valid options are `'row'`, `'col'` or `'none'`.
    sharey
        Whether to share y-axis. Valid options are `'row'`, `'col'` or `'none'`.
    gene_as_title
        Whether to show gene names as titles instead on y-axis.
    legend_loc
        Location of the legend displaying lineages. Only used when `same_plot=True`.
    ncols
        Number of columns of the plot when pl multiple genes. Only used when ``same_plot=True``.
    suptitle
        Suptitle of the figure.
    %(parallel)s
    %(plotting)s
    plot_kwargs
        Keyword arguments for :meth:`cellrank.ul.models.BaseModel.plot`.
    **kwargs
        Keyword arguments for :meth:`cellrank.ul.models.BaseModel.prepare`.

    Returns
    -------
    %(just_plots)s
    """

    if isinstance(genes, str):
        genes = [genes]
    genes = _unique_order_preserving(genes)

    if data_key != "obs":
        _check_collection(
            adata, genes, "var_names", use_raw=kwargs.get("use_raw", False)
        )
    else:
        _check_collection(adata, genes, "obs", use_raw=kwargs.get("use_raw", False))

    ln_key = str(AbsProbKey.BACKWARD if backward else AbsProbKey.FORWARD)
    if ln_key not in adata.obsm:
        raise KeyError(f"Lineages key `{ln_key!r}` not found in `adata.obsm`.")

    if lineages is None:
        lineages = adata.obsm[ln_key].names
    elif isinstance(lineages, str):
        lineages = [lineages]
    elif all(map(lambda ln: ln is None, lineages)):  # no lineage, all the weights are 1
        lineages = [None]
        show_cbar = False
        logg.debug("All lineages are `None`, setting the weights to `1`")
    lineages = _unique_order_preserving(lineages)

    if same_plot:
        gene_as_title = True if gene_as_title is None else gene_as_title
        sharex = "all" if sharex is None else sharex
        sharey = "none" if sharey is None else sharey
        ncols = len(genes) if ncols >= len(genes) else ncols
        nrows = int(np.ceil(len(genes) / ncols))
    else:
        gene_as_title = False if gene_as_title is None else gene_as_title
        sharex = "col" if sharex is None else sharex
        sharey = ("none" if hide_cells else "row") if sharey is None else sharey
        nrows = len(genes)
        ncols = len(lineages)

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        sharex=sharex,
        sharey=sharey,
        figsize=(6 * ncols, 4 * nrows) if figsize is None else figsize,
        constrained_layout=True,
    )
    axes = np.reshape(axes, (-1, ncols))

    _ = adata.obsm[ln_key][[lin for lin in lineages if lin is not None]]

    if isinstance(time_range, (tuple, float, int, type(None))):
        time_range = [time_range] * len(lineages)
    elif len(time_range) != len(lineages):
        raise ValueError(
            f"Expected time ranges to be of length `{len(lineages)}`, found `{len(time_range)}`."
        )

    models = _create_models(model, genes, lineages)
    callbacks = _create_callbacks(adata, callback, genes, lineages)

    kwargs["data_key"] = data_key
    kwargs["backward"] = backward
    kwargs["conf_int"] = conf_int

    plot_kwargs = dict(plot_kwargs)
    if plot_kwargs.get("xlabel", None) is None:
        plot_kwargs["xlabel"] = kwargs.get("time_key", None)

    n_jobs = _get_n_cores(n_jobs, len(genes))
    backend = _get_backend(model, backend)

    start = logg.info(f"Computing trends using `{n_jobs}` core(s)")
    models = parallelize(
        _fit_gene_trends,
        genes,
        unit="gene" if data_key != "obs" else "obs",
        backend=backend,
        n_jobs=n_jobs,
        extractor=lambda modelss: {k: v for m in modelss for k, v in m.items()},
        show_progress_bar=show_progres_bar,
    )(models, callbacks, lineages, time_range, **kwargs)
    logg.info("    Finish", time=start)

    logg.info("Plotting trends")

    cnt = 0
    for row in range(len(axes)):
        for col in range(len(axes[row])):
            if cnt >= len(genes):
                break
            gene = genes[cnt]

            _trends_helper(
                adata,
                models,
                gene=gene,
                lineage_names=lineages,
                ln_key=ln_key,
                same_plot=same_plot,
                hide_cells=hide_cells,
                perc=perc,
                lineage_cmap=lineage_cmap,
                abs_prob_cmap=abs_prob_cmap,
                cell_color=cell_color,
                alpha=cell_alpha,
                lineage_alpha=lineage_alpha,
                size=size,
                lw=lw,
                show_cbar=show_cbar,
                margins=margins,
                sharey=sharey,
                gene_as_title=gene_as_title,
                legend_loc=legend_loc,
                dpi=dpi,
                figsize=figsize,
                fig=fig,
                axes=axes[row, col] if same_plot else axes[cnt],
                show_ylabel=col == 0,
                show_lineage=cnt == 0 or same_plot,
                show_xticks_and_label=((row + 1) * ncols + col >= len(genes))
                if same_plot
                else (cnt == len(axes) - 1),
                **plot_kwargs,
            )
            cnt += 1

    if same_plot and (col != ncols):
        for ax in np.ravel(axes)[cnt:]:
            ax.remove()

    fig.suptitle(suptitle)

    if save is not None:
        save_fig(fig, save)
