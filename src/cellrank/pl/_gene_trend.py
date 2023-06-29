import pathlib
import types
from typing import Any, List, Mapping, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm

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
    _input_model_type,
    _return_model_type,
    _time_range_type,
    _trends_helper,
)

__all__ = ["gene_trends"]


@d.dedent
@_genesymbols
def gene_trends(
    adata: AnnData,
    model: _input_model_type,
    genes: Union[str, Sequence[str]],
    time_key: str,
    lineages: Optional[Union[str, Sequence[str]]] = None,
    backward: bool = False,
    data_key: str = "X",
    time_range: Optional[Union[_time_range_type, List[_time_range_type]]] = None,
    transpose: bool = False,
    callback: _callback_type = None,
    conf_int: Union[bool, float] = True,
    same_plot: bool = False,
    hide_cells: bool = False,
    perc: Optional[Union[Tuple[float, float], Sequence[Tuple[float, float]]]] = None,
    lineage_cmap: Optional[matplotlib.colors.ListedColormap] = None,
    fate_prob_cmap: matplotlib.colors.ListedColormap = cm.viridis,
    cell_color: Optional[str] = None,
    cell_alpha: float = 0.6,
    lineage_alpha: float = 0.2,
    size: float = 15,
    lw: float = 2,
    cbar: bool = True,
    margins: float = 0.015,
    sharex: Optional[Union[str, bool]] = None,
    sharey: Optional[Union[str, bool]] = None,
    gene_as_title: Optional[bool] = None,
    legend_loc: Optional[str] = "best",
    obs_legend_loc: Optional[str] = "best",
    ncols: int = 2,
    suptitle: Optional[str] = None,
    return_models: bool = False,
    n_jobs: Optional[int] = 1,
    backend: Backend_t = DEFAULT_BACKEND,
    show_progress_bar: bool = True,
    figsize: Optional[Tuple[float, float]] = None,
    dpi: Optional[int] = None,
    save: Optional[Union[str, pathlib.Path]] = None,
    return_figure: bool = False,
    plot_kwargs: Mapping[str, Any] = types.MappingProxyType({}),
    **kwargs: Any,
) -> Optional[_return_model_type]:
    """Plot gene expression trends along lineages.

    .. seealso::
        - See :doc:`../../../notebooks/tutorials/estimators/800_gene_trends` on how to
          visualize the gene trends.

    Each lineage is defined via its lineage weights. This function accepts any model based off
    :class:`~cellrank.models.BaseModel` to fit gene expression, where we take the lineage weights
    into account in the loss function.

    Parameters
    ----------
    %(adata)s
    %(model)s
    %(genes)s
    time_key
        Key in :attr:`~anndata.AnnData.obs` where the pseudotime is stored.
    lineages
        Names of the lineages to plot. If :obj:`None`, plot all lineages.
    %(backward)s
    data_key
        Key in :attr:`~anndata.AnnData.layers` or ``'X'`` for :attr:`~anndata.AnnData.X` where the data is stored.
    %(time_range)s
        This can also be specified on per-lineage basis.
    %(gene_symbols)s
    transpose
        If ``same_plot = True``, group the trends by ``lineages`` instead of ``genes``.
        This forces ``hide_cells = True``.
        If ``same_plot = False``, show ``lineages`` in rows and ``genes`` in columns.
    %(model_callback)s
    conf_int
        Whether to compute and show confidence interval. If the ``model`` is :class:`~cellrank.models.GAMR`,
        it can also specify the confidence level, the default is :math:`0.95`.
    same_plot
        Whether to plot all lineages for each gene in the same plot.
    hide_cells
        If :obj:`True`, hide all cells.
    perc
        Percentile for colors. Valid values are in :math:`[0, 100]`.
        This can improve visualization. Can be specified individually for each lineage.
    lineage_cmap
        Categorical colormap to use when coloring in the lineages. If :obj:`None` and ``same_plot = True``,
        use the corresponding colors in :attr:`~anndata.AnnData.uns`, otherwise use ``'black'``.
    fate_prob_cmap
        Continuous colormap to use when visualizing the fate probabilities for each lineage.
        Only used when ``same_plot = False``.
    cell_color
        Key in :attr:`~anndata.AnnData.obs` or :attr:`~anndata.AnnData.var_names` used for coloring the cells.
    cell_alpha
        Alpha channel for cells.
    lineage_alpha
        Alpha channel for lineage confidence intervals.
    size
        Size of the points.
    lw
        Line width of the smoothed values.
    cbar
        Whether to show colorbar. Always shown when percentiles for lineages differ.
        Only used when ``same_plot = False``.
    margins
        Margins around the plot.
    sharex
        Whether to share x-axis. Valid options are ``'row'``, ``'col'`` or ``'none'``.
    sharey
        Whether to share y-axis. Valid options are ``'row'`, ``'col'`` or ``'none'``.
    gene_as_title
        Whether to show gene names as titles instead on y-axis.
    legend_loc
        Location of the legend displaying lineages. Only used when ``same_plot = True``.
    obs_legend_loc
        Location of the legend when ``cell_color`` corresponds to a categorical variable.
    ncols
        Number of columns of the plot when plotting multiple genes. Only used when ``same_plot = True``.
    suptitle
        Suptitle of the figure.
    return_figure
        Whether to return the figure object.
    %(return_models)s
    %(parallel)s
    %(plotting)s
    plot_kwargs
        Keyword arguments for the :meth:`~cellrank.models.BaseModel.plot`.
    kwargs
        Keyword arguments for :meth:`~cellrank.models.BaseModel.prepare`.

    Returns
    -------
    %(plots_or_returns_models)s
    """
    if isinstance(genes, str):
        genes = [genes]
    genes = _unique_order_preserving(genes)

    _check_collection(
        adata,
        genes,
        "obs" if data_key == "obs" else "var_names",
        use_raw=kwargs.get("use_raw", False),
    )

    probs = Lineage.from_adata(adata, backward=backward)
    if lineages is None:
        lineages = probs.names
    elif isinstance(lineages, str):
        lineages = [lineages]
    elif all(ln is None for ln in lineages):  # no lineage, all the weights are 1
        lineages = [None]
        cbar = False
        logg.debug("All lineages are `None`, setting the weights to `1`")
    lineages = _unique_order_preserving(lineages)

    if isinstance(time_range, (tuple, float, int, type(None))):
        time_range = [time_range] * len(lineages)
    elif len(time_range) != len(lineages):
        raise ValueError(f"Expected time ranges to be of length `{len(lineages)}`, found `{len(time_range)}`.")

    kwargs["time_key"] = time_key
    kwargs["data_key"] = data_key
    kwargs["backward"] = backward
    kwargs["conf_int"] = conf_int  # prepare doesnt take or need this
    models = _create_models(model, genes, lineages)

    all_models, models, genes, lineages = _fit_bulk(
        models,
        _create_callbacks(adata, callback, genes, lineages, **kwargs),
        genes,
        lineages,
        time_range,
        return_models=True,
        filter_all_failed=False,
        parallel_kwargs={
            "show_progress_bar": show_progress_bar,
            "n_jobs": _get_n_cores(n_jobs, len(genes)),
            "backend": _get_backend(models, backend),
        },
        **kwargs,
    )

    lineages = sorted(lineages)
    probs = probs[lineages]
    if lineage_cmap is None and not transpose:
        lineage_cmap = probs.colors

    plot_kwargs = dict(plot_kwargs)
    plot_kwargs["obs_legend_loc"] = obs_legend_loc
    if transpose:
        all_models = pd.DataFrame(all_models).T.to_dict()
        models = pd.DataFrame(models).T.to_dict()
        genes, lineages = lineages, genes
        hide_cells = same_plot or hide_cells
    else:
        # information overload otherwise
        plot_kwargs["lineage_probability"] = False
        plot_kwargs["lineage_probability_conf_int"] = False

    tmp = pd.DataFrame(models).T.astype(bool)
    start_rows = np.argmax(tmp.values, axis=0)
    end_rows = tmp.shape[0] - np.argmax(tmp[::-1].values, axis=0) - 1

    if same_plot:
        gene_as_title = True if gene_as_title is None else gene_as_title
        sharex = "all" if sharex is None else sharex
        if sharey is None:
            sharey = "row" if plot_kwargs.get("lineage_probability", False) else "none"
        ncols = len(genes) if ncols >= len(genes) else ncols
        nrows = int(np.ceil(len(genes) / ncols))
    else:
        gene_as_title = False if gene_as_title is None else gene_as_title
        sharex = "col" if sharex is None else sharex
        if sharey is None:
            sharey = "row" if not hide_cells or plot_kwargs.get("lineage_probability", False) else "none"
        nrows = len(genes)
        ncols = len(lineages)

    plot_kwargs = dict(plot_kwargs)
    if plot_kwargs.get("xlabel", None) is None:
        plot_kwargs["xlabel"] = time_key

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        sharex=sharex,
        sharey=sharey,
        figsize=(6 * ncols, 4 * nrows) if figsize is None else figsize,
        tight_layout=True,
        dpi=dpi,
    )
    axes = np.reshape(axes, (nrows, ncols))

    cnt = 0
    plot_kwargs["obs_legend_loc"] = None if same_plot else obs_legend_loc

    logg.info("Plotting trends")
    for row in range(len(axes)):
        for col in range(len(axes[row])):
            if cnt >= len(genes):
                break
            gene = genes[cnt]
            if same_plot and plot_kwargs.get("lineage_probability", False) and transpose:
                lpc = probs[gene].colors[0]
            else:
                lpc = None

            if same_plot:
                plot_kwargs["obs_legend_loc"] = obs_legend_loc if row == 0 and col == len(axes[0]) - 1 else None

            _trends_helper(
                models,
                gene=gene,
                lineage_names=lineages,
                transpose=transpose,
                same_plot=same_plot,
                hide_cells=hide_cells,
                perc=perc,
                lineage_cmap=lineage_cmap,
                fate_prob_cmap=fate_prob_cmap,
                lineage_probability_color=lpc,
                cell_color=cell_color,
                alpha=cell_alpha,
                lineage_alpha=lineage_alpha,
                size=size,
                lw=lw,
                cbar=cbar,
                margins=margins,
                sharey=sharey,
                gene_as_title=gene_as_title,
                legend_loc=legend_loc,
                figsize=figsize,
                fig=fig,
                axes=axes[row, col] if same_plot else axes[cnt],
                show_ylabel=col == 0,
                show_lineage=same_plot or (cnt == start_rows),
                show_xticks_and_label=((row + 1) * ncols + col >= len(genes)) if same_plot else (cnt == end_rows),
                **plot_kwargs,
            )
            cnt += 1  # plot legend on the 1st plot

            if not same_plot:
                plot_kwargs["obs_legend_loc"] = None

    if same_plot and (col != ncols):
        for ax in np.ravel(axes)[cnt:]:
            ax.remove()

    fig.suptitle(suptitle, y=1.05)

    if return_figure:
        return fig

    if save is not None:
        save_fig(fig, save)

    if return_models:
        return all_models
