"""Gene trend module."""

from types import MappingProxyType
from typing import List, Tuple, Union, Mapping, TypeVar, Optional, Sequence
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm

from cellrank import logging as logg
from cellrank.ul._docs import d
from cellrank.pl._utils import (
    _fit_bulk,
    _get_backend,
    _callback_type,
    _create_models,
    _trends_helper,
    _time_range_type,
    _create_callbacks,
    _input_model_type,
    _return_model_type,
)
from cellrank.tl._utils import save_fig, _unique_order_preserving
from cellrank.ul._utils import _get_n_cores, _check_collection
from cellrank.tl._constants import _DEFAULT_BACKEND, AbsProbKey

AnnData = TypeVar("AnnData")


@d.dedent
def gene_trends(
    adata: AnnData,
    model: _input_model_type,
    genes: Union[str, Sequence[str]],
    lineages: Optional[Union[str, Sequence[str]]] = None,
    backward: bool = False,
    data_key: str = "X",
    time_key: str = "latent_time",
    transpose: bool = False,
    time_range: Optional[Union[_time_range_type, List[_time_range_type]]] = None,
    callback: _callback_type = None,
    conf_int: Union[bool, float] = True,
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
    cbar: bool = True,
    margins: float = 0.015,
    sharex: Optional[Union[str, bool]] = None,
    sharey: Optional[Union[str, bool]] = None,
    gene_as_title: Optional[bool] = None,
    legend_loc: Optional[str] = "best",
    ncols: int = 2,
    suptitle: Optional[str] = None,
    return_models: bool = False,
    n_jobs: Optional[int] = 1,
    backend: str = _DEFAULT_BACKEND,
    show_progress_bar: bool = True,
    figsize: Optional[Tuple[float, float]] = None,
    dpi: Optional[int] = None,
    save: Optional[Union[str, Path]] = None,
    plot_kwargs: Mapping = MappingProxyType({}),
    **kwargs,
) -> Optional[_return_model_type]:
    """
    Plot gene expression trends along lineages.

    Each lineage is defined via it's lineage weights which we compute using :func:`cellrank.tl.lineages`. This
    function accepts any model based off :class:`cellrank.ul.models.BaseModel` to fit gene expression,
    where we take the lineage weights into account in the loss function.

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
    time_key
        Key in ``adata.obs`` where the pseudotime is stored.
    %(time_ranges)s
    transpose
        If ``same_plot=True``, group the trends by ``lineages`` instead of ``genes``. This enforces ``hide_cells=True``.
        If ``same_plot=False``, show ``lineages`` in rows and ``genes`` in columns.
    %(model_callback)s
    conf_int
        Whether to compute and show confidence interval. If the :paramref:`model` is :class:`cellrank.ul.models.GAMR`,
        it can also specify the confidence level, the default is `0.95`.
    same_plot
        Whether to plot all lineages for each gene in the same plot.
    hide_cells
        If `True`, hide all cells.
    perc
        Percentile for colors. Valid values are in interval `[0, 100]`.
        This can improve visualization. Can be specified individually for each lineage.
    lineage_cmap
        Categorical colormap to use when coloring in the lineages. If `None` and ``same_plot``,
        use the corresponding colors in ``adata.uns``, otherwise use `'black'`.
    abs_prob_cmap
        Continuous colormap to use when visualizing the absorption probabilities for each lineage.
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
    cbar
        Whether to show colorbar. Always shown when percentiles for lineages differ. Only used when ``same_plot=False``.
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
        Number of columns of the plot when plotting multiple genes. Only used when ``same_plot=True``.
    suptitle
        Suptitle of the figure.
    %(return_models)s
    %(parallel)s
    %(plotting)s
    plot_kwargs
        Keyword arguments for :meth:`cellrank.ul.models.BaseModel.plot`.
    kwargs
        Keyword arguments for :meth:`cellrank.ul.models.BaseModel.prepare`.

    Returns
    -------
    %(plots_or_returns_models)s
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
        cbar = False
        logg.debug("All lineages are `None`, setting the weights to `1`")
    lineages = _unique_order_preserving(lineages)

    tmp = adata.obsm[ln_key][lineages].colors
    if lineage_cmap is None and not transpose:
        lineage_cmap = tmp

    if isinstance(time_range, (tuple, float, int, type(None))):
        time_range = [time_range] * len(lineages)
    elif len(time_range) != len(lineages):
        raise ValueError(
            f"Expected time ranges to be of length `{len(lineages)}`, found `{len(time_range)}`."
        )

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

    plot_kwargs = dict(plot_kwargs)
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
            sharey = (
                "row"
                if not hide_cells or plot_kwargs.get("lineage_probability", False)
                else "none"
            )
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
        constrained_layout=True,
        dpi=dpi,
    )
    axes = np.reshape(axes, (nrows, ncols))

    logg.info("Plotting trends")
    cnt = 0

    for row in range(len(axes)):
        for col in range(len(axes[row])):
            if cnt >= len(genes):
                break
            gene = genes[cnt]
            if (
                same_plot
                and plot_kwargs.get("lineage_probability", False)
                and transpose
            ):
                lpc = adata.obsm[ln_key][gene].colors[0]
            else:
                lpc = None

            _trends_helper(
                models,
                gene=gene,
                lineage_names=lineages,
                transpose=transpose,
                same_plot=same_plot,
                hide_cells=hide_cells,
                perc=perc,
                lineage_cmap=lineage_cmap,
                abs_prob_cmap=abs_prob_cmap,
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
                show_xticks_and_label=((row + 1) * ncols + col >= len(genes))
                if same_plot
                else (cnt == end_rows),
                **plot_kwargs,
            )
            cnt += 1

    if same_plot and (col != ncols):
        for ax in np.ravel(axes)[cnt:]:
            ax.remove()

    fig.suptitle(suptitle, y=1.05)

    if save is not None:
        save_fig(fig, save)

    if return_models:
        return all_models
