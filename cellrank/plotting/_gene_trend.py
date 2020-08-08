# -*- coding: utf-8 -*-
"""Gene trend module."""

import os
from types import MappingProxyType
from typing import List, Tuple, Union, Mapping, TypeVar, Optional, Sequence
from pathlib import Path

import numpy as np

import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from cellrank import logging as logg
from cellrank.utils._docs import d
from cellrank.tools._utils import save_fig, _unique_order_preserving
from cellrank.utils._utils import _get_n_cores, check_collection
from cellrank.plotting._utils import (
    _model_type,
    _get_backend,
    _create_models,
    _trends_helper,
    _fit_gene_trends,
    _time_range_type,
    _maybe_create_dir,
)
from cellrank.tools._constants import AbsProbKey
from cellrank.utils._parallelize import parallelize

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
    conf_int: bool = True,
    same_plot: bool = False,
    hide_cells: bool = False,
    perc: Optional[Union[Tuple[float, float], Sequence[Tuple[float, float]]]] = None,
    lineage_cmap: Optional[matplotlib.colors.ListedColormap] = None,
    abs_prob_cmap: matplotlib.colors.ListedColormap = cm.viridis,
    cell_color: str = "black",
    color: str = "black",
    cell_alpha: float = 0.6,
    lineage_alpha: float = 0.2,
    size: float = 15,
    lw: float = 2,
    show_cbar: bool = True,
    margins: float = 0.015,
    sharey: bool = False,
    ncols: int = 2,
    ext: str = "png",
    suptitle: Optional[str] = None,
    n_jobs: Optional[int] = 1,
    backend: str = "multiprocessing",
    show_progres_bar: bool = True,
    dirname: Optional[str] = None,
    figsize: Optional[Tuple[float, float]] = None,
    dpi: Optional[int] = None,
    save: Optional[Union[str, Path]] = None,
    plot_kwargs: Mapping = MappingProxyType({}),
    **kwargs,
) -> None:
    """
    Plot gene expression trends along lineages.

    Each lineage is defined via it's lineage weights which we compute using :func:`cellrank.tl.lineages`. This
    function accepts any `scikit-learn` model wrapped in :class:`cellrank.ul.models.SKLearnModel`
    to fit gene expression, where we take the lineage weights into account in the loss function.

    .. image:: https://raw.githubusercontent.com/theislab/cellrank/master/resources/images/gene_trends.png
       :width: 400px
       :align: center

    Parameters
    ----------
    %(adata)s
    %(model)s
    genes
        Genes in :paramref:`adata` `.var_names` to plot or in :paramref:`adata` `.raw.var_names`, if `use_raw=True`.
    lineages
        Names of the lineages which to plot. If `None`, plot all lineages.
    %(backward)s
    data_key
        Key in :paramref:`adata` `.layers` or `'X'` for :paramref:`adata` `.X` where the data is stored.
    %(time_range)s
    conf_int
        Whether to compute and show confidence intervals.
    same_plot
        Whether to plot all lineages for each gene in the same plot.
    hide_cells
        If `True`, hide all the cells.
    perc
        Percentile for colors. Valid values are in range `[0, 100]`.
        This can improve visualization. Can be specified separately for each lineage separately.
    lineage_cmap
        Colormap to use when coloring in the lineages. Only used when :paramref:`same_plot` `=True`.
    abs_prob_cmap
        Colormap to use when visualizing the absorption probabilities for each lineage.
        Only used when :paramref:`same_plot` `=False`.
    cell_color
        Color of the cells when not visualizing absorption probabilities.
        Only used when :paramref:`same_plot` `=True`.
    color
        Color for the lineages, when each lineage is on
        separate plot, otherwise according to :paramref:`lineage_cmap`.
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
    margins
        Margins around the plot.
    sharey
        Whether to share y-axis. Only used when :paramref:`same_plot` `=False`.
    ncols
        Number of columns of the plot when plotting multiple genes.
        Only used when :paramref:`same_plot` `=True`.
    %(parallel)s
    suptitle
        Suptitle of the figure. Only used when :paramref:`same_plot` `=True`.
    ext
        Extension to use when saving files, such as `'pdf'`.
        Only used when :paramref:`same_plot` `=False`.
    dirname
        Directory where to save the plots, one per gene in :paramref:`genes`. If `None`, just show the plots.
        Only used when :paramref:`same_plot` `=False`.

        The figures will be saved as :paramref:`dirname` /`{gene}`. :paramref:`ext`.
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
        check_collection(
            adata, genes, "var_names", use_raw=kwargs.get("use_raw", False)
        )
    else:
        check_collection(adata, genes, "obs", use_raw=kwargs.get("use_raw", False))

    nrows = int(np.ceil(len(genes) / ncols))
    fig = None
    axes = [None] * len(genes)

    if same_plot:
        fig, axes = plt.subplots(
            nrows=nrows,
            ncols=ncols,
            sharey=sharey,
            figsize=(15 * ncols, 10 * nrows) if figsize is None else figsize,
        )
        axes = np.ravel(axes)
    elif dirname is not None:
        _maybe_create_dir(dirname)
    elif save is not None:
        logg.warning("No directory specified for saving. Ignoring `save` argument")

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

    _ = adata.obsm[ln_key][[lin for lin in lineages if lin is not None]]

    if isinstance(time_range, (tuple, float, int, type(None))):
        time_range = [time_range] * len(lineages)
    elif len(time_range) != len(lineages):
        raise ValueError(
            f"Expected time ranges to be of length `{len(lineages)}`, found `{len(time_range)}`."
        )

    kwargs["models"] = _create_models(model, genes, lineages)
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
    )(lineages, time_range, **kwargs)
    logg.info("    Finish", time=start)

    logg.debug("Plotting trends")
    for gene, ax in zip(genes, axes):
        f = (
            None
            if (same_plot or dirname is None)
            else os.path.join(dirname, f"{gene}.{ext}")
        )
        _trends_helper(
            adata,
            models,
            gene=gene,
            lineage_names=lineages,
            ln_key=ln_key,
            same_plot=same_plot,
            hide_cells=hide_cells,
            perc=perc,
            cmap=lineage_cmap,
            abs_prob_cmap=abs_prob_cmap,
            cell_color=cell_color,
            color=color,
            alpha=cell_alpha,
            lineage_alpha=lineage_alpha,
            size=size,
            lw=lw,
            show_cbar=show_cbar,
            margins=margins,
            sharey=sharey,
            dpi=dpi,
            figsize=figsize,
            fig=fig,
            ax=ax,
            save=f,
            **plot_kwargs,
        )

    if same_plot:
        for j in range(len(genes), len(axes)):
            axes[j].remove()

        fig.suptitle(suptitle)

        if save is not None:
            save_fig(fig, save)
