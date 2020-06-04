# -*- coding: utf-8 -*-
"""Gene trend module."""

import os
from types import MappingProxyType
from typing import Tuple, Union, Mapping, Optional, Sequence
from pathlib import Path

import numpy as np

import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt

import scanpy as sc
from scanpy import logging as logg
from anndata import AnnData

from cellrank.tools._utils import save_fig
from cellrank.utils._utils import _get_n_cores, _make_unique, check_collection
from cellrank.plotting._utils import (
    _fit,
    _model_type,
    _create_models,
    _trends_helper,
    _is_any_gam_mgcv,
)
from cellrank.tools._constants import LinKey
from cellrank.utils._parallelize import parallelize


def gene_trends(
    adata: AnnData,
    model: _model_type,
    genes: Union[str, Sequence[str]],
    lineages: Optional[Union[str, Sequence[str]]] = None,
    data_key: str = "X",
    final: bool = True,
    start_lineage: Optional[Union[str, Sequence[str]]] = None,
    end_lineage: Optional[Union[str, Sequence[str]]] = None,
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
    figsize: Optional[Tuple[float, float]] = None,
    dpi: Optional[int] = None,
    ncols: int = 2,
    n_jobs: Optional[int] = 1,
    backend: str = "multiprocessing",
    ext: str = "png",
    suptitle: Optional[str] = None,
    save: Optional[Union[str, Path]] = None,
    dirname: Optional[str] = None,
    plot_kwargs: Mapping = MappingProxyType({}),
    show_progres_bar: bool = True,
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

    Params
    ------
    adata : :class:`anndata.AnnData`
        Annotated data object.
    genes
        Genes in :paramref:`adata` `.var_names` to plot.
    model
        Model to fit.

        - If a :class:`dict`, gene and lineage specific models can be specified. Use `'*'` to indicate
        all genes or lineages, for example `{'Map2': {'*': ...}, 'Dcx': {'Alpha': ..., '*': ...}}`.
    lineage_names
        Lineages names for which to show the gene expression.
    data_key
        Key in :paramref:`adata` `.layers` or `'X'` for :paramref:`adata` `.X` where the data is stored.
    final
        Whether to consider cells going to final states or vice versa.
    start_lineage
        Lineage from which to select cells with lowest pseudotime as starting points.
        If specified, the trends start at the earliest pseudotime within that lineage,
        otherwise they start from time `0`.
    end_lineage
        Lineage from which to select cells with highest pseudotime as endpoints.
        If specified, the trends end at the latest pseudotime within that lineage,
        otherwise, it is determined automatically.
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
        Colormap to use when coloring in the lineages.
        Only used when :paramref:`same_plot` `=True`.
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
        Whether to share y-axis.
        Only used when :paramref:`same_plot` `=False`.
    figsize
        Size of the figure.
    dpi
        Dots per inch.
    ncols
        Number of columns of the plot when plotting multiple genes.
        Only used when :paramref:`same_plot` `=True`.
    suptitle
        Suptitle of the figure.
        Only used when :paramref:`same_plot` `=True`.
    n_jobs
        Number of parallel jobs. If `-1`, use all available cores. If `None` or `1`, the execution is sequential.
    backend
        Which backend to use for multiprocessing.
        See :class:`joblib.Parallel` for valid options.
    ext
        Extension to use when saving files, such as `'pdf'`.
        Only used when :paramref:`same_plot` `=False`.
    save
        Filename where to save the plots.
        If `None`, just show the plots.
    dirname
        Directory where to save the plots, one per gene in :paramref:`genes`.
        If `None`, just show the plots.
        Only used when :paramref:`same_plot` `=False`.
        The figures will be saved as :paramref:`dirname` /`{gene}`. :paramref:`ext`.
    plot_kwargs:
        Keyword arguments for :meth:`cellrank.ul.models.Model.plot`.
    kwargs
        Keyword arguments for :meth:`cellrank.ul.models.Model.prepare`.

    Returns
    -------
    None
        Nothings just plots and optionally saves the plots.
    """

    if isinstance(genes, str):
        genes = [genes]
    genes = _make_unique(genes)

    if data_key != "obs":
        check_collection(adata, genes, "var_names")
    else:
        check_collection(adata, genes, "obs")

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
        figdir = sc.settings.figdir
        if figdir is None:
            raise RuntimeError(
                f"Invalid combination: figures directory `cellrank.settings.figdir` is `None`, "
                f"but `dirname={dirname}`."
            )
        if os.path.isabs(dirname):
            if not os.path.isdir(dirname):
                os.makedirs(dirname, exist_ok=True)
        elif not os.path.isdir(os.path.join(figdir, dirname)):
            os.makedirs(os.path.join(figdir, dirname), exist_ok=True)
    elif save is not None:
        logg.warning("No directory specified for saving. Ignoring `save` argument")

    ln_key = str(LinKey.FORWARD if final else LinKey.BACKWARD)
    if ln_key not in adata.obsm:
        raise KeyError(f"Lineages key `{ln_key!r}` not found in `adata.obsm`.")

    if lineages is None:
        lineages = adata.obsm[ln_key].names
    elif isinstance(lineages, str):
        lineages = [lineages]
    elif all(map(lambda ln: ln is None, lineages)):  # no lineage, all the weights are 1
        lineages = [None]
        show_cbar = False
        logg.debug("DEBUG: All lineages are `None`, setting weights to be `1`")
    lineages = _make_unique(lineages)

    for ln in filter(lambda ln: ln is not None, lineages):
        _ = adata.obsm[ln_key][ln]
    n_lineages = len(lineages)

    if isinstance(start_lineage, (str, type(None))):
        start_lineage = [start_lineage] * n_lineages
    if isinstance(end_lineage, (str, type(None))):
        end_lineage = [end_lineage] * n_lineages

    if len(start_lineage) != n_lineages:
        raise ValueError(
            f"Expected the number of start lineages to be the same as number of lineages "
            f"({n_lineages}), found `{len(start_lineage)}`."
        )
    if len(end_lineage) != n_lineages:
        raise ValueError(
            f"Expected the number of end lineages to be the same as number of lineages "
            f"({n_lineages}), found `{len(start_lineage)}`."
        )

    kwargs["models"] = _create_models(model, genes, lineages)
    kwargs["data_key"] = data_key
    kwargs["final"] = final
    kwargs["conf_int"] = conf_int

    plot_kwargs = dict(plot_kwargs)
    if plot_kwargs.get("xlabel", None) is None:
        plot_kwargs["xlabel"] = kwargs.get("time_key", None)

    if _is_any_gam_mgcv(kwargs["models"]):
        logg.debug(
            "DEBUG: Setting backend to multiprocessing because model is `GamMGCV`"
        )
        backend = "multiprocessing"

    n_jobs = _get_n_cores(n_jobs, len(genes))

    start = logg.info(f"Computing trends using `{n_jobs}` core(s)")
    models = parallelize(
        _fit,
        genes,
        unit="gene" if data_key != "obs" else "obs",
        backend=backend,
        n_jobs=n_jobs,
        extractor=lambda modelss: {k: v for m in modelss for k, v in m.items()},
        show_progress_bar=show_progres_bar,
    )(lineages, start_lineage, end_lineage, **kwargs)
    logg.info("    Finish", time=start)

    logg.debug("DEBUG: Plotting trends")
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
