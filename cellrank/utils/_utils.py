# -*- coding: utf-8 -*-
from cellrank.tools._utils import *  # this prevents the circular imports

from multiprocessing import cpu_count
from typing import Iterable, Hashable, Dict

import anndata
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from cellrank.tools._constants import _colors


def check_collection(
    adata: anndata.AnnData,
    needles: Iterable[str],
    attr_name: str,
    key_name: str = "Gene",
) -> None:
    """
    Check if given collection contains all the keys.

    Params
    ------
    adata: :class:`anndata.AnnData`
        Annotated data object.
    needles
        Keys to check.
    attr_name
        Attribute of :paramref:`adata` where the needles are stored.

    Returns
    -------
    None
    """

    haystack = getattr(adata, attr_name)
    for needle in needles:
        if needle not in haystack:
            raise KeyError(f"{key_name} `{needle}` not found in `adata.{attr_name}`.")


def _get_n_cores(n_cores: Optional[int], n_genes: int) -> int:
    """
    Make number of cores a positive integer.

    Params
    ------
    n_cores
        Number of cores to use.
    n_genes.
        Number of genes.

    Returns
    -------
    int
        Positive integer corresponding to how many cores to use.
    """

    if n_cores == 0:
        raise ValueError("Number of cores cannot be `0`.")
    if n_genes == 1 or n_cores is None:
        return 1
    if n_cores < 0:
        return cpu_count() + 1 + n_cores

    return n_cores


def _minmax(
    data: np.ndarray, perc: Optional[Tuple[float, float]] = None
) -> Tuple[float, float]:
    """
    Return minimum and maximum value of the data.

    Params
    ------
    data
        Values for which to return the minimum and maximum.
    perc
        If not `None`, clip the values by the percentiles before getting the result.

    Returns
    -------
    :class:`tuple`
        Minimum and maximum values, respectively.
    """

    if perc is not None:
        data = np.clip(data, *np.percentile(data, sorted(perc)))

    return float(np.nanmin(data)), float(np.nanmax(data))


def _make_unique(collection: Iterable[Hashable]) -> List[Hashable]:
    """
    Make a collection unique while maintaining the order.

    Params
    ------
    collection
        Values to make unique.

    Returns
    -------
    :class:`list`
        The same collection with unique values.
    """

    seen, res = set(), []
    for item in collection:
        if item not in seen:
            seen.add(item)
            res.append(item)

    return res


def _trends_helper(
    adata: anndata.AnnData,
    models: Dict[str, Dict[str, Any]],
    gene: str,
    ln_key: str,
    lineage_names: Optional[Sequence[str]] = None,
    same_plot: bool = False,
    sharey: bool = True,
    cmap=None,
    fig: mpl.figure.Figure = None,
    ax: mpl.axes.Axes = None,
    save: Optional[str] = None,
    **kwargs,
) -> None:
    """
    Plot an expression gene for some lineages.

    Params
    ------
    adata: :class:`anndata.AnnData`
        Annotated data object.
    models
        Gene and lineage specific models can be specified. Use `'*'` to indicate
        all genes or lineages, for example `{'Map2': {'*': ...}, 'Dcx': {'Alpha': ..., '*': ...}}`.
    gene
        Name of the gene in `adata.var_names`.
    fig
        Figure to use, if `None`, create a new one.
    ax
        Ax to use, if `None`, create a new one.
    save
        Filename where to save the plot.
        If `None`, just shows the plots.
    **kwargs
        Keyword arguments for :meth:`cellrank.ul.models.Model.plot`.

    Returns
    -------
    None
        Nothing, just plots the trends.
        Optionally saves the figure based on :paramref:`save`.
    """

    n_lineages = len(lineage_names)
    if same_plot:
        if fig is None and ax is None:
            fig, ax = plt.subplots(
                1,
                figsize=kwargs.get("figsize", None) or (15, 10),
                constrained_layout=True,
            )
        axes = [ax] * len(lineage_names)
    else:
        fig, axes = plt.subplots(
            ncols=n_lineages,
            figsize=kwargs.get("figsize", None) or (6 * n_lineages, 6),
            sharey=sharey,
            constrained_layout=True,
        )
    axes = np.ravel(axes)
    percs = kwargs.pop("perc", None)
    if percs is None or not isinstance(percs[0], (tuple, list)):
        percs = [percs]

    same_perc = False  # we need to show colorbar always if percs differ
    if len(percs) != n_lineages or n_lineages == 1:
        if len(percs) != 1:
            raise ValueError(
                f"Percentile must be a collection of size `1` or `{n_lineages}`, got `{len(percs)}`."
            )
        same_perc = True
        percs = percs * n_lineages

    hide_cells = kwargs.pop("hide_cells", False)
    show_cbar = kwargs.pop("show_cbar", True)
    lineage_color = kwargs.pop("color", "black")

    lc = (
        cmap.colors
        if cmap is not None and hasattr(cmap, "colors")
        else adata.uns.get(f"{_colors(ln_key)}", cm.Set1.colors)
    )

    for i, (name, ax, perc) in enumerate(zip(lineage_names, axes, percs)):
        title = name if name is not None else "No lineage"
        models[gene][name].plot(
            ax=ax,
            fig=fig,
            perc=perc,
            show_cbar=True
            if not same_perc
            else False
            if not show_cbar
            else (i == n_lineages - 1),
            title=title,
            hide_cells=hide_cells or (same_plot and i != n_lineages - 1),
            same_plot=same_plot,
            color=lc[i] if same_plot and name is not None else lineage_color,
            ylabel=gene if not same_plot or name is None else "Expression",
            **kwargs,
        )

    if same_plot and lineage_names != [None]:
        ax.set_title(gene)
        ax.legend()

    if save is not None:
        save_fig(fig, save)
