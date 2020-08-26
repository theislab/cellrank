# -*- coding: utf-8 -*-
"""Utility functions for CellRank pl."""

import os
from copy import copy
from typing import (
    Any,
    Dict,
    Tuple,
    Union,
    Mapping,
    TypeVar,
    Callable,
    Optional,
    Sequence,
)
from pathlib import Path
from collections import defaultdict

import numpy as np
from pandas.core.dtypes.common import is_categorical_dtype

import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from cellrank import logging as logg
from cellrank.ul._docs import d
from cellrank.tl._utils import save_fig, _unique_order_preserving
from cellrank.ul.models import GAMR, BaseModel, SKLearnModel
from cellrank.tl._constants import _DEFAULT_BACKEND, _colors

AnnData = TypeVar("AnnData")


_ERROR_INCOMPLETE_SPEC = (
    "No options were specified for{}`{!r}`. "
    "Consider specifying a fallback model using '*'."
)
_time_range_type = Optional[Union[float, Tuple[Optional[float], Optional[float]]]]
_model_type = Union[BaseModel, Mapping[str, Mapping[str, BaseModel]]]
_callback_type = Optional[Union[Callable, Mapping[str, Mapping[str, Callable]]]]


def _curved_edges(
    G,
    pos,
    radius_fraction: float,
    dist_ratio: float = 0.2,
    bezier_precision: int = 20,
    polarity: str = "directed",
) -> np.ndarray:
    """
    Create curved edges from a graph. Modified from: https://github.com/beyondbeneath/bezier-curved-edges-networkx.

    Parameters
    ----------
    G: :class:`networkx.Graph`
        Graph for which to create curved edges.
    pos
        Mapping of nodes to positions.
    radius_fraction
        Fraction of a unit circle when self loops are present.
    dist_ratio
        Distance of control points of bezier curves.
    bezier_precision
        Number of points in the curves.
    polarity
        Polarity of curves, one of `'random', 'directed' or 'fixed'`.
        If using `'random'`, incoming and outgoing edges may overlap.

    Returns
    -------
    :class:`np.ndarray`
        Array of shape ``(n_edges, bezier_precision, 2)`` containing the curved edges.
    """

    try:
        import bezier
    except ImportError as e:
        raise ImportError("Please install `bezier` as `pip install bezier`.") from e

    # Get nodes into np array
    edges = np.array(G.edges())
    n_edges = edges.shape[0]

    self_loop_mask = edges[:, 0] == edges[:, 1]
    pos_sl = {edge[0]: pos[edge[0]] for edge in edges[self_loop_mask, ...]}

    if polarity == "random":
        # Random polarity of curve
        rnd = np.where(np.random.randint(2, size=n_edges) == 0, -1, 1)
    elif polarity == "directed":
        rnd = np.where(edges[:, 0] > edges[:, 1], -1, 1)
    elif polarity == "fixed":
        # Create a fixed (hashed) polarity column in the case we use fixed polarity
        # This is useful, e.g., for animations
        rnd = np.where(
            np.mod(np.vectorize(hash)(edges[:, 0]) + np.vectorize(hash)(edges[:, 1]), 2)
            == 0,
            -1,
            1,
        )
    else:
        raise ValueError(
            f"Polarity `{polarity!r}` is not a valid option. "
            f"Valid options are: `'random', 'directed' or 'fixed'`."
        )

    # Coordinates (x, y) of both nodes for each edge
    # Note the np.vectorize method doesn't work for all node position dictionaries for some reason
    u, inv = np.unique(edges, return_inverse=True)
    coords = np.array([pos[x] for x in u])[inv].reshape(
        [edges.shape[0], 2, edges.shape[1]]
    )
    coords_node1 = coords[:, 0, :]
    coords_node2 = coords[:, 1, :]

    # Swap node1/node2 allocations to make sure the directionality works correctly
    should_swap = coords_node1[:, 0] > coords_node2[:, 0]
    coords_node1[should_swap], coords_node2[should_swap] = (
        coords_node2[should_swap],
        coords_node1[should_swap],
    )

    # Distance for control points
    dist = dist_ratio * np.sqrt(np.sum((coords_node1 - coords_node2) ** 2, axis=1))

    # Gradients of line connecting node & perpendicular
    m1 = (coords_node2[:, 1] - coords_node1[:, 1]) / (
        coords_node2[:, 0] - coords_node1[:, 0]
    )
    m2 = -1 / m1

    # Temporary points along the line which connects two nodes
    t1 = dist / np.sqrt(1 + m1 ** 2)
    v1 = np.array([np.ones(n_edges), m1])
    coords_node1_displace = coords_node1 + (v1 * t1).T
    coords_node2_displace = coords_node2 - (v1 * t1).T

    # Control points, same distance but along perpendicular line
    # rnd gives the 'polarity' to determine which side of the line the curve should arc
    t2 = dist / np.sqrt(1 + m2 ** 2)
    v2 = np.array([np.ones(len(edges)), m2])
    coords_node1_ctrl = coords_node1_displace + (rnd * v2 * t2).T
    coords_node2_ctrl = coords_node2_displace + (rnd * v2 * t2).T

    # Combine all these four (x,y) columns into a 'node matrix'
    node_matrix = np.array(
        [coords_node1, coords_node1_ctrl, coords_node2_ctrl, coords_node2]
    )

    nums = np.linspace(0, 2 * np.pi, bezier_precision)

    # Create the Bezier curves and store them in a list

    self_loops = []
    for p in pos_sl.values():
        self_loops.append(np.c_[np.cos(nums), np.sin(nums)] * radius_fraction + p)

    curveplots = []
    for i in range(len(edges)):
        nodes = node_matrix[:, i, :].T
        curveplots.append(
            bezier.Curve(nodes, degree=3)
            .evaluate_multi(np.linspace(0, 1, bezier_precision))
            .T
        )

    # Return an array of these curves
    curves = np.array(curveplots)
    if any(self_loop_mask):
        curves[self_loop_mask, ...] = self_loops

    return curves


@d.dedent
def composition(
    adata: AnnData,
    key: str,
    figsize: Optional[Tuple[float, float]] = None,
    dpi: Optional[float] = None,
    save: Optional[Union[str, Path]] = None,
) -> None:
    """
    Plot a pie chart for categorical annotation.

    .. image:: https://raw.githubusercontent.com/theislab/cellrank/master/resources/images/composition.png
       :width: 400px
       :align: center

    Parameters
    ----------
    %(adata)s
    key
        Key in ``adata.obs`` containing categorical observation.
    %(plotting)s

    Returns
    -------
    %(just_plots)s
    """

    if key not in adata.obs:
        raise KeyError(f"Key `{key!r}` not found in `adata.obs`.")
    if not is_categorical_dtype(adata.obs[key]):
        raise TypeError(f"Observation `adata.obs[{key!r}]` is not categorical.")

    cats = adata.obs[key].cat.categories
    colors = adata.uns.get(f"{key}_colors", None)
    x = [np.sum(adata.obs[key] == cl) for cl in cats]
    cats_frac = x / np.sum(x)

    # plot these fractions in a pie plot
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    ax.pie(x=cats_frac, labels=cats, colors=colors)
    ax.set_title(f"composition by {key}")

    if save is not None:
        save_fig(fig, save)

    fig.show()


def _is_any_gam_mgcv(models: Union[BaseModel, Dict[str, Dict[str, BaseModel]]]) -> bool:
    """
    Return whether any models to be fit are from R's mgcv package.

    Parameters
    ----------
    models
        Model(s) used for fitting.

    Returns
    -------
    bool
        `True` if any of the models is from R's mgcv package, else `False`.
    """

    return isinstance(models, GAMR) or (
        isinstance(models, dict)
        and any(isinstance(m, GAMR) for ms in models.values() for m in ms.values())
    )


def _create_models(
    model: _model_type, obs: Sequence[str], lineages: Sequence[Optional[str]]
) -> Dict[str, Dict[str, BaseModel]]:
    """
    Create models for each gene and lineage.

    Parameters
    ----------
    obs
        Sequence of observations, such as genes.
    lineages
        Sequence of genes.

    Returns
    -------
        The created models.
    """

    def process_lineages(
        obs_name: str, lin_names: Union[BaseModel, Dict[Optional[str], Any]]
    ):
        if isinstance(lin_names, BaseModel):
            # sharing the same models for all lineages
            for lin_name in lineages:
                models[obs_name][lin_name] = lin_names
            return
        lin_rest_model = lin_names.get("*", None)  # do not pop

        for lin_name, mod in lin_names.items():
            if lin_name == "*":
                continue
            models[obs_name][lin_name] = copy(mod)

        if lin_rest_model is not None:
            for lin_name in lineages - set(models[obs_name].keys()):
                models[obs_name][lin_name] = copy(lin_rest_model)
        else:
            raise RuntimeError(_ERROR_INCOMPLETE_SPEC.format(" lineage ", obs_name))

    if isinstance(model, BaseModel):
        return {o: {lin: copy(model) for lin in lineages} for o in obs}

    lineages, obs = (
        set(_unique_order_preserving(lineages)),
        set(_unique_order_preserving(obs)),
    )
    models = defaultdict(dict)

    if isinstance(model, dict):
        obs_rest_model = model.pop("*", None)
        for obs_name, lin_names in model.items():
            process_lineages(obs_name, lin_names)

        if obs_rest_model is not None:
            for obs_name in obs - set(model.keys()):
                process_lineages(obs_name, model.get(obs_name, obs_rest_model))
        elif set(model.keys()) != obs:
            raise RuntimeError(_ERROR_INCOMPLETE_SPEC.format(" ", "genes"))
    else:
        raise TypeError(
            f"Class `{type(model).__name__!r}` must be of type `cellrank.ul.BaseModel` or a dictionary of such models."
        )

    return models


def _fit_gene_trends(
    genes: Sequence[str],
    models: _model_type,
    callbacks: _callback_type,
    lineages: Sequence[Optional[str]],
    time_range: Sequence[Union[float, Tuple[float, float]]],
    queue,
    **kwargs,
) -> Dict[str, Dict[str, Any]]:
    """
    Fit model for given genes and lineages.

    Parameters
    ----------
    genes
        Genes for which to fit the models.
    models
        Gene and lineage specific models.
    callbacks
        Gene and lineage specific prepare callbacks.
    lineages
        Lineages for which to fit the models.
    time_range
        Minimum and maximum pseudotimes.
    queue
        Signalling queue in the parent process/thread used to update the progress bar.
    kwargs
        Keyword arguments for :func:`cellrank.ul.models.BaseModel.prepare`.

    Returns
    -------
        The fitted models, optionally containing the confidence interval.
    """

    res = {}
    conf_int = kwargs.pop("conf_int", False)

    for gene in genes:
        res[gene] = {}
        for ln, tr in zip(lineages, time_range):
            cb = callbacks[gene][ln]
            model = cb(
                models[gene][ln], gene=gene, lineage=ln, time_range=tr, **kwargs
            ).fit()
            model.predict()
            if conf_int:
                model.confidence_interval()

            res[gene][ln] = model

        if queue is not None:
            queue.put(1)

    if queue is not None:
        queue.put(None)

    return res


@d.dedent
def _trends_helper(
    adata: AnnData,
    models: Dict[str, Dict[str, Any]],
    gene: str,
    ln_key: str,
    lineage_names: Optional[Sequence[str]] = None,
    same_plot: bool = False,
    sharey: Union[str, bool] = False,
    show_ylabel: bool = True,
    show_lineage: bool = True,
    show_xticks_and_label: bool = True,
    lineage_cmap=None,
    abs_prob_cmap=cm.viridis,
    gene_as_title: bool = False,
    legend_loc: Optional[str] = "best",
    fig: mpl.figure.Figure = None,
    axes: Union[mpl.axes.Axes, Sequence[mpl.axes.Axes]] = None,
    **kwargs,
) -> None:
    """
    Plot an expression gene for some lineages.

    Parameters
    ----------
    %(adata)s
    %(model)s
    gene
        Name of the gene in `adata.var_names``.
    ln_key
        Key in ``adata.obsm`` where to find the lineages.
    lineage_names
        Names of lineages to plot.
    same_plot
        Whether to plot all lineages in the same plot or separately.
    sharey
        Whether the y-axis is being shared.
    show_ylabel
        Whether to show y-label on the y-axis. Usually, only the first column will contain the y-label.
    show_lineage
        Whether to show the lineage as the title. Usually, only first row will contain the lineage names.
    show_xticks_and_label
        Whether to show x-ticks and x-label. Usually, only the last row will show this.
    lineage_cmap
        Colormap to use when coloring the the lineage. If `None`, use colors from ``adata.uns``.
    abs_prob_cmap:
        Colormap to use when coloring in the absorption probabilities, if they are being plotted.
    gene_as_title
        Whether to use the gene names as titles (with lineage names as well) or on the y-axis.
    legend_loc
        Location of the legend. If `None`, don't show any legend.
    fig
        Figure to use.
    ax
        Ax to use.
    **kwargs
        Keyword arguments for :meth:`cellrank.ul.models.BaseModel.plot`.

    Returns
    -------
    %(just_plots)s
    """

    n_lineages = len(lineage_names)
    if same_plot:
        axes = [axes] * len(lineage_names)

    fig.tight_layout()
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

    if same_plot:
        lineage_colors = (
            lineage_cmap.colors
            if lineage_cmap is not None and hasattr(lineage_cmap, "colors")
            else adata.uns.get(f"{_colors(ln_key)}", cm.Set1.colors)
        )
    else:
        lineage_colors = (
            "black" if not mcolors.is_color_like(lineage_cmap) else lineage_cmap
        )

    for i, (name, ax, perc) in enumerate(zip(lineage_names, axes, percs)):
        if same_plot:
            if gene_as_title:
                title = gene
                ylabel = "expression" if show_ylabel else None
            else:
                title = ""
                ylabel = gene
        else:
            if gene_as_title:
                title = None
                ylabel = "expression" if i == 0 else None
            else:
                title = (
                    (name if name is not None else "no lineage") if show_lineage else ""
                )
                ylabel = gene if i == 0 else None

        models[gene][name].plot(
            ax=ax,
            fig=fig,
            perc=perc,
            show_cbar=False,
            title=title,
            hide_cells=hide_cells or (same_plot and i == n_lineages - 1),
            same_plot=same_plot,
            lineage_color=lineage_colors[i]
            if same_plot and name is not None
            else lineage_colors,
            abs_prob_cmap=abs_prob_cmap,
            ylabel=ylabel,
            **kwargs,
        )
        if sharey in ("row", "all", True) and i == 0:
            plt.setp(ax.get_yticklabels(), visible=True)

        if show_xticks_and_label:
            plt.setp(ax.get_xticklabels(), visible=True)
        else:
            ax.set_xlabel(None)

    if not same_plot and same_perc and show_cbar and not hide_cells:
        vmin = np.min([models[gene][ln].w_all for ln in lineage_names])
        vmax = np.max([models[gene][ln].w_all for ln in lineage_names])
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

        for ax in axes:
            [
                c
                for c in ax.get_children()
                if isinstance(c, mpl.collections.PathCollection)
            ][0].set_norm(norm)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="2.5%", pad=0.1)
        _ = mpl.colorbar.ColorbarBase(
            cax, norm=norm, cmap=abs_prob_cmap, label="absorption probability"
        )

    if same_plot and lineage_names != [None] and legend_loc is not None:
        ax.legend(lineage_names, loc=legend_loc)


def _position_legend(ax: mpl.axes.Axes, legend_loc: str, **kwargs) -> mpl.legend.Legend:
    """
    Position legend in- or outside the figure.

    Parameters
    ----------
    ax
        Ax where to position the legend.
    legend_loc
        Position of legend.
    **kwargs
        Keyword arguments for :func:`matplotlib.pyplot.legend`.

    Returns
    -------
    :class: `matplotlib.legend.Legend`
        The created legend.
    """

    if legend_loc == "center center out":
        raise ValueError(
            "Invalid option: `'center center out'`. Doesn't really make sense, does it?"
        )
    if legend_loc == "best":
        return ax.legend(loc="best", **kwargs)

    tmp, loc = legend_loc.split(" "), ""

    if len(tmp) == 1:
        height, rest = tmp[0], []
        width = "right" if height in ("upper", "top", "center") else "left"
    else:
        height, width, *rest = legend_loc.split(" ")
        if rest:
            if len(rest) != 1:
                raise ValueError(
                    f"Expected only 1 additional modifier ('in' or 'out'), found `{list(rest)}`."
                )
            elif rest[0] not in ("in", "out"):
                raise ValueError(
                    f"Invalid modifier `{rest[0]!r}`. Valid options are: `'in', 'out'`."
                )
            if rest[0] == "in":  # ignore in, it's default
                rest = []

    if height in ("upper", "top"):
        y = 1.55 if width == "center" else 1.025
        loc += "upper"
    elif height == "center":
        y = 0.5
        loc += "center"
    elif height in ("lower", "bottom"):
        y = -0.55 if width == "center" else -0.025
        loc += "lower"
    else:
        raise ValueError(
            f"Invalid legend position on y-axis: `{height!r}`. "
            f"Valid options are: `'upper', 'top', 'center', 'lower', 'bottom'`."
        )

    if width == "left":
        x = -0.05
        loc += " right" if rest else " left"
    elif width == "center":
        x = 0.5
        if height != "center":  # causes to be like top center
            loc += " center"
    elif width == "right":
        x = 1.05
        loc += " left" if rest else " right"
    else:
        raise ValueError(
            f"Invalid legend position on x-axis: `{width!r}`. "
            f"Valid options are: `'left', 'center', 'right'`."
        )

    if rest:
        kwargs["bbox_to_anchor"] = (x, y)

    return ax.legend(loc=loc, **kwargs)


def _maybe_create_dir(dirname: Optional[Union[str, Path]]) -> None:
    if dirname is None:
        return

    from cellrank import logging as logg
    from cellrank import settings

    figdir = settings.figdir

    if figdir is None:
        raise RuntimeError(
            f"Figures directory `cellrank.settings.figdir` is `None`, but `dirname={dirname!r}`."
        )
    if os.path.isabs(dirname):
        if not os.path.isdir(dirname):
            os.makedirs(dirname, exist_ok=True)
            logg.debug(f"Creating directory `{dirname!r}`")
    elif not os.path.isdir(os.path.join(figdir, dirname)):
        os.makedirs(os.path.join(figdir, dirname), exist_ok=True)
        logg.debug(f"Creating directory `{dirname!r}`")


def _get_backend(model, backend: str) -> str:
    return _DEFAULT_BACKEND if _is_any_gam_mgcv(model) else backend


@d.dedent
def _create_callbacks(
    adata: AnnData,
    callback: Optional[Callable],
    obs: Sequence[str],
    lineages: Sequence[Optional[str]],
    perform_sanity_check: bool = True,
) -> Dict[str, Dict[str, Callable]]:
    """
    Create models for each gene and lineage.

    Parameters
    ----------
    %(adata)s
    callback
        Gene and lineage specific prepare callbacks.
    obs
        Sequence of observations, such as genes.
    lineages
        Sequence of genes.
    perform_sanity_check
        Whether to check if all callables have the correct signature. This is done by instantiating
        dummy model and running the function. We're assuming that the callback isn't really a pricey operation.
    Returns
    -------
        The created callbacks.
    """

    def process_lineages(
        obs_name: str, lin_names: Optional[Union[Callable, Dict[Optional[str], Any]]]
    ):
        if lin_names is None:
            lin_names = _default_model_callback

        if callable(lin_names):
            # sharing the same models for all lineages
            for lin_name in lineages:
                callbacks[obs_name][lin_name] = lin_names
            return
        lin_rest_callback = (
            lin_names.get("*", _default_model_callback) or _default_model_callback
        )  # do not pop

        for lin_name, cb in lin_names.items():
            if lin_name == "*":
                continue
            callbacks[obs_name][lin_name] = cb

        if callable(lin_rest_callback):
            for lin_name in lineages - set(callbacks[obs_name].keys()):
                callbacks[obs_name][lin_name] = lin_rest_callback
        else:
            raise TypeError(
                f"Expected the callback for the rest of lineages to be `callable`, "
                f"found `{type(lin_rest_callback).__name__!r}`."
            )

    def maybe_sanity_check(callbacks: Dict[str, Dict[str, Callable]]) -> None:
        if not perform_sanity_check:
            return

        from sklearn.svm import SVR

        logg.debug("Performing callback sanity checks")
        for gene in callbacks.keys():
            for lineage, cb in callbacks[gene].items():
                # create the model here because the callback can search the attribute
                dummy_model = SKLearnModel(adata, model=SVR())
                try:
                    model = cb(dummy_model, gene=gene, lineage=lineage)
                    assert model is dummy_model, (
                        "Creation of new models is not allowed. "
                        "Ensure that callback returns the same model."
                    )
                    assert (
                        model.prepared
                    ), "Model is not prepared. Ensure that callback calls `.prepare()`."
                    assert (
                        model._gene == gene
                    ), f"Callback modified the gene from `{gene!r}` to `{model._gene!r}`."
                    assert (
                        model._gene == gene
                    ), f"Callback modified the lineage from `{lineage!r}` to `{model._lineage!r}`."
                except Exception as e:
                    raise RuntimeError(
                        f"Callback validation failed for gene `{gene!r}` and lineage `{lineage!r}`."
                    ) from e

    if callback is None:
        callback = _default_model_callback

    if callable(callback):
        callbacks = {o: {lin: copy(callback) for lin in lineages} for o in obs}
        maybe_sanity_check(callbacks)
        return callbacks

    lineages, obs = (
        set(_unique_order_preserving(lineages)),
        set(_unique_order_preserving(obs)),
    )
    callbacks = defaultdict(dict)

    if isinstance(callback, dict):
        for obs_name, lin_names in callback.items():
            process_lineages(obs_name, lin_names)

        # can be specified as None
        obs_rest_callback = (
            callback.pop("*", _default_model_callback) or _default_model_callback
        )

        if callable(obs_rest_callback):
            for obs_name in obs - set(callback.keys()):
                process_lineages(obs_name, callback.get(obs_name, obs_rest_callback))
        else:
            raise TypeError(
                f"Expected the callback for the rest of genes to be `callable`, "
                f"found `{type(obs_rest_callback).__name__!r}`."
            )
    else:
        raise TypeError(
            f"Class `{type(callback).__name__!r}` must be callable` or a dictionary of callables."
        )

    maybe_sanity_check(callbacks)

    return callbacks


def _default_model_callback(model: BaseModel, **kwargs) -> BaseModel:
    return model.prepare(**kwargs)
