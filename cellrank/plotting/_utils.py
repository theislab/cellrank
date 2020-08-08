# -*- coding: utf-8 -*-
"""Utility functions for CellRank plotting."""

import os
from copy import copy
from typing import (
    Any,
    Dict,
    Tuple,
    Union,
    Mapping,
    TypeVar,
    Iterable,
    Optional,
    Sequence,
)
from pathlib import Path
from collections import defaultdict

import numpy as np
from scipy.sparse import diags
from pandas.core.dtypes.common import is_categorical_dtype

import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from cellrank.utils._docs import d
from cellrank.tools._utils import save_fig, _unique_order_preserving
from cellrank.utils.models import GAMR, BaseModel
from cellrank.tools.kernels import PrecomputedKernel
from cellrank.tools._constants import _DEFAULT_BACKEND, _colors
from cellrank.tools.estimators._cflare import CFLARE

AnnData = TypeVar("AnnData")


_ERROR_INCOMPLETE_SPEC = (
    "No options were specified for{}`{!r}`. "
    "Consider specifying a fallback model using '*'."
)
_time_range_type = Optional[Union[float, Tuple[Optional[float], Optional[float]]]]
_model_type = Union[BaseModel, Mapping[str, Mapping[str, BaseModel]]]


@d.dedent
def lineages(
    adata: AnnData,
    lineages: Optional[Union[str, Iterable[str]]] = None,
    backward: bool = False,
    cluster_key: Optional[str] = None,
    mode: str = "embedding",
    time_key: str = "latent_time",
    cmap: Union[str, mpl.colors.ListedColormap] = cm.viridis,
    **kwargs,
) -> None:
    """
    Plot lineages that were uncovered using :func:`cellrank.tl.lineages`.

    For each lineage, we show all cells in an embedding (default is UMAP, but can be any) and color them by their
    probability of belonging to this lineage. For cells that are already committed, this probability will be one for
    their respective lineage and zero otherwise. For naive cells, these probabilities will be more balanced, reflecting
    the fact that naive cells have the potential to develop towards multiple endpoints.

    .. image:: https://raw.githubusercontent.com/theislab/cellrank/master/resources/images/lineages.png
       :width: 400px
       :align: center

    Parameters
    ----------
    %s(adata)s
    lineages
        Plot only these lineages. If `None`, plot all lineages.
    %(backward)s
    cluster_key
        If given, plot cluster annotations left of the lineage probabilities.
    mode
        Can be either `'embedding'` or `'time'`:

            - `'embedding'` - plot the embedding while coloring in the absorption probabilities.
            - `'time'` - plot the pseudotime on x-axis and the absorption probabilities on y-axis.
    time_key
        Key from `adata.obs` to use as a pseudotime ordering of the cells.
    cmap
        Colormap to use.
    **kwargs
        Keyword arguments for :func:`scvelo.pl.scatter`.

    Returns
    -------
    %s(just_plots)s
    """

    # create a dummy kernel object
    dummy_transition_matrix = diags(np.ones(adata.n_obs))
    k = PrecomputedKernel(dummy_transition_matrix, adata, backward=backward)

    # use this to initialize an MC object
    mc = CFLARE(k, read_from_adata=True, write_to_adata=False)

    # plot using the MC object
    mc.plot_absorption_probabilities(
        lineages=lineages,
        cluster_key=cluster_key,
        mode=mode,
        time_key=time_key,
        cmap=cmap,
        **kwargs,
    )


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
        Polarity of curves, one of `'random'`, `'directed' or `'fixed'`.`
        If using `'random'`, incoming and outgoing edges may overlap.

    Returns
    -------
    :class:`np.ndarray`
        Array of shape (n_edges, :paramref:`bezier_precision`, 2) containing the curved edges.
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
            f"Valid options are: `'random', 'fixed' or 'fixed'`."
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
        Key in :paramref:`adata` `.obs` containing categorical observation.
    %s(plotting)s

    Returns
    -------
    %s(just_plots)s
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
        `True` if any of the models is from R's mgcv package, else `False`.
    """

    return isinstance(models, GAMR) or (
        isinstance(models, dict)
        and any(isinstance(m, GAMR) for ms in models.values() for m in ms.values())
    )


def _create_models(
    model: _model_type, obs: Sequence[str], lineages: Sequence[str]
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

    def process_lineages(obs_name: str, lin_names: Union[BaseModel, Dict[str, Any]]):
        if isinstance(lin_names, BaseModel):
            for lin_name in lineages:
                models[obs_name][lin_name] = lin_names
            return
        lin_rest_model = lin_names.get("*", None)  # do not pop

        for lin_name, mod in lin_names.items():
            if lin_name == "*":
                continue
            models[obs_name][lin_name] = copy(mod)

        if lin_rest_model is not None:
            for lin_name in set(lineages) - set(models[obs_name].keys()):
                models[obs_name][lin_name] = copy(lin_rest_model)
        elif set(models[obs_name].keys()) != lineages:
            raise RuntimeError(_ERROR_INCOMPLETE_SPEC.format(" lineage ", obs_name))

    if isinstance(model, BaseModel):
        return {o: {lin: copy(model) for lin in lineages} for o in obs}

    lineages, obs = _unique_order_preserving(lineages), _unique_order_preserving(obs)
    models = defaultdict(dict)

    if isinstance(model, BaseModel):
        model = {"*": {"*": model}}

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
    lineage_names: Sequence[Optional[str]],
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
    lineage_names
        Lineages for which to fit the models.
    time_range
        Minimum and maximum pseudotimes.
    queue
        Signalling queue in the parent process/thread used to update the progress bar.
    kwargs
        Keyword arguments for :func:`cellrank.utils.models.BaseModel.prepare`.

    Returns
    -------
        The fitted models, optionally containing the confidence interval.
    """

    res = {}
    models = kwargs.pop("models")
    conf_int = kwargs.pop("conf_int", False)

    for gene in genes:
        res[gene] = {}
        for ln, tr in zip(lineage_names, time_range):
            model = models[gene][ln].prepare(gene, ln, time_range=tr, **kwargs).fit()
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
    sharey: bool = True,
    cmap=None,
    fig: mpl.figure.Figure = None,
    ax: mpl.axes.Axes = None,
    save: Optional[Union[str, Path]] = None,
    **kwargs,
) -> None:
    """
    Plot an expression gene for some lineages.

    Parameters
    ----------
    %(adata)s
    %(model)s
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
        Keyword arguments for :meth:`cellrank.ul.models.BaseModel.plot`.

    Returns
    -------
    %(just_plots)s
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
        title = name if name is not None else "no lineage"
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
            ylabel=gene if not same_plot or name is None else "expression",
            **kwargs,
        )

    if same_plot and lineage_names != [None]:
        ax.set_title(gene)
        ax.legend()

    if save is not None:
        save_fig(fig, save)


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

    from cellrank import settings, logging as logg

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
