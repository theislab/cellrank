"""Graph pl module."""

from copy import deepcopy
from types import MappingProxyType
from typing import Dict, Tuple, Union, TypeVar, Callable, Optional, Sequence
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import issparse, spmatrix
from pandas.api.types import is_categorical_dtype

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap
from matplotlib.patches import ArrowStyle, FancyArrowPatch
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable

from cellrank import logging as logg
from cellrank.ul._docs import d
from cellrank.tl._utils import save_fig
from cellrank.ul._utils import _read_graph_data
from cellrank.tl._constants import _colors

AnnData = TypeVar("AnnData")


KEYLOCS = str
KEYS = str
_msg_shown = False


@d.dedent
def graph(
    data: Union[AnnData, np.ndarray, spmatrix],
    graph_key: Optional[str] = None,
    ixs: Optional[Union[range, np.array]] = None,
    layout: Union[str, Dict, Callable] = "umap",
    keys: Sequence[KEYS] = ("incoming",),
    keylocs: Union[KEYLOCS, Sequence[KEYLOCS]] = "uns",
    node_size: float = 400,
    labels: Optional[Union[Sequence[str], Sequence[Sequence[str]]]] = None,
    top_n_edges: Optional[Union[int, Tuple[int, bool, str]]] = None,
    self_loops: bool = True,
    self_loop_radius_frac: Optional[float] = None,
    filter_edges: Optional[Tuple[float, float]] = None,
    edge_reductions: Union[Callable, Sequence[Callable]] = np.sum,
    edge_reductions_restrict_to_ixs: Optional[Union[range, np.ndarray]] = None,
    edge_weight_scale: float = 10,
    edge_width_limit: Optional[float] = None,
    edge_alpha: float = 1.0,
    edge_normalize: bool = False,
    edge_use_curved: bool = True,
    arrows: bool = True,
    font_size: int = 12,
    font_color: str = "black",
    color_nodes: bool = True,
    cat_cmap: ListedColormap = cm.Set3,
    cont_cmap: ListedColormap = cm.viridis,
    legend_loc: Optional[str] = "best",
    title: Optional[Union[str, Sequence[Optional[str]]]] = None,
    figsize: Optional[Tuple[float, float]] = None,
    dpi: Optional[int] = None,
    save: Optional[Union[str, Path]] = None,
    layout_kwargs: Dict = MappingProxyType({}),
) -> None:
    """
    Plot a graph, visualizing incoming and outgoing edges or self-transitions.

    This is a utility function to look in more detail at the transition matrix in areas of interest, e.g. around an
    endpoint of development. This function is meant to visualise a small subset of nodes (~100-500) and the most likely
    transitions between them. Note that limiting edges visualized using ``top_n_edges`` will speed things up,
    as well as reduce the visual clutter.

    Parameters
    ----------
    data
        The graph data to be plotted.
    graph_key
        Key in ``adata.obsp`` or ``adata.uns`` where the graph is stored. Only used
        when ``data`` is :class:`~anndata.Anndata` object.
    ixs
        Subset of indices of the graph to visualize.
    layout
        Layout to use for graph drawing.

        - If :class:`str`, search for embedding in ``adata.obsm['X_{layout}']``.
          Use ``layout_kwargs={'components': [0, 1]}`` to select components.
        - If :class:`dict`, keys should be values in interval ``[0, len(ixs))``
          and values `(x, y)` pairs corresponding to node positions.
    keys
        Keys in ``adata.obs``, ``adata.obsm`` or ``adata.obsp`` to color the nodes.

        - If `'incoming'`, `'outgoing'` or `'self_loops'`, visualize reduction (see ``edge_reductions``)
          for each node based on incoming or outgoing edges, respectively.
    keylocs
        Locations of ``keys``. Can be any attribute of ``data`` if it's :class:`anndata.AnnData` object.
    node_size
        Size of the nodes.
    labels
        Labels of the nodes.
    top_n_edges
        Either top N outgoing edges in descending order or a tuple
        ``(top_n_edges, in_ascending_order, {'incoming', 'outgoing'})``. If `None`, show all edges.
    self_loops
        Whether visualize self transitions and also to consider them in ``top_n_edges``.
    self_loop_radius_frac
        Fraction of a unit circle to visualize self transitions. If `None`, use ``node_size / 1000``.
    filter_edges
        Whether to remove all edges not in `[min, max]` interval.
    edge_reductions
        Aggregation function to use when coloring nodes by edge weights.
    edge_reductions_restrict_to_ixs
        Whether to use the full graph when calculating the ``edge_reductions`` or just use the nodes
        marked by the ``ixs`` and this parameter. If `None`, it's the same as ``ixs``.
    edge_weight_scale
        Number by which to scale the width of the edges. Useful when the weights are small.
    edge_width_limit
        Upper bound for the width of the edges. Useful when weights are unevenly distributed.
    edge_alpha
        Alpha channel value for edges and arrows.
    edge_normalize
        If `True`, normalize edges to `[0, 1]` interval prior to applying any scaling or truncation.
    edge_use_curved
        If `True`, use curved edges. This can improve visualization at a small performance cost.
    arrows
        Whether to show the arrows. Setting this to `False` may dramatically speed things up.
    font_size
        Font size for node labels.
    font_color
        Label color of the nodes.
    color_nodes
        Whether to color the nodes
    cat_cmap
        Categorical colormap used when ``keys`` contain categorical variables.
    cont_cmap
        Continuous colormap used when ``keys`` contain continuous variables.
    legend_loc
        Location of the legend.
    title
        Title of the figure(s), one for each ``key``.
    %(plotting)s
    layout_kwargs
        Additional kwargs for ``layout``.

    Returns
    -------
    %(just_plots)s
    """

    from anndata import AnnData as _AnnData

    import networkx as nx

    def plot_arrows(curves, G, pos, ax, edge_weight_scale):
        for line, (edge, val) in zip(curves, G.edges.items()):
            if edge[0] == edge[1]:
                continue

            mask = (~np.isnan(line)).all(axis=1)
            line = line[mask, :]
            if not len(line):  # can be all NaNs
                continue

            line = line.reshape((-1, 2))
            X, Y = line[:, 0], line[:, 1]

            node_start = pos[edge[0]]
            # reverse
            if np.where(np.isclose(node_start - line, [0, 0]).all(axis=1))[0][0]:
                X, Y = X[::-1], Y[::-1]

            mid = len(X) // 2
            posA, posB = zip(X[mid : mid + 2], Y[mid : mid + 2])

            arrow = FancyArrowPatch(
                posA=posA,
                posB=posB,
                # we clip because too small values
                # cause it to crash
                arrowstyle=ArrowStyle.CurveFilledB(
                    head_length=np.clip(
                        val["weight"] * edge_weight_scale * 4,
                        _min_edge_weight,
                        edge_width_limit,
                    ),
                    head_width=np.clip(
                        val["weight"] * edge_weight_scale * 2,
                        _min_edge_weight,
                        edge_width_limit,
                    ),
                ),
                color="k",
                zorder=float("inf"),
                alpha=edge_alpha,
                linewidth=0,
            )
            ax.add_artist(arrow)

    def normalize_weights():
        weights = np.array([v["weight"] for v in G.edges.values()])
        minn = np.min(weights)
        weights = (weights - minn) / (np.max(weights) - minn)
        for v, w in zip(G.edges.values(), weights):
            v["weight"] = w

    def remove_top_n_edges():
        if top_n_edges is None:
            return

        if isinstance(top_n_edges, (tuple, list)):
            to_keep, ascending, group_by = top_n_edges
        else:
            to_keep, ascending, group_by = top_n_edges, False, "out"

        if group_by not in ("incoming", "outgoing"):
            raise ValueError(
                "Argument `groupby` in `top_n_edges` must be either `'incoming`' or `'outgoing'`."
            )

        source, target = zip(*G.edges)
        weights = [v["weight"] for v in G.edges.values()]
        tmp = pd.DataFrame({"outgoing": source, "incoming": target, "w": weights})

        if not self_loops:
            # remove self loops
            tmp = tmp[tmp["incoming"] != tmp["outgoing"]]

        to_keep = set(
            map(
                tuple,
                tmp.groupby(group_by)
                .apply(
                    lambda g: g.sort_values("w", ascending=ascending).take(
                        range(min(to_keep, len(g)))
                    )
                )[["outgoing", "incoming"]]
                .values,
            )
        )

        for e in list(G.edges):
            if e not in to_keep:
                G.remove_edge(*e)

    def remove_low_weight_edges():
        if filter_edges is None or filter_edges == (None, None):
            return

        minn, maxx = filter_edges
        minn = minn if minn is not None else -np.inf
        maxx = maxx if maxx is not None else np.inf

        for e, attr in list(G.edges.items()):
            if attr["weight"] < minn or attr["weight"] > maxx:
                G.remove_edge(*e)

    _min_edge_weight = 0.00001

    if edge_width_limit is None:
        logg.debug("Not limiting width of edges")
        edge_width_limit = float("inf")

    if self_loop_radius_frac is None:
        self_loop_radius_frac = (
            node_size / 2000 if node_size >= 200 else node_size / 1000
        )
        logg.debug(f"Setting self loop radius fraction to `{self_loop_radius_frac}`")

    if not isinstance(keys, (tuple, list)):
        keys = [keys]

    if not isinstance(keylocs, (tuple, list)):
        keylocs = [keylocs] * len(keys)
    elif len(keylocs) == 1:
        keylocs = keylocs * 3
    elif all(map(lambda k: k in ("incoming", "outgoing", "self_loops"), keys)):
        # don't care about keylocs since they are irrelevant
        logg.debug("Ignoring key locations")
        keylocs = [None] * len(keys)

    if not isinstance(edge_reductions, (tuple, list)):
        edge_reductions = [edge_reductions] * len(keys)
    if not all(map(callable, edge_reductions)):
        raise ValueError("Not all `edge_reductions` functions are callable.")

    if not isinstance(labels, (tuple, list)):
        labels = [labels] * len(keys)
    elif not len(labels):
        labels = [None] * len(keys)
    elif not isinstance(labels[0], (tuple, list)):
        labels = [labels] * len(keys)

    if len(keys) != len(labels):
        raise ValueError(
            f"`Keys` and `labels` must be of the same shape, found `{len(keys)}` and `{len(labels)}`."
        )
    if title is None or isinstance(title, str):
        title = [title] * len(keys)
    if len(title) != len(keys):
        raise ValueError(
            f"`Titles` and `keys` must be of the same shape, found `{len(title)}` and `{len(keys)}`."
        )

    if isinstance(data, _AnnData):
        if graph_key is None:
            raise ValueError(
                "Argument `graph_key` cannot be `None` when `data` is `anndata.Anndata` object."
            )
        gdata = _read_graph_data(data, graph_key)
    elif isinstance(data, (np.ndarray, spmatrix)):
        gdata = data
    else:
        raise TypeError(
            f"Expected argument `data` to be one of `anndata.AnnData`, `numpy.ndarray`, `scipy.sparse.spmatrix`, "
            f"found `{type(data).__name__!r}`."
        )
    is_sparse = issparse(gdata)

    gdata_full = gdata
    if ixs is not None:
        gdata = gdata[ixs, :][:, ixs]
    else:
        ixs = list(range(gdata.shape[0]))

    if edge_reductions_restrict_to_ixs is None:
        edge_reductions_restrict_to_ixs = ixs

    start = logg.info("Creating graph")
    G = (
        nx.from_scipy_sparse_matrix(gdata, create_using=nx.DiGraph)
        if is_sparse
        else nx.from_numpy_array(gdata, create_using=nx.DiGraph)
    )

    remove_low_weight_edges()
    remove_top_n_edges()
    if edge_normalize:
        normalize_weights()
    logg.info("    Finish", time=start)

    # do NOT recreate the graph, for the edge reductions
    # gdata = nx.to_numpy_array(G)

    if figsize is None:
        figsize = (12, 8 * len(keys))

    fig, axes = plt.subplots(nrows=len(keys), ncols=1, figsize=figsize, dpi=dpi)
    if not isinstance(axes, np.ndarray):
        axes = np.array([axes])
    axes = np.ravel(axes)

    if isinstance(layout, str):
        if f"X_{layout}" not in data.obsm:
            raise KeyError(f"Unable to find embedding `'X_{layout}'` in `adata.obsm`.")
        components = layout_kwargs.get("components", [0, 1])
        if len(components) != 2:
            raise ValueError(
                f"Components in `layout_kwargs` must be of length `2`, found `{len(components)}`."
            )
        emb = data.obsm[f"X_{layout}"][:, components]
        pos = {i: emb[ix, :] for i, ix in enumerate(ixs)}
        logg.info(f"Embedding graph using `{layout!r}` layout")
    elif isinstance(layout, dict):
        rng = range(len(ixs))
        for k, v in layout.items():
            if k not in rng:
                raise ValueError(
                    f"Key in `layout` must be in `range(len(ixs))`, found `{k}`."
                )
            if len(v) != 2:
                raise ValueError(
                    f"Value in `layout` must be a `tuple` or a `list` of length 2, found `{len(v)}`."
                )
        pos = layout
        logg.debug("Using precomputed layout")
    elif callable(layout):
        start = logg.info(f"Embedding graph using `{layout.__name__!r}` layout")
        pos = layout(G, **layout_kwargs)
        logg.info("    Finish", time=start)
    else:
        raise TypeError(
            f"Argument `layout` must be either a `string`, "
            f"a `dict` or a `callable`, found `{type(layout)}`."
        )

    curves, lc = None, None
    if edge_use_curved:
        try:
            from cellrank.pl._utils import _curved_edges

            logg.debug("Creating curved edges")
            curves = _curved_edges(G, pos, self_loop_radius_frac, polarity="directed")
            lc = LineCollection(
                curves,
                colors="black",
                linewidths=np.clip(
                    np.ravel([v["weight"] for v in G.edges.values()])
                    * edge_weight_scale,
                    0,
                    edge_width_limit,
                ),
                alpha=edge_alpha,
            )
        except ImportError as e:
            global _msg_shown
            if not _msg_shown:
                print(
                    str(e)[:-1],
                    "in order to use curved edges or specify `edge_use_curved=False`.",
                )
                _msg_shown = True

    for ax, keyloc, title, key, labs, er in zip(
        axes, keylocs, title, keys, labels, edge_reductions
    ):
        label_col = {}  # dummy value

        if key in ("incoming", "outgoing", "self_loops"):
            if key in ("incoming", "outgoing"):
                axis = int(key == "outgoing")
                tmp_data = (
                    gdata_full[ixs, :][:, edge_reductions_restrict_to_ixs]
                    if axis
                    else gdata_full[edge_reductions_restrict_to_ixs, :][:, ixs]
                )
                vals = er(tmp_data, axis=axis)
            else:
                vals = gdata.diagonal() if is_sparse else np.diag(gdata)

            if issparse(vals):
                vals = vals.toarray()
            vals = np.squeeze(np.array(vals))

            assert len(pos) == len(vals), (
                f"Expected reduction values to be of length `{len(pos)}`, "
                f"found `{len(vals)}`."
            )

            node_v = dict(zip(pos.keys(), vals))
        else:
            label_col = getattr(data, keyloc)
            if key in label_col:
                node_v = dict(zip(pos.keys(), label_col[key]))
            else:
                raise RuntimeError(f"Key `{key!r}` not found in `adata.{keyloc}`.")

        if labs is not None:
            if len(labs) != len(pos):
                raise RuntimeError(
                    f"Number of labels ({len(labels)}) and nodes ({len(pos)}) mismatch."
                )
            nx.draw_networkx_labels(
                G,
                pos,
                labels=labs if isinstance(labs, dict) else dict(zip(pos.keys(), labs)),
                ax=ax,
                font_color=font_color,
                font_size=font_size,
            )

        if lc is not None and curves is not None:
            ax.add_collection(deepcopy(lc))  # copying necessary
            if arrows:
                plot_arrows(curves, G, pos, ax, edge_weight_scale)
        else:
            nx.draw_networkx_edges(
                G,
                pos,
                width=[
                    np.clip(
                        v["weight"] * edge_weight_scale,
                        _min_edge_weight,
                        edge_width_limit,
                    )
                    for _, v in G.edges.items()
                ],
                alpha=edge_alpha,
                edge_color="black",
                arrows=True,
                arrowstyle="-|>",
            )

        if key in label_col and is_categorical_dtype(label_col[key]):
            values = label_col[key]
            if keyloc in ("obs", "obsm"):
                values = values[ixs]
            categories = values.cat.categories
            color_key = _colors(key)
            if color_key in data.uns:
                mapper = dict(zip(categories, data.uns[color_key]))
            else:
                mapper = dict(
                    zip(categories, map(cat_cmap.get, range(len(categories))))
                )

            colors = []
            seen = set()

            for v in values:
                colors.append(mapper[v])
                seen.add(v)

            nodes_kwargs = dict(cmap=cat_cmap, node_color=colors)  # noqa
            if legend_loc not in ("none", None):
                x, y = pos[0]
                for label in sorted(seen):
                    ax.plot([x], [y], label=label, color=mapper[label])
                ax.legend(loc=legend_loc)
        else:
            values = list(node_v.values())
            vmin, vmax = np.min(values), np.max(values)
            nodes_kwargs = dict(  # noqa
                cmap=cont_cmap, node_color=values, vmin=vmin, vmax=vmax
            )

            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="1.5%", pad=0.05)
            _ = mpl.colorbar.ColorbarBase(
                cax,
                cmap=cont_cmap,
                norm=mpl.colors.Normalize(vmin=vmin, vmax=vmax),
                ticks=np.linspace(vmin, vmax, 10),
            )

        if color_nodes is False:
            nodes_kwargs = {}

        nx.draw_networkx_nodes(G, pos, node_size=node_size, ax=ax, **nodes_kwargs)

        ax.set_title(key if title is None else title)
        ax.axis("off")

    if save is not None:
        save_fig(fig, save)

    fig.show()
