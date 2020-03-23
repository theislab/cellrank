# -*- coding: utf-8 -*-
from typing import Sequence, Dict, Optional, Tuple, Any, Union, Iterable
from collections import defaultdict
from copy import copy

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from scipy.sparse import csr_matrix

from anndata import AnnData
from pandas.core.dtypes.common import is_categorical_dtype

from cellrank.tools._markov_chain import MarkovChain
from cellrank.plotting._constants import _model_type
from cellrank.utils.models import Model, GamMGCVModel
from cellrank.tools._utils import save_fig
from cellrank.tools.kernels import VelocityKernel

_ERROR_INCOMPLETE_SPEC = (
    "No options were specified for{}`{!r}`. "
    "Consider specifying a fallback model using '*'."
)


def lineages(
    adata: AnnData,
    lineages: Optional[Union[str, Iterable[str]]] = None,
    final: bool = True,
    cluster_key: Optional[str] = None,
    mode: str = "embedding",
    time_key: str = "latent_time",
    color_map: str = "viridis",
    **kwargs,
) -> None:
    """
    Plot lineages that were uncovered using :func:`cellrank.tl.lineages`.

    For each lineage, we show all cells in an embedding (default is UMAP but can be any) and color them by their
    probability of belonging to this lineage. For cells that are already committed, this probability will be one for
    their respective lineage and zero otherwise. For naive cells, these probabilities will be more balanced, reflecting
    the fact that naive cells have the potential to develop towards multiple endpoints.

    .. image:: https://raw.githubusercontent.com/theislab/cellrank/master/resources/images/lineages.png
       :width: 400px
       :align: center

    Params
    ------

    adata : :class:`adata.AnnData`
        Annotated data object.
    lineages
        Only show these lineages. If `None`, plot all lineages.
    final
        Whether to consider cells going to final states or vice versa.
    cluster_key
        If given, plot cluster annotations left of the lineage probabilities.
    mode
        Can be either `'embedding'` or `'time'`.

        - If `'embedding'`, plots the embedding while coloring in the absorption probabilities.
        - If `'time'`, plots the pseudotime on x-axis and the absorption probabilities on y-axis.
    time_key
        Key from `adata.obs` to use as a pseudotime ordering of the cells.
    color_map
        Colormap to use.
    kwargs
        Keyword arguments for :func:`scvelo.pl.scatter`.

    Returns
    -------
    None
        Just plots the lineage probabilities.
    """

    adata_dummy = adata.copy()

    # create a dummy kernel object
    vk = VelocityKernel(adata_dummy, backward=not final)
    vk.transition_matrix = csr_matrix((adata_dummy.n_obs, adata_dummy.n_obs))

    # use this to initialize an MC object
    mc = MarkovChain(vk)

    # plot using the MC object
    mc.plot_lin_probs(
        lineages=lineages,
        cluster_key=cluster_key,
        mode=mode,
        time_key=time_key,
        color_map=color_map,
        **kwargs,
    )


def curved_edges(
    G: nx.Graph,
    pos,
    radius_fraction: float,
    dist_ratio: float = 0.2,
    bezier_precision: int = 20,
    polarity: str = "directed",
) -> np.ndarray:
    """
    Create curved edges from a graph. Modified from: https://github.com/beyondbeneath/bezier-curved-edges-networkx.

    Params
    ------
    G: nx.Graph
        Networkx graph.
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
    except ImportError:
        raise ImportError("Please install `bezier` as `pip install bezier`.")

    # Get nodes into np array
    edges = np.array(G.edges())
    l = edges.shape[0]

    self_loop_mask = edges[:, 0] == edges[:, 1]
    pos_sl = {edge[0]: pos[edge[0]] for edge in edges[self_loop_mask, ...]}

    if polarity == "random":
        # Random polarity of curve
        rnd = np.where(np.random.randint(2, size=l) == 0, -1, 1)
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
        (coords_node2[:, 0] - coords_node1[:, 0])
    )
    m2 = -1 / m1

    # Temporary points along the line which connects two nodes
    t1 = dist / np.sqrt(1 + m1 ** 2)
    v1 = np.array([np.ones(l), m1])
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


def composition(
    adata: AnnData,
    key,
    figsize: Optional[Tuple[float, float]] = None,
    dpi: Optional[float] = None,
    save: Optional[str] = None,
) -> None:
    """
    Utility function to plot pie chart for categorical annotation

    .. image:: https://raw.githubusercontent.com/theislab/cellrank/master/resources/images/composition.png
       :width: 400px
       :align: center

    Params
    ------
    adata
        Annotated data object.
    key
        Key in :paramref:`adata` `.obs` containing categorical observation.
    figsize
        Size of the figure.
    dpi
        Dots per inch.
    save
        Filename where to save the plots.
        If `None`, just shows the plot.

    Returns
    -------
    None
        Nothing, just plots the similarity matrix.
        Optionally saves the figure based on :paramref:`save`.
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
    ax.set_title(f"Composition by {key}")

    if save is not None:
        save_fig(fig, save)

    fig.show()


def _is_any_gam_mgcv(models: Dict[str, Dict[str, Model]]) -> bool:
    """
    Return whether any models to be fit are from R's mgcv package.

    Params
    ------
    models
        Models used for fitting.

    Returns
    -------
        `True` if any of the models is from R's mgcv package, else `False`.
    """

    return any(
        isinstance(m, GamMGCVModel) for ms in models.values() for m in ms.values()
    )


def _create_models(
    model: _model_type, obs: Sequence[str], lineages: Sequence[str]
) -> Dict[str, Dict[str, Model]]:
    """
    Create models for each gene and lineage.

    Params
    ------
    obs
        Sequence of observations, such as genes.
    lineages
        Sequence of genes.

    Returns
    -------
        The created models.
    """

    def process_lineages(obs_name: str, lin_names: Union[Model, Dict[str, Any]]):
        if isinstance(lin_names, Model):
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
        elif set(models[obs_name].keys()) != lineages:
            raise RuntimeError(_ERROR_INCOMPLETE_SPEC.format(" lineage ", obs_name))

    if isinstance(model, Model):
        return {o: {l: copy(model) for l in lineages} for o in obs}

    lineages, obs = set(lineages), set(obs)
    models = defaultdict(dict)

    if isinstance(model, Model):
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
        raise ValueError(
            "Model must be of type `cellrank.ul.Model` or a dictionary of such models."
        )

    return models


def _fit(
    genes: Sequence[str],
    lineage_names: Sequence[Optional[str]],
    start_clusters: Sequence[Optional[str]],
    end_clusters: Sequence[Optional[str]],
    queue,
    **kwargs,
) -> Dict[str, Dict[str, Any]]:
    """
    Fit model for given genes and lineages.

    Params
    ------
    genes
        Genes for which to fit the models.
    lineage_names
        Lineages for which to fit the models.
    start_clusters
        Start clusters for given :paramref:`lineage_names`.
    end_clusters
        End clusters for given :paramref:`lineage_names`.
    queue
        Signalling queue in the parent process/thread used to update the progress bar.
    kwargs
        Keyword arguments for :func:`cellrank.utils.models.Model.prepare`.

    Returns
    -------
        The fitted models, optionally containing the confidence interval.
    """

    res = {}
    models = kwargs.pop("models")
    conf_int = kwargs.pop("conf_int", False)

    for gene in genes:
        res[gene] = {}
        for ln, sc, ec in zip(lineage_names, start_clusters, end_clusters):
            model = (
                models[gene][ln]
                .prepare(gene, ln, start_cluster=sc, end_cluster=ec, **kwargs)
                .fit()
            )
            model.predict()
            if conf_int:
                model.confidence_interval()

            res[gene][ln] = model
        queue.put(1)
    queue.put(None)

    return res
