import enum
import pathlib
import types
from typing import Any, Callable, Literal, Mapping, Optional, Sequence, Tuple, Union

import scvelo as scv

import numpy as np
import pandas as pd
from sklearn.metrics import pairwise_distances

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, LogNorm

from anndata import AnnData
from scanpy._utils import deprecated_arg_names

from cellrank import logging as logg
from cellrank._utils import Lineage
from cellrank._utils._docs import d
from cellrank._utils._enum import ModeEnum
from cellrank._utils._key import Key
from cellrank._utils._lineage import PrimingDegree
from cellrank._utils._utils import _check_collection, _unique_order_preserving, save_fig
from cellrank.pl._utils import _held_karp

__all__ = ["circular_projection"]


class LineageOrder(ModeEnum):
    DEFAULT = enum.auto()
    OPTIMAL = enum.auto()


class LabelRot(ModeEnum):
    DEFAULT = enum.auto()
    BEST = enum.auto()


Metric_T = Union[str, Callable, np.ndarray, pd.DataFrame]
_N = 200


def _get_distances(data: Union[np.ndarray, Lineage], metric: Metric_T) -> np.ndarray:
    if isinstance(data, Lineage):
        data = data.X

    if isinstance(metric, str) or callable(metric):
        metric = pairwise_distances(data.T, metric=metric)
    elif isinstance(metric, (pd.DataFrame, np.ndarray)):
        shape = (data.shape[1], data.shape[1])
        if metric.shape != shape:
            raise ValueError(
                f"Expected an `numpy.array` or `pandas.DataFrame` of shape `{shape}`, found `{metric.shape}`."
            )
    else:
        raise TypeError(
            f"Expected either metric defined by `str`, `callable` or a pairwise distance matrix of type"
            f" `numpy.ndarray` or `pandas.DataFrame`, found `{type(metric).__name__}`."
        )

    return np.asarray(metric, dtype=np.float64)


def _get_optimal_order(data: Lineage, metric: Metric_T) -> Tuple[float, np.ndarray]:
    """Solve the TSP using dynamic programming."""
    return _held_karp(_get_distances(data, metric))


@d.dedent
@deprecated_arg_names({"labeldistance": "label_distance", "labelrot": "label_rot"})
def circular_projection(
    adata: AnnData,
    keys: Union[str, Sequence[str]],
    backward: bool = False,
    lineages: Optional[Union[str, Sequence[str]]] = None,
    early_cells: Optional[Union[Mapping[str, Sequence[str]], Sequence[str]]] = None,
    lineage_order: Optional[Literal["default", "optimal"]] = None,
    metric: Union[str, Callable, np.ndarray, pd.DataFrame] = "correlation",
    normalize_by_mean: bool = True,
    ncols: int = 4,
    space: float = 0.25,
    use_raw: bool = False,
    text_kwargs: Mapping[str, Any] = types.MappingProxyType({}),
    label_distance: float = 1.25,
    label_rot: Union[Literal["default", "best"], float] = "best",
    show_edges: bool = True,
    key_added: Optional[str] = None,
    figsize: Optional[Tuple[float, float]] = None,
    dpi: Optional[int] = None,
    save: Optional[Union[str, pathlib.Path]] = None,
    **kwargs: Any,
) -> None:
    r"""Visualize fate probabilities in a circular embedding :cite:`velten:17,jaitin:14`.

    We arrange all computed terminal states evenly spaced around the unit circle and place
    cells inside the unit circle in a way that reflects their fate probabilities. In other words,
    the more certain we are that a cell will transition towards a given terminal state, the closer
    we place it to that terminal state. Uncommitted cells thus reside in the middle of the circle.

    Parameters
    ----------
    %(adata)s
    keys
        Keys in :attr:`~anndata.AnnData.obs` or :attr:`~anndata.AnnData.var_names`. Additional keys are:

        - ``'kl_divergence'`` - as in :cite:`velten:17`, computes KL-divergence between the fate probabilities
          of a cell and the average fate probabilities. See ``early_cells`` for more information.
        - ``'entropy'`` - as in :cite:`setty:19`, computes entropy over a cells fate probabilities.
    %(backward)s
    lineages
        Lineages to plot. If :obj:`None`, plot all lineages.
    early_cells
        Cell ids or a mask marking early cells used to define the average fate probabilities. If :obj:`None`,
        use all cells. Only used when ``'kl_divergence'`` is in the ``keys``. If a :class:`dict`, key specifies
        a cluster key in :attr:`~anndata.AnnData.obs` and the values specify cluster labels containing early cells.
    lineage_order
        Can be one of the following:

        - :obj:`None` - it will be determined automatically, based on the number of lineages.
        - ``'optimal'`` - order lineages optimally by solving the Traveling salesman problem (TSP).
          Recommended for <= :math:`20` lineages.
        - ``'default'`` - use the order as specified by ``lineages``.
    metric
        Metric to use when constructing pairwise distance matrix when ``lineage_order = 'optimal'``. For available
        options, see :func:`~sklearn.metrics.pairwise_distances`.
    normalize_by_mean
        If :obj:`True`, normalize each lineage by its mean probability, as done in :cite:`velten:17`.
    ncols
        Number of columns when plotting multiple ``keys``.
    space
        Horizontal and vertical space between for :func:`~matplotlib.pyplot.subplots_adjust`.
    use_raw
        Whether to access :attr:`~anndata.AnnData.raw` when there are ``keys`` in :attr:`~anndata.AnnData.var_names`.
    text_kwargs
        Keyword arguments for :meth:`~matplotlib.axes.Axes.text`.
    label_distance
        Distance at which the lineage labels will be drawn.
    label_rot
        How to rotate the labels. Valid options are:

        - ``'best'`` - rotate labels so that they are easily readable.
        - ``'default'`` - use :mod:`matplotlib`'s default.
        - :obj:`None` - same as ``'default'``.
        - :class:`float`, all labels will be rotated by this many degrees.
    show_edges
        Whether to show the edges surrounding the simplex.
    key_added
        Key in :attr:`~anndata.AnnData.obsm` where to add the circular embedding. If :obj:`None`, it will be set to
        ``'X_fate_simplex_{fwd,bwd}'``, based on the ``backward``.
    %(plotting)s
    kwargs
        Keyword arguments for :func:`~scvelo.pl.scatter`.

    Returns
    -------
    %(just_plots)s Also updates ``adata`` with the following fields:

    - :attr:`adata.obsm['{key_added}'] <anndata.AnnData.obsm>` - the circular projection.
    - :attr:`adata.obs['to_{initial,terminal}_states_{method}'] <anndata.AnnData.obs>` - the priming degree,
      if a method is present in the ``keys``.
    """
    if label_distance is not None and label_distance < 0:
        raise ValueError(f"Expected `label_distance` to be positive, found `{label_distance}`.")

    if label_rot is None:
        label_rot = LabelRot.DEFAULT
    label_rot = LabelRot(label_rot)

    suffix = "bwd" if backward else "fwd"
    if key_added is None:
        key_added = "X_fate_simplex_" + suffix

    if isinstance(keys, str):
        keys = (keys,)

    keys = _unique_order_preserving(keys)
    keys_ = _check_collection(adata, keys, "obs", key_name="Observation", raise_exc=False) + _check_collection(
        adata, keys, "var_names", key_name="Gene", raise_exc=False, use_raw=use_raw
    )
    haystack = set(PrimingDegree)
    keys = keys_ + [k for k in keys if k in haystack]
    keys = _unique_order_preserving(keys)

    if not len(keys):
        raise ValueError("No valid keys have been selected.")

    probs = Lineage.from_adata(adata, backward=backward)
    if isinstance(lineages, str):
        lineages = [lineages]
    elif lineages is None:
        lineages = probs.names
    lineages = _unique_order_preserving(lineages)
    probs = probs[lineages]
    if probs.nlin < 3:
        raise ValueError(f"Expected at least `3` lineages, found `{probs.nlin}`.")

    X = probs.X.copy()
    if normalize_by_mean:
        X /= np.mean(X, axis=0)[None, :]
        X /= X.sum(1)[:, None]
        # this happens when cells for sel. lineages sum to 1 (or when the lineage average is 0, which is unlikely)
        X = np.nan_to_num(X, nan=1.0 / probs.nlin, copy=False)

    if lineage_order is None:
        lineage_order = LineageOrder.OPTIMAL if 3 < probs.nlin <= 20 else LineageOrder.DEFAULT
        logg.debug(f"Set ordering to `{lineage_order}`")
    lineage_order = LineageOrder(lineage_order)

    if lineage_order == LineageOrder.OPTIMAL:
        logg.info(f"Solving TSP for `{probs.nlin}` states")
        _, order = _get_optimal_order(X, metric=metric)
    else:
        order = np.arange(probs.nlin)

    probs = probs[:, order]
    X = X[:, order]

    angle_vec = np.linspace(0, 2 * np.pi, probs.nlin, endpoint=False)
    angle_vec_sin = np.cos(angle_vec)
    angle_vec_cos = np.sin(angle_vec)

    x = np.sum(X * angle_vec_sin, axis=1)
    y = np.sum(X * angle_vec_cos, axis=1)
    adata.obsm[key_added] = np.c_[x, y]

    nrows = int(np.ceil(len(keys) / ncols))
    fig, ax = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(ncols * 5, nrows * 5) if figsize is None else figsize,
        dpi=dpi,
    )

    fig.subplots_adjust(wspace=space, hspace=space)
    axes = np.ravel([ax])

    text_kwargs = dict(text_kwargs)
    text_kwargs["ha"] = "center"
    text_kwargs["va"] = "center"

    _i, lineage_key = 0, Key.obsm.fate_probs(backward)
    for _i, (k, ax) in enumerate(zip(keys, axes)):
        set_lognorm, colorbar = False, kwargs.pop("colorbar", True)
        try:
            _ = PrimingDegree(k)
            logg.debug(f"Calculating priming degree using `method={k}`")
            val = probs.priming_degree(method=k, early_cells=early_cells)
            k = f"{lineage_key}_{k}"
            adata.obs[k] = val
        except ValueError:
            # TODO(michalk8): parse the exception
            pass

        scv.pl.scatter(
            adata,
            basis=key_added,
            color=k,
            show=False,
            ax=ax,
            use_raw=use_raw,
            norm=LogNorm() if set_lognorm else None,
            colorbar=colorbar,
            **kwargs,
        )
        if colorbar and set_lognorm:
            cbar = ax.collections[0].colorbar
            cax = cbar.locator.axis
            ticks = cax.minor.locator.tick_values(cbar.vmin, cbar.vmax)
            ticks = [ticks[0], ticks[len(ticks) // 2 + 1], ticks[-1]]
            cbar.set_ticks(ticks)
            cbar.set_ticklabels([f"{t:.2f}" for t in ticks])
            cbar.update_ticks()

        patches, texts = ax.pie(
            np.ones_like(angle_vec),
            labeldistance=label_distance,
            rotatelabels=True,
            labels=probs.names[::-1],
            startangle=-360 / len(angle_vec) / 2,
            counterclock=False,
            textprops=text_kwargs,
        )

        for patch in patches:
            patch.set_visible(False)

        # clockwise
        for color, text in zip(probs.colors[::-1], texts):
            if isinstance(label_rot, (int, float, np.number)):
                text.set_rotation(label_rot)
            elif label_rot == LabelRot.BEST:
                rot = text.get_rotation()
                text.set_rotation(rot + 90 + (1 - rot // 180) * 180)
            elif label_rot != LabelRot.DEFAULT:
                raise NotImplementedError(f"Label rotation `{label_rot}` is not yet implemented.")
            text.set_color(color)

        if not show_edges:
            continue

        for i, color in enumerate(probs.colors):
            next = (i + 1) % probs.nlin
            x = 1.04 * np.linspace(angle_vec_sin[i], angle_vec_sin[next], _N)
            y = 1.04 * np.linspace(angle_vec_cos[i], angle_vec_cos[next], _N)
            points = np.array([x, y]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)

            cmap = LinearSegmentedColormap.from_list("fate_prob_cmap", [color, probs.colors[next]], N=_N)
            lc = LineCollection(segments, cmap=cmap, zorder=-1)
            lc.set_array(np.linspace(0, 1, _N))
            lc.set_linewidth(2)
            ax.add_collection(lc)

    for j in range(_i + 1, len(axes)):
        axes[j].remove()

    if save is not None:
        save_fig(fig, save)
