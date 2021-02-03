from types import MappingProxyType
from typing import Any, Tuple, Union, Mapping, Callable, Optional, Sequence
from pathlib import Path

import scvelo as scv
from anndata import AnnData

import numpy as np
import pandas as pd
from sklearn.metrics import pairwise_distances

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, LinearSegmentedColormap
from matplotlib.collections import LineCollection

from cellrank import logging as logg
from cellrank.tl import Lineage
from cellrank.ul._docs import d
from cellrank.pl._utils import _held_karp
from cellrank.tl._utils import save_fig, _unique_order_preserving
from cellrank.ul._utils import _check_collection
from cellrank.tl._constants import ModeEnum, AbsProbKey


class LineageOrder(ModeEnum):  # noqa: D101
    DEFAULT = "default"
    OPTIMAL = "optimal"


class SpecialKey(ModeEnum):  # noqa: D101
    DEGREE = "priming_degree"
    DIRECTION = "priming_direction"


class LabelRot(ModeEnum):  # noqa: D101
    DEFAULT = "default"
    BEST = "best"


Metric_T = Union[str, Callable, np.ndarray, pd.DataFrame]
_N = 200


def _get_distances(data: Union[np.ndarray, Lineage], metric: Metric_T) -> np.ndarray:
    if isinstance(data, Lineage):
        data = data.X

    if isinstance(metric, str) or callable(metric):
        metric = pairwise_distances(data.T, metric=metric)
    elif isinstance(metric, np.ndarray):
        shape = (data.shape[1], data.shape[1])
        if metric.shape != shape:
            raise ValueError(
                f"Expected a `numpy.ndarray` of shape `{shape}`, found `{metric.shape}`."
            )
    elif isinstance(metric, pd.DataFrame):
        if np.any(metric.index != data.names):
            raise ValueError(
                f"Expected `pandas.DataFrame` to have the following index `{list(data.names)}`, "
                f"found `{list(metric.columns)}`."
            )
        if np.any(metric.columns != data.names):
            raise ValueError(
                f"Expected `pandas.DataFrame` to have the following columns `{list(data.names)}`, "
                f"found `{list(metric.columns)}`."
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


def _priming_degree(probs: Union[np.ndarray, Lineage]) -> np.ndarray:
    if isinstance(probs, Lineage):
        probs = probs.X
    return np.sum(probs * np.log2(probs / np.mean(probs, axis=0)), axis=1)


def _priming_direction(probs: Union[np.ndarray, Lineage]) -> np.ndarray:
    if isinstance(probs, Lineage):
        probs = probs.X
    return np.argmax(probs / np.sum(probs, axis=0), axis=1)


@d.dedent
def circular_projection(
    adata: AnnData,
    keys: Union[str, Sequence[str]],
    backward: bool = False,
    lineages: Optional[Union[str, Sequence[str]]] = None,
    lineage_order: Optional[str] = None,
    metric: Union[str, Callable, np.ndarray, pd.DataFrame] = "correlation",
    normalize_by_mean: bool = True,
    ncols: int = 4,
    space: float = 0.25,
    use_raw: bool = False,
    text_kwargs: Mapping[str, Any] = MappingProxyType({}),
    labeldistance: float = 1.25,
    labelrot: Union[str, float] = "best",
    show_edges: bool = True,
    key_added: Optional[str] = None,
    figsize: Optional[Tuple[float, float]] = None,
    dpi: Optional[int] = None,
    save: Optional[Union[str, Path]] = None,
    **kwargs,
):
    r"""
    Plot absorption probabilities on a circular embedding as done in [Velten17]_.

    Parameters
    ----------
    %(adata)s
    keys
        Keys in :attr:`anndata.AnnData.obs` or :attr:`anndata.AnnData.var_names`. Extra available keys:

            - `'priming_degree'`: the entropy relative to the "average" cell, as described in [Velten17]_.
              It follows  the formula :math:`\sum_{j} p_{i, j} \log \frac{p_{ij}}{\overline{p_j}}` where :math:`i`
              corresponds to a cell, :math:`j` corresponds to a lineage and :math:`\overline{p_j}` is the average
              probability for lineage :math:`j`.
            - `'priming_direction'`: defined as :math:`d_i = argmax_j \frac{p_{ij}}{\overline{p_j}}` in [Velten17]_
              where :math:`i`, :math:`j` and :math:`\overline{p_j}` are defined as above.

    %(backward)s
    lineages
        Lineages to plot. If `None`, plot all lineages.
    lineage_order
        Can be one of the following:

            - `None`: it will determined automatically, based on the number of lineages.
            - `'optimal'`: order the lineages optimally by solving the Travelling salesman problem (TSP).
              Recommended for <= `20` lineages.
            - `'default'`: use the order as specified in ``lineages``.

    metric
        Metric to use when constructing pairwise distance matrix when ``lineage_order = 'optimal'``. For available
        options, see :func:`sklearn.metrics.pairwise_distances`.
    normalize_by_mean
        If `True`, normalize each lineage by its mean probability, as done in [Velten17]_.
    ncols
        Number of columns when plotting multiple ``keys``.
    space
        Horizontal and vertical space between for :func:`matplotlib.pyplot.subplots_adjust`.
    use_raw
        Whether to access :attr:`anndata.AnnData.raw` when there are ``keys`` in :attr:`anndata.AnnData.var_names`.
    text_kwargs
        Keyword arguments for :func:`matplotlib.pyplot.text`.
    labeldistance
        Distance at which the lineage labels will be drawn.
    labelrot
        How to rotate the labels. Valid options are:

            - `'best'`: rotate labels so that they are easily readable.
            - `'default'`: use :mod:`matplotlib`'s default.
            - `None`: same as `'default'`.

        If a :class:`float`, all labels will be rotated by this many degrees.
    show_edges
        Whether to show the edges surrounding the simplex.
    key_added
        Key in :attr:`anndata.AnnData.obsm` where to add the circular embedding. If `None`, it will be set to
        `'X_fate_simplex_{fwd,bwd}'`, based on ``backward``.
    %(plotting)s
    kwargs
        Keyword arguments for :func:`scvelo.pl.scatter`.

    Returns
    -------
    %(just_plots)s
        Also updates ``adata`` with the following fields:

            - :attr:`anndata.AnnData.obsm` ``[{key_added}]``: the circular projection.
            - :attr:`anndata.AnnData.obs` ``['priming_degree_{fwd,bwd}']``: the priming degree, if in ``keys``.
            - :attr:`anndata.AnnData.obs` ``['priming_direction_{fwd,bwd}']``: the priming direction, if in ``keys``.
    """
    if labeldistance is not None and labeldistance < 0:
        raise ValueError(f"Expected `delta` to be positive, found `{labeldistance}`.")

    if labelrot is None:
        labelrot = LabelRot.DEFAULT
    if isinstance(labelrot, str):
        labelrot = LabelRot(labelrot)

    suffix = "bwd" if backward else "fwd"
    if key_added is None:
        key_added = "X_fate_simplex_" + suffix

    if isinstance(keys, str):
        keys = (keys,)

    keys = _unique_order_preserving(keys)
    keys_ = _check_collection(
        adata, keys, "obs", key_name="Observation", raise_exc=False
    ) + _check_collection(
        adata, keys, "var_names", key_name="Gene", raise_exc=False, use_raw=use_raw
    )
    haystack = {s.s for s in SpecialKey}
    keys = keys_ + [k for k in keys if k in haystack]
    keys = _unique_order_preserving(keys)

    if not len(keys):
        raise ValueError("No valid keys have been selected.")

    lineage_key = str(AbsProbKey.BACKWARD if backward else AbsProbKey.FORWARD)
    if lineage_key not in adata.obsm:
        raise KeyError(f"Lineages key `{lineage_key!r}` not found in `adata.obsm`.")

    probs = adata.obsm[lineage_key]

    if isinstance(lineages, str):
        lineages = (lineages,)
    elif lineages is None:
        lineages = probs.names

    probs = adata.obsm[lineage_key][lineages]
    n_lin = probs.shape[1]
    if n_lin <= 2:
        raise ValueError(f"Expected at least `3` lineages, found `{n_lin}`")

    X = probs.X.copy()
    if normalize_by_mean:
        X /= np.mean(X, axis=0)[None, :]
        X /= X.sum(1)[:, None]
        # this happens when cells for sel. lineages sum to 1 (or when the lineage average is 0, which is unlikely)
        X = np.nan_to_num(X, nan=1.0 / n_lin, copy=False)

    if lineage_order is None:
        lineage_order = LineageOrder.OPTIMAL if n_lin <= 15 else LineageOrder.DEFAULT
        logg.debug(f"Set ordering to `{lineage_order}`")
    lineage_order = LineageOrder(lineage_order)

    if lineage_order == LineageOrder.OPTIMAL:
        logg.info(f"Solving TSP for `{n_lin}` states")
        _, order = _get_optimal_order(X, metric=metric)
    else:
        order = np.arange(n_lin)

    probs = probs[:, order]
    X = X[:, order]

    angle_vec = np.linspace(0, 2 * np.pi, n_lin, endpoint=False)
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

    _i = 0
    for _i, (k, ax) in enumerate(zip(keys, axes)):

        set_lognorm, colorbar = False, kwargs.pop("colorbar", True)
        try:
            k = SpecialKey(k)
            logg.debug(f"Calculating `{k}`")
            if k == SpecialKey.DEGREE:
                val = _priming_degree(probs)
                set_lognorm = True
            elif k == SpecialKey.DIRECTION:
                val = pd.Series(
                    probs.names[_priming_direction(probs)], dtype="category"
                )
                color_mapper = dict(zip(probs.names, probs.colors))
                adata.uns[f"{k.s}_{suffix}_colors"] = [
                    color_mapper[c] for c in val.cat.categories
                ]
                val = val.values
            else:
                raise NotImplementedError(f"Special key `{k}` is not implemented.")
            k = f"{k}_{suffix}"
            if k in adata.obs:
                logg.debug(f"Overwriting `adata.obs[{k!r}]`")
            adata.obs[k] = val
        except ValueError:
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
            labeldistance=labeldistance,
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
            if isinstance(labelrot, (int, float)):
                text.set_rotation(labelrot)
            elif labelrot == LabelRot.BEST:
                rot = text.get_rotation()
                text.set_rotation(rot + 90 + (1 - rot // 180) * 180)
            elif labelrot != LabelRot.DEFAULT:
                raise NotImplementedError(
                    f"Label rotation `{labelrot}` is not yet implemented."
                )
            text.set_color(color)

        if not show_edges:
            continue

        for i, color in enumerate(probs.colors):
            next = (i + 1) % n_lin
            x = 1.04 * np.linspace(angle_vec_sin[i], angle_vec_sin[next], _N)
            y = 1.04 * np.linspace(angle_vec_cos[i], angle_vec_cos[next], _N)
            points = np.array([x, y]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)

            cmap = LinearSegmentedColormap.from_list(
                "abs_prob_cmap", [color, probs.colors[next]], N=_N
            )
            lc = LineCollection(segments, cmap=cmap, zorder=-1)
            lc.set_array(np.linspace(0, 1, _N))
            lc.set_linewidth(2)
            ax.add_collection(lc)

    for j in range(_i + 1, len(axes)):
        axes[j].remove()

    if save is not None:
        save_fig(fig, save)
