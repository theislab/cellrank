from typing import Any, Tuple, Union, Iterable, Optional, Sequence

from copy import copy
from pathlib import Path

from anndata import AnnData
from cellrank import logging as logg
from cellrank._key import Key
from cellrank.ul._docs import d
from cellrank.pl._utils import _position_legend, _get_categorical_colors
from cellrank.tl._utils import save_fig, _unique_order_preserving

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from seaborn import stripplot
from matplotlib.cm import ScalarMappable
from matplotlib.axes import Axes
from matplotlib.colors import Normalize, to_hex


@d.dedent
def log_odds(
    adata: AnnData,
    lineage_1: str,
    lineage_2: Optional[str] = None,
    time_key: str = "exp_time",
    backward: bool = False,
    keys: Optional[Union[str, Sequence[str]]] = None,
    threshold: Optional[Union[float, Sequence]] = None,
    threshold_color: str = "red",
    layer: Optional[str] = None,
    use_raw: bool = False,
    size: float = 2.0,
    cmap: str = "viridis",
    alpha: Optional[float] = 0.8,
    ncols: Optional[int] = None,
    fontsize: Optional[Union[float, str]] = None,
    xticks_step_size: Optional[int] = 1,
    legend_loc: Optional[str] = "best",
    jitter: Union[bool, float] = True,
    seed: Optional[int] = None,
    figsize: Optional[Tuple[float, float]] = None,
    dpi: Optional[int] = None,
    save: Optional[Union[str, Path]] = None,
    show: bool = True,
    **kwargs: Any,
) -> Optional[Union[Axes, Sequence[Axes]]]:
    """
    Plot log-odds ratio between lineages.

    Log-odds are plotted as a function of the experimental time.

    Parameters
    ----------
    %(adata)s
    lineage_1
        The first lineage for which to compute the log-odds.
    lineage_2
        The second lineage for which to compute the log-odds. If `None`, use the rest of the lineages.
    time_key
        Key in :attr:`anndata.AnnData.obs` containing the experimental time.
    %(backward)s
    keys
        Key in :attr:`anndata.AnnData.obs` or :attr:`anndata.AnnData.var_names`.
    threshold
        Visualize whether total expression per cell is greater than ``threshold``.
        If a :class:`typing.Sequence`, it should be the same length as ``keys``.
    threshold_color
        Color to use when plotting thresholded expression values.
    layer
        Which layer to use to get expression values. If `None` or `'X'`, use :attr:`anndata.AnnData.X`.
    use_raw
        Whether to access :attr:`anndata.AnnData.raw`. If `True`, ``layer`` is ignored.
    size
        Size of the dots.
    cmap
        Colormap to use for continuous variables in ``keys``.
    alpha
        Alpha values for the dots.
    ncols
        Number of columns.
    fontsize
        Size of the font for the title, x- and y-label.
    xticks_step_size
        Show only every n-th ticks on x-axis. If `None`, don't show any ticks.
    legend_loc
        Position of the legend. If `None`, do not show the legend.
    jitter
        Amount of jitter to apply along x-axis.
    seed
        Seed for ``jitter`` to ensure reproducibility.
    %(plotting)s
    show
        If `False`, return :class:`matplotlib.pyplot.Axes` or a sequence of them.
    kwargs
        Keyword arguments for :func:`seaborn.stripplot`.

    Returns
    -------
    The axes object(s), if ``show = False``.
    %(just_plots)s
    """
    from cellrank.tl.kernels._utils import _ensure_numeric_ordered

    def decorate(
        ax: Axes, *, title: Optional[str] = None, show_ylabel: bool = True
    ) -> None:
        ax.set_xlabel(time_key, fontsize=fontsize)
        ax.set_title(title, fontdict={"fontsize": fontsize})
        ax.set_ylabel(ylabel if show_ylabel else "", fontsize=fontsize)

        if xticks_step_size is None:
            ax.set_xticks([])
        else:
            step = max(1, xticks_step_size)
            ax.set_xticks(np.arange(0, n_cats, step))
            ax.set_xticklabels(df[time_key].cat.categories[::step])

    def cont_palette(values: np.ndarray) -> Tuple[np.ndarray, ScalarMappable]:
        cm = copy(plt.get_cmap(cmap))
        cm.set_bad("grey")
        sm = ScalarMappable(
            cmap=cm, norm=Normalize(vmin=np.nanmin(values), vmax=np.nanmax(values))
        )
        return np.array([to_hex(v) for v in (sm.to_rgba(values))]), sm

    def get_data(
        key: str,
        thresh: Optional[float] = None,
    ) -> Tuple[
        Optional[str], Optional[np.ndarray], Optional[np.ndarray], ScalarMappable
    ]:
        try:
            _, palette = _get_categorical_colors(adata, key)
            df[key] = adata.obs[key].values[mask]
            df[key] = df[key].cat.remove_unused_categories()
            try:
                # seaborn doesn't like numeric categories
                df[key] = df[key].astype(float)
                palette = {float(k): v for k, v in palette.items()}
            except ValueError:
                pass
            # otherwise seaborn plots all
            palette = {k: palette[k] for k in df[key].unique()}
            hue, thresh_mask, sm = key, None, None
        except TypeError:
            palette, hue, thresh_mask, sm = (
                cont_palette(adata.obs[key].values[mask])[0],
                None,
                None,
                None,
            )
        except KeyError:
            try:
                # fmt: off
                if thresh is None:
                    values = adata.raw.obs_vector(key) if use_raw else adata.obs_vector(key, layer=layer)
                    palette, sm = cont_palette(values)
                    hue, thresh_mask = None, None
                else:
                    if use_raw:
                        values = np.asarray(adata.raw[:, key].X[mask].sum(1)).squeeze()
                    elif layer not in (None, "X"):
                        values = np.asarray(adata[:, key].layers[layer][mask].sum(1)).squeeze()
                    else:
                        values = np.asarray(adata[:, key].X[mask].sum(1)).squeeze()
                    thresh_mask = values > thresh
                    hue, palette, sm = None, None, None
                # fmt: on
            except KeyError as e:
                raise e from None

        return hue, palette, thresh_mask, sm

    np.random.seed(seed)
    _ = kwargs.pop("orient", None)
    if use_raw and adata.raw is None:
        logg.warning("No raw attribute set. Setting `use_raw=False`")
        use_raw = False

    # define log-odds
    lineage_key = Key.obsm.abs_probs(backward)
    if lineage_key not in adata.obsm:
        raise KeyError(f"Lineages key `{lineage_key!r}` not found in `adata.obsm`.")
    time = _ensure_numeric_ordered(adata, time_key)
    order = time.cat.categories[:: -1 if backward else 1]

    fate1 = adata.obsm[lineage_key][lineage_1].X.squeeze(-1)
    if lineage_2 is None:
        fate2 = 1 - fate1
        ylabel = rf"$\log{{\frac{{{lineage_1}}}{{rest}}}}$"
    else:
        fate2 = adata.obsm[lineage_key][lineage_2].X.squeeze(-1)
        ylabel = rf"$\log{{\frac{{{lineage_1}}}{{{lineage_2}}}}}$"

    # fmt: off
    df = pd.DataFrame(
        {
            "log_odds": np.log(np.divide(fate1, fate2, where=fate2 != 0, out=np.zeros_like(fate1)) + 1e-12),
            time_key: time,
        }
    )
    mask = (fate1 != 0) & (fate2 != 0)
    df = df[mask]
    n_cats = len(df[time_key].cat.categories)
    # fmt: on

    if keys is None:
        if figsize is None:
            figsize = np.array([n_cats, n_cats * 4 / 6]) / 2

        fig, ax = plt.subplots(figsize=figsize, dpi=dpi, tight_layout=True)
        ax = stripplot(
            x=time_key,
            y="log_odds",
            data=df,
            order=order,
            jitter=jitter,
            color="k",
            size=size,
            ax=ax,
            **kwargs,
        )
        decorate(ax)
        if save is not None:
            save_fig(fig, save)
        return None if show else ax

    if isinstance(keys, str):
        keys = (keys,)
    if not len(keys):
        raise ValueError("No keys have been selected.")
    keys = _unique_order_preserving(keys)

    if not isinstance(threshold, Iterable):
        threshold = (threshold,) * len(keys)
    if len(threshold) != len(keys):
        raise ValueError(
            f"Expected `threshold` to be of length `{len(keys)}`, found `{len(threshold)}`."
        )

    ncols = max(len(keys) if ncols is None else ncols, 1)
    nrows = int(np.ceil(len(keys) / ncols))
    if figsize is None:
        figsize = np.array([n_cats * ncols, n_cats * nrows * 4 / 6]) / 2

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=figsize,
        dpi=dpi,
        constrained_layout=True,
        sharey="all",
    )
    axes = np.ravel([axes])

    i = 0
    for i, (key, ax, thresh) in enumerate(zip(keys, axes, threshold)):
        hue, palette, thresh_mask, sm = get_data(key, thresh)
        show_ylabel = i % ncols == 0

        ax = stripplot(
            x=time_key,
            y="log_odds",
            data=df if thresh_mask is None else df[~thresh_mask],
            hue=hue,
            order=order,
            jitter=jitter,
            color="black",
            palette=palette,
            size=size,
            alpha=alpha if alpha is not None else None if thresh_mask is None else 0.8,
            ax=ax,
            **kwargs,
        )
        if thresh_mask is not None:
            stripplot(
                x=time_key,
                y="log_odds",
                data=df if thresh_mask is None else df[thresh_mask],
                hue=hue,
                order=order,
                jitter=jitter,
                color=threshold_color,
                palette=palette,
                size=size * 2,
                alpha=0.9,
                ax=ax,
                **kwargs,
            )
            key = rf"${key} > {thresh}$"
        if sm is not None:
            cax = ax.inset_axes([1.02, 0, 0.025, 1], transform=ax.transAxes)
            fig.colorbar(sm, ax=ax, cax=cax)
        else:
            if legend_loc in (None, "none"):
                legend = ax.get_legend()
                if legend is not None:
                    legend.remove()
            else:
                handles, labels = ax.get_legend_handles_labels()
                if len(handles):
                    _position_legend(
                        ax, legend_loc=legend_loc, handles=handles, labels=labels
                    )

        decorate(ax, title=key, show_ylabel=show_ylabel)

    for ax in axes[i + 1 :]:
        ax.remove()
    axes = axes[: i + 1]

    if save is not None:
        save_fig(fig, save)

    return None if show else axes[0] if len(axes) == 1 else axes
