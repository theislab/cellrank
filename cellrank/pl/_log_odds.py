from typing import Any, Tuple, Union, Optional, Sequence
from pathlib import Path

from anndata import AnnData
from cellrank import logging as logg
from cellrank.ul._docs import d
from cellrank.pl._utils import _get_categorical_colors
from cellrank.tl._utils import save_fig, _unique_order_preserving
from cellrank.tl._constants import AbsProbKey

import numpy as np
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt
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
    threshold: Optional[float] = None,
    threshold_color: str = "red",
    layer: Optional[str] = None,
    use_raw: bool = False,
    size: float = 2,
    cmap: str = "viridis",
    alpha: Optional[float] = 0.8,
    ncols: Optional[int] = None,
    figsize: Optional[Tuple[float, float]] = None,
    dpi: Optional[int] = None,
    save: Optional[Union[str, Path]] = None,
    **kwargs: Any,
) -> None:
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
        Visualize whether total expression per cell is greater or equal than ``threshold``.
    threshold_color
        Color to use when plotting thresholded expression values. Only used when ``threshold!=None``.
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
    %(plotting)s
    kwargs
        Keyword arguments for :func:`seaborn.stripplot`.

    Returns
    -------
    %(just_plots)s
    """
    from cellrank.tl.kernels._utils import _ensure_numeric_ordered

    def decorate(
        ax: Axes, *, title: Optional[str] = None, show_ylabel: bool = True
    ) -> None:
        ax.set_xlabel(time_key)
        ax.set_title(title)
        ax.set_ylabel(ylabel if show_ylabel else "")

    def cont_palette(values: np.ndarray) -> np.ndarray:
        sm = ScalarMappable(
            cmap=cmap, norm=Normalize(vmin=np.nanmin(values), vmax=np.nanmax(values))
        )
        return np.array([to_hex(v) for v in (sm.to_rgba(values))])

    def get_data(
        key: str,
    ) -> Tuple[Optional[str], Optional[np.ndarray], Optional[np.ndarray]]:
        try:
            _, palette = _get_categorical_colors(adata, key)
            df[key] = adata.obs[key].values[mask]
            hue, thresh_mask = key, None
        except TypeError:
            palette, hue, thresh_mask = (
                cont_palette(adata.obs[key].values[mask]),
                None,
                None,
            )
        except KeyError:
            try:
                # fmt: off
                if threshold is None:
                    values = (adata.raw if use_raw else adata).obs_vector(key, layer=layer)
                    palette = cont_palette(values)
                    hue, thresh_mask = None, None
                else:
                    if use_raw:
                        values = np.asarray(adata.raw[:, key].X[mask].sum(1)).squeeze()
                    elif layer not in (None, "X"):
                        values = np.asarray(adata[:, key].layers[layer][mask].sum(1)).squeeze()
                    else:
                        values = np.asarray(adata[:, key].X[mask].sum(1)).squeeze()
                    thresh_mask = values >= threshold
                    hue, palette = None, None
                # fmt: on
            except KeyError as e:
                raise e from None

        return hue, palette, thresh_mask

    if use_raw and adata.raw is None:
        logg.warning("No raw attribute set. Setting `use_raw=False`")
        use_raw = False

    # define log-odds
    ln_key = str(AbsProbKey.BACKWARD if backward else AbsProbKey.FORWARD)
    if ln_key not in adata.obsm:
        raise KeyError(f"Lineages key `{ln_key!r}` not found in `adata.obsm`.")
    time = _ensure_numeric_ordered(adata, time_key, convert=False)
    order = time.cat.categories[:: -1 if backward else 1]
    jitter = np.nanmax(np.abs(time.cat.categories[:-1] - time.cat.categories[1:])) / 2.5

    fate1 = adata.obsm[ln_key][lineage_1].X.squeeze(-1)
    if lineage_2 is None:
        fate2 = 1 - fate1
        ylabel = rf"$\log{{\frac{{{lineage_1}}}{{rest}}}}$"
    else:
        fate2 = adata.obsm[ln_key][lineage_2].X.squeeze(-1)
        ylabel = rf"$\log{{\frac{{{lineage_1}}}{{{lineage_2}}}}}$"

    # fmt: off
    df = pd.DataFrame(
        {
            "log_odds": np.log(np.divide(fate1, fate2, where=fate2 != 0, out=np.zeros_like(fate1) + 1e-12)),
            time_key: time,
        }
    )
    mask = (fate1 != 0) & (fate2 != 0)
    df = df[mask]
    # fmt: on

    if keys is None:
        if figsize is None:
            n_cats = len(time.cat.categories)
            figsize = np.array([n_cats, n_cats * 4 / 6]) / 4

        fig, ax = plt.subplots(figsize=figsize, dpi=dpi, tight_layout=True)
        ax = sns.stripplot(
            time_key,
            "log_odds",
            data=df,
            order=order,
            jitter=0.4,
            color="k",
            size=size,
            ax=ax,
            **kwargs,
        )
        decorate(ax)
        return

    if isinstance(keys, str):
        keys = (keys,)
    if not len(keys):
        raise ValueError("No keys have been selected.")
    keys = _unique_order_preserving(keys)

    ncols = max(len(keys) if ncols is None else ncols, 1)
    nrows = int(np.ceil(len(keys) / ncols))
    if figsize is None:
        n_cats = len(time.cat.categories)
        figsize = np.array([n_cats * ncols, n_cats * nrows * 4 / 6]) / 4

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=figsize,
        dpi=dpi,
        tight_layout=True,
        sharey="all",
    )
    axes = np.ravel([axes])

    i = 0
    for i, (key, ax) in enumerate(zip(keys, axes)):
        hue, palette, thresh_mask = get_data(key)
        show_ylabel = i % ncols == 0

        ax = sns.stripplot(
            time_key,
            "log_odds",
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
            ax = sns.stripplot(
                time_key,
                "log_odds",
                data=df if thresh_mask is None else df[thresh_mask],
                hue=hue,
                order=order,
                jitter=jitter,
                color=threshold_color,
                palette=palette,
                size=size * 1.5,
                alpha=0.9,
                ax=ax,
                **kwargs,
            )
            key = rf"${key} \geq {threshold}$"
        decorate(ax, title=key, show_ylabel=show_ylabel)

    for ax in axes[i + 1 :]:
        ax.remove()

    if save is not None:
        save_fig(fig, save)
