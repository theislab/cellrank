"""Private plotting utilities formerly imported from scVelo."""

from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from anndata import AnnData
from matplotlib.axes import Axes

__all__: list[str] = []


def _quiver_autoscale(X_emb: np.ndarray, V_emb: np.ndarray) -> float:
    """Compute auto-scale factor for quiver arrow plots.

    Creates a temporary quiver plot to determine the appropriate scale
    factor, then removes it.

    Parameters
    ----------
    X_emb
        Embedding coordinates, shape ``(n_cells, 2)``.
    V_emb
        Embedding velocity vectors, shape ``(n_cells, 2)``.

    Returns
    -------
    The quiver scale normalized by the embedding scale factor.
    """
    scale_factor = np.abs(X_emb).max()
    fig, ax = plt.subplots()
    Q = ax.quiver(
        X_emb[:, 0] / scale_factor,
        X_emb[:, 1] / scale_factor,
        V_emb[:, 0],
        V_emb[:, 1],
        angles="xy",
        scale_units="xy",
        scale=None,
    )
    Q._init()
    fig.clf()
    plt.close(fig)
    return Q.scale / scale_factor


def _default_size(adata: AnnData) -> float:
    """Compute default scatter point size based on the number of cells.

    Parameters
    ----------
    adata
        Annotated data matrix.

    Returns
    -------
    Point size scaled inversely with cell count.
    """
    return 1.2e5 / adata.n_obs


def _plot_outline(
    x: np.ndarray,
    y: np.ndarray,
    *,
    outline_color: tuple[str, str] = ("black", "white"),
    outline_width: tuple[float, float] = (0.3, 0.05),
    kwargs: dict[str, Any] | None = None,
    ax: Axes | None = None,
    **scatter_kwargs: Any,
) -> None:
    """Draw outlined scatter points using a double-scatter technique.

    Draws two layers: a large background scatter (``outline_color[0]``) and
    a medium gap scatter (``outline_color[1]``).  The caller typically draws
    the actual data scatter on top.

    The outline sizes are computed from the base ``s`` value following scVelo's
    quadratic formula so that ``outline_width`` controls the visual thickness.

    Parameters
    ----------
    x
        X coordinates.
    y
        Y coordinates.
    outline_color
        Tuple of ``(background_color, gap_color)``.
    outline_width
        ``(bg_frac, gap_frac)`` relative to the dot radius.
    kwargs
        Base keyword arguments for all scatter calls (e.g. ``s``, ``alpha``).
    ax
        Matplotlib axes to draw on. If :obj:`None`, uses current axes.
    **scatter_kwargs
        Additional keyword arguments forwarded to all scatter calls (e.g. ``zorder``).
    """
    if ax is None:
        ax = plt.gca()
    if kwargs is None:
        kwargs = {}

    bg_color, gap_color = outline_color
    bg_frac, gap_frac = outline_width

    s = kwargs.pop("s", 20)
    point = np.sqrt(s)
    gap_size = (2 * point * gap_frac + point) ** 2
    bg_size = (2 * point * bg_frac + np.sqrt(gap_size)) ** 2

    base = {**kwargs, **scatter_kwargs, "edgecolors": "none", "marker": "."}

    ax.scatter(x, y, s=bg_size, color=bg_color, **base)
    ax.scatter(x, y, s=gap_size, color=gap_color, **base)
