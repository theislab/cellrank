"""IO module."""

from typing import Union, TypeVar, Callable, Optional
from pathlib import Path

from scvelo import read as scv_read

from matplotlib.colors import is_color_like

from cellrank import logging as logg
from cellrank.ul._docs import d
from cellrank.tl._colors import _create_categorical_colors
from cellrank.tl._lineage import Lineage
from cellrank.tl._constants import (
    Direction,
    AbsProbKey,
    DirectionPlot,
    _colors,
    _lin_names,
)

AnnData = TypeVar("AnnData")


@d.dedent
def read(
    path: Union[Path, str],
    key: Optional[str] = None,
    read_callback: Callable = scv_read,
    **kwargs,
) -> AnnData:
    """
    Read file and return :class:`anndata.AnnData` object.

    Parameters
    ----------
    path
        Path to the annotated data object.
    key
        Key in ``adata.obsm`` where the :class:`cellrank.tl.Lineage` is stored.
        If `None`, it is determined automatically.
    read_callback
        Function that actually reads the :class:`anndata.AnnData` object, such as
        :func:`scvelo.read` (default) or :func:`scanpy.read`.
    kwargs
        Keyword arguments for ``read_callback``.

    Returns
    -------
    %(adata)s
    """

    def maybe_create_lineage(
        direction: Union[str, Direction], pretty_name: Optional[str] = None
    ):
        if isinstance(direction, Direction):
            lin_key = str(
                AbsProbKey.FORWARD
                if direction == Direction.FORWARD
                else AbsProbKey.BACKWARD
            )
        else:
            lin_key = direction

        pretty_name = "" if pretty_name is None else (pretty_name + " ")
        names_key, colors_key = _lin_names(lin_key), _colors(lin_key)

        if lin_key in adata.obsm.keys():
            n_cells, n_lineages = adata.obsm[lin_key].shape
            logg.info(f"Creating {pretty_name}`Lineage` from `adata.obsm[{lin_key!r}]`")

            if names_key not in adata.uns.keys():
                logg.warning(
                    f"    Lineage names not found in `adata.uns[{names_key!r}]`, creating new names"
                )
                names = [f"Lineage {i}" for i in range(n_lineages)]
            elif len(adata.uns[names_key]) != n_lineages:
                logg.warning(
                    f"    Lineage names are don't have the required length ({n_lineages}), creating new names"
                )
                names = [f"Lineage {i}" for i in range(n_lineages)]
            else:
                logg.info("    Successfully loaded names")
                names = adata.uns[names_key]

            if colors_key not in adata.uns.keys():
                logg.warning(
                    f"    Lineage colors not found in `adata.uns[{colors_key!r}]`, creating new colors"
                )
                colors = _create_categorical_colors(n_lineages)
            elif len(adata.uns[colors_key]) != n_lineages or not all(
                map(lambda c: is_color_like(c), adata.uns[colors_key])
            ):
                logg.warning(
                    f"    Lineage colors don't have the required length ({n_lineages}) "
                    f"or are not color-like, creating new colors"
                )
                colors = _create_categorical_colors(n_lineages)
            else:
                logg.info("    Successfully loaded colors")
                colors = adata.uns[colors_key]

            adata.obsm[lin_key] = Lineage(
                adata.obsm[lin_key], names=names, colors=colors
            )
            adata.uns[colors_key] = colors
            adata.uns[names_key] = names
        else:
            logg.debug(
                f"Unable to load {pretty_name}`Lineage` from `adata.obsm[{lin_key!r}]`"
            )

    adata = read_callback(path, **kwargs)

    if key is None:
        maybe_create_lineage(Direction.FORWARD, pretty_name=DirectionPlot.FORWARD.s)
        maybe_create_lineage(Direction.BACKWARD, pretty_name=DirectionPlot.BACKWARD.s)
    else:
        maybe_create_lineage(key)

    return adata
