# -*- coding: utf-8 -*-
"""IO module."""

from typing import Union, Callable
from pathlib import Path

from matplotlib.colors import is_color_like

import scvelo as scv
import anndata
from scanpy import logging as logg

from cellrank.tools._colors import _create_categorical_colors
from cellrank.tools._lineage import Lineage
from cellrank.tools._constants import LinKey, Direction, _colors, _lin_names


def read(
    path: Union[Path, str], read_callback: Callable = scv.read, **kwargs
) -> anndata.AnnData:
    """
    Read file and return :class:`anndata.AnnData` object.

    Params
    ------
    path
        Path to the annotated data object.
    read_callback
        Function that actually reads the :class:`anndata.AnnData` object, such as
        :func:`scvelo.read` (default) or :func:`scanpy.read`.
    kwargs
        Keyword arguments for :paramref:`read_callback`.

    Returns
    -------
    :class:`anndata.AnnData`
        The annotated data object.
    """

    def maybe_create_lineage(direction: Direction):
        lin_key = str(
            LinKey.FORWARD if direction == Direction.FORWARD else LinKey.BACKWARD
        )
        names_key, colors_key = _lin_names(lin_key), _colors(lin_key)
        if lin_key in adata.obsm.keys():
            n_cells, n_lineages = adata.obsm[lin_key].shape
            logg.info(
                f"Creating {'forward' if direction == Direction.FORWARD else 'backward'} `Lineage` object"
            )

            if names_key not in adata.uns.keys():
                logg.warning(
                    f"Lineage names not found in `adata.uns[{names_key!r}]`, creating dummy names"
                )
                names = [f"Lineage {i}" for i in range(n_lineages)]
            elif len(adata.uns[names_key]) != n_lineages:
                logg.warning(
                    f"Lineage names are don't have the required length ({n_lineages}), creating dummy names"
                )
                names = [f"Lineage {i}" for i in range(n_lineages)]
            else:
                logg.info("Succesfully loaded names")
                names = adata.uns[names_key]

            if colors_key not in adata.uns.keys():
                logg.warning(
                    f"Lineage colors not found in `adata.uns[{colors_key!r}]`, creating new colors"
                )
                colors = _create_categorical_colors(n_lineages)
            elif len(adata.uns[colors_key]) != n_lineages or not all(
                map(lambda c: is_color_like(c), adata.uns[colors_key])
            ):
                logg.warning(
                    f"Lineage colors don't have the required length ({n_lineages}) "
                    f"or are not color-like, creating new colors"
                )
                colors = _create_categorical_colors(n_lineages)
            else:
                logg.info("Succesfully loaded colors")
                colors = adata.uns[colors_key]

            adata.obsm[lin_key] = Lineage(
                adata.obsm[lin_key], names=names, colors=colors
            )
            adata.uns[colors_key] = colors
            adata.uns[names_key] = names
        else:
            logg.debug(
                f"DEBUG: Unable to load {'forward' if direction == Direction.FORWARD else 'backward'} "
                f"`Lineage` from `adata.obsm[{lin_key!r}]`"
            )

    adata = read_callback(path, **kwargs)

    maybe_create_lineage(Direction.FORWARD)
    maybe_create_lineage(Direction.BACKWARD)

    return adata
