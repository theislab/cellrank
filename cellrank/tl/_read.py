from typing import Any, Union, Callable, Optional

from pathlib import Path

from scvelo import read as scv_read
from anndata import AnnData
from cellrank import logging as logg
from cellrank._key import Key
from cellrank.ul._docs import d
from cellrank.tl._utils import _deprecate
from cellrank.tl._colors import _create_categorical_colors
from cellrank.tl._lineage import Lineage

from matplotlib.colors import is_color_like


@d.dedent
@_deprecate(version="2.0")
def read(
    path: Union[Path, str],
    read_callback: Callable = scv_read,
    **kwargs: Any,
) -> AnnData:
    """
    Read file and return :class:`anndata.AnnData` object.

    Parameters
    ----------
    path
        Path to the annotated data object.
    read_callback
        Function that actually reads the :class:`anndata.AnnData` object, such as
        :func:`scvelo.read` (default) or :func:`scanpy.read`.
    kwargs
        Keyword arguments for ``read_callback``.

    Returns
    -------
    %(adata)s
    """

    def maybe_create_lineage(backward: bool, pretty_name: Optional[str] = None) -> None:
        lin_key = Key.obsm.abs_probs(backward)
        pretty_name = "" if pretty_name is None else (pretty_name + " ")
        names_key = Key.obs.term_states(backward)
        colors_key = Key.uns.colors(names_key)

        if lin_key in adata.obsm.keys():
            n_cells, n_lineages = adata.obsm[lin_key].shape
            logg.info(f"Creating {pretty_name}`Lineage` from `adata.obsm[{lin_key!r}]`")

            if names_key not in adata.obs:
                logg.warning(
                    f"    Lineage names not found in `adata.uns[{names_key!r}]`, creating new names"
                )
                names = [f"Lineage {i}" for i in range(n_lineages)]
            elif len(adata.obs[names_key].cat.categories) != n_lineages:
                logg.warning(
                    f"    Lineage names are don't have the required length ({n_lineages}), creating new names"
                )
                names = [f"Lineage {i}" for i in range(n_lineages)]
            else:
                logg.info("    Successfully loaded names")
                names = list(adata.obs[names_key].cat.categories)

            if colors_key not in adata.uns:
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

    maybe_create_lineage(False, pretty_name="forward")
    maybe_create_lineage(True, pretty_name="backward")

    return adata
