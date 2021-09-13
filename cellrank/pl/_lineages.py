from typing import Any, Union, Optional, Sequence
from typing_extensions import Literal

import cellrank.logging as logg
from cellrank._key import Key
from cellrank.ul._docs import d
from cellrank.pl._utils import AnnData
from cellrank.tl.estimators.terminal_states._cflare import CFLARE


@d.dedent
def lineages(
    adata: AnnData,
    lineages: Optional[Union[str, Sequence[str]]] = None,
    backward: bool = False,
    color: Optional[str] = None,
    mode: Literal["embedding", "time"] = "embedding",
    time_key: str = "latent_time",
    **kwargs: Any,
) -> None:
    """
    Plot lineages that were uncovered using :func:`cellrank.tl.lineages`.

    For each lineage, we show all cells in an embedding (default is UMAP) and color them by their probability of
    belonging to this lineage. For cells that are already committed, this probability will be one for  their respective
    lineage and zero otherwise. For naive cells, these probabilities will be more balanced, reflecting
    the fact that naive cells have the potential to develop towards multiple endpoints.

    Parameters
    ----------
    %(adata)s
    lineages
        Plot only these lineages. If `None`, plot all lineages.
    %(backward)s
    color
        If given, plot cluster annotations left of the lineage probabilities.
    %(time_mode)s
    time_key
        Key in ``adata.obs`` where the pseudotime is stored.
    %(basis)s
    kwargs
        Keyword arguments for :meth:`cellrank.tl.estimators.GPCCA.plot_absorption_probabilities`.

    Returns
    -------
    %(just_plots)s
    """

    # use CFLARE (no macrostates required)
    mc = CFLARE.from_adata(adata, obsp_key=Key.uns.kernel(backward))
    if mc.absorption_probabilities is None:
        raise RuntimeError(
            f"Compute absorption probabilities first as `cellrank.tl.lineages(..., backward={backward})`."
        )

    # plot using the MC object
    color = kwargs.pop("cluster_key", color)
    mc.plot_absorption_probabilities(
        states=lineages,
        color=color,
        mode=mode,
        time_key=time_key,
        **kwargs,
    )


@d.dedent
def lineage_drivers(
    adata: AnnData,
    lineage: str,
    backward: bool = False,
    n_genes: int = 8,
    ncols: Optional[int] = None,
    use_raw: bool = False,
    title_fmt: str = "{gene} qval={qval:.4e}",
    **kwargs: Any,
) -> None:
    """
    Plot lineage drivers that were uncovered using :func:`cellrank.tl.lineage_drivers`.

    Parameters
    ----------
    %(adata)s
    %(backward)s
    %(plot_lineage_drivers.parameters)s

    Returns
    -------
    %(just_plots)s
    """

    mc = CFLARE.from_adata(adata, obsp_key=Key.uns.kernel(backward))

    if use_raw and adata.raw is None:
        logg.warning("No raw attribute set. Using `use_raw=False`")
        use_raw = False

    key = Key.varm.lineage_drivers(backward)
    haystack = (adata.raw if use_raw else adata).varm

    if key not in haystack:
        raise RuntimeError(
            "Unable to find lineage drivers in "
            f"`{'adata.raw.varm' if use_raw else 'adata.varm'}[{key!r}]`. "
            f"Compute lineage drivers first as `cellrank.tl.lineage_drivers(lineages={lineage!r}, "
            f"use_raw={use_raw}, backward={backward}).`"
        )

    mc._lineage_drivers = haystack[key]
    mc.plot_lineage_drivers(
        lineage,
        n_genes=n_genes,
        use_raw=use_raw,
        ncols=ncols,
        title_fmt=title_fmt,
        **kwargs,
    )
