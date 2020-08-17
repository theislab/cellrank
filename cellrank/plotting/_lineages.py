# -*- coding: utf-8 -*-
"""Module for plotting lineage-related stuff."""
from typing import Union, Iterable, Optional

import pandas as pd
from scipy.sparse import eye as speye

import matplotlib as mpl
from matplotlib import cm as cm

import cellrank.logging as logg
from cellrank.tools import CFLARE
from cellrank.utils._docs import d
from cellrank.tools.kernels import PrecomputedKernel
from cellrank.plotting._utils import AnnData
from cellrank.tools._constants import DirPrefix
from cellrank.tools.estimators._constants import A, P


@d.dedent
def lineages(
    adata: AnnData,
    lineages: Optional[Union[str, Iterable[str]]] = None,
    backward: bool = False,
    cluster_key: Optional[str] = None,
    mode: str = "embedding",
    time_key: str = "latent_time",
    cmap: Union[str, mpl.colors.ListedColormap] = cm.viridis,
    **kwargs,
) -> None:
    """
    Plot lineages that were uncovered using :func:`cellrank.tl.lineages`.

    For each lineage, we show all cells in an embedding (default is UMAP, but can be any) and color them by their
    probability of belonging to this lineage. For cells that are already committed, this probability will be one for
    their respective lineage and zero otherwise. For naive cells, these probabilities will be more balanced, reflecting
    the fact that naive cells have the potential to develop towards multiple endpoints.

    .. image:: https://raw.githubusercontent.com/theislab/cellrank/master/resources/images/lineages.png
       :width: 400px
       :align: center

    Parameters
    ----------
    %s(adata)s
    lineages
        Plot only these lineages. If `None`, plot all lineages.
    %(backward)s
    cluster_key
        If given, plot cluster annotations left of the lineage probabilities.
    mode
        Can be either `'embedding'` or `'time'`:

            - `'embedding'` - plot the embedding while coloring in the absorption probabilities.
            - `'time'` - plot the pseudotime on x-axis and the absorption probabilities on y-axis.
    time_key
        Key from `adata.obs` to use as a pseudotime ordering of the cells.
    cmap
        Colormap to use.
    **kwargs
        Keyword arguments for :func:`scvelo.pl.scatter`.

    Returns
    -------
    %s(just_plots)s
    """

    # create a dummy kernel object
    pk = PrecomputedKernel(
        speye(adata.n_obs, adata.n_obs, format="csr"), adata=adata, backward=backward
    )

    # use this to initialize an MC object
    mc = CFLARE(pk, read_from_adata=True, write_to_adata=False)
    if mc._get(P.ABS_PROBS) is None:
        raise RuntimeError(
            "Compute absorption probabilities first as `cellrank.tl.lineages()`."
        )

    # plot using the MC object
    mc.plot_absorption_probabilities(
        lineages=lineages,
        cluster_key=cluster_key,
        mode=mode,
        time_key=time_key,
        cmap=cmap,
        **kwargs,
    )


@d.dedent
def lineage_drivers(
    adata: AnnData,
    lineage: str,
    backward: bool = False,
    n_genes: int = 10,
    use_raw: bool = False,
    **kwargs,
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
    %s(just_plots)s
    """

    # create a dummy kernel object
    pk = PrecomputedKernel(
        speye(adata.n_obs, adata.n_obs, format="csr"), adata=adata, backward=backward
    )

    # use this to initialize an MC object
    mc = CFLARE(pk, read_from_adata=True, write_to_adata=False)

    needle = f"{DirPrefix.BACKWARD if backward else DirPrefix.FORWARD} {lineage}"

    if use_raw and adata.raw is None:
        logg.warning("No raw attribute set. Using `adata.var` instead")
        use_raw = False

    haystack = adata.raw.var if use_raw else adata.var

    if needle not in haystack:
        raise RuntimeError(
            f"Unable to find lineage drivers in "
            f"`{'adata.raw.var' if use_raw else 'adata.var'}[{needle!r}]`. "
            f"Try computing lineage drivers as `cellrank.tl.lineage_drivers(lineages={lineage!r}, "
            f"use_raw={use_raw}).`"
        )

    drivers = pd.DataFrame(haystack[needle])
    drivers.columns = [lineage]
    mc._set(A.LIN_DRIVERS, drivers)

    mc.plot_lineage_drivers(lineage, n_genes=n_genes, **kwargs)
