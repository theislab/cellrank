# -*- coding: utf-8 -*-
"""Module for pl lineage-related stuff."""
from typing import Union, Iterable, Optional

import pandas as pd

import matplotlib as mpl
from matplotlib import cm as cm

import cellrank.logging as logg
from cellrank.ul._docs import d
from cellrank.pl._utils import AnnData
from cellrank.tl._constants import DirPrefix
from cellrank.tl.estimators import GPCCA
from cellrank.tl.estimators._constants import A, P
from cellrank.tl.kernels._precomputed_kernel import DummyKernel


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
    %(adata)s
    lineages
        Plot only these lineages. If `None`, plot all lineages.
    %(backward)s
    cluster_key
        If given, plot cluster annotations left of the lineage probabilities.
    %(time_mode)s
    time_key
        Key from ``adata.obs`` to use as a pseudotime ordering of the cells.
    cmap
        Colormap to use.
    **kwargs
        Keyword arguments for :func:`scvelo.pl.scatter`.

    Returns
    -------
    %(just_plots)s
    """

    pk = DummyKernel(adata, backward=backward)
    mc = GPCCA(pk, read_from_adata=True, write_to_adata=False)
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
    %(just_plots)s
    """

    pk = DummyKernel(adata, backward=backward)
    mc = GPCCA(pk, read_from_adata=True, write_to_adata=False)

    if use_raw and adata.raw is None:
        logg.warning("No raw attribute set. Using `adata.var` instead")
        use_raw = False

    needle = f"{DirPrefix.BACKWARD if backward else DirPrefix.FORWARD} {lineage}"

    haystack = adata.raw.var if use_raw else adata.var

    if needle not in haystack:
        raise RuntimeError(
            f"Unable to find lineage drivers in "
            f"`{'adata.raw.var' if use_raw else 'adata.var'}[{needle!r}]`. "
            f"Compute lineage drivers first as `cellrank.tl.lineage_drivers(lineages={lineage!r}, "
            f"use_raw={use_raw}, backward={backward}).`"
        )

    drivers = pd.DataFrame(haystack[needle])
    drivers.columns = [lineage]
    mc._set(A.LIN_DRIVERS, drivers)

    mc.plot_lineage_drivers(lineage, n_genes=n_genes, **kwargs)
