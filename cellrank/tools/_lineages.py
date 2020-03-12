# -*- coding: utf-8 -*-
from typing import Optional, Sequence

from anndata import AnnData
from scanpy import logging as logg

from cellrank.tools import MarkovChain
from cellrank.tools._constants import LinKey, RcKey, _transition, Direction
from cellrank.tools.kernels import VelocityKernel


def lineages(
    adata: AnnData,
    final: bool = True,
    keys: Optional[Sequence[str]] = None,
    copy: bool = False,
) -> Optional[AnnData]:
    """
    Computes probabilistic lineage assignment using RNA velocity.

    For each cell i in {1, ..., N} and start/endpoint j in {1, ..., M}, the probability is computed that cell i
    is either going to j (end point) or coming from j (start point). Mathematically, this computes absorption
    probabilities to approximate recurrent classes using an RNA velocity based Markov Chain.

    Note that absorption probabilities have been used in the single cell context to infer lineage probabilities e.g.
    in [Setty19]_ or [Weinreb18]_ and we took inspiration from there.

    Before running this function, compute start/endpoints using :func:`cellrank.tl.root_final`.

    Parameters
    --------
    adata : :class:`anndata.AnnData`
        Annotated data object
    final
        If `True`, computes final cells, i.e. end points.
        Otherwise, computes root cells, i.e. starting points.
    keys
        Determines which end/start-points to use by passing their names. Further, start/end-points can be combined.
        If e.g. the endpoints are ['Neuronal_1', 'Neuronal_1', 'Astrocytes', 'OPC'], then passing
        keys=['Neuronal_1, Neuronal_2', 'OPC'] means that the two neuronal endpoints are treated as one and
        Astrocytes are excluded.
    copy
        Whether to update the existing AnnData object or to return a copy.

    Returns
    --------
    :class:`anndata.AnnData` or :class:`NoneType`
        Depending on :paramref:`copy`, either updates the existing :paramref:`adata` object or returns a copy.
    """

    # Set the keys and print info
    adata = adata.copy() if copy else adata

    if final:
        direction = Direction.FORWARD
        lin_key = LinKey.FORWARD
        rc_key = RcKey.FORWARD
    else:
        direction = Direction.BACKWARD
        lin_key = LinKey.BACKWARD
        rc_key = RcKey.BACKWARD

    transition_key = _transition(direction)
    if transition_key not in adata.uns.keys():
        raise ValueError("Please run `cellrank.tl.root_final` first.")

    start = logg.info(f"Computing lineage probabilities towards `{rc_key}`")

    # get the transition matrix from the AnnData object and initialise MC object
    # TODO: should kernel object read the transition matrix from AnnData?
    vk = VelocityKernel(adata, backward=not final)
    vk.transition_matrix = adata.uns[transition_key]["T"]
    mc = MarkovChain(vk)

    # compute the absoprtion probabilities
    mc.compute_lin_probs(keys=keys)

    logg.info(f"Added key `{lin_key!r}` to `adata.obsm`\n" f"    Finish", time=start)

    return adata if copy else None
