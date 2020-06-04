# -*- coding: utf-8 -*-
"""Lineages module."""

from typing import Optional, Sequence

from scanpy import logging as logg
from anndata import AnnData

from cellrank.tools.kernels import VelocityKernel
from cellrank.tools._constants import LinKey, StateKey, Direction, _transition
from cellrank.tools.estimators import GPCCA, CFLARE
from cellrank.tools.estimators._base_estimator import BaseEstimator


def lineages(
    adata: AnnData,
    estimator: type(BaseEstimator) = GPCCA,
    final: bool = True,
    keys: Optional[Sequence[str]] = None,
    copy: bool = False,
    method: str = "krylov",
    **kwargs,
) -> Optional[AnnData]:
    """
    Compute probabilistic lineage assignment using RNA velocity.

    For each cell i in {1, ..., N} and start/endpoint j in {1, ..., M}, the probability is computed that cell i
    is either going to j (end point) or coming from j (start point). Mathematically, this computes absorption
    probabilities to approximate recurrent classes using an RNA velocity based Markov chain.

    Note that absorption probabilities have been used in the single cell context to infer lineage probabilities e.g.
    in [Setty19]_ or [Weinreb18]_ and we took inspiration from there.

    Before running this function, compute start/endpoints using :func:`cellrank.tl.root_final`.

    Parameters
    --------
    adata : :class:`anndata.AnnData`
        Annotated data object
    estimator
        Estimator to use to compute the lineage probabilities.
    final
        If `True`, computes final cells, i.e. end points. Otherwise, computes root cells, i.e. starting points.
    keys
        Determines which end/start-points to use by passing their names. Further, start/end-points can be combined.
        If e.g. the endpoints are ['Neuronal_1', 'Neuronal_1', 'Astrocytes', 'OPC'], then passing
        keys=['Neuronal_1, Neuronal_2', 'OPC'] means that the two neuronal endpoints are treated as one and
        Astrocytes are excluded.
    method
        Method to compute Schur vectors when :paramref:`estimator` is :class:`cellrank.tl.GPCCA`.
    copy
        Whether to update the existing AnnData object or to return a copy.
    kwargs
        Keyword arguments for :meth:`cellrank.tl.estimators.BaseEstimator.compute_metastable_states`.

    Returns
    --------
    :class:`anndata.AnnData` or :class:`NoneType`
        Depending on :paramref:`copy`, either updates the existing :paramref:`adata` object or returns a copy.
    """

    if not isinstance(estimator, type):
        raise TypeError(
            f"Expected estimator to be a class, found `{type(estimator).__name__}`."
        )

    if not issubclass(estimator, BaseEstimator):
        raise TypeError(
            f"Expected estimator to be a subclass of `BaseEstimator`, found `{type(estimator).__name__}`"
        )

    # Set the keys and print info
    adata = adata.copy() if copy else adata

    if final:
        direction = Direction.FORWARD
        lin_key = LinKey.FORWARD
        rc_key = StateKey.FORWARD
    else:
        direction = Direction.BACKWARD
        lin_key = LinKey.BACKWARD
        rc_key = StateKey.BACKWARD

    transition_key = _transition(direction)
    if transition_key not in adata.uns.keys():
        key = "final" if final else "root"
        raise ValueError(f"Compute {key} cells first as `cellrank.tl.find_{key}`.")

    start = logg.info(f"Computing lineage probabilities towards `{rc_key}`")

    # get the transition matrix from the AnnData object and initialise MC object
    vk = VelocityKernel(adata, backward=not final)
    vk.transition_matrix = adata.uns[transition_key]["T"]
    mc = estimator(vk, read_from_adata=False)

    # compute the absorption probabilities
    if isinstance(mc, CFLARE):
        mc.compute_metastable_states(**kwargs)
        mc.compute_lin_probs(keys=keys)
    elif isinstance(mc, GPCCA):
        n_states = kwargs["n_states"]
        if n_states == 1:
            mc.compute_eig()
        mc.compute_schur(n_components=n_states, method=method)
        mc.compute_metastable_states(**kwargs)
        mc.set_main_states()
    else:
        raise NotImplementedError(
            f"Pipeline not implemented for `{type(bytes).__name__}`"
        )

    logg.info(f"Added key `{lin_key!r}` to `adata.obsm`\n    Finish", time=start)

    return adata if copy else None
