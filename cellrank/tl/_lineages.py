# -*- coding: utf-8 -*-
"""Lineages module."""

from typing import TypeVar, Optional

from cellrank import logging as logg
from cellrank.ul._docs import d
from cellrank.tl.kernels import PrecomputedKernel
from cellrank.tl._constants import AbsProbKey, FinalStatesKey, FinalStatesPlot
from cellrank.tl.estimators import GPCCA
from cellrank.tl.estimators._constants import P

AnnData = TypeVar("AnnData")


@d.dedent
def lineages(
    adata: AnnData,
    backward: bool = False,
    copy: bool = False,
    return_estimator: bool = False,
    **kwargs,
) -> Optional[AnnData]:
    """
    Compute probabilistic lineage assignment using RNA velocity.

    For each cell `i` in :math:`{1, ..., N}` and %(initial_or_terminal)s state `j` in :math:`{1, ..., M}`,
    the probability is computed that cell `i` is either going to %(terminal)s state `j` (``backward=False``)
    or is coming from %(initial)s state `j` (``backward=True``).
    We provide two estimators for computing these probabilities:

    For the estimator :class:`cellrank.tl.estimators.GPCCA`, we perform Generalized Perron Cluster Cluster Analysis
    [GPCCA18]_. Cells are mapped to a simplex where each corner represents a %(initial_or_terminal) state, and the
    position of a cell in the simplex determines its probability of going to a %(terminal)s states or
    coming from %(initial)s states.

    For the estimator :class:`cellrank.tl.estimators.CFLARE`, we compute absorption probabilities towards
    the %(initial_or_terminal)s states of the Markov chain.
    For related approaches in the single cell context that utilise absorption probabilities to map cells to lineages,
    see [Setty19]_ or [Weinreb18]_.

    Before running this function, compute %(initial_or_terminal) states using :func:`cellrank.tl.initial_states` or
    :func:`cellrank.tl.terminal_states`, respectively.

    Parameters
    ----------
    %(adata)s
    %(backward)s
    copy
        Whether to update the existing ``adata`` object or to return a copy.
    return_estimator
        Whether to return the estimator. Only available when ``copy=False``.
    **kwargs
        Keyword arguments for :meth:`cellrank.tl.estimators.BaseEstimator.compute_absorption_probabilities`.

    Returns
    -------
    :class:`anndata.AnnData`, :class:`cellrank.tl.estimators.BaseEstimator` or :obj:`None`
        Depending on ``copy`` and ``return_estimato``, either updates the existing ``adata`` object,
        returns its copy or returns the estimator.
    """

    if backward:
        lin_key = AbsProbKey.BACKWARD
        fs_key = FinalStatesKey.BACKWARD
        fs_key_pretty = FinalStatesPlot.BACKWARD
    else:
        lin_key = AbsProbKey.FORWARD
        fs_key = FinalStatesKey.FORWARD
        fs_key_pretty = FinalStatesPlot.FORWARD
    try:
        pk = PrecomputedKernel(adata=adata, backward=backward)
    except KeyError as e:
        raise RuntimeError(
            "Compute transition matrix first a `cellrank.tl.transition_matrix()`."
        ) from e

    start = logg.info(f"Computing lineage probabilities towards {fs_key_pretty.s}")
    mc = GPCCA(
        pk, read_from_adata=True, inplace=not copy
    )  # GPCCA is more general than CFLARE
    if mc._get(P.FIN) is None:
        raise RuntimeError(f"Compute the states first as `cellrank.tl.{fs_key.s}()`.")

    # compute the absorption probabilities
    mc.compute_absorption_probabilities(**kwargs)

    logg.info(f"Adding lineages to `adata.obsm[{lin_key.s!r}]`\n    Finish", time=start)

    return mc.adata if copy else mc if return_estimator else None
