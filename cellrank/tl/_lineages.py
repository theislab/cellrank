# -*- coding: utf-8 -*-
"""Lineages module."""

from typing import TypeVar, Optional

from cellrank import logging as logg
from cellrank.ul._docs import d
from cellrank.tl._utils import _check_estimator_type
from cellrank.tl.kernels import PrecomputedKernel
from cellrank.tl._constants import AbsProbKey, FinalStatesKey
from cellrank.tl.estimators import GPCCA
from cellrank.tl.estimators._constants import P
from cellrank.tl.estimators._base_estimator import BaseEstimator

AnnData = TypeVar("AnnData")


@d.dedent
def lineages(
    adata: AnnData,
    estimator: BaseEstimator = GPCCA,
    backward: bool = False,
    copy: bool = False,
    return_estimator: bool = False,
    **kwargs,
) -> Optional[AnnData]:
    """
    Compute probabilistic lineage assignment using RNA velocity.

    For each cell `i` in :math:`{1, ..., N}` and %(root_or_final)s state `j` in :math:`{1, ..., M}`,
    the probability is computed that cell `i` is either going to %(final)s state `j` (``backward=False``)
    or is coming from %(root)s state `j` (`backward=True`).
    We provide two estimators for computing these probabilities:

    For the estimator :class:`cellrank.tl.estimators.GPCCA`, we perform Generalized Perron Cluster Cluster Analysis
    [GPCCA18]_. Cells are mapped to a simplex where each corner represents a %(root_or_final) state, and the position
    of a cell in the simplex determines its probability of going to a %(final)s states or coming from %(root)s states.

    For the estimator :class:`cellrank.tl.estimators.CFLARE`, we compute absorption probabilities towards
    the %(root_or_final)s states of the Markov chain.
    For related approaches in the single cell context that utilise absorption probabilities to map cells to lineages,
    see [Setty19]_ or [Weinreb18]_.

    Before running this function, compute %(root_or_final) states using :func:`cellrank.tl.root_states` or
    :func:`cellrank.tl.final_states`, respectively.

    Parameters
    ----------
    %(adata)s
    estimator
        Estimator class to use to compute the lineage probabilities.
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
        Depending on ``copy``, either updates the existing ``adata`` object or returns a copy or returns the estimator.
    """

    _check_estimator_type(estimator)
    adata = adata.copy() if copy else adata

    if backward:
        lin_key = AbsProbKey.BACKWARD
        rc_key = FinalStatesKey.BACKWARD
    else:
        lin_key = AbsProbKey.FORWARD
        rc_key = FinalStatesKey.FORWARD

    try:
        pk = PrecomputedKernel(adata=adata, backward=backward)
    except KeyError as e:
        raise RuntimeError(
            "Compute transition matrix first a `cellrank.tl.transition_matrix()`."
        ) from e

    start = logg.info(f"Computing lineage probabilities towards `{rc_key}`")
    mc = estimator(pk, read_from_adata=True)
    if mc._get(P.FIN) is None:
        raise RuntimeError(
            f"Compute the states first as `cellrank.tl.{'root' if backward else 'final'}_states()`."
        )

    # compute the absorption probabilities
    mc.compute_absorption_probabilities(**kwargs)

    logg.info(f"Adding lineages to `adata.obsm[{lin_key.s!r}]`\n    Finish", time=start)

    return adata if copy else mc if return_estimator else None
