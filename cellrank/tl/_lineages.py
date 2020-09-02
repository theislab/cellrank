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

    This function computes the absorption probabilities of a Markov chain towards %(initial_or_terminal) states
    uncovered by :func:`cellrank.tl.initial_states` or :func:`cellrank.tl.terminal_states`.

    It's also possible to calculate mean and variance of time until absorption all or just a subset
    of %(initial_or_terminal)s states

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
        Depending on ``copy`` and ``return_estimator``, either updates the existing ``adata`` object,
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
            f"Compute {'backward' if backward else 'forward'} transition matrix first as "
            f"`cellrank.tl.transition_matrix(..., backward={backward})`."
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
