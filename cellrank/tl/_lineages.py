"""Lineages module."""

from typing import Union, TypeVar, Optional, Sequence

import pandas as pd

from cellrank import logging as logg
from cellrank.ul._docs import d
from cellrank.tl._utils import TestMethod
from cellrank.tl.kernels import PrecomputedKernel
from cellrank.tl._constants import AbsProbKey, TermStatesKey, TerminalStatesPlot
from cellrank.tl.estimators import GPCCA
from cellrank.tl.estimators._constants import P
from cellrank.tl.kernels._precomputed_kernel import DummyKernel

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

    This function computes the absorption probabilities of a Markov chain towards the %(initial_or_terminal) states
    uncovered by :func:`cellrank.tl.initial_states` or :func:`cellrank.tl.terminal_states` using a highly efficient
    implementation that scales to large cell numbers.

    It's also possible to calculate mean and variance of the time until absorption for all or just a subset
    of the %(initial_or_terminal)s states. This can be seen as a pseudotemporal measure, either towards any terminal
    population of the state change trajectory, or towards specific ones.

    Parameters
    ----------
    %(adata)s
    %(backward)s
    copy
        Whether to update the existing ``adata`` object or to return a copy.
    return_estimator
        Whether to return the estimator. Only available when ``copy=False``.
    kwargs
        Keyword arguments for :meth:`cellrank.tl.estimators.BaseEstimator.compute_absorption_probabilities`.

    Returns
    -------
    :class:`anndata.AnnData`, :class:`cellrank.tl.estimators.BaseEstimator` or :obj:`None`
        Depending on ``copy`` and ``return_estimator``, either updates the existing ``adata`` object,
        returns its copy or returns the estimator.
    """

    if backward:
        lin_key = AbsProbKey.BACKWARD
        fs_key = TermStatesKey.BACKWARD
        fs_key_pretty = TerminalStatesPlot.BACKWARD
    else:
        lin_key = AbsProbKey.FORWARD
        fs_key = TermStatesKey.FORWARD
        fs_key_pretty = TerminalStatesPlot.FORWARD

    try:
        pk = PrecomputedKernel(adata=adata, backward=backward)
    except KeyError as e:
        raise RuntimeError(
            f"Compute transition matrix first as `cellrank.tl.transition_matrix(..., backward={backward})`."
        ) from e

    start = logg.info(f"Computing lineage probabilities towards {fs_key_pretty.s}")
    mc = GPCCA(
        pk, read_from_adata=True, inplace=not copy
    )  # GPCCA is more general than CFLARE, in terms of what is saves
    if mc._get(P.TERM) is None:
        raise RuntimeError(
            f"Compute the states first as `cellrank.tl.{fs_key.s}(..., backward={backward})`."
        )

    # compute the absorption probabilities
    mc.compute_absorption_probabilities(**kwargs)

    logg.info(f"Adding lineages to `adata.obsm[{lin_key.s!r}]`\n    Finish", time=start)

    return mc.adata if copy else mc if return_estimator else None


@d.dedent
def lineage_drivers(
    adata: AnnData,
    backward: bool = False,
    lineages: Optional[Union[Sequence, str]] = None,
    method: str = TestMethod.FISCHER.s,
    cluster_key: Optional[str] = None,
    clusters: Optional[Union[Sequence, str]] = None,
    layer: str = "X",
    use_raw: bool = False,
    confidence_level: float = 0.95,
    n_perms: int = 1000,
    seed: Optional[int] = None,
    return_drivers: bool = True,
    **kwargs,
) -> Optional[pd.DataFrame]:
    """
    %(lineage_drivers.full_desc)s

    Parameters
    ----------
    %(adata)s
    %(backward)s
    %(lineage_drivers.parameters)s

    Returns
    -------
    %(lineage_drivers.returns)s

    References
    ----------
    %(lineage_drivers.references)s
    """  # noqa: D400

    # create dummy kernel and estimator
    pk = DummyKernel(adata, backward=backward)
    g = GPCCA(pk, read_from_adata=True, write_to_adata=False)
    if g._get(P.ABS_PROBS) is None:
        raise RuntimeError(
            f"Compute absorption probabilities first as `cellrank.tl.lineages(..., backward={backward})`."
        )

    # call the underlying function to compute and store the lineage drivers
    return g.compute_lineage_drivers(
        method=method,
        lineages=lineages,
        cluster_key=cluster_key,
        clusters=clusters,
        layer=layer,
        use_raw=use_raw,
        confidence_level=confidence_level,
        n_perms=n_perms,
        seed=seed,
        return_drivers=return_drivers,
        **kwargs,
    )
