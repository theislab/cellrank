from typing import Any, Optional

from anndata import AnnData
from cellrank._key import Key
from cellrank.ul._docs import d
from cellrank.tl._utils import _deprecate
from cellrank.tl.estimators import CFLARE


@_deprecate(version="2.0")
@d.dedent
def lineages(
    adata: AnnData,
    backward: bool = False,
    copy: bool = False,
    return_estimator: bool = False,
    **kwargs: Any,
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
        Whether to return the estimator. Only available when ``copy = False``.
    kwargs
        Keyword arguments for :meth:`cellrank.tl.estimators.BaseEstimator.compute_absorption_probabilities`.

    Returns
    -------
    Depending on ``copy`` and ``return_estimator``, either updates the existing ``adata`` object,
    returns its copy or returns the estimator.
    """
    if copy:
        adata = adata.copy()

    try:
        mc = CFLARE.from_adata(adata, obsp_key=Key.uns.kernel(backward))
        assert mc.adata is adata
    except KeyError as e:
        raise RuntimeError(
            f"Compute transition matrix first as `cellrank.tl.transition_matrix(..., backward={backward})`."
        ) from e

    if mc.terminal_states is None:
        fs_key = Key.obs.term_states(backward)
        raise RuntimeError(
            f"Compute the states first as `cellrank.tl.{fs_key}(..., backward={backward})`."
        )

    # compute the absorption probabilities
    mc.compute_absorption_probabilities(**kwargs)

    return mc.adata if copy else mc if return_estimator else None
