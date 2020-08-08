# -*- coding: utf-8 -*-
"""Module used for finding root and final states."""

from typing import Union, TypeVar, Optional

from cellrank.utils._docs import d, _root, _final, inject_docs
from cellrank.tools._utils import (
    _check_estimator_type,
    _info_if_obs_keys_categorical_present,
)
from cellrank.tools.estimators import GPCCA, CFLARE
from cellrank.tools._transition_matrix import transition_matrix
from cellrank.tools.estimators._constants import P
from cellrank.tools.estimators._base_estimator import BaseEstimator

AnnData = TypeVar("AnnData")


_find_docs = """\
Compute {direction} states based on RNA velocity, see [Manno18]_. The tool models dynamic cellular
processes as a Markov chain, where the transition matrix is computed based on the velocity vectors of each
individual cell. Based on this Markov chain, we provide two estimators to compute {direction} states, both of which
are based on spectral methods.

For the estimator :class:`cellrank.tl.GPCCA`, cells are fuzzily clustered into metastable states,
using Generalized Perron Cluster Cluster Analysis [GPCCA18]_.
In short, this coarse-grains the Markov chain into a set of macrostates representing the slow
time-scale dynamics, i.e. transitions between these macrostates are rare. The most stable ones of these will represent
{direction}, while the others will represent transient, metastable states.

For the estimator :class:`cellrank.tl.CFLARE`, cells are filtered into transient/recurrent cells using the
left eigenvectors of the transition matrix and clustered into distinct groups of {direction} states using the right
eigenvectors of the transition matrix of the Markov chain.

Parameters
----------
%(adata)s
estimator
    Estimator to use to compute the {direction} states.
mode
    How to compute transition probabilities. Options are "stochastic" (propagate uncertainty analytically),
    "deterministic" (don't propagate uncertainty) and "sampling" (sample from velocity distribution).
n_states
    If you know how many {direction} states you are expecting, you can provide this number.
    Otherwise, an `eigengap` heuristic is used.
cluster_key
    Key from :paramref:`adata` `.obs` where cluster annotations are stored.
    These are used to give names to the {direction} states.
weight_connectivities
    Weight given to a transition matrix computed on the basis of the KNN connectivities. Must be in `[0, 1]`.
    This can help in situations where we have noisy velocities and want to give some weight to
    transcriptomic similarity.
show_plots
    Whether to show plots of the spectrum and eigenvectors in the embedding.
copy
    Whether to update the existing :paramref:`adata` object or to return a copy.
return_estimator
    Whether to return the estimator. Only available when :paramref:`copy=False`.
**kwargs
    Keyword arguments for :meth:`cellrank.tl.BaseEstimator.fit`, such as `n_cells`.

Returns
-------
:class:`anndata.AnnData`, :class:`cellrank.tools.estimators.BaseEstimator` or :class:`NoneType`
    Depending on :paramref:`copy`, either updates the existing :paramref:`adata` object, returns a copy or
    returns the estimator.

    Marked cells can be found in :paramref:`adata` `.obs[`{key_added!r}]``.
"""


def _root_final(
    adata: AnnData,
    estimator: type(BaseEstimator) = GPCCA,
    backward: bool = False,
    mode: str = "deterministic",
    n_states: Optional[int] = None,
    cluster_key: Optional[str] = None,
    weight_connectivities: float = None,
    show_plots: bool = False,
    copy: bool = False,
    return_estimator: bool = False,
    **kwargs,
) -> Optional[Union[AnnData, BaseEstimator]]:

    _check_estimator_type(estimator)
    adata = adata.copy() if copy else adata

    # compute kernel object
    kernel = transition_matrix(
        adata,
        backward=backward,
        mode=mode,
        weight_connectivities=weight_connectivities,
    )
    # create estimator object
    mc = estimator(kernel, read_from_adata=False)

    if cluster_key is None:
        _info_if_obs_keys_categorical_present(
            adata,
            keys=["louvain", "clusters"],
            msg_fmt="Found categorical observation in `adata.obs[{!r}]`. "
            "Consider specifying it as `cluster_key`.",
        )

    mc.fit(
        n_lineages=n_states,
        cluster_key=cluster_key,
        compute_absorption_probabilities=False,
        **kwargs,
    )

    if show_plots:
        mc.plot_spectrum(real_only=True)
        if isinstance(mc, CFLARE):
            mc.plot_eigendecomposition(abs_value=True, perc=[0, 98], use=n_states)
            mc.plot_final_states(discrete=True, same_plot=False)
        elif isinstance(mc, GPCCA):
            n_states = len(mc._get(P.META).cat.categories)
            if n_states > 1:
                mc.plot_schur()
            mc.plot_final_states(discrete=True, same_plot=False)
            if n_states > 1:
                mc.plot_coarse_T()
        else:
            raise NotImplementedError(
                f"Pipeline not implemented for `{type(mc).__name__!r}.`"
            )

    return adata if copy else mc if return_estimator else None


@d.dedent
@inject_docs(root=_find_docs.format(direction=_root, key_added=f"{_root}_states",))
def root_states(
    adata: AnnData,
    estimator: type(BaseEstimator) = GPCCA,
    mode: str = "deterministic",
    n_states: Optional[int] = None,
    cluster_key: Optional[str] = None,
    weight_connectivities: float = None,
    show_plots: bool = False,
    copy: bool = False,
    return_estimator: bool = False,
    **kwargs,
) -> Optional[AnnData]:
    """
    Find %(root)s states of a dynamic process of single cells.

    {root}
    """

    return _root_final(
        adata,
        estimator=estimator,
        mode=mode,
        backward=True,
        n_states=n_states,
        cluster_key=cluster_key,
        weight_connectivities=weight_connectivities,
        show_plots=show_plots,
        copy=copy,
        return_estimator=return_estimator,
        **kwargs,
    )


@d.dedent
@inject_docs(final=_find_docs.format(direction=_final, key_added=f"{_final}_states",))
def final_states(
    adata: AnnData,
    estimator: type(BaseEstimator) = GPCCA,
    mode: str = "deterministic",
    n_states: Optional[int] = None,
    cluster_key: Optional[str] = None,
    weight_connectivities: float = None,
    show_plots: bool = False,
    copy: bool = False,
    return_estimator: bool = False,
    **kwargs,
) -> Optional[AnnData]:
    """
    Find %(final)s states of a dynamic process of single cells.

    {final}
    """

    return _root_final(
        adata,
        estimator=estimator,
        mode=mode,
        backward=False,
        n_states=n_states,
        cluster_key=cluster_key,
        weight_connectivities=weight_connectivities,
        show_plots=show_plots,
        copy=copy,
        return_estimator=return_estimator,
        **kwargs,
    )
