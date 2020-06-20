# -*- coding: utf-8 -*-
"""Module used for finding root and final states."""

from typing import Union, Optional

from scanpy import logging as logg
from anndata import AnnData

from cellrank.utils._docs import inject_docs
from cellrank.tools._utils import _info_if_obs_keys_categorical_present
from cellrank.tools._constants import StateKey
from cellrank.tools.estimators import GPCCA, CFLARE
from cellrank.tools._transition_matrix import transition_matrix
from cellrank.tools.estimators._base_estimator import BaseEstimator

_find_docs = """\
Compute {cells} states based on RNA velocity, see [Manno18]_.The tool models dynamic cellular
processes as a Markov chain, where the transition matrix is computed based on the velocity vectors of each
individual cell. Based on this Markov chain, we provide two estimators to compute {cells} states, both of which
are based on spectral methods.

For the estimator GPCCA, cells are fuzzily clustered into metastable states, using Generalized Perron Cluster Cluster
Analysis [GPCCA18]_. In short, this coarse-grains the Markov chain into a set of macrostates representing the slow
time-scale dynamics, i.e. transitions between these macrostates are rare. The most stable ones of these will represent
{cells}, while the others will represent transient, metastable states.

For the estimator CFLARE, cells are filtered into transient/recurrent cells using the left eigenvectors of the
transition matrix and clustered into distinct groups of {cells} states using the right eigenvectors of the transition
matrix of the Markov chain.

Params
------
adata : :class:`adata.AnnData`
    Annotated data object.
estimator
    Estimator to use to compute the {cells} states.
n_states
    If you know how many {direction} states you are expecting, you can provide this number.
    Otherwise, an `eigen-gap` heuristic is used.
cluster_key
    Key from `adata.obs` where cluster annotations are stored. These are used to give names to the {direction} states.
weight_connectivities
    Weight given to a transition matrix computed on the basis of the KNN connectivities. Should be in `[0, 1]`. This
    can help in situations where we have noisy velocities and want to give some weight to transcriptomic similarity.
use_velocity_uncertainty
    Whether to use velocity uncertainty. Uncertainties are computed independently per gene using the neighborhood graph.
    They are then propagated into cosine similarities and finally used as a scaling factor in the softmax which
    transforms cosine similarities to probabilities, i.e. transitions we are uncertain about are down-weighted.
method
    Method to use when computing the Schur decomposition. Only needed when :paramref:`estimator`
    is :class`:cellrank.tl.GPCCA:.
    Valid options are: `'krylov'`, `'brandts'`.
show_plots
    Whether to show plots of the spectrum and eigenvectors in the embedding.
copy
    Whether to update the existing :paramref:`adata` object or to return a copy.
return_estimator
    Whether to return the estimator. Only available when :paramref:`copy=False`.
kwargs
    Keyword arguments for :meth:`cellrank.tl.estimators.BaseEstimator.compute_metastable_states`.

Returns
-------

:class:`anndata.AnnData`, :class:`cellrank.tools.estimators.BaseEstimator` or :class:`NoneType`
    Depending on :paramref:`copy`, either updates the existing :paramref:`adata` object or returns a copy or
    returns the estimator.
    Marked cells can be found in :paramref:`adata` `.obs` under `{key_added!r}`.
"""


def _root_final(
    adata: AnnData,
    estimator: type(BaseEstimator) = GPCCA,
    final: bool = True,
    n_states: Optional[int] = None,
    cluster_key: Optional[str] = None,
    weight_connectivities: float = None,
    use_velocity_uncertainty: bool = False,
    method: str = "krylov",
    show_plots: bool = False,
    copy: bool = False,
    return_estimator: bool = False,
    **kwargs,
) -> Optional[Union[AnnData, BaseEstimator]]:
    """Compute root or final states of  Markov Chain."""

    key = StateKey.FORWARD if final else StateKey.BACKWARD
    logg.info(f"Computing `{key}`")
    adata = adata.copy() if copy else adata

    # compute kernel object
    kernel = transition_matrix(
        adata,
        backward=not final,
        weight_connectivities=weight_connectivities,
        scale_by_variances=use_velocity_uncertainty,
    )
    # create MarkovChain object
    mc = estimator(kernel, read_from_adata=False)

    if cluster_key is None:
        _info_if_obs_keys_categorical_present(
            adata,
            keys=["louvain", "clusters"],
            msg_fmt="Found categorical observation in `adata.obs[{!r}]`. "
            "Consider specifying it as `cluster_key`.",
        )

    # run the computation
    if isinstance(mc, CFLARE):
        kwargs["use"] = n_states

        mc.compute_eig()
        mc.compute_metastable_states(cluster_key=cluster_key, **kwargs)

        if show_plots:
            mc.plot_spectrum(real_only=True)
            mc.plot_eig_embedding(abs_value=True, perc=[0, 98], use=n_states)
            mc.plot_eig_embedding(left=False, use=n_states)

    elif isinstance(mc, GPCCA):
        if n_states is None or n_states == 1:
            mc.compute_eig()
            if n_states is None:
                n_states = mc.eigendecomposition["eigengap"] + 1
        if n_states > 1:
            mc.compute_schur(n_states + 1, method=method)

        try:
            mc.compute_metastable_states(
                n_states=n_states, cluster_key=cluster_key, **kwargs
            )
        except ValueError:
            logg.warning(
                f"Computing {n_states} metastable states cuts through a block of complex conjugates. "
                f"Increasing `n_states` to {n_states + 1}"
            )
            mc.compute_metastable_states(
                n_states=n_states + 1, cluster_key=cluster_key, **kwargs
            )
        mc.set_main_states()  # write to adata

        if show_plots:
            mc.plot_spectrum(real_only=True)
            if n_states > 1:
                mc.plot_schur_embedding()
            mc.plot_metastable_states(same_plot=False)
            if n_states > 1:
                mc.plot_coarse_T()
    else:
        raise NotImplementedError(
            f"Pipeline not implemented for `{type(bytes).__name__}`"
        )

    return adata if copy else mc if return_estimator else None


@inject_docs(
    root=_find_docs.format(cells="root", direction="start", key_added="root_states")
)
def root_states(
    adata: AnnData,
    estimator: type(BaseEstimator) = GPCCA,
    n_states: Optional[int] = None,
    weight_connectivities: float = None,
    method: str = "krylov",
    show_plots: bool = False,
    copy: bool = False,
    return_estimator: bool = False,
    **kwargs,
) -> Optional[AnnData]:
    """
    Find root states of a dynamic process of single cells.

    {root}
    """

    return _root_final(
        adata,
        estimator=estimator,
        final=False,
        n_states=n_states,
        weight_connectivities=weight_connectivities,
        method=method,
        show_plots=show_plots,
        copy=copy,
        return_estimator=return_estimator,
        **kwargs,
    )


@inject_docs(
    final=_find_docs.format(cells="final", direction="end", key_added="final_states")
)
def final_states(
    adata: AnnData,
    estimator: type(BaseEstimator) = GPCCA,
    n_states: Optional[int] = None,
    weight_connectivities: float = None,
    method: str = "krylov",
    show_plots: bool = False,
    copy: bool = False,
    return_estimator: bool = False,
    **kwargs,
) -> Optional[AnnData]:
    """
    Find final states of a dynamic process of single cells.

    {final}
    """

    return _root_final(
        adata,
        estimator=estimator,
        final=True,
        n_states=n_states,
        weight_connectivities=weight_connectivities,
        method=method,
        show_plots=show_plots,
        copy=copy,
        return_estimator=return_estimator,
        **kwargs,
    )
