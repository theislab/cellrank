# -*- coding: utf-8 -*-
"""Module used for finding initial and terminal states."""
from types import MappingProxyType
from typing import Union, Mapping, TypeVar, Optional

from cellrank import logging as logg
from cellrank.ul._docs import d, _initial, _terminal, inject_docs
from cellrank.tl._utils import (
    _check_estimator_type,
    _info_if_obs_keys_categorical_present,
)
from cellrank.tl.kernels import PrecomputedKernel
from cellrank.tl._constants import TermStatesKey
from cellrank.tl.estimators import GPCCA, CFLARE
from cellrank.tl._transition_matrix import transition_matrix
from cellrank.tl.estimators._constants import P
from cellrank.tl.kernels._velocity_kernel import BackwardMode, VelocityMode
from cellrank.tl.estimators._base_estimator import BaseEstimator

AnnData = TypeVar("AnnData")


_docstring = """\
Find {direction} states of a dynamic process of single cells based on RNA velocity [Manno18]_.

The function models dynamic cellular processes as a Markov chain, where the transition matrix is computed based
on the velocity vectors of each individual cell. Based on this Markov chain, we provide two estimators
to compute {direction} states, both of which are based on spectral methods.

For the estimator :class:`cellrank.tl.estimators.GPCCA`, cells are fuzzily clustered into macrostates,
using Generalized Perron Cluster Cluster Analysis [GPCCA18]_.
In short, this coarse-grains the Markov chain into a set of macrostates representing the slow
time-scale dynamics, i.e. transitions between these macrostates are rare. The most stable ones of these will represent
{direction}, while the others represent intermediate macrostates.

For the estimator :class:`cellrank.tl.estimators.CFLARE`, cells are filtered into transient/recurrent cells using the
left eigenvectors of the transition matrix and clustered into distinct groups of {direction} states using the right
eigenvectors of the transition matrix of the Markov chain.

Parameters
----------
%(adata)s
estimator
    Estimator class to use to compute the {direction} states.
%(velocity_mode)s{bwd_mode}
n_states
    If you know how many {direction} states you are expecting, you can provide this number.
    Otherwise, an `eigengap` heuristic is used.
cluster_key
    Key from ``adata.obs`` where cluster annotations are stored. These are used to give names to the {direction} states.
key
    Key in ``adata.obsp`` where the transition matrix is saved.
    If not found, compute a new one using :func:`cellrank.tl.transition_matrix`.
weight_connectivities
    Weight given to a transition matrix computed on the basis of the KNN connectivities. Must be in `[0, 1]`.
    This can help in situations where we have noisy velocities and want to give some weight to transcriptomic
    similarity.
show_plots
    Whether to show plots of the spectrum and eigenvectors in the embedding.
%(n_jobs)s
copy
    Whether to update the existing ``adata`` object or to return a copy.
return_estimator
    Whether to return the estimator. Only available when ``copy=False``.
fit_kwargs
    Keyword arguments for :meth:`cellrank.tl.BaseEstimator.fit`, such as ``n_cells``.
**kwargs
    Keyword arguments for :func:`cellrank.tl.transition_matrix`, such as ``weight_connectivities`` or ``softmax_scale``.

Returns
-------
:class:`anndata.AnnData`, :class:`cellrank.tl.estimators.BaseEstimator` or :obj:`None`
    Depending on ``copy`` and ``return_estimator``, either updates the existing ``adata`` object,
    returns its copy or returns the estimator.

    Marked cells are added to ``adata.obs[{key_added!r}]``.
"""


def _initial_terminal(
    adata: AnnData,
    estimator: type(BaseEstimator) = GPCCA,
    backward: bool = False,
    mode: str = VelocityMode.DETERMINISTIC.s,
    backward_mode: str = BackwardMode.TRANSPOSE.s,
    n_states: Optional[int] = None,
    cluster_key: Optional[str] = None,
    key: Optional[str] = None,
    show_plots: bool = False,
    copy: bool = False,
    return_estimator: bool = False,
    fit_kwargs: Mapping = MappingProxyType({}),
    **kwargs,
) -> Optional[Union[AnnData, BaseEstimator]]:

    _check_estimator_type(estimator)

    try:
        kernel = PrecomputedKernel(key, adata=adata, backward=backward)
        write_to_adata = False  # no need to write
        logg.info("Using precomputed transition matrix")
    except KeyError:
        # compute kernel object
        kernel = transition_matrix(
            adata,
            backward=backward,
            mode=mode,
            backward_mode=backward_mode,
            **kwargs,
        )
        write_to_adata = True

    # create estimator object
    mc = estimator(
        kernel,
        read_from_adata=False,
        inplace=not copy,
        key=key,
        write_to_adata=write_to_adata,
    )

    if cluster_key is None:
        _info_if_obs_keys_categorical_present(
            adata,
            keys=["louvain", "leiden", "clusters"],
            msg_fmt="Found categorical observation in `adata.obs[{!r}]`. Consider specifying it as `cluster_key`.",
        )

    mc.fit(
        n_lineages=n_states,
        cluster_key=cluster_key,
        compute_absorption_probabilities=False,
        **fit_kwargs,
    )

    if show_plots:
        mc.plot_spectrum(real_only=True)
        if isinstance(mc, CFLARE):
            mc.plot_eigendecomposition(abs_value=True, perc=[0, 98], use=n_states)
            mc.plot_terminal_states(discrete=True, same_plot=False)
        elif isinstance(mc, GPCCA):
            n_states = len(mc._get(P.MACRO).cat.categories)
            if n_states > 1:
                mc.plot_schur()
            mc.plot_terminal_states(discrete=True, same_plot=False)
            if n_states > 1:
                mc.plot_coarse_T()
        else:
            raise NotImplementedError(
                f"Pipeline not implemented for `{type(mc).__name__!r}.`"
            )

    return mc.adata if copy else mc if return_estimator else None


@inject_docs(m=VelocityMode, b=BackwardMode)
@d.dedent
@inject_docs(
    __doc__=_docstring.format(
        direction=_initial,
        key_added=TermStatesKey.BACKWARD.s,
        bwd_mode="\n%(velocity_backward_mode_high_lvl)s",
    )
)
def initial_states(
    adata: AnnData,
    estimator: type(BaseEstimator) = GPCCA,
    mode: str = VelocityMode.DETERMINISTIC.s,
    backward_mode: str = BackwardMode.TRANSPOSE.s,
    n_states: Optional[int] = None,
    cluster_key: Optional[str] = None,
    key: Optional[str] = None,
    show_plots: bool = False,
    copy: bool = False,
    return_estimator: bool = False,
    fit_kwargs: Mapping = MappingProxyType({}),
    **kwargs,
) -> Optional[AnnData]:  # noqa

    return _initial_terminal(
        adata,
        estimator=estimator,
        mode=mode,
        backward_mode=backward_mode,
        backward=True,
        n_states=n_states,
        cluster_key=cluster_key,
        key=key,
        show_plots=show_plots,
        copy=copy,
        return_estimator=return_estimator,
        fit_kwargs=fit_kwargs,
        **kwargs,
    )


@inject_docs(m=VelocityMode, b=BackwardMode)
@d.dedent
@inject_docs(
    __doc__=_docstring.format(
        direction=_terminal, key_added=TermStatesKey.FORWARD.s, bwd_mode=""
    )
)
def terminal_states(
    adata: AnnData,
    estimator: type(BaseEstimator) = GPCCA,
    mode: str = VelocityMode.DETERMINISTIC.s,
    n_states: Optional[int] = None,
    cluster_key: Optional[str] = None,
    key: Optional[str] = None,
    show_plots: bool = False,
    copy: bool = False,
    return_estimator: bool = False,
    fit_kwargs: Mapping = MappingProxyType({}),
    **kwargs,
) -> Optional[AnnData]:  # noqa

    return _initial_terminal(
        adata,
        estimator=estimator,
        mode=mode,
        backward=False,
        n_states=n_states,
        cluster_key=cluster_key,
        key=key,
        show_plots=show_plots,
        copy=copy,
        return_estimator=return_estimator,
        fit_kwargs=fit_kwargs,
        **kwargs,
    )
