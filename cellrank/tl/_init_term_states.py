from typing import Any, Union, Mapping, Optional, Sequence

from types import MappingProxyType

from anndata import AnnData
from cellrank import logging as logg
from cellrank.ul._docs import d, _initial, _terminal, inject_docs
from cellrank.tl._utils import (
    _deprecate,
    _check_estimator_type,
    _info_if_obs_keys_categorical_present,
)
from cellrank.tl.kernels import PrecomputedKernel
from cellrank.tl.estimators import GPCCA, CFLARE, TermStatesEstimator
from cellrank.tl._transition_matrix import transition_matrix
from cellrank.tl.kernels._velocity_kernel import BackwardMode, VelocityMode
from cellrank.tl.estimators._base_estimator import BaseEstimator

_docstring = """\
Find {direction} states of a dynamic process of single cells based on RNA velocity :cite:`manno:18`.

The function models dynamic cellular processes as a Markov chain, where the transition matrix is computed based
on the velocity vector of each individual cell. Based on this Markov chain, we provide two estimators
to compute {direction} states, both of which are based on spectral methods.

For the estimator :class:`cellrank.tl.estimators.GPCCA`, cells are fuzzily clustered into macrostates,
using Generalized Perron Cluster Cluster Analysis :cite:`reuter:18`. In short, this coarse-grains the Markov chain into
a set of macrostates representing the slow time-scale dynamics, i.e. transitions between these macrostates are rare.
The most stable ones of these will represent {direction}, while the others represent intermediate macrostates.

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
force_recompute
    Whether to always recompute the transition matrix even if one exists.
show_plots
    Whether to show plots of the spectrum and eigenvectors in the embedding.
%(n_jobs)s
copy
    Whether to update the existing ``adata`` object or to return a copy.
return_estimator
    Whether to return the estimator. Only available when ``copy=False``.
fit_kwargs
    Keyword arguments for :meth:`cellrank.tl.BaseEstimator.fit`, such as ``n_cells``.
kwargs
    Keyword arguments for :func:`cellrank.tl.transition_matrix`, such as ``weight_connectivities`` or ``softmax_scale``.

Returns
-------
Depending on ``copy`` and ``return_estimator``, either updates the existing ``adata`` object,
returns its copy or returns the estimator.

Marked cells are added to ``adata.obs[{key_added!r}]``.
"""


def _fit(
    estim: TermStatesEstimator,
    n_lineages: Optional[int] = None,
    keys: Optional[Sequence[str]] = None,
    cluster_key: Optional[str] = None,
    compute_absorption_probabilities: bool = True,
    **kwargs: Any,
) -> TermStatesEstimator:
    if isinstance(estim, CFLARE):
        estim.compute_eigendecomposition(k=20 if n_lineages is None else n_lineages + 1)
        if n_lineages is None:
            n_lineages = estim.eigendecomposition["eigengap"] + 1

        estim.compute_terminal_states(
            use=n_lineages,
            cluster_key=cluster_key,
            n_clusters_kmeans=n_lineages,
            method=kwargs.pop("method", "kmeans"),
            **kwargs,
        )
    elif isinstance(estim, GPCCA):
        if n_lineages is None or n_lineages == 1:
            estim.compute_eigendecomposition()
            if n_lineages is None:
                n_lineages = estim.eigendecomposition["eigengap"] + 1

        if n_lineages > 1:
            estim.compute_schur(n_lineages, method=kwargs.pop("method", "krylov"))

        try:
            estim.compute_macrostates(
                n_states=n_lineages, cluster_key=cluster_key, **kwargs
            )
        except ValueError:
            logg.warning(
                f"Computing `{n_lineages}` macrostates cuts through a block of complex conjugates. "
                f"Increasing `n_lineages` to {n_lineages + 1}"
            )
            estim.compute_macrostates(
                n_states=n_lineages + 1, cluster_key=cluster_key, **kwargs
            )

        fs_kwargs = {"n_cells": kwargs["n_cells"]} if "n_cells" in kwargs else {}

        if n_lineages is None:
            estim.compute_terminal_states(method="eigengap", **fs_kwargs)
        else:
            estim.set_terminal_states_from_macrostates(**fs_kwargs)
    else:
        raise NotImplementedError(type(estim))

    if compute_absorption_probabilities:
        estim.compute_absorption_probabilities(keys=keys)

    return estim


def _initial_terminal(
    adata: AnnData,
    estimator: type(BaseEstimator) = GPCCA,
    backward: bool = False,
    mode: str = VelocityMode.DETERMINISTIC,
    backward_mode: str = BackwardMode.TRANSPOSE,
    n_states: Optional[int] = None,
    cluster_key: Optional[str] = None,
    key: Optional[str] = None,
    force_recompute: bool = False,
    show_plots: bool = False,
    copy: bool = False,
    return_estimator: bool = False,
    fit_kwargs: Mapping = MappingProxyType({}),
    **kwargs,
) -> Optional[Union[AnnData, BaseEstimator]]:
    _check_estimator_type(estimator)

    if copy:
        adata = adata.copy()
    try:
        if force_recompute:
            raise KeyError("Forcing transition matrix recomputation.")
        kernel = PrecomputedKernel(key, adata=adata, backward=backward)
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

    mc = estimator(kernel)

    if cluster_key is None:
        _info_if_obs_keys_categorical_present(
            adata,
            keys=["leiden", "louvain", "cluster", "clusters"],
            msg_fmt="Found categorical observation in `adata.obs[{!r}]`. Consider specifying it as `cluster_key`.",
        )

    _fit(
        mc,
        n_lineages=n_states,
        cluster_key=cluster_key,
        compute_absorption_probabilities=False,
        **fit_kwargs,
    )

    if show_plots:
        mc.plot_spectrum(real_only=True)
        if isinstance(mc, CFLARE):
            mc.plot_terminal_states(discrete=True, same_plot=False)
        elif isinstance(mc, GPCCA):
            n_states = len(mc.macrostates.cat.categories)
            mc.plot_terminal_states(discrete=True, same_plot=False)
            if n_states > 1:
                mc.plot_coarse_T()
        else:
            raise NotImplementedError(
                f"Pipeline not implemented for `{type(mc).__name__!r}.`"
            )

    return mc.adata if copy else mc if return_estimator else None


@_deprecate(version="2.0")
@inject_docs(m=VelocityMode, b=BackwardMode)
@d.dedent
@inject_docs(
    __doc__=_docstring.format(
        direction=_initial,
        key_added="initial_states",
        bwd_mode="\n%(velocity_backward_mode_high_lvl)s",
    )
)
def initial_states(  # noqa: D103
    adata: AnnData,
    estimator: type(BaseEstimator) = GPCCA,
    mode: str = VelocityMode.DETERMINISTIC,
    backward_mode: str = BackwardMode.TRANSPOSE,
    n_states: Optional[int] = None,
    cluster_key: Optional[str] = None,
    key: Optional[str] = None,
    force_recompute: bool = False,
    show_plots: bool = False,
    copy: bool = False,
    return_estimator: bool = False,
    fit_kwargs: Mapping = MappingProxyType({}),
    **kwargs,
) -> Optional[Union[AnnData, BaseEstimator]]:

    return _initial_terminal(
        adata,
        estimator=estimator,
        mode=mode,
        backward_mode=backward_mode,
        backward=True,
        n_states=n_states,
        cluster_key=cluster_key,
        key=key,
        force_recompute=force_recompute,
        show_plots=show_plots,
        copy=copy,
        return_estimator=return_estimator,
        fit_kwargs=fit_kwargs,
        **kwargs,
    )


@_deprecate(version="2.0")
@inject_docs(m=VelocityMode, b=BackwardMode)
@d.dedent
@inject_docs(
    __doc__=_docstring.format(
        direction=_terminal, key_added="terminal_states", bwd_mode=""
    )
)
def terminal_states(  # noqa: D103
    adata: AnnData,
    estimator: type(BaseEstimator) = GPCCA,
    mode: str = VelocityMode.DETERMINISTIC,
    n_states: Optional[int] = None,
    cluster_key: Optional[str] = None,
    key: Optional[str] = None,
    force_recompute: bool = False,
    show_plots: bool = False,
    copy: bool = False,
    return_estimator: bool = False,
    fit_kwargs: Mapping = MappingProxyType({}),
    **kwargs,
) -> Optional[Union[AnnData, BaseEstimator]]:

    return _initial_terminal(
        adata,
        estimator=estimator,
        mode=mode,
        backward=False,
        n_states=n_states,
        cluster_key=cluster_key,
        key=key,
        force_recompute=force_recompute,
        show_plots=show_plots,
        copy=copy,
        return_estimator=return_estimator,
        fit_kwargs=fit_kwargs,
        **kwargs,
    )
