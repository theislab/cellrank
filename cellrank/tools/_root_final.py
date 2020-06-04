# -*- coding: utf-8 -*-
"""Module used for finding root and final states."""

from typing import Union, Optional

from scanpy import logging as logg
from anndata import AnnData

from cellrank.utils._docs import inject_docs
from cellrank.tools._constants import StateKey
from cellrank.tools.estimators import GPCCA, CFLARE
from cellrank.tools._transition_matrix import transition_matrix
from cellrank.tools.estimators._base_estimator import BaseEstimator

_find_docs = """\
Compute {cells} cells based on RNA velocity, see [Manno18]_.The tool models dynamic cellular
processes as a Markov chain, where the transition matrix is computed based on the velocity vectors of each
individual cell. The spectrum of the transition matrix can be used to query approximate recurrent classes of the
Markov chain, which represent groups of {cells} cells.

Cells are filtered into transient/recurrent cells using the left eigenvectors of the transition matrix and clustered
into distinct groups of {cells} cells using the right eigenvectors of the transition matrix of the Markov chain.

Params
------
adata : :class:`adata.AnnData`
    Annotated data object.
estimator
    Estimator to use to compute the lineage probabilities.
cluster_key
    The tool can match computed {direction}points against pre-computed clusters to annotate the {direction}points.
    For this, provide a key from :paramref:`adata` `.obs` where cluster labels have been computed.
weight_connectivities
    Weight given to a transition matrix computed on the basis of the KNN connectivities. Should be in `[0, 1]`. This
    can help in situations where we have noisy velocities and want to give some weight to transcriptomic similarity.
percentile
    When making a distinction between transient and recurrent cells, a percentile is used for filtering. Choose
    this value according to the percentage of transient cells you expect to see in your data.
    E.g. :paramref:`percentile` `=98` means you are expecting 98% of your cells to be transient
    and 2% to be recurrent {direction}points.
n_matches_min
    Parameter used to remove some noise. If `n_matches_min = L`, required that at least L of the nearest neighbors of
    cells *i* belong to the same {direction}point, otherwise, *i* is not considered a {direction}point itself.
n_start_end
    If you know how many {direction}points you are expecting, you can provide this number.
    Otherwise, an eigen-gap heuristic is used.
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
    cluster_key: Optional[str] = None,
    weight_connectivities: float = None,
    percentile: int = 98,
    n_matches_min: Optional[int] = 1,
    n_start_end: Optional[int] = None,
    show_plots: bool = False,
    copy: bool = False,
    return_estimator: bool = False,
    **kwargs,
) -> Optional[Union[AnnData, BaseEstimator]]:

    key = StateKey.FORWARD if final else StateKey.BACKWARD
    logg.info(f"Computing `{key}`")
    adata = adata.copy() if copy else adata

    # compute kernel object
    kernel = transition_matrix(
        adata, backward=not final, weight_connectivities=weight_connectivities
    )

    # create MarkovChain object
    mc = estimator(kernel, read_from_adata=False)

    # run the computation
    mc.compute_eig()

    if isinstance(mc, CFLARE):
        mc.compute_metastable_states(
            percentile=percentile,
            n_matches_min=n_matches_min,
            use=n_start_end,
            n_clusters_kmeans=n_start_end,
            cluster_key=cluster_key,
        )

        if show_plots:
            mc.plot_spectrum(real_only=True)
            mc.plot_eig_embedding(abs_value=True, perc=[0, 98], use=n_start_end)
            mc.plot_eig_embedding(left=False, use=n_start_end)
    elif isinstance(mc, GPCCA):
        mc.compute_metastable_states(**kwargs)
        mc.set_main_states()

        if show_plots:
            pass
    else:
        raise NotImplementedError(
            f"Pipeline not implemented for `{type(bytes).__name__}`"
        )

    return adata if copy else mc if return_estimator else None


@inject_docs(
    root=_find_docs.format(cells="root", direction="start", key_added="root_cells")
)
def find_root(
    adata: AnnData,
    estimator: type(BaseEstimator) = GPCCA,
    cluster_key: Optional[str] = None,
    weight_connectivities: float = None,
    percentile: int = 98,
    n_start_end: Optional[int] = None,
    show_plots: bool = False,
    copy: bool = False,
    return_estimator: bool = False,
    **kwargs,
) -> Optional[AnnData]:
    """
    Find root cells of a dynamic process in single cells.

    {root}
    """

    return _root_final(
        adata,
        estimator=estimator,
        final=False,
        cluster_key=cluster_key,
        weight_connectivities=weight_connectivities,
        percentile=percentile,
        n_start_end=n_start_end,
        show_plots=show_plots,
        copy=copy,
        return_estimator=return_estimator,
    )


@inject_docs(
    final=_find_docs.format(cells="final", direction="end", key_added="final_cells")
)
def find_final(
    adata: AnnData,
    estimator: type(BaseEstimator) = GPCCA,
    cluster_key: Optional[str] = None,
    weight_connectivities: float = None,
    percentile: int = 98,
    n_start_end: Optional[int] = None,
    show_plots: bool = False,
    copy: bool = False,
    return_estimator: bool = False,
    **kwargs,
) -> Optional[AnnData]:
    """
    Find final cells of a dynamic process in single cells.

    {final}
    """

    return _root_final(
        adata,
        estimator=estimator,
        final=True,
        cluster_key=cluster_key,
        weight_connectivities=weight_connectivities,
        percentile=percentile,
        n_start_end=n_start_end,
        show_plots=show_plots,
        copy=copy,
        return_estimator=return_estimator,
    )
