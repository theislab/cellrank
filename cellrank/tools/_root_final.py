# -*- coding: utf-8 -*-
from anndata import AnnData
from typing import Optional
from scanpy import logging as logg

from cellrank.tools._markov_chain import MarkovChain
from cellrank.tools._constants import RcKey
from cellrank.tools._transition_matrix import transition_matrix


def root_final(
    adata: AnnData,
    final: bool = True,
    cluster_key: Optional[str] = None,
    weight_connectivities: float = 0.0,
    percentile: int = 98,
    n_start_end: Optional[int] = None,
    show_plots=False,
    copy: bool = False,
) -> Optional[AnnData]:
    """
    Root cells and final cells of a dynamic process in single cells.

    Computes root and final cells based on RNA velocity, see [Manno18]_. The tool models dynamic cellular
    processes as a Markov chain, where the transition matrix is computed based on the velocity vectors of each
    individual cell. The spectrum of the transition matrix can be used to query approximate recurrent classes of the
    Markov chain, which represent groups of root- or final cells.

    Cells are filtered into transient/recurrent cells using the left eigenvectors of the transition matrix and clustered
    into distinct groups of root or final cells using the right eigenvectors of the transition matrix of the Markov
    Chain.

    Params
    ------
    adata : :class:`adata.AnnData`
        Annotated data object.
    final
        If true, computes final cells, i.e. end points. Otherwise, computes root cells, i.e. starting points.
    cluster_key
        The tool can match computed start/end-points against pre-computed clusters to annotate the start/end-points.
        For this, provide a key from :paramref:`adata` `.obs` where cluster labels have been computed.
    weight_connectivities
        Weight given to a transition matrix computed on the basis of the KNN connectivities. Should be in `[0, 1]`. This
        can help in situations where we have noisy velocities and want to give some weight to transcriptomic similarity.
    percentile
        When making a distinction between transient and recurrent cells, a percentile is used for filtering. Choose
        this value according to the percentage of transient cells you expect to see in your data.
        E.g. :paramref:`percentile` `=98` means you are expecting 98% of your cells to be transient
        and 2% to be recurrent (start/end-points).
    n_start_end
        If you know how many start/end-points you are expecting, you can provide this number.
        Otherwise, an eigen-gap heuristic is used.
    show_plots
        Whether to show plots of the spectrum and eigenvectors in the embedding.
    copy
        Whether to update the existing :paramref:`adata` object or to return a copy.

    Returns
    -------
    :class:`anndata.AnnData` or :class:`NoneType`
        Depending on :paramref:`copy`, either updates the existing :paramref:`adata` object or returns a copy.
        Root or final cells can be found in :paramref:`adata` `.obs` under either `['root_cells']` or `['final_cells']`.
    """

    key = RcKey.FORWARD if final else RcKey.BACKWARD
    start = logg.info(f"Computing `{key}`")
    adata = adata.copy() if copy else adata

    # compute kernel object
    kernel = transition_matrix(
        adata, backward=not final, weight_connectivities=weight_connectivities
    )

    # create MarkovChain object
    mc = MarkovChain(kernel)

    # run the computation
    mc.compute_eig()
    mc.compute_approx_rcs(
        percentile=percentile,
        use=n_start_end,
        n_clusters_kmeans=n_start_end,
        cluster_key=cluster_key,
    )

    if show_plots:
        mc.plot_real_spectrum()
        mc.plot_eig_embedding(abs_value=True, perc=[0, 98], use=n_start_end)
        mc.plot_eig_embedding(left=False, use=n_start_end)

    logg.info(f"Added key `{key!r}` to `adata.obs`" f"    Finish", time=start)

    return adata if copy else None
