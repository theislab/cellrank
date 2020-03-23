# -*- coding: utf-8 -*-
from anndata import AnnData
from typing import Optional
from scanpy import logging as logg

from cellrank.tools._markov_chain import MarkovChain
from cellrank.tools._constants import RcKey
from cellrank.tools._transition_matrix import transition_matrix
from cellrank.utils._docs import inject_docs


_find_docs = """\
Computes {cells} cells based on RNA velocity, see [Manno18]_. The tool models dynamic cellular
processes as a Markov chain, where the transition matrix is computed based on the velocity vectors of each
individual cell. The spectrum of the transition matrix can be used to query approximate recurrent classes of the
Markov chain, which represent groups of {cells} cells.

Cells are filtered into transient/recurrent cells using the left eigenvectors of the transition matrix and clustered
into distinct groups of {cells} cells using the right eigenvectors of the transition matrix of the Markov chain.

Params
------
adata : :class:`adata.AnnData`
    Annotated data object.
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
    cells i belong to the same start or endpoint, otherwise, i is not considered a start/endpoint itself.
n_start_end
    If you know how many {direction}points you are expecting, you can provide this number.
    Otherwise, an eigen-gap heuristic is used.
show_plots
    Whether to show plots of the spectrum and eigenvectors in the embedding.
copy
    Whether to update the existing :paramref:`adata` object or to return a copy.

Returns
-------
:class:`anndata.AnnData` or :class:`NoneType`
    Depending on :paramref:`copy`, either updates the existing :paramref:`adata` object or returns a copy.
    Marked cells can be found in :paramref:`adata` `.obs` under `{key_added!r}`.
"""


def _root_final(
    adata: AnnData,
    final: bool = True,
    cluster_key: Optional[str] = None,
    weight_connectivities: float = None,
    percentile: int = 98,
    n_matches_min: Optional[int] = 1,
    n_start_end: Optional[int] = None,
    show_plots=False,
    copy: bool = False,
) -> Optional[AnnData]:

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
        n_matches_min=n_matches_min,
        use=n_start_end,
        n_clusters_kmeans=n_start_end,
        cluster_key=cluster_key,
    )

    if show_plots:
        mc.plot_real_spectrum()
        mc.plot_eig_embedding(abs_value=True, perc=[0, 98], use=n_start_end)
        mc.plot_eig_embedding(left=False, use=n_start_end)

    logg.info(f"Added key `{key!r}` to `adata.obs`\n    Finish", time=start)

    return adata if copy else None


@inject_docs(
    root=_find_docs.format(cells="root", direction="start", key_added="root_cells")
)
def find_root(
    adata: AnnData,
    cluster_key: Optional[str] = None,
    weight_connectivities: float = None,
    percentile: int = 98,
    n_start_end: Optional[int] = None,
    show_plots=False,
    copy: bool = False,
) -> Optional[AnnData]:
    """
    Root cells of a dynamic process in single cells.

    {root}
    """

    return _root_final(
        adata,
        final=False,
        cluster_key=cluster_key,
        weight_connectivities=weight_connectivities,
        percentile=percentile,
        n_start_end=n_start_end,
        show_plots=show_plots,
        copy=copy,
    )


@inject_docs(
    final=_find_docs.format(cells="final", direction="end", key_added="final_cells")
)
def find_final(
    adata: AnnData,
    cluster_key: Optional[str] = None,
    weight_connectivities: float = None,
    percentile: int = 98,
    n_start_end: Optional[int] = None,
    show_plots=False,
    copy: bool = False,
) -> Optional[AnnData]:
    """
    Final cells of a dynamic process in single cells.

    {final}
    """

    return _root_final(
        adata,
        final=True,
        cluster_key=cluster_key,
        weight_connectivities=weight_connectivities,
        percentile=percentile,
        n_start_end=n_start_end,
        show_plots=show_plots,
        copy=copy,
    )
