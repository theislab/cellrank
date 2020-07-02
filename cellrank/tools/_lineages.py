# -*- coding: utf-8 -*-
"""Lineages module."""

from typing import Optional, Sequence

from scanpy import logging as logg
from anndata import AnnData

from cellrank.tools._utils import _info_if_obs_keys_categorical_present
from cellrank.utils._utils import _read_graph_data
from cellrank.tools.kernels import VelocityKernel
from cellrank.tools._constants import LinKey, StateKey, Direction, _transition
from cellrank.tools.estimators import GPCCA, CFLARE
from cellrank.tools.estimators._base_estimator import BaseEstimator


def lineages(
    adata: AnnData,
    estimator: type(BaseEstimator) = GPCCA,
    final: bool = True,
    cluster_key: Optional[str] = None,
    keys: Optional[Sequence[str]] = None,
    n_lineages: Optional[int] = None,
    method: str = "krylov",
    copy: bool = False,
    return_estimator: bool = False,
    **kwargs,
) -> Optional[AnnData]:
    """
    Compute probabilistic lineage assignment using RNA velocity.

    For each cell `i` in {1, ..., N} and root/final state j in {1, ..., M}, the probability is computed that cell `i`
    is either going to final state `j` (`final=True`) or coming from root state `j` (`final=False`). We provide two
    estimators for computing these probabilities:

    For the estimator :classL`cellrank.tl.GPCCA`, we perform Generalized Perron Cluster Cluster Analysis [GPCCA18]_.
    Cells are mapped to a simplex where each corner represents a final/root state, and the position of a cell in the
    simplex determines its probability of going to a final states/coming from a root state.

    For the estimator :class:`cellrank.tl.CFLARE`, we compute absorption probabilities towards the root/final states
    of the Markov chain.
    For related approaches in the single cell context that utilise absorption probabilities to map cells to lineages,
    see [Setty19]_ or [Weinreb18]_.

    Before running this function, compute root/final states using :func:`cellrank.tl.root_states` or
    :func:`cellrank.tl.final_states`, respectively.

    Parameters
    ----------
    adata : :class:`anndata.AnnData`
        Annotated data object
    estimator
        Estimator to use to compute the lineage probabilities.
    final
        If `True`, computes final states. Otherwise, computes root states.
    cluster_key
        Match computed {direction} states against pre-computed clusters to annotate the {direction} states.
        For this, provide a key from :paramref:`adata` `.obs` where cluster labels have been computed.
    keys
        Determines which root/final states to use by passing their names. Further, root/final states can be combined.
        If e.g. the final states are ['Neuronal_1', 'Neuronal_1', 'Astrocytes', 'OPC'], then passing
        keys=['Neuronal_1, Neuronal_2', 'OPC'] means that the two neuronal final states are treated as one and the
        'Astrocyte' state is excluded.
    n_lineages
        Number of lineages when :paramref:`estimator` `=GPCCA`. If `None`, it will be based on `eigengap`.
    method
        Method to use when computing the Schur decomposition. Only needed when :paramref:`estimator`
        is :class`:cellrank.tl.GPCCA:.
        Valid options are: `'krylov'`, `'brandts'`.
    copy
        Whether to update the existing AnnData object or to return a copy.
    return_estimator
        Whether to return the estimator. Only available when :paramref:`copy=False`.
    kwargs
        Keyword arguments for :meth:`cellrank.tl.estimators.BaseEstimator.compute_metastable_states`.

    Returns
    -------
    :class:`anndata.AnnData`, :class:`cellrank.tools.estimators.BaseEstimator` or :class:`NoneType`
        Depending on :paramref:`copy`, either updates the existing :paramref:`adata` object
        or returns a copy or returns the estimator.
    """

    if not isinstance(estimator, type):
        raise TypeError(
            f"Expected estimator to be a class, found `{type(estimator).__name__}`."
        )

    if not issubclass(estimator, BaseEstimator):
        raise TypeError(
            f"Expected estimator to be a subclass of `BaseEstimator`, found `{type(estimator).__name__}`"
        )

    # Set the keys and print info
    adata = adata.copy() if copy else adata

    if final:
        direction = Direction.FORWARD
        lin_key = LinKey.FORWARD
        rc_key = StateKey.FORWARD
    else:
        direction = Direction.BACKWARD
        lin_key = LinKey.BACKWARD
        rc_key = StateKey.BACKWARD

    transition_key = _transition(direction)
    vk = VelocityKernel(adata, backward=not final)

    try:
        vk._transition_matrix = _read_graph_data(adata, transition_key)
    except KeyError as e:
        key = "final" if final else "root"
        raise KeyError(
            f"Compute {key} states first as `cellrank.tl.find_{key}`."
        ) from e

    start = logg.info(f"Computing lineage probabilities towards `{rc_key}`")
    mc = estimator(vk, read_from_adata=False)

    if cluster_key is None:
        _info_if_obs_keys_categorical_present(
            adata,
            keys=["louvain", "clusters"],
            msg_fmt="Found categorical observation in `adata.obs[{!r}]`. "
            "Consider specifying it as `cluster_key`.",
        )

    # compute the absorption probabilities
    if isinstance(mc, CFLARE):
        mc.compute_eig()
        mc.compute_metastable_states(cluster_key=cluster_key, **kwargs)
        mc.compute_lin_probs(keys=keys)
    elif isinstance(mc, GPCCA):
        if n_lineages is None or n_lineages == 1:
            mc.compute_eig()
            if n_lineages is None:
                n_lineages = mc.eigendecomposition["eigengap"] + 1

        if n_lineages > 1:
            mc.compute_schur(n_lineages + 1, method=method)

        try:
            mc.compute_metastable_states(
                n_states=n_lineages, cluster_key=cluster_key, **kwargs
            )
        except ValueError:
            logg.warning(
                f"Computing {n_lineages} metastable states cuts through a block of complex conjugates. "
                f"Increasing `n_lineages` to {n_lineages + 1}"
            )
            mc.compute_metastable_states(
                n_states=n_lineages + 1, cluster_key=cluster_key, **kwargs
            )
        mc.set_main_states(names=keys)
    else:
        raise NotImplementedError(
            f"Pipeline not implemented for `{type(bytes).__name__}`"
        )

    logg.info(f"Added key `{lin_key!r}` to `adata.obsm`\n    Finish", time=start)

    return adata if copy else mc if return_estimator else None
