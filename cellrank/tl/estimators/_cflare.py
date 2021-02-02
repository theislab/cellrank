"""Clustering and Filtering of Left and Right Eigenvectors based on Markov chains."""
from typing import List, Tuple, Union, Optional, Sequence

import numpy as np
from pandas import Series
from scipy.stats import zscore

from cellrank import logging as logg
from cellrank.ul._docs import d, inject_docs
from cellrank.tl._utils import (
    _cluster_X,
    _filter_cells,
    _complex_warning,
    _get_connectivities,
)
from cellrank.tl.estimators._constants import A, P
from cellrank.tl.estimators._decomposition import Eigen
from cellrank.tl.estimators._base_estimator import BaseEstimator


@d.dedent
class CFLARE(BaseEstimator, Eigen):
    """
    Clustering and Filtering of Left and Right Eigenvectors based on Markov chains.

    This is one of the two main classes of CellRank. We model cellular development as a Markov chain (MC), where each
    measured cell is represented by a state in the MC. We assume that transition probabilities between these states
    have already been computed using either the :class:`cellrank.tl.kernels.Kernel` class directly or the
    :func:`cellrank.tl.transition_matrix` high level function.

    The MC is time-homogeneous, i.e. the transition probabilities don't change over time. Further, it's
    discrete, as every state in the MC is given by a measured cell state. The state space is finite, as is the number
    of measured cells and we consider discrete time-increments.

    Parameters
    ----------
    %(base_estimator.parameters)s
    """

    @d.dedent
    @inject_docs(fs=P.TERM, fsp=P.TERM_PROBS)
    def compute_terminal_states(
        self,
        use: Optional[Union[int, Tuple[int], List[int], range]] = None,
        percentile: Optional[int] = 98,
        method: str = "kmeans",
        cluster_key: Optional[str] = None,
        n_clusters_kmeans: Optional[int] = None,
        n_neighbors: int = 20,
        resolution: float = 0.1,
        n_matches_min: Optional[int] = 0,
        n_neighbors_filtering: int = 15,
        basis: Optional[str] = None,
        n_comps: int = 5,
        scale: bool = False,
        en_cutoff: Optional[float] = 0.7,
        p_thresh: float = 1e-15,
    ) -> None:
        """
        Find approximate recurrent classes of the Markov chain.

        Filter to obtain recurrent states in left eigenvectors.
        Cluster to obtain approximate recurrent classes in right eigenvectors.

        Parameters
        ----------
        use
            Which or how many first eigenvectors to use as features for clustering/filtering.
            If `None`, use the `eigengap` statistic.
        percentile
            Threshold used for filtering out cells which are most likely transient states. Cells which are in the
            lower ``percentile`` percent of each eigenvector will be removed from the data matrix.
        method
            Method to be used for clustering. Must be one of `'louvain'`, `'leiden'` or `'kmeans'`.
        cluster_key
            If a key to cluster labels is given, :paramref:`{fs}` will get associated with these for naming and colors.
        n_clusters_kmeans
            If `None`, this is set to ``use + 1``.
        n_neighbors
            If we use `'louvain'` or `'leiden'` for clustering cells, we need to build a KNN graph.
            This is the :math:`K` parameter for that, the number of neighbors for each cell.
        resolution
            Resolution parameter for `'louvain'` or `'leiden'` clustering. Should be chosen relatively small.
        n_matches_min
            Filters out cells which don't have at least n_matches_min neighbors from the same class.
            This filters out some cells which are transient but have been misassigned.
        n_neighbors_filtering
            Parameter for filtering cells. Cells are filtered out if they don't have at least ``n_matches_min``
            neighbors among their ``n_neighbors_filtering`` nearest cells.
        basis
            Key from :paramref`adata` ``.obsm`` to be used as additional features for the clustering.
        n_comps
            Number of embedding components to be use when ``basis`` is not `None`.
        scale
            Scale to z-scores. Consider using this if appending embedding to features.
        %(en_cutoff_p_thresh)s

        Returns
        -------
        None
            Nothing, but updates the following fields:

                - :paramref:`{fsp}`
                - :paramref:`{fs}`
        """

        def _compute_macrostates_prob() -> Series:
            """Compute a global score of being an approximate recurrent class."""

            # get the truncated eigendecomposition
            V, evals = eig["V_l"].real[:, use], eig["D"].real[use]

            # shift and scale
            V_pos = np.abs(V)
            V_shifted = V_pos - np.min(V_pos, axis=0)
            V_scaled = V_shifted / np.max(V_shifted, axis=0)

            # check the ranges are correct
            assert np.allclose(np.min(V_scaled, axis=0), 0), "Lower limit it not zero."
            assert np.allclose(np.max(V_scaled, axis=0), 1), "Upper limit is not one."

            # further scale by the eigenvalues
            V_eigs = V_scaled / evals

            # sum over cols and scale
            c_ = np.sum(V_eigs, axis=1)
            c = c_ / np.max(c_)

            return Series(c, index=self.adata.obs_names)

        def check_use(use) -> List[int]:
            if method not in ["kmeans", "louvain", "leiden"]:
                raise ValueError(
                    f"Invalid method `{method!r}`. Valid options are `'louvain'`, `'leiden'` or `'kmeans'`."
                )

            if use is None:
                use = eig["eigengap"] + 1  # add one b/c indexing starts at 0
            if isinstance(use, int):
                use = list(range(use))
            elif not isinstance(use, (tuple, list, range)):
                raise TypeError(
                    f"Argument `use` must be either `int`, `tuple`, `list` or `range`, "
                    f"found `{type(use).__name__!r}`."
                )
            else:
                if not all(map(lambda u: isinstance(u, int), use)):
                    raise TypeError("Not all values in `use` argument are integers.")
            use = list(use)

            if len(use) == 0:
                raise ValueError(
                    f"Number of eigenvector must be larger than `0`, found `{len(use)}`."
                )

            muse = max(use)
            if muse >= eig["V_l"].shape[1] or muse >= eig["V_r"].shape[1]:
                raise ValueError(
                    f"Maximum specified eigenvector `{muse}` is larger "
                    f'than the number of computed eigenvectors `{eig["V_l"].shape[1]}`. '
                    f"Use `.compute_eigendecomposition(k={muse})` to recompute the eigendecomposition."
                )

            return use

        eig = self._get(P.EIG)
        if eig is None:
            raise RuntimeError(
                "Compute eigendecomposition first as `.compute_eigendecomposition()`."
            )
        use = check_use(use)

        start = logg.info("Computing approximate recurrent classes")
        # we check for complex values only in the left, that's okay because the complex pattern
        # will be identical for left and right
        V_l, V_r = eig["V_l"][:, use], eig["V_r"].real[:, use]
        V_l = _complex_warning(V_l, use, use_imag=False)

        # compute a rc probability
        logg.debug("Computing probabilities of approximate recurrent classes")
        self._set(A.TERM_PROBS, _compute_macrostates_prob())

        # retrieve embedding and concatenate
        if basis is not None:
            bkey = f"X_{basis}"
            if bkey not in self.adata.obsm.keys():
                raise KeyError(f"Basis key `{bkey!r}` not found in `adata.obsm`")

            X_em = self.adata.obsm[bkey][:, :n_comps]
            X = np.concatenate([V_r, X_em], axis=1)
        else:
            logg.debug("Basis is `None`. Setting X equal to the right eigenvectors")
            X = V_r

        # filter out cells which are in the lowest q percentile in abs value in each eigenvector
        if percentile is not None:
            logg.debug("Filtering out cells according to percentile")
            if percentile < 0 or percentile > 100:
                raise ValueError(
                    f"Percentile must be in interval `[0, 100]`, found `{percentile}`."
                )
            cutoffs = np.percentile(np.abs(V_l), percentile, axis=0)
            ixs = np.sum(np.abs(V_l) < cutoffs, axis=1) < V_l.shape[1]
            X = X[ixs, :]

        # scale
        if scale:
            X = zscore(X, axis=0)

        # cluster X
        if method == "kmeans" and n_clusters_kmeans is None:
            n_clusters_kmeans = len(use) + (percentile is None)
            if X.shape[0] < n_clusters_kmeans:
                raise ValueError(
                    f"Filtering resulted in only {X.shape[0]} cell(s), insufficient to cluster into "
                    f"`{n_clusters_kmeans}` clusters. Consider decreasing the value of `percentile`."
                )

        logg.debug(
            f"Using `{use}` eigenvectors, basis `{basis!r}` and method `{method!r}` for clustering"
        )
        labels = _cluster_X(
            X,
            method=method,
            n_clusters=n_clusters_kmeans,
            n_neighbors=n_neighbors,
            resolution=resolution,
        )

        # fill in the labels in case we filtered out cells before
        if percentile is not None:
            rc_labels = np.repeat(None, self.adata.n_obs)
            rc_labels[ixs] = labels
        else:
            rc_labels = labels

        rc_labels = Series(rc_labels, index=self.adata.obs_names, dtype="category")
        rc_labels.cat.categories = list(rc_labels.cat.categories.astype("str"))

        # filtering to get rid of some of the left over transient states
        if n_matches_min > 0:
            logg.debug(f"Filtering according to `n_matches_min={n_matches_min}`")
            distances = _get_connectivities(
                self.adata, mode="distances", n_neighbors=n_neighbors_filtering
            )
            rc_labels = _filter_cells(
                distances, rc_labels=rc_labels, n_matches_min=n_matches_min
            )

        self.set_terminal_states(
            labels=rc_labels,
            cluster_key=cluster_key,
            en_cutoff=en_cutoff,
            p_thresh=p_thresh,
            add_to_existing=False,
            time=start,
        )

    def _fit_terminal_states(
        self,
        n_lineages: Optional[int] = None,
        keys: Optional[Sequence[str]] = None,
        cluster_key: Optional[str] = None,
        compute_absorption_probabilities: bool = True,
        **kwargs,
    ) -> None:
        self.compute_eigendecomposition(k=20 if n_lineages is None else n_lineages + 1)
        if n_lineages is None:
            n_lineages = self._get(P.EIG)["eigengap"] + 1

        self.compute_terminal_states(
            use=n_lineages,
            cluster_key=cluster_key,
            n_clusters_kmeans=n_lineages,
            method=kwargs.pop("method", "kmeans"),
            **kwargs,
        )

    @d.dedent  # because of fit
    @d.dedent
    @inject_docs(fs=P.TERM, fsp=P.TERM_PROBS, ap=P.ABS_PROBS, dp=P.DIFF_POT)
    def fit(
        self,
        n_lineages: Optional[int],
        keys: Optional[Sequence[str]] = None,
        cluster_key: Optional[str] = None,
        compute_absorption_probabilities: bool = True,
        **kwargs,
    ):
        """
        Run the pipeline, computing the %(initial_or_terminal)s states and optionally the absorption probabilities.

        It is equivalent to running::

            compute_eigendecomposition(...)
            compute_terminal_states(...)
            compute_absorption_probabilities(...)

        Parameters
        ----------
        %(fit)s
        kwargs
            Keyword arguments for :meth:`compute_terminal_states`, such as ``n_cells``.

        Returns
        -------
        None
            Nothing, just makes available the following fields:

                - :paramref:`{fsp}`
                - :paramref:`{fs}`
                - :paramref:`{ap}`
                - :paramref:`{dp}`
        """
        super().fit(
            n_lineages=n_lineages,
            keys=keys,
            cluster_key=cluster_key,
            compute_absorption_probabilities=compute_absorption_probabilities,
            **kwargs,
        )
