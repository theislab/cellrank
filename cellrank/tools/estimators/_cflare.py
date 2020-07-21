# -*- coding: utf-8 -*-
"""Clustering and Filtering of Left and Right Eigenvectors based on Markov chains."""
from typing import List, Tuple, Union, Optional

import numpy as np
from pandas import Series
from cellrank import logging as logg
from scipy.stats import zscore
from cellrank.utils._docs import d, inject_docs
from cellrank.tools._utils import (
    _cluster_X,
    _filter_cells,
    _complex_warning,
    _get_connectivities,
)
from cellrank.tools.estimators._constants import A, F, P
from cellrank.tools.estimators._decomposition import Eigen
from cellrank.tools.estimators._base_estimator import BaseEstimator


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

    @inject_docs(fin_states=P.FIN, fin_states_probs=P.FIN_PROBS)
    def compute_final_states(
        self,
        use: Optional[Union[int, Tuple[int], List[int], range]] = None,
        percentile: Optional[int] = 98,
        method: str = "kmeans",
        cluster_key: Optional[str] = None,
        n_clusters_kmeans: Optional[int] = None,
        n_neighbors_louvain: int = 20,
        resolution_louvain: float = 0.1,
        n_matches_min: Optional[int] = 0,
        n_neighbors_filtering: int = 15,
        basis: Optional[str] = None,
        n_comps: int = 5,
        scale: bool = False,
        en_cutoff: Optional[float] = 0.7,
        p_thresh: float = 1e-15,
    ) -> None:
        """
        Find approximate recurrent classes in the Markov chain.

        Filter to obtain recurrent states in left eigenvectors.
        Cluster to obtain approximate recurrent classes in right eigenvectors.

        Parameters
        ----------
        use
            Which or how many first eigenvectors to use as features for clustering/filtering.
            If `None`, use `eigengap` statistic.
        percentile
            Threshold used for filtering out cells which are most likely transient states.
            Cells which are in the lower :paramref:`percentile` percent
            of each eigenvector will be removed from the data matrix.
        method
            Method to be used for clustering. Must be one of `['louvain', 'kmeans']`.
        cluster_key
            If a key to cluster labels is given, :paramref:`metastable_states` will ge associated
            with these for naming and colors.
        n_clusters_kmeans
            If `None`, this is set to :paramref:`use` `+ 1`.
        n_neighbors_louvain
            If we use `'louvain'` for clustering cells, we need to build a KNN graph.
            This is the K parameter for that, the number of neighbors for each cell.
        resolution_louvain
            Resolution parameter from the `louvain` algorithm. Should be chosen relatively small.
        n_matches_min
            Filters out cells which don't have at least n_matches_min neighbors from the same class.
            This filters out some cells which are transient but have been misassigned.
        n_neighbors_filtering
            Parameter for filtering cells. Cells are filtered out if they don't have at
            least :paramref:`n_matches_min` neighbors.
            among their n_neighbors_filtering nearest cells.
        basis
            Key from :paramref`adata` `.obsm` to be used as additional features for the clustering.
        n_comps
            Number of embedding components to be use.
        scale
            Scale to z-scores. Consider using if appending embedding to features.
        en_cutoff
            If :paramref:`cluster_key` is given, this parameter determines when an approximate recurrent class will
            be labelled as *'Unknown'*, based on the entropy of the distribution of cells over transcriptomic clusters.
        p_thresh
            If cell cycle scores were provided, a *Wilcoxon rank-sum test* is conducted to identify cell-cycle driven
            start- or endpoints.
            If the test returns a positive statistic and a p-value smaller than :paramref:`p_thresh`,
            a warning will be issued.

        Returns
        -------
        None
            Nothing, but updates the following fields:

                - :paramref:`{fin_states}`
                - :paramref:`{fin_states_probs}`
        """

        def compute_metastable_states_prob() -> Series:
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
            if method not in ["kmeans", "louvain"]:
                raise ValueError(
                    f"Invalid method `{method!r}`. Valid options are `'kmeans', 'louvain'`."
                )

            if use is None:
                use = eig["eigengap"] + 1  # add one b/c indexing starts at 0
            if isinstance(use, int):
                use = list(range(use))
            elif not isinstance(use, (tuple, list, range)):
                raise TypeError(
                    f"Argument `use` must be either `int`, `tuple`, `list` or `range`, "
                    f"found `{type(use).__name__}`."
                )
            else:
                if not all(map(lambda u: isinstance(u, int), use)):
                    raise TypeError("Not all values in `use` argument are integers.")
            use = list(use)

            if len(use) == 0:
                raise ValueError(
                    f"Number of eigenvector must be larger than `0`,  found `{len(use)}`."
                )

            muse = max(use)
            if muse >= eig["V_l"].shape[1] or muse >= eig["V_r"].shape[1]:
                raise ValueError(
                    f"Maximum specified eigenvector ({muse}) is larger "
                    f'than the number of computed eigenvectors ({eig["V_l"].shape[1]}). '
                    f"Use `.{F.COMPUTE.fmt(P.EIG)}(k={muse})` to recompute the eigendecomposition."
                )

            return use

        eig = self._get(P.EIG)
        if eig is None:
            raise RuntimeError(
                f"Compute eigendecomposition first as `.{F.COMPUTE.fmt(P.EIG)}()`"
            )
        use = check_use(use)

        start = logg.info("Computing approximate recurrent classes")
        # we check for complex values only in the left, that's okay because the complex pattern
        # will be identical for left and right
        V_l, V_r = eig["V_l"][:, use], eig["V_r"].real[:, use]
        V_l = _complex_warning(V_l, use, use_imag=False)

        # compute a rc probability
        logg.debug("Computing probabilities of approximate recurrent classes")
        self._set(A.FIN_PROBS, compute_metastable_states_prob())

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
            # TODO: @Marius - why the percentile?
            n_clusters_kmeans = len(use) + (percentile is None)
            if X.shape[0] < n_clusters_kmeans:
                raise ValueError(
                    f"Filtering resulted in only {X.shape[0]} cell(s), insufficient to cluster into "
                    f"{n_clusters_kmeans} clusters. Consider decreasing the value of `percentile`. "
                )

        logg.debug(
            f"Using `{use}` eigenvectors, basis `{basis!r}` and method `{method!r}` for clustering"
        )
        labels = _cluster_X(
            X,
            method=method,
            n_clusters=n_clusters_kmeans,
            n_neighbors=n_neighbors_louvain,
            resolution=resolution_louvain,
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
            logg.debug("Filtering according to `n_matches_min`")
            distances = _get_connectivities(
                self.adata, mode="distances", n_neighbors=n_neighbors_filtering
            )
            rc_labels = _filter_cells(
                distances, rc_labels=rc_labels, n_matches_min=n_matches_min
            )

        self.set_final_states(
            labels=rc_labels,
            cluster_key=cluster_key,
            en_cutoff=en_cutoff,
            p_thresh=p_thresh,
            add_to_existing=False,
            time=start,
        )

    # TODO: @Marius, you've written: "this won't be able to deal with combined states", is it fixed?
    def _get_restriction_to_main(self) -> Tuple[Series, np.ndarray]:
        """
        Restrict the categorical of metastable states.

        This restricts the categorical Series object where we store metastable states to the set of those states
        that we computed lineage probabilities for. This is a utility function - it is needed because in CFLARE,
        we currently have no possibility to conveniently restrict the metastable states to a core set of main states,
        other than by computing lineage probabilities

        Returns
        -------
        :class:`pandas.Series`, :class:`numpy.ndararay`
            The restricted categorical annotations and matching colors.
        """

        # get the names of the main states, remove 'rest' if present
        main_names = self._get(P.ABS_PROBS).names
        main_names = main_names[main_names != "rest"]

        # get the metastable annotations & colors
        cats_main = self._get(P.FIN).copy()
        colors_main = np.array(self._get(A.FIN_COLORS).copy())

        # restrict both colors and categories
        mask = np.in1d(cats_main.cat.categories, main_names)
        colors_main = colors_main[mask]
        cats_main.cat.remove_categories(cats_main.cat.categories[~mask], inplace=True)

        return cats_main, colors_main
