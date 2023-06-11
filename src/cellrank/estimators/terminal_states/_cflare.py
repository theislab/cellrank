from typing import Any, Dict, List, Literal, Optional, Sequence, Union

import numpy as np
import pandas as pd
import scipy.stats as st

from anndata import AnnData

from cellrank import logging as logg
from cellrank._utils._docs import d
from cellrank._utils._utils import (
    _cluster_X,
    _complex_warning,
    _filter_cells,
    _get_connectivities,
)
from cellrank.estimators.mixins import EigenMixin, LinDriversMixin
from cellrank.estimators.terminal_states._term_states_estimator import (
    TermStatesEstimator,
)
from cellrank.kernels._utils import _get_basis

__all__ = ["CFLARE"]


@d.dedent
class CFLARE(TermStatesEstimator, LinDriversMixin, EigenMixin):
    """Clustering and filtering of left and right eigenvectors (CFLARE).

    This estimator computes initial and terminal states, as well as fate probabilities.
    It uses the left eigenvectors of the transition matrix to filter to a set of (approximately)
    recurrent cells and the right eigenvectors to cluster this set of cells into discrete groups.
    CFLARE exists mostly for legacy reasons, we recommend using the :class:`cellrank.estimators.GPCCA`
    estimator, which is mathematically more principled.

    Parameters
    ----------
    %(base_estimator.parameters)s
    """

    @d.dedent
    def fit(self, k: int = 20, **kwargs: Any) -> "TermStatesEstimator":
        """Prepare self for terminal states prediction.

        Parameters
        ----------
        k
            Number of eigenvectors to compute.
        kwargs
            Keyword arguments for :meth:`compute_eigendecomposition`.

        Returns
        -------
        Returns self and updates the following fields:

        - :attr:`eigendecomposition` - %(eigen.summary)s
        """
        return self.compute_eigendecomposition(k=k, only_evals=False, **kwargs)

    @d.dedent
    def predict(
        self,
        use: Optional[Union[int, Sequence[int]]] = None,
        percentile: Optional[int] = 98,
        method: Literal["leiden", "kmeans"] = "leiden",
        cluster_key: Optional[str] = None,
        n_clusters_kmeans: Optional[int] = None,
        n_neighbors: int = 20,
        resolution: float = 0.1,
        n_matches_min: int = 0,
        n_neighbors_filtering: int = 15,
        basis: Optional[str] = None,
        n_comps: int = 5,
        scale: Optional[bool] = None,
    ) -> "CFLARE":
        """Find approximate recurrent classes of the Markov chain.

        Filter to obtain recurrent states from left eigenvectors.
        Cluster to obtain approximate recurrent classes from right eigenvectors.

        Parameters
        ----------
        use
            Which or how many first eigenvectors to use as features for filtering and clustering.
            If :obj:`None`, use the *eigengap* statistic.
        percentile
            Threshold used for filtering out cells which are most likely transient states. Cells which are in the
            lower ``percentile`` percent of each eigenvector will be removed from the data matrix.
        method
            Method to be used for clustering. Valid option are:

            - ``'kmeans'`` - :class:`~sklearn.cluster.KMeans`.
            - ``'leiden'`` - :func:`~scanpy.tl.leiden`.
        cluster_key
            Key in :attr:`~anndata.AnnData.obs` in order to associate names and colors with :attr:`terminal_states`.
        n_clusters_kmeans
            Number of clusters when ``method = 'kmeans'``. If :obj:`None`, this is set to ``use + 1``.
        n_neighbors
            Number of neighbors in a kNN graph. This is the :math:`K` parameter for that,
            the number of neighbors for each cell. Only used when ``method = 'leiden'``.
        resolution
            Resolution parameter for :func:`~scanpy.tl.leiden`. Should be chosen relatively small.
        n_matches_min
            Filters out cells which don't have at least ``n_matches_min`` neighbors from the same category.
            This filters out some cells which are transient but have been misassigned.
        n_neighbors_filtering
            Parameter for filtering cells. Cells are filtered out if they don't have at least ``n_matches_min``
            neighbors among the ``n_neighbors_filtering`` nearest cells.
        basis
            Key from :attr:`~anndata.AnnData.obsm` as additional features for clustering.
            If `None`, use only the right eigenvectors.
        n_comps
            Number of embedding components to be used when ``basis != None``.
        scale
            Scale the values to z-scores. If :obj:`None`, scale the values if ``basis != None``.

        Returns
        -------
        Returns self and just updates the following fields:

        - :attr:`terminal_states` - %(tse_term_states.summary)s
        - :attr:`terminal_states_probabilities` - %(tse_term_states_probs.summary)s
        """

        def convert_use(use: Optional[Union[int, Sequence[int], np.ndarray]]) -> List[int]:
            if method not in ["kmeans", "leiden"]:
                raise ValueError(f"Invalid method `{method!r}`. Valid options are `leiden` or `kmeans`.")

            if use is None:
                use = eig["eigengap"] + 1  # add one b/c indexing starts at 0
            if isinstance(use, int):
                use = list(range(use))
            elif not isinstance(use, (np.ndarray, Sequence)):
                raise TypeError(f"Expected `use` to be `int` or a `Sequence`, found `{type(use).__name__}`.")
            use = list(use)

            if not use:
                raise ValueError("No eigenvectors have been selected.")

            muse, mevecs = max(use), max(eig["V_l"].shape[1], eig["V_r"].shape[1])
            if muse >= mevecs:
                raise ValueError(
                    f"Maximum specified eigenvector `{muse}` is larger than the number of eigenvectors `{mevecs}`. "
                    f"Use `.compute_eigendecomposition(k={muse})` to recompute the eigendecomposition."
                )

            return use

        eig = self.eigendecomposition
        if eig is None:
            raise RuntimeError("Compute eigendecomposition first as `.compute_eigendecomposition()`.")
        use = convert_use(use)

        start = logg.info("Computing terminal states")
        # we check for complex values only in the left, that's okay because the complex pattern
        # will be identical for left and right
        V_l, V_r = eig["V_l"][:, use], eig["V_r"].real[:, use]
        V_l = _complex_warning(V_l, use, use_imag=False)

        # retrieve embedding and concatenate
        if basis is not None:
            X_em = _get_basis(self.adata, basis=basis)[:, :n_comps]
            X = np.concatenate([V_r, X_em], axis=1)
        else:
            X = V_r

        # filter out cells which are in the lowest q percentile in abs value in each eigenvector
        if percentile is not None:
            logg.debug("Filtering out cells according to percentile")
            if not (0 <= percentile <= 100):
                raise ValueError(f"Expected `percentile` to be in interval `[0, 100]`, found `{percentile}`.")
            cutoffs = np.percentile(np.abs(V_l), percentile, axis=0)
            ixs = np.any(cutoffs <= np.abs(V_l), axis=1)
            X = X[ixs, :]

        if scale is None:
            scale = basis is not None
        if scale:
            X = st.zscore(X, axis=0)

        # cluster X
        if method == "kmeans" and n_clusters_kmeans is None:
            n_clusters_kmeans = len(use) + (percentile is None)
            if X.shape[0] < n_clusters_kmeans:
                raise ValueError(
                    f"Filtering resulted in only {X.shape[0]} cell(s), insufficient to cluster into "
                    f"`{n_clusters_kmeans}` clusters. Consider decreasing the value of `percentile`."
                )

        # fmt: off
        logg.debug(f"Using `{use}` eigenvectors, basis `{basis!r}` and method `{method!r}` for clustering")
        clusters = _cluster_X(
            X,
            method=method,
            n_clusters=n_clusters_kmeans,
            n_neighbors=n_neighbors,
            resolution=resolution,
        )

        # fill in the labels in case we filtered out cells before
        if percentile is not None:
            labels = np.repeat(None, len(self))
            labels[ixs] = clusters
        else:
            labels = clusters
        labels = pd.Series(labels, index=self.adata.obs_names, dtype="category")
        labels = labels.cat.rename_categories({c: str(c) for c in labels.cat.categories})

        # filtering to get rid of some of the leftover transient states
        if n_matches_min > 0:
            logg.debug(f"Filtering according to `n_matches_min={n_matches_min}`")
            distances = _get_connectivities(self.adata, mode="distances", n_neighbors=n_neighbors_filtering)
            labels = _filter_cells(distances, rc_labels=labels, n_matches_min=n_matches_min)
        # fmt: on

        return self.set_terminal_states(
            states=labels,
            cluster_key=cluster_key,
            probs=self._compute_term_states_probs(eig, use),
            params=self._create_params(),
            time=start,
        )

    def _compute_term_states_probs(self, eig: Dict[str, Any], use: List[int]) -> pd.Series:
        # get the truncated eigendecomposition
        V, evals = eig["V_l"].real[:, use], eig["D"].real[use]

        # shift and scale
        V_pos = np.abs(V)
        V_shifted = V_pos - np.min(V_pos, axis=0)
        V_scaled = V_shifted / np.max(V_shifted, axis=0)

        # check the ranges are correct
        np.testing.assert_allclose(np.min(V_scaled, axis=0), 0)
        np.testing.assert_allclose(np.max(V_scaled, axis=0), 1)

        # further scale by the eigenvalues
        V_eigs = V_scaled / evals

        # sum over cols and scale
        c = np.sum(V_eigs, axis=1)
        c /= np.max(c)

        return pd.Series(c, index=self.adata.obs_names)

    def _read_from_adata(self, adata: AnnData, **kwargs: Any) -> bool:
        ok = super()._read_from_adata(adata, **kwargs)
        ok = ok and self._read_eigendecomposition(adata, allow_missing=False)
        return ok and self._read_fate_probabilities(adata) and self._read_absorption_times(adata)
