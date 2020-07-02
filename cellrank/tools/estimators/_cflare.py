# -*- coding: utf-8 -*-
"""Clustering Left and Right Eigenvectors (CFLARE) module."""
from copy import copy, deepcopy
from typing import Any, Dict, List, Tuple, Union, Iterable, Optional, Sequence

import numpy as np
import scipy
from pandas import Series
from scipy.stats import zscore, entropy
from scipy.linalg import solve as solve
from scipy.sparse import issparse
from scipy.sparse.linalg import gmres

import matplotlib as mpl
import matplotlib.cm as cm

import scvelo as scv
from scanpy import logging as logg
from anndata import AnnData

from cellrank.tools._utils import (
    _cluster_X,
    _normalize,
    _filter_cells,
    _process_series,
    _complex_warning,
    _get_connectivities,
)
from cellrank.tools._colors import _convert_to_hex_colors, _create_categorical_colors
from cellrank.tools._lineage import Lineage
from cellrank.tools._constants import _dp, _probs, _colors, _lin_names
from cellrank.tools.kernels._kernel import KernelExpression
from cellrank.tools.estimators._base_estimator import BaseEstimator

EPS = np.finfo(np.float64).eps


class CFLARE(BaseEstimator):
    """
    Clustering and Filtering of Left and Right Eigenvectors based on Markov chains.

    This is one of the two main classes of CellRank. We model cellular development as a Markov chain (MC), where each
    measured cell is represented by a state in the MC. We assume that transition probabilities between these states
    have already been computed using either the :class:`cellrank.tl.kernels.Kernel` class directly or the
    :func:`cellrank.tl.transition_matrix` high level function.

    The MC is time-homogeneous, i.e. the transition probabilities don't change over time. Further, it's
    discrete, as every state in the MC is given by a measured cell state. The state space is finite, as is the number
    of measured cells and we consider discrete time-increments.

    Params
    ------
    kernel
        Kernel object that stores a transition matrix.
    adata : :class:`anndata.AnnData`
        Optional annotated data object. If given, pre-computed lineages can be read in from this.
        Otherwise, read the object from the specified :paramref:`kernel`.
    inplace
        Whether to modify :paramref:`adata` object inplace or make a copy.
    read_from_adata
        Whether to read available attributes in :paramref:`adata`, if present.
    g2m_key
        Key from :paramref:`adata` `.obs`. Can be used to detect cell-cycle driven start- or endpoints.
    s_key
        Key from :paramref:`adata` `.obs`. Can be used to detect cell-cycle driven start- or endpoints.
    key_added
        Key in :paramref:`adata` where to store the final transition matrix.
    """

    def __init__(
        self,
        kernel: KernelExpression,
        adata: Optional[AnnData] = None,
        inplace: bool = True,
        read_from_adata: bool = True,
        g2m_key: Optional[str] = "G2M_score",
        s_key: Optional[str] = "S_score",
        key_added: Optional[str] = None,
    ):

        self._meta_states, self._meta_states_colors, self._meta_states_probs = (
            None,
            None,
            None,
        )

        super().__init__(
            kernel,
            adata,
            inplace=inplace,
            read_from_adata=read_from_adata,
            g2m_key=g2m_key,
            s_key=s_key,
            key_added=key_added,
        )

    def _read_from_adata(
        self, g2m_key: Optional[str] = None, s_key: Optional[str] = None, **kwargs
    ) -> None:
        if f"eig_{self._direction}" in self._adata.uns.keys():
            self._eig = self._adata.uns[f"eig_{self._direction}"]
        else:
            logg.debug(
                f"DEBUG: `eig_{self._direction}` not found. Setting `.eig` to `None`"
            )

        if self._rc_key in self._adata.obs.keys():
            self._meta_states = self._adata.obs[self._rc_key]
        else:
            logg.debug(
                f"DEBUG: `{self._rc_key}` not found in `adata.obs`. Setting `.metastable_states` to `None`"
            )

        if _colors(self._rc_key) in self._adata.uns.keys():
            self._meta_states_colors = self._adata.uns[_colors(self._rc_key)]
        else:
            logg.debug(
                f"DEBUG: `{_colors(self._rc_key)}` not found in `adata.uns`. "
                f"Setting `.metastable_states_colors`to `None`"
            )

        if self._lin_key in self._adata.obsm.keys():
            lineages = range(self._adata.obsm[self._lin_key].shape[1])
            colors = _create_categorical_colors(len(lineages))
            self._lin_probs = Lineage(
                self._adata.obsm[self._lin_key],
                names=[f"Lineage {i + 1}" for i in lineages],
                colors=colors,
            )
            self._adata.obsm[self._lin_key] = self._lin_probs
        else:
            logg.debug(
                f"DEBUG: `{self._lin_key}` not found in `adata.obsm`. Setting `.lin_probs` to `None`"
            )

        if _dp(self._lin_key) in self._adata.obs.keys():
            self._dp = self._adata.obs[_dp(self._lin_key)]
        else:
            logg.debug(
                f"DEBUG: `{_dp(self._lin_key)}` not found in `adata.obs`. Setting `.diff_potential` to `None`"
            )

        if g2m_key and g2m_key in self._adata.obs.keys():
            self._G2M_score = self._adata.obs[g2m_key]
        else:
            logg.debug(
                f"DEBUG: `{g2m_key}` not found in `adata.obs`. Setting `.G2M_score` to `None`"
            )

        if s_key and s_key in self._adata.obs.keys():
            self._S_score = self._adata.obs[s_key]
        else:
            logg.debug(
                f"DEBUG: `{s_key}` not found in `adata.obs`. Setting `.S_score` to `None`"
            )

        if _probs(self._rc_key) in self._adata.obs.keys():
            self._meta_states_probs = self._adata.obs[_probs(self._rc_key)]
        else:
            logg.debug(
                f"DEBUG: `{_probs(self._rc_key)}` not found in `adata.obs`. "
                f"Setting `.metastable_states_probs` to `None`"
            )

        if self._lin_probs is not None:
            if _lin_names(self._lin_key) in self._adata.uns.keys():
                self._lin_probs = Lineage(
                    np.array(self._lin_probs),
                    names=self._adata.uns[_lin_names(self._lin_key)],
                    colors=self._lin_probs.colors,
                )
                self._adata.obsm[self._lin_key] = self._lin_probs
            else:
                logg.debug(
                    f"DEBUG: `{_lin_names(self._lin_key)}` not found in `adata.uns`. "
                    f"Using default names"
                )

            if _colors(self._lin_key) in self._adata.uns.keys():
                self._lin_probs = Lineage(
                    np.array(self._lin_probs),
                    names=self._lin_probs.names,
                    colors=self._adata.uns[_colors(self._lin_key)],
                )
                self._adata.obsm[self._lin_key] = self._lin_probs
            else:
                logg.debug(
                    f"DEBUG: `{_colors(self._lin_key)}` not found in `adata.uns`. "
                    f"Using default colors"
                )

    def compute_eig(self, k: int = 20, which: str = "LR", alpha: float = 1) -> None:
        """
        Compute eigendecomposition of the transition matrix.

        Uses a sparse implementation, if possible, and only computes the top k eigenvectors
        to speed up the computation. Computes both left and right eigenvectors.

        Params
        ------
        k
            Number of eigenvalues/vectors to compute.
        which
            Eigenvalues are in general complex. `'LR'` - largest real part, `'LM'` - largest magnitude.
        alpha
            Used to compute the `eigengap`. :paramref:`alpha` is the weight given
            to the deviation of an eigenvalue from one.

        Returns
        -------
        None
            Nothing, but updates the following fields: :paramref:`eigendecomposition`.
        """

        self._compute_eig(k, which=which, alpha=alpha, only_evals=False)

    def plot_eig_embedding(
        self,
        left: bool = True,
        use: Optional[Union[int, Tuple[int], List[int]]] = None,
        abs_value: bool = False,
        use_imag: bool = False,
        cluster_key: Optional[str] = None,
        **kwargs,
    ) -> None:
        """
        Plot eigenvectors in an embedding.

        Params
        ------
        left
            Whether to use left or right eigenvectors.
        use
            Which or how many eigenvectors to be plotted. If `None`, it will be chosen by `eigengap`.
        abs_value
            Whether to take the absolute value before plotting.
        use_imag
            Whether to show real or imaginary part for complex eigenvectors
        cluster_key
            Key from :paramref:`adata` `.obs` to plot cluster annotations.
        kwargs
            Keyword arguments for :func:`scvelo.pl.scatter`.

        Returns
        -------
        None
            Nothing, just plots the eigenvectors.
        """

        if self._eig is None:
            raise RuntimeError("Compute eigendecomposition first as `.compute_eig()`")

        # set the direction and get the vectors
        side = "left" if left else "right"
        D, V = self._eig["D"], self._eig[f"V_{side[0]}"]

        # if irreducible, first rigth e-vec should be const.
        if side == "right":
            # quick check for irreducibility:
            if np.sum(np.isclose(D, 1, rtol=1e2 * EPS, atol=1e2 * EPS)) == 1:
                V[:, 0] = 1.0

        if use is None:
            use = self._eig["eigengap"] + 1  # add one because first e-vec has index 0

        self._plot_vectors(
            V,
            "eigen",
            abs_value=abs_value,
            cluster_key=cluster_key,
            use=use,
            use_imag=use_imag,
            D=D,
            **kwargs,
        )

    def set_metastable_states(
        self,
        labels: Union[Series, Dict[Any, Any]],
        cluster_key: Optional[str] = None,
        en_cutoff: Optional[float] = None,
        p_thresh: Optional[float] = None,
        add_to_existing: bool = False,
    ):
        """
        Set the approximate recurrent classes, if they are known a priori.

        Params
        ------
        categories
            Either a categorical :class:`pandas.Series` with index as cell names, where `NaN` marks marks a cell
            belonging to a transient state or a :class:`dict`, where each key is the name of the recurrent class and
            values are list of cell names.
        cluster_key
            If a key to cluster labels is given, :paramref"`metastable_states` will ge associated
            with these for naming and colors.
        en_cutoff
            If :paramref:`cluster_key` is given, this parameter determines when an approximate recurrent class will
            be labelled as *'Unknown'*, based on the entropy of the distribution of cells over transcriptomic clusters.
        p_thresh
            If cell cycle scores were provided, a *Wilcoxon rank-sum test* is conducted to identify cell-cycle driven
            start- or endpoints.
            If the test returns a positive statistic and a p-value smaller than :paramref:`p_thresh`,
            a warning will be issued.
        add_to_existing
            Whether to add thses categories to existing ones. Cells already belonging to recurrent classes will be
            updated if there's an overlap.
            Throws an error if previous approximate recurrent classes have not been calculated.

        Returns
        -------
        None
            Nothing, but updates the following fields: :paramref:`approx_recurrent_classes`.
        """

        self._set_categorical_labels(
            attr_key="_meta_states",
            pretty_attr_key="approx_recurrent_classes",
            cat_key=self._rc_key,
            add_to_existing_error_msg="Compute approximate recurrent classes first as `.compute_metastable_states()`.",
            categories=labels,
            cluster_key=cluster_key,
            en_cutoff=en_cutoff,
            p_thresh=p_thresh,
            add_to_existing=add_to_existing,
        )

    def compute_metastable_states(
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

        Params
        ------
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
            Nothing, but updates the following fields: :paramref:`approx_recurrent_classes`.
        """

        if self._eig is None:
            raise RuntimeError("Compute eigendecomposition first as `.compute_eig()`")

        start = logg.info("Computing approximate recurrent classes")

        if method not in ["kmeans", "louvain"]:
            raise ValueError(
                f"Invalid method `{method!r}`. Valid options are `'kmeans', 'louvain'`."
            )

        if use is None:
            use = self._eig["eigengap"] + 1  # add one b/c indexing starts at 0
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
        if muse >= self._eig["V_l"].shape[1] or muse >= self._eig["V_r"].shape[1]:
            raise ValueError(
                f"Maximum specified eigenvector ({muse}) is larger "
                f'than the number of computed eigenvectors ({self._eig["V_l"].shape[1]}). '
                f"Use `.compute_eig(k={muse})` to recompute the eigendecomposition."
            )

        logg.debug("DEBUG: Retrieving eigendecomposition")
        # we check for complex values only in the left, that's okay because the complex pattern
        # will be identical for left and right
        V_l, V_r = self._eig["V_l"][:, use], self._eig["V_r"].real[:, use]
        V_l = _complex_warning(V_l, use, use_imag=False)

        # compute a rc probability
        logg.debug("DEBUG: Computing probabilities of approximate recurrent classes")
        probs = self._compute_metastable_states_prob(use)
        self._meta_states_probs = probs
        self._adata.obs[_probs(self._rc_key)] = probs

        # retrieve embedding and concatenate
        if basis is not None:
            if f"X_{basis}" not in self._adata.obsm.keys():
                raise KeyError(f"Compute basis `{basis!r}` first.")
            X_em = self._adata.obsm[f"X_{basis}"][:, :n_comps]
            X = np.concatenate([V_r, X_em], axis=1)
        else:
            logg.debug("DEBUG: Basis is `None`. Setting X equal to right eigenvectors")
            X = V_r

        # filter out cells which are in the lowest q percentile in abs value in each eigenvector
        if percentile is not None:
            logg.debug("DEBUG: Filtering out cells according to percentile")
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
        if method == "kmeans":
            if n_clusters_kmeans is None:
                if percentile is not None:
                    n_clusters_kmeans = len(use)
                else:
                    n_clusters_kmeans = len(use) + 1

            if X.shape[0] < n_clusters_kmeans:
                raise ValueError(
                    f"Filtering resulted in only {X.shape[0]} cell(s), insufficient to cluster into "
                    f"{n_clusters_kmeans} clusters. Consider decreasing the value of `percentile`. "
                )

        logg.debug(
            f"DEBUG: Using `{use}` eigenvectors, basis `{basis!r}` and method `{method!r}` for clustering"
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
            rc_labels = np.repeat(None, self._adata.n_obs)
            rc_labels[ixs] = labels
        else:
            rc_labels = labels
        rc_labels = Series(rc_labels, index=self._adata.obs_names, dtype="category")
        rc_labels.cat.categories = list(rc_labels.cat.categories.astype("str"))

        # filtering to get rid of some of the left over transient states
        if n_matches_min > 0:
            logg.debug("DEBUG: Filtering according to `n_matches_min`")
            distances = _get_connectivities(
                self._adata, mode="distances", n_neighbors=n_neighbors_filtering
            )
            rc_labels = _filter_cells(
                distances, rc_labels=rc_labels, n_matches_min=n_matches_min
            )

        self.set_metastable_states(
            labels=rc_labels,
            cluster_key=cluster_key,
            en_cutoff=en_cutoff,
            p_thresh=p_thresh,
            add_to_existing=False,
        )

        logg.info(
            f"Adding `adata.obs[{_probs(self._rc_key)!r}]`\n"
            f"       `adata.obs[{self._rc_key!r}]`\n"
            f"       `.approx_recurrent_classes_probabilities`\n"
            f"       `.approx_recurrent_classes`\n"
            f"    Finish",
            time=start,
        )

    def plot_metastable_states(
        self, cluster_key: Optional[str] = None, **kwargs
    ) -> None:
        """
        Plot the approximate recurrent classes in a given embedding.

        Params
        ------
        cluster_key
            Key from `.obs` to plot clusters.
        kwargs
            Keyword arguments for :func:`scvelo.pl.scatter`.

        Returns
        -------
        None
            Nothing, just plots the approximate recurrent classes.
        """

        if self._meta_states is None:
            raise RuntimeError(
                "Compute approximate recurrent classes first as `.compute_metastable_states()`"
            )

        self._adata.obs[self._rc_key] = self._meta_states

        # check whether the length of the color array matches the number of clusters
        color_key = _colors(self._rc_key)
        if color_key in self._adata.uns and len(self._adata.uns[color_key]) != len(
            self._meta_states.cat.categories
        ):
            del self._adata.uns[_colors(self._rc_key)]
            self._meta_states_colors = None

        color = self._rc_key if cluster_key is None else [cluster_key, self._rc_key]
        scv.pl.scatter(self._adata, color=color, **kwargs)

        if color_key in self._adata.uns:
            self._meta_states_colors = self._adata.uns[color_key]

    def compute_lin_probs(
        self,
        keys: Optional[Sequence[str]] = None,
        check_irred: bool = False,
        norm_by_frequ: bool = False,
        use_iterative_solver: Optional[bool] = None,
        tol: float = 1e-5,
        use_initialization: bool = True,
    ) -> None:
        """
        Compute absorption probabilities for a Markov chain.

        For each cell, this computes the probability of it reaching any of the approximate recurrent classes.
        This also computes the entropy over absorption probabilities, which is a measure of cell plasticity, see
        [Setty19]_.

        Params
        ------
        keys
            Comma separated sequence of keys defining the recurrent classes.
        check_irred
            Check whether the matrix restricted to the given transient states is irreducible.
        norm_by_frequ
            Divide absorption probabilities for `rc_i` by `|rc_i|`.
        use_iterative_solver
            Whether to use an iterative solver for the linear system. Makes sense for large problems (> 20k cells)
            with few final states (<50).
        tol
            Convergence tolerance for the iterative solver. The default is fine for most cases, only consider
            decreasing this for severely ill-conditioned matrices.
        use_initialization
            Only relevant when using an iterative solver. In that case, the solution of absorbing states from the same
            recurrent class can be used as initialization to the iterative solver.

        Returns
        -------
        None
            Nothing, but updates the following fields: :paramref:`lineage_probabilities`, :paramref:`diff_potential`.
        """

        if self._meta_states is None:
            raise RuntimeError(
                "Compute approximate recurrent classes first as `.compute_metastable_states()`"
            )
        if keys is not None:
            keys = sorted(set(keys))

        # Note: There are three relevant data structures here
        # - self.metastable_states: pd.Series which contains annotations for approx rcs. Associated colors in
        #   self.metastable_states_colors
        # - self.lin_probs: Linage object which contains the lineage probabilities with associated names and colors
        # -_metastable_states: pd.Series, temporary copy of self.approx rcs used in the context of this function.
        #   In this copy, some metastable_states may be removed or combined with others
        start = logg.info("Computing absorption probabilities")

        # get the transition matrix
        t = self._T
        n_cells = t.shape[0]
        if not issparse(t):
            logg.warning(
                "Attempting to solve a potentially large linear system with dense transition matrix"
            )

        # colors are created in `compute_metastable_states`, this is just in case
        self._check_and_create_colors()

        # process the current annotations according to `keys`
        metastable_states_, colors_ = _process_series(
            series=self._meta_states, keys=keys, colors=self._meta_states_colors
        )

        #  create empty lineage object
        if self._lin_probs is not None:
            logg.debug("DEBUG: Overwriting `.lin_probs`")
        self._lin_probs = Lineage(
            np.empty((1, len(colors_))),
            names=metastable_states_.cat.categories,
            colors=colors_,
        )

        # warn in case only one state is left
        keys = list(metastable_states_.cat.categories)
        if len(keys) == 1:
            logg.warning(
                "There is only one recurrent class, all cells will have probability 1 of going there"
            )

        # create arrays of all recurrent and transient indices
        mask = np.repeat(False, len(metastable_states_))
        rec_start_indices = [0]
        for cat in metastable_states_.cat.categories:
            mask = np.logical_or(mask, metastable_states_ == cat)
            rec_start_indices.append(np.sum(mask))
        rec_indices, trans_indices = np.where(mask)[0], np.where(~mask)[0]

        # create Q (restriction transient-transient), S (restriction transient-recurrent) and I (Q-sized identity)
        q = t[trans_indices, :][:, trans_indices]
        s = t[trans_indices, :][:, rec_indices]

        if check_irred:
            if self._is_irreducible is None:
                self.compute_partition()
            if not self._is_irreducible:
                logg.warning("Restriction Q is not irreducible")

        # determine whether it makes sense you use a iterative solver
        if use_iterative_solver is None:
            if issparse(t) and n_cells > 1e4 and s.shape[1] < 100:
                use_iterative_solver = True
            elif issparse(t) and n_cells > 5 * 1e5:
                use_iterative_solver = True
            else:
                use_iterative_solver = False
        solver = "an interative" if use_iterative_solver else "a direct"
        logg.debug(
            f"DEBUG: Found {n_cells} cells and {s.shape[1]} absorbing states. Using {solver} solver"
        )

        if not use_iterative_solver:
            logg.debug("DEBUG: Densifying matrices for direct solver ")
            q = q.toarray() if issparse(q) else q
            s = s.toarray() if issparse(s) else s
            eye = np.eye(len(trans_indices))
        else:
            if not issparse(q):
                q = scipy.sparse.csr_matrix(q)
            if not issparse(s):
                s = scipy.sparse.csr_matrix(s)
            eye = scipy.sparse.eye(len(trans_indices))

        # compute abs probs. Since we don't expect sparse solution, dense computation is faster.
        if use_iterative_solver:

            def flex_solve(M, B, solver=scipy.sparse.linalg.gmres, init_indices=None):
                # solve a series of linear problems using an iterative solver
                x_list, info_list = [], []
                x = None
                for ix, b in enumerate(B.T):
                    # use the previous problem solution as initialisation
                    if init_indices is not None:
                        x0 = None if ix in init_indices else x
                    else:
                        x0 = None
                    x, info = solver(M, b.toarray().flatten(), tol=tol, x0=x0)
                    x_list.append(x[:, None])
                    info_list.append(info)

                return np.concatenate(x_list, axis=1), info_list

            logg.debug("DEBUG: Solving the linear system using GMRES")
            init_indices = rec_start_indices if use_initialization else None
            abs_states, info = flex_solve(eye - q, s, gmres, init_indices=init_indices)
            if not (all(con == 0 for con in info)):
                logg.warning("Some linear solves did not converge for GMRES")
        else:
            logg.debug("DEBUG: Solving the linear system using direct factorization")
            abs_states = solve(eye - q, s)

        # aggregate to class level by summing over columns belonging to the same metastable_states
        approx_rc_red = metastable_states_[mask]
        rec_classes_red = {
            key: np.where(approx_rc_red == key)[0]
            for key in approx_rc_red.cat.categories
        }
        _abs_classes = np.concatenate(
            [
                np.array(abs_states[:, rec_classes_red[key]].sum(1)).flatten()[:, None]
                for key in approx_rc_red.cat.categories
            ],
            axis=1,
        )

        if norm_by_frequ:
            logg.debug("DEBUG: Normalizing by frequency")
            _abs_classes /= [len(value) for value in rec_classes_red.values()]
        _abs_classes = _normalize(_abs_classes)

        # for recurrent states, set their self-absorption probability to one
        abs_classes = np.zeros((self._n_states, len(rec_classes_red)))
        rec_classes_full = {
            cl: np.where(metastable_states_ == cl)
            for cl in metastable_states_.cat.categories
        }
        for col, cl_indices in enumerate(rec_classes_full.values()):
            abs_classes[trans_indices, col] = _abs_classes[:, col]
            abs_classes[cl_indices, col] = 1

        self._dp = entropy(abs_classes.T)
        self._lin_probs = Lineage(
            abs_classes,
            names=list(self._lin_probs.names),
            colors=list(self._lin_probs.colors),
        )

        self._adata.obsm[self._lin_key] = self._lin_probs
        self._adata.obs[_dp(self._lin_key)] = self._dp
        self._adata.uns[_lin_names(self._lin_key)] = self._lin_probs.names
        self._adata.uns[_colors(self._lin_key)] = self._lin_probs.colors

        logg.info("    Finish", time=start)

    def plot_lin_probs(
        self,
        lineages: Optional[Union[str, Iterable[str]]] = None,
        cluster_key: Optional[str] = None,
        mode: str = "embedding",
        time_key: str = "latent_time",
        show_dp: bool = False,
        same_plot: bool = False,
        title: Optional[str] = None,
        cmap: Union[str, mpl.colors.ListedColormap] = cm.viridis,
        **kwargs,
    ) -> None:
        """
        Plot the absorption probabilities in the given embedding.

        Params
        ------
        lineages
            Only show these lineages. If `None`, plot all lineages.
        cluster_key
            Key from :paramref`adata: `.obs` for plotting cluster labels.
        mode
            Can be either `'embedding'` or `'time'`.

            - If `'embedding'`, plot the embedding while coloring in the absorption probabilities.
            - If `'time'`, plos the pseudotime on x-axis and the absorption probabilities on y-axis.
        time_key
            Key from `adata.obs` to use as a pseudotime ordering of the cells.
        show_dp
            Whether to show :paramref:`diff_potential` if present.
        same_plot
            Whether to plot the lineages on the same plot using color gradients when :paramref:`mode='embedding'`.
        title
            Either `None`, in which case titles are "to/from final/root state X",
            or an array of titles, one per lineage.
        cmap
            Colormap to use.
        kwargs
            Keyword arguments for :func:`scvelo.pl.scatter`.

        Returns
        -------
        None
            Nothing, just plots the absorption probabilities.
        """

        self._plot_probabilities(
            attr="_lin_probs",
            error_msg="Compute lineage probabilities first as `.compute_lin_probs()`.",
            lineages=lineages,
            cluster_key=cluster_key,
            mode=mode,
            time_key=time_key,
            show_dp=show_dp,
            title=title,
            same_plot=same_plot,
            color_map=cmap,
            **kwargs,
        )

    def _compute_metastable_states_prob(
        self, use: Union[Tuple[int], List[int], range]
    ) -> np.ndarray:
        """Compute a global score of being an approximate recurrent class."""

        if self._eig is None:
            raise RuntimeError("Compute eigendecomposition first as `.compute_eig()`.")

        # get the truncated eigendecomposition
        V, evals = self._eig["V_l"].real[:, use], self._eig["D"].real[use]

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

        return c

    def _check_and_create_colors(self):
        n_cats = len(self._meta_states.cat.categories)
        color_key = _colors(self._rc_key)

        if self._meta_states_colors is None:
            if color_key in self._adata.uns and n_cats == len(
                self._adata.uns[color_key]
            ):
                logg.debug("DEBUG: Loading colors from `.adata` object")
                self._meta_states_colors = _convert_to_hex_colors(
                    self._adata.uns[color_key]
                )
            else:
                self._meta_states_colors = _create_categorical_colors(n_cats)
                self._adata.uns[color_key] = self._meta_states_colors
        elif len(self._meta_states_colors) != n_cats:
            self._meta_states_colors = _create_categorical_colors(n_cats)
            self._adata.uns[color_key] = self._meta_states_colors

    def _get_restriction_to_main(self) -> Tuple[Series, np.ndarray]:
        """
        Restrict the categorical of metastable states.

        This restricts the categorical Series object where we store metastable states to the set of those states
        that we computed lineage probabilities for. This is a utility function - it is needed because in CFLARE,
        we currently have no possibility to conveniently restrict the metastable states to a core set of main states,
        other than by computing lineage probabilities
        #TODO this won't be able to deal with combined states

        Returns
        -------
        :class:`pandas.Series`, :class:`numpy.ndararay`
            The restricted categorical annotations and matching colors.
        """

        # get the names of the main states, remove 'rest' if present
        main_names = self.lineage_probabilities.names
        main_names = main_names[main_names != "rest"]

        # get the metastable annotations & colors
        cats_main = self.metastable_states.copy()
        colors_main = np.array(self._meta_states_colors.copy())

        # restrict both colors and categories
        mask = np.in1d(cats_main.cat.categories, main_names)
        colors_main = colors_main[mask]
        cats_main.cat.remove_categories(cats_main.cat.categories[~mask], inplace=True)

        return cats_main, colors_main

    def copy(self) -> "CFLARE":
        """
        Returns
        -------
        :class:`cellrank.tl.CFLARE`
            A copy of itself.
        """  # noqa

        kernel = copy(self.kernel)  # doesn't copy the adata object
        c = CFLARE(kernel, self.adata.copy(), inplace=False, read_from_adata=False)

        c._is_irreducible = self.irreducible
        c._rec_classes = copy(self._rec_classes)
        c._trans_classes = copy(self._trans_classes)

        c._eig = deepcopy(self.eigendecomposition)
        c._lin_probs = copy(self.lineage_probabilities)
        c._dp = copy(self.diff_potential)

        c._meta_states = copy(self._meta_states)
        c._meta_states_probs = copy(self._meta_states_probs)
        c._meta_states_colors = copy(self._meta_states_colors)

        c._G2M_score = copy(self._G2M_score)
        c._S_score = copy(self._S_score)

        c._g2m_key = self._g2m_key
        c._s_key = self._s_key
        c._key_added = self._key_added

        return c

    @property
    def metastable_states(self) -> Series:
        """
        Returns
        -------
        :class:`pandas.Series`
            The approximate recurrent classes, where `NaN` marks cells which are transient.
        """  # noqa
        return self._meta_states

    @property
    def metastable_states_probabilities(self) -> Series:
        """
        Returns
        -------
        :class:`pandas.Series`
            Probabilities of cells belonging to the approximate recurrent classes.
        """  # noqa
        return self._meta_states_probs
