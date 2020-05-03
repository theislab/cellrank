# -*- coding: utf-8 -*-
from pathlib import Path

from cellrank.tools.kernels._kernel import KernelExpression
from typing import Optional, Tuple, Sequence, List, Any, Union, Dict

import numpy as np
import scvelo as scv
import matplotlib.pyplot as plt

from anndata import AnnData
from copy import copy, deepcopy
from pandas import Series
from scanpy import logging as logg
from scipy.linalg import solve
from scipy.stats import zscore, entropy


from cellrank.tools.estimators._base_estimator import BaseEstimator
from cellrank.tools._lineage import Lineage
from cellrank.tools._constants import _probs, _colors, _lin_names
from cellrank.tools._utils import (
    _process_series,
    _complex_warning,
    _cluster_X,
    _get_connectivities,
    _normalize,
    _filter_cells,
    _make_cat,
    _convert_to_hex_colors,
    _create_categorical_colors,
    partition,
    save_fig,
)


class MarkovChain(BaseEstimator):
    """
    Class modelling cellular development as a Markov chain.

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
        self._is_irreducible = None
        self._rec_classes = None
        self._trans_classes = None

        self._approx_rcs, self._approx_rcs_colors, self._approx_rcs_probs = (
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
            self._approx_rcs = self._adata.obs[self._rc_key]
        else:
            logg.debug(
                f"DEBUG: `{self._rc_key}` not found in `adata.obs`. Setting `.approx_rcs` to `None`"
            )

        if _colors(self._rc_key) in self._adata.uns.keys():
            self._approx_rcs_colors = self._adata.uns[_colors(self._rc_key)]
        else:
            logg.debug(
                f"DEBUG: `{_colors(self._rc_key)}` not found in `adata.uns`. "
                f"Setting `.approx_rcs_colors`to `None`"
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

        if f"{self._lin_key}_dp" in self._adata.obs.keys():
            self._dp = self._adata.obs[f"{self._lin_key}_dp"]
        else:
            logg.debug(
                f"DEBUG: `{self._lin_key}_dp` not found in `adata.obs`. Setting `.dp` to `None`"
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
            self._approx_rcs_probs = self._adata.obs[_probs(self._rc_key)]
        else:
            logg.debug(
                f"DEBUG: `{_probs(self._rc_key)}` not found in `adata.obs`. "
                f"Setting `.approx_rcs_probs` to `None`"
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

    def compute_partition(self) -> None:
        """
        Computes communication classes for the Markov chain.

        Returns
        -------
        None
            Nothing, but updates the following fields:
                - :paramref:`recurrent_classes`
                - :paramref:`transient_classes`
                - :paramref:`irreducible`
        """

        start = logg.info("Computing communication classes")

        rec_classes, trans_classes = partition(self._T)

        self._is_irreducible = len(rec_classes) == 1 and len(trans_classes) == 0

        if not self._is_irreducible:
            self._trans_classes = _make_cat(
                trans_classes, self._n_states, self._adata.obs_names
            )
            self._rec_classes = _make_cat(
                rec_classes, self._n_states, self._adata.obs_names
            )
            self._adata.obs[f"{self._rc_key}_rec_classes"] = self._rec_classes
            self._adata.obs[f"{self._rc_key}_trans_classes"] = self._trans_classes
            logg.info(
                f"Found `{(len(rec_classes))}` recurrent and `{len(trans_classes)}` transient classes\n"
                f"Adding `.recurrent_classes`\n"
                f"       `.transient_classes`\n"
                f"       `.irreducible`\n"
                f"    Finish",
                time=start,
            )
        else:
            logg.warning(
                "The transition matrix is irreducible - cannot further partition it\n    Finish",
                time=start,
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

    def plot_eig(
        self,
        dpi: int = 100,
        figsize: Optional[Tuple[float, float]] = (5, 5),
        legend_loc: Optional[str] = None,
        title: Optional[str] = None,
        save: Optional[Union[str, Path]] = None,
    ) -> None:
        """
        Plot the top eigenvalues in complex plane.

        Params
        ------
        dpi
            Dots per inch.
        figsize
            Size of the figure.
        save
            Filename where to save the plots. If `None`, just shows the plot.

        Returns
        -------
        None
            Nothing, just plots the spectrum in complex plane.
        """

        if self._eig is None:
            logg.warning(
                "No eigendecomposition found, computing with default parameters"
            )
            self.compute_eig()
        D = self._eig["D"]
        params = self._eig["params"]

        # create fiture and axes
        fig, ax = plt.subplots(nrows=1, ncols=1, dpi=dpi, figsize=figsize)

        # get the original data ranges
        lam_x, lam_y = D.real, D.imag
        x_min, x_max = np.min(lam_x), np.max(lam_x)
        y_min, y_max = np.min(lam_y), np.max(lam_y)
        x_range, y_range = x_max - x_min, y_max - y_min
        final_range = np.max([x_range, y_range]) + 0.05

        # define a function to make the data limits rectangular
        adapt_range = lambda min_, max_, range_: (
            min_ + (max_ - min_) / 2 - range_ / 2,
            min_ + (max_ - min_) / 2 + range_ / 2,
        )
        x_min_, x_max_ = adapt_range(x_min, x_max, final_range)
        y_min_, y_max_ = adapt_range(y_min, y_max, final_range)

        # plot the data and the unit circle
        ax.scatter(D.real, D.imag, marker=".", label="Eigenvalue")
        t = np.linspace(0, 2 * np.pi, 500)
        x_circle, y_circle = np.sin(t), np.cos(t)
        ax.plot(x_circle, y_circle, "k-", label="Unit circle")

        # set labels, ranges and legend
        ax.set_xlabel("Im($\lambda$)")
        ax.set_ylabel("Re($\lambda$)")
        ax.set_xlim(x_min_, x_max_)
        ax.set_ylim(y_min_, y_max_)
        key = "real part" if params["which"] == "LR" else "magnitude"
        if title is None:
            title_str = f"top {params['k']} eigenvalues according to their {key}"
        else:
            title_str = title
        ax.set_title(title_str)
        ax.legend(loc=legend_loc)

        if save is not None:
            save_fig(fig, save)

        fig.show()

    def plot_real_spectrum(
        self,
        dpi: int = 100,
        figsize: Optional[Tuple[float, float]] = None,
        legend_loc: Optional[str] = None,
        save: Optional[Union[str, Path]] = None,
    ) -> None:
        """
        Plot the real part of the top eigenvalues.

        Params
        ------
        dpi
            Dots per inch.
        figsize
            Size of the figure.
        save
            Filename where to save the plots. If `None`, just shows the plot.

        Returns
        -------
        None
            Nothing, just plots the spectrum.
        """

        if self._eig is None:
            logg.warning(
                "No eigendecomposition found, computing with default parameters"
            )
            self.compute_eig()

        # Obtain the eigendecomposition, create the color code
        D, params = self._eig["D"], self._eig["params"]
        D_real, D_imag = D.real, D.imag
        ixs = np.arange(len(D))
        mask = D_imag == 0

        # plot the top eigenvalues
        fig, ax = plt.subplots(nrows=1, ncols=1, dpi=dpi, figsize=figsize)
        ax.scatter(ixs[mask], D_real[mask], marker="o", label="Real eigenvalue")
        ax.scatter(ixs[~mask], D_real[~mask], marker="o", label="Complex eigenvalue")

        # add dashed line for the eigengap, ticks, labels, title and legend
        ax.axvline(self._eig["eigengap"], label="Eigengap", ls="--")
        ax.set_xticks(range(len(D)))
        ax.set_xlabel("index")
        ax.set_ylabel("Re($\lambda_i$)")
        key = "real part" if params["which"] == "LR" else "magnitude"
        ax.set_title(
            f"Real part of top {params['k']} eigenvalues according to their {key}"
        )
        ax.legend(loc=legend_loc)

        if save is not None:
            save_fig(fig, save)

        fig.show()

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

    def set_approx_rcs(
        self,
        rc_labels: Union[Series, Dict[Any, Any]],
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
            If a key to cluster labels is given, `approx_rcs` will ge associated with these for naming and colors.
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
            attr_key="_approx_rcs",
            pretty_attr_key="approx_recurrent_classes",
            cat_key=self._rc_key,
            add_to_existing_error_msg="Compute approximate recurrent classes first as `.compute_approx_rcs()`.",
            categories=rc_labels,
            cluster_key=cluster_key,
            en_cutoff=en_cutoff,
            p_thresh=p_thresh,
            add_to_existing=add_to_existing,
        )

    def compute_approx_rcs(
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
            If a key to cluster labels is given, `approx_rcs` will ge associated with these for naming and colors.
        n_clusters_kmeans
            If `None`, this is set to :paramref:`use` `+ 1`.
        n_neighbors_louvain
            If we use `'louvain'` for clustering cells, we need to build a KNN graph.
            This is the K parameter for that, the number of neighbors for each cell.
        resolution_louvain
            Resolution parameter from the `louvain` algorithm. Should be chosen relatively small.
        n_matches_min
            Filters out cells which don't have at leas n_matches_min neighbors from the same class.
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
            If the test returns a positive statistic and a p-value smaller than :paramref:`p_thresh`, a warning will be issued.

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
        probs = self._compute_approx_rcs_prob(use)
        self._approx_rcs_probs = probs
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
        logg.debug(
            f"DEBUG: Using `{use}` eigenvectors, basis `{basis!r}` and method `{method!r}` for clustering"
        )
        labels = _cluster_X(
            X,
            method=method,
            n_clusters_kmeans=n_clusters_kmeans,
            percentile=percentile,
            use=use,
            n_neighbors_louvain=n_neighbors_louvain,
            resolution_louvain=resolution_louvain,
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

        self.set_approx_rcs(
            rc_labels=rc_labels,
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

    def plot_approx_rcs(self, cluster_key: Optional[str] = None, **kwargs) -> None:
        """
        Plots the approximate recurrent classes in a given embedding.

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

        if self._approx_rcs is None:
            raise RuntimeError(
                "Compute approximate recurrent classes first as `.compute_approx_rcs()`"
            )

        self._adata.obs[self._rc_key] = self._approx_rcs

        # check whether the length of the color array matches the number of clusters
        color_key = _colors(self._rc_key)
        if color_key in self._adata.uns and len(self._adata.uns[color_key]) != len(
            self._approx_rcs.cat.categories
        ):
            del self._adata.uns[_colors(self._rc_key)]
            self._approx_rcs_colors = None

        color = self._rc_key if cluster_key is None else [cluster_key, self._rc_key]
        scv.pl.scatter(self._adata, color=color, **kwargs)

        if color_key in self._adata.uns:
            self._approx_rcs_colors = self._adata.uns[color_key]

    def compute_lin_probs(
        self,
        keys: Optional[Sequence[str]] = None,
        check_irred: bool = False,
        norm_by_frequ: bool = False,
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

        Returns
        -------
        None
            Nothing, but updates the following fields: :paramref:`lineage_probabilities`, :paramref:`diff_potential`.
        """

        if self._approx_rcs is None:
            raise RuntimeError(
                "Compute approximate recurrent classes first as `.compute_approx_rcs()`"
            )
        if keys is not None:
            keys = sorted(set(keys))

        # Note: There are three relevant data structures here
        # - self.approx_rcs: pd.Series which contains annotations for approx rcs. Associated colors in
        #   self.approx_rcs_colors
        # - self.lin_probs: Linage object which contains the lineage probabilities with associated names and colors
        # - _approx_rcs: pd.Series, temporary copy of self.approx rcs used in the context of this function. In this
        #   copy, some approx_rcs may be removed or combined with others
        start = logg.info("Computing absorption probabilities")

        # we don't expect the abs. probs. to be sparse, therefore, make T dense. See scipy docs about sparse lin solve.
        t = self._T.A if self._is_sparse else self._T

        # colors are created in `compute_approx_rcs`, this is just in case
        self._check_and_create_colors()

        # process the current annotations according to `keys`
        approx_rcs_, colors_ = _process_series(
            series=self._approx_rcs, keys=keys, colors=self._approx_rcs_colors
        )

        #  create empty lineage object
        if self._lin_probs is not None:
            logg.debug("DEBUG: Overwriting `.lin_probs`")
        self._lin_probs = Lineage(
            np.empty((1, len(colors_))),
            names=approx_rcs_.cat.categories,
            colors=colors_,
        )

        # warn in case only one state is left
        keys = list(approx_rcs_.cat.categories)
        if len(keys) == 1:
            logg.warning(
                "There is only one recurrent class, all cells will have probability 1 of going there"
            )

        # create arrays of all recurrent and transient indices
        mask = np.repeat(False, len(approx_rcs_))
        for cat in approx_rcs_.cat.categories:
            mask = np.logical_or(mask, approx_rcs_ == cat)
        rec_indices, trans_indices = np.where(mask)[0], np.where(~mask)[0]

        # create Q (restriction transient-transient), S (restriction transient-recurrent) and I (Q-sized identity)
        q = t[trans_indices, :][:, trans_indices]
        s = t[trans_indices, :][:, rec_indices]
        eye = np.eye(len(trans_indices))

        if check_irred:
            if self._is_irreducible is None:
                self.compute_partition()
            if not self._is_irreducible:
                logg.warning("Restriction Q is not irreducible")

        # compute abs probs. Since we don't expect sparse solution, dense computation is faster.
        logg.debug("DEBUG: Solving the linear system to find absorption probabilities")
        abs_states = solve(eye - q, s)

        # aggregate to class level by summing over columns belonging to the same approx_rcs
        approx_rc_red = approx_rcs_[mask]
        rec_classes_red = {
            key: np.where(approx_rc_red == key)[0]
            for key in approx_rc_red.cat.categories
        }
        _abs_classes = np.concatenate(
            [
                np.sum(abs_states[:, rec_classes_red[key]], axis=1)[:, None]
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
            cl: np.where(approx_rcs_ == cl) for cl in approx_rcs_.cat.categories
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
        self._adata.obs[f"{self._lin_key}_dp"] = self._dp
        self._adata.uns[_lin_names(self._lin_key)] = self._lin_probs.names
        self._adata.uns[_colors(self._lin_key)] = self._lin_probs.colors

        logg.info("    Finish", time=start)

    def _compute_approx_rcs_prob(
        self, use: Union[Tuple[int], List[int], range]
    ) -> np.ndarray:
        """
        Utility function which computes a global score of being an approximate recurrent class.
        """

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
        n_cats = len(self._approx_rcs.cat.categories)
        if self._approx_rcs_colors is None:
            color_key = _colors(self._rc_key)
            if color_key in self._adata.uns and n_cats == len(
                self._adata.uns[color_key]
            ):
                logg.debug("DEBUG: Loading colors from `.adata` object")
                self._approx_rcs_colors = _convert_to_hex_colors(
                    self._adata.uns[color_key]
                )
            else:
                self._approx_rcs_colors = _create_categorical_colors(n_cats)
                self._adata.uns[_colors(self._rc_key)] = self._approx_rcs_colors
        elif len(self._approx_rcs_colors) != n_cats:
            self._approx_rcs_colors = _create_categorical_colors(n_cats)
            self._adata.uns[_colors(self._rc_key)] = self._approx_rcs_colors

    def copy(self) -> "MarkovChain":
        """
        Return a copy of itself.
        """

        kernel = copy(self.kernel)  # doesn't copy the adata object
        mc = MarkovChain(
            kernel, self.adata.copy(), inplace=False, read_from_adata=False
        )

        mc._is_irreducible = self.irreducible
        mc._rec_classes = copy(self._rec_classes)
        mc._trans_classes = copy(self._trans_classes)

        mc._eig = deepcopy(self.eigendecomposition)
        mc._lin_probs = copy(self.lineage_probabilities)
        mc._dp = copy(self.diff_potential)

        mc._approx_rcs = copy(self.approx_recurrent_classes)
        mc._approx_rcs_probs = copy(self.approx_recurrent_classes_probabilities)
        mc._approx_rcs_colors = copy(self._approx_rcs_colors)

        mc._G2M_score = copy(self._G2M_score)
        mc._S_score = copy(self._S_score)

        mc._g2m_key = self._g2m_key
        mc._s_key = self._s_key
        mc._key_added = self._key_added

        return mc

    @property
    def irreducible(self) -> Optional[bool]:
        """
        Whether the Markov chain is irreducible or not.
        """
        return self._is_irreducible

    @property
    def recurrent_classes(self) -> Optional[List[List[Any]]]:
        """
        The recurrent classes of the Markov chain.
        """
        return self._rec_classes

    @property
    def transient_classes(self) -> Optional[List[List[Any]]]:
        """
        The recurrent classes of the Markov chain.
        """
        return self._trans_classes

    @property
    def approx_recurrent_classes(self) -> Series:
        """
        The approximate recurrent classes, where `NaN` marks cells which are transient.
        """
        return self._approx_rcs

    @property
    def approx_recurrent_classes_probabilities(self):
        """
        Probabilities of cells belonging to the approximate recurrent classes.
        """
        return self._approx_rcs_probs
