# -*- coding: utf-8 -*-
from cellrank.tools.kernels._kernel import KernelExpression
from typing import Optional, Tuple, Sequence, List, Any, Union, Dict, Iterable
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import scvelo as scv

from anndata import AnnData
from itertools import combinations
from pandas import Series, DataFrame, to_numeric
from pandas.api.types import is_categorical_dtype
from scanpy import logging as logg
from scipy.linalg import solve
from scipy.sparse import issparse
from scipy.sparse.linalg import eigs
from scipy.stats import zscore, entropy, ranksums


from cellrank.tools._lineage import Lineage
from cellrank.tools._constants import (
    Direction,
    RcKey,
    LinKey,
    Prefix,
    _probs,
    _colors,
    _lin_names,
)
from cellrank.tools._utils import (
    _complex_warning,
    _cluster_X,
    _get_connectivities,
    _normalize,
    _eigengap,
    _filter_cells,
    _make_cat,
    _create_colors,
    _convert_to_hex_colors,
    _vec_mat_corr,
    _create_categorical_colors,
    _compute_mean_color,
    _convert_to_categorical_series,
    partition,
    save_fig,
)


class MarkovChain:
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
        if adata is not None:
            self._adata = adata if inplace else adata.copy()
        else:
            logg.debug("DEBUG: Loading `adata` object from kernel.")
            self._adata = kernel.adata if inplace else kernel.adata.copy()

        if kernel.backward:
            self._direction: Direction = Direction.BACKWARD
            self._rc_key: str = str(RcKey.BACKWARD)
            self._lin_key: str = str(LinKey.BACKWARD)
            self._prefix: str = str(Prefix.BACKWARD)
        else:
            self._direction: Direction = Direction.FORWARD
            self._rc_key: str = str(RcKey.FORWARD)
            self._lin_key: str = str(LinKey.FORWARD)
            self._prefix: str = str(Prefix.FORWARD)

        # import transition matrix and parameters
        if kernel.transition_matrix is None:
            logg.debug("DEBUG: Computing transition matrix using default parameters.")
            kernel.compute_transition_matrix()
        kernel.write_to_adata(key_added=key_added)

        self._kernel = kernel
        self._T = kernel.transition_matrix
        self._is_sparse = issparse(self._T)
        self._n_states = self._T.shape[0]
        if self._n_states != self._adata.n_obs:
            raise ValueError("Number of states and number of cells differ.")

        self._is_irreducible = None
        self._rec_classes = None
        self._trans_classes = None

        # read eig, approx_rcs and lin_probs from adata if present
        self._eig, self._approx_rcs, self._approx_rcs_colors, self._lin_probs, self._dp, self._G2M_score, self._S_score, self._approx_rcs_probs = (
            [None] * 8
        )

        if read_from_adata:
            logg.debug(
                "DEBUG: Reading `eig`, `approx_rcs` and `lin_probs` from `adata` object"
            )
            self._read_from_adata(g2m_key, s_key)

    def _read_from_adata(
        self, g2m_key: Optional[str] = None, s_key: Optional[str] = None
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
        Compute eigendecomposition of transition matrix.

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

        logg.info("Computing eigendecomposition of transition matrix")
        if self._is_sparse:
            logg.debug(f"DEBUG: Computing top `{k}` eigenvalues for sparse matrix")
            D, V_l = eigs(self._T.T, k=k, which=which)
            _, V_r = eigs(self._T, k=k, which=which)
        else:
            logg.warning(
                "This transition matrix is not sparse, computing full eigendecomposition"
            )
            D, V_l = np.linalg.eig(self._T.T)
            _, V_r = np.linalg.eig(self._T)

        # Sort the eigenvalues and eigenvectors and take the real part
        logg.debug("DEBUG: Sorting eigenvalues by their real part")
        p = np.flip(np.argsort(D.real))
        D, V_l, V_r = D[p], V_l[:, p], V_r[:, p]
        e_gap = _eigengap(D.real, alpha)

        # write to class and AnnData object
        if self._eig is not None:
            logg.debug("DEBUG: Overwriting `.eig`")
        else:
            logg.debug(f"DEBUG: Adding `.eig` and `adata.uns['eig_{self._direction}']`")

        eig_dict = {
            "D": D,
            "V_l": V_l,
            "V_r": V_r,
            "eigengap": e_gap,
            "params": {"which": which, "k": k, "alpha": alpha},
        }
        self._eig = eig_dict
        self._adata.uns[f"eig_{self._direction}"] = eig_dict

    def plot_eig(
        self,
        dpi: int = 100,
        figsize: Optional[Tuple[float, float]] = (5, 5),
        legend_loc: Optional[str] = None,
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
        ax.set_title(f"Top {params['k']} eigenvalues according to their {key}")
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
        use: Optional[Union[int, tuple, list]] = None,
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
            Which or how many eigenvectors to be plotted. If None, it will be chosen by `eigengap`.
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

        # set the direction
        side = "left" if left else "right"

        # get the eigendecomposition
        D, V = self._eig["D"], self._eig[f"V_{side[0]}"]

        # check whether dimensions are consistent
        if self._adata.n_obs != V.shape[0]:
            raise ValueError(
                "Cell number inconsistent with dimensions of eigendecomposition."
            )

        # subset eigenvectors
        if use is None:
            use = self._eig["eigengap"] + 1  # add one because first e-vec has index 0
        if isinstance(use, int):
            if use > 15:
                logg.warning(
                    f"Too many eigenvectors ({use}) at once. Showing the first `15`"
                )
                use = 15
            use = list(range(use))
        elif not isinstance(use, (tuple, list, range)):
            raise TypeError(
                f"Argument `use` must be either `int`, `tuple`, `list` or `range`,"
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

        # retrieve the eigendecomposition, check for imaginary values and issue warnings
        D, V = D[use], V[:, use]
        V_ = _complex_warning(V, use, use_imag=use_imag)

        # take absolute value
        if abs_value:
            V_ = np.abs(V_)

        if cluster_key is not None:
            color = [cluster_key] + [v for v in V_.T]
        else:
            color = [v for v in V_.T]

        # actual plotting with scvelo
        logg.debug(f"DEBUG: Showing `{use}` {side} eigenvectors")
        scv.pl.scatter(
            self._adata,
            color=color,
            title=[f"$\lambda_{i}$={d:.02f}" for i, d in zip(use, D)],
            **kwargs,
        )

    def set_approx_rcs(
        self,
        rc_labels: Union[Series, Dict[Any, Any]],
        cluster_key: Optional[str] = None,
        en_cutoff: Optional[float] = None,
        p_thresh: Optional[float] = None,
    ):
        """
        Set the approximate recurrent classes, if they are known a priori.

        Params
        ------
        rc_labels
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
            If the test returns a positive statistic and a p-value smaller than :paramref:`p_thresh`, a warning will be issued.

        Returns
        -------
        None
            Nothing, but updates the following fields: :paramref:`approx_recurrent_classes`.
        """

        if isinstance(rc_labels, dict):
            rc_labels = _convert_to_categorical_series(
                rc_labels, list(self.adata.obs_names)
            )
        if not is_categorical_dtype(rc_labels):
            raise TypeError(
                f"Approximate recurrent classes must be `categorical`, found `{type(rc_labels).__name__}`."
            )

        if cluster_key is not None:
            logg.debug(f"DEBUG: Creating colors based on `{cluster_key}`")
            approx_rcs_names, self._approx_rcs_colors = self._get_lin_names_colors(
                rc_labels, cluster_key, en_cutoff
            )
            rc_labels.cat.categories = approx_rcs_names
        else:
            self._approx_rcs_colors = _create_categorical_colors(
                len(rc_labels.cat.categories)
            )

        if p_thresh is not None:
            self._detect_cc_stages(rc_labels, p_thresh=p_thresh)

        # write to class and adata
        if self._approx_rcs is not None:
            logg.debug("DEBUG: Overwriting `.approx_rcs`")

        self._approx_rcs = rc_labels
        self._adata.obs[self._rc_key] = self._approx_rcs
        self._adata.uns[_colors(self._rc_key)] = self._approx_rcs_colors

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
        self._compute_approx_rcs_prob(use)

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
            rc_labels, cluster_key=cluster_key, en_cutoff=en_cutoff, p_thresh=p_thresh
        )

        logg.info(
            f"Adding `adata.uns[{_colors(self._rc_key)!r}]`\n"
            f"       `adata.obs[{_probs(self._rc_key)!r}]`\n"
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
        n_cats = len(self._approx_rcs.cat.categories)
        if self._approx_rcs_colors is None:
            color_key = _colors(self._rc_key)
            if color_key in self._adata.uns and n_cats == len(
                self._adata.uns[color_key]
            ):
                logg.debug("DEBUG: Loading colors from `.adata` object")
                self._approx_rcs_colors = self._adata.uns[color_key]
            else:
                self._approx_rcs_colors = _create_categorical_colors(n_cats)
                self._adata.uns[_colors(self._rc_key)] = self._approx_rcs_colors
        elif len(self._approx_rcs_colors) != n_cats:
            self._approx_rcs_colors = _create_categorical_colors(n_cats)
            self._adata.uns[_colors(self._rc_key)] = self._approx_rcs_colors

        # initialize empty lineage object with names and colors from recurrent classes
        if self._lin_probs is not None:
            logg.debug("DEBUG: Overwriting `.lin_probs`")
        rc_names = list(self._approx_rcs.cat.categories)
        self._lin_probs = Lineage(
            np.empty((1, len(rc_names))), names=rc_names, colors=self._approx_rcs_colors
        )

        # if keys are given, remove some rc's and combine others
        if keys is None:
            approx_rcs_ = self._approx_rcs
        else:
            logg.debug(f"DEBUG: Combining recurrent classes according to `{keys}`")
            approx_rcs_ = self._prep_rc_classes(keys)
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

        # for recurrent states, take the maximum prob among all transient states
        abs_classes = np.zeros((self._n_states, len(rec_classes_red)))
        rec_classes_full = {
            cl: np.where(approx_rcs_ == cl) for cl in approx_rcs_.cat.categories
        }
        for col, cl_indices in enumerate(rec_classes_full.values()):
            abs_classes[trans_indices, col] = _abs_classes[:, col]
            abs_classes[cl_indices, col] = np.max(_abs_classes[:, col])

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

    def plot_lin_probs(
        self,
        lineages: Optional[Union[str, Iterable[str]]] = None,
        cluster_key: Optional[str] = None,
        mode: str = "embedding",
        time_key: str = "latent_time",
        cmap: Union[str, matplotlib.colors.ListedColormap] = cm.viridis,
        **kwargs,
    ) -> None:
        """
        Plots the absorption probabilities in the given embedding.

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
        cmap
            Colormap to use.
        kwargs
            Keyword arguments for :func:`scvelo.pl.scatter`.

        Returns
        -------
        None
            Nothing, just plots the absorption probabilities.
        """

        if self._lin_probs is None:
            raise RuntimeError(
                "Compute lineage probabilities first as `.compute_lin_probs()`."
            )
        if isinstance(lineages, str):
            lineages = [lineages]

        if lineages is None:
            lineages = self._lin_probs.names
            A = self._lin_probs.X
        else:
            for lineage in lineages:
                if lineage not in self._lin_probs.names:
                    raise ValueError(
                        f"Invalid lineage name `{lineages!r}`. Valid options are `{list(self._lin_probs.names)}`."
                    )
            A = self._lin_probs[lineages].X

        if mode == "time":
            if time_key not in self._adata.obs.keys():
                raise KeyError(f"Time key `{time_key}` not in `adata.obs`.")
            t = self._adata.obs[time_key]
            cluster_key = None

        rc_titles = [f"{self._prefix} {rc}" for rc in lineages] + [
            "Differentiation Potential"
        ]

        if cluster_key is not None:
            color = [cluster_key] + [a for a in A.T] + [self._dp]
            titles = [cluster_key] + rc_titles
        else:
            color = [a for a in A.T] + [self._dp]
            titles = rc_titles

        if mode == "embedding":
            scv.pl.scatter(
                self._adata, color=color, title=titles, color_map=cmap, **kwargs
            )
        elif mode == "time":
            xlabel, ylabel = (
                list(np.repeat(time_key, len(titles))),
                list(np.repeat("probability", len(titles) - 1)) + ["entropy"],
            )
            scv.pl.scatter(
                self._adata,
                x=t,
                color_map=cmap,
                y=[a for a in A.T] + [self._dp],
                title=titles,
                xlabel=time_key,
                ylabel=ylabel,
                **kwargs,
            )
        else:
            raise ValueError(
                f"Invalid mode `{mode!r}`. Valid options are: `'embedding', 'time'`."
            )

    def compute_lineage_drivers(
        self,
        lin_names: Optional[Sequence] = None,
        cluster_key: Optional[str] = "louvain",
        clusters: Optional[Sequence] = None,
        layer: str = "X",
        use_raw: bool = True,
        inplace: bool = True,
    ):
        """
        Compute driver genes per lineage.

        Correlates gene expression with lineage probabilities, for a given lineage and set of clusters.
        Often, it makes sense to restrict this to a set of clusters which are relevant for the lineage under consideration.

        Params
        --------
        lin_keys
            Either a set of lineage names from :paramref:`lineage_probabilities` `.names` or None,
            in which case all lineages are considered.
        cluster_key
            Key from :paramref:`adata` `.obs` to obtain cluster annotations.
            These are considered for :paramref:`clusters`.
        clusters
            Restrict the correlations to these clusters.
        layer
            Key from :paramref:`adata` `.layers`.
        use_raw
            Whether or not to use :paramref:`adata` `.raw` to correlate gene expression.
            If using a layer other than `.X`, this must be set to `False`.

        Returns
        --------
        :class:`pandas.DataFrame` or :class:`NoneType`
            Writes to :paramref:`adata` `.var` or :paramref:`adata` `.raw.var`,
            depending on the value of :paramref:`use_raw`.
            For each lineage specified, a key is added to `.var` and correlations are saved there.

            Returns `None` if :paramref:`inplace` `=True`, otherwise a dataframe.
        """

        # check that lineage probs have been computed
        if self._lin_probs is None:
            raise RuntimeError(
                "Compute lineage probabilities first as `.compute_lin_probs()`."
            )

        # check all lin_keys exist in self.lin_names
        if lin_names is not None:
            _ = self._lin_probs[lin_names]
        else:
            lin_names = self._lin_probs.names

        # check the cluster key exists in adata.obs and check that all clusters exist
        if cluster_key is not None and cluster_key not in self._adata.obs.keys():
            raise KeyError(f"Key `{cluster_key!r}` not found in `adata.obs`.")

        if clusters is not None:
            all_clusters = np.array(self._adata.obs[cluster_key].cat.categories)
            cluster_mask = np.array([name not in all_clusters for name in clusters])
            if any(cluster_mask):
                raise KeyError(
                    f"Clusters `{list(np.array(clusters)[cluster_mask])}` not found in "
                    f"`adata.obs[{cluster_key!r}]`."
                )

            subset_mask = np.in1d(self._adata.obs[cluster_key], clusters)
            adata_comp = self._adata[subset_mask].copy()
            lin_probs = self._lin_probs[subset_mask, :]
        else:
            adata_comp = self._adata.copy()
            lin_probs = self._lin_probs

        # check that the layer exists, and that use raw is only used with layer X
        if layer != "X":
            if layer not in self._adata.layers:
                raise KeyError(f"Layer `{layer!r}` not found in `adata.layers`.")
            if use_raw:
                raise ValueError("For `use_raw=True`, layer must be 'X'.")
            data = adata_comp.layers[layer]
            var_names = adata_comp.var_names
        else:
            if use_raw and self._adata.raw is None:
                raise AttributeError("No raw attribute set")
            data = adata_comp.raw.X if use_raw else adata_comp.X
            var_names = adata_comp.raw.var_names if use_raw else adata_comp.var_names

        start = logg.info(
            f"Computing correlations for lineages `{lin_names}` restricted to clusters `{clusters}` in "
            f"layer `{layer}` with `use_raw={use_raw}`"
        )

        # loop over lineages
        lin_corrs = {}
        for lineage in lin_names:
            y = lin_probs[:, lineage].X.squeeze()
            correlations = _vec_mat_corr(data, y)

            if inplace:
                if use_raw:
                    self._adata.raw.var[f"{self._prefix} {lineage} corr"] = correlations
                else:
                    self._adata.var[f"{self._prefix} {lineage} corr"] = correlations
            else:
                lin_corrs[lineage] = correlations

        if not inplace:
            return DataFrame(lin_corrs, index=var_names)

        field = "raw.var" if use_raw else "var"
        logg.info(
            f"Adding gene correlations to `.adata.{field}`\n    Finish", time=start
        )

    def _get_lin_names_colors(
        self, rc_labels: Series, cluster_key: str, en_cutoff: Optional[float]
    ) -> Tuple[Series, List[Any]]:
        """
        Utility function to map approx_rcs to ts_clusters.

        Params
        ------
        rc_labels
            Raw labels for the approx_rcs. These are labels like [0, 1, 2, 3]. This function is meant to associate
            these uninformative labels with pre-computed clusters based on transcriptomic similarities, using e.g.
            louvain, leiden or k-means clustering.
        cluster_key
            String to a key from :paramref:`adata` `.obs`. This is how the pre-computed cluster labels are accessed.
        en_cutoff
            Threshold for defining when to associate an approximate rc with a pre-computed cluster. Imagine a
            very fine pre-computed louvain clustering and an approximate rc spanning two louvain clusters, with approx.
            50% of the cells coming from either cluster. What name should this approx. rc get? I thins case, we
            should call it 'Unknown', since we don't know. The entropy threshold decides when we call an approximate rc
            'Unknown'

        Returns
        -------
        :class:`pandas.Series`, :class:`list`
            Array-like names and colors for the approximate RC's.
            These link back to the underlying pre-computed clusters.
        """

        if cluster_key not in self._adata.obs:
            raise KeyError(f"Cluster key `{cluster_key!r}` not found in `.adata.obs`.")

        # create dataframe to store approx_rc associtation with ts_clusters
        ts_labels = self._adata.obs[cluster_key]
        approx_rc_clusters = rc_labels.cat.categories
        ts_clusters = ts_labels.cat.categories
        rc_df = DataFrame(None, index=approx_rc_clusters, columns=ts_clusters)

        # populate the df - compute the overlap
        for cl in approx_rc_clusters:
            row = [
                np.sum(ts_labels.loc[np.array(rc_labels == cl)] == key)
                for key in ts_clusters
            ]
            rc_df.loc[cl] = row
        rc_df = rc_df.apply(to_numeric)

        # label endpoints and add uncertainty through entropy
        rc_df["entropy"] = entropy(rc_df.T)
        rc_df["name"] = rc_df.T.idxmax()

        # add cluster colors
        # clusters_colors = self.lin_probs.colors  # _convert_to_hex_colors(self._adata.uns[_colors(cluster_key)])
        clusters_colors = _convert_to_hex_colors(self._adata.uns[_colors(cluster_key)])
        names = rc_df["name"]
        lin_colors = []
        for name in names:
            mask = ts_clusters == name
            color = np.array(clusters_colors)[mask][0]
            lin_colors.append(color)
        rc_df["color"] = lin_colors

        # make rc_labels unique
        names = Series(rc_df["name"], dtype="category")
        lin_colors = Series(lin_colors, dtype="category")  # colors must be hex
        frequ = {key: np.sum(names == key) for key in names.cat.categories}

        # Create unique names by adding suffixes "..._1, ..._2" etc and unique colors by shifting the original color
        names_new = np.array(names.copy())
        lin_colors_new = np.array(lin_colors.copy())

        for key, value in frequ.items():
            if value == 1:
                continue  # already unique, skip

            # deal with non-unique names
            suffix = list(np.arange(1, value + 1).astype("str"))
            unique_names = [f"{key}_{rep}" for rep in suffix]
            names_new[names == key] = unique_names

            color = rc_df[rc_df["name"] == key]["color"].values[0]
            shifted_colors = _create_colors(color, value, saturation_range=None)
            lin_colors_new[lin_colors == color] = shifted_colors

        rc_df["name"] = names_new
        rc_df["color"] = lin_colors_new

        # issue a warning for rcs with high entropy
        if en_cutoff is not None:
            critical_rcs = list(rc_df.loc[rc_df["entropy"] > en_cutoff, "name"].values)
            if len(critical_rcs) > 0:
                logg.warning(
                    f"The following groups of {self._rc_key} contain many different cell types: `{critical_rcs}`"
                )

        return rc_df["name"], list(rc_df["color"])

    def _prep_rc_classes(self, keys: Sequence[str]) -> Tuple[Series, np.ndarray]:
        """
        Utility function to remove and combine rcs.

        This function takes a list of strings, like ['endpt_1, endpt_2', 'endpt_3'], and compares this with the
        names of the approximate recurrent classes to figure out which ones to leave out and which ones to combine. If
        e.g. the approx rcs had names ['endpt_1', 'endpt_2, 'endpt_3', 'endpt_4'], then this function would combine
        endpoints 1 and 2 and kick out 4.

        Params
        ------
        keys
            Sequence of strings that defines how absorption probabilities should be computed.

        Returns
        -------
        :class:`pandas.Series`
            Categorical updated annotation. Each cell is assigned to either `NaN`
            or one of updated approximate recurrent classes.
        """

        if self._approx_rcs is None:
            raise RuntimeError(
                "Compute approximate recurrent classes first as `.compute_approx_rcs()`"
            )

        # initialize a copy of the approx_rcs Series
        approx_rcs_temp = self._approx_rcs.copy()

        # define a set of keys
        keys_ = {
            tuple((key.strip() for key in rc.strip(" ,").split(","))) for rc in keys
        }

        overlap = [set(ks) for ks in keys_]
        for c1, c2 in combinations(overlap, 2):
            overlap = c1 & c2
            if overlap:
                raise ValueError(f"Found overlapping keys: `{list(overlap)}`.")

        # remove the unused categories, both in approx_rcs_temp as well as in the lineage object
        remaining_cat = [b for a in keys_ for b in a]
        removed_cat = list(set(approx_rcs_temp.cat.categories) - set(remaining_cat))
        approx_rcs_temp.cat.remove_categories(removed_cat, inplace=True)
        original_colors = list(self._lin_probs[remaining_cat].colors)
        original_len = len(original_colors)

        # loop over all indiv. or combined rc's
        lin_colors = {}
        for cat in keys_:
            # if there are more than two keys in this category, combine them
            if len(cat) > 1:
                new_cat_name = " or ".join(cat)
                mask = np.repeat(False, len(approx_rcs_temp))
                for key in cat:
                    mask = np.logical_or(mask, approx_rcs_temp == key)
                    remaining_cat.remove(key)
                approx_rcs_temp.cat.add_categories(new_cat_name, inplace=True)
                remaining_cat.append(new_cat_name)
                approx_rcs_temp[mask] = new_cat_name

                # apply the same to the colors array. We just append new colors at the end
                color_mask = np.in1d(approx_rcs_temp.cat.categories[:original_len], cat)
                colors_merge = np.array(original_colors)[:original_len][color_mask]
                lin_colors[new_cat_name] = _compute_mean_color(colors_merge)
            else:
                lin_colors[cat[0]] = self._lin_probs[cat].colors[0]

        # Since we have just appended colors at the end, we must now delete the unused ones
        approx_rcs_temp.cat.remove_unused_categories(inplace=True)
        approx_rcs_temp.cat.reorder_categories(remaining_cat, inplace=True)

        self._lin_probs = Lineage(
            np.empty((1, len(lin_colors))),
            names=approx_rcs_temp.cat.categories,
            colors=[lin_colors[c] for c in approx_rcs_temp.cat.categories],
        )

        return approx_rcs_temp

    def _detect_cc_stages(self, rc_labels: Series, p_thresh: float = 1e-15) -> None:
        """
        Utility function to detect cell-cycle driven start or endpoints.
        """

        # initialise the groups (start or end clusters) and scores
        groups = rc_labels.cat.categories
        scores = []
        if self._G2M_score is not None:
            scores.append(self._G2M_score)
        if self._S_score is not None:
            scores.append(self._S_score)

        # loop over groups and scores
        for group in groups:
            flag = False
            mask = rc_labels == group
            for score in scores:
                a, b = score[mask], score[~mask]
                result = ranksums(a, b)
                if result.statistic > 0 and result.pvalue < p_thresh:
                    flag = True
            if flag:
                logg.warning(f"Group `{group}` appears to be cell-cycle driven")

    def _compute_approx_rcs_prob(self, use: Union[Tuple[int], List[int], range]):
        """
        Utility function which computes a global score of being an approximate recurrent class
        """

        if self._eig is None:
            raise RuntimeError("Compute eigendecomposition first as `.compute_eig()`")

        # get the truncated eigendecomposition
        V, evals = self._eig["V_l"].real[:, use], self._eig["D"].real[use]

        # shift and scale
        V_pos = np.abs(V)
        V_shifted = V_pos - np.min(V_pos, axis=0)
        V_scaled = V_shifted / np.max(V_shifted, axis=0)

        # check the ranges are correct
        assert np.allclose(np.min(V_scaled, axis=0), 0), "Lower limit it not zero"
        assert np.allclose(np.max(V_scaled, axis=0), 1), "Upper limit is not one"

        # further scale by the eigenvalues
        V_eigs = V_scaled / evals

        # sum over cols and scale
        c_ = np.sum(V_eigs, axis=1)
        c = c_ / np.max(c_)

        self._approx_rcs_probs = c
        self._adata.obs[_probs(self._rc_key)] = c

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
    def eigendecomposition(self) -> Optional[Dict[str, np.ndarray]]:
        """
        A dictionary with following fields:

        - `'D'` eigenvalues of left eigenvectors
        - `'V_l'` left eigenvectors
        - `'V_r'` right eigenvectors
        """
        return self._eig

    @property
    def lineage_probabilities(self) -> Lineage:
        """
        A `numpy`-like array with names and colors, where
        each column represents one lineage.
        """
        return self._lin_probs

    @property
    def diff_potential(self) -> np.ndarray:
        """
        Differentiation potential for each lineage.
        """
        return self._dp

    @property
    def approx_recurrent_classes(self) -> DataFrame:
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

    @property
    def adata(self) -> AnnData:
        """
        Returns
        -------
        :class:`anndata.AnnData`
            The underlying annotated data object.
        """
        return self._adata

    @property
    def kernel(self) -> KernelExpression:
        """
        The underlying kernel expression.
        """
        return self._kernel

    def __len__(self) -> int:
        return self._n_states

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}[n={len(self)}, kernel={repr(self._kernel)}]"

    def __str__(self) -> str:
        return f"{self.__class__.__name__}[n={len(self)}, kernel={str(self._kernel)}]"
