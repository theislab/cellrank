# -*- coding: utf-8 -*-
"""Generalized Perron Cluster Cluster Analysis (GPCCA) module."""

from copy import copy, deepcopy
from types import MappingProxyType
from typing import Any, Dict, List, Tuple, Union, Mapping, Iterable, Optional
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import entropy

import seaborn as sns
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import scvelo as scv
from scanpy import logging as logg
from anndata import AnnData

from cellrank.tools._utils import (
    save_fig,
    _eigengap,
    _fuzzy_to_discrete,
    _convert_lineage_name,
    _generate_random_keys,
    _series_from_one_hot_matrix,
)
from cellrank.tools._colors import _get_black_or_white
from cellrank.tools._lineage import Lineage
from cellrank.tools._constants import Lin, MetaKey, _dp, _probs, _colors, _lin_names
from cellrank.tools.kernels._kernel import KernelExpression
from cellrank.tools.estimators._base_estimator import BaseEstimator
from cellrank._vendor.msmtools.analysis.dense.gpcca import GPCCA as _GPPCA

# whether to remove overlapping cells from both states, or assign them to the most likely clusters
REMOVE_OVERLAP = False


class GPCCA(BaseEstimator):
    """
    Generalized Perron Cluster Cluster Analysis [GPCCA18]_.

    Params
    ------
    kernel
        Kernel object that stores a transition matrix.
    adata : :class:`anndata.AnnData`
        Optional annotated data object. If given, precomputed lineages can be read in from this.
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
        super().__init__(
            kernel,
            adata,
            inplace=inplace,
            read_from_adata=read_from_adata,
            g2m_key=g2m_key,
            s_key=s_key,
            key_added=key_added,
        )
        self._meta_key = (
            str(MetaKey.BACKWARD) if kernel.backward else str(MetaKey.FORWARD)
        )

        self._gpcca: Optional[GPCCA] = None
        self._schur_vectors = None
        self._schur_matrix = None
        self._coarse_T = None
        self._coarse_init_dist = None
        self._coarse_stat_dist = None

        self._meta_states = None
        self._meta_states_colors = None
        self._meta_lin_probs = None

        self._main_states = None
        self._main_states_probabilities = None

    def compute_eig(
        self,
        k: int = 20,
        which: str = "LM",
        alpha: float = 1,
        ncv: Optional[int] = None,
    ) -> None:
        """
        Compute eigendecomposition of the transition matrix.

        Uses a sparse implementation, if possible, and only computes the top k eigenvalues.

        Params
        ------
        k
            Number of eigenvalues/vectors to compute.
        which
            Eigenvalues are in general complex. `'LR'` - largest real part, `'LM'` - largest magnitude.
        alpha
            Used to compute the `eigengap`. paramref:`alpha` is the weight given
            to the deviation of an eigenvalue from one.
        ncv
            Number of Lanczos vectors generated.

        Returns
        -------
        None
            Nothing, but updates the following fields: paramref:`eigendecomposition`.
        """
        self._compute_eig(k=k, which=which, alpha=alpha, only_evals=True, ncv=ncv)

    def _read_from_adata(
        self, g2m_key: Optional[str] = None, s_key: Optional[str] = None, **kwargs
    ) -> None:
        pass

    def plot_schur_embedding(
        self,
        use: Optional[Union[int, tuple, list]] = None,
        abs_value: bool = False,
        cluster_key: Optional[str] = None,
        **kwargs,
    ) -> None:
        """
        Plot Schur vectors in an embedding.

        Params
        ------
        use
            Which or how many Schur vectors to be plotted. If `None`, all will be chosen.
        abs_value
            Whether to take the absolute value before plotting.
        cluster_key
            Key from :paramref:`adata` `.obs` to plot cluster annotations.
        kwargs
            Keyword arguments for :func:`scvelo.pl.scatter`.

        Returns
        -------
        None
            Nothing, just plots the Schur vectors.
        """

        if self.schur_vectors is None:
            raise RuntimeError(
                "Compute Schur vectors as `.compute_schur()` or `.compute_metastable_states()` with `n_states` > 1."
            )

        self._plot_vectors(
            self.schur_vectors,
            "schur",
            abs_value=abs_value,
            use=use,
            cluster_key=cluster_key,
            **kwargs,
        )

    def plot_schur_matrix(
        self,
        title: Optional[str] = "schur matrix",
        cmap: str = "viridis",
        figsize: Optional[Tuple[float, float]] = None,
        dpi: Optional[float] = 80,
        save: Optional[Union[str, Path]] = None,
        **kwargs,
    ):
        """
        Plot the Schur matrix.

        title
            Title of the figure.
        cmap
            Colormap to use.
        figsize
            Size of the figure.
        dpi
            Dots per inch.
        save
            Filename where to save the plots. If `None`, just shows the plot.

        Returns
        -------
        None
            Nothing, just plots the Schur matrix.
        """

        if self._schur_matrix is None:
            raise RuntimeError(
                "Compute Schur matrix first as `.compute_schur()` or "
                "`.compute_metastable_states()` with `n_states` > 1."
            )

        fig, ax = plt.subplots(
            figsize=self._schur_matrix.shape if figsize is None else figsize, dpi=dpi
        )

        divider = make_axes_locatable(ax)  # square=True make the colorbar a bit bigger
        cbar_ax = divider.append_axes("right", size="2.5%", pad=0.05)

        mask = np.zeros_like(self._schur_matrix, dtype=np.bool)
        mask[np.tril_indices_from(mask, k=-1)] = True
        mask[~np.isclose(self._schur_matrix, 0.0)] = False

        vmin, vmax = (
            np.min(self._schur_matrix[~mask]),
            np.max(self._schur_matrix[~mask]),
        )

        kwargs["fmt"] = kwargs.get("fmt", "0.2f")
        sns.heatmap(
            self._schur_matrix,
            cmap=cmap,
            square=True,
            annot=True,
            vmin=vmin,
            vmax=vmax,
            cbar_ax=cbar_ax,
            mask=mask,
            xticklabels=[],
            yticklabels=[],
            ax=ax,
            **kwargs,
        )

        ax.set_title(title)

        if save is not None:
            save_fig(fig, path=save)

    def compute_schur(
        self,
        n_components: int = 10,
        initial_distribution: Optional[np.ndarray] = None,
        method: str = "krylov",
        which: str = "LM",
        alpha: float = 1,
    ):
        """
        Compute the Schur decomposition.

        Params
        ------
        n_components
            Number of vectors to compute.
        initial_distribution
            Input probability distribution over all cells. If `None`, uniform is chosen.
        method
            Method for calculating the Schur vectors. Valid options are: `'krylov'` or `'brandts'`.
            For benefits of each method, see :class:`msmtools.analysis.dense.gpcca.GPCCA`. The former is
            an iterative procedure that computes a partial, sorted Schur decomposition for large, sparse
            matrices whereas the latter computes a full sorted Schur decomposition of a dense matrix.
        which
            Eigenvalues are in general complex. `'LR'` - largest real part, `'LM'` - largest magnitude.
        alpha
            Used to compute the `eigengap`. paramref:`alpha` is the weight given
            to the deviation of an eigenvalue from one.

        Returns
        -------
        None
            Nothing, but updates the following fields:

                - :paramref:`schur_vectors`
        """

        if n_components < 2:
            raise ValueError(
                f"Number of components must be `>=2`, found `{n_components}`."
            )

        self._gpcca = _GPPCA(self._T, eta=initial_distribution, z=which, method=method)
        try:
            self._gpcca._do_schur_helper(n_components)
        except ValueError:
            logg.warning(
                f"Using {n_components} components would split a block of complex conjugates. "
                f"Increasing `n_components` to {n_components + 1}"
            )
            self._gpcca._do_schur_helper(n_components + 1)

        self._write_eig_to_adata(
            {
                "D": self._gpcca.eigenvalues,
                "eigengap": _eigengap(self._gpcca.eigenvalues, alpha),
                "params": {
                    "which": which,
                    "k": len(self._gpcca.eigenvalues),
                    "alpha": alpha,
                },
            }
        )
        # make it available for plotting
        self._schur_vectors = self._gpcca.X
        self._schur_matrix = self._gpcca.R

    def _compute_meta_states_1_state(
        self,
        n_cells: int,
        cluster_key: Optional[str],
        en_cutoff: Optional[float],
        p_thresh: float,
    ) -> None:
        start = logg.info("Computing metastable states")
        logg.warning("For `n_states=1`, stationary distribution is computed")

        self._compute_eig(only_evals=False, which="LM")
        stationary_dist = self.eigendecomposition["stationary_dist"]

        self._assign_metastable_states(
            memberships=stationary_dist[:, None],
            n_cells=n_cells,
            cluster_key=cluster_key,
            p_thresh=p_thresh,
            en_cutoff=en_cutoff,
        )

        (
            self._lin_probs,
            self._schur_vectors,
            self._coarse_T,
            self._coarse_init_dist,
            self._coarse_stat_dist,
            self._schur_matrix,
        ) = [None] * 6

        logg.info("Adding `.metastable_states`\n    Finish", time=start)

    def _get_n_states_from_minchi(
        self, n_states: Union[Tuple[int, int], List[int], Dict[str, int]]
    ) -> int:
        if not isinstance(n_states, (dict, tuple, list)):
            raise TypeError(
                f"Expected `n_states` to be either `dict`, `tuple` or a `list`, "
                f"found `{type(n_states).__name__}`."
            )
        if len(n_states) != 2:
            raise ValueError(
                f"Expected `n_states` to be of size `2`, found `{len(n_states)}`."
            )

        minn, maxx = (
            (n_states["min"], n_states["max"])
            if isinstance(n_states, dict)
            else n_states
        )
        if minn <= 1:
            raise ValueError(f"Minimum value must be > 1, found `{minn}`.")
        elif minn == 2:
            logg.warning(
                "In most cases, 2 clusters will always be optimal. "
                "If you really expect 2 clusters, use `n_states=2`.\nSetting the minimum to 3"
            )
            minn = 3

        logg.debug(f"DEBUG: Calculating minChi within interval [{minn}, {maxx}]")
        return int(np.arange(minn, maxx)[np.argmax(self._gpcca.minChi(minn, maxx))])

    def compute_metastable_states(
        self,
        n_states: Optional[
            Union[int, Tuple[int, int], List[int], Dict[str, int]]
        ] = None,
        use_min_chi: bool = False,
        n_cells: Optional[int] = 30,
        cluster_key: str = None,
        en_cutoff: Optional[float] = 0.7,
        p_thresh: float = 1e-15,
    ):
        """
        Compute the metastable states.

        Params
        ------
        n_states
            Number of metastable states. If `None`, use the `eigengap` heuristic.
        use_min_chi
            Whether to use :meth:`msmtools.analysis.dense.gpcca.GPCCA.minChi` to calculate the number of metastable
            states. If `True`, :paramref:`n_states` corresponds to an interval `[min, max]` inside of which
            the potentially optimal number of metastable states is searched.
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

        Returns
        -------
        None
            Nothing, but updates the following fields:

                - :paramref:`schur_vectors`
                - :paramref:`coarse_T`
                - :paramref:`coarse_stationary_distribution`
        """

        was_from_eigengap = False
        if n_states is None:
            if self.eigendecomposition is None:
                raise RuntimeError(
                    "Compute eigendecomposition first as `.compute_eig()` or `.compute_schur()`."
                )
            was_from_eigengap = True
            n_states = self.eigendecomposition["eigengap"] + 1
            logg.info(f"Using `{n_states}` states based on eigengap")

        if n_states == 1:
            self._compute_meta_states_1_state(
                n_cells=n_cells,
                cluster_key=cluster_key,
                p_thresh=p_thresh,
                en_cutoff=en_cutoff,
            )
            return

        if self._gpcca is None:
            if was_from_eigengap:
                logg.warning(
                    f"Number of states ({n_states}) was automatically determined by `eigengap` "
                    "but no Schur decomposition was found. Computing with default parameters"
                )
                self.compute_schur(n_states + 1)
            else:
                raise RuntimeError(
                    "Compute Schur decomposition first as `.compute_schur()`."
                )

        if use_min_chi:
            n_states = self._get_n_states_from_minchi(n_states)
        elif not isinstance(n_states, int):
            raise ValueError(
                f"Expected `n_states` to be an integer when `use_min_chi=False`, found `{type(n_states).__name__}`."
            )

        if self._gpcca.X.shape[1] < n_states:
            logg.warning(
                f"Requested more metastable states ({n_states}) than available "
                f"Schur vectors ({self._gpcca.X.shape[1]}). "
                f"Recomputing the decomposition"
            )

        start = logg.info("Computing metastable states")
        try:
            self._gpcca = self._gpcca.optimize(m=n_states)
        except ValueError as e:
            if n_states != self._schur_vectors.shape[1]:
                raise e
            logg.warning(
                f"Unable to perform the optimization using `{self._schur_vectors.shape[1]}` Schur vectors. "
                f"Recomputing the decomposition"
            )
            self.compute_schur(
                n_states + 1,
                initial_distribution=self._gpcca.eta,
                method=self._gpcca.method,
                which=self._gpcca.z,
                alpha=self.eigendecomposition["params"]["alpha"],
            )
            self._gpcca = self._gpcca.optimize(m=n_states)

        self._assign_metastable_states(
            memberships=self._gpcca.memberships,
            n_cells=n_cells,
            cluster_key=cluster_key,
            p_thresh=p_thresh,
            en_cutoff=en_cutoff,
        )

        self._lin_probs = None

        # cache the results and make sure we don't overwrite
        self._schur_vectors = self._gpcca.X
        self._schur_matrix = self._gpcca.R

        names = self._meta_lin_probs.names
        self._coarse_T = pd.DataFrame(
            self._gpcca.coarse_grained_transition_matrix, index=names, columns=names,
        )
        self._coarse_init_dist = pd.Series(
            self._gpcca.coarse_grained_input_distribution, index=names
        )
        # careful here, in case computing the stat. dist failed
        if self._gpcca.coarse_grained_stationary_probability is not None:
            self._coarse_stat_dist = pd.Series(
                self._gpcca.coarse_grained_stationary_probability, index=names,
            )
        else:
            logg.warning("No stationary distribution found in GPCCA object")

        logg.info(
            "Adding `.schur_vectors`\n"
            "       `.metastable_states`\n"
            "       `.coarse_T`\n"
            "       `.coarse_stationary_distribution`\n"
            "    Finish",
            time=start,
        )

    def plot_metastable_states(
        self,
        discrete: bool = False,
        lineages: Optional[Union[str, Iterable[str]]] = None,
        cluster_key: Optional[str] = None,
        mode: str = "embedding",
        time_key: str = "latent_time",
        same_plot: bool = True,
        cmap: Union[str, mpl.colors.ListedColormap] = cm.viridis,
        title: Optional[str] = None,
        **kwargs,
    ) -> None:
        """
        Plot the absorption probabilities of metastable states in the given embedding.

        Params
        ------
        discrete
            Whether to plot the top cells from each linages or the probabilities.
        lineages
            Only show these lineages. If `None`, plot all lineages.
        cluster_key
            Key from :paramref`adata: `.obs` for plotting cluster labels.
        mode
            Can be either `'embedding'` or `'time'`.

            - If `'embedding'`, plot the embedding while coloring in the absorption probabilities.
            - If `'time'`, plot the pseudotime on x-axis and the absorption probabilities on y-axis.
        time_key
            Key from `adata.obs` to use as a pseudotime ordering of the cells.
        same_plot
            Whether to plot the lineages on the same plot using color gradients when :paramref:`mode='embedding'`.
        cmap
            Colormap to use.
        title
            Either `None`, in which case titles are "to/from final/root state X",
            or an array of titles, one per lineage.
        kwargs
            Keyword arguments for :func:`scvelo.pl.scatter`.

        Returns
        -------
        None
            Nothing, just plots the metastable states.
        """

        attr = "_meta_lin_probs"
        error_msg = "Compute metastable states first as `.compute_metastable_states()`."

        if not discrete:
            self._plot_probabilities(
                attr=attr,
                error_msg=error_msg,
                lineages=lineages,
                cluster_key=cluster_key,
                mode=mode,
                time_key=time_key,
                show_dp=False,
                title=title,
                same_plot=same_plot,
                color_map=cmap,
                **kwargs,
            )
        else:
            self._plot_states(
                attr=attr,
                error_msg=error_msg,
                same_plot=same_plot,
                title=title,
                **kwargs,
            )

    def plot_main_states(
        self,
        discrete: bool = False,
        lineages: Optional[Union[str, Iterable[str]]] = None,
        cluster_key: Optional[str] = None,
        mode: str = "embedding",
        time_key: str = "latent_time",
        same_plot: bool = True,
        show_dp: bool = False,
        title: Optional[str] = None,
        cmap: Union[str, mpl.colors.ListedColormap] = cm.viridis,
        **kwargs,
    ) -> None:
        """
        Plot the absorption probabilities in the given embedding.

        Params
        ------
        discrete
            Whether to plot the top cells from each linages or the probabilities.
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
        same_plot
            Whether to plot the lineages on the same plot using color gradients when :paramref:`method='embedding'`.
        show_dp
            Whether to show :paramref:`diff_potential` when :paramref:`n_cells` `=None`.
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
            Nothing, just plots the main states.
        """

        attr = "_lin_probs"
        error_msg = (
            "Compute main states first as `.compute_main_states()` "
            "or set them manually as `.set_main_states()`."
        )

        if not discrete:
            self._plot_probabilities(
                attr=attr,
                error_msg=error_msg,
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
        else:
            self._plot_states(
                attr=attr,
                error_msg=error_msg,
                same_plot=same_plot,
                title=title,
                **kwargs,
            )

    def _set_main_states(self, n_cells: int, write_to_adata: bool = True) -> None:
        probs = self._lin_probs[[n for n in self._lin_probs.names if n != "rest"]]
        a_discrete, _ = _fuzzy_to_discrete(
            a_fuzzy=probs,
            n_most_likely=n_cells,
            remove_overlap=REMOVE_OVERLAP,
            raise_threshold=0.2,
            check_row_sums=False,
        )
        self._main_states = _series_from_one_hot_matrix(
            a=a_discrete, index=self.adata.obs_names, names=probs.names
        )

        if write_to_adata:
            self.adata.obs[self._rc_key] = self._main_states
            self.adata.uns[_colors(self._rc_key)] = probs[
                list(self._main_states.cat.categories)
            ].colors

    def set_main_states(
        self,
        names: Optional[Union[Iterable[str], str]] = None,
        n_cells: int = 30,
        redistribute: bool = True,
        **kwargs,
    ):
        """
        Manually select the main states from the metastable states.

        Params
        ------
        names
            Names of the main states. Multiple states can be combined using `','`, such as `['Alpha, Beta', 'Epsilon']`.
        n_cells
            Number of most likely cells from each main state to select. If `None`, on main states will be selected,
            only the lineage probabilities may be redistributed.
        redistribute
            Whether to redistribute the probability mass of unselected lineages or create a `'rest'` lineage.
        kwargs
            Keyword arguments for :meth:`cellrank.tl.Lineage.reduce` when redistributing the probability mass.

        Returns
        -------
        None
            Nothing, but updates the following fields:

                - :paramref:`lineage_probabilities`
                - :paramref:`diff_potential`
        """
        if not isinstance(n_cells, int):
            raise TypeError(
                f"Expected `n_cells` to be of type `int`, found `{type(n_cells).__name__}`."
            )

        if n_cells <= 0:
            raise ValueError(f"Expected `n_cells` to be positive, found `{n_cells}`.")

        if names is None:
            names = self._meta_lin_probs.names
            redistribute = False

        if isinstance(names, str):
            names = [names]

        if len(names) == 1 and redistribute:
            logg.warning(
                "Redistributing the mass to only 1 state will create a constant vector of ones. Not redistributing"
            )
            redistribute = False

        names = list(names)
        kwargs["return_weights"] = False

        if redistribute:
            self._lin_probs = self._meta_lin_probs[names + [Lin.OTHERS]]
            self._lin_probs = self._lin_probs.reduce(
                [" or ".join(_convert_lineage_name(name)) for name in names], **kwargs
            )
        else:
            self._lin_probs = self._meta_lin_probs[names + [Lin.REST]]

        self._set_main_states(n_cells)
        self._dp = entropy(self._lin_probs.X.T)

        # compute the aggregated probability of being a root/final state (no matter which)
        scaled_probs = self._lin_probs[
            [n for n in self._lin_probs.names if n != "rest"]
        ].X
        scaled_probs /= scaled_probs.max(0)
        self._main_states_probabilities = scaled_probs.max(1)

        # write to adata
        self.adata.obs[_dp(self._lin_key)] = self._dp
        self.adata.obs[_probs(self._rc_key)] = self._main_states_probabilities

        self.adata.uns[_lin_names(self._lin_key)] = self._lin_probs.names
        self.adata.uns[_colors(self._lin_key)] = self._lin_probs.colors

        self.adata.obsm[self._lin_key] = self._lin_probs

        logg.info("Adding `.lineage_probabilities\n       `.diff_potential`")

    def compute_main_states(
        self,
        method: str = "eigengap",
        redistribute: bool = True,
        alpha: Optional[float] = 1,
        min_self_prob: Optional[float] = None,
        n_main_states: Optional[int] = None,
        n_cells: int = 30,
        **kwargs,
    ):
        """
        Automatically select the main states from metastable states.

        Params
        ------
        method
            One of the following:

            - `'eigengap'` - select the number of states based on the eigengap of the transition matrix
            - `'eigengap_coarse'`- select the number of states based on the eigengap of the diagonal
                of the coarse-grained transition matrix
            - `'min_self_prob'`- select states which have the given minimum probability on the diagonal
                of the coarse-grained transition matrix
            - `'top_n'`- select top :paramref:`n_main_states` based on the probability on the diagonal
                of the coarse-grained transition matrix
        redistribute
            Whether to redistribute the probability mass of unselected lineages or create a `'rest'` lineage.
        alpha
            Used when :paramref:`method` `='eigengap'` or `='eigengap_coarse`.
        min_self_prob
            Used when :paramref:`method` `='min_self_prob'`.
        n_main_states
            Used when :paramref:`method` `='top_n'`.
        n_cells
            Number of most likely cells from each main state to select.
        kwargs
            Keyword arguments for :meth:`cellrank.tl.Lineage.reduce` when redistributing the mass.

        Returns
        -------
        None
            Nothing, just updates the following fields:

                - :paramref:`lineage_probabilities`
                - :paramref:`diff_potential`
        """
        if len(self.metastable_states.cat.categories) == 1:
            logg.warning(
                "Found only one metastable state. Making it the single main state. "
            )
            self.set_main_states(None, redistribute=False, n_cells=n_cells, **kwargs)
            return

        if method == "eigengap":
            if self.eigendecomposition is None:
                raise RuntimeError(
                    "Compute eigendecomposition first as `.compute_eig()`."
                )
            n_main_states = _eigengap(self.eigendecomposition["D"], alpha=alpha) + 1
        elif method == "eigengap_coarse":
            if self._coarse_T is None:
                raise RuntimeError(
                    "Compute metastable states first as `.compute_metastable_states()`."
                )
            n_main_states = _eigengap(
                np.sort(np.diag(self._coarse_T)[::-1]), alpha=alpha
            )
        elif method == "top_n":
            if n_main_states is None:
                raise ValueError(
                    "Argument `n_main_states` must not be `None` for `method='top_n'`."
                )
            elif n_main_states <= 0:
                raise ValueError(
                    f"Expected `n_main_states` to be positive, found `{n_main_states}`."
                )
        elif method == "min_self_prob":
            if min_self_prob is None:
                raise ValueError(
                    "Argument `min_self_prob` must not be `None` for `method='min_self_prob'`."
                )
            self_probs = pd.Series(
                np.diag(self._coarse_T), index=self._coarse_T.columns
            )
            names = self_probs[self_probs.values >= min_self_prob].index
            self.set_main_states(
                names, redistribute=redistribute, n_cells=n_cells, **kwargs
            )
            return
        else:
            raise ValueError(
                f"Invalid method `{method!r}`. Valid options are `'eigengap', 'eigengap_coarse', "
                f"'top_n' and 'min_self_prob'`."
            )

        names = self._coarse_T.columns[np.argsort(np.diag(self._coarse_T))][
            -n_main_states:
        ]
        self.set_main_states(
            names, redistribute=redistribute, n_cells=n_cells, **kwargs
        )

    def _assign_metastable_states(
        self,
        memberships: np.ndarray,
        n_cells: Optional[int] = 30,
        cluster_key: str = "clusters",
        en_cutoff: Optional[float] = 0.7,
        p_thresh: float = 1e-15,
        check_row_sums: bool = True,
    ) -> None:
        """
        Map a fuzzy clustering to pre-computed annotations to get names and colors.

        Given the fuzzy clustering we have computed, we would like to select the most likely cells from each state
        and use these to give each state a name and a color by comparing with pre-computed, categorical cluster
        annotations.

        Params
        --------
        memberships
            Fuzzy clustering
        n_cells
            Number of cells to be used to represent each state.
        cluster_key
            Key from `adata.obs` to get reference cluster annotations.
        en_cutoff
            Threshold to decide when we we want to warn the user about an uncertain name mapping. This happens when
            one fuzzy state overlaps with several reference clusters, and the most likely cells are distributed almost
            evenly across the reference clusters.
        p_thresh
            Only used to detect cell cycle stages. These have to be present in `adata.obs` as `G2M_score` and `S_score`.
        check_row_sums
            Check whether rows in `memberships` sum to 1.

        Returns
        --------
        None
            Writes a lineage object which mapped names and colors. Also creates a categorical :class:`pandas.Series`
            `.metastable_states`, where the top `n_cells` cells represent each fuzzy state.
        """

        if n_cells is None:
            max_assignment = np.argmax(memberships, axis=1)
            _meta_assignment = pd.Series(
                index=self.adata.obs_names, data=max_assignment, dtype="category"
            )
            # sometimes, the assignment can have a missing category and the Lineage creation therefore fails
            # keep it as ints when `n_cells != None`
            _meta_assignment.cat.set_categories(
                list(range(memberships.shape[1])), inplace=True
            )

            logg.debug(
                "DEBUG: Setting the metastable states using metastable assignment"
            )
            metastable_states = _meta_assignment.astype(str).astype("category").copy()
            not_enough_cells = []
        else:
            if n_cells <= 0:
                raise ValueError(
                    f"Expected `n_cells` to be positive, found `{n_cells}`."
                )

            logg.debug(
                "DEBUG: Setting the metastable states using metastable memberships"
            )

            # select the most likely cells from each metastable state
            a_discrete, not_enough_cells = _fuzzy_to_discrete(
                a_fuzzy=memberships,
                n_most_likely=n_cells,
                remove_overlap=REMOVE_OVERLAP,
                raise_threshold=0.2,
                check_row_sums=check_row_sums,
            )
            metastable_states = _series_from_one_hot_matrix(
                a=a_discrete, index=self.adata.obs_names
            )
            not_enough_cells = not_enough_cells.astype("str")

        # _set_categorical_labels creates the names, we still need to remap the group names
        orig_cats = metastable_states.cat.categories
        self._set_categorical_labels(
            attr_key="_meta_states",
            pretty_attr_key="metastable_states",
            cat_key=self._meta_key,
            add_to_existing_error_msg="Compute metastable states first as `.compute_metastable_states()`.",
            categories=metastable_states,
            cluster_key=cluster_key,
            en_cutoff=en_cutoff,
            p_thresh=p_thresh,
            add_to_existing=False,
        )
        name_mapper = dict(zip(orig_cats, self.metastable_states.cat.categories))
        _print_insufficient_number_of_cells(
            [name_mapper.get(n, n) for n in not_enough_cells], n_cells
        )

        logg.debug(
            "DEBUG: Setting metastable lineage probabilities based on GPCCA membership vectors"
        )
        self._meta_lin_probs = Lineage(
            memberships,
            names=list(metastable_states.cat.categories),
            colors=self._meta_states_colors,
        )

    def _plot_states(
        self,
        attr: str,
        error_msg: str,
        same_plot: bool = True,
        title: Optional[Union[str, List[str]]] = None,
        **kwargs,
    ):
        """
        Plot the main states for each uncovered lineage.

        Params
        ------
        n_cells
            Number of most likely cells per lineage.
        same_plot
            Whether to plot the lineages on the same plot or separately.
        title
            The title of the plot.
        kwargs
            Keyword arguments for :func:`scvelo.pl.scatter`.

        Returns
        -------
        None
            Nothing, just plots the metastable or main states.
        """

        def cleanup():
            for key in to_clean:
                try:
                    del self.adata.obs[key]
                except KeyError:
                    pass
                try:
                    del self.adata.uns[f"{key}_colors"]
                except KeyError:
                    pass

        probs = getattr(self, attr, None)
        if probs is None:
            raise RuntimeError(error_msg)

        probs = probs[[n for n in probs.names if n != "rest"]]

        if attr == "_meta_lin_probs":
            _main_states = self._meta_states
        elif attr == "_lin_probs":
            _main_states = self._main_states
        else:
            raise RuntimeError(f"Invalid attribute name: `{attr!r}`.")

        to_clean = []
        try:
            if same_plot:
                key = _generate_random_keys(self.adata, "obs")[0]
                to_clean = [key]
                self.adata.obs[key] = _main_states
                self.adata.uns[f"{key}_colors"] = probs.colors

                if title is None:
                    title = (
                        "metastable states (backward)"
                        if self.kernel.backward
                        else "metastable states (forward)"
                    )
                scv.pl.scatter(self.adata, title=title, color=key, **kwargs)
            else:
                keys = _generate_random_keys(
                    self.adata, "obs", len(_main_states.cat.categories)
                )

                to_clean = keys

                for key, cat in zip(keys, _main_states.cat.categories):
                    d = _main_states.copy()
                    d[_main_states != cat] = None
                    d.cat.set_categories([cat], inplace=True)

                    self.adata.obs[key] = d
                    self.adata.uns[f"{key}_colors"] = probs[cat].colors

                scv.pl.scatter(
                    self.adata,
                    color=keys,
                    title=list(_main_states.cat.categories) if title is None else title,
                    **kwargs,
                )
        except Exception as e:
            raise e
        finally:
            cleanup()

    def plot_coarse_T(
        self,
        show_stationary_dist: bool = True,
        show_initial_dist: bool = False,
        cmap: Union[str, mcolors.ListedColormap] = "viridis",
        xtick_rotation: float = 45,
        annotate: bool = True,
        show_cbar: bool = True,
        title: Optional[str] = None,
        figsize: Tuple[float, float] = (8, 8),
        dpi: float = 80,
        save: Optional[Union[Path, str]] = None,
        text_kwargs: Mapping[str, Any] = MappingProxyType({}),
        **kwargs,
    ) -> None:
        """
        Plot the coarse-grained transition matrix of the metastable states.

        Params
        ------
        show_stationary_dist
            Whether to show the stationary distribution, if present.
        show_initial_dist
            Whether to show the initial distribution.
        cmap
            Colormap to use.
        xtick_rotation
            Rotation of ticks on the x-axis.
        annotate
            Whether to display the text on each cell.
        show_cbar
            Whether to show colorbar.
        title
            Title of the figure.
        figsize
            Size of the figure.
        dpi
            Dots per inch.
        save
            Filename where to save the plots.
            If `None`, just show the plots.
        text_kwargs
            Keyword arguments for :func:`matplotlib.pyplot.text`.
        kwargs
            Keyword arguments for :func:`matplotlib.pyplot.imshow`.

        Returns
        -------
        None
            Nothing, just plots and optionally saves the plot.
        """

        def stylize_dist(
            ax, data: np.ndarray, xticks_labels: Union[List[str], Tuple[str]] = ()
        ):
            _ = ax.imshow(data, aspect="auto", cmap=cmap, norm=norm)
            for spine in ax.spines.values():
                spine.set_visible(False)

            ax.set_xticklabels(xticks_labels)
            if xticks_labels:
                ax.set_xticks(np.arange(data.shape[1]))
                plt.setp(
                    ax.get_xticklabels(),
                    rotation=xtick_rotation,
                    ha="right",
                    rotation_mode="anchor",
                )
            else:
                ax.tick_params(
                    which="both", top=False, right=False, bottom=False, left=False
                )

            ax.set_yticks([])

        def annotate_heatmap(im, valfmt: str = "{x:.2f}"):
            # modified from matplotlib's site

            data = im.get_array()
            kw = {"ha": "center", "va": "center"}
            kw.update(**text_kwargs)

            # Get the formatter in case a string is supplied
            if isinstance(valfmt, str):
                valfmt = mpl.ticker.StrMethodFormatter(valfmt)

            # Loop over the data and create a `Text` for each "pixel".
            # Change the text's color depending on the data.
            texts = []
            for i in range(data.shape[0]):
                for j in range(data.shape[1]):
                    kw.update(color=_get_black_or_white(im.norm(data[i, j]), cmap))
                    text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
                    texts.append(text)

        def annotate_dist_ax(
            ax, data: np.ndarray, is_vertical: bool, valfmt: str = "{x:.2f}"
        ):
            if ax is None:
                return

            if isinstance(valfmt, str):
                valfmt = mpl.ticker.StrMethodFormatter(valfmt)

            kw = {"ha": "center", "va": "center"}
            kw.update(**text_kwargs)
            if is_vertical:
                kw["rotation"] = -90

            for i, val in enumerate(data):
                kw.update(color=_get_black_or_white(im.norm(val), cmap))
                ax.text(
                    0 if is_vertical else i,
                    i if is_vertical else 0,
                    valfmt(val, None),
                    **kw,
                )

        if self.coarse_T is None:
            raise RuntimeError(
                "Compute coarse transition matrix first as `.compute_metastable_states()` with `n_states` > 1."
            )

        if show_stationary_dist and self.coarse_stationary_distribution is None:
            logg.warning("Coarse stationary distribution is `None`, not plotting")
            show_stationary_dist = False
        if show_initial_dist and self._coarse_init_dist is None:
            logg.warning("Coarse initial distribution is `None`, not plotting")
            show_initial_dist = False

        hrs, wrs = [1], [1]
        if show_initial_dist:
            wrs += [0.05]
        if show_stationary_dist:
            hrs += [0.05]
        if show_cbar:
            wrs += [
                0.025
            ] * 2  # dirty trick so that ylabel doesn't overlap with colorbar

        fig = plt.figure(constrained_layout=False, figsize=figsize, dpi=dpi)
        gs = plt.GridSpec(
            1 + show_stationary_dist,
            1 + show_initial_dist + (show_cbar * 2),
            height_ratios=hrs,
            width_ratios=wrs,
            wspace=0.10,
            hspace=0.05,
        )
        if isinstance(cmap, str):
            cmap = plt.get_cmap(cmap)

        ax = fig.add_subplot(gs[0, 0])
        cax = fig.add_subplot(gs[:1, -1])
        init_ax, stat_ax = None, None

        labels = list(self._coarse_T.columns)

        tmp = self.coarse_T
        if show_initial_dist:
            tmp = np.c_[tmp, self.coarse_stationary_distribution]
        if show_initial_dist:
            tmp = np.c_[tmp, self._coarse_init_dist]
        norm = mpl.colors.Normalize(vmin=np.nanmin(tmp), vmax=np.nanmax(tmp))

        if show_stationary_dist:
            stat_ax = fig.add_subplot(gs[1, 0])
            stylize_dist(
                stat_ax,
                np.array(self.coarse_stationary_distribution).reshape(1, -1),
                xticks_labels=labels,
            )
            stat_ax.set_xlabel("stationary distribution")  # , ha="right", x=1)

        if show_initial_dist:
            init_ax = fig.add_subplot(gs[0, 1])
            stylize_dist(init_ax, np.array(self._coarse_init_dist).reshape(-1, 1))

            init_ax.yaxis.set_label_position("right")
            init_ax.set_ylabel("initial distribution", rotation=-90, va="bottom")

        im = ax.imshow(self.coarse_T, aspect="auto", cmap=cmap, **kwargs)
        ax.set_title("coarse-grained transition matrix" if title is None else title)

        if show_cbar:
            _ = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)

        ax.set_yticks(np.arange(self.coarse_T.shape[0]))
        ax.set_yticklabels(labels)

        ax.tick_params(
            top=False,
            bottom=not show_stationary_dist,
            labeltop=False,
            labelbottom=not show_stationary_dist,
        )

        for spine in ax.spines.values():
            spine.set_visible(False)

        if not show_stationary_dist:
            ax.set_xticks(np.arange(self.coarse_T.shape[1]))
            ax.set_xticklabels(labels)
            plt.setp(
                ax.get_xticklabels(),
                rotation=xtick_rotation,
                ha="right",
                rotation_mode="anchor",
            )
        else:
            ax.set_xticks([])

        ax.set_yticks(np.arange(self.coarse_T.shape[0] + 1) - 0.5, minor=True)
        ax.tick_params(
            which="minor", bottom=not show_stationary_dist, left=False, top=False
        )

        if annotate:
            annotate_heatmap(im)
            annotate_dist_ax(
                stat_ax, self.coarse_stationary_distribution.values, is_vertical=False
            )
            annotate_dist_ax(init_ax, self._coarse_init_dist, is_vertical=True)

        if save:
            save_fig(fig, save)

        fig.show()

    def compute_gdpt(
        self, n_components: int = 10, key_added: str = "gdpt_pseudotime", **kwargs
    ):
        """
        Compute generalized DPT making use of the real Schur decomposition.

        Params
        ------
        n_components
            Number of real Schur vectors to consider.
        key_added
            Key in :paramref:`adata` `.obs` where to save the pseudotime.
        kwargs
            Keyword arguments for :math:`compute_schur` if Schur decomposition if not found.

        Returns
        -------
        None
            Nothing, just updates :paramref:`adata` .obs[:paramref:`key_added`] with the computed pseudotime.
        """

        def _get_dpt_row(e_vals: np.ndarray, e_vecs: np.ndarray, i: int):
            row = sum(
                (
                    np.abs(e_vals[eval_ix])
                    / (1 - np.abs(e_vals[eval_ix]))
                    * (e_vecs[i, eval_ix] - e_vecs[:, eval_ix])
                )
                ** 2
                # account for float32 precision
                for eval_ix in range(0, e_vals.size)
                if np.abs(e_vals[eval_ix]) < 0.9994
            )

            return np.sqrt(row)

        if "iroot" not in self.adata.uns.keys():
            raise KeyError("Key `'iroot'` not found in `adata.uns`.")

        iroot = self.adata.uns["iroot"]
        if isinstance(iroot, str):
            iroot = np.where(self.adata.obs_names == iroot)[0]
            if not len(iroot):
                raise ValueError(
                    f"Unable to find cell with name `{self.adata.uns['iroot']!r}` in `adata.obs_names`."
                )
            iroot = iroot[0]

        if n_components < 2:
            raise ValueError(
                f"Expected number of components >= 2, found `{n_components}`."
            )

        if self._schur_vectors is None:
            logg.warning("No Schur decomposition found. Computing")
            self.compute_schur(n_components, **kwargs)
        elif self._schur_matrix.shape[1] < n_components:
            logg.warning(
                f"Requested `{n_components}` components, but only `{self._schur_matrix.shape[1]}` were found. "
                f"Recomputing using default values"
            )
            self.compute_schur(n_components)
        else:
            logg.debug("DEBUG: Using cached Schur decomposition")

        start = logg.info(
            f"Computing Generalized Diffusion Pseudotime using n_components = {n_components}"
        )

        Q, eigenvalues = (
            self._schur_vectors,
            self.eigendecomposition["D"],
        )
        # may have to remove some values if too many converged
        Q, eigenvalues = Q[:, :n_components], eigenvalues[:n_components]

        D = _get_dpt_row(eigenvalues, Q, i=iroot)
        pseudotime = D / np.max(D[np.isfinite(D)])
        self.adata.obs[key_added] = pseudotime

        logg.info(f"Adding `{key_added!r}` to `adata.obs`\n    Finish", time=start)

    def copy(self) -> "GPCCA":
        """
        Return a copy of itself.

        Returns
        -------
        :class:`cellrank.tl.GPCCA`
            A copy of itself.
        """

        kernel = copy(self.kernel)  # doesn't copy the adata object
        g = GPCCA(kernel, self.adata.copy(), inplace=False, read_from_adata=False)

        g._eig = deepcopy(self.eigendecomposition)

        g._lin_probs = copy(self.lineage_probabilities)
        g._dp = copy(self.diff_potential)

        g._gpcca = deepcopy(self._gpcca)

        g._schur_vectors = copy(self.schur_vectors)
        g._schur_matrix = copy(self._schur_matrix)
        g._coarse_T = copy(self.coarse_T)

        g._meta_states = copy(self._meta_states)
        g._meta_states_colors = copy(self._meta_states_colors)
        g._meta_lin_probs = copy(self._meta_lin_probs)

        g._main_states = copy(self.main_states)
        g._main_states_probabilities = copy(self._main_states_probabilities)

        g._coarse_stat_dist = copy(self.coarse_stationary_distribution)
        g._coarse_init_dist = copy(self._coarse_init_dist)

        g._G2M_score = copy(self._G2M_score)
        g._S_score = copy(self._S_score)

        g._g2m_key = self._g2m_key
        g._s_key = self._s_key
        g._key_added = self._key_added

        g._is_irreducible = self.irreducible
        g._rec_classes = copy(self._rec_classes)
        g._trans_classes = copy(self._trans_classes)

        return g

    @property
    def schur_vectors(self) -> np.ndarray:
        """Schur vectors."""
        return self._schur_vectors

    @property
    def coarse_T(self) -> pd.DataFrame:
        """Coarse-grained transition matrix between metastable states."""
        return self._coarse_T

    @property
    def metastable_states(self) -> pd.Series:
        """Metastable states."""
        return self._meta_states

    @property
    def coarse_stationary_distribution(self) -> pd.Series:
        """Coarse-grained stationary distribution of metastable states."""
        return self._coarse_stat_dist

    @property
    def main_states(self) -> pd.Series:
        """Subset and/or combination of main states."""
        return self._main_states

    @property
    def main_states_probabilities(self) -> pd.Series:
        """Upper bound of of becoming a main states."""
        return self._main_states_probabilities


def _print_insufficient_number_of_cells(groups: Iterable[Any], n_cells: int):
    if groups:
        logg.debug(
            f"DEBUG: The following groups have less than requested number of cells ({n_cells}): "
            f"`{', '.join(sorted(map(str, groups)))}`"
        )
