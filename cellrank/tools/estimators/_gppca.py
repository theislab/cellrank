# -*- coding: utf-8 -*-
import os
from copy import copy, deepcopy
from types import MappingProxyType
from typing import Any, Dict, List, Tuple, Union, Mapping, Iterable, Optional
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import entropy
from pandas.api.types import is_categorical_dtype

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
    _convert_lineage_name,
    _generate_random_keys,
)
from cellrank.tools._colors import _get_black_or_white
from cellrank.tools._lineage import Lineage
from cellrank.tools._constants import Lin, MetaKey, _dp, _probs, _colors, _lin_names
from msmtools.analysis.dense.gpcca import GPCCA as _GPPCA
from cellrank.tools.kernels._kernel import KernelExpression
from cellrank.tools.estimators._base_estimator import BaseEstimator


class GPCCA(BaseEstimator):
    """
    Generalized Perron Cluster Cluster Analysis [GPCCA18]/.

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
        super().__init__(
            kernel,
            adata,
            inplace=inplace,
            read_from_adata=read_from_adata,
            g2m_key=g2m_key,
            s_key=s_key,
            key_added=key_added,
        )
        if kernel.backward:
            self._meta_key = str(MetaKey.BACKWARD)
        else:
            self._meta_key = str(MetaKey.FORWARD)

        self._gpcca: GPCCA = None
        self._schur_vectors = None
        self._coarse_T = None
        self._coarse_init_dist = None
        self._coarse_stat_dist = None
        self._gppca_overlap = None

        self._meta_states = None
        self._meta_states_colors = None
        self._meta_lin_probs = None

        self._main_states = None
        self._main_states_probabilities = None
        self._n_cells = None  # serves as a cache for plotting

    def compute_eig(
        self,
        k: int = 20,
        which: str = "LR",
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
                "Compute Schur matrix first as `.compute_schur()` or `.compute_metastable_states()` with `n_states` > 1."
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

        # make it available for plotting
        self._schur_vectors = self._gpcca.X
        self._schur_matrix = self._gpcca.R

    def compute_metastable_states(
        self,
        n_states: Union[int, Tuple[int, int], List[int], Dict[str, int]],
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
            Number of metastable states.
        use_min_chi
            Whether to use :meth:`msmtools.analysis.dense.gpcca.GPCCA.minChi` to calculate the number of metastable states.
            If `True`, :paramref:`n_states` corresponds to an interval `[min, max]` inside of which
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
            Nothings, but updates the following fields:

                - :paramref:`schur_vectors`
                - :paramref:`coarse_T`
                - :paramref:`coarse_stationary_distribution`
        """

        if n_states == 1:
            if self.eigendecomposition is None:
                raise RuntimeError(
                    "Compute eigendecomposition first as `.compute_eig()`."
                )

            start = logg.info("Computing metastable states")
            logg.warning("For `n_states=1`, stationary distribution is computed")

            k = self.eigendecomposition["params"]["k"]
            which = self.eigendecomposition["params"]["which"]
            alpha = self.eigendecomposition["params"]["alpha"]

            self._compute_eig(k=k, which=which, alpha=alpha, only_evals=False)
            stationary_dist = self.eigendecomposition["stationary_dist"]

            self._assign_metastable_states(
                stationary_dist[:, None],
                np.zeros_like(stationary_dist),
                n_cells,
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
        else:
            if self._gpcca is None:
                raise RuntimeError(
                    "Compute Schur decomposition first as `.compute_schur()`."
                )

            if use_min_chi:
                if not isinstance(n_states, (dict, tuple, list)):
                    raise TypeError(
                        f"Expected `n_states` to be either `dict`, `tuple` or a `list`, found `{type(n_states).__name__}`."
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

                logg.debug(
                    f"DEBUG: Calculating minChi within interval [{minn}, {maxx}]"
                )
                n_states = int(
                    np.arange(minn, maxx)[np.argmax(self._gpcca.minChi(minn, maxx))]
                )
            elif not isinstance(n_states, int):
                raise ValueError(
                    f"Expected `n_states` to be an integer when `use_min_chi=False`, found `{type(n_states).__name__}`."
                )

            if self._gpcca.X.shape[1] < n_states:
                logg.warning(
                    f"Requested more metastable states ({n_states}) than available Schur vectors ({self._gpcca.X.shape[1]}). Recomputing the decomposition"
                )

            start = logg.info("Computing metastable states")

            self._gpcca = self._gpcca.optimize(m=n_states)

            self._assign_metastable_states(
                self._gpcca.memberships,
                self._gpcca.metastable_assignment,
                n_cells,
                cluster_key=cluster_key,
                p_thresh=p_thresh,
                en_cutoff=en_cutoff,
            )

            self._lin_probs = None
            self._schur_vectors = self._gpcca.schur_vectors
            self._schur_matrix = self._gpcca.R  # gpcca.schur_matrix

            names = self._meta_lin_probs.names
            self._coarse_T = pd.DataFrame(
                self._gpcca.coarse_grained_transition_matrix,
                index=names,
                columns=names,
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
        n_cells: Optional[int] = None,
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
        Plots the absorption probabilities of metastable states in the given embedding.

        Params
        ------
        n_cells
            Number of cells to sample from the distribution. If `None`, don't sample, justs plot the distribution.
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

        if n_cells is None:
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
                n_cells=n_cells,
                same_plot=same_plot,
                title=title,
                **kwargs,
            )

    def plot_main_states(
        self,
        n_cells: Optional[int] = None,
        lineages: Optional[Union[str, Iterable[str]]] = None,
        cluster_key: Optional[str] = None,
        mode: str = "embedding",
        time_key: str = "latent_time",
        same_plot: bool = False,
        show_dp: bool = False,
        title: Optional[str] = None,
        cmap: Union[str, mpl.colors.ListedColormap] = cm.viridis,
        **kwargs,
    ) -> None:
        """
        Plots the absorption probabilities in the given embedding.

        Params
        ------
        n_cells
            Number of cells to sample from the distribution. If `None`, don't sample, just plots the distribution.
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
        error_msg = "Compute main states first as `.compute_main_states()` or set them manually as `.set_main_states()`."

        if n_cells is None:
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
                n_cells=n_cells,
                same_plot=same_plot,
                title=title,
                **kwargs,
            )

    def _set_main_states(
        self, n_cells: Optional[int], write_to_adata: bool = True
    ) -> None:
        self._n_cells = n_cells

        if n_cells is None:
            return None

        probs = self._lin_probs[[n for n in self._lin_probs.names if n != "rest"]]
        self._n_cells = n_cells
        self._main_states = self._select_cells(n_cells, memberships=probs)

        if write_to_adata:
            self.adata.obs[self._rc_key] = self._main_states
            self.adata.uns[_colors(self._rc_key)] = probs[
                list(self._main_states.cat.categories)
            ].colors

    def set_main_states(
        self,
        names: Optional[Union[Iterable[str], str]] = None,
        n_cells: Optional[int] = 30,
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
            Number of most likely cells from each main state to select.
        redistribute
            Whether to redistribute the probability mass of unselected lineages or create a `'rest'` lineage.
        kwargs
            Keyword arguments for :meth:`cellrank.tl.Lineage.reduce` when redistributing the mass.

        Returns
        -------
        None
            Nothing, but updates the following fields:

                - :paramref:`lineage_probabilities`
                - :paramref:`diff_potential`
        """

        if names is None:
            names = self._meta_lin_probs.names
            redistribute = False

        if isinstance(names, str):
            names = [names]

        if len(names) == 1 and redistribute:
            logg.warning(
                "Redistributing the mass only to 1 state will create a constant vector of 1s. Not redistributing"
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

        # compute the aggregated probability of being a final state (no matter which)
        aggregated_state_probability = self._lin_probs[
            [n for n in self._lin_probs.names if n != "rest"]
        ].X.max(axis=1)
        aggregated_state_probability /= np.max(aggregated_state_probability)

        # write to adata
        self.adata.obs[_dp(self._lin_key)] = self._dp
        self.adata.obs[_probs(self._rc_key)] = aggregated_state_probability
        self.adata.obsm[self._lin_key] = self._lin_probs
        self.adata.uns[_lin_names(self._lin_key)] = self._lin_probs.names
        self.adata.uns[_colors(self._lin_key)] = self._lin_probs.colors

        logg.info("Adding `.lineage_probabilities\n       `.diff_potential`")

    def compute_main_states(
        self,
        method: str = "eigengap",
        redistribute: bool = True,
        alpha: Optional[float] = 1,
        min_self_prob: Optional[float] = None,
        n_main_states: Optional[int] = None,
        n_cells: Optional[int] = 30,
        **kwargs,
    ):
        """
        Automatically select the main states from metastable states.

        Params
        ------
        method
            One of the following:

            - `'eigengap'` - select the number of states based on the eigengap of the transition matrix
            - `'eigengap_coarse'`- select the number of states based on the eigengap of the diagonal of the coarse-grained transition matrix
            - `'min_self_prob'`- select states which have the given minimum probability on the diagonal of the coarse-grained transition matrix
            - `'top_n'`- select top :paramref:`n_main_states` based on the probability on the diagonal of the coarse-grained transition matrix
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
            Nothings, just updates the following fields:

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
                f"Invalid method `{method!r}`. Valid options are `'eigengap', 'eigengap_coarse', 'top_n' and 'min_self_prob'`."
            )

        names = self._coarse_T.columns[np.argsort(np.diag(self._coarse_T))][
            -n_main_states:
        ]
        self.set_main_states(
            names, redistribute=redistribute, n_cells=n_cells, **kwargs
        )

    def _select_cells(
        self, n_cells: int, memberships: Union[np.ndarray, Lineage],
    ) -> Tuple[pd.Series, Union[np.ndarray, Lineage], List[int]]:
        def set_categories(group):
            category = group["assignment"].iloc[0]

            if n_cells > len(group):
                logg.warning(
                    f"Number of requested cells ({n_cells}) exceeds the number of available cells ({len(group)}) for cluster `{category}`. Selecting all cells"
                )
                return pd.Series(group["assignment"])

            ixs = np.argpartition(group["values"], -n_cells)[-n_cells:]

            res = pd.Series([np.nan] * len(group), dtype="category", index=group.index)
            res.cat.add_categories(category, inplace=True)
            res.iloc[ixs.values] = category

            return res

        if isinstance(memberships, Lineage):
            # same is in msmtools, just update it
            self._meta_assignment = pd.Series(
                memberships.names[np.argmax(memberships.X, axis=1)],
                dtype="category",
                index=self.adata.obs_names,
            )

        membership_assignments = np.array(
            memberships[np.arange(memberships.shape[0]), self._meta_assignment.values]
        ).squeeze()  # squeeze because of Lineage

        meta_assignment = pd.DataFrame(
            {"assignment": self._meta_assignment, "values": membership_assignments}
        )
        meta_assignment["assignment"] = (
            meta_assignment["assignment"].astype(str).astype("category")
        )

        metastable_states = meta_assignment.groupby("assignment").apply(set_categories)
        # when there's only 1 state, we get back strange DataFrame
        metastable_states = (
            metastable_states.T.iloc[:, 0]
            if memberships.shape[1] == 1
            else metastable_states.droplevel(0)
        )
        metastable_states = metastable_states.sort_index().astype("category")
        metastable_states.index = self.adata.obs.index

        return metastable_states

    def _assign_metastable_states(
        self,
        memberships: np.ndarray,
        metastable_assignment: np.ndarray,
        n_cells: Optional[int],
        cluster_key: str,
        p_thresh,
        en_cutoff,
    ) -> None:
        # keep it as int for indexing when `n_cells!=None`
        self._meta_assignment = pd.Series(
            index=self.adata.obs_names, data=metastable_assignment, dtype="category",
        )
        # sometimes, the assignment can have a missing category and the Lineage creation therefore fails
        self._meta_assignment.cat.set_categories(
            range(memberships.shape[1]), inplace=True
        )

        if n_cells is None:
            logg.debug(
                "DEBUG: Setting the metastable states using metastable assignment"
            )
            metastable_states = (
                self._meta_assignment.astype(str).astype("category").copy()
            )
        else:
            logg.debug(
                "DEBUG: Setting the metastable states using metastable memberships"
            )
            metastable_states = self._select_cells(n_cells, memberships)

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

        # this makes indexing easier for Lineage class
        self._meta_assignment.cat.rename_categories(
            metastable_states.cat.categories, inplace=True
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
        n_cells: int,
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
            Nothings, just plots the main states.
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

        if n_cells <= 0:
            raise ValueError(f"Expected `n_cells` to be positive, found `{n_cells}`.")

        if attr == "_meta_lin_probs":  # plotting meta_lin_probs
            _main_states = self._select_cells(n_cells, memberships=probs)
        elif attr == "_lin_probs":
            if n_cells == self._n_cells:
                logg.debug("DEBUG: Using cached main states")
            else:
                self._set_main_states(n_cells, write_to_adata=False)
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
        save: Optional[Union[os.PathLike, str]] = None,
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
            Nothings just plots and optionally saves the plot.
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
            kw = dict(ha="center", va="center")
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

            kw = dict(ha="center", va="center")
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
                f"Compute coarse transition matrix first as `.compute_metastable_states()` with `n_states` > 1."
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
        ax.set_title("Coarse-grained transition matrix" if title is None else title)

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

    def copy(self) -> "GPCCA":
        """
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

        g._schur_vectors = copy(self.schur_vectors)
        g._coarse_T = copy(self.coarse_T)

        self._gppca_overlap = deepcopy(self._gppca_overlap)

        g._meta_states = copy(self._meta_states)
        g._meta_states_colors = copy(self._meta_states_colors)
        g._meta_lin_probs = copy(self._meta_lin_probs)

        g._main_states = copy(self.main_states)

        g._n_cells = self._n_cells

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
        """
        Returns
        -------
        :class:`pandas.DataFrame`
            Coarse-grained transition matrix between metastable states.
        """
        return self._coarse_T

    @property
    def metastable_states(self) -> pd.Series:
        """
        Returns
        -------
        :class:`pandas.Series`
            Metastable states
        """
        return self._meta_states

    @property
    def coarse_stationary_distribution(self) -> pd.Series:
        """
        Returns
        -------
        :class:`pandas.Series`
            Coarse-grained stationary distribution of metastable states
        """
        return self._coarse_stat_dist

    @property
    def main_states(self) -> pd.Series:
        """Subset and/or combination of main states."""
        return self._main_states

    @property
    def main_states_probabilities(self) -> pd.Series:
        """
        Returns
        -------
        :class:`pandas.Series`
            Upper bound of of becoming a main states.
        """
        return self._main_states_probabilities
