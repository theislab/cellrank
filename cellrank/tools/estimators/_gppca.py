# -*- coding: utf-8 -*-
from typing import Optional, List, Tuple, Dict, Union, Mapping, Any, Iterable
from types import MappingProxyType
from anndata import AnnData
from msmtools.analysis.dense.gpcca import GPCCA as _GPPCA
from scanpy import logging as logg
from scipy.stats import entropy
from copy import copy, deepcopy
from pathlib import Path
from mpl_toolkits.axes_grid1 import make_axes_locatable

from cellrank.tools._lineage import Lineage
from cellrank.tools._constants import Lin, MetaKey, _colors, _lin_names, _dp, _probs
from cellrank.tools.estimators._base_estimator import BaseEstimator
from cellrank.tools._utils import (
    _eigengap,
    _get_black_or_white,
    _convert_lineage_name,
    generate_random_keys,
    save_fig,
)
from cellrank.tools.kernels._kernel import KernelExpression

import os
import numpy as np
import pandas as pd
import scvelo as scv
import seaborn as sns
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


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

    def compute_eig(self, k: int = 20, which: str = "LR", alpha: float = 1) -> None:
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

        Returns
        -------
        None
            Nothing, but updates the following fields: paramref:`eigendecomposition`.
        """
        self._compute_eig(k=k, which=which, alpha=alpha, only_evals=True)

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
                "Compute Schur vectors as `.compute_metastable_states()` first."
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
        upper_triangular_only: bool = True,
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
        upper_triangular_only
            Whether to show only the upper triangular matrix, including diagonal. Schur matrix should be an
            upper triangular matrix, but there can be small numerical imprecision.
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
                "Compute Schur matrix as `.compute_metastable_states()` first."
            )

        fig, ax = plt.subplots(
            figsize=self._schur_matrix.shape if figsize is None else figsize, dpi=dpi
        )

        divider = make_axes_locatable(ax)  # square=True make the colorbar a bit bigger
        cbar_ax = divider.append_axes("right", size="2.5%", pad=0.05)

        if upper_triangular_only:
            mask = np.zeros_like(self._schur_matrix, dtype=np.bool)
            mask[np.tril_indices_from(mask, k=-1)] = True
            if not np.allclose(self._schur_matrix[mask], 0.0):
                logg.warning(
                    "Elements of the lower triangular matrix are not close to `0`. "
                    "Consider visualizing the deviance `upper_triangular_only=False`."
                )

            vmin, vmax = (
                np.min(self._schur_matrix[~mask]),
                np.max(self._schur_matrix[~mask]),
            )
        else:
            mask = None
            vmin, vmax = np.min(self._schur_matrix), np.max(self._schur_matrix)

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

    def compute_metastable_states(
        self,
        n_states: Union[int, Tuple[int, int], List[int], Dict[str, int]],
        initial_distribution: Optional[np.ndarray] = None,
        use_min_chi: bool = False,
        method: str = "krylov",
        which: str = "LM",
        n_cells: Optional[int] = 30,
        cluster_key: str = None,
        en_cutoff: Optional[float] = 0.7,
        p_thresh: float = 1e-15,
    ):
        """
        Calculate the metastable states.

        Params
        ------
        n_states
            Number of metastable states.
        initial_distribution
            Input probability distribution over all cells. If `None`, uniform is chosen.
        use_min_chi
            Whether to use :meth:`msmtools.analysis.dense.gpcca.GPCCA.minChi` to calculate the number of metastable states.
            If `True`, :paramref:`n_states` corresponds to an interval `[min, max]` inside of which
            the potentially optimal number of metastable states is searched.
        method
            Method for calculating the Schur vectors. Valid options are: `'krylov'`, `'brandts'` and `'scipy'`.
            For benefits of each method, see :class:`msmtoots.analysis.dense.gpcca.GPCCA`.
        which
            Eigenvalues are in general complex. `'LR'` - largest real part, `'LM'` - largest magnitude.
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

        gpcca = _GPPCA(self._T, eta=initial_distribution, z=which, method=method)

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

            logg.debug(f"DEBUG: Calculating minChi within interval [{minn}, {maxx}]")
            n_states = int(np.arange(minn, maxx)[np.argmax(gpcca.minChi(minn, maxx))])
        elif not isinstance(n_states, int):
            raise ValueError(
                f"Expected `n_states` to be integer when `use_min_chi=False`, found `{type(n_states).__name__}`."
            )

        start = logg.info("Computing metastable states")

        gpcca = gpcca.optimize(m=n_states)

        # when `n_cells!=None` and the overlap is high, we're skipping some metastable states
        valid_ixs = self._assign_metastable_states(
            gpcca.memberships,
            gpcca.metastable_assignment,
            n_cells,
            cluster_key=cluster_key,
            p_thresh=p_thresh,
            en_cutoff=en_cutoff,
        )
        logg.debug(
            f"Selected `{len(valid_ixs)}` out of `{n_states}` due to an overlapp caused by `n_cells={n_cells}`"
        )

        self._lin_probs = None
        self._schur_vectors = gpcca.schur_vectors
        self._schur_matrix = gpcca.R  # gpcca.schur_matrix

        names = self._meta_lin_probs.names
        self._coarse_T = pd.DataFrame(
            gpcca.coarse_grained_transition_matrix[valid_ixs, :][:, valid_ixs],
            index=names,
            columns=names,
        )
        self._coarse_init_dist = pd.Series(
            gpcca.coarse_grained_input_distribution[valid_ixs], index=names
        )
        self._coarse_stat_dist = pd.Series(
            gpcca.coarse_grained_stationary_probability[valid_ixs], index=names
        )

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
        max_avail_cells = n_cells * probs.shape[1]
        if max_avail_cells > self.adata.n_obs:
            raise ValueError(
                f"Total number of requested cells ({max_avail_cells}) "
                f"exceeds the total number of cells ({self.adata.n_obs})."
            )

        self._n_cells = n_cells
        self._main_states, _, _ = self._select_cells(n_cells, memberships=probs)

        if write_to_adata:
            self.adata.obs[self._rc_key] = self._main_states
            self.adata.uns[_colors(self._rc_key)] = probs[
                list(self._main_states.cat.categories)
            ].colors

    def set_main_states(
        self,
        names: Optional[Iterable[str]] = None,
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

        # write to adata
        self.adata.obs[_dp(self._lin_key)] = self._dp
        self.adata.obs[_probs(self._rc_key)] = self._lin_probs[
            [n for n in self._lin_probs.names if n != "rest"]
        ].X.max(axis=1)
        self.adata.obsm[self._lin_key] = self._lin_probs
        self.adata.uns[_lin_names(self._lin_key)] = self._lin_probs.names
        self.adata.uns[_colors(self._lin_key)] = self._lin_probs.colors

        logg.info("Adding `.lineage_probabilities\n       `.diff_potential`")

    def compute_main_states(
        self,
        method: str = "eigengap",
        redistribute: bool = False,
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
        self, n_cells: int, memberships: Union[np.ndarray, Lineage]
    ) -> Tuple[pd.Series, Union[np.ndarray, Lineage], List[int]]:
        # in this case, we need to be a bit careful because fuzzy clusters can largely overlap
        metastable_states = pd.Series(index=self.adata.obs_names, dtype="category")
        overlaps, cols, valid_ixs = {}, [], []

        if isinstance(memberships, Lineage):
            names = memberships.names
            # we always retain the same shape, causes problems
            memberships = memberships.X
        else:
            names = list(map(str, range(memberships.shape[1])))

        for ix, (name, col) in enumerate(zip(names, memberships.T)):
            p = np.argpartition(col, -n_cells)[-n_cells:]

            # handle the case of overlapping cells (fuzzy clustering)
            if len(metastable_states.cat.categories) > 0:
                current_labels = metastable_states.iloc[p]
                overlap = {
                    cl: np.sum(current_labels == cl)
                    for cl in current_labels.cat.categories
                    if np.sum(current_labels == cl) > 0
                }
                overlaps[name] = overlap
                if any(np.fromiter(overlap.values(), dtype=float) / n_cells > 0.8):
                    logg.warning(
                        "Found overlapping clusters with overlap > 80%. Skipping"
                    )
                    continue

            self._gppca_overlap = overlaps

            metastable_states.cat.add_categories(name, inplace=True)
            metastable_states.iloc[p] = name

            cols.append(col[:, None])
            valid_ixs.append(ix)

        # aggregate the non-overlapping columns together
        _memberships = np.concatenate(cols, axis=1)

        return metastable_states, _memberships, valid_ixs

    def _assign_metastable_states(
        self,
        memberships: np.ndarray,
        metastable_assignment: np.ndarray,
        n_cells: Optional[int],
        cluster_key: str,
        p_thresh,
        en_cutoff,
    ) -> List[int]:
        if n_cells is None:
            logg.debug(
                "DEBUG: Setting the metastable states using metastable assignment"
            )
            metastable_states = pd.Series(
                index=self._adata.obs_names,
                data=map(str, metastable_assignment),
                dtype="category",
            )
            # sometimes, the assignment can have a missing category and the Lineage creation therefore fails
            metastable_states.cat.set_categories(
                map(str, range(memberships.shape[1])), inplace=True
            )
            _memberships = memberships
            valid_ixs = (
                metastable_states.cat.categories,
                list(range(memberships.shape[1])),
            )
        else:
            logg.debug(
                "DEBUG: Setting the metastable states using metastable memberships"
            )
            metastable_states, _memberships, valid_ixs = self._select_cells(
                n_cells, memberships
            )

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

        logg.debug(
            "DEBUG: Setting metastable lineage probabilities based on GPCCA membership vectors"
        )
        self._meta_lin_probs = Lineage(
            _memberships,
            names=list(self._meta_states.cat.categories),
            colors=self._meta_states_colors,
        )

        return valid_ixs

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
            _main_states, *_ = self._select_cells(n_cells, probs)
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
                key = generate_random_keys(self.adata, "obs")[0]
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
                keys = generate_random_keys(
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
                f"Compute coarse transition matrix first as `.compute_metastable_states()`."
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
