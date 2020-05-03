# -*- coding: utf-8 -*-
from typing import Optional, List, Tuple, Dict, Union, Mapping, Any, Iterable
from types import MappingProxyType
from anndata import AnnData
from msmtools.analysis.dense.gpcca import GPCCA as _GPPCA
from scanpy import logging as logg
from scipy.stats import entropy
from copy import copy, deepcopy

from cellrank.tools._lineage import Lineage
from cellrank.tools._constants import Lin, RcKey, _colors, _lin_names
from cellrank.tools.estimators._base_estimator import BaseEstimator
from cellrank.tools._utils import save_fig, _eigengap, generate_random_keys
from cellrank.tools.kernels._kernel import KernelExpression

import os
import numpy as np
import pandas as pd
import scvelo as scv
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


class GPCCA(BaseEstimator):
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
            self._ms_key = str(RcKey.BACKWARD)
        else:
            self._ms_key = str(RcKey.FORWARD)

        self._schur_vectors = None
        self._coarse_T = None
        self._coarse_init_dist = None
        self._coarse_stat_dist = None
        self._gppca_overlap = None

        self._meta_states = None
        self._meta_states_colors = None
        self._meta_lin_probs = None

        self._main_states = None
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
            Used to compute the `eigengap`. gref:`alpha` is the weight given
            to the deviation of an eigenvalue from one.

        Returns
        -------
        None
            Nothing, but updates the following fields: gref:`eigendecomposition`.
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
        Params
        ------
        n_states:
        initial_distribution
        use_min_chi
        method
        which
            Eigenvalues are in general complex. `'LR'` - largest real part, `'LM'` - largest magnitude.
        cluster_key
        en_cutoff
        p_thresh

        Returns
        -------
        None
            Nothings, but updates the following fields:
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
                (n_states["n_min"], n_states["n_max"])
                if isinstance(n_states, dict)
                else n_states
            )
            if minn <= 1:
                raise ValueError(f"Minimum value must be > 1, found `{minn}`.")
            elif minn == 2:
                logg.warning(
                    "In most cases, 2 clusters will always be optimal. "
                    "If you really expect 2 clusters, use `n_clusters=2`.\nSetting minimum to 3"
                )
                minn = 3

            logg.debug(f"DEBUG: Calculating minChi within interval [{minn}, {maxx}]")
            n_states = np.arange(minn, maxx)[np.argmax(gpcca.minChi(minn, maxx))]

        start = logg.info("Computing metastable states")

        gpcca = gpcca.optimize(m=n_states)

        self._assign_metastable_states(
            gpcca.memberships,
            gpcca.metastable_assignment,
            n_cells,
            cluster_key=cluster_key,
            p_thresh=p_thresh,
            en_cutoff=en_cutoff,
        )
        self._lin_probs = None

        self._schur_vectors = gpcca.schur_vectors
        self._coarse_T = pd.DataFrame(
            gpcca.coarse_grained_transition_matrix,
            index=self._meta_lin_probs.names,
            columns=self._meta_lin_probs.names,
        )
        self._coarse_init_dist = pd.Series(
            gpcca.coarse_grained_input_distribution, index=self._meta_lin_probs.names
        )
        self._coarse_stat_dist = pd.Series(
            gpcca.coarse_grained_stationary_probability,
            index=self._meta_lin_probs.names,
        )

        logg.info(
            "Adding `.schur_vectors`\n"
            "       `.metastable_states`\n"
            "       `.coarse_transition_matrix`\n"
            "       `.coarse_stationary_distribution`\n"
            "    Finish",
            time=start,
        )

    def plot_meta_lin_probs(
        self,
        lineages: Optional[Union[str, Iterable[str]]] = None,
        cluster_key: Optional[str] = None,
        mode: str = "embedding",
        time_key: str = "latent_time",
        same_plot: bool = False,
        color_map: Union[str, mpl.colors.ListedColormap] = cm.viridis,
        **kwargs,
    ) -> None:
        """
        Plots the absorption probabilities of metastable states in the given embedding.

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
        same_plot
            Whether to plot the lineages on the same plot using color gradients when :paramref:`mode='embedding'`.
        color_map
            Colormap to use.
        kwargs
            Keyword arguments for :func:`scvelo.pl.scatter`.

        Returns
        -------
        None
            Nothing, just plots the absorption probabilities.
        """

        self._plot_probabilities(
            attr="_meta_lin_probs",
            error_msg="Compute metastable lineage probabilities first as `.compute_metastable_states()`.",
            lineages=lineages,
            cluster_key=cluster_key,
            mode=mode,
            time_key=time_key,
            same_plot=same_plot,
            color_map=color_map,
            **kwargs,
        )

    def set_main_states(self, names: Iterable[str], mode: str = "normalize"):
        """
        Params
        ------
        Returns
        -------
        Nothing, but updates the following fields:
        """
        names = list(names)
        if mode == "normalize":
            names += [Lin.NORM]
        elif mode == "rest":
            names += [Lin.REST]
        else:
            raise ValueError(
                f"Invalid mode `{mode!r}`. Valid options are `'normalize', 'rest'`."
            )

        self._n_cells = None  # invalidate cache
        self._lin_probs = self._meta_lin_probs[names]
        self._dp = entropy(self._lin_probs.X.T)

        # write to adata
        self.adata.obs[f"{self._lin_key}_dp"] = self._dp
        self.adata.uns[_lin_names(self._lin_key)] = self._lin_probs.names
        self.adata.uns[_colors(self._lin_key)] = self._lin_probs.colors

        logg.info(
            "Adding `.lineage_probabilities" "       `.diff_potential``" "    Finish"
        )

    def compute_main_states(
        self,
        method: str = "eigengap",
        mode: str = "normalize",
        alpha: Optional[float] = 1,
        min_self_prob: Optional[float] = None,
        n_main_states: Optional[int] = None,
    ):
        """
        """

        if method == "eigengap":
            if self.eigendecomposition is None:
                raise RuntimeError(
                    "Compute eigendecomposition first as `.compute_eig()`."
                )
            n_main_states = _eigengap(self.eigendecomposition["D"], alpha=alpha)
        elif method == "eigengap_coarse":
            if self._coarse_T is None:
                raise RuntimeError(
                    "Compute metastable states first as `.compute_metastable_states()`."
                )
            n_main_states = _eigengap(
                np.sort(np.diag(self._coarse_T)[::-1]), alpha=alpha
            )
        elif method == "top_k":
            if n_main_states is None:
                raise ValueError(
                    "Argument `n_main_states` must not be `None` for `method='top_k'`."
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
            self.set_main_states(names, mode)
            return
        else:
            raise ValueError(
                f"Invalid method `{method!r}`. Valid options are `'eigengap', 'eigengap_coarse', 'top_k' and 'min_self_prob'`."
            )

        names = self._coarse_T.columns[np.argsort(np.diag(self._coarse_T))][
            -n_main_states:
        ]
        self.set_main_states(names, mode)

    def _select_cells(
        self, n_cells: int, memberships: Union[np.ndarray, Lineage]
    ) -> Tuple[pd.Series, Union[np.ndarray, Lineage]]:
        # in this case, we need to be a bit careful because fuzzy clusters can largely overlap
        metastable_states = pd.Series(index=self.adata.obs_names, dtype="category")
        overlaps, cols = {}, []
        if isinstance(memberships, Lineage):
            names = memberships.names
            memberships = (
                memberships.X
            )  # we always retain the same shape, causes problems
        else:
            names = map(str, range(memberships.shape[1]))

        for name, col in zip(names, memberships.T):
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

        # aggregate the non-overlapping columns together
        _memberships = np.concatenate(cols, axis=1)

        return metastable_states, _memberships

    def _assign_metastable_states(
        self,
        memberships: np.ndarray,
        metastable_assignment: np.ndarray,
        n_cells: Optional[int],
        cluster_key: str,
        p_thresh,
        en_cutoff,
    ):
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
        else:
            logg.debug(
                "DEBUG: Setting the metastable states using metastable memberships"
            )
            metastable_states, _memberships = self._select_cells(n_cells, memberships)

        self._set_categorical_labels(
            attr_key="_meta_states",
            pretty_attr_key="metastable_states",
            cat_key=self._rc_key,
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

    def plot_main_states(self, n_cells: int, same_plot: bool = True, **kwargs):
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

        if self._lin_probs is None:
            raise RuntimeError(
                "Compute main states as `.compute_main_states()` or set them manually."
            )
        if n_cells <= 0:
            raise ValueError(f"Expected `n_cells` to be positive, found `{n_cells}`.")

        if n_cells != self._n_cells:
            self._main_states, _ = self._select_cells(n_cells, self._lin_probs)
            self._n_cells = n_cells
        else:
            logg.debug("DEBUG: Using cached main states")

        if n_cells * len(self._main_states.cat.categories) > self.adata.n_obs:
            raise ValueError(
                f"Total number of requested cells ({n_cells * len(self._main_states.cat.categories)}) "
                f"exceeds the total number of cells ({self.adata.n_obs})."
            )

        to_clean = []
        try:
            if same_plot:
                key = generate_random_keys(self.adata, "obs")
                self.adata.obs[key] = self._main_states
                self.adata.uns[f"{key}_colors"] = self._lin_probs.colors
                to_clean = [key]

                title = "from root cells" if self.kernel.backward else "to final cells"
                scv.pl.scatter(self.adata, title=title, color=key, **kwargs)
            else:
                prefix = "from" if self.kernel.backward else "to"
                titles = []
                keys = generate_random_keys(
                    self.adata, "obs", len(self._main_states.cat.categories)
                )
                to_clean = keys

                for key, cat in zip(keys, self._main_states.cat.categories):
                    d = self._main_states.copy()
                    d[self._main_states != cat] = None
                    d.cat.set_categories([cat], inplace=True)

                    titles.append(f"{prefix} {cat}")

                    self.adata.obs[key] = d
                    self.adata.uns[f"{key}_colors"] = self._lin_probs[cat].colors

                scv.pl.scatter(self.adata, color=keys, title=titles, **kwargs)
        except Exception as e:
            raise e
        finally:
            cleanup()

    def plot_coarse_T(
        self,
        show_stationary_dist: bool = True,
        show_initial_dist: bool = False,
        cmap: mcolors.ListedColormap = cm.viridis,
        xtick_rotation: float = 45,
        annotate: bool = True,
        show_cbar: bool = True,
        figsize: Tuple[float, float] = (8, 8),
        dpi: float = 80,
        save: Optional[Union[os.PathLike, str]] = None,
        text_kwargs: Mapping[str, Any] = MappingProxyType({}),
        **kwargs,
    ) -> None:
        """
        Plot the coarse-grained transition matrix between metastable states.

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
        figsize
            Size of the figure.
        dpi
            Dots per inch.
        save
            Filename where to save the plots.
            If `None`, just show the plots.
        text_kwargs
            Keyword arguments for `text`.
        kwargs
            Keyword arguments for `imshow`.

        Returns
        -------
        None
            Nothings just plots and optionally saves the plots.
        """

        def stylize_dist(
            ax, data: np.ndarray, xticks_labels: Union[List[str], Tuple[str]] = ()
        ):
            _ = ax.imshow(data, aspect="auto", cmap=cmap)
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

        def annotate_heatmap(
            im,
            valfmt: str = "{x:.2f}",
            textcolors: Tuple[str, str] = ("black", "white"),
            threshold: Optional[float] = None,
        ):
            # modified from matplotlib's site

            data = im.get_array()
            # Normalize the threshold to the images color range.
            threshold = (
                im.norm(np.median(data)) if threshold is None else im.norm(threshold)
            )
            kw = dict(horizontalalignment="center", verticalalignment="center")
            kw.update(**text_kwargs)

            # Get the formatter in case a string is supplied
            if isinstance(valfmt, str):
                valfmt = mpl.ticker.StrMethodFormatter(valfmt)

            # Loop over the data and create a `Text` for each "pixel".
            # Change the text's color depending on the data.
            texts = []
            for i in range(data.shape[0]):
                for j in range(data.shape[1]):
                    kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
                    text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
                    texts.append(text)

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

        ax = fig.add_subplot(gs[0, 0])
        cax = fig.add_subplot(gs[:1, -1])

        labels = list(self._coarse_T.columns)
        tmp = self.coarse_T

        if show_stationary_dist:
            stat_ax = fig.add_subplot(gs[1, 0])
            tmp = np.c_[tmp, self.coarse_stationary_distribution]

            stylize_dist(
                stat_ax,
                np.array(self.coarse_stationary_distribution).reshape(1, -1),
                xticks_labels=labels,
            )
            stat_ax.set_xlabel("Stationary Distribution")
        if show_initial_dist:
            init_ax = fig.add_subplot(gs[0, 1])
            tmp = np.c_[tmp, self._coarse_init_dist]

            stylize_dist(init_ax, np.array(self._coarse_init_dist).reshape(-1, 1))
            init_ax.yaxis.set_label_position("right")
            init_ax.set_ylabel("Initial Distribution", rotation=-90, va="bottom")

        im = ax.imshow(self.coarse_T, aspect="auto", cmap=cmap, **kwargs)
        ax.set_title("Coarse-grained Transition Matrix")

        if show_cbar:
            norm = mpl.colors.Normalize(vmin=np.nanmin(tmp), vmax=np.nanmax(tmp))
            _ = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)
            cax.set_ylabel("Value", rotation=-90, va="bottom")

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

        if save:
            save_fig(fig, save)

        fig.show()

    def copy(self) -> "GPCCA":
        """
        Return a copy of itself.
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
        return self._schur_vectors

    @property
    def coarse_T(self) -> pd.DataFrame:
        return self._coarse_T

    @property
    def metastable_states(self) -> pd.Series:
        return self._meta_states

    @property
    def coarse_stationary_distribution(self) -> pd.Series:
        return self._coarse_stat_dist

    @property
    def main_states(self) -> np.ndarray:
        return self._main_states
