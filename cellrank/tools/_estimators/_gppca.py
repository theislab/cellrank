# -*- coding: utf-8 -*-
from typing import Optional, List, Tuple, Dict, Union, Mapping, Any
from types import MappingProxyType
from anndata import AnnData
from msmtools.analysis.dense.gpcca import GPCCA as _GPPCA
from scanpy import logging as logg
from scipy.stats import entropy

from cellrank.tools._lineage import Lineage
from cellrank.tools._constants import RcKey, _colors, _lin_names
from cellrank.tools._estimators._base_estimator import BaseEstimator
from cellrank.tools._utils import save_fig
from cellrank.tools.kernels._kernel import KernelExpression

import os
import numpy as np
import pandas as pd
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

        self._gpcca = None
        self._schur_vectors = None
        self._coarse_T = None
        self._coarse_init_dist = None
        self._coarse_stat_dist = None

        self._main_states = None
        self._main_states_colors = None

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
        n_cells: Optional[int] = None,
        cluster_key: str = "louvain",
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
            n_states = np.arange(minn, maxx)[np.argmax(self._gpcca.minChi(minn, maxx))]

        start = logg.info("Computing metastable states")

        self._gpcca: _GPPCA = _GPPCA(
            self._T, eta=initial_distribution, z=which, method=method
        ).optimize(m=n_states)

        self._schur_vectors = self._gpcca.schur_vectors
        self._coarse_T = self._gpcca.coarse_grained_transition_matrix
        self._coarse_init_dist = self._gpcca.coarse_grained_input_distribution
        self._coarse_stat_dist = self._gpcca.coarse_grained_stationary_probability

        self._assign_main_states(
            self._gpcca.memberships,
            n_cells,
            cluster_key=cluster_key,
            p_thresh=p_thresh,
            en_cutoff=en_cutoff,
        )

        logg.info(
            f"Adding `adata.obs[{self._rc_key!r}]`\n"
            f"       `.main_states`\n"
            f"       `.coarse_transition_matrix`\n"
            f"    Finish",
            time=start,
        )

    def _assign_main_states(
        self, memberships, n_cells: Optional[int], cluster_key: str, p_thresh, en_cutoff
    ):
        if n_cells is None:
            logg.debug("DEBUG: Setting the main states using metastable assignment")
            main_states = pd.Series(
                index=self._adata.obs_names,
                data=map(str, self._gpcca.metastable_assignment),
                dtype="category",
            )
        else:
            logg.debug("DEBUG: Setting the main states using metastable memberships")
            # in this case, we need to be a bit careful because fuzzy clusters can largely overlap
            main_states = pd.Series(index=self.adata.obs_names, dtype="category")
            overlaps, cols = {}, []

            for i, col in enumerate(memberships.T):
                p = np.flip(np.argsort(col))[:n_cells]

                # handle the case of overlapping cells (fuzzy clustering)
                i = str(i)
                if len(main_states.cat.categories) > 0:
                    current_labels = main_states.iloc[p]
                    overlap = {
                        cl: np.sum(current_labels == cl)
                        for cl in current_labels.cat.categories
                        if np.sum(current_labels == cl) > 0
                    }
                    overlaps[i] = overlap
                    if any(np.fromiter(overlap.values(), dtype=float) / n_cells > 0.8):
                        logg.warning(
                            "Found overlapping clusters with overlap > 80%. Skipping"
                        )
                        continue

                self.eigendecomposition["gpcca_overlap"] = overlaps

                main_states.cat.add_categories(i, inplace=True)
                main_states.iloc[p] = i
                cols.append(col[:, None])

            # aggregate the non-overlapping columns together
            _memberships = np.concatenate(cols, axis=1)

        self.set_main_states(
            main_states,
            cluster_key=cluster_key,
            en_cutoff=en_cutoff,
            p_thresh=p_thresh,
            add_to_existing=False,
        )

        logg.debug(
            "DEBUG: Setting lineage probabilities based on GPCCA membership vectors"
        )
        self._lin_probs = Lineage(
            memberships,
            names=list(self._main_states.cat.categories),
            colors=self._main_states_colors,
        )
        self._dp = entropy(self._lin_probs.T)

        self._adata.obsm[self._lin_key] = self._lin_probs
        self._adata.obs[f"{self._lin_key}_dp"] = self._dp
        self._adata.uns[_lin_names(self._lin_key)] = self._lin_probs.names
        self._adata.uns[_colors(self._lin_key)] = self._lin_probs.colors

    def set_main_states(
        self,
        states: Union[pd.Series, Dict[Any, Any]],
        cluster_key: Optional[str] = None,
        en_cutoff: Optional[float] = None,
        p_thresh: Optional[float] = None,
        add_to_existing: bool = False,
    ):
        self._set_categorical_labels(
            attr_key="_main_states",
            pretty_attr_key="main_states",
            cat_key=self._rc_key,
            add_to_existing_error_msg="Compute main states first as `.compute_metastable_states()`.",
            categories=states,
            cluster_key=cluster_key,
            en_cutoff=en_cutoff,
            p_thresh=p_thresh,
            add_to_existing=add_to_existing,
        )

    def plot_metastable_states(
        self, n_cells: Optional[int] = None, same_plot: bool = True, **kwargs
    ):
        if self.schur_vectors is None:
            raise RuntimeError("Compute Schur vectors as `.metastable_state()` first.")

    def plot_coarse_T(
        self,
        show_stationary_dist: bool = True,
        show_initial_dist: bool = False,  # TODO @Marius: do we even want this?
        cmap: mcolors.ListedColormap = cm.viridis,
        xtick_rotation: float = 45,
        annotate: bool = False,
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

        labels = list(self._lin_probs.names)
        tmp = self.coarse_T

        if show_stationary_dist:
            stat_ax = fig.add_subplot(gs[1, 0])
            tmp = np.c_[tmp, self.coarse_stationary_distribution]

            stylize_dist(
                stat_ax,
                self.coarse_stationary_distribution.reshape(1, -1),
                xticks_labels=labels,
            )
            stat_ax.set_xlabel("Stationary Distribution")
        if show_initial_dist:
            init_ax = fig.add_subplot(gs[0, 1])
            tmp = np.c_[tmp, self._coarse_init_dist]

            stylize_dist(init_ax, self._coarse_init_dist.reshape(-1, 1))
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
        raise NotImplementedError()

    @property
    def schur_vectors(self):
        return self._schur_vectors

    @property
    def coarse_T(self):
        return self._coarse_T

    @property
    def coarse_stationary_distribution(self):
        return self._coarse_stat_dist

    @property
    def main_states(self):
        return self._main_states
