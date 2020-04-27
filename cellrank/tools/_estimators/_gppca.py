# -*- coding: utf-8 -*-
from typing import Optional, List, Tuple, Dict, Union
from anndata import AnnData
from msmtools.analysis.dense.gpcca import GPCCA as _GPPCA
from scanpy import logging as logg

from cellrank.tools._estimators._base_estimator import BaseEstimator
from cellrank.tools._utils import save_fig
from cellrank.tools.kernels._kernel import KernelExpression

import os
import numpy as np
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
        self._gpcca = None
        self._schur_vectors = None
        self._coarse_T = None
        self._coarse_init_dist = None
        self._coarse_stat_dist = None

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
        Plt Schur vectors in an embedding.

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
            raise RuntimeError("Compute Schur vectors as `.metastable_state()` first.")

        self._plot_vectors(
            self.schur_vectors,
            "schur",
            abs_value=abs_value,
            use=use,
            cluster_key=cluster_key,
            **kwargs,
        )

    def metastable_states(
        self,
        n_states: Union[int, Tuple[int, int], List[int], Dict[str, int]],
        initial_distribution: Optional[np.ndarray] = None,
        use_min_chi: bool = False,
        method: str = "krylov",
        which: str = "LM",
        cluster_key: Optional[str] = "louvain",
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

        Returns
        -------
        None
            Nothings, but updates the following fields:

            - TODO
        """
        if (
            use_min_chi
        ):  # TODO: @Marius - this is the cleanest option I could thought of
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
            logg.debug(f"DEBUG: Calculating min Chi in interval [{minn}, {maxx}]")
            n_states = np.arange(minn, maxx)[np.argmax(self._gpcca.minChi(minn, maxx))]

        start = logg.info("Computing metastable states")

        self._gpcca: _GPPCA = _GPPCA(
            self._T, eta=initial_distribution, z=which, method=method
        ).optimize(m=n_states)

        self._schur_vectors = self._gpcca.schur_vectors
        self._coarse_T = self._gpcca.coarse_grained_transition_matrix
        self._coarse_init_dist = self._gpcca.coarse_grained_input_distribution
        self._coarse_stat_dist = self._gpcca.coarse_grained_stationary_probability

        logg.info("Adding `...`\n" "    Finish", time=start)

    def plot_metastable_states(
        self, n_cells: Optional[int] = None, same_plot: bool = True, **kwargs
    ):
        raise NotImplementedError()

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
        **kwargs,
    ) -> None:
        """
        TODO

        Params
        ------
        show_stationary_dist
        show_initial_dist
        cmap
        xtick_rotation
        annotate
        show_cbar
        figsize
        dpi
        save
        kwargs

        Returns
        -------
        None
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
            **textkw,
        ):
            # modified from matplotlib's site

            data = im.get_array()
            # Normalize the threshold to the images color range.
            threshold = (
                im.norm(np.median(data)) if threshold is None else im.norm(threshold)
            )
            kw = dict(horizontalalignment="center", verticalalignment="center")
            kw.update(textkw)

            # Get the formatter in case a string is supplied
            if isinstance(valfmt, str):
                valfmt = mpl.ticker.StrMethodFormatter(valfmt)

            # Loop over the data and create a `Text` for each "pixel".
            # Change the text's color depending on the data.
            texts = []
            for i in range(data.shape[0]):
                for j in range(data.shape[1]):
                    # TODO: @Marius do we want to change the color based on thresh?
                    kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
                    text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
                    texts.append(text)

        if self.coarse_T is None:
            raise RuntimeError(
                f"Compute coarse transition matrix first as `.metastable_states()`."
            )

        if show_stationary_dist and self.coarse_stat_dist is None:
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

        labels = [
            str(i) for i in range(self.coarse_T.shape[0])
        ]  # TODO: names of metastable states
        tmp = self.coarse_T

        if show_stationary_dist:
            stat_ax = fig.add_subplot(gs[1, 0])
            tmp = np.c_[tmp, self.coarse_stat_dist]

            stylize_dist(
                stat_ax, self.coarse_stat_dist.reshape(1, -1), xticks_labels=labels
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

    def set_main_states(self, keys: Optional[List[str]] = None, **kwargs):
        raise NotImplementedError()

    def copy(self) -> "GPCCA":
        raise NotImplementedError()

    @property
    def schur_vectors(self):
        return self._schur_vectors

    @property
    def coarse_T(self):
        return self._coarse_T

    @property
    def coarse_stat_dist(self):
        return self._coarse_stat_dist
