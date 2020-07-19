# -*- coding: utf-8 -*-
"""Generalized Perron Cluster Cluster Analysis [GPCCA18]_."""

from types import MappingProxyType
from typing import Any, Dict, List, Tuple, Union, Mapping, Iterable, Optional
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

from cellrank import logging as logg
from cellrank.tools import Lineage
from cellrank.tools._utils import (
    save_fig,
    _eigengap,
    _fuzzy_to_discrete,
    _convert_lineage_name,
    _series_from_one_hot_matrix,
)
from cellrank.tools._colors import _get_black_or_white
from cellrank.tools._constants import Lin
from cellrank.tools.estimators._utils import (
    Metadata,
    _print_insufficient_number_of_cells,
)
from cellrank.tools.estimators._property import MetaStates
from cellrank.tools.estimators._constants import A, F, P
from cellrank.tools.estimators._decomposition import Eigen, Schur
from cellrank.tools.estimators._base_estimator import BaseEstimator


class GPCCA(BaseEstimator, MetaStates, Schur, Eigen):
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

    __prop_metadata__ = [
        Metadata(
            attr=A.COARSE_T, prop=P.COARSE_T, compute_fmt=F.NO_FUNC, dtype=pd.DataFrame
        ),
        Metadata(attr=A.FIN_ABS_PROBS, prop=P.NO_PROPERTY, dtype=Lineage),
        Metadata(attr=A.COARSE_INIT_D, prop=P.COARSE_INIT_D, dtype=pd.Series),
        Metadata(attr=A.COARSE_STAT_D, prop=P.COARSE_STAT_D, dtype=pd.Series),
    ]

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
            if self._get(P.EIG) is None:
                raise RuntimeError(
                    f"Compute eigendecomposition first as `.{F.COMPUTE.fmt(P.EIG)}()` "
                    f"or `.{F.COMPUTE.fmt(P.SCHUR)}()`."
                )
            was_from_eigengap = True
            n_states = self._get(P.EIG)["eigengap"] + 1
            logg.info(f"Using `{n_states}` states based on eigengap")

        if n_states == 1:
            self._compute_meta_for_one_state(
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
                self._get(F.COMPUTE.fmt(P.SCHUR))(n_states + 1)
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
            if n_states != self._get(P.SCHUR).shape[1]:
                raise e
            logg.warning(
                f"Unable to perform the optimization using `{self._get(P.SCHUR).shape[1]}` Schur vectors. "
                f"Recomputing the decomposition"
            )

            self._get(F.COMPUTE.fmt(P.SCHUR))(
                n_states + 1,
                initial_distribution=self._gpcca.eta,
                method=self._gpcca.method,
                which=self._gpcca.z,
                alpha=self._get(P.EIG)["params"]["alpha"],
            )
            self._gpcca = self._gpcca.optimize(m=n_states)

        self._set_meta_states(
            memberships=self._gpcca.memberships,
            n_cells=n_cells,
            cluster_key=cluster_key,
            p_thresh=p_thresh,
            en_cutoff=en_cutoff,
        )

        # cache the results and make sure we don't overwrite
        self._set(A.SCHUR, self._gpcca.X)
        self._set(A.SCHUR_MAT, self._gpcca.R)

        names = self._get(P.META_PROBS).names

        self._set(
            A.COARSE_T,
            pd.DataFrame(
                self._gpcca.coarse_grained_transition_matrix,
                index=names,
                columns=names,
            ),
        )
        self._set(
            A.COARSE_INIT_D,
            pd.Series(self._gpcca.coarse_grained_input_distribution, index=names),
        )

        # careful here, in case computing the stat. dist failed
        if self._gpcca.coarse_grained_stationary_probability is not None:
            self._set(
                A.COARSE_STAT_D,
                pd.Series(
                    self._gpcca.coarse_grained_stationary_probability, index=names,
                ),
            )
        else:
            logg.warning("No stationary distribution found in GPCCA object")

        logg.info(
            f"Adding `.{P.SCHUR}`\n"
            f"       `.{P.COARSE_T}`\n"
            f"       `.{P.COARSE_STAT_D}`\n"
            f"    Finish",
            time=start,
        )

    def set_final_states_from_metastable_states(
        self,
        names: Optional[Union[Iterable[str], str]] = None,
        redistribute: bool = True,
        n_cells: int = 30,
        **kwargs,
    ):
        """
        Manually select the main states from the metastable states.

        Params
        ------
        names
            Names of the main states. Multiple states can be combined using `','`, such as `['Alpha, Beta', 'Epsilon']`.
        redistribute
            Whether to redistribute the probability mass of unselected lineages or create a `'rest'` lineage.
        n_cells
            Number of most likely cells from each main state to select. If `None`, on main states will be selected,
            only the lineage probabilities may be redistributed.
        kwargs
            Keyword arguments for :meth:`cellrank.tl.Lineage.reduce` when redistributing the probability mass.

        Returns
        -------
        None
            Nothing, just sets the final states.
        """

        if not isinstance(n_cells, int):
            raise TypeError(
                f"Expected `n_cells` to be of type `int`, found `{type(n_cells).__name__}`."
            )

        if n_cells <= 0:
            raise ValueError(f"Expected `n_cells` to be positive, found `{n_cells}`.")

        probs = self._get(P.META_PROBS)

        if names is None:
            names = probs.names
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
            meta_states_probs = probs[names + [Lin.OTHERS]]
            meta_states_probs = meta_states_probs.reduce(
                [" or ".join(_convert_lineage_name(name)) for name in names], **kwargs
            )
        else:
            meta_states_probs = probs[names + [Lin.REST]]

        # compute the aggregated probability of being a root/final state (no matter which)
        scaled_probs = meta_states_probs[
            [n for n in meta_states_probs.names if n != "rest"]
        ].X.copy()
        scaled_probs /= scaled_probs.max(0)

        self._set(A.FIN, self._create_states(meta_states_probs, n_cells))
        self._set(
            A.FIN_PROBS, pd.Series(scaled_probs.max(1), index=self.adata.obs_names)
        )
        self._set(
            A.FIN_COLORS,
            meta_states_probs[list(self._get(P.FIN).cat.categories)].colors,
        )

        # TODO: do we even want this?
        self._set(A.FIN_ABS_PROBS, meta_states_probs)

        self._write_final_states()

    def compute_final_states(
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
            One of following:

                - `'eigengap'` - select the number of states based on the eigengap of the transition matrix
                - `'eigengap_coarse'` - select the number of states based on the eigengap of the diagonal of the coarse-grained transition matrix
                - `'min_self_prob'` - select states which have the given minimum probability on the diagonal of the coarse-grained transition matrix
                - `'top_n'` - select top :paramref:`n_main_states` based on the probability on the diagonal of the coarse-grained transition matrix
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

                - :paramref:`final_states_probabilities`
                - :paramref:`final_states`
        """  # noqa
        if len(self._get(P.META).cat.categories) == 1:
            logg.warning(
                "Found only one metastable state. Making it the single main state. "
            )
            self.set_final_states_from_metastable_states(
                None, redistribute=False, n_cells=n_cells, **kwargs
            )
            return

        coarse_T = self._get(P.COARSE_T)

        if method == "eigengap":
            if self._get(P.EIG) is None:
                raise RuntimeError(
                    "Compute eigendecomposition first as `.compute_eig()`."
                )
            n_main_states = _eigengap(self._get(P.EIG)["D"], alpha=alpha) + 1
        elif method == "eigengap_coarse":
            if coarse_T is None:
                raise RuntimeError(
                    "Compute metastable states first as `.compute_metastable_states()`."
                )
            n_main_states = _eigengap(np.sort(np.diag(coarse_T)[::-1]), alpha=alpha)
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
            self_probs = pd.Series(np.diag(coarse_T), index=coarse_T.columns)
            names = self_probs[self_probs.values >= min_self_prob].index
            self.set_final_states_from_metastable_states(
                names, redistribute=redistribute, n_cells=n_cells, **kwargs
            )
            return
        else:
            raise ValueError(
                f"Invalid method `{method!r}`. Valid options are `'eigengap', 'eigengap_coarse', "
                f"'top_n' and 'min_self_prob'`."
            )

        names = coarse_T.columns[np.argsort(np.diag(coarse_T))][-n_main_states:]
        self.set_final_states_from_metastable_states(
            names, redistribute=redistribute, n_cells=n_cells, **kwargs
        )

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
            Keyword arguments for :meth:`cellrank.tl.GPCCA.compute_schur` if Schur decomposition is not found.

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

        if self._get(P.SCHUR) is None:
            logg.warning("No Schur decomposition found. Computing")
            self._get(F.COMPUTE.fmt(P.SCHUR))(n_components, **kwargs)
        elif self._get(P.SCHUR_MAT).shape[1] < n_components:
            logg.warning(
                f"Requested `{n_components}` components, but only `{self._get(P.SCHUR_MAT).shape[1]}` were found. "
                f"Recomputing using default values"
            )
            self._get(F.COMPUTE.fmt(P.SCHUR))(n_components)
        else:
            logg.debug("Using cached Schur decomposition")

        start = logg.info(
            f"Computing Generalized Diffusion Pseudotime using n_components = {n_components}"
        )

        Q, eigenvalues = (
            self._get(P.SCHUR),
            self._get(P.EIG)["D"],
        )
        # may have to remove some values if too many converged
        Q, eigenvalues = Q[:, :n_components], eigenvalues[:n_components]
        D = _get_dpt_row(eigenvalues, Q, i=iroot)
        pseudotime = D / np.max(D[np.isfinite(D)])

        self.adata.obs[key_added] = pseudotime

        logg.info(f"Adding `{key_added!r}` to `adata.obs`\n    Finish", time=start)

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
        dpi: int = 80,
        save: Optional[Union[Path, str]] = None,
        text_kwargs: Mapping[str, Any] = MappingProxyType({}),
        **kwargs,
    ) -> None:
        """
        Plot the coarse-grained transition matrix of the metastable states.

        .. image:: https://raw.githubusercontent.com/theislab/cellrank/master/resources/images/coarse_T.png
           :alt: image of coarse transition matrix
           :width: 400px
           :align: center

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

        coarse_T = self._get(P.COARSE_T)
        coarse_stat_d = self._get(P.COARSE_STAT_D)
        coarse_init_d = self._get(P.COARSE_INIT_D)

        if coarse_T is None:
            raise RuntimeError(
                f"Compute coarse-grained transition matrix first as `.{F.COMPUTE.fmt(P.META)}()` with `n_states` > 1."
            )

        if show_stationary_dist and coarse_stat_d is None:
            logg.warning("Coarse stationary distribution is `None`, not plotting")
            show_stationary_dist = False
        if show_initial_dist and coarse_init_d is None:
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

        labels = list(self.coarse_T.columns)

        tmp = coarse_T
        if show_initial_dist:
            tmp = np.c_[tmp, coarse_stat_d]
        if show_initial_dist:
            tmp = np.c_[tmp, coarse_init_d]
        norm = mpl.colors.Normalize(vmin=np.nanmin(tmp), vmax=np.nanmax(tmp))

        if show_stationary_dist:
            stat_ax = fig.add_subplot(gs[1, 0])
            stylize_dist(
                stat_ax, np.array(coarse_stat_d).reshape(1, -1), xticks_labels=labels,
            )
            stat_ax.set_xlabel("stationary distribution")  # , ha="right", x=1)

        if show_initial_dist:
            init_ax = fig.add_subplot(gs[0, 1])
            stylize_dist(init_ax, np.array(coarse_init_d).reshape(-1, 1))

            init_ax.yaxis.set_label_position("right")
            init_ax.set_ylabel("initial distribution", rotation=-90, va="bottom")

        im = ax.imshow(coarse_T, aspect="auto", cmap=cmap, **kwargs)
        ax.set_title("coarse-grained transition matrix" if title is None else title)

        if show_cbar:
            _ = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)

        ax.set_yticks(np.arange(coarse_T.shape[0]))
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
            ax.set_xticks(np.arange(coarse_T.shape[1]))
            ax.set_xticklabels(labels)
            plt.setp(
                ax.get_xticklabels(),
                rotation=xtick_rotation,
                ha="right",
                rotation_mode="anchor",
            )
        else:
            ax.set_xticks([])

        ax.set_yticks(np.arange(coarse_T.shape[0] + 1) - 0.5, minor=True)
        ax.tick_params(
            which="minor", bottom=not show_stationary_dist, left=False, top=False
        )

        if annotate:
            annotate_heatmap(im)
            annotate_dist_ax(stat_ax, coarse_stat_d.values, is_vertical=False)
            annotate_dist_ax(init_ax, coarse_init_d, is_vertical=True)

        if save:
            save_fig(fig, save)

        fig.show()

    def _compute_meta_for_one_state(
        self,
        n_cells: int,
        cluster_key: Optional[str],
        en_cutoff: Optional[float],
        p_thresh: float,
    ) -> None:
        start = logg.info("Computing metastable states")
        logg.warning("For `n_states=1`, stationary distribution is computed")

        self._get(F.COMPUTE.fmt(P.EIG))(only_evals=False, which="LM")
        stationary_dist = self._get(P.EIG)["stationary_dist"]

        self._set_meta_states(
            memberships=stationary_dist[:, None],
            n_cells=n_cells,
            cluster_key=cluster_key,
            p_thresh=p_thresh,
            en_cutoff=en_cutoff,
        )

        # reset all the things
        for key in (
            A.META_PROBS,
            A.ABS_RPOBS,
            A.SCHUR,
            A.SCHUR_MAT,
            A.COARSE_T,
            A.COARSE_STAT_D,
            A.COARSE_STAT_D,
        ):
            self._set(key.s, None)

        logg.info(f"Adding `.{P.META}`\n    Finish", time=start)

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

        logg.debug(f"Calculating minChi within interval `[{minn}, {maxx}]`")

        return int(np.arange(minn, maxx)[np.argmax(self._gpcca.minChi(minn, maxx))])

    def _set_meta_states(
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
            Fuzzy clustering.
        n_cells
            Number of cells to be used to represent each state.
        cluster_key
            Key from :paramref:`adata` `.obs` to get reference cluster annotations.
        en_cutoff
            Threshold to decide when we we want to warn the user about an uncertain name mapping. This happens when
            one fuzzy state overlaps with several reference clusters, and the most likely cells are distributed almost
            evenly across the reference clusters.
        p_thresh
            Only used to detect cell cycle stages. These have to be present in :paramref:`adata` `.obs`
            as `G2M_score` and `S_score`.
        check_row_sums
            Check whether rows in `memberships` sum to 1.

        Returns
        --------
        None
            Writes a lineage object which mapped names and colors. Also creates a categorical :class:`pandas.Series`
            `.metastable_states`, where the top `n_cells` cells represent each fuzzy state.
        """

        if n_cells is None:
            logg.debug("Setting the metastable states using metastable assignment")

            max_assignment = np.argmax(memberships, axis=1)
            _meta_assignment = pd.Series(
                index=self.adata.obs_names, data=max_assignment, dtype="category"
            )
            # sometimes, the assignment can have a missing category and the Lineage creation therefore fails
            # keep it as ints when `n_cells != None`
            _meta_assignment.cat.set_categories(
                list(range(memberships.shape[1])), inplace=True
            )

            metastable_states = _meta_assignment.astype(str).astype("category").copy()
            not_enough_cells = []
        else:
            logg.debug("Setting the metastable states using metastable memberships")

            # select the most likely cells from each metastable state
            metastable_states, not_enough_cells = self._create_states(
                memberships,
                n_cells=n_cells,
                check_row_sums=check_row_sums,
                return_not_enough_cells=True,
            )
            not_enough_cells = not_enough_cells.astype("str")

        # _set_categorical_labels creates the names, we still need to remap the group names
        orig_cats = metastable_states.cat.categories
        self._set_categorical_labels(
            attr_key=A.META.v,
            color_key=A.META_COLORS.v,
            pretty_attr_key=P.META.v,
            add_to_existing_error_msg=f"Compute metastable states first as `.{F.COMPUTE.fmt(P.META)}()`.",
            categories=metastable_states,
            cluster_key=cluster_key,
            en_cutoff=en_cutoff,
            p_thresh=p_thresh,
            add_to_existing=False,
        )

        name_mapper = dict(zip(orig_cats, self._get(P.META).cat.categories))
        _print_insufficient_number_of_cells(
            [name_mapper.get(n, n) for n in not_enough_cells], n_cells
        )

        logg.debug(
            "Setting metastable lineage probabilities based on GPCCA membership vectors"
        )

        self._set(
            A.META_PROBS,
            Lineage(
                memberships,
                names=list(metastable_states.cat.categories),
                colors=self._get(A.META_COLORS),
            ),
        )

    def _create_states(
        self,
        probs: Union[np.ndarray, Lineage],
        n_cells: int,
        check_row_sums: bool = False,
        return_not_enough_cells: bool = False,
    ) -> pd.Series:
        if n_cells <= 0:
            raise ValueError(f"Expected `n_cells` to be positive, found `{n_cells}`.")

        if isinstance(probs, Lineage):
            probs = probs[[n for n in probs.names if n != "rest"]]

        a_discrete, not_enough_cells = _fuzzy_to_discrete(
            a_fuzzy=probs,
            n_most_likely=n_cells,
            remove_overlap=False,
            raise_threshold=0.2,
            check_row_sums=check_row_sums,
        )

        states = _series_from_one_hot_matrix(
            a=a_discrete,
            index=self.adata.obs_names,
            names=probs.names if isinstance(probs, Lineage) else None,
        )

        return (states, not_enough_cells) if return_not_enough_cells else states

    # TODO: docrep
    def copy(self) -> "GPCCA":
        """Return a copy of self."""
        raise NotImplementedError()

    # TODO: just call super + extra
    # or save meta states as well (NYI) and handle it here
    def __getstate__(self):
        pass

    # TODO: call super + extra
    def __setstate__(self, state):
        pass
