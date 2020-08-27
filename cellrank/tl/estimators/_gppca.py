# -*- coding: utf-8 -*-
"""Generalized Perron Cluster Cluster Analysis [GPCCA18]_."""

from types import MappingProxyType
from typing import Any, Dict, List, Tuple, Union, Mapping, Iterable, Optional, Sequence
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

from cellrank import logging as logg
from cellrank.ul._docs import d, inject_docs
from cellrank.tl._utils import (
    save_fig,
    _eigengap,
    _fuzzy_to_discrete,
    _series_from_one_hot_matrix,
)
from cellrank.tl._colors import _get_black_or_white
from cellrank.tl._lineage import Lineage
from cellrank.tl.estimators._utils import Metadata, _print_insufficient_number_of_cells
from cellrank.tl.estimators._property import MetaStates
from cellrank.tl.estimators._constants import A, F, P
from cellrank.tl.estimators._decomposition import Eigen, Schur
from cellrank.tl.estimators._base_estimator import BaseEstimator


@d.dedent
class GPCCA(BaseEstimator, MetaStates, Schur, Eigen):
    """
    Generalized Perron Cluster Cluster Analysis [GPCCA18]_.

    Parameters
    ----------
    %(base_estimator.parameters)s
    """

    __prop_metadata__ = [
        Metadata(
            attr=A.COARSE_T,
            prop=P.COARSE_T,
            compute_fmt=F.NO_FUNC,
            plot_fmt=F.NO_FUNC,
            dtype=pd.DataFrame,
            doc="Coarse-grained transition matrix.",
        ),
        Metadata(attr=A.FIN_ABS_PROBS, prop=P.NO_PROPERTY, dtype=Lineage),
        Metadata(attr=A.COARSE_INIT_D, prop=P.COARSE_INIT_D, dtype=pd.Series),
        Metadata(attr=A.COARSE_STAT_D, prop=P.COARSE_STAT_D, dtype=pd.Series),
    ]

    @inject_docs(
        ms=P.META,
        msp=P.META_PROBS,
        schur=P.SCHUR.s,
        coarse_T=P.COARSE_T,
        coarse_stat=P.COARSE_STAT_D,
    )
    @d.dedent
    def compute_metastable_states(
        self,
        n_states: Optional[
            Union[int, Tuple[int, int], List[int], Dict[str, int]]
        ] = None,
        n_cells: Optional[int] = 30,
        use_min_chi: bool = False,
        cluster_key: str = None,
        en_cutoff: Optional[float] = 0.7,
        p_thresh: float = 1e-15,
    ):
        """
        Compute the metastable states.

        Parameters
        ----------
        n_states
            Number of metastable states. If `None`, use the `eigengap` heuristic.
        %(n_cells)s
        use_min_chi
            Whether to use :meth:`msmtools.analysis.dense.gpcca.GPCCA.minChi` to calculate the number of metastable
            states. If `True`, ``n_states`` corresponds to an interval `[min, max]` inside of which
            the potentially optimal number of metastable states is searched.
        cluster_key
            If a key to cluster labels is given, names and colors of the states will be associated with the clusters.
        en_cutoff
            If ``cluster_key`` is given, this parameter determines when an approximate recurrent class will
            be labelled as *'Unknown'*, based on the entropy of the distribution of cells over transcriptomic clusters.
        p_thresh
            If cell cycle scores were provided, a *Wilcoxon rank-sum test* is conducted to identify cell-cycle driven
            start- or endpoints.
            If the test returns a positive statistic and a p-value smaller than ``p_thresh``, a warning will be issued.

        Returns
        -------
        None
            Nothing, but updates the following fields:

                - :paramref:`{msp}`
                - :paramref:`{ms}`
                - :paramref:`{schur}`
                - :paramref:`{coarse_T}`
                - :paramref:`{coarse_stat}`
        """

        was_from_eigengap = False
        if n_states is None:
            if self._get(P.EIG) is None:
                raise RuntimeError(
                    "Compute eigendecomposition first as `.compute_eigendecomposition()` or `.compute_schur()`."
                )
            was_from_eigengap = True
            n_states = self._get(P.EIG)["eigengap"] + 1
            logg.info(f"Using `{n_states}` states based on eigengap")

        if n_states <= 0:
            raise ValueError(
                f"Expected `n_states` to be positive or `None`, found `{n_states}`."
            )

        if self._invalid_n_states is not None and n_states in self._invalid_n_states:
            logg.warning(
                f"Unable to compute metastable states with `n_states={n_states}` because it will "
                f"split the conjugate eigenvalues. Increasing `n_states` to `{n_states + 1}`"
            )
            n_states += 1  # cannot force recomputation of Schur decomposition
            assert n_states not in self._invalid_n_states, "Sanity check failed."

        if n_states == 1:
            self._compute_meta_for_one_state(
                n_cells=n_cells,
                cluster_key=cluster_key,
                p_thresh=p_thresh,
                en_cutoff=en_cutoff,
            )
            return

        if self._gpcca is None:
            if not was_from_eigengap:
                raise RuntimeError(
                    "Compute Schur decomposition first as `.compute_schur()`."
                )

            logg.warning(
                f"Number of states `{n_states}` was automatically determined by `eigengap` "
                "but no Schur decomposition was found. Computing with default parameters"
            )
            # this cannot fail if splitting occurs
            # if it were to split, it's automatically increased in `compute_schur`
            self.compute_schur(n_states + 1)

        if use_min_chi:
            n_states = self._get_n_states_from_minchi(n_states)
        elif not isinstance(n_states, int):
            raise ValueError(
                f"Expected `n_states` to be an integer when `use_min_chi=False`, found `{type(n_states).__name__!r}`."
            )

        if self._gpcca.X.shape[1] < n_states:
            logg.warning(
                f"Requested more metastable states `{n_states}` than available "
                f"Schur vectors `{self._gpcca.X.shape[1]}`. Recomputing the decomposition"
            )

        start = logg.info("Computing metastable states")
        try:
            self._gpcca = self._gpcca.optimize(m=n_states)
        except ValueError as e:
            # this is the following cage - we have 4 Schur vectors, user requests 5 states, but it splits the conj. ev.
            # in the try block, schur decomposition with 5 vectors is computed, but it fails (no way of knowing)
            # so in this case, we increate it by 1
            n_states += 1
            logg.warning(f"{e}\nIncreasing `n_states` to `{n_states}`")
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
            logg.info(
                f"Adding `.{P.META_PROBS}`\n"
                f"       `.{P.META}`\n"
                f"       `.{P.SCHUR}`\n"
                f"       `.{P.COARSE_T}`\n"
                f"       `.{P.COARSE_STAT_D}`\n"
                f"    Finish",
                time=start,
            )
        else:
            logg.warning("No stationary distribution found in GPCCA object")
            logg.info(
                f"Adding `.{P.META_PROBS}`\n"
                f"       `.{P.META}`\n"
                f"       `.{P.SCHUR}`\n"
                f"       `.{P.COARSE_T}`\n"
                f"    Finish",
                time=start,
            )

    @d.dedent
    @inject_docs(fs=P.FIN, fsp=P.FIN_PROBS)
    def set_final_states_from_metastable_states(
        self, names: Optional[Union[Iterable[str], str]] = None, n_cells: int = 30,
    ):
        """
        Manually select the main states from the metastable states.

        Parameters
        ----------
        names
            Names of the main states. Multiple states can be combined using `','`, such as `['Alpha, Beta', 'Epsilon']`.
        %(n_cells)s

        Returns
        -------
        None
            Nothing, just updates the following fields:

                - :paramref:`{fsp}`
                - :paramref:`{fs}`
        """

        if not isinstance(n_cells, int):
            raise TypeError(
                f"Expected `n_cells` to be of type `int`, found `{type(n_cells).__name__}`."
            )

        if n_cells <= 0:
            raise ValueError(f"Expected `n_cells` to be positive, found `{n_cells}`.")

        probs = self._get(P.META_PROBS)

        if self._get(P.META_PROBS) is None:
            raise RuntimeError(
                "Compute metastable_states first as `.compute_metastable_states()`."
            )
        elif probs.shape[1] == 1:
            self._set(A.FIN, self._create_states(probs, n_cells=n_cells))
            self._set(A.FIN_COLORS, self._get(A.META_COLORS))
            self._set(A.FIN_PROBS, probs / probs.max())
            self._set(A.FIN_ABS_PROBS, probs)
            self._write_final_states()
            return

        if names is None:
            names = probs.names

        if isinstance(names, str):
            names = [names]

        meta_states_probs = probs[list(names)]

        # compute the aggregated probability of being a root/final state (no matter which)
        scaled_probs = meta_states_probs[
            [n for n in meta_states_probs.names if n != "rest"]
        ].copy()
        scaled_probs /= scaled_probs.max(0)

        self._set(A.FIN, self._create_states(meta_states_probs, n_cells))
        self._set(
            A.FIN_PROBS, pd.Series(scaled_probs.X.max(1), index=self.adata.obs_names)
        )
        self._set(
            A.FIN_COLORS,
            meta_states_probs[list(self._get(P.FIN).cat.categories)].colors,
        )

        self._set(A.FIN_ABS_PROBS, scaled_probs)
        self._write_final_states()

    @inject_docs(fs=P.FIN, fsp=P.FIN_PROBS)
    @d.dedent
    def compute_final_states(
        self,
        method: str = "eigengap",
        n_cells: int = 30,
        alpha: Optional[float] = 1,
        min_self_prob: Optional[float] = None,
        n_main_states: Optional[int] = None,
    ):
        """
        Automatically select the main states from metastable states.

        Parameters
        ----------
        method
            One of following:

                - `'eigengap'` - select the number of states based on the eigengap of the transition matrix.
                - `'eigengap_coarse'` - select the number of states based on the eigengap of the diagonal \
                    of the coarse-grained transition matrix.
                - `'top_n'` - select top :paramref:`n_main_states` based on the probability of the diagonal \
                    of the coarse-grained transition matrix.
                - `'min_self_prob'` - select states which have the given minimum probability of the diagonal \
                    of the coarse-grained transition matrix.
        %(n_cells)s
        alpha
            Weight given to the deviation of an eigenvalue from one. Used when ``method='eigengap'``
            or ``method='eigengap_coarse'``.
        min_self_prob
            Used when ``method='min_self_prob'``.
        n_main_states
            Used when ``method='top_n'``.

        Returns
        -------
        None
            Nothing, just updates the following fields:

                - :paramref:`{fsp}`
                - :paramref:`{fs}`
        """  # noqa

        if len(self._get(P.META).cat.categories) == 1:
            logg.warning(
                "Found only one metastable state. Making it the single main state"
            )
            self.set_final_states_from_metastable_states(None, n_cells=n_cells)
            return

        coarse_T = self._get(P.COARSE_T)

        if method == "eigengap":
            if self._get(P.EIG) is None:
                raise RuntimeError(
                    "Compute eigendecomposition first as `.compute_eigendecomposition()`."
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
                    "Argument `n_main_states` must be != `None` for `method='top_n'`."
                )
            elif n_main_states <= 0:
                raise ValueError(
                    f"Expected `n_main_states` to be positive, found `{n_main_states}`."
                )
        elif method == "min_self_prob":
            if min_self_prob is None:
                raise ValueError(
                    "Argument `min_self_prob` must be != `None` for `method='min_self_prob'`."
                )
            self_probs = pd.Series(np.diag(coarse_T), index=coarse_T.columns)
            names = self_probs[self_probs.values >= min_self_prob].index
            self.set_final_states_from_metastable_states(names, n_cells=n_cells)
            return
        else:
            raise ValueError(
                f"Invalid method `{method!r}`. Valid options are `'eigengap', 'eigengap_coarse', "
                f"'top_n' and 'min_self_prob'`."
            )

        names = coarse_T.columns[np.argsort(np.diag(coarse_T))][-n_main_states:]
        self.set_final_states_from_metastable_states(names, n_cells=n_cells)

    def compute_gdpt(
        self, n_components: int = 10, key_added: str = "gdpt_pseudotime", **kwargs
    ):
        """
        Compute generalized Diffusion pseudotime from [Haghverdi16]_ making use of the real Schur decomposition.

        Parameters
        ----------
        n_components
            Number of real Schur vectors to consider.
        key_added
            Key in :paramref:`adata` ``.obs`` where to save the pseudotime.
        **kwargs
            Keyword arguments for :meth:`cellrank.tl.GPCCA.compute_schur` if Schur decomposition is not found.

        Returns
        -------
        None
            Nothing, just updates :paramref:`adata` ``.obs[`key_added]`` with the computed pseudotime.
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
            self.compute_schur(n_components, **kwargs)
        elif self._get(P.SCHUR_MAT).shape[1] < n_components:
            logg.warning(
                f"Requested `{n_components}` components, but only `{self._get(P.SCHUR_MAT).shape[1]}` were found. "
                f"Recomputing using default values"
            )
            self.compute_schur(n_components)
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

    @d.dedent
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

        Parameters
        ----------
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
        %(plotting)s
        text_kwargs
            Keyword arguments for :func:`matplotlib.pyplot.text`.
        **kwargs
            Keyword arguments for :func:`matplotlib.pyplot.imshow`.

        Returns
        -------
        %(just_plots)s
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
                "Compute coarse-grained transition matrix first as `.compute_metastable_states()` with `n_states > 1`."
            )

        if show_stationary_dist and coarse_stat_d is None:
            logg.warning("Coarse stationary distribution is `None`, ignoring")
            show_stationary_dist = False
        if show_initial_dist and coarse_init_d is None:
            logg.warning("Coarse initial distribution is `None`, ignoring")
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
            stat_ax.set_xlabel("stationary distribution")

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

        eig = self._get(P.EIG)
        if (
            eig is not None
            and "stationary_dist" in eig
            and eig["params"]["which"] == "LM"
        ):
            stationary_dist = eig["stationary_dist"]
        else:
            self.compute_eigendecomposition(only_evals=False, which="LM")
            stationary_dist = self._get(P.EIG)["stationary_dist"]

        self._set_meta_states(
            memberships=stationary_dist[:, None],
            n_cells=n_cells,
            cluster_key=cluster_key,
            p_thresh=p_thresh,
            en_cutoff=en_cutoff,
        )
        self._set(
            A.META_PROBS,
            Lineage(
                stationary_dist,
                names=list(self._get(A.META).cat.categories),
                colors=self._get(A.META_COLORS),
            ),
        )

        # reset all the things
        for key in (
            A.ABS_RPOBS,
            A.SCHUR,
            A.SCHUR_MAT,
            A.COARSE_T,
            A.COARSE_STAT_D,
            A.COARSE_STAT_D,
        ):
            self._set(key.s, None)

        logg.info(
            f"Adding `.{P.META_PROBS}`\n" f"        .{P.META}\n" f"    Finish",
            time=start,
        )

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
                "If you really expect 2 clusters, use `n_states=2`. Setting the minimum to 3"
            )
            minn = 3

        logg.debug(f"Calculating minChi within interval `[{minn}, {maxx}]`")

        return int(np.arange(minn, maxx)[np.argmax(self._gpcca.minChi(minn, maxx))])

    @d.dedent
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

        Parameters
        --------
        memberships
            Fuzzy clustering.
        %(n_cells)s
        cluster_key
            Key from :paramref:`adata` ``.obs`` to get reference cluster annotations.
        en_cutoff
            Threshold to decide when we we want to warn the user about an uncertain name mapping. This happens when
            one fuzzy state overlaps with several reference clusters, and the most likely cells are distributed almost
            evenly across the reference clusters.
        p_thresh
            Only used to detect cell cycle stages. These have to be present in
            :paramref:`adata` ``.obs`` as `'G2M_score'` and `'S_score'`.
        check_row_sums
            Check whether rows in `memberships` sum to `1`.

        Returns
        -------
        None
            Writes a :class:`cellrank.tl.Lineage` object which mapped names and colors.
            Also writes a categorical :class:`pandas.Series`, where top ``n_cells`` cells represent each fuzzy state.
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
            add_to_existing_error_msg="Compute metastable states first as `.compute_metastable_states()`.",
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
            membership=a_discrete,
            index=self.adata.obs_names,
            names=probs.names if isinstance(probs, Lineage) else None,
        )

        return (states, not_enough_cells) if return_not_enough_cells else states

    def _fit_final_states(
        self,
        n_lineages: Optional[int] = None,
        cluster_key: Optional[str] = None,
        method: str = "krylov",
        **kwargs,
    ) -> None:
        if n_lineages is None or n_lineages == 1:
            self.compute_eigendecomposition()
            if n_lineages is None:
                n_lineages = self.eigendecomposition["eigengap"] + 1

        if n_lineages > 1:
            self.compute_schur(n_lineages + 1, method=method)

        try:
            self.compute_metastable_states(
                n_states=n_lineages, cluster_key=cluster_key, **kwargs
            )
        except ValueError:
            logg.warning(
                f"Computing `{n_lineages}` metastable states cuts through a block of complex conjugates. "
                f"Increasing `n_lineages` to {n_lineages + 1}"
            )
            self.compute_metastable_states(
                n_states=n_lineages + 1, cluster_key=cluster_key, **kwargs
            )

        fs_kwargs = {"n_cells": kwargs["n_cells"]} if "n_cells" in kwargs else {}

        if n_lineages is None:
            self.compute_final_states(method="eigengap", **fs_kwargs)
        else:
            self.set_final_states_from_metastable_states(**fs_kwargs)

    @d.dedent  # because of fit
    @d.dedent
    @inject_docs(
        ms=P.META,
        msp=P.META_PROBS,
        fs=P.FIN,
        fsp=P.FIN_PROBS,
        ap=P.ABS_PROBS,
        dp=P.DIFF_POT,
    )
    def fit(
        self,
        n_lineages: Optional[int] = None,
        cluster_key: Optional[str] = None,
        keys: Optional[Sequence[str]] = None,
        method: str = "krylov",
        compute_absorption_probabilities: bool = True,
        **kwargs,
    ):
        """
        Run the pipeline, computing the metastable states, %(final)s states and optionally the absorption probabilities.

        It is equivalent to running::

            compute_eigendecomposition(...)  # if needed
            compute_schur(...)
            compute_metastable_states(...)

            compute_final_states(...)   # if n_lineages=None
            set_final_states_from_metastable_states(...)   # otherwise

            compute_absorption_probabilities(...)  # optional

        Parameters
        ----------
        %(fit)s
        method
            Method to use when computing the Schur decomposition. Valid options are: `'krylov'` or `'brandts'`.
        compute_absorption_probabilities
            Whether to compute absorption probabilities or only final states.
        **kwargs
            Keyword arguments for :meth:`cellrank.tl.GPCCA.compute_metastable_states`.

        Returns
        -------
        None
            Nothing, just makes available the following fields:

                - :paramref:`{msp}`
                - :paramref:`{ms}`
                - :paramref:`{fsp}`
                - :paramref:`{fs}`
                - :paramref:`{ap}`
                - :paramref:`{dp}`
        """

        super().fit(
            n_lineages=n_lineages,
            cluster_key=cluster_key,
            keys=keys,
            method=method,
            compute_absorption_probabilities=compute_absorption_probabilities,
            **kwargs,
        )
