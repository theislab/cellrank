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
from cellrank.tl._constants import TermStatesKey, _probs, _colors, _lin_names
from cellrank.tl.estimators._utils import Metadata, _print_insufficient_number_of_cells
from cellrank.tl.estimators._property import Macrostates
from cellrank.tl.estimators._constants import A, F, P
from cellrank.tl.estimators._decomposition import Eigen, Schur
from cellrank.tl.estimators._base_estimator import BaseEstimator


@d.dedent
class GPCCA(BaseEstimator, Macrostates, Schur, Eigen):
    """
    Generalized Perron Cluster Cluster Analysis [GPCCA18]_ as implemented in `pyGPCCA <https://pygpcca.readthedocs.io/en/latest/>`_.

    Coarse-grains a discrete Markov chain into a set of macrostates and computes coarse-grained transition probabilities
    among the macrostates. Each macrostate corresponds to an area of the state space, i.e. to a subset of cells. The
    assignment is soft, i.e. each cell is assigned to every macrostate with a certain weight, where weights sum to
    one per cell. Macrostates are computed by maximizing the 'crispness' which can be thought of as a measure for
    minimal overlap between macrostates in a certain inner-product sense. Once the macrostates have been computed,
    we project the large transition matrix onto a coarse-grained transition matrix among the macrostates via
    a Galerkin projection. This projection is based on invariant subspaces of the original transition matrix which
    are obtained using the real Schur decomposition [GPCCA18]_.

    Parameters
    ----------
    %(base_estimator.parameters)s
    """  # noqa: E501

    __prop_metadata__ = [
        Metadata(
            attr=A.COARSE_T,
            prop=P.COARSE_T,
            compute_fmt=F.NO_FUNC,
            plot_fmt=F.NO_FUNC,
            dtype=pd.DataFrame,
            doc="Coarse-grained transition matrix.",
        ),
        Metadata(attr=A.TERM_ABS_PROBS, prop=P.NO_PROPERTY, dtype=Lineage),
        Metadata(attr=A.COARSE_INIT_D, prop=P.COARSE_INIT_D, dtype=pd.Series),
        Metadata(attr=A.COARSE_STAT_D, prop=P.COARSE_STAT_D, dtype=pd.Series),
    ]

    def _read_from_adata(self) -> None:
        super()._read_from_adata()
        self._reconstruct_lineage(
            A.TERM_ABS_PROBS,
            self._term_abs_prob_key,
        )

    @inject_docs(
        ms=P.MACRO,
        msp=P.MACRO_MEMBER,
        schur=P.SCHUR.s,
        coarse_T=P.COARSE_T,
        coarse_stat=P.COARSE_STAT_D,
    )
    @d.dedent
    def compute_macrostates(
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
        Compute the macrostates.

        Parameters
        ----------
        n_states
            Number of macrostates. If `None`, use the `eigengap` heuristic.
        %(n_cells)s
        use_min_chi
            Whether to use :meth:`pygpcca.GPCCA.minChi` to calculate the number of macrostates.
            If `True`, ``n_states`` corresponds to a closed interval `[min, max]` inside of which the potentially
            optimal number of macrostates is searched.
        cluster_key
            If a key to cluster labels is given, names and colors of the states will be associated with the clusters.
        %(en_cutoff_p_thresh)s

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

        if use_min_chi:
            n_states = self._get_n_states_from_minchi(n_states)

        if n_states is None:
            if self._get(P.EIG) is None:
                raise RuntimeError(
                    "Compute eigendecomposition first as `.compute_eigendecomposition()` or `.compute_schur()`."
                )
            was_from_eigengap = True
            n_states = self._get(P.EIG)["eigengap"] + 1
            logg.info(f"Using `{n_states}` states based on eigengap")
        elif not isinstance(n_states, int):
            raise ValueError(
                f"Expected `n_states` to be an integer when `use_min_chi=False`, "
                f"found `{type(n_states).__name__!r}`."
            )

        if n_states <= 0:
            raise ValueError(
                f"Expected `n_states` to be positive or `None`, found `{n_states}`."
            )

        n_states = self._check_states_validity(n_states)
        if n_states == 1:
            self._compute_one_macrostate(
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

        # pre-computed X
        if self._gpcca._p_X.shape[1] < n_states:
            logg.warning(
                f"Requested more macrostates `{n_states}` than available "
                f"Schur vectors `{self._gpcca.X.shape[1]}`. Recomputing the decomposition"
            )

        start = logg.info(f"Computing `{n_states}` macrostates")
        try:
            self._gpcca = self._gpcca.optimize(m=n_states)
        except ValueError as e:
            # this is the following case - we have 4 Schur vectors, user requests 5 states, but it splits the conj. ev.
            # in the try block, Schur decomposition with 5 vectors is computed, but it fails (no way of knowing)
            # so in this case, we increase it by 1
            n_states += 1
            logg.warning(f"{e}\nIncreasing `n_states` to `{n_states}`")
            self._gpcca = self._gpcca.optimize(m=n_states)

        self._set_macrostates(
            memberships=self._gpcca.memberships,
            n_cells=n_cells,
            cluster_key=cluster_key,
            p_thresh=p_thresh,
            en_cutoff=en_cutoff,
        )

        # cache the results and make sure we don't overwrite
        self._set(A.SCHUR, self._gpcca._p_X)
        self._set(A.SCHUR_MAT, self._gpcca._p_R)

        names = self._get(P.MACRO_MEMBER).names

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
                    self._gpcca.coarse_grained_stationary_probability,
                    index=names,
                ),
            )
            logg.info(
                f"Adding `.{P.MACRO_MEMBER}`\n"
                f"       `.{P.MACRO}`\n"
                f"       `.{P.SCHUR}`\n"
                f"       `.{P.COARSE_T}`\n"
                f"       `.{P.COARSE_STAT_D}`\n"
                f"    Finish",
                time=start,
            )
        else:
            logg.warning("No stationary distribution found in GPCCA object")
            logg.info(
                f"Adding `.{P.MACRO_MEMBER}`\n"
                f"       `.{P.MACRO}`\n"
                f"       `.{P.SCHUR}`\n"
                f"       `.{P.COARSE_T}`\n"
                f"    Finish",
                time=start,
            )

    @d.dedent
    @inject_docs(fs=P.TERM, fsp=P.TERM_PROBS)
    def set_terminal_states_from_macrostates(
        self,
        names: Optional[Union[Sequence[str], Mapping[str, str], str]] = None,
        n_cells: int = 30,
    ):
        """
        Manually select terminal states from macrostates.

        Parameters
        ----------
        names
            Names of the macrostates to be marked as terminal. Multiple states can be combined using `','`,
            such as ``["Alpha, Beta", "Epsilon"]``.  If a :class:`dict`, keys correspond to the names
            of the macrostates and the values to the new names.  If `None`, select all macrostates.
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

        probs = self._get(P.MACRO_MEMBER)
        if probs is None:
            raise RuntimeError("Compute macrostates first as `.compute_macrostates()`.")

        rename = True
        if names is None:
            names = probs.names
            rename = False
        if isinstance(names, str):
            names = [names]
            rename = False
        if not isinstance(names, dict):
            names = {n: n for n in names}
            rename = False

        if not len(names):
            raise ValueError("No macrostates have been selected.")

        if not all(isinstance(old, str) for old in names.keys()):
            raise TypeError("Not all new names are strings.")

        if not all(isinstance(new, (str, int)) for new in names.values()):
            raise TypeError("Not all macrostates names are strings or integers.")

        # this also checks that the names are correct before renaming
        macrostates_probs = probs[list(names.keys())]

        # we do this also here because if `rename_terminal_states` fails
        # invalid states would've been written to this object and nothing to adata
        new_names = {k: str(v) for k, v in names.items()}
        # TODO: seems wrong
        names_after_renaming = [new_names.get(n, n) for n in probs.names]
        if len(set(names_after_renaming)) != probs.shape[1]:
            raise ValueError(
                f"After renaming, the names will not be unique: `{names_after_renaming}`."
            )

        if probs.shape[1] == 1:
            self._set(A.TERM, self._create_states(probs, n_cells=n_cells))
            self._set(A.TERM_COLORS, self._get(A.MACRO_COLORS))
            self._set(A.TERM_PROBS, probs / probs.max())
            self._set(A.TERM_ABS_PROBS, probs)
            if rename:
                # access lineage renames join states, e.g. 'Alpha, Beta' becomes 'Alpha or Beta' + whitespace stripping
                self.rename_terminal_states(
                    dict(zip(self._get(P.TERM).cat.categories, names.values()))
                )

            self._write_terminal_states()
            return

        # compute the aggregated probability of being a initial/terminal state (no matter which)
        scaled_probs = macrostates_probs.copy()
        scaled_probs /= scaled_probs.max(0)

        self._set(A.TERM, self._create_states(macrostates_probs, n_cells=n_cells))
        self._set(
            A.TERM_PROBS, pd.Series(scaled_probs.X.max(1), index=self.adata.obs_names)
        )
        self._set(
            A.TERM_COLORS,
            macrostates_probs[list(self._get(P.TERM).cat.categories)].colors,
        )
        self._set(A.TERM_ABS_PROBS, scaled_probs)
        if rename:
            self.rename_terminal_states(
                dict(zip(self._get(P.TERM).cat.categories, names.values()))
            )

        self._write_terminal_states()

    @inject_docs(fs=P.TERM, fsp=P.TERM_PROBS)
    @d.dedent
    def compute_terminal_states(
        self,
        method: str = "stability",
        n_cells: int = 30,
        alpha: Optional[float] = 1,
        stability_threshold: float = 0.96,
        n_states: Optional[int] = None,
    ):
        """
        Automatically select terminal states from macrostates.

        Parameters
        ----------
        method
            One of following:

                - `'eigengap'` - select the number of states based on the `eigengap` of the transition matrix.
                - `'eigengap_coarse'` - select the number of states based on the `eigengap` of the diagonal
                  of the coarse-grained transition matrix.
                - `'top_n'` - select top ``n_states`` based on the probability of the diagonal
                  of the coarse-grained transition matrix.
                - `'stability'` - select states which have a stability index >= ``stability_threshold``. The stability
                  index is given by the diagonal elements of the coarse-grained transition matrix.
        %(n_cells)s
        alpha
            Weight given to the deviation of an eigenvalue from one. Used when ``method='eigengap'``
            or ``method='eigengap_coarse'``.
        stability_threshold
            Threshold used when ``method='stability'``.
        n_states
            Numer of states used when ``method='top_n'``.

        Returns
        -------
        None
            Nothing, just updates the following fields:

                - :paramref:`{fsp}`
                - :paramref:`{fs}`
        """

        if len(self._get(P.MACRO).cat.categories) == 1:
            logg.warning("Found only one macrostate. Making it the single main state")
            self.set_terminal_states_from_macrostates(None, n_cells=n_cells)
            return

        coarse_T = self._get(P.COARSE_T)

        if method == "eigengap":
            if self._get(P.EIG) is None:
                raise RuntimeError(
                    "Compute eigendecomposition first as `.compute_eigendecomposition()`."
                )
            n_states = _eigengap(self._get(P.EIG)["D"], alpha=alpha) + 1
        elif method == "eigengap_coarse":
            if coarse_T is None:
                raise RuntimeError(
                    "Compute macrostates first as `.compute_macrostates()`."
                )
            n_states = _eigengap(np.sort(np.diag(coarse_T)[::-1]), alpha=alpha)
        elif method == "top_n":
            if n_states is None:
                raise ValueError(
                    "Argument `n_states` must be != `None` for `method='top_n'`."
                )
            elif n_states <= 0:
                raise ValueError(
                    f"Expected `n_states` to be positive, found `{n_states}`."
                )
        elif method == "stability":
            if stability_threshold is None:
                raise ValueError(
                    "Argument `stability_threshold` must be != `None` for `method='stability'`."
                )
            self_probs = pd.Series(np.diag(coarse_T), index=coarse_T.columns)
            names = self_probs[self_probs.values >= stability_threshold].index
            self.set_terminal_states_from_macrostates(names, n_cells=n_cells)
            return
        else:
            raise ValueError(
                f"Invalid method `{method!r}`. Valid options are `'eigengap', 'eigengap_coarse', "
                f"'top_n' and 'min_self_prob'`."
            )

        names = coarse_T.columns[np.argsort(np.diag(coarse_T))][-n_states:]
        self.set_terminal_states_from_macrostates(names, n_cells=n_cells)

    def compute_gdpt(
        self, n_components: int = 10, key_added: str = "gdpt_pseudotime", **kwargs
    ):
        """
        Compute generalized Diffusion pseudotime from [Haghverdi16]_ using the real Schur decomposition.

        Parameters
        ----------
        n_components
            Number of real Schur vectors to consider.
        key_added
            Key in :paramref:`adata` ``.obs`` where to save the pseudotime.
        kwargs
            Keyword arguments for :meth:`cellrank.tl.GPCCA.compute_schur` if Schur decomposition is not found.

        Returns
        -------
        None
            Nothing, just updates :paramref:`adata` ``.obs[key_added]`` with the computed pseudotime.
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
            f"Computing Generalized Diffusion Pseudotime using `n_components={n_components}`"
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
        Plot the coarse-grained transition matrix between macrostates.

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
        kwargs
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

            if xticks_labels is not None:
                ax.set_xticklabels(xticks_labels)
                ax.set_xticks(np.arange(data.shape[1]))
                plt.setp(
                    ax.get_xticklabels(),
                    rotation=xtick_rotation,
                    ha="right",
                    rotation_mode="anchor",
                )
            else:
                ax.set_xticks([])
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

        def annotate_dist_ax(ax, data: np.ndarray, valfmt: str = "{x:.2f}"):
            if ax is None:
                return

            if isinstance(valfmt, str):
                valfmt = mpl.ticker.StrMethodFormatter(valfmt)

            kw = {"ha": "center", "va": "center"}
            kw.update(**text_kwargs)

            for i, val in enumerate(data):
                kw.update(color=_get_black_or_white(im.norm(val), cmap))
                ax.text(
                    i,
                    0,
                    valfmt(val, None),
                    **kw,
                )

        coarse_T = self._get(P.COARSE_T)
        coarse_stat_d = self._get(P.COARSE_STAT_D)
        coarse_init_d = self._get(P.COARSE_INIT_D)

        if coarse_T is None:
            raise RuntimeError(
                "Compute coarse-grained transition matrix first as `.compute_macrostates()` with `n_states > 1`."
            )

        if show_stationary_dist and coarse_stat_d is None:
            logg.warning("Coarse stationary distribution is `None`, ignoring")
            show_stationary_dist = False
        if show_initial_dist and coarse_init_d is None:
            logg.warning("Coarse initial distribution is `None`, ignoring")
            show_initial_dist = False

        hrs, wrs = [1], [1]
        if show_stationary_dist:
            hrs += [0.05]
        if show_initial_dist:
            hrs += [0.05]
        if show_cbar:
            wrs += [0.025]

        dont_show_dist = not show_initial_dist and not show_stationary_dist

        fig = plt.figure(constrained_layout=False, figsize=figsize, dpi=dpi)
        gs = plt.GridSpec(
            1 + show_stationary_dist + show_initial_dist,
            1 + show_cbar,
            height_ratios=hrs,
            width_ratios=wrs,
            wspace=0.05,
            hspace=0.05,
        )
        if isinstance(cmap, str):
            cmap = plt.get_cmap(cmap)

        ax = fig.add_subplot(gs[0, 0])
        cax = fig.add_subplot(gs[:1, -1]) if show_cbar else None
        init_ax, stat_ax = None, None

        labels = list(self.coarse_T.columns)

        tmp = coarse_T
        if show_initial_dist:
            tmp = np.c_[tmp, coarse_stat_d]
        if show_initial_dist:
            tmp = np.c_[tmp, coarse_init_d]

        minn, maxx = np.nanmin(tmp), np.nanmax(tmp)
        norm = mpl.colors.Normalize(vmin=minn, vmax=maxx)

        if show_stationary_dist:
            stat_ax = fig.add_subplot(gs[1, 0])
            stylize_dist(
                stat_ax,
                np.array(coarse_stat_d).reshape(1, -1),
                xticks_labels=labels if not show_initial_dist else None,
            )
            stat_ax.yaxis.set_label_position("right")
            stat_ax.set_ylabel("stationary dist", rotation=0, ha="left", va="center")

        if show_initial_dist:
            init_ax = fig.add_subplot(gs[show_stationary_dist + show_initial_dist, 0])
            stylize_dist(
                init_ax, np.array(coarse_init_d).reshape(1, -1), xticks_labels=labels
            )

            init_ax.yaxis.set_label_position("right")
            init_ax.set_ylabel("initial dist", rotation=0, ha="left", va="center")

        im = ax.imshow(coarse_T, aspect="auto", cmap=cmap, norm=norm, **kwargs)
        ax.set_title("coarse-grained transition matrix" if title is None else title)

        if cax is not None:
            _ = mpl.colorbar.ColorbarBase(
                cax,
                cmap=cmap,
                norm=norm,
                ticks=np.linspace(minn, maxx, 10),
                format="%0.3f",
            )

        ax.set_yticks(np.arange(coarse_T.shape[0]))
        ax.set_yticklabels(labels)

        ax.tick_params(
            top=False,
            bottom=dont_show_dist,
            labeltop=False,
            labelbottom=dont_show_dist,
        )

        for spine in ax.spines.values():
            spine.set_visible(False)

        if dont_show_dist:
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
        ax.tick_params(which="minor", bottom=dont_show_dist, left=False, top=False)

        if annotate:
            annotate_heatmap(im)
            annotate_dist_ax(stat_ax, coarse_stat_d.values)
            annotate_dist_ax(init_ax, coarse_init_d)

        if save:
            save_fig(fig, save)

        fig.show()

    def _compute_one_macrostate(
        self,
        n_cells: int,
        cluster_key: Optional[str],
        en_cutoff: Optional[float],
        p_thresh: float,
    ) -> None:
        start = logg.warning("For 1 macrostate, stationary distribution is computed")

        eig = self._get(P.EIG)
        if (
            eig is not None
            and "stationary_dist" in eig
            and eig["params"]["which"] == "LR"
        ):
            stationary_dist = eig["stationary_dist"]
        else:
            self.compute_eigendecomposition(only_evals=False, which="LR")
            stationary_dist = self._get(P.EIG)["stationary_dist"]

        self._set_macrostates(
            memberships=stationary_dist[:, None],
            n_cells=n_cells,
            cluster_key=cluster_key,
            p_thresh=p_thresh,
            en_cutoff=en_cutoff,
        )
        self._set(
            A.MACRO_MEMBER,
            Lineage(
                stationary_dist,
                names=list(self._get(A.MACRO).cat.categories),
                colors=self._get(A.MACRO_COLORS),
            ),
        )

        # reset all the things
        for key in (
            A.ABS_PROBS,
            A.SCHUR,
            A.SCHUR_MAT,
            A.COARSE_T,
            A.COARSE_STAT_D,
            A.COARSE_STAT_D,
        ):
            self._set(key.s, None)

        logg.info(
            f"Adding `.{P.MACRO_MEMBER}`\n        `.{P.MACRO}`\n    Finish",
            time=start,
        )

    def _get_n_states_from_minchi(
        self, n_states: Union[Tuple[int, int], List[int], Dict[str, int]]
    ) -> int:
        if self._gpcca is None:
            raise RuntimeError(
                "Compute Schur decomposition first as `.compute_schur()` when `use_min_chi=True`."
            )

        if not isinstance(n_states, (dict, tuple, list)):
            raise TypeError(
                f"Expected `n_states` to be either `dict`, `tuple` or a `list`, "
                f"found `{type(n_states).__name__}`."
            )
        if len(n_states) != 2:
            raise ValueError(
                f"Expected `n_states` to be of size `2`, found `{len(n_states)}`."
            )

        if isinstance(n_states, dict):
            if "min" not in n_states or "max" not in n_states:
                raise KeyError(
                    f"Expected the dictionary to have `'min'` and `'max'` keys, "
                    f"found `{tuple(n_states.keys())}`."
                )
            minn, maxx = n_states["min"], n_states["max"]
        else:
            minn, maxx = n_states

        if minn > maxx:
            logg.debug(f"Swapping minimum and maximum because `{minn}` > `{maxx}`")
            minn, maxx = maxx, minn

        if minn <= 1:
            raise ValueError(f"Minimum value must be > `1`, found `{minn}`.")
        elif minn == 2:
            logg.warning(
                "In most cases, 2 clusters will always be optimal. "
                "If you really expect 2 clusters, use `n_states=2` and `use_minchi=False`. Setting minimum to `3`"
            )
            minn = 3

        if minn >= maxx:
            maxx = minn + 1
            logg.debug(
                f"Setting maximum to `{maxx}` because it was <= than minimum `{minn}`"
            )

        logg.info(f"Calculating minChi within interval `[{minn}, {maxx}]`")

        return int(np.arange(minn, maxx + 1)[np.argmax(self._gpcca.minChi(minn, maxx))])

    @d.dedent
    def _set_macrostates(
        self,
        memberships: np.ndarray,
        n_cells: Optional[int] = 30,
        cluster_key: str = "clusters",
        en_cutoff: Optional[float] = 0.7,
        p_thresh: float = 1e-15,
        check_row_sums: bool = True,
    ) -> None:
        """
        Map fuzzy clustering to pre-computed annotations to get names and colors.

        Given the fuzzy clustering, we would like to select the most likely cells from each state and use these to
        give each state a name and a color by comparing with pre-computed, categorical cluster annotations.

        Parameters
        ----------
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
            Only used to detect cell cycle stages. These have to be present in :paramref:`adata` ``.obs`` as
            `'G2M_score'` and `'S_score'`.
        check_row_sums
            Check whether rows in `memberships` sum to `1`.

        Returns
        -------
        None
            Writes a :class:`cellrank.tl.Lineage` object which mapped names and colors.
            Also writes a categorical :class:`pandas.Series`, where top ``n_cells`` cells represent each fuzzy state.
        """

        if n_cells is None:
            logg.debug("Setting the macrostates using macrostate assignment")

            max_assignment = np.argmax(memberships, axis=1)
            _macro_assignment = pd.Series(
                index=self.adata.obs_names, data=max_assignment, dtype="category"
            )
            # sometimes, the assignment can have a missing category and the Lineage creation therefore fails
            # keep it as ints when `n_cells != None`
            _macro_assignment.cat.set_categories(
                list(range(memberships.shape[1])), inplace=True
            )

            macrostates = _macro_assignment.astype(str).astype("category").copy()
            not_enough_cells = []
        else:
            logg.debug("Setting the macrostates using macrostates memberships")

            # select the most likely cells from each macrostate
            macrostates, not_enough_cells = self._create_states(
                memberships,
                n_cells=n_cells,
                check_row_sums=check_row_sums,
                return_not_enough_cells=True,
            )
            not_enough_cells = not_enough_cells.astype("str")

        # _set_categorical_labels creates the names, we still need to remap the group names
        orig_cats = macrostates.cat.categories
        self._set_categorical_labels(
            attr_key=A.MACRO.v,
            color_key=A.MACRO_COLORS.v,
            pretty_attr_key=P.MACRO.v,
            add_to_existing_error_msg="Compute macrostates first as `.compute_macrostates()`.",
            categories=macrostates,
            cluster_key=cluster_key,
            en_cutoff=en_cutoff,
            p_thresh=p_thresh,
            add_to_existing=False,
        )

        name_mapper = dict(zip(orig_cats, self._get(P.MACRO).cat.categories))
        _print_insufficient_number_of_cells(
            [name_mapper.get(n, n) for n in not_enough_cells], n_cells
        )
        logg.debug("Setting macrostates memberships based on GPCCA membership vectors")

        self._set(
            A.MACRO_MEMBER,
            Lineage(
                memberships,
                names=list(macrostates.cat.categories),
                colors=self._get(A.MACRO_COLORS),
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

    def _check_states_validity(self, n_states: int) -> int:
        if self._invalid_n_states is not None and n_states in self._invalid_n_states:
            logg.warning(
                f"Unable to compute macrostates with `n_states={n_states}` because it will "
                f"split the conjugate eigenvalues. Increasing `n_states` to `{n_states + 1}`"
            )
            n_states += 1  # cannot force recomputation of the Schur decomposition
            assert n_states not in self._invalid_n_states, "Sanity check failed."

        return n_states

    def _fit_terminal_states(
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
            self.compute_macrostates(
                n_states=n_lineages, cluster_key=cluster_key, **kwargs
            )
        except ValueError:
            logg.warning(
                f"Computing `{n_lineages}` macrostates cuts through a block of complex conjugates. "
                f"Increasing `n_lineages` to {n_lineages + 1}"
            )
            self.compute_macrostates(
                n_states=n_lineages + 1, cluster_key=cluster_key, **kwargs
            )

        fs_kwargs = {"n_cells": kwargs["n_cells"]} if "n_cells" in kwargs else {}

        if n_lineages is None:
            self.compute_terminal_states(method="eigengap", **fs_kwargs)
        else:
            self.set_terminal_states_from_macrostates(**fs_kwargs)

    @d.dedent  # because of fit
    @d.dedent
    @inject_docs(
        ms=P.MACRO,
        msp=P.MACRO_MEMBER,
        fs=P.TERM,
        fsp=P.TERM_PROBS,
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
        Run the pipeline, computing the macrostates, %(initial_or_terminal)s states \
        and optionally the absorption probabilities.

        It is equivalent to running::

            if n_lineages is None or n_lineages == 1:
                compute_eigendecomposition(...)  # get the stationary distribution
            if n_lineages > 1:
                compute_schur(...)

            compute_macrostates(...)

            if n_lineages is None:
                compute_terminal_states(...)
            else:
                set_terminal_states_from_macrostates(...)

            if compute_absorption_probabilities:
                compute_absorption_probabilities(...)

        Parameters
        ----------
        %(fit)s
        method
            Method to use when computing the Schur decomposition. Valid options are: `'krylov'` or `'brandts'`.
        compute_absorption_probabilities
            Whether to compute the absorption probabilities or only the %(initial_or_terminal)s states.
        kwargs
            Keyword arguments for :meth:`cellrank.tl.estimators.GPCCA.compute_macrostates`.

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

    @d.dedent
    def _compute_initial_states(self, n_states: int = 1, n_cells: int = 30) -> None:
        """
        Compute initial states from macrostates.

        Parameters
        ----------
        n_states
            Number of initial states.
        %(n_cells)s

        Returns
        -------
        %(set_initial_states_from_macrostates.returns)s
        """

        if n_states <= 0:
            raise ValueError(f"Expected `n_states` to be positive, found `{n_states}`.")

        if n_cells <= 0:
            raise ValueError(f"Expected `n_cells` to be positive, found `{n_cells}`.")

        probs = self._get(P.MACRO_MEMBER)

        if probs is None:
            raise RuntimeError("Compute macrostates first as `.compute_macrostates()`.")

        if n_states > probs.shape[1]:
            raise ValueError(
                f"Requested `{n_states}` initial states, but only `{probs.shape[1]}` macrostates have been computed."
            )

        if probs.shape[1] == 1:
            self._set_initial_states_from_macrostates(n_cells=n_cells)
            return

        stat_dist = self._get(P.COARSE_STAT_D)
        if stat_dist is None:
            raise RuntimeError("No coarse-grained stationary distribution found.")

        self._set_initial_states_from_macrostates(
            stat_dist[np.argsort(stat_dist)][:n_states].index, n_cells=n_cells
        )

    @d.get_sections(base="set_initial_states_from_macrostates", sections=["Returns"])
    @d.dedent
    @inject_docs(
        key=TermStatesKey.BACKWARD.s, probs_key=_probs(TermStatesKey.BACKWARD.s)
    )
    def _set_initial_states_from_macrostates(
        self,
        names: Optional[Union[Iterable[str], str]] = None,
        n_cells: int = 30,
    ) -> None:
        """
        Manually select initial states from macrostates.

        Note that no check is performed to ensure initial and terminal states are distinct.

        Parameters
        ----------
        names
            Names of the macrostates to be marked as initial states. Multiple states can be combined using `','`,
            such as `["Alpha, Beta", "Epsilon"]`.
        %(n_cells)s

        Returns
        -------
        None
            Nothing, just writes to :paramref:`adata`:

                - ``.obs[{key!r}]`` - probability of being an initial state.
                - ``.obs[{probs_key!r}]`` - top ``n_cells`` from each initial state.
        """

        if not isinstance(n_cells, int):
            raise TypeError(
                f"Expected `n_cells` to be of type `int`, found `{type(n_cells).__name__!r}`."
            )

        if n_cells <= 0:
            raise ValueError(f"Expected `n_cells` to be positive, found `{n_cells}`.")

        probs = self._get(P.MACRO_MEMBER)

        if probs is None:
            raise RuntimeError("Compute macrostates first as `.compute_macrostates()`.")
        elif probs.shape[1] == 1:
            categorical = self._create_states(probs, n_cells=n_cells)
            scaled = probs / probs.max()
        else:
            if names is None:
                names = probs.names
            if isinstance(names, str):
                names = [names]

            probs = probs[list(names)]
            categorical = self._create_states(probs, n_cells=n_cells)
            probs /= probs.max(0)

            # compute the aggregated probability of being a initial/terminal state (no matter which)
            scaled = probs.X.max(1)

        self._write_initial_states(membership=probs, probs=scaled, cats=categorical)

    def _write_initial_states(
        self, membership: Lineage, probs: pd.Series, cats: pd.Series, time=None
    ) -> None:
        key = TermStatesKey.BACKWARD.s

        self.adata.obs[key] = cats
        self.adata.obs[_probs(key)] = probs

        self.adata.uns[_colors(key)] = membership.colors
        self.adata.uns[_lin_names(key)] = membership.names

        logg.info(
            f"Adding `adata.obs[{_probs(key)!r}]`\n       `adata.obs[{key!r}]`\n",
            time=time,
        )
