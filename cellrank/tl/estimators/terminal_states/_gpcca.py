from typing import Any, Dict, Tuple, Union, Mapping, Optional, Sequence
from typing_extensions import Literal

from enum import auto
from types import MappingProxyType
from pathlib import Path
from datetime import datetime

from anndata import AnnData
from cellrank import logging as logg
from cellrank._key import Key
from cellrank.tl._enum import ModeEnum
from cellrank.ul._docs import d
from cellrank.tl._utils import (
    save_fig,
    _eigengap,
    _fuzzy_to_discrete,
    _series_from_one_hot_matrix,
)
from cellrank.tl._colors import _get_black_or_white, _create_categorical_colors
from cellrank.tl._lineage import Lineage
from cellrank.tl.estimators._utils import SafeGetter
from cellrank.tl.estimators.mixins import EigenMixin, SchurMixin, LinDriversMixin
from cellrank.tl.kernels._base_kernel import KernelExpression
from cellrank.tl.estimators.mixins._utils import logger, shadow, register_plotter
from cellrank.tl.estimators.terminal_states._term_states_estimator import (
    TermStatesEstimator,
)

import numpy as np
import pandas as pd
from scipy.sparse import spmatrix
from pandas.api.types import infer_dtype, is_categorical_dtype

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.colors import Normalize, ListedColormap
from matplotlib.ticker import StrMethodFormatter
from matplotlib.colorbar import ColorbarBase


class TermStatesMethod(ModeEnum):  # noqa: D101
    EIGENGAP = auto()
    EIGENGAP_COARSE = auto()
    TOP_N = auto()
    STABILITY = auto()


@d.dedent
class GPCCA(TermStatesEstimator, LinDriversMixin, SchurMixin, EigenMixin):
    """
    Generalized Perron Cluster Cluster Analysis :cite:`reuter:18` as implemented in \
    `pyGPCCA <https://pygpcca.readthedocs.io/en/latest/>`_.

    Coarse-grains a discrete Markov chain into a set of macrostates and computes coarse-grained transition probabilities
    among the macrostates. Each macrostate corresponds to an area of the state space, i.e. to a subset of cells. The
    assignment is soft, i.e. each cell is assigned to every macrostate with a certain weight, where weights sum to
    one per cell. Macrostates are computed by maximizing the 'crispness' which can be thought of as a measure for
    minimal overlap between macrostates in a certain inner-product sense. Once the macrostates have been computed,
    we project the large transition matrix onto a coarse-grained transition matrix among the macrostates via
    a Galerkin projection. This projection is based on invariant subspaces of the original transition matrix which
    are obtained using the real Schur decomposition :cite:`reuter:18`.

    Parameters
    ----------
    %(base_estimator.parameters)s
    """

    def __init__(
        self,
        obj: Union[AnnData, np.ndarray, spmatrix, KernelExpression],
        obsp_key: Optional[str] = None,
        **kwargs: Any,
    ):
        super().__init__(obj=obj, obsp_key=obsp_key, **kwargs)

        self._coarse_init_dist: Optional[pd.Series] = None
        self._coarse_stat_dist: Optional[pd.Series] = None
        self._coarse_tmat: Optional[pd.DataFrame] = None

        self._macrostates: Optional[pd.Series] = None
        self._macrostates_memberships: Optional[Lineage] = None
        self._macrostates_colors: Optional[np.ndarray] = None

        self._term_states_memberships: Optional[Lineage] = None

    @property
    @d.get_summary(base="gpcca_macro")
    def macrostates(self) -> Optional[pd.Series]:
        """Macrostates of the transition matrix."""
        return self._macrostates

    @property
    @d.get_summary(base="gpcca_macro_memberships")
    def macrostates_memberships(self) -> Optional[Lineage]:
        """Macrostate membership matrix.

        Soft assignment of microstates (cells) to macrostates.
        """
        return self._macrostates_memberships

    @property
    @d.get_summary(base="gpcca_term_states_memberships")
    def terminal_states_memberships(self) -> Optional[Lineage]:
        """Terminal state membership matrix.

        Soft assignment of cells to terminal states.
        """
        return self._term_states_memberships

    @property
    @d.get_summary(base="gpcca_coarse_init")
    def coarse_initial_distribution(self) -> Optional[pd.Series]:
        """Coarse-grained initial distribution."""
        return self._coarse_init_dist

    @property
    @d.get_summary(base="gpcca_coarse_stat")
    def coarse_stationary_distribution(self) -> Optional[pd.Series]:
        """Coarse-grained stationary distribution."""
        return self._coarse_stat_dist

    @property
    @d.get_summary(base="gpcca_coarse_tmat")
    def coarse_T(self) -> Optional[pd.DataFrame]:
        """Coarse-grained transition matrix."""
        return self._coarse_tmat

    @d.get_sections(base="gpcca_compute_macro", sections=["Parameters", "Returns"])
    @d.dedent
    def compute_macrostates(
        self,
        n_states: Optional[Union[int, Sequence[int]]] = None,
        n_cells: Optional[int] = 30,
        cluster_key: Optional[str] = None,
        **kwargs: Any,
    ) -> None:
        """
        Compute the macrostates.

        Parameters
        ----------
        n_states
            Number of macrostates. If a :class:`typing.Sequence`, use the *minChi* criterion :cite:`reuter:18`.
            If `None`, use the *eigengap* heuristic.
        %(n_cells)s
        cluster_key
            If a key to cluster labels is given, names and colors of the states will be associated with the clusters.
        kwargs
            Keyword arguments for :meth:`compute_schur`.

        Returns
        -------
        Nothing, just updates the following fields:

            - :attr:`macrostates` - %(gpcca_macro.summary)s
            - :attr:`macrostates_memberships` - %(gpcca_macro_memberships.summary)s
            - :attr:`coarse_T` - %(gpcca_coarse_tmat.summary)s
            - :attr:`coarse_initial_distribution` - %(gpcca_coarse_init.summary)s
            - :attr:`coarse_stationary_distribution` - %(gpcca_coarse_stat.summary)s
            - :attr:`schur_vectors` - %(schur_vectors.summary)s
            - :attr:`schur_matrix` - %(schur_matrix.summary)s
            - :attr:`eigendecomposition` - %(eigen.summary)s
        """

        n_states = self._n_states(n_states)
        if n_states == 1:
            self._compute_one_macrostate(
                n_cells=n_cells,
                cluster_key=cluster_key,
            )
            return

        if self._gpcca is None or kwargs:
            self.compute_schur(n_states, **kwargs)
        n_states = self._validate_n_states(n_states)

        if self._gpcca._p_X.shape[1] < n_states:
            # precomputed X
            logg.warning(
                f"Requested more macrostates `{n_states}` than available "
                f"Schur vectors `{self._gpcca._p_X.shape[1]}`. Recomputing the decomposition"
            )

        start = logg.info(f"Computing `{n_states}` macrostates")
        try:
            self._gpcca = self._gpcca.optimize(m=n_states)
        except ValueError as e:
            if "will split complex conjugate eigenvalues" not in str(e):
                raise
            # this is the following case - we have 4 Schur vectors, user requests 5 states, but it splits the conj. ev.
            # in the try block, Schur decomposition with 5 vectors is computed, but it fails (no way of knowing)
            # so in this case, we increase it by 1
            logg.warning(
                f"Unable to compute macrostates with `n_states={n_states}` because it will "
                f"split complex conjugate eigenvalues. Using `n_states={n_states + 1}`"
            )
            self._gpcca = self._gpcca.optimize(m=n_states + 1)

        self._set_macrostates(
            memberships=self._gpcca.memberships,
            n_cells=n_cells,
            cluster_key=cluster_key,
            params=self._create_params(),
            time=start,
        )

    @d.dedent
    def predict(
        self,
        method: Literal[
            "stability", "top_n", "eigengap", "eigengap_coarse"
        ] = TermStatesMethod.STABILITY,
        n_cells: int = 30,
        alpha: Optional[float] = 1,
        stability_threshold: float = 0.96,
        n_states: Optional[int] = None,
    ) -> None:
        """
        Automatically select terminal states from macrostates.

        Parameters
        ----------
        method
            How to select the terminal states. Valid option are:

                - `'eigengap'` - select the number of states based on the *eigengap* of :attr:`transition_matrix`.
                - `'eigengap_coarse'` - select the number of states based on the *eigengap* of the diagonal
                  of :attr:`coarse_T`.
                - `'top_n'` - select top ``n_states`` based on the probability of the diagonal of :attr:`coarse_T`.
                - `'stability'` - select states which have a stability >= ``stability_threshold``.
                  The stability is given by the diagonal elements of :attr:`coarse_T`.
        %(n_cells)s
        alpha
            Weight given to the deviation of an eigenvalue from one.
            Only used when ``method = 'eigengap'`` or ``method = 'eigengap_coarse'``.
        stability_threshold
            Threshold used when ``method = 'stability'``.
        n_states
            Number of states used when ``method = 'top_n'``.

        Returns
        -------
        Nothing, just updates the following fields:

            - :attr:`terminal_states` - %(tse_term_states.summary)s
            - :attr:`terminal_states_memberships` - %(gpcca_term_states_memberships.summary)s
            - :attr:`terminal_states_probabilities` - %(tse_term_states_probs.summary)s
        """
        if self.macrostates is None:
            raise RuntimeError("Compute macrostates first as `.compute_macrostates()`.")

        # fmt: off
        if len(self._macrostates.cat.categories) == 1:
            logg.warning("Found only one macrostate. Making it the single terminal state")
            self.set_terminal_states_from_macrostates(None, n_cells=n_cells, params=self._create_params())
            return

        method = TermStatesMethod(method)
        eig = self.eigendecomposition
        coarse_T = self.coarse_T

        if method == TermStatesMethod.EIGENGAP:
            if eig is None:
                raise RuntimeError("Compute eigendecomposition first as `.compute_eigendecomposition()`.")
            n_states = _eigengap(eig["D"], alpha=alpha) + 1
        elif method == TermStatesMethod.EIGENGAP_COARSE:
            if coarse_T is None:
                raise RuntimeError("Compute macrostates first as `.compute_macrostates()`.")
            n_states = _eigengap(np.sort(np.diag(coarse_T)[::-1]), alpha=alpha)
        elif method == TermStatesMethod.TOP_N:
            if n_states is None:
                raise ValueError("Expected `n_states != None` for `method='top_n'`.")
            elif n_states <= 0:
                raise ValueError(f"Expected `n_states` to be positive, found `{n_states}`.")
        elif method == TermStatesMethod.STABILITY:
            if stability_threshold is None:
                raise ValueError("Expected `stability_threshold != None` for `method='stability'`.")
            stability = pd.Series(np.diag(coarse_T), index=coarse_T.columns)
            names = stability[stability.values >= stability_threshold].index
            self.set_terminal_states_from_macrostates(names, n_cells=n_cells, params=self._create_params())
            return
        else:
            raise NotImplementedError(f"Method `{method}` is not yet implemented.")
        # fmt: on

        names = coarse_T.columns[np.argsort(np.diag(coarse_T))][-n_states:]
        self.set_terminal_states_from_macrostates(
            names, n_cells=n_cells, params=self._create_params()
        )

        return

    @d.dedent
    def set_terminal_states_from_macrostates(
        self,
        names: Optional[Union[str, Sequence[str], Mapping[str, str]]] = None,
        n_cells: int = 30,
        **kwargs: Any,
    ) -> None:
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
        Nothing, just updates the following fields:

            - :attr:`terminal_states` - %(tse_term_states.summary)s
            - :attr:`terminal_states_probabilities` - %(tse_term_states_probs.summary)s
            - :attr:`terminal_states_probabilities_memberships` - %(gpcca_term_states_memberships.summary)s
        """
        if n_cells <= 0:
            raise ValueError(f"Expected `n_cells` to be positive, found `{n_cells}`.")

        memberships = self.macrostates_memberships
        if memberships is None:
            raise RuntimeError("Compute macrostates first as `.compute_macrostates()`.")

        rename = True
        if names is None:
            names = memberships.names
            rename = False
        if isinstance(names, str):
            names = [names]
            rename = False
        if not isinstance(names, dict):
            names = {n: n for n in names}
            rename = False
        if not len(names):
            raise ValueError("No macrostates have been selected.")

        # we do this also here because if `rename_terminal_states` fails
        # invalid states would've been written to this object and nothing to adata
        names = {str(k): str(v) for k, v in names.items()}
        names_after_renaming = {names.get(n, n) for n in memberships.names}
        if len(names_after_renaming) != memberships.shape[1]:
            raise ValueError(
                f"After renaming, terminal state names will no longer be unique: `{names_after_renaming}`."
            )

        # this also checks that the names are correct before renaming
        is_singleton = memberships.shape[1] == 1
        memberships = memberships[list(names.keys())].copy()

        states = self._create_states(memberships, n_cells=n_cells, check_row_sums=False)
        if is_singleton:
            colors = self._macrostates_colors.copy()
            probs = memberships.X.squeeze() / memberships.X.max()
        else:
            colors = memberships[list(states.cat.categories)].colors
            probs = (memberships.X / memberships.X.max(0)).max(1)
        probs = pd.Series(probs, index=self.adata.obs_names)

        self._write_terminal_states(
            states, colors, probs, memberships, params=kwargs.pop("params", {})
        )
        if rename:
            # TODO(michalk8): in a future PR, remove this behavior in Lineage
            # access lineage renames join states, e.g. 'Alpha, Beta' becomes 'Alpha or Beta' + whitespace stripping
            self.rename_terminal_states(
                dict(zip(self.terminal_states.cat.categories, names.values()))
            )

    @d.dedent
    def rename_terminal_states(self, new_names: Mapping[str, str]) -> None:
        """
        %(tse_rename_term_states.full_desc)s

        Parameters
        ----------
        %(tse_rename_term_states.parameters)s

        Returns
        -------
        %(tse_rename_term_states.returns)s
            - :attr:`terminal_states_memberships` - %(gpcca_term_states_memberships.summary)s
        """  # noqa: D400
        term_states_memberships = self.terminal_states_memberships
        super().rename_terminal_states(new_names)

        # fmt: off
        new_names = {str(k): str(v) for k, v in new_names.items()}
        term_states_memberships.names = [new_names.get(n, n) for n in term_states_memberships.names]
        self._set("_term_states_memberships", value=term_states_memberships, shadow_only=True)
        # fmt: on

        with self._shadow:
            key = Key.obsm.memberships(Key.obs.macrostates(self.backward))
            self._set(obj=self.adata.obsm, key=key, value=term_states_memberships)

    @d.dedent
    def fit(
        self,
        n_states: Optional[Union[int, Sequence[int]]] = None,
        n_cells: Optional[int] = 30,
        cluster_key: Optional[str] = None,
        **kwargs: Any,
    ) -> "GPCCA":
        """
        Prepare self for terminal states prediction.

        Parameters
        ----------
        %(gpcca_compute_macro.parameters)s

        Returns
        -------
        %(gpcca_compute_macro.returns)s
        """
        if n_states is None:
            self.compute_eigendecomposition()
            n_states = self.eigendecomposition["eigengap"] + 1
        if isinstance(n_states, int) and n_states == 1:
            self.compute_eigendecomposition()

        self.compute_macrostates(n_states=n_states, cluster_key=cluster_key, **kwargs)

        return self

    @d.dedent
    def plot_coarse_T(
        self,
        show_stationary_dist: bool = True,
        show_initial_dist: bool = False,
        cmap: Union[str, ListedColormap] = "viridis",
        xtick_rotation: float = 45,
        annotate: bool = True,
        show_cbar: bool = True,
        title: Optional[str] = None,
        figsize: Tuple[float, float] = (8, 8),
        dpi: int = 80,
        save: Optional[Union[Path, str]] = None,
        text_kwargs: Mapping[str, Any] = MappingProxyType({}),
        **kwargs: Any,
    ) -> None:
        """
        Plot the coarse-grained transition matrix between macrostates.

        Parameters
        ----------
        show_stationary_dist
            Whether to show :attr:`coarse_stationary_distribution`, if present.
        show_initial_dist
            Whether to show :attr:`coarse_initial_distribution`.
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
            ax: Axes, data: np.ndarray, xticks_labels: Sequence[str] = ()
        ) -> None:
            _ = ax.imshow(data, aspect="auto", cmap=cmap, norm=norm)
            for spine in ax.spines.values():
                spine.set_visible(False)

            if xticks_labels is not None:
                ax.set_xticks(np.arange(data.shape[1]))
                ax.set_xticklabels(xticks_labels)
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

        def annotate_heatmap(im, valfmt: str = "{x:.2f}") -> None:
            # modified from matplotlib's site

            data = im.get_array()
            kw = {"ha": "center", "va": "center"}
            kw.update(**text_kwargs)

            # Get the formatter in case a string is supplied
            if isinstance(valfmt, str):
                valfmt = StrMethodFormatter(valfmt)

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
                valfmt = StrMethodFormatter(valfmt)

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

        coarse_T = self.coarse_T
        coarse_init_d = self.coarse_initial_distribution
        coarse_stat_d = self.coarse_stationary_distribution

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
        norm = Normalize(vmin=minn, vmax=maxx)

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
            _ = ColorbarBase(
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
            if show_stationary_dist:
                annotate_dist_ax(stat_ax, coarse_stat_d.values)
            if show_initial_dist:
                annotate_dist_ax(init_ax, coarse_init_d.values)

        if save:
            save_fig(fig, save)

    @d.dedent
    def plot_macrostate_composition(
        self,
        key: str,
        width: float = 0.8,
        title: Optional[str] = None,
        labelrot: float = 45,
        legend_loc: Optional[str] = "upper right out",
        figsize: Optional[Tuple[float, float]] = None,
        dpi: Optional[int] = None,
        save: Optional[Union[str, Path]] = None,
        show: bool = True,
    ) -> Optional[Axes]:
        """
        Plot stacked histogram of macrostates over categorical annotations.

        Parameters
        ----------
        %(adata)s
        key
            Key from :attr:`anndata.AnnData.obs` containing categorical annotations.
        width
            Bar width in `[0, 1]`.
        title
            Title of the figure. If `None`, create one automatically.
        labelrot
            Rotation of labels on x-axis.
        legend_loc
            Position of the legend. If `None`, don't show legend.
        %(plotting)s
        show
            If `False`, return :class:`matplotlib.pyplot.Axes`.

        Returns
        -------
        The axes object, if ``show = False``.
        %(just_plots)s
        """
        from cellrank.pl._utils import _position_legend

        macrostates = self.macrostates
        if macrostates is None:
            raise RuntimeError("Compute macrostates first as `.compute_macrostates()`.")
        if key not in self.adata.obs:
            raise KeyError(f"Data not found in `adata.obs[{key!r}]`.")
        if not is_categorical_dtype(self.adata.obs[key]):
            raise TypeError(
                f"Expected `adata.obs[{key!r}]` to be `categorical`, "
                f"found `{infer_dtype(self.adata.obs[key])}`."
            )

        mask = ~macrostates.isnull()
        df = (
            pd.DataFrame({"macrostates": macrostates, key: self.adata.obs[key]})[mask]
            .groupby([key, "macrostates"])
            .size()
        )
        try:
            cats_colors = self.adata.uns[f"{key}_colors"]
        except KeyError:
            cats_colors = _create_categorical_colors(
                len(self.adata.obs[key].cat.categories)
            )
        cat_color_mapper = dict(zip(self.adata.obs[key].cat.categories, cats_colors))
        x_indices = np.arange(len(macrostates.cat.categories))
        bottom = np.zeros_like(x_indices, dtype=np.float32)

        width = min(1, max(0, width))
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi, tight_layout=True)
        for cat, color in cat_color_mapper.items():
            frequencies = df.loc[cat]
            # do not add to legend if category is missing
            if np.sum(frequencies) > 0:
                ax.bar(
                    x_indices,
                    frequencies,
                    width,
                    label=cat,
                    color=color,
                    bottom=bottom,
                    ec="black",
                    lw=0.5,
                )
                bottom += np.array(frequencies)

        ax.set_xticks(x_indices)
        ax.set_xticklabels(
            # assuming at least 1 category
            frequencies.index,
            rotation=labelrot,
            ha="center" if labelrot in (0, 90) else "right",
        )
        y_max = bottom.max()
        ax.set_ylim([0, y_max + 0.05 * y_max])
        ax.set_yticks(np.linspace(0, y_max, 5))
        ax.margins(0.05)

        ax.set_xlabel("macrostate")
        ax.set_ylabel("frequency")
        if title is None:
            title = f"distribution over {key}"
        ax.set_title(title)
        if legend_loc not in (None, "none"):
            _position_legend(ax, legend_loc=legend_loc)

        if save is not None:
            save_fig(fig, save)

        if not show:
            return ax

    def _n_states(self, n_states: Optional[Union[int, Sequence[int]]]) -> int:
        if n_states is None:
            if self.eigendecomposition is None:
                raise RuntimeError(
                    "Compute eigendecomposition first as `.compute_eigendecomposition()` or "
                    "supply `n_states != None`."
                )
            return self.eigendecomposition["eigengap"] + 1

        # fmt: off
        if isinstance(n_states, int):
            if n_states <= 0:
                raise ValueError(f"Expected `n_states` to be positive, found `{n_states}`.")
            return n_states

        if self._gpcca is None:
            raise RuntimeError("Compute Schur decomposition first as `.compute_schur()`.")

        if not isinstance(n_states, Sequence):
            raise TypeError(f"Expected `n_states` to be a `Sequence`, found `{type(n_states).__name__!r}`.")
        if len(n_states) != 2:
            raise ValueError(f"Expected `n_states` to be of size `2`, found `{len(n_states)}`.")

        minn, maxx = sorted(n_states)
        if minn <= 1:
            minn = 2
            logg.warning(f"Minimum value must be larger than `1`, found `{minn}`. Setting `min={minn}`")
        if minn == 2:
            minn = 3
            logg.warning(
                f"In most cases, 2 clusters will always be optimal. "
                f"If you really expect 2 clusters, use `n_states=2`. Setting `min={minn}`"
            )
        # fmt: on
        maxx = max(minn + 1, maxx)

        logg.info(f"Calculating minChi criterion in interval `[{minn}, {maxx}]`")

        return int(np.arange(minn, maxx + 1)[np.argmax(self._gpcca.minChi(minn, maxx))])

    def _create_states(
        self,
        probs: Union[np.ndarray, Lineage],
        n_cells: int,
        check_row_sums: bool = False,
        return_not_enough_cells: bool = False,
    ) -> Union[pd.Series, Tuple[pd.Series, np.ndarray]]:
        if n_cells <= 0:
            raise ValueError(f"Expected `n_cells` to be positive, found `{n_cells}`.")

        discrete, not_enough_cells = _fuzzy_to_discrete(
            a_fuzzy=probs,
            n_most_likely=n_cells,
            remove_overlap=False,
            raise_threshold=0.2,
            check_row_sums=check_row_sums,
        )

        states = _series_from_one_hot_matrix(
            membership=discrete,
            index=self.adata.obs_names,
            names=probs.names if isinstance(probs, Lineage) else None,
        )

        return (states, not_enough_cells) if return_not_enough_cells else states

    def _validate_n_states(self, n_states: int) -> int:
        if self._invalid_n_states is not None and n_states in self._invalid_n_states:
            logg.warning(
                f"Unable to compute macrostates with `n_states={n_states}` because it will "
                f"split complex conjugate eigenvalues. Using `n_states={n_states + 1}`"
            )
            n_states += 1  # cannot force recomputation of the Schur decomposition
            assert n_states not in self._invalid_n_states, "Sanity check failed."

        return n_states

    def _compute_one_macrostate(
        self,
        n_cells: Optional[int],
        cluster_key: Optional[str],
    ) -> None:
        start = logg.info("For 1 macrostate, stationary distribution is computed")

        eig = self.eigendecomposition
        if (
            eig is not None
            and "stationary_dist" in eig
            and eig["params"]["which"] == "LR"
        ):
            stationary_dist = eig["stationary_dist"]
        else:
            self.compute_eigendecomposition(only_evals=False, which="LR")
            stationary_dist = self.eigendecomposition["stationary_dist"]

        self._set_macrostates(
            memberships=stationary_dist[:, None],
            n_cells=n_cells,
            cluster_key=cluster_key,
            check_row_sums=False,
            time=start,
        )

    @d.dedent
    def _set_macrostates(
        self,
        memberships: np.ndarray,
        n_cells: Optional[int] = 30,
        cluster_key: str = "clusters",
        check_row_sums: bool = True,
        time: Optional[datetime] = None,
        params: Dict[str, Any] = MappingProxyType({}),
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
            Key from :attr:`anndata.AnnData.obs` to get reference cluster annotations.
        check_row_sums
            Check whether rows in `memberships` sum to `1`.
        time
            Start time of macrostates computation.
        params
            Parameters used in macrostates computation.

        Returns
        -------
        Nothing, just updates the field as described in :meth:`compute_macrostates`.
        """

        if n_cells is None:
            # fmt: off
            logg.debug("Setting the macrostates using macrostate assignment")
            assignment = pd.Series(np.argmax(memberships, axis=1).astype(str), dtype="category")
            # sometimes, a category can be missing
            assignment = assignment.cat.reorder_categories([str(i) for i in range(memberships.shape[1])])
            not_enough_cells = []
            # fmt: on
        else:
            logg.debug("Setting the macrostates using macrostates memberships")

            # select the most likely cells from each macrostate
            assignment, not_enough_cells = self._create_states(
                memberships,
                n_cells=n_cells,
                check_row_sums=check_row_sums,
                return_not_enough_cells=True,
            )

        # remove previous fields
        self._write_terminal_states(None, None, None, None, log=False)

        # fmt: off
        assignment, colors = self._set_categorical_labels(assignment, cluster_key=cluster_key)
        memberships = Lineage(memberships, names=list(assignment.cat.categories), colors=colors)
        # fmt: on

        groups = pd.DataFrame(assignment).groupby(0).size()
        groups = groups[groups != n_cells].to_dict()
        if len(groups):
            logg.warning(
                f"The following terminal states have different number "
                f"of cells than requested ({n_cells}): {groups}"
            )

        self._write_macrostates(
            assignment, colors, memberships, time=time, params=params
        )

    @logger
    @shadow
    def _write_macrostates(
        self,
        macrostates: pd.Series,
        colors: np.ndarray,
        memberships: Lineage,
        params: Dict[str, Any] = MappingProxyType({}),
    ) -> str:
        # fmt: off
        names = list(macrostates.cat.categories)

        key = Key.obs.macrostates(self.backward)
        self._set("_macrostates", obj=self.adata.obs, key=key, value=macrostates, shadow_only=True)
        ckey = Key.uns.colors(key)
        self._set("_macrostates_colors", obj=self.adata.uns, key=ckey, value=colors, shadow_only=True)
        mkey = Key.obsm.memberships(key)
        self._set("_macrostates_memberships", obj=self.adata.obsm, key=mkey, value=memberships, shadow_only=True)
        self.params[key] = dict(params)

        if len(names) > 1:
            # not using stationary distribution
            g = self._gpcca
            tmat = pd.DataFrame(g.coarse_grained_transition_matrix, index=names, columns=names)
            init_dist = pd.Series(g.coarse_grained_input_distribution, index=names)
            stat_dist = pd.Series(g.coarse_grained_stationary_probability, index=names)
            dists = pd.DataFrame({"coarse_init_dist": init_dist})
            if stat_dist is not None:
                dists["coarse_stat_dist"] = pd.Series(stat_dist, index=names)

            key = Key.obsm.schur_vectors(self.backward)
            self._set("_schur_vectors", obj=self.adata.obsm, key=key, value=g._p_X, shadow_only=True)
            key = Key.uns.schur_matrix(self.backward)
            self._set("_schur_matrix", obj=self.adata.uns, key=key, value=g._p_R, shadow_only=True)
            self._set("_coarse_tmat", value=tmat, shadow_only=True)
            self._set("_coarse_init_dist", value=init_dist, shadow_only=True)
            self._set("_coarse_stat_dist", value=stat_dist, shadow_only=True)
            self._set(obj=self.adata.uns, key=Key.uns.coarse(self.backward), value=AnnData(tmat, obs=dists))
        else:
            for attr in ["_schur_vectors", "_schur_matrix", "_coarse_tmat", "_coarse_init_dist", "_coarse_stat_dist"]:
                self._set(attr, value=None, shadow_only=True)
            self._set(obj=self.adata.uns, key=Key.uns.coarse(self.backward), value=None)
        # fmt: on

        return (
            "Adding `.macrostates`\n"
            "       `.macrostates_memberships`\n"
            "       `.coarse_T`\n"
            "       `.coarse_initial_distribution\n"
            "       `.coarse_stationary_distribution`\n"
            "       `.schur_vectors`\n"
            "       `.schur_matrix`\n"
            "       `.eigendecomposition`\n"
            "    Finish"
        )

    @logger
    @shadow
    def _write_terminal_states(
        self,
        states: Optional[pd.Series],
        colors: Optional[np.ndarray],
        probs: Optional[pd.Series] = None,
        memberships: Optional[Lineage] = None,
        params: Dict[str, Any] = MappingProxyType({}),
    ) -> str:
        # fmt: off
        msg = super()._write_terminal_states(states, colors, probs, params=params, log=False)
        msg = "\n".join(msg.split("\n")[:-1])
        msg += "\n       `.terminal_states_memberships\n    Finish`"

        self._write_absorption_probabilities(None, None, log=False)
        key = Key.obsm.memberships(Key.obs.term_states(self.backward))
        self._set("_term_states_memberships", obj=self.adata.obsm, key=key, value=memberships)
        # fmt: on

        return msg

    def _read_from_adata(self, adata: AnnData, **kwargs: Any) -> bool:
        _ = self._read_eigendecomposition(adata, allow_missing=True)
        # TODO(michalk8): reintroduce in 2.0
        ok = self._read_schur_decomposition(adata, allow_missing=True)
        if not ok:
            return False

        # fmt: off
        with SafeGetter(self, allowed=KeyError) as sg:
            key = Key.obs.macrostates(self.backward)
            # TODO(michalk8): in the future, be more stringent and ensure the categories match the macro memberships
            self._get("_macrostates", self.adata.obs, key=key, where="obs", dtype=pd.Series)
            ckey = Key.uns.colors(key)
            self._get("_macrostates_colors", self.adata.uns, key=ckey, where="uns", dtype=(list, tuple, np.ndarray))
            mkey = Key.obsm.memberships(key)
            self._get("_macrostates_memberships", self.adata.obsm, key=mkey, where="obsm", dtype=(Lineage, np.ndarray))
            self._ensure_lineage_object("_macrostates_memberships", kind="macrostates")

            self._macrostates_colors = self.macrostates_memberships.colors.copy()
            self.params[key] = self._read_params(key)

            # TODO(michalk8): allow missing?
            tmat: AnnData = self.adata.uns[Key.uns.coarse(self.backward)]
            if not isinstance(tmat, AnnData):
                raise TypeError(f"Expected coarse-grained transition matrix to be stored "
                                f"as `anndata.AnnData`, found `{type(tmat).__name__}`.")
            tmat = tmat.copy()
            names = tmat.obs_names

            self._coarse_tmat = pd.DataFrame(tmat.X, index=names, columns=names)
            self._coarse_init_dist = tmat.obs["coarse_init_dist"]
            self._coarse_stat_dist = tmat.obs.get("coarse_stat_dist", None)

            self._set(obj=self._shadow_adata.uns, key=Key.uns.coarse(self.backward), value=tmat)

        # TODO(michalk8): reintroduce this in 2.0 - this is done for high-level plotting of init/term states only
        # if not sg.ok:
        #    return False

        if not super()._read_from_adata(adata, **kwargs):
            return False

        with SafeGetter(self, allowed=KeyError) as sg:
            key = Key.obsm.memberships(Key.obs.term_states(self.backward))
            self._get("_term_states_memberships", self.adata.obsm, key=key, where="obsm", dtype=(np.ndarray, Lineage))
            self._ensure_lineage_object("_term_states_memberships", kind="term_states")
        # fmt: on

        return sg.ok and self._read_absorption_probabilities(adata)

    plot_macrostates = register_plotter(
        discrete="macrostates", continuous="macrostates_memberships"
    )
    plot_terminal_states = register_plotter(
        discrete="terminal_states", continuous="terminal_states_memberships"
    )

    @d.dedent
    def _compute_initial_states(self, n_states: int = 1, n_cells: int = 30) -> None:
        """
        Compute initial states from macrostates using :attr:`coarse_stationary_distribution`.

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

        probs = self.macrostates_memberships
        if probs is None:
            raise RuntimeError("Compute macrostates first as `.compute_macrostates()`.")

        if n_states > probs.shape[1]:
            raise ValueError(
                f"Requested `{n_states}` initial states, but only `{probs.shape[1]}` macrostates have been computed."
            )

        if probs.shape[1] == 1:
            self._set_initial_states_from_macrostates(n_cells=n_cells)
            return

        stat_dist = self.coarse_stationary_distribution
        if stat_dist is None:
            raise RuntimeError("No coarse-grained stationary distribution found.")

        self._set_initial_states_from_macrostates(
            stat_dist[np.argsort(stat_dist)][:n_states].index, n_cells=n_cells
        )

    @d.get_sections(base="set_initial_states_from_macrostates", sections=["Returns"])
    @d.dedent
    def _set_initial_states_from_macrostates(
        self,
        names: Optional[Union[str, Sequence[str]]] = None,
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
        Nothing, just modifies :attr:`anndata.AnnData.obs`. The actual keys depend of :attr:`backward`.
        """

        if not isinstance(n_cells, int):
            raise TypeError(
                f"Expected `n_cells` to be of type `int`, found `{type(n_cells).__name__}`."
            )

        if n_cells <= 0:
            raise ValueError(f"Expected `n_cells` to be positive, found `{n_cells}`.")

        probs = self.macrostates_memberships
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
        key = Key.obs.term_states(not self.backward)

        self.adata.obs[key] = cats
        self.adata.obs[Key.obs.probs(key)] = probs

        self.adata.uns[Key.uns.colors(key)] = membership.colors

        logg.info(
            f"Adding `adata.obs[{key!r}]`\n       `adata.obs[{Key.uns.colors(key)!r}]`\n",
            time=time,
        )
