import collections
import datetime
import enum
import pathlib
import types
from pathlib import Path
from typing import Any, Dict, List, Literal, Mapping, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
import scipy.sparse as sp
from pandas.api.types import infer_dtype

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.axes import Axes
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.ticker import StrMethodFormatter

from anndata import AnnData

from cellrank import logging as logg
from cellrank._utils._colors import _create_categorical_colors, _get_black_or_white
from cellrank._utils._docs import d, inject_docs
from cellrank._utils._enum import ModeEnum
from cellrank._utils._key import Key
from cellrank._utils._lineage import Lineage
from cellrank._utils._utils import (
    _eigengap,
    _fuzzy_to_discrete,
    _series_from_one_hot_matrix,
    save_fig,
)
from cellrank.estimators.mixins import EigenMixin, LinDriversMixin, SchurMixin
from cellrank.estimators.mixins._utils import SafeGetter, StatesHolder, logger, shadow
from cellrank.estimators.terminal_states._term_states_estimator import (
    TermStatesEstimator,
)
from cellrank.kernels._base_kernel import KernelExpression

__all__ = ["GPCCA"]


class TermStatesMethod(ModeEnum):
    EIGENGAP = enum.auto()
    EIGENGAP_COARSE = enum.auto()
    TOP_N = enum.auto()
    STABILITY = enum.auto()


class CoarseTOrder(ModeEnum):
    STABILITY = enum.auto()  # diagonal
    INCOMING = enum.auto()
    OUTGOING = enum.auto()
    STAT_DIST = enum.auto()


@d.dedent
class GPCCA(TermStatesEstimator, LinDriversMixin, SchurMixin, EigenMixin):
    """Generalized Perron Cluster Cluster Analysis (GPCCA) :cite:`reuter:18,reuter:19`.

    .. seealso::
        - See :doc:`../../../notebooks/tutorials/estimators/600_initial_terminal` on how to compute the
          :attr:`initial <initial_states>` and :attr:`terminal <terminal_states>` states.
        - See :doc:`../../../notebooks/tutorials/estimators/700_fate_probabilities` on how to compute the
          :attr:`fate_probabilities` and :attr:`lineage_drivers`.

    This is our main and recommended estimator implemented in `pyGPCCA <https://pygpcca.readthedocs.io/en/latest/>`_ .
    Use it to compute macrostates, automatically and semi-automatically classify these as initial, intermediate and
    terminal states, compute fate probabilities towards macrostates, uncover driver genes, and much more. To compute and
    classify macrostates, we run the GPCCA algorithm under the hood, which returns a soft assignment of cells
    to macrostates, as well as a coarse-grained transition matrix among the set of macrostates
    :cite:`reuter:18,reuter:19`. This estimator allows you to inject prior knowledge where available
    to guide the identification of initial, intermediate and terminal states.

    Parameters
    ----------
    %(base_estimator.parameters)s
    """

    def __init__(
        self,
        object: Union[str, bool, np.ndarray, sp.spmatrix, AnnData, KernelExpression],
        **kwargs: Any,
    ):
        super().__init__(object=object, **kwargs)
        self._macrostates = StatesHolder()
        self._coarse_init_dist: Optional[pd.Series] = None
        self._coarse_stat_dist: Optional[pd.Series] = None
        self._coarse_tmat: Optional[pd.DataFrame] = None
        self._tsi: Optional[AnnData] = None

    @property
    @d.get_summary(base="gpcca_macro")
    def macrostates(self) -> Optional[pd.Series]:
        """Macrostates of the transition matrix."""
        return self._macrostates.assignment

    @property
    @d.get_summary(base="gpcca_macro_memberships")
    def macrostates_memberships(self) -> Optional[Lineage]:
        """Macrostate memberships.

        Soft assignment of microstates (cells) to macrostates.
        """
        return self._macrostates.memberships

    @property
    @d.get_summary(base="gpcca_init_states_memberships")
    def initial_states_memberships(self) -> Optional[Lineage]:
        """Initial states memberships.

        Soft assignment of cells to initial states.
        """
        return self._init_states.memberships

    @property
    @d.get_summary(base="gpcca_term_states_memberships")
    def terminal_states_memberships(self) -> Optional[Lineage]:
        """Terminal states memberships.

        Soft assignment of cells to terminal states.
        """
        return self._term_states.memberships

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
    ) -> "GPCCA":
        """Compute the macrostates.

        Parameters
        ----------
        n_states
            Number of macrostates to compute. If a :class:`~typing.Sequence`, use the *minChi*
            criterion :cite:`reuter:18`. If :obj:`None`, use the `eigengap <https://en.wikipedia.org/wiki/Eigengap>`__
            heuristic.
        %(n_cells)s
        cluster_key
            If a key to cluster labels is given, names and colors of the states will be associated with the clusters.
        kwargs
            Keyword arguments for :meth:`compute_schur`.

        Returns
        -------
        Returns self and updates the following fields:

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
            return self

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
        return self

    def predict(self, *args: Any, **kwargs: Any) -> "GPCCA":
        """Alias for :meth:`predict_terminal_states`.

        Parameters
        ----------
        args
            Positional arguments.
        kwargs
            Keyword arguments.

        Returns
        -------
        Same as :meth:`predict_terminal_states`.
        """
        return self.predict_terminal_states(*args, **kwargs)

    @d.dedent
    def predict_terminal_states(
        self,
        method: Literal["stability", "top_n", "eigengap", "eigengap_coarse"] = TermStatesMethod.STABILITY,
        n_cells: int = 30,
        alpha: Optional[float] = 1,
        stability_threshold: float = 0.96,
        n_states: Optional[int] = None,
        allow_overlap: bool = False,
    ) -> "GPCCA":
        """Automatically select terminal states from macrostates.

        Parameters
        ----------
        method
            How to select the terminal states. Valid option are:

            - ``'eigengap'`` - select the number of states based on the
              `eigengap <https://en.wikipedia.org/wiki/Eigengap>`__ of :attr:`transition_matrix`.
            - ``'eigengap_coarse'`` - select the number of states based on the *eigengap* of the diagonal
              of :attr:`coarse_T`.
            - ``'top_n'`` - select top ``n_states`` based on the probability of the diagonal of :attr:`coarse_T`.
            - ``'stability'`` - select states which have a stability >= ``stability_threshold``.
              The stability is given by the diagonal elements of :attr:`coarse_T`.
        %(n_cells)s
        alpha
            Weight given to the deviation of an eigenvalue from one.
            Only used when ``method = 'eigengap'`` or ``method = 'eigengap_coarse'``.
        stability_threshold
            Threshold used when ``method = 'stability'``.
        n_states
            Number of states used when ``method = 'top_n'``.
        %(allow_overlap)s

        Returns
        -------
        Returns self and updates the following fields:

        - :attr:`terminal_states` - %(tse_term_states.summary)s
        - :attr:`terminal_states_probabilities` - %(tse_term_states_probs.summary)s
        - :attr:`terminal_states_memberships` - %(gpcca_term_states_memberships.summary)s
        """
        if self.macrostates is None:
            raise RuntimeError("Compute macrostates first as `.compute_macrostates()`.")

        if len(self.macrostates.cat.categories) == 1:
            logg.warning("Found only one macrostate, making it the singular terminal state")
            return self.set_terminal_states(
                states=None,
                n_cells=n_cells,
                allow_overlap=allow_overlap,
                params=self._create_params(),
            )

        method = TermStatesMethod(method)
        eig = self.eigendecomposition
        coarse_T = self.coarse_T

        # fmt: off
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
            if n_states <= 0:
                raise ValueError(f"Expected `n_states` to be positive, found `{n_states}`.")
        elif method == TermStatesMethod.STABILITY:
            if stability_threshold is None:
                raise ValueError("Expected `stability_threshold != None` for `method='stability'`.")
            stability = pd.Series(np.diag(coarse_T), index=coarse_T.columns)
            names = stability[stability.values >= stability_threshold].index
            return self.set_terminal_states(
                names,
                n_cells=n_cells,
                allow_overlap=allow_overlap,
                params=self._create_params()
            )
        else:
            raise NotImplementedError(f"Method `{method}` is not yet implemented.")
        # fmt: on

        names = coarse_T.columns[np.argsort(np.diag(coarse_T))][-n_states:]
        return self.set_terminal_states(
            names,
            n_cells=n_cells,
            allow_overlap=allow_overlap,
            params=self._create_params(),
        )

    @d.dedent
    def predict_initial_states(self, n_states: int = 1, n_cells: int = 30, allow_overlap: bool = False) -> "GPCCA":
        """Compute initial states from macrostates using :attr:`coarse_stationary_distribution`.

        Parameters
        ----------
        n_states
            Number of initial states.
        %(n_cells)s
        %(allow_overlap)s

        Returns
        -------
        Returns self and updates the following fields:

        - :attr:`initial_states` - %(tse_init_states.summary)s
        - :attr:`initial_states_probabilities` - %(tse_init_states_probs.summary)s
        - :attr:`initial_states_memberships` - %(gpcca_init_states_memberships.summary)s
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
            return self.set_initial_states(states=None, n_cells=n_cells, allow_overlap=allow_overlap)

        stat_dist = self.coarse_stationary_distribution
        if stat_dist is None:
            raise RuntimeError("No coarse-grained stationary distribution found.")

        states = list(stat_dist.iloc[np.argsort(stat_dist)][:n_states].index)
        return self.set_initial_states(states, n_cells=n_cells, allow_overlap=allow_overlap)

    @d.dedent
    def set_terminal_states(
        self,
        states: Optional[Union[str, Sequence[str], Dict[str, Sequence[str]], pd.Series]] = None,
        n_cells: int = 30,
        allow_overlap: bool = False,
        cluster_key: Optional[str] = None,
        **kwargs: Any,
    ) -> "GPCCA":
        """Set the :attr:`terminal_states`.

        Parameters
        ----------
        states
            Which states to select. Valid options are:

            - :class:`str`, :class:`~typing.Sequence` - subset of :attr:`macrostates`. Multiple states can be
              combined using ``','``, such as ``['Alpha, Beta', 'Epsilon']``.
            - :class:`dict` - keys correspond to terminal states and values to cell IDs in
              :attr:`~anndata.AnnData.obs_names`.
            - :class:`~pandas.Series` - categorical series where each category corresponds to a macrostate.
              `NaN` values mark cells that should not be marked as :attr:`terminal_states`.
            - :obj:`None` - select all :attr:`macrostates`.
        %(n_cells)s
        %(allow_overlap)s
        cluster_key
            Key in :attr:`~anndata.AnnData.obs` to associate names and colors with :attr:`terminal_states`.
            Each state will be given the name and color corresponding to the cluster it mostly overlaps with.
            Only used when ``states`` is a :class:`dict` or :class:`~pandas.Series`.
        kwargs
            Additional keyword arguments.

        Returns
        -------
        Returns self and updates the following fields:

        - :attr:`terminal_states` - %(tse_term_states.summary)s
        - :attr:`terminal_states_probabilities` - %(tse_term_states_probs.summary)s
        - :attr:`terminal_states_memberships` - %(gpcca_term_states_memberships.summary)s
        """
        if isinstance(states, (dict, pd.Series)):
            return super().set_terminal_states(states, cluster_key=cluster_key, allow_overlap=allow_overlap, **kwargs)

        if n_cells <= 0:
            raise ValueError(f"Expected `n_cells` to be positive, found `{n_cells}`.")

        memberships = self.macrostates_memberships
        if memberships is None:
            raise RuntimeError("Compute macrostates first as `.compute_macrostates()`.")

        if states is None:
            states = memberships.names
        elif isinstance(states, str):
            states = [states]
        if not len(states):  # unset the states
            raise ValueError("No macrostates have been selected.")

        is_singleton = memberships.shape[1] == 1
        memberships = memberships[list(states)].copy()

        states = self._create_states(memberships, n_cells=n_cells, check_row_sums=False)
        if is_singleton:
            colors = self._macrostates.colors.copy()
            probs = memberships.X.squeeze() / memberships.X.max()
        else:
            colors = memberships[list(states.cat.categories)].colors
            probs = (memberships.X / memberships.X.max(axis=0)).max(axis=1)
        probs = pd.Series(probs, index=self.adata.obs_names)

        self._write_states(
            "terminal",
            states=states,
            colors=colors,
            probs=probs,
            memberships=memberships,
            params=kwargs.pop("params", {}),
            allow_overlap=allow_overlap,
            **kwargs,
        )
        return self

    @d.dedent
    def set_initial_states(
        self,
        states: Optional[Union[str, Sequence[str], Dict[str, Sequence[str]], pd.Series]] = None,
        n_cells: int = 30,
        allow_overlap: bool = False,
        cluster_key: Optional[str] = None,
        **kwargs: Any,
    ) -> "GPCCA":
        """Set the :attr:`initial_states`.

        Parameters
        ----------
        states
            Which states to select. Valid options are:

            - :class:`str`, :class:`~typing.Sequence` - subset of :attr:`macrostates`. Multiple states can be
              combined using ``','``, such as ``['Alpha, Beta', 'Epsilon']``.
            - :class:`dict` - keys correspond to initial states and values to cell IDs in
              :attr:`~anndata.AnnData.obs_names`.
            - :class:`~pandas.Series` - categorical series where each category corresponds to a macrostate.
              `NaN` values mark cells that should not be marked as :attr:`initial_states`.
            - :obj:`None` - select all :attr:`macrostates`.
        %(n_cells)s
        %(allow_overlap)s
        cluster_key
            Key in :attr:`~anndata.AnnData.obs` to associate names and colors with :attr:`initial_states`.
            Each state will be given the name and color corresponding to the cluster it mostly overlaps with.
            Only used when ``states`` is a :class:`dict` or :class:`~pandas.Series`.
        kwargs
            Additional keyword arguments.

        Returns
        -------
        Returns self and updates the following fields:

        - :attr:`initial_states` - %(tse_init_states.summary)s
        - :attr:`initial_states_probabilities` - %(tse_init_states_probs.summary)s
        - :attr:`initial_states_memberships` - %(gpcca_init_states_memberships.summary)s
        """
        if isinstance(states, (pd.Series, dict)):
            return super().set_initial_states(states, cluster_key=cluster_key, allow_overlap=allow_overlap, **kwargs)

        if n_cells <= 0:
            raise ValueError(f"Expected `n_cells` to be positive, found `{n_cells}`.")

        memberships = self.macrostates_memberships
        if memberships is None:
            raise RuntimeError("Compute macrostates first as `.compute_macrostates()`.")

        if states is None:
            states = memberships.names
        elif isinstance(states, str):
            states = [states]
        if not len(states):  # unset the states
            raise ValueError("No macrostates have been selected.")

        is_singleton = memberships.shape[1] == 1
        memberships = memberships[list(states)].copy()

        states = self._create_states(memberships, n_cells=n_cells, check_row_sums=False)
        if is_singleton:
            colors = self._macrostates.colors.copy()
            probs = memberships.X.squeeze() / memberships.X.max()
        else:
            colors = memberships[list(states.cat.categories)].colors
            probs = (memberships.X / memberships.X.max(axis=0)).max(axis=1)
        probs = pd.Series(probs, index=self.adata.obs_names)

        self._write_states(
            "initial",
            states=states,
            colors=colors,
            probs=probs,
            memberships=memberships,
            params=kwargs.pop("params", {}),
            allow_overlap=allow_overlap,
            **kwargs,
        )
        return self

    # TODO: Add definition/link to paper.
    def tsi(
        self,
        n_macrostates: int,
        terminal_states: Optional[List[str]] = None,
        cluster_key: Optional[str] = None,
        **kwargs: Any,
    ) -> float:
        """Compute terminal state identification (TSI) score.

        Parameters
        ----------
        n_macrostates
            Maximum number of macrostates to consider.
        terminal_states
            List of terminal states.
        cluster_key
            Key in :attr:`~anndata.AnnData.obs` defining cluster labels including terminal states.
        kwargs
            Keyword arguments passed to :meth:`compute_macrostates` function.

        Returns
        -------
        Returns TSI score.
        """
        tsi_precomputed = (self._tsi is not None) and (self._tsi[:, "number_of_macrostates"].X.max() >= n_macrostates)
        if terminal_states is not None:
            tsi_precomputed = tsi_precomputed and (set(self._tsi.uns["terminal_states"]) == set(terminal_states))
        if cluster_key is not None:
            tsi_precomputed = tsi_precomputed and (self._tsi.uns["cluster_key"] == cluster_key)

        if not tsi_precomputed:
            if terminal_states is None:
                raise RuntimeError("`terminal_states` needs to be specified to compute TSI.")
            if cluster_key is None:
                raise RuntimeError("`cluster_key` needs to be specified to compute TSI.")

            # create a new GPCCA object to avoid unsetting attributes
            # that depend on the macrostates, e.g. the terminal states
            g = self.copy(deep=True)
            macrostates = {}
            for n_states in range(n_macrostates, 0, -1):
                g = g.compute_macrostates(n_states=n_states, cluster_key=cluster_key, **kwargs)
                macrostates[n_states] = g.macrostates.cat.categories

            max_terminal_states = len(terminal_states)

            tsi_df = collections.defaultdict(list)
            for n_states, states in macrostates.items():
                n_terminal_states = (
                    states.str.replace(r"(_).*", "", regex=True).drop_duplicates().isin(terminal_states).sum()
                )
                tsi_df["number_of_macrostates"].append(n_states)
                tsi_df["identified_terminal_states"].append(n_terminal_states)

                tsi_df["optimal_identification"].append(min(n_states, max_terminal_states))

            tsi_df = AnnData(pd.DataFrame(tsi_df), uns={"terminal_states": terminal_states, "cluster_key": cluster_key})
            self._tsi = tsi_df

        tsi_df = self._tsi.to_df()
        row_mask = tsi_df["number_of_macrostates"] <= n_macrostates
        optimal_score = tsi_df.loc[row_mask, "optimal_identification"].sum()

        return tsi_df.loc[row_mask, "identified_terminal_states"].sum() / optimal_score

    @d.dedent
    def plot_tsi(
        self,
        n_macrostates: Optional[int] = None,
        x_offset: Tuple[float, float] = (0.2, 0.2),
        y_offset: Tuple[float, float] = (0.1, 0.1),
        figsize: Tuple[float, float] = (6, 4),
        dpi: Optional[int] = None,
        save: Optional[Union[str, Path]] = None,
        **kwargs: Any,
    ) -> Axes:
        """Plot terminal state identificiation (TSI).

        Requires computing TSI with :meth:`tsi` first.

        Parameters
        ----------
        n_macrostates
            Maximum number of macrostates to consider. Defaults to using all.
        x_offset
            Offset of x-axis.
        y_offset
            Offset of y-axis.
        %(plotting)s
        kwargs
            Keyword arguments for :func:`~seaborn.lineplot`.

        Returns
        -------
        Plot TSI of the kernel and an optimal identification strategy.
        """
        if self._tsi is None:
            raise RuntimeError("Compute TSI with `tsi` first as `.tsi()`.")

        tsi_df = self._tsi.to_df()
        if n_macrostates is not None:
            tsi_df = tsi_df.loc[tsi_df["number_of_macrostates"] <= n_macrostates, :]

        optimal_identification = tsi_df[["number_of_macrostates", "optimal_identification"]]
        optimal_identification = optimal_identification.rename(
            columns={"optimal_identification": "identified_terminal_states"}
        )
        optimal_identification["method"] = "Optimal identification"
        optimal_identification["line_style"] = "--"

        df = tsi_df[["number_of_macrostates", "identified_terminal_states"]]
        df["method"] = self.kernel.__class__.__name__
        df["line_style"] = "-"

        df = pd.concat([df, optimal_identification])

        fig, ax = plt.subplots(figsize=figsize, dpi=dpi, tight_layout=True)
        sns.lineplot(
            data=df,
            x="number_of_macrostates",
            y="identified_terminal_states",
            hue="method",
            style="line_style",
            drawstyle="steps-post",
            ax=ax,
            **kwargs,
        )

        ax.set_xticks(df["number_of_macrostates"].unique().astype(int))
        # Plot is generated from large to small values on the x-axis
        for label_id, label in enumerate(ax.xaxis.get_ticklabels()[::-1]):
            if ((label_id + 1) % 5 != 0) and label_id != 0:
                label.set_visible(False)
        ax.set_yticks(df["identified_terminal_states"].unique())

        x_min = df["number_of_macrostates"].min() - x_offset[0]
        x_max = df["number_of_macrostates"].max() + x_offset[1]
        y_min = df["identified_terminal_states"].min() - y_offset[0]
        y_max = df["identified_terminal_states"].max() + y_offset[1]
        ax.set(
            xlim=[x_min, x_max],
            ylim=[y_min, y_max],
            xlabel="Number of macrostates",
            ylabel="Identified terminal states",
        )

        ax.get_legend().remove()

        n_methods = len(df["method"].unique())
        handles, labels = ax.get_legend_handles_labels()
        handles[n_methods].set_linestyle("--")
        handles = handles[: (n_methods + 1)]
        labels = labels[: (n_methods + 1)]
        labels[0] = "Method"
        fig.legend(handles=handles, labels=labels, loc="lower center", ncol=(n_methods + 1), bbox_to_anchor=(0.5, -0.1))

        if save is not None:
            save_fig(fig=fig, path=save)

        return ax

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
        elif isinstance(n_states, int) and n_states == 1:
            self.compute_eigendecomposition()

        n = n_states if isinstance(n_states, int) else max(n_states)
        # call explicitly since `compute_macrostates` doesn't handle the case
        # when `minChi` is used for `n_states` and `self._gpcca` is uninitialized
        _ = self.compute_schur(n, **kwargs)
        return self.compute_macrostates(n_states=n_states, cluster_key=cluster_key, n_cells=n_cells)

    @d.dedent
    @inject_docs(o=CoarseTOrder)
    def plot_coarse_T(
        self,
        show_stationary_dist: bool = True,
        show_initial_dist: bool = False,
        order: Optional[Literal["stability", "incoming", "outgoing", "stat_dist"]] = "stability",
        cmap: Union[str, ListedColormap] = "viridis",
        xtick_rotation: float = 45,
        annotate: bool = True,
        show_cbar: bool = True,
        title: Optional[str] = None,
        figsize: Tuple[float, float] = (8, 8),
        dpi: int = 80,
        save: Optional[Union[str, pathlib.Path]] = None,
        text_kwargs: Mapping[str, Any] = types.MappingProxyType({}),
        **kwargs: Any,
    ) -> None:
        """Plot the coarse-grained transition matrix.

        Parameters
        ----------
        show_stationary_dist
            Whether to show the :attr:`coarse_stationary_distribution`, if present.
        show_initial_dist
            Whether to show the :attr:`coarse_initial_distribution`.
        order
            How to order the coarse-grained transition matrix. Valid options are:

            - ``{o.STABILITY!r}`` - order by the values on the diagonal.
            - ``{o.INCOMING!r}`` - order by the incoming mass, excluding the diagonal.
            - ``{o.OUTGOING!r}`` - order by the outgoing mass, excluding the diagonal.
            - ``{o.STAT_DIST!r}`` - order by coarse stationary distribution. If not present, use ``{o.STABILITY!r}``.
        cmap
            Colormap to use.
        xtick_rotation
            Rotation of ticks on the x-axis.
        annotate
            Whether to display the text on each cell.
        show_cbar
            Whether to show the colorbar.
        title
            Title of the figure.
        %(plotting)s
        text_kwargs
            Keyword arguments for :meth:`~matplotlib.axes.Axes.text`.
        kwargs
            Keyword arguments for :meth:`~matplotlib.axes.Axes.imshow`.

        Returns
        -------
        %(just_plots)s
        """

        def order_matrix(
            order: Optional[CoarseTOrder],
        ) -> Tuple[pd.DataFrame, Optional[pd.Series], Optional[pd.Series]]:
            coarse_T = self.coarse_T
            init_d = self.coarse_initial_distribution
            stat_d = self.coarse_stationary_distribution

            if order is None:
                return coarse_T, init_d, stat_d

            order = CoarseTOrder(order)
            if order == CoarseTOrder.STAT_DIST and stat_d is None:
                order = CoarseTOrder.STABILITY
                logg.warning(
                    f"Unable to order by `{CoarseTOrder.STAT_DIST}`, no coarse stationary distribution. "
                    f"Using `order={order}`"
                )

            if order == CoarseTOrder.INCOMING:
                values = (coarse_T.sum(0) - np.diag(coarse_T)).argsort(kind="stable")
                names = values.index[values]
            elif order == CoarseTOrder.OUTGOING:
                values = (coarse_T.sum(1) - np.diag(coarse_T)).argsort(kind="stable")
                names = values.index[values]
            elif order == CoarseTOrder.STABILITY:
                names = coarse_T.index[np.argsort(np.diag(coarse_T), kind="stable")]
            elif order == CoarseTOrder.STAT_DIST:
                names = stat_d.index[stat_d.argsort(kind="stable")]
            else:
                raise NotImplementedError(f"Order `{order}` is not yet implemented.")

            coarse_T = coarse_T.loc[names][names]
            if init_d is not None:
                init_d = init_d[names]
            if stat_d is not None:
                stat_d = stat_d[names]

            return coarse_T, init_d, stat_d

        def stylize_dist(ax: Axes, data: np.ndarray, xticks_labels: Sequence[str] = ()) -> None:
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
                ax.tick_params(which="both", top=False, right=False, bottom=False, left=False)

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

        def annotate_dist_ax(ax: Axes, data: np.ndarray, valfmt: str = "{x:.2f}") -> None:
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

        if self.coarse_T is None:
            raise RuntimeError(
                "Compute coarse-grained transition matrix first as `.compute_macrostates()` with `n_states > 1`."
            )

        coarse_T, coarse_init_d, coarse_stat_d = order_matrix(order)
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

        labels = list(coarse_T.columns)

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
            stylize_dist(init_ax, np.array(coarse_init_d).reshape(1, -1), xticks_labels=labels)

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
        save: Optional[Union[str, pathlib.Path]] = None,
        show: bool = True,
    ) -> Optional[Axes]:
        """Plot histogram of macrostates over categorical annotations.

        Parameters
        ----------
        %(adata)s
        key
            Key from :attr:`~anndata.AnnData.obs` containing categorical annotations.
        width
            Bar width in :math:`[0, 1]`.
        title
            Title of the figure. If :obj:`None`, create one automatically.
        labelrot
            Rotation of labels on x-axis.
        legend_loc
            Position of the legend. If :obj:`None`, don't show the legend.
        %(plotting)s
        show
            If `False`, return the :class:`~matplotlib.axes.Axes` object.

        Returns
        -------
        If ``show = True``, nothing, just plots, otherwise returns the axes object.
        Optionally saves it based on ``save``.
        """
        from cellrank.pl._utils import _position_legend

        macrostates = self.macrostates
        if macrostates is None:
            raise RuntimeError("Compute macrostates first as `.compute_macrostates()`.")
        if key not in self.adata.obs:
            raise KeyError(f"Data not found in `adata.obs[{key!r}]`.")
        if not isinstance(self.adata.obs[key].dtype, pd.CategoricalDtype):
            raise TypeError(
                f"Expected `adata.obs[{key!r}]` to be `categorical`, " f"found `{infer_dtype(self.adata.obs[key])}`."
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
            cats_colors = _create_categorical_colors(len(self.adata.obs[key].cat.categories))
        cat_color_mapper = dict(zip(self.adata.obs[key].cat.categories, cats_colors))
        x_indices = np.arange(len(macrostates.cat.categories))
        bottom = np.zeros_like(x_indices, dtype=float)

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
        if eig is not None and "stationary_dist" in eig and eig["params"]["which"] == "LR":
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
        time: Optional[datetime.datetime] = None,
        params: Dict[str, Any] = types.MappingProxyType({}),
    ) -> None:
        """Map fuzzy clustering to pre-computed annotations to get names and colors.

        Given the fuzzy clustering, we would like to select the most likely cells from each state and use these to
        give each state a name and a color by comparing with pre-computed, categorical cluster annotations.

        Parameters
        ----------
        memberships
            Fuzzy clustering.
        %(n_cells)s
        cluster_key
            Key from :attr:`~anndata.AnnData.obs` to get reference cluster annotations.
        check_row_sums
            Check whether rows in `memberships` sum to :math:`1`.
        time
            Start time of macrostates computation.
        params
            Parameters used in macrostates computation.

        Returns
        -------
        Nothing, just updates the fields as described in :meth:`compute_macrostates`.
        """
        if n_cells is None:
            # fmt: off
            logg.debug("Setting the macrostates using macrostate assignment")
            assignment = pd.Series(np.argmax(memberships, axis=1).astype(str), dtype="category")
            # sometimes, a category can be missing
            assignment = assignment.cat.reorder_categories([str(i) for i in range(memberships.shape[1])])
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
        # fmt: off
        self._write_states("initial", states=None, colors=None, probs=None, memberships=None, log=False)
        self._write_states("terminal", states=None, colors=None, probs=None, memberships=None, log=False)

        assignment, colors = self._set_categorical_labels(assignment, cluster_key=cluster_key)
        memberships = Lineage(memberships, names=list(assignment.cat.categories), colors=colors)
        # fmt: on

        groups = assignment.value_counts()
        groups = groups[groups != n_cells].to_dict()
        if len(groups):
            logg.warning(
                f"The following terminal states have different number " f"of cells than requested ({n_cells}): {groups}"
            )

        self._write_macrostates(
            macrostates=assignment,
            colors=colors,
            memberships=memberships,
            time=time,
            params=params,
        )

    @logger
    @shadow
    def _write_macrostates(
        self,
        macrostates: pd.Series,
        colors: np.ndarray,
        memberships: Lineage,
        params: Dict[str, Any] = types.MappingProxyType({}),
    ) -> str:
        # fmt: off
        key = Key.obs.macrostates(self.backward)
        self._set(obj=self.adata.obs, key=key, value=macrostates)
        ckey = Key.uns.colors(key)
        self._set(obj=self.adata.uns, key=ckey, value=colors)
        mkey = Key.obsm.memberships(key)
        self._set(obj=self.adata.obsm, key=mkey, value=memberships)
        self._macrostates = self._macrostates.set(assignment=macrostates, colors=colors, memberships=memberships)
        self.params[key] = dict(params)

        names = list(macrostates.cat.categories)
        if len(names) > 1:
            # not using stationary distribution
            g = self._gpcca
            tmat = pd.DataFrame(g.coarse_grained_transition_matrix, index=names, columns=names)
            init_dist = pd.Series(g.coarse_grained_input_distribution, index=names)
            if g.coarse_grained_stationary_probability is None:
                stat_dist = None
            else:
                stat_dist = pd.Series(g.coarse_grained_stationary_probability, index=names)
            dists = pd.DataFrame({"coarse_init_dist": init_dist}, index=names)
            if stat_dist is not None:
                dists["coarse_stat_dist"] = stat_dist

            key = Key.obsm.schur_vectors(self.backward)
            self._set("_schur_vectors", obj=self.adata.obsm, key=key, value=g._p_X, shadow_only=True)
            key = Key.uns.schur_matrix(self.backward)
            self._set("_schur_matrix", obj=self.adata.uns, key=key, value=g._p_R, shadow_only=True)
            self._set("_coarse_tmat", value=tmat, shadow_only=True)
            self._set("_coarse_init_dist", value=init_dist, shadow_only=True)
            self._set("_coarse_stat_dist", value=stat_dist, shadow_only=True)
            self._set(
                obj=self.adata.uns, key=Key.uns.coarse(self.backward),
                value=AnnData(tmat, obs=dists),
            )
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
    def _write_states(
        self,
        which: Literal["initial", "terminal"],
        states: Optional[pd.Series],
        colors: Optional[np.ndarray],
        probs: Optional[pd.Series] = None,
        memberships: Optional[Lineage] = None,
        params: Dict[str, Any] = types.MappingProxyType({}),
        allow_overlap: bool = False,
    ) -> str:
        msg = super()._write_states(
            which,
            states=states,
            colors=colors,
            probs=probs,
            params=params,
            allow_overlap=allow_overlap,
            log=False,
        )
        # fmt: off
        msg = "\n".join(msg.split("\n")[:-1])
        msg += f"\n       `.{which}_states_memberships\n    Finish`"

        self._write_fate_probabilities(None, log=False)
        # TODO(michalk8): CFLARE doesn't remove the downstream properties
        self._write_absorption_times(None, log=False)

        backward = which == "initial"
        key = Key.obsm.memberships(Key.obs.term_states(self.backward, bwd=backward))
        self._set(obj=self.adata.obsm, key=key, value=memberships)
        if backward:
            self._init_states = self._init_states.set(memberships=memberships)
        else:
            self._term_states = self._term_states.set(memberships=memberships)
        # fmt: on

        return msg

    def _read_from_adata(self, adata: AnnData, **kwargs: Any) -> bool:
        # design choice: no need to eigen/Schur-decomposition
        _ = self._read_eigendecomposition(adata, allow_missing=True)
        ok = self._read_schur_decomposition(adata, allow_missing=True)
        if not ok:
            return False

        # fmt: off
        with SafeGetter(self, allowed=KeyError) as sg:
            key = Key.obs.macrostates(self.backward)
            # TODO(michalk8): in the future, be more stringent and ensure the categories match the macro memberships
            assignment = self._get(obj=adata.obs, key=key, shadow_attr="obs", dtype=pd.Series)
            ckey = Key.uns.colors(key)
            colors = self._get(obj=adata.uns, key=ckey, shadow_attr="uns", dtype=(list, tuple, np.ndarray))
            mkey = Key.obsm.memberships(key)
            memberships = self._get(obj=adata.obsm, key=mkey, shadow_attr="obsm", dtype=(Lineage, np.ndarray))
            memberships = self._ensure_lineage_object(memberships, backward=self.backward, kind="macrostates")
            self._macrostates = StatesHolder(assignment=assignment, colors=colors, memberships=memberships)
            self.params[key] = self._read_params(key)

            tmat = adata.uns[Key.uns.coarse(self.backward)].copy()
            if not isinstance(tmat, AnnData):
                raise TypeError(
                    f"Expected coarse-grained transition matrix to be stored "
                    f"as `AnnData`, found `{type(tmat).__name__}`."
                )

            self._coarse_tmat = pd.DataFrame(tmat.X, index=tmat.obs_names, columns=tmat.obs_names)
            self._coarse_init_dist = tmat.obs["coarse_init_dist"]
            self._coarse_stat_dist = tmat.obs.get("coarse_stat_dist", None)

            self._set(obj=self._shadow_adata.uns, key=Key.uns.coarse(self.backward), value=tmat)

        if not sg.ok:
            return False

        if not super()._read_from_adata(adata, **kwargs):
            return False

        # status is based on `backward=False` by design
        for backward in [True, False]:
            with SafeGetter(self, allowed=KeyError) as sg:
                key = Key.obsm.memberships(Key.obs.term_states(self.backward, bwd=backward))
                memberships = self._get(obj=adata.obsm, key=key, shadow_attr="obsm", dtype=(np.ndarray, Lineage))
                memberships = self._ensure_lineage_object(memberships, backward=backward, kind="term_states")
                if backward:
                    self._init_states = self._init_states.set(memberships=memberships)
                else:
                    self._term_states = self._term_states.set(memberships=memberships)
                self._set(obj=self._shadow_adata.obsm, key=key, value=memberships)
        # fmt: on

        fate_prob_ok = self._read_fate_probabilities(adata)
        abs_time_ok = self._read_absorption_times(adata)
        return sg.ok and fate_prob_ok and abs_time_ok
