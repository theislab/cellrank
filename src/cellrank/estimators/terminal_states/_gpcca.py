import collections
import enum
import logging
import pathlib
import time as _time
import types
import warnings
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any, Literal

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
import seaborn as sns
from anndata import AnnData
from matplotlib.axes import Axes
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.ticker import StrMethodFormatter
from pandas.api.types import infer_dtype

from cellrank._utils._colors import _create_categorical_colors, _get_black_or_white
from cellrank._utils._docs import d, inject_docs
from cellrank._utils._enum import ModeEnum
from cellrank._utils._key import Key
from cellrank._utils._lineage import Lineage
from cellrank._utils._utils import (
    _aggregate_proportions,
    _check_proportions,
    _eigengap,
    _fuzzy_to_discrete,
    _obsm_proportion_weights,
    _obsm_proportions,
    _resolve_composition_key,
    _series_from_one_hot_matrix,
    save_fig,
)
from cellrank.estimators.mixins import EigenMixin, LinDriversMixin, SchurMixin
from cellrank.estimators.mixins._utils import SafeGetter, StatesHolder, log_writer, shadow
from cellrank.estimators.terminal_states._term_states_estimator import (
    TermStatesEstimator,
)
from cellrank.kernels._base_kernel import KernelExpression

logger = logging.getLogger(__name__)
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


class StatesAggregation(ModeEnum):
    TOP_N = enum.auto()
    UNION = enum.auto()


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
        object: str | bool | np.ndarray | sp.spmatrix | AnnData | KernelExpression,
        **kwargs: Any,
    ):
        super().__init__(object=object, **kwargs)
        self._macrostates = StatesHolder()
        self._coarse_init_dist: pd.Series | None = None
        self._coarse_stat_dist: pd.Series | None = None
        self._coarse_tmat: pd.DataFrame | None = None
        self._tsi: AnnData | None = None

    @property
    @d.get_summary(base="gpcca_macro")
    def macrostates(self) -> pd.Series | None:
        """Macrostates of the transition matrix."""
        return self._macrostates.assignment

    @property
    @d.get_summary(base="gpcca_macro_memberships")
    def macrostates_memberships(self) -> Lineage | None:
        """Macrostate memberships.

        Soft assignment of microstates (cells) to macrostates.
        """
        return self._macrostates.memberships

    @property
    @d.get_summary(base="gpcca_init_states_memberships")
    def initial_states_memberships(self) -> Lineage | None:
        """Initial states memberships.

        Soft assignment of cells to initial states.
        """
        return self._init_states.memberships

    @property
    @d.get_summary(base="gpcca_term_states_memberships")
    def terminal_states_memberships(self) -> Lineage | None:
        """Terminal states memberships.

        Soft assignment of cells to terminal states.
        """
        return self._term_states.memberships

    @property
    @d.get_summary(base="gpcca_coarse_init")
    def coarse_initial_distribution(self) -> pd.Series | None:
        """Coarse-grained initial distribution."""
        return self._coarse_init_dist

    @property
    @d.get_summary(base="gpcca_coarse_stat")
    def coarse_stationary_distribution(self) -> pd.Series | None:
        """Coarse-grained stationary distribution."""
        return self._coarse_stat_dist

    @property
    @d.get_summary(base="gpcca_coarse_tmat")
    def coarse_T(self) -> pd.DataFrame | None:
        """Coarse-grained transition matrix."""
        return self._coarse_tmat

    @d.get_sections(base="gpcca_compute_macro", sections=["Parameters", "Returns"])
    @d.dedent
    def compute_macrostates(
        self,
        n_states: int | Sequence[int] | None = None,
        n_cells: int | None = 30,
        cluster_key: str | Mapping[str, str] | None = None,
        weight_key: str | None = None,
        optimizer: str = "CG",
        n_starts: int = 10,
        perturbation_scale: float = 0.1,
        seed: int | None = None,
        **kwargs: Any,
    ) -> "GPCCA":
        r"""Compute the macrostates.

        Parameters
        ----------
        n_states
            Number of macrostates to compute. If a :class:`~typing.Sequence`, use the *minChi*
            criterion :cite:`reuter:18`. If :obj:`None`, use the `eigengap <https://en.wikipedia.org/wiki/Eigengap>`__
            heuristic.
        %(n_cells)s
        cluster_key
            If given, names and colors of the states will be associated with these reference annotations.
            Either:

            - a :class:`str` (or ``{"obs": <column>}``), a categorical column in
              :attr:`~anndata.AnnData.obs`. Each macrostate is named after the dominant category
              among its most-likely observations (unchanged behavior).
            - ``{"obsm": <key>}``, where :attr:`adata.obsm[key] <anndata.AnnData.obsm>` is a
              :class:`~pandas.DataFrame` whose columns are the categories and whose rows are
              per-observation proportions summing to :math:`1` (e.g. cell-type or condition
              fractions of aggregated samples). For each macrostate, the proportion rows of its
              most-likely observations are summed (see ``weight_key`` for the weighting) and the
              macrostate is named after the category with the largest resulting total.
        weight_key
            Only used when ``cluster_key`` points to :attr:`~anndata.AnnData.obsm`. Controls how
            each observation's proportion row is weighted when summing proportions to pick a
            macrostate's dominant category:

            - :obj:`None` (default) - **every observation contributes equally**, with weight
              :math:`1`, irrespective of any per-observation quantity.
            - a key in :attr:`~anndata.AnnData.obs` - each observation's proportion row is scaled
              by ``adata.obs[weight_key]`` before summing, so observations with larger weights
              count proportionally more. The weights can be any per-observation quantity; a common
              choice is the number of cells per aggregated sample, which makes naming reflect
              cell-level rather than sample-level dominance.
        optimizer
            Optimization method for the GPCCA rotation matrix. Valid options are
            ``'CG'``, ``'L-BFGS-B'``, ``'BFGS'``, and ``'Nelder-Mead'``.
            ``'CG'`` (conjugate gradient) uses an analytical Jacobian and scales
            as :math:`O(n \cdot m^2)` per iteration, where *n* is the number of
            cells and *m* is the number of macrostates.
            ``'Nelder-Mead'`` is derivative-free and was the original default.
            It optimizes over :math:`O(m^2)` rotation parameters and becomes
            impractical for ``n_states > 15``.
        n_starts
            Number of optimization runs. The first run uses the deterministic
            ISA initialization; subsequent runs perturb the initial rotation
            matrix on the SO(k) manifold and keep the result with the best
            crispness. Set to ``1`` to disable perturbation (fully
            deterministic).
        perturbation_scale
            Angular scale for the rotation perturbation (only used when
            ``n_starts > 1``). Recommended range is 0.05-0.2; larger values
            risk producing degenerate solutions.
        seed
            Random seed for reproducibility of the perturbations
            (only used when ``n_starts > 1``).
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
                weight_key=weight_key,
            )
            return self

        if self._gpcca is None or kwargs:
            self.compute_schur(n_states, **kwargs)
        n_states = self._validate_n_states(n_states)

        if self._gpcca._p_X.shape[1] < n_states:
            # precomputed X
            logger.warning(
                "Requested more macrostates `%s` than available Schur vectors `%s`. Recomputing the decomposition",
                n_states,
                self._gpcca._p_X.shape[1],
            )

        if n_states > 15 and optimizer == "Nelder-Mead":
            logger.warning(
                "Using Nelder-Mead with `n_states=%s` may be very slow. "
                "Consider using optimizer='CG' (the default) or another gradient-based optimizer",
                n_states,
            )

        _start = _time.perf_counter()
        logger.info("Computing %d macrostates", n_states)
        try:
            self._gpcca = self._gpcca.optimize(
                m=n_states,
                method=optimizer,
                n_starts=n_starts,
                perturbation_scale=perturbation_scale,
                seed=seed,
            )
        except ValueError as e:
            if "will split complex conjugate eigenvalues" not in str(e):
                raise
            # this is the following case - we have 4 Schur vectors, user requests 5 states, but it splits the conj. ev.
            # in the try block, Schur decomposition with 5 vectors is computed, but it fails (no way of knowing)
            # so in this case, we increase it by 1
            logger.warning(
                "Unable to compute macrostates with `n_states=%s` because it will "
                "split complex conjugate eigenvalues. Using `n_states=%s`",
                n_states,
                n_states + 1,
            )
            self._gpcca = self._gpcca.optimize(
                m=n_states + 1,
                method=optimizer,
                n_starts=n_starts,
                perturbation_scale=perturbation_scale,
                seed=seed,
            )

        self._set_macrostates(
            memberships=self._gpcca.memberships,
            n_cells=n_cells,
            cluster_key=cluster_key,
            weight_key=weight_key,
            params=self._create_params(),
            time=_start,
        )
        return self

    def predict(self, *args: Any, **kwargs: Any) -> "GPCCA":
        """Alias for :meth:`predict_terminal_states`.

        .. deprecated:: 2.1
            Will be removed in CellRank 3.0. Use :meth:`predict_terminal_states` directly.

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
        warnings.warn(
            "`GPCCA.predict()` is deprecated and will be removed in CellRank 3.0. "
            "Use `GPCCA.predict_terminal_states()` directly.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.predict_terminal_states(*args, **kwargs)

    @d.dedent
    def predict_terminal_states(
        self,
        method: Literal["stability", "top_n", "eigengap", "eigengap_coarse"] = TermStatesMethod.STABILITY,
        n_cells: int = 30,
        alpha: float | None = 1,
        stability_threshold: float = 0.96,
        n_states: int | None = None,
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
            logger.warning("Found only one macrostate, making it the singular terminal state")
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
            if not len(names):
                raise ValueError(
                    f"No macrostate has a stability `>= {stability_threshold}` "
                    f"(maximum stability is `{stability.max():.4f}`), so no terminal states could be selected. "
                    f"Decrease `stability_threshold`, use a different `method` (e.g. `'top_n'`), "
                    f"or recompute the macrostates with `.compute_macrostates()`."
                )
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
            raise RuntimeError(
                "Unable to predict initial states because no coarse-grained stationary distribution is "
                "available. It could not be computed because the stationary distribution of the transition "
                "matrix could not be calculated, which typically happens for very large or (nearly) reducible "
                "transition matrices where the underlying solver fails to converge. Initial state prediction "
                "relies on this distribution, so it cannot be done automatically here. Set the initial states "
                "manually using `.set_initial_states()` instead."
            )

        states = list(stat_dist.iloc[np.argsort(stat_dist)][:n_states].index)
        return self.set_initial_states(states, n_cells=n_cells, allow_overlap=allow_overlap)

    @d.dedent
    def set_terminal_states(
        self,
        states: str | Sequence[str] | dict[str, Sequence[str]] | pd.Series | None = None,
        n_cells: int = 30,
        allow_overlap: bool = False,
        cluster_key: str | Mapping[str, str] | None = None,
        weight_key: str | None = None,
        agg: Literal["top_n", "union"] = StatesAggregation.TOP_N,
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
            Reference annotations to associate names and colors with :attr:`terminal_states`; see
            :meth:`compute_macrostates`. Only used when ``states`` is a :class:`dict` or :class:`~pandas.Series`.
        weight_key
            Per-observation weights, only used when ``cluster_key`` points to
            :attr:`~anndata.AnnData.obsm`; see :meth:`compute_macrostates`.
        %(agg)s
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
            return super().set_terminal_states(
                states, cluster_key=cluster_key, weight_key=weight_key, allow_overlap=allow_overlap, **kwargs
            )

        if n_cells <= 0:
            raise ValueError(f"Expected `n_cells` to be positive, found `{n_cells}`.")
        agg = StatesAggregation(agg)

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
        macrostates_memberships, selected = memberships, list(states)
        memberships = memberships[selected].copy()

        if agg == StatesAggregation.UNION and not is_singleton:
            states = self._create_aggregated_states(
                macrostates_memberships, states=selected, group_labels=list(memberships.names), n_cells=n_cells
            )
        else:
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
        states: str | Sequence[str] | dict[str, Sequence[str]] | pd.Series | None = None,
        n_cells: int = 30,
        allow_overlap: bool = False,
        cluster_key: str | Mapping[str, str] | None = None,
        weight_key: str | None = None,
        agg: Literal["top_n", "union"] = StatesAggregation.TOP_N,
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
            Reference annotations to associate names and colors with :attr:`initial_states`; see
            :meth:`compute_macrostates`. Only used when ``states`` is a :class:`dict` or :class:`~pandas.Series`.
        weight_key
            Per-observation weights, only used when ``cluster_key`` points to
            :attr:`~anndata.AnnData.obsm`; see :meth:`compute_macrostates`.
        %(agg)s
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
            return super().set_initial_states(
                states, cluster_key=cluster_key, weight_key=weight_key, allow_overlap=allow_overlap, **kwargs
            )

        if n_cells <= 0:
            raise ValueError(f"Expected `n_cells` to be positive, found `{n_cells}`.")
        agg = StatesAggregation(agg)

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
        macrostates_memberships, selected = memberships, list(states)
        memberships = memberships[selected].copy()

        if agg == StatesAggregation.UNION and not is_singleton:
            states = self._create_aggregated_states(
                macrostates_memberships, states=selected, group_labels=list(memberships.names), n_cells=n_cells
            )
        else:
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
        terminal_states: list[str] | None = None,
        cluster_key: str | None = None,
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
        if not (cluster_key is None or isinstance(cluster_key, str)):
            # `cluster_key` is round-tripped through `.uns` for cache validity, so it must be serializable;
            # the `{"obsm": ...}` proportion form supported by `compute_macrostates` is not allowed here.
            raise TypeError(f"Expected `cluster_key` to be a `str` or `None`, found `{type(cluster_key).__name__!r}`.")

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
            comp_failed = {}
            for n_states in range(n_macrostates, 0, -1):
                try:
                    g = g.compute_macrostates(n_states=n_states, cluster_key=cluster_key, **kwargs)
                    macrostates[n_states] = g.macrostates.cat.categories
                    comp_failed[n_states] = False
                except ValueError:
                    comp_failed[n_states] = True
                    macrostates[n_states] = pd.Series(np.nan, index=self.adata.obs_names)

            max_terminal_states = len(terminal_states)

            tsi_df = collections.defaultdict(list)
            for n_states, states in macrostates.items():
                if comp_failed[n_states] and (n_states == 1):
                    n_terminal_states = 0
                elif comp_failed[n_states] and (n_states > 1):
                    n_terminal_states = np.nan
                    tsi_df["identified_terminal_states"][-1]
                else:
                    n_terminal_states = (
                        states.str.replace(r"(_).*", "", regex=True).drop_duplicates().isin(terminal_states).sum()
                    )
                tsi_df["number_of_macrostates"].append(n_states)
                tsi_df["identified_terminal_states"].append(n_terminal_states)

                tsi_df["optimal_identification"].append(min(n_states, max_terminal_states))

            tsi_df = pd.DataFrame(tsi_df)
            tsi_df.sort_values(by="number_of_macrostates", inplace=True)
            tsi_df["identified_terminal_states"] = tsi_df["identified_terminal_states"].ffill()
            tsi_df.sort_values(by="number_of_macrostates", ascending=False, inplace=True)
            tsi_df.index = tsi_df.index.astype(str)
            tsi_df = AnnData(tsi_df, uns={"terminal_states": terminal_states, "cluster_key": cluster_key})
            self._tsi = tsi_df

        tsi_df = self._tsi.to_df()
        row_mask = tsi_df["number_of_macrostates"] <= n_macrostates
        optimal_score = tsi_df.loc[row_mask, "optimal_identification"].sum()

        return tsi_df.loc[row_mask, "identified_terminal_states"].sum() / optimal_score

    @d.dedent
    def plot_tsi(
        self,
        n_macrostates: int | None = None,
        x_offset: tuple[float, float] = (0.2, 0.2),
        y_offset: tuple[float, float] = (0.1, 0.1),
        figsize: tuple[float, float] = (6, 4),
        dpi: int | None = None,
        save: str | Path | None = None,
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
        n_states: int | Sequence[int] | None = None,
        n_cells: int | None = 30,
        cluster_key: str | Mapping[str, str] | None = None,
        weight_key: str | None = None,
        optimizer: str = "CG",
        n_starts: int = 10,
        perturbation_scale: float = 0.1,
        seed: int | None = None,
        **kwargs: Any,
    ) -> "GPCCA":
        """
        Prepare self for terminal states prediction.

        .. deprecated:: 2.1
            Will be removed in CellRank 3.0. Use :meth:`compute_schur` and
            :meth:`compute_macrostates` directly.

        Parameters
        ----------
        %(gpcca_compute_macro.parameters)s

        Returns
        -------
        %(gpcca_compute_macro.returns)s
        """
        warnings.warn(
            "`GPCCA.fit()` is deprecated and will be removed in CellRank 3.0. "
            "Use `GPCCA.compute_schur()` and `GPCCA.compute_macrostates()` directly.",
            DeprecationWarning,
            stacklevel=2,
        )
        if n_states is None:
            self.compute_eigendecomposition()
            n_states = self.eigendecomposition["eigengap"] + 1
        elif isinstance(n_states, int) and n_states == 1:
            self.compute_eigendecomposition()

        n = n_states if isinstance(n_states, int) else max(n_states)
        # call explicitly since `compute_macrostates` doesn't handle the case
        # when `minChi` is used for `n_states` and `self._gpcca` is uninitialized
        _ = self.compute_schur(n, **kwargs)
        return self.compute_macrostates(
            n_states=n_states,
            cluster_key=cluster_key,
            weight_key=weight_key,
            n_cells=n_cells,
            optimizer=optimizer,
            n_starts=n_starts,
            perturbation_scale=perturbation_scale,
            seed=seed,
        )

    @d.dedent
    @inject_docs(o=CoarseTOrder)
    def plot_coarse_T(
        self,
        show_stationary_dist: bool = True,
        show_initial_dist: bool = False,
        order: Literal["stability", "incoming", "outgoing", "stat_dist"] | None = "stability",
        cmap: str | ListedColormap = "viridis",
        xtick_rotation: float = 45,
        annotate: bool = True,
        show_cbar: bool = True,
        title: str | None = None,
        figsize: tuple[float, float] = (8, 8),
        dpi: int = 80,
        save: str | pathlib.Path | None = None,
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
            order: CoarseTOrder | None,
        ) -> tuple[pd.DataFrame, pd.Series | None, pd.Series | None]:
            coarse_T = self.coarse_T
            init_d = self.coarse_initial_distribution
            stat_d = self.coarse_stationary_distribution

            if order is None:
                return coarse_T, init_d, stat_d

            order = CoarseTOrder(order)
            if order == CoarseTOrder.STAT_DIST and stat_d is None:
                order = CoarseTOrder.STABILITY
                logger.warning(
                    "Unable to order by `%s`, no coarse stationary distribution. Using `order=%s`",
                    CoarseTOrder.STAT_DIST,
                    order,
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
            logger.warning("Coarse stationary distribution is `None`, ignoring")
            show_stationary_dist = False
        if show_initial_dist and coarse_init_d is None:
            logger.warning("Coarse initial distribution is `None`, ignoring")
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
        key: str | Mapping[str, str],
        weight_key: str | None = None,
        width: float = 0.8,
        title: str | None = None,
        labelrot: float = 45,
        legend_loc: str | None = "upper right out",
        figsize: tuple[float, float] | None = None,
        dpi: int | None = None,
        save: str | pathlib.Path | None = None,
        show: bool = True,
    ) -> Axes | None:
        """Plot histogram of macrostates over categorical annotations.

        The annotation can either be a hard categorical assignment, with one label per
        observation, or a soft assignment, with a distribution over categories per
        observation (e.g. cell-type fractions of aggregated samples).

        Parameters
        ----------
        %(adata)s
        key
            Source of the categorical annotation. Either:

            - a :class:`str`, interpreted as a categorical column in
              :attr:`~anndata.AnnData.obs`. Each macrostate's bar stacks the *counts* of
              the categories among its most-likely observations.
            - a :class:`dict` of the form ``{"obs": <column>}`` (equivalent to passing a
              :class:`str`) or ``{"obsm": <key>}``. In the latter case,
              :attr:`adata.obsm[key] <anndata.AnnData.obsm>` must be a
              :class:`~pandas.DataFrame` whose columns are the categories and whose rows
              are per-observation proportions summing to :math:`1` (e.g. cell-type
              fractions of aggregated samples). The proportion rows are summed per
              macrostate, so the bar semantics are identical to the categorical case --
              each observation contributes :math:`1`, split across categories -- only the
              observations are now samples rather than cells.
        weight_key
            Only used when ``key`` points to :attr:`~anndata.AnnData.obsm`. Key from
            :attr:`~anndata.AnnData.obs` with per-observation weights, e.g. the number of
            cells per aggregated sample, so the bars reflect cell-level rather than
            sample-level frequencies. If :obj:`None`, each observation contributes equally.
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
        from cellrank.pl._utils import _position_legend  # circular import

        macrostates = self.macrostates
        if macrostates is None:
            raise RuntimeError("Compute macrostates first as `.compute_macrostates()`.")

        source, name = _resolve_composition_key(key)
        if weight_key is not None and source != "obsm":
            raise ValueError("`weight_key` is only supported when `key` points to `adata.obsm`.")

        mask = ~macrostates.isnull()
        if source == "obs":
            composition, categories = self._obs_composition(name, macrostates, mask)
        else:
            composition, categories = self._obsm_composition(name, weight_key, macrostates, mask)

        try:
            cats_colors = self.adata.uns[f"{name}_colors"]
        except KeyError:
            cats_colors = _create_categorical_colors(len(categories))
        if len(cats_colors) < len(categories):
            cats_colors = _create_categorical_colors(len(categories))
        cat_color_mapper = dict(zip(categories, cats_colors))
        n_states = len(macrostates.cat.categories)
        x_indices = np.arange(n_states)
        bottom = np.zeros_like(x_indices, dtype=float)

        if figsize is None:
            # scale width with the number of macrostates and reserve room for an
            # outside legend, otherwise the bars get squeezed into a narrow strip
            fig_width = max(6.0, 1.1 * n_states)
            if legend_loc not in (None, "none") and legend_loc.endswith("out"):
                fig_width += 0.18 * max(len(str(cat)) for cat in categories) + 1.5
            figsize = (fig_width, 5.0)

        width = min(1, max(0, width))
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi, tight_layout=True)
        for cat, color in cat_color_mapper.items():
            frequencies = composition.loc[cat].to_numpy()
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
                bottom += frequencies

        ax.set_xticks(x_indices)
        ax.set_xticklabels(
            macrostates.cat.categories,
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
            title = f"distribution over {name}"
        ax.set_title(title)
        if legend_loc not in (None, "none"):
            _position_legend(ax, legend_loc=legend_loc)

        if save is not None:
            save_fig(fig, save)

        if not show:
            return ax

    def _obs_composition(self, key: str, macrostates: pd.Series, mask: pd.Series) -> tuple[pd.DataFrame, pd.Index]:
        """Count each :attr:`~anndata.AnnData.obs` category per macrostate."""
        if key not in self.adata.obs:
            raise KeyError(f"Data not found in `adata.obs[{key!r}]`.")
        if not isinstance(self.adata.obs[key].dtype, pd.CategoricalDtype):
            raise TypeError(
                f"Expected `adata.obs[{key!r}]` to be `categorical`, found `{infer_dtype(self.adata.obs[key])}`."
            )

        categories = self.adata.obs[key].cat.categories
        composition = (
            pd.DataFrame({"macrostates": macrostates, key: self.adata.obs[key]})[mask]
            .groupby([key, "macrostates"], observed=False)
            .size()
            .unstack("macrostates")
            .reindex(index=categories, columns=macrostates.cat.categories, fill_value=0)
        )
        return composition, categories

    def _obsm_composition(
        self, key: str, weight_key: str | None, macrostates: pd.Series, mask: pd.Series
    ) -> tuple[pd.DataFrame, pd.Index]:
        """Sum the (weighted) proportions per macrostate from an :attr:`~anndata.AnnData.obsm` frame.

        Each observation contributes its proportion row, so an unweighted bar height equals the number
        of observations in the macrostate -- the direct analog of the cell counts in :meth:`_obs_composition`.
        """
        fractions = _obsm_proportions(self.adata, key)

        assigned = macrostates[mask]
        fractions = fractions.loc[assigned.index]
        _check_proportions(fractions, key)
        weights = _obsm_proportion_weights(self.adata, weight_key, assigned.index)

        categories = fractions.columns
        # `_aggregate_proportions` returns macrostates x categories; transpose for the stacked bars
        composition = _aggregate_proportions(fractions, assigned, weights).T
        return composition, categories

    def _n_states(self, n_states: int | Sequence[int] | None) -> int:
        if n_states is None:
            if self.eigendecomposition is None:
                raise RuntimeError(
                    "Compute eigendecomposition first as `.compute_eigendecomposition()` or supply `n_states != None`."
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
            logger.warning("Minimum value must be larger than `1`, found `%s`. Setting `min=%s`", minn, minn)
        if minn == 2:
            minn = 3
            logger.warning(
                "In most cases, 2 clusters will always be optimal. "
                "If you really expect 2 clusters, use `n_states=2`. Setting `min=%s`",
                minn,
            )
        # fmt: on
        maxx = max(minn + 1, maxx)

        logger.info("Calculating minChi criterion in interval `[%s, %s]`", minn, maxx)

        return int(np.arange(minn, maxx + 1)[np.argmax(self._gpcca.minChi(minn, maxx))])

    def _create_states(
        self,
        probs: np.ndarray | Lineage,
        n_cells: int,
        check_row_sums: bool = False,
        return_not_enough_cells: bool = False,
    ) -> pd.Series | tuple[pd.Series, np.ndarray]:
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

    def _create_aggregated_states(
        self,
        memberships: Lineage,
        states: Sequence[str],
        group_labels: Sequence[str],
        n_cells: int,
    ) -> pd.Series:
        """Select the ``n_cells`` most likely cells of *each* macrostate and union them per (aggregated) state.

        In contrast to selecting the most likely cells of the *combined* membership, this guarantees that every
        constituent macrostate is represented, even if one of them dominates the combined membership.

        Parameters
        ----------
        memberships
            Memberships of the individual macrostates.
        states
            Selected states, where each entry may combine several macrostates, e.g. ``'Alpha, Beta'``.
        group_labels
            Canonical label of each entry in ``states`` (i.e. the names of the aggregated ``memberships``).
        n_cells
            Number of most likely cells to select from each individual macrostate.

        Returns
        -------
        Categorical series assigning the selected cells to the aggregated states.
        """
        constituents: list[str] = []
        mapping: dict[str, str] = {}
        for state, label in zip(states, group_labels):
            for name in sorted({n.strip(" ") for n in str(state).strip(" ,").split(",")}):
                constituents.append(name)
                mapping[name] = label

        per_macrostate = self._create_states(memberships[constituents], n_cells=n_cells, check_row_sums=False)
        aggregated = per_macrostate.map(mapping)
        # preserve the canonical order and drop states that ended up with no cells
        present = [label for label in group_labels if label in set(aggregated.dropna())]
        return pd.Series(
            pd.Categorical(aggregated.to_numpy(), categories=present),
            index=per_macrostate.index,
        )

    def _validate_n_states(self, n_states: int) -> int:
        if self._invalid_n_states is not None and n_states in self._invalid_n_states:
            logger.warning(
                "Unable to compute macrostates with `n_states=%s` because it will "
                "split complex conjugate eigenvalues. Using `n_states=%s`",
                n_states,
                n_states + 1,
            )
            n_states += 1  # cannot force recomputation of the Schur decomposition
            assert n_states not in self._invalid_n_states, "Sanity check failed."

        return n_states

    def _compute_one_macrostate(
        self,
        n_cells: int | None,
        cluster_key: str | Mapping[str, str] | None,
        weight_key: str | None = None,
    ) -> None:
        _start = _time.perf_counter()
        logger.info("For 1 macrostate, stationary distribution is computed")

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
            weight_key=weight_key,
            check_row_sums=False,
            time=_start,
        )

    @d.dedent
    def _set_macrostates(
        self,
        memberships: np.ndarray,
        n_cells: int | None = 30,
        cluster_key: str | Mapping[str, str] = "clusters",
        weight_key: str | None = None,
        check_row_sums: bool = True,
        time: float | None = None,
        params: dict[str, Any] = types.MappingProxyType({}),
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
            Reference annotations to name and color the states; see :meth:`compute_macrostates`.
        weight_key
            Per-observation weights, only used when ``cluster_key`` points to
            :attr:`~anndata.AnnData.obsm`; see :meth:`compute_macrostates`.
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
            logger.debug("Setting the macrostates using macrostate assignment")
            assignment = pd.Series(np.argmax(memberships, axis=1).astype(str), dtype="category")
            # sometimes, a category can be missing
            assignment = assignment.cat.reorder_categories([str(i) for i in range(memberships.shape[1])])
            # fmt: on
        else:
            logger.debug("Setting the macrostates using macrostates memberships")

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

        assignment, colors = self._set_categorical_labels(assignment, cluster_key=cluster_key, weight_key=weight_key)
        memberships = Lineage(memberships, names=list(assignment.cat.categories), colors=colors)
        # fmt: on

        groups = assignment.value_counts()
        groups = groups[groups != n_cells].to_dict()
        if len(groups):
            logger.warning(
                "The following terminal states have different number of cells than requested (%s): %s",
                n_cells,
                groups,
            )

        self._write_macrostates(
            macrostates=assignment,
            colors=colors,
            memberships=memberships,
            time=time,
            params=params,
        )

    @log_writer
    @shadow
    def _write_macrostates(
        self,
        macrostates: pd.Series,
        colors: np.ndarray,
        memberships: Lineage,
        params: dict[str, Any] = types.MappingProxyType({}),
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
                logger.warning(
                    "Unable to compute the coarse-grained stationary distribution because the stationary "
                    "distribution of the transition matrix could not be calculated. This typically happens "
                    "for very large or (nearly) reducible transition matrices, where the underlying solver "
                    "fails to converge. Automatic initial state prediction via `.predict_initial_states()` "
                    "will not be available; set the initial states manually using `.set_initial_states()` "
                    "instead"
                )
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
                value=AnnData(tmat, obs=dists), copy=False,
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

    @log_writer
    @shadow
    def _write_states(
        self,
        which: Literal["initial", "terminal"],
        states: pd.Series | None,
        colors: np.ndarray | None,
        probs: pd.Series | None = None,
        memberships: Lineage | None = None,
        params: dict[str, Any] = types.MappingProxyType({}),
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

            self._set(obj=self._shadow_adata.uns, key=Key.uns.coarse(self.backward), value=tmat, copy=False)

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
