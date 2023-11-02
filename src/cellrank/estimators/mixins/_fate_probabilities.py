import types
from typing import Any, Dict, Literal, Mapping, NamedTuple, Optional, Sequence, Union

import numpy as np
import pandas as pd
import scipy.sparse as sp
from pandas.api.types import infer_dtype

from anndata import AnnData

from cellrank import logging as logg
from cellrank._utils._docs import d
from cellrank._utils._enum import DEFAULT_BACKEND, Backend_t
from cellrank._utils._key import Key
from cellrank._utils._lineage import Lineage
from cellrank._utils._linear_solver import _solve_lin_system
from cellrank._utils._utils import (
    _calculate_lineage_absorption_time_means,
    _get_cat_and_null_indices,
    _pairwise,
    _process_series,
)
from cellrank.estimators.mixins._utils import (
    BaseProtocol,
    PlotMode,
    SafeGetter,
    StatesHolder,
    logger,
    shadow,
)

__all__ = ["FateProbsMixin"]


class RecTransStates(NamedTuple):
    q: Union[np.ndarray, sp.spmatrix]
    s: Union[np.ndarray, sp.spmatrix]
    trans_indices: np.ndarray
    rec_indices: np.ndarray
    name_to_ixs: Dict[str, np.ndarray]
    term_states: pd.Series
    term_states_colors: Sequence[Any]


class FateProbsProtocol(BaseProtocol):
    _term_states: StatesHolder

    @property
    def transition_matrix(self) -> Union[np.ndarray, sp.spmatrix]:
        ...

    @property
    def terminal_states(self) -> pd.Series:
        ...

    @property
    def fate_probabilities(self) -> Optional[Lineage]:
        ...

    @property
    def absorption_times(self) -> Optional[pd.DataFrame]:
        ...

    @property
    def priming_degree(self) -> Optional[pd.Series]:
        ...

    def __len__(self) -> int:
        ...

    def _compute_fate_probabilities(
        self,
        q: Union[np.ndarray, sp.spmatrix],
        s: Union[np.ndarray, sp.spmatrix],
        trans_indices: np.ndarray,
        term_states: np.ndim,
        solver: str,
        use_petsc: bool,
        n_jobs: Optional[int],
        backend: str,
        tol: float,
        show_progress_bar: bool,
        preconditioner: str,
    ) -> np.ndarray:
        ...

    def _rec_trans_states(
        self,
        keys: Optional[Sequence[str]],
        *,
        ctx: Literal["fate_probs", "time_to_absorption"],
    ) -> RecTransStates:
        ...

    def _ensure_lineage_object(
        self,
        obj: Union[str, np.ndarray, Lineage],
        *,
        kind: Literal["macrostates", "term_states", "fate_probs"],
        backward: bool,
        **kwargs: Any,
    ) -> Lineage:
        ...

    def _write_fate_probabilities(
        self,
        fate_probs: Optional[Lineage],
    ) -> str:
        ...

    def _write_absorption_times(
        self,
        abs_times: Optional[pd.DataFrame],
        params: Mapping[str, Any] = types.MappingProxyType({}),
    ) -> str:
        ...

    def _write_lineage_priming(
        self,
        priming_degree: Optional[pd.Series],
    ) -> str:
        ...

    def _plot_continuous(
        self,
        _data: Lineage,
        _colors: Optional[np.ndarray] = None,
        _title: Optional[str] = None,
        states: Optional[Union[str, Sequence[str]]] = None,
        color: Optional[str] = None,
        mode: Literal["embedding", "time"] = PlotMode.EMBEDDING,
        time_key: Optional[str] = None,
        title: Optional[Union[str, Sequence[str]]] = None,
        same_plot: bool = True,
        cmap: str = "viridis",
        **kwargs: Any,
    ) -> None:
        ...


class FateProbsMixin:
    """Mixin that supports computation of fate probabilities and mean times to absorption."""

    def __init__(self, **kwargs: Any):
        super().__init__(**kwargs)

        self._fate_probabilities: Optional[Lineage] = None
        self._absorption_times: Optional[pd.DataFrame] = None
        self._priming_degree: Optional[pd.Series] = None

    @property
    @d.get_summary(base="fate_probs")
    def fate_probabilities(self) -> Optional[Lineage]:
        """Fate probabilities.

        Informally, given a (finite, discrete) Markov chain with a set of transient states :math:`T` and
        a set of absorbing states :math:`A`, the absorption probability for cell :math:`i` from :math:`T`
        to reach cell :math:`j` from :math:`R` is the probability that a random walk initialized in :math:`i`
        will reach absorbing state :math:`j`.

        In our context, states correspond to cells, in particular, absorbing states correspond to cells
        in :attr:`terminal_states`.
        """
        return self._fate_probabilities

    @property
    @d.get_summary(base="abs_times")
    def absorption_times(self) -> Optional[pd.DataFrame]:
        """Mean and variance of the time until absorption.

        Related to conditional mean first passage times. Corresponds to the expectation of the time until absorption,
        depending on initialization, and the variance.
        """
        return self._absorption_times

    @property
    @d.get_summary(base="priming_degree")
    def priming_degree(self) -> Optional[pd.Series]:
        """Priming degree.

        Given a cell :math:`i` and a set of terminal states, this quantifies how committed vs. naive cell :math:`i` is,
        i.e. its degree of pluripotency. Low values correspond to naive cells (high degree of pluripotency), high values
        correspond to committed cells (low degree of pluripotency).
        """
        return self._priming_degree

    @d.dedent
    def compute_fate_probabilities(
        self: FateProbsProtocol,
        keys: Optional[Sequence[str]] = None,
        solver: Union[str, Literal["direct", "gmres", "lgmres", "bicgstab", "gcrotmk"]] = "gmres",
        use_petsc: bool = True,
        n_jobs: Optional[int] = None,
        backend: Backend_t = DEFAULT_BACKEND,
        show_progress_bar: bool = True,
        tol: float = 1e-6,
        preconditioner: Optional[str] = None,
    ) -> None:
        """Compute fate probabilities.

        For each cell, this computes the probability of being absorbed in any of the :attr:`terminal_states`. In
        particular, this corresponds to the probability that a random walk initialized in transient cell :math:`i`
        will reach any cell from a fixed transient state before reaching a cell from any other transient state.

        Parameters
        ----------
        keys
            Terminal states for which to compute the fate probabilities.
            If :obj:`None`, use all states defined in :attr:`terminal_states`.
        %(absorption_utils)s

        Returns
        -------
        Nothing, just updates the following fields:

        - :attr:`fate_probabilities` - %(fate_probs.summary)s
        """
        start = logg.info("Computing fate probabilities")
        data = self._rec_trans_states(keys, ctx="fate_probs")
        fate_probs = self._compute_fate_probabilities(
            data.q,
            data.s,
            trans_indices=data.trans_indices,
            term_states=data.term_states,
            solver=solver,
            use_petsc=use_petsc,
            n_jobs=n_jobs,
            backend=backend,
            tol=tol,
            show_progress_bar=show_progress_bar,
            preconditioner=preconditioner,
        )
        fate_probs = Lineage(
            fate_probs,
            names=list(data.term_states.cat.categories),
            colors=data.term_states_colors,
        )

        params = self._create_params(remove=["use_petsc", "n_jobs", "backend", "show_progress_bar"])
        self._write_fate_probabilities(
            fate_probs,
            params=params,
            time=start,
        )

    @d.dedent
    def plot_fate_probabilities(
        self: FateProbsProtocol,
        states: Optional[Union[str, Sequence[str]]] = None,
        color: Optional[str] = None,
        mode: Literal["embedding", "time"] = PlotMode.EMBEDDING,
        time_key: Optional[str] = None,
        same_plot: bool = True,
        title: Optional[Union[str, Sequence[str]]] = None,
        cmap: str = "viridis",
        **kwargs: Any,
    ) -> None:
        """Plot fate probabilities.

        Parameters
        ----------
        states
            Subset of the macrostates to show. If :obj:`None`, plot all macrostates.
        color
            Key in :attr:`~anndata.AnnData.obs` or :attr:`anndata.AnnData.var` used to color the observations.
        mode
            Whether to plot the probabilities in an embedding or along the pseudotime.
        time_key
            Key in :attr:`~anndata.AnnData.obs` where pseudotime is stored. Only used when ``mode = 'time'``.
        title
            Title of the plot.
        same_plot
            Whether to plot the data on the same plot or not. Only use when ``mode = 'embedding'``.
            If `True` and ``discrete = False``, ``color`` is ignored.
        cmap
            Colormap for continuous annotations.
        kwargs
            Keyword arguments for :func:`~scvelo.pl.scatter`.

        Returns
        -------
        %(just_plots)s
        """
        if self.fate_probabilities is None:
            raise RuntimeError("Compute fate probabilities first as `.compute_fate_probabilities()`.")

        return self._plot_continuous(
            _data=self.fate_probabilities,
            _colors=self.fate_probabilities.colors,
            _title="fate probabilities",
            states=states,
            color=color,
            mode=mode,
            time_key=time_key,
            same_plot=same_plot,
            title=title,
            cmap=cmap,
            **kwargs,
        )

    @d.dedent
    def compute_absorption_times(
        self: FateProbsProtocol,
        keys: Optional[Sequence[str]] = None,
        calculate_variance: bool = False,
        solver: Union[str, Literal["direct", "gmres", "lgmres", "bicgstab", "gcrotmk"]] = "gmres",
        use_petsc: bool = True,
        n_jobs: Optional[int] = None,
        backend: Backend_t = DEFAULT_BACKEND,
        show_progress_bar: Optional[bool] = None,
        tol: float = 1e-6,
        preconditioner: Optional[str] = None,
    ) -> None:
        """Compute the mean time to absorption and optionally its variance.

        Parameters
        ----------
        keys
            Terminal states for which to compute the fate probabilities.
            If :obj:`None`, use all states defined in :attr:`terminal_states`.
        calculate_variance
            Whether to calculate the variance.
        %(absorption_utils)s

        Returns
        -------
        Nothing, just updates the following fields:

        - :attr:`absorption_times` - %(abs_times.summary)s
        """
        start = logg.info("Computing absorption times")
        if show_progress_bar is None:
            # prevent from displaying too many progress bars
            show_progress_bar = not calculate_variance

        data = self._rec_trans_states(keys, ctx="time_to_absorption")
        abs_times = _calculate_lineage_absorption_time_means(
            data.q,
            data.s,
            calculate_variance=calculate_variance,
            trans_indices=data.trans_indices,
            ixs=data.name_to_ixs,
            solver=solver,
            use_petsc=use_petsc,
            n_jobs=n_jobs,
            backend=backend,
            tol=tol,
            show_progress_bar=show_progress_bar,
            preconditioner=preconditioner,
            index=self.adata.obs_names,
        )

        params = self._create_params(remove=["use_petsc", "n_jobs", "backend", "show_progress_bar"])
        self._write_absorption_times(
            abs_times,
            params=params,
            time=start,
        )

    @d.dedent
    def compute_lineage_priming(
        self: FateProbsProtocol,
        method: Literal["kl_divergence", "entropy"] = "kl_divergence",
        early_cells: Optional[Union[Mapping[str, Sequence[str]], Sequence[str]]] = None,
    ) -> pd.Series:
        """%(lin_pd.full_desc)s

        Parameters
        ----------
        %(lin_pd.parameters)s
            If a :class:`dict`, the key specifies a cluster key in :attr:`~anndata.AnnData.obs` and the values
            specify cluster labels containing early cells.

        Returns
        -------
        Returns the priming degree and updates the following fields:

        - :attr:`priming_degree` - %(priming_degree.summary)s
        """  # noqa: D400
        fate_probs = self.fate_probabilities
        if fate_probs is None:
            raise RuntimeError("Compute fate probabilities first as `.compute_fate_probabilities()`.")
        if isinstance(early_cells, dict):
            if len(early_cells) != 1:
                raise ValueError(f"Expected a dictionary with only 1 key, found `{len(early_cells)}`.")
            key = next(iter(early_cells.keys()))
            if key not in self.adata.obs:
                raise KeyError(f"Unable to find clusters in `adata.obs[{key!r}]`.")
            if not isinstance(self.adata.obs[key].dtype, pd.CategoricalDtype):
                raise TypeError(
                    f"Expected `adata.obs[{key!r}]` to be categorical, " f"found `{infer_dtype(self.adata.obs[key])}`."
                )
            early_cells = self.adata.obs[key].isin(early_cells[key])
        elif early_cells is not None:
            early_cells = np.asarray(early_cells)
            if not np.issubdtype(early_cells.dtype, np.bool_):
                early_cells = np.isin(self.adata.obs_names, early_cells)

        values = pd.Series(fate_probs.priming_degree(method, early_cells), index=self.adata.obs_names)
        self._write_lineage_priming(values)

        return values

    def _rec_trans_states(
        self: FateProbsProtocol,
        keys: Optional[Sequence[str]] = None,
        *,
        ctx: Literal["fate_probs", "time_to_absorption"],
    ) -> RecTransStates:
        if self.terminal_states is None:
            raise RuntimeError("Compute terminal states first as `.predict_terminal_states()`.")
        if keys is not None:
            keys = sorted(set(keys))

        # get the transition matrix
        if not sp.issparse(self.transition_matrix):
            logg.warning("Attempting to solve a potentially large linear system with dense transition matrix")

        # process the current annotations according to `keys`
        term_states, colors = _process_series(series=self.terminal_states, keys=keys, cols=self._term_states.colors)
        # warn in case only one state is left
        keys = list(term_states.cat.categories)
        if ctx == "fate_probs" and len(keys) == 1:
            logg.warning("There is only `1` terminal state, all cells will have probability `1` of going there")

        # get indices corresponding to recurrent and transient states
        rec_indices, trans_indices, name_to_ixs = _get_cat_and_null_indices(term_states)
        if not len(trans_indices):
            raise RuntimeError("Markov chain is irreducible.")

        # create Q (restriction transient-transient), S (restriction transient-recurrent)
        q = self.transition_matrix[trans_indices, :][:, trans_indices]
        s = self.transition_matrix[trans_indices, :][:, rec_indices]

        # take individual solutions and piece them together to get absorption probabilities towards the classes
        macro_ix_helper = np.cumsum([0] + [len(indices) for indices in name_to_ixs.values()])
        # `s` can be sparse or dense, ensure the correct shape
        s = np.concatenate(
            [s[:, np.arange(a, b)].sum(axis=1).reshape(-1, 1) for a, b in _pairwise(macro_ix_helper)],
            axis=1,
        )

        return RecTransStates(
            q=q,
            s=s,
            trans_indices=trans_indices,
            rec_indices=rec_indices,
            name_to_ixs=name_to_ixs,
            term_states=term_states,
            term_states_colors=colors,
        )

    def _compute_fate_probabilities(
        self: FateProbsProtocol,
        q: Union[np.ndarray, sp.spmatrix],
        s: Union[np.ndarray, sp.spmatrix],
        trans_indices: np.ndarray,
        term_states: np.ndim,
        solver: str,
        use_petsc: bool,
        n_jobs: Optional[int],
        backend: str,
        tol: float,
        show_progress_bar: bool,
        preconditioner: str,
    ) -> np.ndarray:
        _abs_classes = _solve_lin_system(
            q,
            s,
            solver=solver,
            use_petsc=use_petsc,
            n_jobs=n_jobs,
            backend=backend,
            tol=tol,
            use_eye=True,
            show_progress_bar=show_progress_bar,
            preconditioner=preconditioner,
        )
        abs_classes = np.zeros(shape=(len(self), len(term_states.cat.categories)), dtype=np.float64)
        for col, rec_class in enumerate(term_states.cat.categories):
            rec_indices = np.where(term_states == rec_class)[0]
            abs_classes[trans_indices, col] = _abs_classes[:, col]
            abs_classes[rec_indices, col] = 1.0

        mask = abs_classes >= 0
        if not np.all(mask):
            raise ValueError(
                f"`{np.sum(~mask)}` value(s) are negative. Try decreasing the tolerance as `tol=...`, "
                f"specifying a preconditioner as `preconditioner=...` or "
                f"use a direct solver as `solver='direct'` if the matrix is small."
            )
        mask = np.isclose(abs_classes.sum(1), 1.0, rtol=1e-3)
        if not np.all(mask):
            raise ValueError(
                f"`{np.sum(~mask)}` value(s) do not sum to 1 (rtol=1e-3). Try decreasing the tolerance as `tol=...`, "
                f"specifying a preconditioner as `preconditioner=...` or "
                f"use a direct solver as `solver='direct'` if the matrix is small."
            )

        return abs_classes

    @logger
    @shadow
    def _write_fate_probabilities(
        self: FateProbsProtocol,
        fate_probs: Optional[Lineage],
        params: Mapping[str, Any] = types.MappingProxyType({}),
    ) -> str:
        # fmt: off
        key = Key.obsm.fate_probs(self.backward)
        self._set("_fate_probabilities", self.adata.obsm, key=key, value=fate_probs)
        self._write_lineage_priming(None, log=False)
        self.params[key] = dict(params)
        # fmt: on

        return f"Adding `adata.obsm[{key!r}]`\n" f"       `.fate_probabilities`\n" f"    Finish"

    @logger
    @shadow
    def _write_absorption_times(
        self: FateProbsProtocol,
        abs_times: Optional[pd.DataFrame],
        params: Mapping[str, Any] = types.MappingProxyType({}),
    ) -> str:
        key = Key.obsm.abs_times(self.backward)
        self._set("_absorption_times", self.adata.obsm, key=key, value=abs_times)
        self.params[key] = dict(params)

        return f"Adding `adata.obsm[{key!r}]`\n       `.absorption_times`\n    Finish"

    @logger
    @shadow
    def _write_lineage_priming(self: FateProbsProtocol, priming_degree: Optional[pd.Series]) -> str:
        self._priming_degree = priming_degree
        key = Key.obs.priming_degree(self.backward)
        self._set("_priming_degree", self.adata.obs, key=key, value=priming_degree)

        # fmt: off
        return (
            f"Adding `adata.obs[{key!r}]`\n"
            f"       `.priming_degree`\n"
            f"    Finish"
        )
        # fmt: on

    def _read_fate_probabilities(self: FateProbsProtocol, adata: AnnData) -> bool:
        with SafeGetter(self, allowed=KeyError) as sg:
            key1 = Key.obsm.fate_probs(self.backward)
            fate_probs = self._get(
                obj=adata.obsm,
                key=key1,
                shadow_attr="obsm",
                dtype=(np.ndarray, Lineage),
            )
            self._fate_probabilities = self._ensure_lineage_object(
                fate_probs, backward=self.backward, kind="fate_probs"
            )
            key = Key.obs.priming_degree(self.backward)
            self._priming_degree = self._get(
                obj=adata.obs,
                key=key,
                shadow_attr="obs",
                dtype=pd.Series,
                allow_missing=True,
            )
            self.params[key1] = self._read_params(key1)

        return sg.ok

    def _read_absorption_times(self: FateProbsProtocol, adata: AnnData) -> bool:
        with SafeGetter(self, allowed=KeyError) as sg:
            key = Key.obsm.abs_times(self.backward)
            self._absorption_times = self._get(
                obj=adata.obsm,
                key=key,
                shadow_attr="obsm",
                dtype=pd.DataFrame,
                allow_missing=True,
            )
            self.params[key] = self._read_params(key)

        return sg.ok

    def _ensure_lineage_object(
        self: FateProbsProtocol,
        obj: Union[str, np.ndarray, Lineage],
        *,
        kind: Literal["macrostates", "term_states", "fate_probs"],
        backward: bool,
        **kwargs: Any,
    ) -> Lineage:
        if isinstance(obj, str):
            obj = getattr(self, obj)
        if isinstance(obj, Lineage):
            return obj

        try:
            return Lineage.from_adata(
                self.adata,
                backward=backward,
                kind=kind,
                estimator_backward=self.backward,
                copy=True,
                **kwargs,
            )
        except Exception as e:  # noqa: BLE001
            raise RuntimeError(f"Unable to reconstruct `.fate_probabilities`. Reason: `{e}`.") from None
