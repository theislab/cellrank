from typing import Any, Dict, Tuple, Union, Mapping, Optional, Sequence, NamedTuple
from typing_extensions import Literal

from types import MappingProxyType

from anndata import AnnData
from cellrank import logging as logg
from cellrank._utils._key import Key
from cellrank._utils._docs import d
from cellrank._utils._enum import _DEFAULT_BACKEND, Backend_t
from cellrank._utils._utils import (
    _pairwise,
    _process_series,
    _get_cat_and_null_indices,
    _calculate_lineage_absorption_time_means,
)
from cellrank._utils._lineage import Lineage
from cellrank._utils._linear_solver import _solve_lin_system
from cellrank.estimators.mixins._utils import (
    SafeGetter,
    BaseProtocol,
    StatesHolder,
    logger,
    shadow,
    register_plotter,
)

import numpy as np
import pandas as pd
from scipy.sparse import issparse, spmatrix
from pandas.api.types import infer_dtype, is_categorical_dtype

__all__ = ["AbsProbsMixin"]


class RecTransStates(NamedTuple):
    q: Union[np.ndarray, spmatrix]
    s: Union[np.ndarray, spmatrix]
    trans_indices: np.ndarray
    rec_indices: np.ndarray
    name_to_ixs: Dict[str, np.ndarray]
    term_states: pd.Series
    term_states_colors: Sequence[Any]


class AbsProbsProtocol(BaseProtocol):
    _term_states: StatesHolder

    @property
    def transition_matrix(self) -> Union[np.ndarray, spmatrix]:
        ...

    @property
    def terminal_states(self) -> pd.Series:
        ...

    @property
    def absorption_probabilities(self) -> Optional[Lineage]:
        ...

    @property
    def absorption_times(self) -> Optional[pd.DataFrame]:
        ...

    @property
    def priming_degree(self) -> Optional[pd.Series]:
        ...

    def __len__(self) -> int:
        ...

    def _compute_absorption_probabilities(
        self,
        q: Union[np.ndarray, spmatrix],
        s: Union[np.ndarray, spmatrix],
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
        ctx: Literal["abs_probs", "time_to_absorption"],
    ) -> RecTransStates:
        ...

    def _ensure_lineage_object(self, attr: str, **kwargs: Any) -> None:
        ...

    def _write_absorption_probabilities(
        self,
        abs_probs: Optional[Lineage],
    ) -> str:
        ...

    def _write_absorption_times(
        self,
        abs_times: Optional[pd.DataFrame],
        params: Mapping[str, Any] = MappingProxyType({}),
    ) -> str:
        ...

    def _write_lineage_priming(
        self,
        priming_degree: Optional[pd.Series],
    ) -> str:
        ...


def _normalize_abs_times(
    keys: Sequence[str], time_to_absorption: Any = None
) -> Dict[Tuple[str, ...], Literal["mean", "var"]]:
    if time_to_absorption is None:
        return {}

    if isinstance(time_to_absorption, (str, tuple)):
        time_to_absorption = [time_to_absorption]
    if not isinstance(time_to_absorption, dict):
        time_to_absorption = {ln: "mean" for ln in time_to_absorption}

    res = {}
    for ln, moment in time_to_absorption.items():
        if moment not in ("mean", "var"):
            raise ValueError(
                f"Moment must be either `'mean'` or `'var'`, found `{moment!r}` in `{ln}`."
            )

        seen = set()
        if isinstance(ln, str):
            ln = tuple(keys) if ln == "all" else (ln,)
        sorted_ln = tuple(sorted(ln))  # preserve the user order

        if sorted_ln not in seen:
            seen.add(sorted_ln)
            for lin in ln:
                if lin not in keys:
                    raise ValueError(
                        f"Invalid absorbing state `{lin!r}` in `{ln}`. "
                        f"Valid options are `{list(keys)}`."
                    )
            res[tuple(ln)] = moment

    return res


class AbsProbsMixin:
    """Mixin that supports computation of absorption probabilities and mean times to absorption."""

    def __init__(self, **kwargs: Any):
        super().__init__(**kwargs)

        self._absorption_probabilities: Optional[Lineage] = None
        self._absorption_times: Optional[pd.DataFrame] = None
        self._priming_degree: Optional[pd.Series] = None

    @property
    @d.get_summary(base="abs_probs")
    def absorption_probabilities(self) -> Optional[Lineage]:
        """Absorption probabilities.

        Informally, given a (finite, discrete) Markov chain with a set of transient states :math:`T` and
        a set of absorbing states :math:`A`, the absorption probability for cell :math:`i` from :math:`T`
        to reach cell :math:`j` from :math:`R` is the probability that a random walk initialized in :math:`i`
        will reach absorbing state :math:`j`.

        In our context, states correspond to cells, in particular, absorbing states correspond to cells in terminal
        states.
        """
        return self._absorption_probabilities

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
    def compute_absorption_probabilities(
        self: AbsProbsProtocol,
        keys: Optional[Sequence[str]] = None,
        solver: Union[
            str, Literal["direct", "gmres", "lgmres", "bicgstab", "gcrotmk"]
        ] = "gmres",
        use_petsc: bool = True,
        n_jobs: Optional[int] = None,
        backend: Backend_t = _DEFAULT_BACKEND,
        show_progress_bar: bool = True,
        tol: float = 1e-6,
        preconditioner: Optional[str] = None,
    ) -> None:
        """
        Compute absorption probabilities.

        For each cell, this computes the probability of being absorbed in any of the :attr:`terminal_states`. In
        particular, this corresponds to the probability that a random walk initialized in transient cell :math:`i`
        will reach any cell from a fixed transient state before reaching a cell from any other transient state.

        Parameters
        ----------
        keys
            Terminal states for which to compute the absorption probabilities.
            If `None`, use all states defined in :attr:`terminal_states`.
        %(absorption_utils)s

        Returns
        -------
        Nothing, just updates the following field:

            - :attr:`absorption_probabilities` - %(abs_probs.summary)s
        """
        start = logg.info("Computing absorption probabilities")
        data = self._rec_trans_states(keys, ctx="abs_probs")
        abs_probs = self._compute_absorption_probabilities(
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
        abs_probs = Lineage(
            abs_probs,
            names=list(data.term_states.cat.categories),
            colors=data.term_states_colors,
        )

        params = self._create_params(
            remove=["use_petsc", "n_jobs", "backend", "show_progress_bar"]
        )
        self._write_absorption_probabilities(
            abs_probs,
            params=params,
            time=start,
        )

    @d.dedent
    def compute_absorption_times(
        self: AbsProbsProtocol,
        keys: Optional[Sequence[str]] = None,
        calculate_variance: bool = False,
        solver: Union[
            str, Literal["direct", "gmres", "lgmres", "bicgstab", "gcrotmk"]
        ] = "gmres",
        use_petsc: bool = True,
        n_jobs: Optional[int] = None,
        backend: Backend_t = _DEFAULT_BACKEND,
        show_progress_bar: Optional[bool] = None,
        tol: float = 1e-6,
        preconditioner: Optional[str] = None,
    ) -> None:
        """
        Compute the mean time to absorption and optionally its variance.

        Parameters
        ----------
        keys
            Terminal states for which to compute the absorption probabilities.
            If `None`, use all states defined in :attr:`terminal_states`.
        calculate_variance
            Whether to calculate the variance.
        %(absorption_utils)s

        Returns
        -------
        Nothing, just updates the following field:

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

        params = self._create_params(
            remove=["use_petsc", "n_jobs", "backend", "show_progress_bar"]
        )
        self._write_absorption_times(
            abs_times,
            params=params,
            time=start,
        )

    @d.dedent
    def compute_lineage_priming(
        self: AbsProbsProtocol,
        method: Literal["kl_divergence", "entropy"] = "kl_divergence",
        early_cells: Optional[Union[Mapping[str, Sequence[str]], Sequence[str]]] = None,
    ) -> pd.Series:
        """
        %(lin_pd.full_desc)s

        Parameters
        ----------
        %(lin_pd.parameters)s
            If a :class:`dict`, the key specifies a cluster key in :attr:`anndata.AnnData.obs` and the values
            specify cluster labels containing early cells.

        Returns
        -------
        %(lin_pd.returns)s

        Also updates the following field:

            - :attr:`priming_degree` - %(priming_degree.summary)s
        """  # noqa: D400
        abs_probs = self.absorption_probabilities
        if abs_probs is None:
            raise RuntimeError(
                "Compute absorption probabilities first as `.compute_absorption_probabilities()`."
            )
        if isinstance(early_cells, dict):
            if len(early_cells) != 1:
                raise ValueError(
                    f"Expected a dictionary with only 1 key, found `{len(early_cells)}`."
                )
            key = next(iter(early_cells.keys()))
            if key not in self.adata.obs:
                raise KeyError(f"Unable to find clusters in `adata.obs[{key!r}]`.")
            if not is_categorical_dtype(self.adata.obs[key]):
                raise TypeError(
                    f"Expected `adata.obs[{key!r}]` to be categorical, "
                    f"found `{infer_dtype(self.adata.obs[key])}`."
                )
            early_cells = self.adata.obs[key].isin(early_cells[key])
        elif early_cells is not None:
            early_cells = np.asarray(early_cells)
            if not np.issubdtype(early_cells.dtype, np.bool_):
                early_cells = np.isin(self.adata.obs_names, early_cells)

        values = pd.Series(
            abs_probs.priming_degree(method, early_cells), index=self.adata.obs_names
        )
        self._write_lineage_priming(values)

        return values

    def _rec_trans_states(
        self: AbsProbsProtocol,
        keys: Optional[Sequence[str]] = None,
        *,
        ctx: Literal["abs_probs", "time_to_absorption"],
    ) -> RecTransStates:
        if self.terminal_states is None:
            raise RuntimeError(
                "Compute terminal states first as `.compute_terminal_states()`."
            )
        if keys is not None:
            keys = sorted(set(keys))

        # get the transition matrix
        if not issparse(self.transition_matrix):
            logg.warning(
                "Attempting to solve a potentially large linear system with dense transition matrix"
            )

        # process the current annotations according to `keys`
        term_states, colors = _process_series(
            series=self.terminal_states, keys=keys, colors=self._term_states_colors
        )
        # warn in case only one state is left
        keys = list(term_states.cat.categories)
        if ctx == "abs_probs" and len(keys) == 1:
            logg.warning(
                "There is only `1` terminal state, all cells will have probability `1` of going there"
            )

        # get indices corresponding to recurrent and transient states
        rec_indices, trans_indices, name_to_ixs = _get_cat_and_null_indices(term_states)
        if not len(trans_indices):
            raise RuntimeError("Markov chain is irreducible.")

        # create Q (restriction transient-transient), S (restriction transient-recurrent)
        q = self.transition_matrix[trans_indices, :][:, trans_indices]
        s = self.transition_matrix[trans_indices, :][:, rec_indices]

        # take individual solutions and piece them together to get absorption probabilities towards the classes
        macro_ix_helper = np.cumsum(
            [0] + [len(indices) for indices in name_to_ixs.values()]
        )
        # `s` can be sparse or dense, ensure the correct shape
        s = np.concatenate(
            [
                s[:, np.arange(a, b)].sum(axis=1).reshape(-1, 1)
                for a, b in _pairwise(macro_ix_helper)
            ],
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

    def _compute_absorption_probabilities(
        self: AbsProbsProtocol,
        q: Union[np.ndarray, spmatrix],
        s: Union[np.ndarray, spmatrix],
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
        abs_classes = np.zeros(
            shape=(len(self), len(term_states.cat.categories)), dtype=np.float64
        )
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
    def _write_absorption_probabilities(
        self: AbsProbsProtocol,
        abs_probs: Optional[Lineage],
        params: Mapping[str, Any] = MappingProxyType({}),
    ) -> str:
        # fmt: off
        key = Key.obsm.abs_probs(self.backward)
        self._set("_absorption_probabilities", self.adata.obsm, key=key, value=abs_probs)
        self._write_lineage_priming(None, log=False)
        self.params[key] = dict(params)
        # fmt: on

        return (
            f"Adding `adata.obsm[{key!r}]`\n"
            f"       `.absorption_probabilities`\n"
            f"    Finish"
        )

    @logger
    @shadow
    def _write_absorption_times(
        self: AbsProbsProtocol,
        abs_times: Optional[pd.DataFrame],
        params: Mapping[str, Any] = MappingProxyType({}),
    ) -> str:
        key = Key.obsm.abs_times(self.backward)
        self._set("_absorption_times", self.adata.obsm, key=key, value=abs_times)
        self.params[key] = dict(params)

        return (
            f"Adding `adata.obsm[{key!r}]`\n"
            f"       `.absorption_times`\n"
            f"    Finish"
        )

    @logger
    @shadow
    def _write_lineage_priming(
        self: AbsProbsProtocol, priming_degree: Optional[pd.Series]
    ) -> str:
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

    def _read_absorption_probabilities(self: AbsProbsProtocol, adata: AnnData) -> bool:
        # fmt: off
        with SafeGetter(self, allowed=KeyError) as sg:
            key1 = Key.obsm.abs_probs(self.backward)
            self._get("_absorption_probabilities", self.adata.obsm, key=key1, where="obsm", dtype=(np.ndarray, Lineage))
            self._ensure_lineage_object("_absorption_probabilities", kind="abs_probs")
            key = Key.obs.priming_degree(self.backward)
            self._get("_priming_degree", self.adata.obs, key=key, where="obs", dtype=pd.Series, allow_missing=True)
            self.params[key1] = self._read_params(key1)
        # fmt: on

        return sg.ok

    def _read_absorption_times(self: AbsProbsProtocol, adata: AnnData) -> bool:
        # fmt: off
        with SafeGetter(self, allowed=KeyError) as sg:
            key = Key.obsm.abs_times(self.backward)
            self._get("_absorption_times", self.adata.obsm, key=key, where="obsm", dtype=pd.DataFrame,
                      allow_missing=True)
            self.params[key] = self._read_params(key)
        # fmt: on

        return sg.ok

    def _ensure_lineage_object(
        self: AbsProbsProtocol, attr: str, **kwargs: Any
    ) -> None:
        if not isinstance(getattr(self, attr), Lineage):
            try:
                lineage = Lineage.from_adata(
                    self.adata, backward=self.backward, copy=True, **kwargs
                )
                setattr(self, attr, lineage)
            except Exception as e:  # noqa: B902
                raise RuntimeError(
                    f"Unable to reconstruct `.absorption_probabilities`. Reason: `{e}`."
                ) from None

    plot_absorption_probabilities = register_plotter(
        fwd_attr="absorption_probabilities"
    )
