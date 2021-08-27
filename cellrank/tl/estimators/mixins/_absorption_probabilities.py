from typing import Any, Dict, Tuple, Union, Mapping, Optional, Sequence
from typing_extensions import Literal, Protocol

from types import MappingProxyType

from anndata import AnnData
from cellrank import logging as logg
from cellrank.tl._key import Key
from cellrank.ul._docs import d
from cellrank.tl._utils import (
    _pairwise,
    _process_series,
    _get_cat_and_null_indices,
    _calculate_lineage_absorption_time_means,
)
from cellrank.tl._lineage import Lineage
from cellrank.tl._linear_solver import _solve_lin_system
from cellrank.tl.estimators._utils import SafeGetter
from cellrank.tl.estimators.mixins._utils import logger, shadow, register_plotter

import numpy as np
import pandas as pd
from scipy.sparse import issparse, csr_matrix


class AbsProbsProtocol(Protocol):
    _term_states_colors: np.ndarray
    priming_degree: pd.Series

    @property
    def adata(self) -> AnnData:
        ...

    @property
    def backward(self) -> bool:
        ...

    @property
    def params(self) -> Dict[str, Any]:
        ...

    @property
    def transition_matrix(self) -> Union[np.ndarray, csr_matrix]:
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

    def _write_absorption_probabilities(
        self,
        abs_probs: Optional[Lineage],
        abs_times: Optional[pd.DataFrame],
    ) -> str:
        ...

    def _write_lineage_priming(
        self,
        priming_degree: Optional[pd.Series],
    ) -> str:
        ...

    # TODO: type
    def _compute_absorption_probabilities(
        self,
        q,
        s,
        trans_indices,
        term_states,
        solver,
        use_petsc,
        n_jobs,
        backend,
        tol,
        show_progress_bar,
        preconditioner,
    ) -> np.ndarray:
        ...

    def _set(
        self,
        attr: str,
        obj: Optional[Union[pd.DataFrame, Mapping[str, Any]]] = None,
        key: Optional[str] = None,
        value: Optional[Union[np.ndarray, pd.Series, pd.DataFrame, Lineage]] = None,
    ):
        ...

    def __len__(self) -> int:
        ...


def _normalize_abs_times(
    keys: Sequence[str], time_to_absorption: Any = None
) -> Dict[Tuple[str, ...], str]:
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
                f"Moment must be either `'mean'` or `'var'`, found `{moment!r}` for `{ln!r}`."
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
                        f"Invalid absorbing state `{lin!r}` for `{ln}`. "
                        f"Valid options are `{list(keys)}`."
                    )
            res[tuple(ln)] = moment

    return res


class AbsProbsMixin:
    def __init__(self, **kwargs: Any):
        super().__init__(**kwargs)

        self._absorption_probabilities: Optional[Lineage] = None
        self._absorption_times: Optional[pd.DataFrame] = None
        self._priming_degree: Optional[pd.Series] = None

    def compute_absorption_probabilities(
        self: AbsProbsProtocol,
        keys: Optional[Sequence[str]] = None,
        # check_irreducibility: bool = False,  # TODO
        solver: str = "gmres",
        use_petsc: bool = True,
        time_to_absorption: Optional[
            Union[
                str,
                Sequence[Union[str, Sequence[str]]],
                Dict[Union[str, Sequence[str]], str],
            ]
        ] = None,
        n_jobs: Optional[int] = None,
        backend: str = "loky",
        show_progress_bar: bool = True,
        tol: float = 1e-6,
        preconditioner: Optional[str] = None,
    ) -> None:
        """
        Compute absorption probabilities.

        For each cell, this computes the probability of being absorbed in the by :attr:`terminal_states`.

        Parameters
        ----------
        keys
            Keys defining the recurrent classes.
        solver
            Solver to use for the linear problem. Options are `'direct', 'gmres', 'lgmres', 'bicgstab' or 'gcrotmk'`
            when ``use_petsc = False`` or one of :class:`petsc4py.PETSc.KPS.Type` otherwise.

            Information on the :mod:`scipy` iterative solvers can be found in :func:`scipy.sparse.linalg` or for
            :mod:`petsc4py` solver `here <https://www.mcs.anl.gov/petsc/documentation/linearsolvertable.html>`__.
        use_petsc
            Whether to use solvers from :mod:`petsc4py` or :mod:`scipy`. Recommended for large problems.
            If no installation is found, defaults to :func:`scipy.sparse.linalg.gmres`.
        time_to_absorption
            Whether to compute mean time to absorption and its variance to specific absorbing states.

            If a :class:`dict`, can be specified as ``{{'Alpha': 'var', ...}}`` to also compute variance.
            In case when states are a :class:`tuple`, time to absorption will be computed to the subset of these states,
            such as ``[('Alpha', 'Beta'), ...]`` or ``{{('Alpha', 'Beta'): 'mean', ...}}``.
            Can be specified as ``'all'`` to compute it to any absorbing state in ``keys``, which is more efficient
            than listing all absorbing states.

            It might be beneficial to disable the progress bar as ``show_progress_bar = False``.
        n_jobs
            Number of parallel jobs to use when using an iterative solver. When ``use_petsc = True`` or for
            quickly-solvable problems, we recommend higher number (>=8) of jobs in order to fully saturate the cores.
        backend
            Which backend to use for multiprocessing. See :class:`joblib.Parallel` for valid options.
        show_progress_bar
            Whether to show progress bar. Only used when the solver isn't a direct one.
        tol
            Convergence tolerance for the iterative solver. The default is fine for most cases, only consider
            decreasing this for severely ill-conditioned matrices.
        preconditioner
            Preconditioner to use, only available when ``use_petsc = True``. For valid options, see
            `here <https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCType.html#PCType>`__ or the values
            of `petsc4py.PETSc.PC.Type`.
            We recommended `'ilu'` preconditioner for badly conditioned problems.

        Returns
        -------
        Nothing, just updates the following fields:

            - :attr:`absorption_probabilities` - probabilities of being absorbed into the :attr:`terminal_states`.
            - :attr:`absorption_times` - mean times until absorption to subsets of absorbing states and
              optionally their variances saved as ``'{{lineage}} mean'`` and ``'{{lineage}} var'``, respectively.
        """
        if self.terminal_states is None:
            raise RuntimeError(
                "Compute `.terminal_states` first as `.compute_terminal_states()`."
            )
        if keys is not None:
            keys = sorted(set(keys))

        start = logg.info("Computing absorption probabilities")

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
        if len(keys) == 1:
            logg.warning(
                "There is only 1 terminal state, all cells will have probability 1 of going there"
            )

        # get indices corresponding to recurrent and transient states
        rec_indices, trans_indices, lookup_dict = _get_cat_and_null_indices(term_states)
        if not len(trans_indices):
            raise RuntimeError("Markov chain is irreducible.")

        # create Q (restriction transient-transient), S (restriction transient-recurrent)
        q = self.transition_matrix[trans_indices, :][:, trans_indices]
        s = self.transition_matrix[trans_indices, :][:, rec_indices]

        # take individual solutions and piece them together to get absorption probabilities towards the classes
        # fmt: off
        macro_ix_helper = np.cumsum([0] + [len(indices) for indices in lookup_dict.values()])
        s = np.concatenate([s[:, np.arange(a, b)].sum(axis=1) for a, b in _pairwise(macro_ix_helper)], axis=1)
        # fmt: on

        abs_probs = self._compute_absorption_probabilities(
            q=q,
            s=s,
            trans_indices=trans_indices,
            term_states=term_states,
            solver=solver,
            use_petsc=use_petsc,
            n_jobs=n_jobs,
            backend=backend,
            tol=tol,
            show_progress_bar=show_progress_bar,
            preconditioner=preconditioner,
        )
        abs_probs = Lineage(abs_probs, names=keys, colors=colors)

        abs_times = None
        if time_to_absorption is not None:
            lineages = _normalize_abs_times(keys, time_to_absorption=time_to_absorption)
            abs_times = _calculate_lineage_absorption_time_means(
                q,
                s,
                trans_indices=trans_indices,
                ixs=lookup_dict,
                lineages=lineages,
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
        self._write_absorption_probabilities(
            abs_probs, abs_times, params=params, time=start
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
                raise KeyError(f"Unable to find clustering in `adata.obs[{key!r}]`.")
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

    # TODO(michalk8): type
    def _compute_absorption_probabilities(
        self: AbsProbsProtocol,
        q,
        s,
        trans_indices,
        term_states,
        solver,
        use_petsc,
        n_jobs,
        backend,
        tol,
        show_progress_bar,
        preconditioner,
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
        abs_classes = np.zeros((len(self), len(term_states.cat.categories)))
        rec_classes_full = {
            cl: np.where(term_states == cl)[0] for cl in term_states.cat.categories
        }
        for col, cl_indices in enumerate(rec_classes_full.values()):
            abs_classes[trans_indices, col] = _abs_classes[:, col]
            abs_classes[cl_indices, col] = 1

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
        abs_times: Optional[pd.DataFrame],
        params: Mapping[str, Any] = MappingProxyType({}),
    ) -> str:
        # fmt: off
        key1 = Key.obsm.abs_probs(self.backward)
        self._set("_absorption_probabilities", self.adata.obsm, key=key1, value=abs_probs)
        key2 = Key.obsm.abs_times(self.backward)
        self._set("_absorption_times", self.adata.obsm, key=key2, value=abs_times)
        self._write_lineage_priming(None, log=False)
        self.params[key1] = dict(params)
        # fmt: on

        if abs_times is None:
            return (
                f"Adding `adata.obsm[{key1!r}]`\n"
                f"       `.absorption_probabilities`\n"
                "    Finish"
            )
        return (
            f"Adding `adata.obsm[{key1!r}]`\n"
            f"       `adata.obsm[{key2!r}]`\n"
            f"       `.absorption_probabilities`\n"
            f"       `.absorption_times`\n"
            "    Finish"
        )

    @logger
    @shadow
    def _write_lineage_priming(
        self: AbsProbsProtocol, priming_degree: Optional[pd.Series]
    ) -> str:
        self._priming_degree = priming_degree
        key = Key.obs.priming_degree(self.backward)
        self._set("_priming_degree", self.adata.obs, key=key, value=priming_degree)

        return f"Adding `adata.obs[{key!r}]`\n       `.priming_degree`\n    Finish"

    def _read_absorption_probabilities(
        self: AbsProbsProtocol, anndata: AnnData
    ) -> bool:
        # fmt: off
        with SafeGetter(self, allowed=KeyError) as sg:
            key1 = Key.obsm.abs_probs(self.backward)
            self._get("_absorption_probabilities", self.adata.obsm, key=key1, where="obsm", dtype=Lineage)
            key = Key.obsm.abs_times(self.backward)
            self._get("_absorption_times", self.adata.obsm, key=key, where="obsm", dtype=pd.DataFrame,
                      allow_missing=True)
            key = Key.obs.priming_degree(self.backward)
            self._get("_priming_degree", self.adata.obs, key=key, where="obs", dtype=pd.Series, allow_missing=True)
            self.params[key1] = self._read_params(key1)
        # fmt: on

        return sg.ok

    plot_absorption_probabilities = register_plotter(
        continuous="absorption_probabilities"
    )

    @property
    def absorption_probabilities(self) -> Optional[Lineage]:
        """TODO."""
        return self._absorption_probabilities

    @property
    def absorption_times(self) -> Optional[pd.DataFrame]:
        """TODO."""
        return self._absorption_times

    @property
    def priming_degree(self) -> Optional[pd.Series]:
        """TODO."""
        return self._priming_degree
