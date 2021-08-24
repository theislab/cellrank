from typing import Any, Dict, Tuple, Union, Optional, Protocol, Sequence

from datetime import datetime

from anndata import AnnData
from cellrank import logging as logg
from cellrank.tl import Lineage
from cellrank.ul._docs import inject_docs
from cellrank.tl._utils import (
    _pairwise,
    _process_series,
    _get_cat_and_null_indices,
    _calculate_lineage_absorption_time_means,
)
from cellrank.tl._constants import _colors, _lin_names
from cellrank.tl._linear_solver import _solve_lin_system
from cellrank.tl._estimators.mixins._constants import Key

import numpy as np
import pandas as pd
from scipy.sparse import issparse, csr_matrix


class AbsProbsProtocol(Protocol):
    @property
    def transition_matrix(self) -> Union[np.ndarray, csr_matrix]:
        ...

    @property
    def adata(self) -> AnnData:
        ...

    @property
    def obs_names(self) -> pd.Index:
        ...

    @property
    def terminal_states(self) -> pd.Series:
        ...

    @property
    def backward(self) -> bool:
        ...

    @property
    def absorption_probabilities(self) -> Optional[Lineage]:
        ...

    @property
    def absorption_times(self) -> Optional[pd.DataFrame]:
        ...

    def _write_absorption_probabilities(self, time: datetime) -> None:
        ...

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

    def _terminal_states_colors(self) -> pd.Series:
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


class AbsProbMixin:
    def __init__(self):
        self._absorption_probabilities: Optional[Lineage] = None
        self._absorption_times: Optional[pd.DataFrame] = None

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
        check_irreducibility
            Check whether the transition matrix is irreducible.
        solver
            Solver to use for the linear problem. Options are `'direct', 'gmres', 'lgmres', 'bicgstab' or 'gcrotmk'`
            when ``use_petsc=False`` or one of :class:`petsc4py.PETSc.KPS.Type` otherwise.

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

            It might be beneficial to disable the progress bar as ``show_progress_bar=False``, because many linear
            systems are being solved.
        n_jobs
            Number of parallel jobs to use when using an iterative solver. When ``use_petsc=True`` or for
            quickly-solvable problems, we recommend higher number (>=8) of jobs in order to fully saturate the cores.
        backend
            Which backend to use for multiprocessing. See :class:`joblib.Parallel` for valid options.
        show_progress_bar
            Whether to show progress bar. Only used when the solver isn't a direct one.
        tol
            Convergence tolerance for the iterative solver. The default is fine for most cases, only consider
            decreasing this for severely ill-conditioned matrices.
        preconditioner
            Preconditioner to use, only available when ``use_petsc=True``. For available values, see
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
            series=self.terminal_states, keys=keys, colors=self._terminal_states_colors
        )
        # warn in case only one state is left
        keys = list(term_states.cat.categories)
        if len(keys) == 1:
            logg.warning(
                "There is only 1 recurrent class, all cells will have probability 1 of going there"
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

        abs_classes = self._compute_absorption_probabilities(
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
                index=self.obs_names,
            )

        self._absorption_probabilities = Lineage(abs_classes, names=keys, colors=colors)
        self._absorption_times = abs_times
        self._write_absorption_probabilities(time=start)

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

    def _write_absorption_probabilities(self: AbsProbsProtocol, time: datetime) -> None:
        key = Key.obsm.abs_probs(self.backward)
        self.adata.obsm[key] = self.absorption_probabilities
        self.adata.uns[Key.uns.names(key)] = self.absorption_probabilities.names
        self.adata.uns[Key.uns.colors(key)] = self.absorption_probabilities.colors

        # extra_msg = ""  # TODO
        if self.lineage_absorption_times is not None:
            self.adata.obsm[Key.obsm.abs_times(key)] = self.absorption_times

        logg.info(
            f"Adding `adata.obsm[{key!r}]`\n"
            f"       `.absorption_probabilities`\n"
            "    Finish",
            time=time,
        )

    @property
    def absorption_probabilities(self) -> Optional[Lineage]:
        """TODO."""
        return self._absorption_probabilities

    @property
    def absorption_times(self) -> Optional[pd.DataFrame]:
        """TODO."""
        return self._absorption_times
