"""Abstract base class for all kernel-holding estimators."""
import pickle
from abc import ABC, abstractmethod
from sys import version_info
from copy import copy, deepcopy
from typing import (
    Any,
    Dict,
    Tuple,
    Union,
    Mapping,
    TypeVar,
    Iterable,
    Optional,
    Sequence,
)
from pathlib import Path
from datetime import datetime

from typing_extensions import Literal

import scanpy as sc
import scvelo as scv
from cellrank import logging as logg
from cellrank.ul._docs import d, inject_docs
from cellrank.tl._utils import (
    TestMethod,
    save_fig,
    _pairwise,
    _irreducible,
    _process_series,
    _correlation_test,
    _get_cat_and_null_indices,
    _merge_categorical_series,
    _convert_to_categorical_series,
    _calculate_lineage_absorption_time_means,
)
from cellrank.tl._colors import (
    _map_names_and_colors,
    _convert_to_hex_colors,
    _create_categorical_colors,
)
from cellrank.tl.kernels import PrecomputedKernel
from cellrank.tl._lineage import Lineage
from cellrank.tl._constants import (
    Macro,
    DirPrefix,
    AbsProbKey,
    PrettyEnum,
    TermStatesKey,
    _pd,
    _probs,
    _colors,
    _lin_names,
)
from cellrank.tl._linear_solver import _solve_lin_system
from cellrank.tl.estimators._property import Partitioner, LineageEstimatorMixin
from cellrank.tl.kernels._base_kernel import KernelExpression
from cellrank.tl.estimators._constants import A, P

import numpy as np
import pandas as pd
from pandas import Series
from scipy.stats import ranksums
from scipy.sparse import spmatrix
from pandas.api.types import infer_dtype, is_categorical_dtype

import matplotlib.pyplot as plt
from matplotlib import rc_context, patheffects
from matplotlib.axes import Axes
from matplotlib.colors import is_color_like
from matplotlib.patches import ArrowStyle

AnnData = TypeVar("AnnData")
_COMP_TERM_STATES_MSG = (
    "Compute terminal states first as `.compute_terminal_states()` or"
    "set them manually as `.set_terminal_states()`."
)


@d.get_sections(base="base_estimator", sections=["Parameters"])
@d.dedent
class BaseEstimator(LineageEstimatorMixin, Partitioner, ABC):
    """
    Base class for all estimators.

    Parameters
    ----------
    obj
        Either a :class:`cellrank.tl.kernels.Kernel` object, an :class:`anndata.AnnData` object which
        stores the transition matrix in ``.obsp`` attribute or :mod:`numpy` or :mod:`scipy` array.
    inplace
        Whether to modify :attr:`adata` object inplace or make a copy.
    read_from_adata
        Whether to read available attributes in :attr:`adata`, if present.
    obsp_key
        Key in ``obj.obsp`` when ``obj`` is an :class:`anndata.AnnData` object.
    g2m_key
        Key in :attr:`adata` ``.obs``. Can be used to detect cell-cycle driven start- or endpoints.
    s_key
        Key in :attr:`adata` ``.obs``. Can be used to detect cell-cycle driven start- or endpoints.
    write_to_adata
        Whether to write the transition matrix to :attr:`adata` ``.obsp`` and the parameters to
        :attr:`adata` ``.uns``.
    %(write_to_adata.parameters)s
        Only used when ``write_to_adata=True``.
    """

    def __init__(
        self,
        obj: Union[KernelExpression, AnnData, spmatrix, np.ndarray],
        inplace: bool = True,
        read_from_adata: bool = False,
        obsp_key: Optional[str] = None,
        g2m_key: Optional[str] = "G2M_score",
        s_key: Optional[str] = "S_score",
        write_to_adata: bool = True,
        key: Optional[str] = None,
    ):
        from anndata import AnnData

        super().__init__(obj, obsp_key=obsp_key, key=key, write_to_adata=write_to_adata)

        if isinstance(obj, (KernelExpression, AnnData)) and not inplace:
            self.kernel._adata = self.adata.copy()

        if self.kernel.backward:
            self._term_key = TermStatesKey.BACKWARD.s
            self._abs_prob_key = AbsProbKey.BACKWARD.s
            self._term_abs_prob_key = Macro.BACKWARD.s
        else:
            self._term_key = TermStatesKey.FORWARD.s
            self._abs_prob_key = AbsProbKey.FORWARD.s
            self._term_abs_prob_key = Macro.FORWARD.s

        self._key_added = key
        self._g2m_key = g2m_key
        self._s_key = s_key

        self._G2M_score = None
        self._S_score = None

        self._absorption_time_mean = None
        self._absorption_time_var = None

        if read_from_adata:
            self._read_from_adata()

    def __init_subclass__(cls, **kwargs: Any):
        super().__init_subclass__()

    def _read_from_adata(self) -> None:
        def _check_term_states():
            term_states = self._get(A.TERM)
            colors = self._get(A.TERM_COLORS)

            if term_states is None or colors is None:
                return
            if not is_categorical_dtype(term_states):
                raise TypeError(
                    f"Expected `adata.obs[{self._term_key!r}]` to be `categorical`, "
                    f"found `{infer_dtype(term_states)}`."
                )
            if len(term_states.cat.categories) != len(colors):
                raise ValueError(
                    f"Expected `adata.uns[{_colors(self._term_key)}]` to have "
                    f"`{len(term_states.cat.categories)}` colors, found `{len(colors)}`."
                )

        self._set_or_debug(f"eig_{self._direction}", self.adata.uns, "_eig")

        self._set_or_debug(self._g2m_key, self.adata.obs, "_G2M_score")
        self._set_or_debug(self._s_key, self.adata.obs, "_S_score")

        self._set_or_debug(self._term_key, self.adata.obs, A.TERM.s)
        self._set_or_debug(_probs(self._term_key), self.adata.obs, A.TERM_PROBS)
        self._set_or_debug(_colors(self._term_key), self.adata.uns, A.TERM_COLORS)
        _check_term_states()

        self._reconstruct_lineage(A.ABS_PROBS, self._abs_prob_key)
        self._set_or_debug(_pd(self._abs_prob_key), self.adata.obs, A.PRIME_DEG)

    def _reconstruct_lineage(self, attr: PrettyEnum, obsm_key: str):
        self._set_or_debug(obsm_key, self.adata.obsm, attr)
        probs = self._get(attr)

        if probs is not None:
            names = self._set_or_debug(_lin_names(self._term_key), self.adata.uns)
            colors = self._set_or_debug(_colors(self._term_key), self.adata.uns)
            if len(names) != len(colors):
                raise ValueError(
                    f"Expected names and colors to be of same length, found `{names}` and `{colors}`."
                )

            if len(names) != probs.shape[1]:
                if isinstance(probs, Lineage):
                    names = probs.names
                else:
                    logg.warning(
                        f"Expected lineage names to be of length `{probs.shape[1]}`, found `{len(names)}`. "
                        f"Creating new names"
                    )
                    names = [f"Lineage {i}" for i in range(probs.shape[1])]
            if len(colors) != probs.shape[1] or not all(
                is_color_like(c) for c in colors
            ):
                if isinstance(probs, Lineage):
                    colors = probs.colors
                else:
                    logg.warning(
                        f"Expected lineage colors to be of length `{probs.shape[1]}`, found `{len(colors)}`. "
                        f"Creating new colors"
                    )
                    colors = _create_categorical_colors(probs.shape[1])

            colors = _convert_to_hex_colors(colors)
            self._set(attr, Lineage(probs, names=names, colors=colors))

            # ensure that the cont. lineage names match with the categorical values
            term_states = self._get(A.TERM)
            if term_states is not None:
                direction = "initial" if self.kernel.backward else "terminal"
                if len(names) != len(term_states.cat.categories):
                    self._set(attr, None)
                    logg.warning(
                        f"Expected to find `{len(names)}` {direction} states "
                        f"(from adata.uns[{_lin_names(self._term_key)!r}]), "
                        f"found `{len(term_states.cat.categories)}` (from adata.obsm[{obsm_key!r}]). Skipping`.{attr}`"
                    )
                    return
                if tuple(term_states.cat.categories) != tuple(names):
                    self._set(attr, None)
                    logg.warning(
                        f"Expected to find `{names}` {direction} states "
                        f"(from adata.uns[{_lin_names(self._term_key)!r}]), "
                        f"found `{term_states.cat.categories}` (from adata.obs[{obsm_key!r}). Skipping `{attr}`"
                    )
                    return

            self.adata.obsm[obsm_key] = self._get(attr)
            self.adata.uns[_lin_names(self._term_key)] = names
            self.adata.uns[_colors(self._term_key)] = colors

    @d.dedent
    @inject_docs(fs=P.TERM.s, fsp=P.TERM_PROBS.s)
    def set_terminal_states(
        self,
        labels: Union[Series, Dict[str, Sequence[Any]]],
        cluster_key: Optional[str] = None,
        en_cutoff: Optional[float] = None,
        p_thresh: Optional[float] = None,
        add_to_existing: bool = False,
        **kwargs,
    ) -> None:
        """
        Manually define terminal states.

        Parameters
        ----------
        labels
            Defines the terminal states. Valid options are:

                - categorical :class:`pandas.Series` where each category corresponds to one terminal state.
                  `NaN` entries denote cells that do not belong to any terminal state, i.e. these are either initial or
                  transient cells.
                - :class:`dict` where keys are terminal states and values are lists of cell barcodes corresponding to
                  annotations in :attr:`adata` ``.obs_names``.
                  If only 1 key is provided, values should correspond to terminal state clusters if a categorical
                  :class:`pandas.Series` can be found in :attr:`adata` ``.obs``.

        cluster_key
            Key from :attr:`adata.obs` where categorical cluster labels are stored. These are used to associate names
            and colors with each terminal state. Each terminal state will be given the name and color corresponding to
            the cluster it mostly overlaps with.
        %(en_cutoff_p_thresh)s
        add_to_existing
            Whether the new terminal states should be added to pre-existing ones. Cells already assigned to a terminal
            state will be re-assigned to the new terminal state if there's a conflict between old and new annotations.
            This throws an error if no previous annotations corresponding to terminal states have been found.

        Returns
        -------
        None
            Nothing, but updates the following fields:

                - :attr:`{fsp}`
                - :attr:`{fs}`
        """

        self._set_categorical_labels(
            attr_key=A.TERM.v,
            color_key=A.TERM_COLORS.v,
            pretty_attr_key=P.TERM.v,
            add_to_existing_error_msg=_COMP_TERM_STATES_MSG,
            categories=labels,
            cluster_key=cluster_key,
            en_cutoff=en_cutoff,
            p_thresh=p_thresh,
            add_to_existing=add_to_existing,
        )
        self._write_terminal_states(time=kwargs.get("time", None))

    @inject_docs(ts=P.TERM.s)
    def rename_terminal_states(
        self, new_names: Mapping[str, str], update_adata: bool = True
    ) -> None:
        """
        Rename the names of :attr:`{ts}`.

        Parameters
        ----------
        new_names
            Mapping where keys are the old names and the values are the new names. New names must be unique.
        update_adata
            Whether to update underlying :attr:`adata` object as well or not.

        Returns
        -------
        None
            Nothing, just updates the names of :attr:`{ts}`.
        """

        term_states = self._get(P.TERM)

        if term_states is None:
            raise RuntimeError(_COMP_TERM_STATES_MSG)

        if not isinstance(new_names, Mapping):
            raise TypeError(f"Expected a `Mapping` type, found `{type(new_names)!r}`.")
        if not len(new_names):
            return

        new_names = {k: str(v) for k, v in new_names.items()}

        mask = np.isin(list(new_names.keys()), term_states.cat.categories)
        if not np.all(mask):
            raise ValueError(
                f"Invalid old terminal states names: `{np.array(list(new_names.keys()))[~mask]}`."
            )

        names_after_renaming = [new_names.get(n, n) for n in term_states.cat.categories]
        if len(set(names_after_renaming)) != len(term_states.cat.categories):
            raise ValueError(
                f"After renaming, the names will not be unique: `{names_after_renaming}`."
            )

        term_states.cat.rename_categories(new_names, inplace=True)

        memberships = (
            self._get(A.TERM_ABS_PROBS) if hasattr(self, A.TERM_ABS_PROBS.s) else None
        )
        if memberships is not None:  # GPCCA
            memberships.names = [new_names.get(n, n) for n in memberships.names]
            self._set(A.TERM_ABS_PROBS, memberships)

        # we can be just computing it and it's not yet saved in adata
        if (
            update_adata
            and self._term_key in self.adata.obs
            and _lin_names(self._term_key) in self.adata.uns
        ):
            self.adata.obs[self._term_key].cat.rename_categories(
                new_names, inplace=True
            )
            self.adata.uns[_lin_names(self._term_key)] = np.array(
                self.adata.obs[self._term_key].cat.categories
            )

    @inject_docs(
        abs_prob=P.ABS_PROBS,
        fs=P.TERM.s,
        lat=P.LIN_ABS_TIMES,
    )
    def compute_absorption_probabilities(
        self,
        keys: Optional[Sequence[str]] = None,
        check_irreducibility: bool = False,
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
        Compute absorption probabilities of a Markov chain.

        For each cell, this computes the probability of it reaching any of the approximate recurrent classes defined
        by :attr:`{fs}`.

        Parameters
        ----------
        keys
            Keys defining the recurrent classes.
        check_irreducibility:
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
            Whether to show progress bar when the solver isn't a direct one.
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
        None
            Nothing, but updates the following fields:

                - :attr:`{abs_prob}` - probabilities of being absorbed into the terminal states.
                - :attr:`{lat}` - mean times until absorption to subset absorbing states and optionally
                  their variances saved as ``'{{lineage}} mean'`` and ``'{{lineage}} var'``, respectively,
                  for each subset of absorbing states specified in ``time_to_absorption``.
        """

        if self._get(P.TERM) is None:
            raise RuntimeError(_COMP_TERM_STATES_MSG)
        if keys is not None:
            keys = sorted(set(keys))

        start = logg.info("Computing absorption probabilities")

        # get the transition matrix
        t = self.transition_matrix
        if not self.issparse:
            logg.warning(
                "Attempting to solve a potentially large linear system with dense transition matrix"
            )

        # process the current annotations according to `keys`
        terminal_states_, colors_ = _process_series(
            series=self._get(P.TERM), keys=keys, colors=self._get(A.TERM_COLORS)
        )
        # warn in case only one state is left
        keys = list(terminal_states_.cat.categories)
        if len(keys) == 1:
            logg.warning(
                "There is only 1 recurrent class, all cells will have probability 1 of going there"
            )

        lin_abs_times = {}
        if time_to_absorption is not None:
            if isinstance(time_to_absorption, (str, tuple)):
                time_to_absorption = [time_to_absorption]
            if not isinstance(time_to_absorption, dict):
                time_to_absorption = {ln: "mean" for ln in time_to_absorption}
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
                                f"Invalid absorbing state `{lin!r}` in `{ln}`. "
                                f"Valid options are `{list(terminal_states_.cat.categories)}`."
                            )
                    lin_abs_times[tuple(ln)] = moment

            # define the dimensions of this problem
        n_cells = t.shape[0]
        n_macrostates = len(terminal_states_.cat.categories)

        # get indices corresponding to recurrent and transient states
        rec_indices, trans_indices, lookup_dict = _get_cat_and_null_indices(
            terminal_states_
        )
        if not len(trans_indices):
            raise RuntimeError("Cannot proceed - Markov chain is irreducible.")

        # create Q (restriction transient-transient), S (restriction transient-recurrent)
        q = t[trans_indices, :][:, trans_indices]
        s = t[trans_indices, :][:, rec_indices]

        # take individual solutions and piece them together to get absorption probabilities towards the classes
        macro_ix_helper = np.cumsum(
            [0] + [len(indices) for indices in lookup_dict.values()]
        )
        s = np.concatenate(
            [s[:, np.arange(a, b)].sum(axis=1) for a, b in _pairwise(macro_ix_helper)],
            axis=1,
        )

        # check for irreducibility
        if check_irreducibility:
            if self.is_irreducible is None:
                self._is_irreducible = _irreducible(self.transition_matrix)
            else:
                if not self.is_irreducible:
                    logg.warning("Transition matrix is not irreducible")
                else:
                    logg.debug("Transition matrix is irreducible")

        logg.debug(f"Found `{n_cells}` cells and `{s.shape[1]}` absorbing states")

        # solve the linear system of equations
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

        if time_to_absorption is not None:
            abs_time_means = _calculate_lineage_absorption_time_means(
                q,
                t[trans_indices, :][:, rec_indices],
                trans_indices,
                n=t.shape[0],
                ixs=lookup_dict,
                lineages=lin_abs_times,
                solver=solver,
                use_petsc=use_petsc,
                n_jobs=n_jobs,
                backend=backend,
                tol=tol,
                show_progress_bar=show_progress_bar,
                preconditioner=preconditioner,
            )
            abs_time_means.index = self.adata.obs_names
        else:
            abs_time_means = None

        # for recurrent states, set their self-absorption probability to one
        abs_classes = np.zeros((len(self), n_macrostates))
        rec_classes_full = {
            cl: np.where(terminal_states_ == cl)[0]
            for cl in terminal_states_.cat.categories
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

        self._set(
            A.ABS_PROBS,
            Lineage(
                abs_classes,
                names=terminal_states_.cat.categories,
                colors=colors_,
            ),
        )

        extra_msg = ""
        if abs_time_means is not None:
            self._set(A.LIN_ABS_TIMES, abs_time_means)
            extra_msg = f"       `.{P.LIN_ABS_TIMES}`\n"

        self._write_absorption_probabilities(time=start, extra_msg=extra_msg)

    @d.dedent
    def compute_lineage_priming(
        self,
        method: Literal["kl_divergence", "entropy"] = "kl_divergence",
        early_cells: Optional[Union[Mapping[str, Sequence[str]], Sequence[str]]] = None,
    ) -> pd.Series:
        """
        %(lin_pd.full_desc)s

        Parameters
        ----------
        %(lin_pd.parameters)s
            Cell ids or a mask marking early cells. If `None`, use all cells. Only used when ``method='kl_divergence'``.
            If a :class:`dict`, the key species a cluster key in :attr:`anndata.AnnData.obs` and the values
            specify cluster labels containing early cells.

        Returns
        -------
        %(lin_pd.returns)s
        """  # noqa: D400
        abs_probs: Optional[Lineage] = self._get(P.ABS_PROBS)
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

        self._set(A.PRIME_DEG, values)
        self.adata.obs[_pd(self._abs_prob_key)] = values

        logg.info(
            f"Adding `adata.obs[{_pd(self._abs_prob_key)!r}]`\n"
            f"       `.{P.PRIME_DEG}`"
        )

        return values

    @d.get_sections(
        base="lineage_drivers", sections=["Parameters", "Returns", "References"]
    )
    @d.get_full_description(base="lineage_drivers")
    @d.dedent
    @inject_docs(lin_drivers=P.LIN_DRIVERS, tm=TestMethod)
    def compute_lineage_drivers(
        self,
        lineages: Optional[Union[str, Sequence]] = None,
        method: str = TestMethod.FISCHER.s,
        cluster_key: Optional[str] = None,
        clusters: Optional[Union[str, Sequence]] = None,
        layer: str = "X",
        use_raw: bool = False,
        confidence_level: float = 0.95,
        n_perms: int = 1000,
        seed: Optional[int] = None,
        return_drivers: bool = True,
        **kwargs,
    ) -> Optional[pd.DataFrame]:
        """
        Compute driver genes per lineage.

        Correlates gene expression with lineage probabilities, for a given lineage and set of clusters.
        Often, it makes sense to restrict this to a set of clusters which are relevant for the specified lineages.

        Parameters
        ----------
        lineages
            Either a set of lineage names from :attr:`absorption_probabilities` `.names` or `None`,
            in which case all lineages are considered.
        method
            Mode to use when calculating p-values and confidence intervals. Valid options are:

                - {tm.FISCHER.s!r} - use Fischer transformation :cite:`fischer:21`.
                - {tm.PERM_TEST.s!r} - use permutation test.
        cluster_key
            Key from :attr:`adata` ``.obs`` to obtain cluster annotations. These are considered for ``clusters``.
        clusters
            Restrict the correlations to these clusters.
        layer
            Key from :attr:`adata` ``.layers``.
        use_raw
            Whether or not to use :attr:`adata` ``.raw`` to correlate gene expression.
            If using a layer other than ``.X``, this must be set to `False`.
        confidence_level
            Confidence level for the confidence interval calculation. Must be in `[0, 1]`.
        n_perms
            Number of permutations to use when ``method={tm.PERM_TEST.s!r}``.
        seed
            Random seed when ``method={tm.PERM_TEST.s!r}``.
        return_drivers
            Whether to return the drivers. This also contains the lower and upper ``confidence_level`` confidence
            interval bounds.
        %(parallel)s

        Returns
        -------
        %(correlation_test.returns)s
        Only if ``return_drivers=True``.

        Otherwise, updates :attr:`adata` ``.var`` or :attr:`adata` ``.raw.var``, depending ``use_raw`` with:

            - ``'{{direction}} {{lineage}} corr'`` - the potential lineage drivers.
            - ``'{{direction}} {{lineage}} qval'`` - the corrected p-values.

        Also updates the following fields:

            - :attr:`{lin_drivers}` - same as the returned values.
        """

        # check that lineage probs have been computed
        method = TestMethod(method)
        abs_probs = self._get(P.ABS_PROBS)
        prefix = DirPrefix.BACKWARD if self.kernel.backward else DirPrefix.FORWARD

        if abs_probs is None:
            raise RuntimeError(
                "Compute absorption probabilities first as `.compute_absorption_probabilities()`."
            )
        elif abs_probs.shape[1] == 1:
            logg.warning(
                "There is only 1 lineage present. Using the stationary distribution instead"
            )
            abs_probs = Lineage(
                self._get(P.TERM_PROBS).values,
                names=abs_probs.names,
                colors=abs_probs.colors,
            )

        # check all lin_keys exist in self.lin_names
        if isinstance(lineages, str):
            lineages = [lineages]
        if lineages is not None:
            _ = abs_probs[lineages]
        else:
            lineages = abs_probs.names

        if not len(lineages):
            raise ValueError("No lineages have been selected.")

        # use `cluster_key` and clusters to subset the data
        if clusters is not None:
            if cluster_key not in self.adata.obs.keys():
                raise KeyError(f"Key `{cluster_key!r}` not found in `adata.obs`.")
            if isinstance(clusters, str):
                clusters = [clusters]

            all_clusters = np.array(self.adata.obs[cluster_key].cat.categories)
            cluster_mask = np.array([name not in all_clusters for name in clusters])

            if any(cluster_mask):
                raise KeyError(
                    f"Clusters `{list(np.array(clusters)[cluster_mask])}` not found in "
                    f"`adata.obs[{cluster_key!r}]`."
                )
            subset_mask = np.in1d(self.adata.obs[cluster_key], clusters)
            adata_comp = self.adata[subset_mask]
            lin_probs = abs_probs[subset_mask, :]
        else:
            adata_comp = self.adata
            lin_probs = abs_probs

        # check that the layer exists, and that use raw is only used with layer X
        if layer != "X":
            if layer not in self.adata.layers:
                raise KeyError(f"Layer `{layer!r}` not found in `adata.layers`.")
            if use_raw:
                raise ValueError("For `use_raw=True`, layer must be 'X'.")
            data = adata_comp.layers[layer]
            var_names = adata_comp.var_names
        else:
            if use_raw and self.adata.raw is None:
                logg.warning("No raw attribute set. Using `.X` instead")
                use_raw = False
            data = adata_comp.raw.X if use_raw else adata_comp.X
            var_names = adata_comp.raw.var_names if use_raw else adata_comp.var_names

        start = logg.info(
            f"Computing correlations for lineages `{lineages}` restricted to clusters `{clusters}` in "
            f"layer `{layer}` with `use_raw={use_raw}`"
        )

        lin_probs = lin_probs[lineages]
        drivers = _correlation_test(
            data,
            lin_probs,
            gene_names=var_names,
            method=method,
            n_perms=n_perms,
            seed=seed,
            confidence_level=confidence_level,
            **kwargs,
        )
        self._set(A.LIN_DRIVERS, drivers)

        corrs, qvals = [f"{lin} corr" for lin in lin_probs.names], [
            f"{lin} qval" for lin in lin_probs.names
        ]
        if use_raw:
            self.adata.raw.var[[f"{prefix} {col}" for col in corrs]] = drivers[corrs]
            self.adata.raw.var[[f"{prefix} {col}" for col in qvals]] = drivers[qvals]
        else:
            self.adata.var[[f"{prefix} {col}" for col in corrs]] = drivers[corrs]
            self.adata.var[[f"{prefix} {col}" for col in qvals]] = drivers[qvals]

        field = "raw.var" if use_raw else "var"
        keys_added = [
            f"`adata.{field}['{prefix} {lin} corr']`" for lin in lin_probs.names
        ]

        logg.info(
            f"Adding `.{P.LIN_DRIVERS}`\n       "
            + "\n       ".join(keys_added)
            + "\n    Finish",
            time=start,
        )

        if return_drivers:
            return drivers

    @d.get_sections(base="plot_lineage_drivers", sections=["Parameters"])
    @d.dedent
    def plot_lineage_drivers(
        self,
        lineage: str,
        n_genes: int = 8,
        ncols: Optional[int] = None,
        use_raw: bool = False,
        title_fmt: str = "{gene} qval={qval:.4e}",
        figsize: Optional[Tuple[float, float]] = None,
        dpi: Optional[int] = None,
        save: Optional[Union[str, Path]] = None,
        **kwargs,
    ) -> None:
        """
        Plot lineage drivers discovered by :meth:`compute_lineage_drivers`.

        Parameters
        ----------
        lineage
            Lineage for which to plot the driver genes.
        n_genes
            Top most correlated genes to plot.
        ncols
            Number of columns.
        use_raw
            Whether to look in :attr:`adata` ``.raw.var`` or :attr:`adata` ``.var``.
        title_fmt
            Title format. Possible keywords include `{gene}`, `{qval}`, `{corr}` for gene name,
            q-value and correlation, respectively.
        %(plotting)s
        kwargs
            Keyword arguments for :func:`scvelo.pl.scatter`.

        Returns
        -------
        %(just_plots)s
        """

        def prepare_format(
            gene: str,
            *,
            pval: Optional[float],
            qval: Optional[float],
            corr: Optional[float],
        ) -> Dict[str, Any]:
            kwargs = {}
            if "{gene" in title_fmt:
                kwargs["gene"] = gene
            if "{pval" in title_fmt:
                kwargs["pval"] = float(pval) if pval is not None else np.nan
            if "{qval" in title_fmt:
                kwargs["qval"] = float(qval) if qval is not None else np.nan
            if "{corr" in title_fmt:
                kwargs["corr"] = float(corr) if corr is not None else np.nan

            return kwargs

        lin_drivers = self._get(P.LIN_DRIVERS)
        if lin_drivers is None:
            raise RuntimeError(
                f"Compute `.{P.LIN_DRIVERS}` first as `.compute_lineage_drivers()`."
            )

        key = f"{lineage} corr"
        if key not in lin_drivers:
            raise KeyError(
                f"Lineage `{key!r}` not found in `{list(lin_drivers.columns)}`."
            )

        if n_genes <= 0:
            raise ValueError(f"Expected `n_genes` to be positive, found `{n_genes}`.")

        kwargs.pop("save", None)
        genes = lin_drivers.sort_values(by=key, ascending=False).head(n_genes)

        ncols = 4 if ncols is None else ncols
        nrows = int(np.ceil(len(genes) / ncols))

        fig, axes = plt.subplots(
            ncols=ncols,
            nrows=nrows,
            dpi=dpi,
            figsize=(ncols * 6, nrows * 4) if figsize is None else figsize,
        )
        axes = np.ravel([axes])

        _i = 0
        for _i, (gene, ax) in enumerate(zip(genes.index, axes)):
            data = genes.loc[gene]
            scv.pl.scatter(
                self.adata,
                color=gene,
                ncols=ncols,
                use_raw=use_raw,
                ax=ax,
                show=False,
                title=title_fmt.format(
                    **prepare_format(
                        gene,
                        pval=data.get(f"{lineage} pval", None),
                        qval=data.get(f"{lineage} qval", None),
                        corr=data.get(f"{lineage} corr", None),
                    )
                ),
                **kwargs,
            )

        for j in range(_i + 1, len(axes)):
            axes[j].remove()

        if save is not None:
            save_fig(fig, save)

    @d.dedent
    def plot_lineage_drivers_correlation(
        self,
        lineage_x: str,
        lineage_y: str,
        color: Optional[str] = None,
        gene_sets: Optional[Dict[str, Iterable]] = None,
        gene_sets_colors: Optional[Iterable] = None,
        use_raw: bool = False,
        cmap: str = "RdYlBu_r",
        fontsize: int = 12,
        adjust_text: bool = False,
        legend_loc: Optional[str] = "best",
        figsize: Optional[Tuple[float, float]] = (4, 4),
        dpi: Optional[int] = None,
        save: Optional[Union[str, Path]] = None,
        show: bool = True,
        **kwargs: Any,
    ) -> Optional[Axes]:
        """
        Show scatter plot of gene-correlations between two lineages.

        Optionally, you can pass a :class:`dict` of gene names that will be annotated in the plot.

        Parameters
        ----------
        lineage_x
            Name of the lineage on the x-axis.
        lineage_y
            Name of the lineage on the y-axis.
        color
            Key in :attr:`adata` ``.var``.
        gene_sets
            Gene sets annotations of the form `{'gene_set_name': ['gene_1', 'gene_2'], ...}`.
        gene_sets_colors
            List of colors where each entry corresponds to a gene set from ``genes_sets``.
            If `None` and keys in ``gene_sets`` correspond to lineage names, use the lineage colors.
            Otherwise, use default colors.
        use_raw
            Whether to access :attr:`adata` ``.raw.var`` or :attr:`adata` ``.var``.
        cmap
            Colormap to use.
        fontsize
            Size of the text when plotting ``gene_sets``.
        adjust_text
            Whether to automatically adjust text in order to reduce overlap.
        legend_loc
            Position of the legend. If `None`, don't show the legend.
            Only used when ``gene_sets!=None``.
        %(plotting)s
        show
            If `False`, return :class:`matplotlib.pyplot.Axes`.
        kwargs
            Keyword arguments for :func:`scanpy.pl.scatter`.

        Returns
        -------
        :class:`matplotlib.pyplot.Axes`
            The axis object if ``show=False``.
        %(just_plots)s

        Notes
        -----
        This plot is based on the following
        `notebook <https://github.com/theislab/gastrulation_analysis/blob/main/6_cellrank.ipynb>`_ by Maren Büttner.
        """
        from cellrank.pl._utils import _position_legend

        if use_raw and self.adata.raw is None:
            logg.warning("No raw attribute set. Setting `use_raw=False`")
            use_raw = False
        adata = self.adata.raw if use_raw else self.adata

        # silent assumption: `.compute_lineage_drivers()` always writes to AnnData
        key1, key2 = f"to {lineage_x} corr", f"to {lineage_y} corr"
        if key1 not in adata.var or key2 not in adata.var:
            haystack = "adata.raw.var" if use_raw else "adata.var"
            raise RuntimeError(
                f"Unable to find correlations in `{haystack}[{key1!r}]` or `{haystack}[{key2!r}]`."
                f"Compute `.{P.LIN_DRIVERS}` first as "
                f"`.compute_lineage_drivers([{lineage_x!r}, {lineage_y!r}], use_raw={use_raw})`."
            )
        # produce the actual scatter plot
        ctx = {"figure.figsize": figsize, "figure.dpi": dpi}
        for key in list(ctx.keys()):
            if ctx[key] is None:
                del ctx[key]
        with rc_context(ctx):
            ax = sc.pl.scatter(
                adata.to_adata() if hasattr(adata, "to_adata") else adata,
                x=key1,
                y=key2,
                color_map=cmap,
                use_raw=False,
                show=False,
                color=color,
                **kwargs,
            )
        fig = ax.figure

        # add some lines to highlight the origin
        xmin, xmax = np.nanmin(adata.var[key1]), np.nanmax(adata.var[key1])
        ymin, ymax = np.nanmin(adata.var[key2]), np.nanmax(adata.var[key2])
        ax.hlines(0, xmin=xmin, xmax=xmax, color="grey", alpha=0.5, zorder=-1)
        ax.vlines(0, ymin=ymin, ymax=ymax, color="grey", alpha=0.5, zorder=-1)

        # annotate the passed set of genes
        if gene_sets is not None:
            if gene_sets_colors is None:
                try:
                    # fmt: off
                    sets = list(gene_sets.keys())
                    gene_sets_colors = self.adata.obsm[self._abs_prob_key][sets].colors
                    # fmt: on
                except KeyError:
                    logg.warning(
                        "Unable to determine gene sets colors from lineages. Using default colors"
                    )
                    gene_sets_colors = _create_categorical_colors(len(gene_sets))
            if len(gene_sets_colors) != len(gene_sets):
                raise ValueError(
                    f"Expected `gene_sets_colors` to be of length `{len(gene_sets)}`, "
                    f"found `{len(gene_sets_colors)}`."
                )

            path_effect = [
                patheffects.Stroke(linewidth=2, foreground="w", alpha=0.8),
                patheffects.Normal(),
            ]
            annots = []
            for (key, values), color in zip(gene_sets.items(), gene_sets_colors):
                arrowprops = (
                    {
                        "arrowstyle": ArrowStyle.CurveFilledB(
                            head_width=0.1, head_length=0.2
                        ),
                        "fc": color,
                        "ec": color,
                    }
                    if adjust_text
                    else None
                )
                values = set(values) & set(adata.var_names)
                for value in values:
                    x = adata.var.loc[value, key1]
                    y = adata.var.loc[value, key2]

                    annot = ax.annotate(
                        value,
                        xy=(x, y),
                        va="top",
                        ha="left",
                        path_effects=path_effect,
                        arrowprops=arrowprops,
                        size=fontsize,
                        c=color,
                    )
                    annots.append(annot)
                if values:
                    ax.scatter([], [], color=color, label=key)

            if adjust_text:
                try:
                    import adjustText

                    start = logg.info("Adjusting text position")
                    adjustText.adjust_text(
                        annots,
                        x=adata.var[key1].values,
                        y=adata.var[key2].values,
                        ax=ax,
                    )
                    logg.info("    Finish", time=start)
                except ImportError:
                    logg.error(
                        "Please install `adjustText` first as `pip install adjustText`"
                    )
            if len(annots) and legend_loc not in (None, "none"):
                _position_legend(ax, legend_loc=legend_loc)

        if save is not None:
            save_fig(fig, path=save)

        if not show:
            return ax

    def _detect_cc_stages(self, rc_labels: Series, p_thresh: float = 1e-15) -> None:
        """
        Detect cell-cycle driven start or endpoints.

        Parameters
        ----------
        rc_labels
            Approximate recurrent classes.
        p_thresh
            P-value threshold for the rank-sum test for the group to be considered cell-cycle driven.

        Returns
        -------
        None
            Nothing, but warns if a group is cell-cycle driven.
        """

        # initialize the groups (start or end clusters) and scores
        groups = rc_labels.cat.categories
        scores = []
        if self._G2M_score is not None:
            scores.append(self._G2M_score)
        if self._S_score is not None:
            scores.append(self._S_score)

        # loop over groups and scores
        for group in groups:
            mask = rc_labels == group
            for score in scores:
                a, b = score[mask], score[~mask]
                statistic, pvalue = ranksums(a, b)
                if statistic > 0 and pvalue < p_thresh:
                    logg.warning(f"Group `{group!r}` appears to be cell-cycle driven")
                    break

    def _set_categorical_labels(
        self,
        attr_key: str,
        color_key: str,
        pretty_attr_key: str,
        categories: Union[Series, Dict[Any, Any]],
        add_to_existing_error_msg: Optional[str] = None,
        cluster_key: Optional[str] = None,
        en_cutoff: Optional[float] = None,
        p_thresh: Optional[float] = None,
        add_to_existing: bool = False,
    ) -> None:
        if isinstance(categories, dict):
            if len(categories) == 1 and is_categorical_dtype(
                self.adata.obs.get(next(iter(categories.keys())), None)
            ):
                key = next(iter(categories.keys()))
                if isinstance(categories[key], str) or not isinstance(
                    categories[key], Iterable
                ):
                    vals = (categories[key],)
                else:
                    vals = categories[key]

                clusters = self.adata.obs[key]
                categories = {
                    cat: self.adata[clusters == cat].obs_names for cat in vals
                }

            categories = _convert_to_categorical_series(
                categories, list(self.adata.obs_names)
            )
        if not is_categorical_dtype(categories):
            raise TypeError(
                f"Object must be `categorical`, found `{infer_dtype(categories)}`."
            )

        if add_to_existing:
            if getattr(self, attr_key) is None:
                raise RuntimeError(add_to_existing_error_msg)
            categories = _merge_categorical_series(
                getattr(self, attr_key), categories, inplace=False
            )

        if cluster_key is not None:
            logg.debug(f"Creating colors based on `{cluster_key}`")

            # check that we can load the reference series from adata
            if cluster_key not in self.adata.obs:
                raise KeyError(
                    f"Cluster key `{cluster_key!r}` not found in `adata.obs`."
                )
            series_query, series_reference = categories, self.adata.obs[cluster_key]

            # load the reference colors if they exist
            if _colors(cluster_key) in self.adata.uns.keys():
                colors_reference = _convert_to_hex_colors(
                    self.adata.uns[_colors(cluster_key)]
                )
            else:
                colors_reference = _create_categorical_colors(
                    len(series_reference.cat.categories)
                )

            approx_rcs_names, colors = _map_names_and_colors(
                series_reference=series_reference,
                series_query=series_query,
                colors_reference=colors_reference,
                en_cutoff=en_cutoff,
            )
            setattr(self, color_key, colors)
            # if approx_rcs_names is categorical, the info is take from .cat.categories
            categories.cat.categories = approx_rcs_names
        else:
            setattr(
                self,
                color_key,
                _create_categorical_colors(len(categories.cat.categories)),
            )

        if p_thresh is not None:
            self._detect_cc_stages(categories, p_thresh=p_thresh)

        # write to class and adata
        if getattr(self, attr_key) is not None:
            logg.debug(f"Overwriting `.{pretty_attr_key}`")

        setattr(self, attr_key, categories)

    def _write_terminal_states(self, time=None) -> None:
        self.adata.obs[self._term_key] = self._get(P.TERM)
        self.adata.obs[_probs(self._term_key)] = self._get(P.TERM_PROBS)

        self.adata.uns[_colors(self._term_key)] = self._get(A.TERM_COLORS)
        self.adata.uns[_lin_names(self._term_key)] = np.array(
            self._get(P.TERM).cat.categories
        )

        logg.info(
            f"Adding `adata.obs[{_probs(self._term_key)!r}]`\n"
            f"       `adata.obs[{self._term_key!r}]`\n"
            f"       `.{P.TERM_PROBS}`\n"
            f"       `.{P.TERM}`\n"
            "    Finish",
            time=time,
        )

    @abstractmethod
    def _fit_terminal_states(self, *args: Any, **kwargs: Any) -> None:
        """
        High level API helper method called inside :meth:`fit` that should compute the terminal states.

        This method would usually call :meth:`compute_terminal_states` after all the functions that required beforehand.

        Parameters
        ----------
        args
            Positional arguments.
        kwargs
            Keyword arguments.

        Returns
        -------
        None
            Nothing, just sets the terminal states.

        See also
        --------
        See :meth:`cellrank.tl.estimators.GPCCA._fit_terminal_states` for an example implementation.
        """

    @d.dedent
    @inject_docs(fs=P.TERM, fsp=P.TERM_PROBS, ap=P.ABS_PROBS, pd=P.PRIME_DEG)
    def fit(
        self,
        keys: Optional[Sequence] = None,
        compute_absorption_probabilities: bool = True,
        **kwargs,
    ) -> None:
        """
        Run the pipeline.

        Parameters
        ----------
        keys
            States for which to compute absorption probabilities.
        compute_absorption_probabilities
            Whether to compute absorption probabilities or just %(initial_or_terminal)s states.
        kwargs
            Keyword arguments.

        Returns
        -------
        None
            Nothing, just makes available the following fields:

                - :attr:`{fsp}`
                - :attr:`{fs}`
                - :attr:`{ap}`
                - :attr:`{pd}`
        """
        self._fit_terminal_states(**kwargs)
        if compute_absorption_probabilities:
            self.compute_absorption_probabilities(keys=keys)

    def _write_absorption_probabilities(
        self, time: datetime, extra_msg: str = ""
    ) -> None:
        self.adata.obsm[self._abs_prob_key] = self._get(P.ABS_PROBS)

        abs_prob = self._get(P.ABS_PROBS)

        self.adata.uns[_lin_names(self._abs_prob_key)] = abs_prob.names
        self.adata.uns[_colors(self._abs_prob_key)] = abs_prob.colors

        logg.info(
            f"Adding `adata.obsm[{self._abs_prob_key!r}]`\n"
            f"{extra_msg}"
            f"       `.{P.ABS_PROBS}`\n"
            "    Finish",
            time=time,
        )

    def _set(self, n: Union[str, PrettyEnum], v: Any) -> None:
        setattr(self, n.s if isinstance(n, PrettyEnum) else n, v)

    def _get(self, n: Union[str, PrettyEnum]) -> Any:
        return getattr(self, n.s if isinstance(n, PrettyEnum) else n)

    def _set_or_debug(
        self, needle: str, haystack, attr: Optional[Union[str, PrettyEnum]] = None
    ) -> Optional[Any]:
        if isinstance(attr, PrettyEnum):
            attr = attr.s
        if needle in haystack:
            if attr is None:
                return haystack[needle]
            setattr(self, attr, haystack[needle])
        elif attr is not None:
            logg.debug(f"Unable to set attribute `.{attr}`, skipping")

    def copy(self) -> "BaseEstimator":
        """Return a copy of self, including the underlying :attr:`adata` object."""
        k = deepcopy(self.kernel)  # ensure we copy the adata object
        res = type(self)(k, read_from_adata=False)
        for k, v in self.__dict__.items():
            if isinstance(v, dict):
                res.__dict__[k] = deepcopy(v)
            elif k != "_kernel":
                res.__dict__[k] = copy(v)

        return res

    def __copy__(self) -> "BaseEstimator":
        return self.copy()

    @d.dedent
    def write(self, fname: Union[str, Path], ext: Optional[str] = "pickle") -> None:
        """
        %(pickleable.full_desc)s

        Parameters
        ----------
        %(pickleable.parameters)s

        Returns
        -------
        %(pickleable.returns)s
        """  # noqa

        fname = str(fname)
        if ext is not None:
            if not ext.startswith("."):
                ext = "." + ext
            if not fname.endswith(ext):
                fname += ext

        logg.debug(f"Writing to `{fname}`")

        with open(fname, "wb") as fout:
            if version_info[:2] > (3, 6):
                pickle.dump(self, fout)
            else:
                # we need to use PrecomputedKernel because Python3.6 can't pickle Enums
                # and they are present in VelocityKernel
                logg.warning("Saving kernel as `cellrank.tl.kernels.PrecomputedKernel`")
                orig_kernel = self.kernel
                self._kernel = PrecomputedKernel(self.kernel)
                try:
                    pickle.dump(self, fout)
                except Exception as e:  # noqa: B902
                    raise e
                finally:
                    self._kernel = orig_kernel
