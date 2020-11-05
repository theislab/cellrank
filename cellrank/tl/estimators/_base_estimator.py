# -*- coding: utf-8 -*-
"""Abstract base class for all kernel-holding estimators."""
import pickle
from abc import ABC, abstractmethod
from sys import version_info
from copy import copy, deepcopy
from math import ceil
from typing import Any, Dict, Union, Mapping, TypeVar, Optional, Sequence
from pathlib import Path
from datetime import datetime

import scvelo as scv

import numpy as np
import pandas as pd
from pandas import Series
from scipy.stats import ranksums
from scipy.sparse import spmatrix
from pandas.api.types import infer_dtype, is_categorical_dtype

from matplotlib.colors import is_color_like

from cellrank import logging as logg
from cellrank.ul._docs import d, inject_docs
from cellrank.tl._utils import (
    TestMethod,
    _pairwise,
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
    _dp,
    _probs,
    _colors,
    _lin_names,
)
from cellrank.tl._linear_solver import _solve_lin_system
from cellrank.tl.estimators._property import Partitioner, LineageEstimatorMixin
from cellrank.tl.kernels._base_kernel import KernelExpression
from cellrank.tl.estimators._constants import A, P

AnnData = TypeVar("AnnData")
_COMP_TERM_STATES_MSG = (
    "Compute terminal states first as `.compute_terminal_states()` or"
    "set them manually as `.set_terminal_states()`."
)


@d.get_sectionsf("base_estimator", sections=["Parameters"])
@d.dedent
class BaseEstimator(LineageEstimatorMixin, Partitioner, ABC):
    """
    Base class for all estimators.

    Parameters
    ----------
    obj
        Either a :class:`cellrank.tl.Kernel` object, an :class:`anndata.AnnData` object which
        stores the transition matrix in ``.obsp`` attribute or :mod:`numpy` or :mod:`scipy` array.
    inplace
        Whether to modify :paramref:`adata` object inplace or make a copy.
    read_from_adata
        Whether to read available attributes in :paramref:`adata`, if present.
    obsp_key
        Key in ``obj.obsp`` when ``obj`` is an :class:`anndata.AnnData` object.
    g2m_key
        Key in :paramref:`adata` ``.obs``. Can be used to detect cell-cycle driven start- or endpoints.
    s_key
        Key in :paramref:`adata` ``.obs``. Can be used to detect cell-cycle driven start- or endpoints.
    write_to_adata
        Whether to write the transition matrix to :paramref:`adata` ``.obsp`` and the parameters to
        :paramref:`adata` ``.uns``.
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

    def _read_from_adata(self) -> None:
        self._set_or_debug(f"eig_{self._direction}", self.adata.uns, "_eig")

        self._set_or_debug(self._g2m_key, self.adata.obs, "_G2M_score")
        self._set_or_debug(self._s_key, self.adata.obs, "_S_score")

        self._set_or_debug(self._term_key, self.adata.obs, A.TERM.s)
        self._set_or_debug(_probs(self._term_key), self.adata.obs, A.TERM_PROBS)
        self._set_or_debug(_colors(self._term_key), self.adata.uns, A.TERM_COLORS)

        self._reconstruct_lineage(A.ABS_PROBS, self._abs_prob_key)
        self._set_or_debug(_dp(self._abs_prob_key), self.adata.obs, A.DIFF_POT)

    def _reconstruct_lineage(self, attr: PrettyEnum, obsm_key: str):

        self._set_or_debug(obsm_key, self.adata.obsm, attr)
        names = self._set_or_debug(_lin_names(self._term_key), self.adata.uns)
        colors = self._set_or_debug(_colors(self._term_key), self.adata.uns)

        # choosing this instead of property because GPCCA doesn't have property for FIN_ABS_PROBS
        probs = self._get(attr)

        if probs is not None:
            if len(names) != probs.shape[1]:
                logg.debug(
                    f"Expected lineage names to be of length `{probs.shape[1]}`, found `{len(names)}`. "
                    f"Creating new names"
                )
                names = [f"Lineage {i}" for i in range(probs.shape[1])]
            if len(colors) != probs.shape[1] or not all(
                map(lambda c: isinstance(c, str) and is_color_like(c), colors)
            ):
                logg.debug(
                    f"Expected lineage colors to be of length `{probs.shape[1]}`, found `{len(names)}`. "
                    f"Creating new colors"
                )
                colors = _create_categorical_colors(probs.shape[1])
            self._set(attr, Lineage(probs, names=names, colors=colors))

            self.adata.obsm[obsm_key] = self._get(attr)
            self.adata.uns[_lin_names(self._term_key)] = names
            self.adata.uns[_colors(self._term_key)] = colors

    @d.dedent
    @inject_docs(fs=P.TERM.s, fsp=P.TERM_PROBS.s)
    def set_terminal_states(
        self,
        labels: Union[Series, Dict[str, Any]],
        cluster_key: Optional[str] = None,
        en_cutoff: Optional[float] = None,
        p_thresh: Optional[float] = None,
        add_to_existing: bool = False,
        **kwargs,
    ) -> None:
        """
        Set the approximate recurrent classes, if they are known a priori.

        Parameters
        ----------
        labels
            Either a categorical :class:`pandas.Series` with index as cell names, where `NaN` marks marks a cell
            belonging to a transient state or a :class:`dict`, where each key is the name of the recurrent class and
            values are list of cell names.
        cluster_key
            If a key to cluster labels is given, :paramref:`{fs}` will ge associated with these for naming and colors.
        %(en_cutoff_p_thresh)s
        add_to_existing
            Whether to add thses categories to existing ones. Cells already belonging to recurrent classes will be
            updated if there's an overlap.
            Throws an error if previous approximate recurrent classes have not been calculated.

        Returns
        -------
        None
            Nothing, but updates the following fields:

                - :paramref:`{fsp}`
                - :paramref:`{fs}`
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
        Rename the names of :paramref:`{ts}`.

        Parameters
        ----------
        new_names
            Mapping where keys are the old names and the values are the new names. New names must be unique.
        update_adata
            Whether to update underlying :paramref:`adata` object as well or not.

        Returns
        -------
        None
            Nothing, just updates the names of :paramref:`{ts}`.
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
        diff_pot=P.DIFF_POT,
        lat=P.LIN_ABS_TIMES,
    )
    def compute_absorption_probabilities(
        self,
        keys: Optional[Sequence[str]] = None,
        check_irred: bool = False,
        solver: Optional[str] = None,
        use_petsc: Optional[bool] = None,
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
        tol: float = 1e-5,
        preconditioner: Optional[str] = None,
    ) -> None:
        """
        Compute absorption probabilities of a Markov chain.

        For each cell, this computes the probability of it reaching any of the approximate recurrent classes defined
        by :paramref:`{fs}`. This also computes the entropy over absorption probabilities, which is a measure of cell
        plasticity, see [Setty19]_.

        Parameters
        ----------
        keys
            Keys defining the recurrent classes.
        check_irred
            Check whether the transition matrix is irreducible.
        solver
            Solver to use for the linear problem. Options are `'direct', 'gmres', 'lgmres', 'bicgstab' or 'gcrotmk'`
            when ``use_petsc=False`` or one of :class:`petsc4py.PETSc.KPS.Type` otherwise.

            Information on the :mod:`scipy` iterative solvers can be found in :func:`scipy.sparse.linalg` or for
            :mod:`petsc4py` solver `here <https://www.mcs.anl.gov/petsc/documentation/linearsolvertable.html>`__.

            If is `None`, the solver is chosen automatically, depending on the problem size.
        use_petsc
            Whether to use solvers from :mod:`petsc4py` or :mod:`scipy`. Recommended for large problems.
            If `None`, it is determined automatically. If no installation is found, defaults to
            :func:`scipy.sparse.linalg.gmres`.
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

                - :paramref:`{abs_prob}` - probabilities of being absorbed into the terminal states.
                - :paramref:`{diff_pot}` - differentiation potential of cells.
                - :paramref:`{lat}` - mean times until absorption to subset absorbing states and optionally
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

        #  create empty lineage object
        if self._get(P.ABS_PROBS) is not None:
            logg.debug(f"Overwriting `.{P.ABS_PROBS}`")

        self._set(
            A.ABS_PROBS,
            Lineage(
                np.empty((1, len(colors_))),
                names=terminal_states_.cat.categories,
                colors=colors_,
            ),
        )

        # get indices corresponding to recurrent and transient states
        rec_indices, trans_indices, lookup_dict = _get_cat_and_null_indices(
            terminal_states_
        )
        if not len(trans_indices):
            raise RuntimeError("Cannot proceed - Markov chain is irreducible.")

        # create Q (restriction transient-transient), S (restriction transient-recurrent)
        q = t[trans_indices, :][:, trans_indices]
        s = t[trans_indices, :][:, rec_indices]

        if check_irred:
            if self.is_irreducible is None:
                self.compute_partition()
            if not self.is_irreducible:
                logg.warning("The transition matrix is not irreducible")

        # determine whether it makes sense you use a iterative solver
        if solver is None:
            solver = (
                "gmres"
                if self.issparse
                and (n_cells >= 3e5 or (n_cells >= 1e4 and s.shape[1] <= 100))
                else "direct"
            )
        if use_petsc is None:
            use_petsc = n_cells >= 3e5

        logg.debug(f"Found `{n_cells}` cells and `{s.shape[1]}` absorbing states")

        # solve the linear system of equations
        mat_x = _solve_lin_system(
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

        # take individual solutions and piece them together to get absorption probabilities towards the classes
        macro_ix_helper = np.cumsum(
            [0] + [len(indices) for indices in lookup_dict.values()]
        )
        _abs_classes = np.concatenate(
            [
                mat_x[:, np.arange(a, b)].sum(1)[:, None]
                for a, b in _pairwise(macro_ix_helper)
            ],
            axis=1,
        )

        # for recurrent states, set their self-absorption probability to one
        abs_classes = np.zeros((len(self), n_macrostates))
        rec_classes_full = {
            cl: np.where(terminal_states_ == cl)[0]
            for cl in terminal_states_.cat.categories
        }
        for col, cl_indices in enumerate(rec_classes_full.values()):
            abs_classes[trans_indices, col] = _abs_classes[:, col]
            abs_classes[cl_indices, col] = 1

        self._set(
            A.ABS_PROBS,
            Lineage(
                abs_classes,
                names=self._get(P.ABS_PROBS).names,
                colors=self._get(P.ABS_PROBS).colors,
            ),
        )
        self._set(
            A.DIFF_POT,
            pd.Series(
                self._get(P.ABS_PROBS).entropy(axis=1).X.squeeze(axis=1),
                index=self.adata.obs.index,
            ),
        )

        extra_msg = ""
        if abs_time_means is not None:
            self._set(A.LIN_ABS_TIMES, abs_time_means)
            extra_msg = f"       `.{P.LIN_ABS_TIMES}`\n"

        self._write_absorption_probabilities(time=start, extra_msg=extra_msg)

    @d.get_sectionsf(
        "lineage_drivers", sections=["Parameters", "Returns", "References"]
    )
    @d.get_full_descriptionf("lineage_drivers")
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
            Either a set of lineage names from :paramref:`absorption_probabilities` `.names` or `None`,
            in which case all lineages are considered.
        method
            Mode to use when calculating p-values and confidence intervals. Can be one of:

                - {tm.FISCHER.s!r} - use Fischer transformation [Fischer21]_.
                - {tm.PERM_TEST.s!r} - use permutation test.
        cluster_key
            Key from :paramref:`adata` ``.obs`` to obtain cluster annotations. These are considered for ``clusters``.
        clusters
            Restrict the correlations to these clusters.
        layer
            Key from :paramref:`adata` ``.layers``.
        use_raw
            Whether or not to use :paramref:`adata` ``.raw`` to correlate gene expression.
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

        None
            Updates :paramref:`adata` ``.var`` or :paramref:`adata` ``.raw.var``, depending ``use_raw`` with:

                - ``'{{direction}} {{lineage}} corr'`` - the potential lineage drivers.
                - ``'{{direction}} {{lineage}} qval'`` - the corrected p-values.

            Updates the following fields:

                - :paramref:`{lin_drivers}` - same as the returned values.

        References
        ----------
        .. [Fischer21] Fisher, R. A. (1921),
            *On the “probable error” of a coefficient of correlation deduced from a small sample.*,
            `Metron 1 3–32 <http://hdl.handle.net/2440/15169>`__.
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

        drivers = _correlation_test(
            data,
            lin_probs[lineages],
            gene_names=var_names,
            method=method,
            n_perms=n_perms,
            seed=seed,
            confidence_level=confidence_level,
            **kwargs,
        )
        self._set(A.LIN_DRIVERS, drivers)

        corrs, qvals = [f"{lin} corr" for lin in lineages], [
            f"{lin} qval" for lin in lineages
        ]
        if use_raw:
            self.adata.raw.var[[f"{prefix} {col}" for col in corrs]] = drivers[corrs]
            self.adata.raw.var[[f"{prefix} {col}" for col in qvals]] = drivers[qvals]
        else:
            self.adata.var[[f"{prefix} {col}" for col in corrs]] = drivers[corrs]
            self.adata.var[[f"{prefix} {col}" for col in qvals]] = drivers[qvals]

        field = "raw.var" if use_raw else "var"
        keys_added = [f"`adata.{field}['{prefix} {lin} corr']`" for lin in lineages]

        logg.info(
            f"Adding `.{P.LIN_DRIVERS}`\n       "
            + "\n       ".join(keys_added)
            + "\n    Finish",
            time=start,
        )

        if return_drivers:
            return drivers

    @d.get_sectionsf("plot_lineage_drivers", sections=["Parameters"])
    @d.dedent
    def plot_lineage_drivers(
        self,
        lineage: str,
        n_genes: int = 8,
        use_raw: bool = False,
        **kwargs,
    ) -> None:
        """
        Plot lineage drivers discovered by :meth:`compute_lineage_drivers`.

        Parameters
        ----------
        lineage
            Lineage for which to plot the driver genes.
        n_genes
            Number of genes to plot.
        use_raw
            Whether to look in :paramref:`adata` ``.raw.var`` or :paramref:`adata` ``.var``.
        **kwargs
            Keyword arguments for :func:`scvelo.pl.scatter`.

        Returns
        -------
        %(just_plots)s
        """

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

        cmap = kwargs.pop("cmap", "viridis")
        ncols = kwargs.pop("ncols", None)

        geness = lin_drivers.sort_values(by=key, ascending=False).head(n_genes)
        kwargs.setdefault(
            "title",
            [
                f"{g} q-value={q}"
                for g, q in dict(geness[f"{lineage} qval"].map("{:.4e}".format)).items()
            ],
        )

        if ncols is None and len(geness) >= 8:
            ncols = 4

        # TODO: scvelo can handle only < 20 plots, see https://github.com/theislab/scvelo/issues/252
        geness = filter(len, np.array_split(geness.index, int(ceil(len(geness) / 8))))

        for genes in geness:
            scv.pl.scatter(
                self.adata,
                color=genes,
                cmap=cmap,
                ncols=ncols,
                use_raw=use_raw,
                **kwargs,
            )

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

        extra_msg = ""
        if getattr(self, A.TERM_ABS_PROBS.s, None) is not None and hasattr(
            self, "_term_abs_prob_key"
        ):
            # checking for None because terminal states can be set using `set_terminal_states`
            # without the probabilities in GPCCA
            self.adata.obsm[self._term_abs_prob_key] = self._get(A.TERM_ABS_PROBS)
            extra_msg = f"       `adata.obsm[{self._term_abs_prob_key!r}]`\n"

        logg.info(
            f"Adding `adata.obs[{_probs(self._term_key)!r}]`\n"
            f"       `adata.obs[{self._term_key!r}]`\n"
            f"{extra_msg}"
            f"       `.{P.TERM_PROBS}`\n"
            f"       `.{P.TERM}`",
            time=time,
        )

    @abstractmethod
    def _fit_terminal_states(self, *args, **kwargs):
        pass

    @d.dedent
    @inject_docs(fs=P.TERM, fsp=P.TERM_PROBS, ap=P.ABS_PROBS, dp=P.DIFF_POT)
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
        **kwargs
            Keyword arguments.

        Returns
        -------
        None
            Nothing, just makes available the following fields:

                - :paramref:`{fsp}`
                - :paramref:`{fs}`
                - :paramref:`{ap}`
                - :paramref:`{dp}`
        """
        self._fit_terminal_states(**kwargs)
        if compute_absorption_probabilities:
            self.compute_absorption_probabilities(keys=keys)

    def _write_absorption_probabilities(
        self, time: datetime, extra_msg: str = ""
    ) -> None:
        self.adata.obsm[self._abs_prob_key] = self._get(P.ABS_PROBS)

        abs_prob = self._get(P.ABS_PROBS)
        self.adata.obs[_dp(self._abs_prob_key)] = self._get(P.DIFF_POT)

        self.adata.uns[_lin_names(self._abs_prob_key)] = abs_prob.names
        self.adata.uns[_colors(self._abs_prob_key)] = abs_prob.colors

        logg.info(
            f"Adding `adata.obsm[{self._abs_prob_key!r}]`\n"
            f"       `adata.obs[{_dp(self._abs_prob_key)!r}]`\n"
            f"{extra_msg}"
            f"       `.{P.ABS_PROBS}`\n"
            f"       `.{P.DIFF_POT}`\n"
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
        """Return a copy of self, including the underlying :paramref:`adata` object."""
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
                except Exception as e:
                    raise e
                finally:
                    self._kernel = orig_kernel
