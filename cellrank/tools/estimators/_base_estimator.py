# -*- coding: utf-8 -*-
"""Abstract base class for all kernel-holding estimators."""

import pickle
from abc import ABC, abstractmethod
from copy import copy, deepcopy
from typing import Any, Dict, Union, TypeVar, Optional, Sequence
from pathlib import Path

import numpy as np
import pandas as pd
from pandas import Series
from scipy.stats import ranksums
from scipy.sparse import spmatrix
from pandas.api.types import infer_dtype, is_categorical_dtype

from matplotlib.colors import is_color_like

from cellrank import logging as logg
from cellrank.tools import Lineage
from cellrank.utils._docs import d, inject_docs
from cellrank.tools._utils import (
    _pairwise,
    _vec_mat_corr,
    _min_max_scale,
    _process_series,
    _get_cat_and_null_indices,
    _merge_categorical_series,
    _convert_to_categorical_series,
    _calculate_absorption_time_moments,
)
from cellrank.tools._colors import (
    _map_names_and_colors,
    _convert_to_hex_colors,
    _create_categorical_colors,
)
from cellrank.tools._constants import (
    DirPrefix,
    AbsProbKey,
    PrettyEnum,
    FinalStatesKey,
    _dp,
    _probs,
    _colors,
    _lin_names,
)
from cellrank.tools._linear_solver import _solve_lin_system
from cellrank.tools.estimators._property import Partitioner, LineageEstimatorMixin
from cellrank.tools.kernels._base_kernel import KernelExpression
from cellrank.tools.estimators._constants import A, P

AnnData = TypeVar("AnnData")


@d.get_sectionsf("base_estimator", sections=["Parameters"])
class BaseEstimator(LineageEstimatorMixin, Partitioner, ABC):
    """
    Base class for all estimators.

    Parameters
    ----------
    obj
        Either a :class:`cellrank.tl.Kernel` object, an :class:`anndata.AnnData` object which
        stores the transition matrix in `.obsp` attribute or :mod:`numpy` or :mod:`scipy` array.
    inplace
        Whether to modify :paramref:`adata` object inplace or make a copy.
    read_from_adata
        Whether to read available attributes in :paramref:`adata`, if present.
    obsp_key
        Key in :paramref:`obj` `.obsp` when :paramref:`obj` is an :class:`anndata.AnnData` object.
    g2m_key
        Key from :paramref:`adata` `.obs`. Can be used to detect cell-cycle driven start- or endpoints.
    s_key
        Key from :paramref:`adata` `.obs`. Can be used to detect cell-cycle driven start- or endpoints.
    write_to_adata
        Whether to write the transition matrix to :paramref:`adata` `.obsp`.
    key_added
        Key in :paramref:`adata` where to store the final transition matrix.
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
        key_added: Optional[str] = None,
    ):
        from anndata import AnnData

        super().__init__(
            obj, obsp_key=obsp_key, key_added=key_added, write_to_adata=write_to_adata
        )

        if isinstance(obj, (KernelExpression, AnnData)) and not inplace:
            self.kernel._adata = self.adata.copy()

        if self.kernel.backward:
            self._fs_key: str = str(FinalStatesKey.BACKWARD)
            self._abs_prob_key: str = str(AbsProbKey.BACKWARD)
        else:
            self._fs_key: str = str(FinalStatesKey.FORWARD)
            self._abs_prob_key: str = str(AbsProbKey.FORWARD)

        self._key_added = key_added
        self._g2m_key = g2m_key
        self._s_key = s_key

        self._G2M_score = None
        self._S_score = None

        if read_from_adata:
            self._read_from_adata()

    def _read_from_adata(self) -> None:
        self._set_or_debug(f"eig_{self._direction}", self.adata.uns, "_eig")

        self._set_or_debug(self._g2m_key, self.adata.obs, "_G2M_score")
        self._set_or_debug(self._s_key, self.adata.obs, "_S_score")

        self._set_or_debug(self._fs_key, self.adata.obs, A.FIN.s)
        self._set_or_debug(_probs(self._fs_key), self.adata.obs, A.FIN_PROBS.s)
        self._set_or_debug(_colors(self._fs_key), self.adata.uns, A.FIN_COLORS.s)

        self._set_or_debug(self._abs_prob_key, self.adata.obsm, A.ABS_RPOBS.s)
        self._set_or_debug(_dp(self._abs_prob_key), self.adata.obs, A.DIFF_POT.s)
        names = self._set_or_debug(_lin_names(self._abs_prob_key), self.adata.uns)
        colors = self._set_or_debug(_colors(self._abs_prob_key), self.adata.uns)

        abs_probs = self._get(P.ABS_PROBS)
        if abs_probs is not None:
            if len(names) != abs_probs.shape[1]:
                logg.debug(
                    f"Expected lineage names to be of length `{abs_probs.shape[1]}`, found `{len(names)}`. "
                    f"Creating new names"
                )
                names = [f"Lineage {i}" for i in range(abs_probs.shape[1])]
            if len(colors) != abs_probs.shape[1] or not all(
                map(lambda c: isinstance(c, str) and is_color_like(c), colors)
            ):
                logg.debug(
                    f"Expected lineage colors to be of length `{abs_probs.shape[1]}`, found `{len(names)}`. "
                    f"Creating new colors"
                )
                colors = _create_categorical_colors(abs_probs.shape[1])
            self._set(A.ABS_RPOBS, Lineage(abs_probs, names=names, colors=colors))

            self.adata.obsm[self._abs_prob_key] = self._get(P.ABS_PROBS)
            self.adata.uns[_lin_names(self._abs_prob_key)] = names
            self.adata.uns[_colors(self._abs_prob_key)] = colors

    @inject_docs(fs=P.FIN.s, fsp=P.FIN_PROBS.s)
    def set_final_states(
        self,
        labels: Union[Series, Dict[Any, Any]],
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
            If a key to cluster labels is given, :paramref:`{fs}` will ge associated
            with these for naming and colors.
        en_cutoff
            If :paramref:`cluster_key` is given, this parameter determines when an approximate recurrent class will
            be labelled as *'Unknown'*, based on the entropy of the distribution of cells over transcriptomic clusters.
        p_thresh
            If cell cycle scores were provided, a *Wilcoxon rank-sum test* is conducted to identify cell-cycle driven
            start- or endpoints.
            If the test returns a positive statistic and a p-value smaller than :paramref:`p_thresh`,
            a warning will be issued.
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
            attr_key=A.FIN.v,
            color_key=A.FIN_COLORS.v,
            pretty_attr_key=P.FIN.v,
            add_to_existing_error_msg="Compute final states first as `.compute_final_states()`.",
            categories=labels,
            cluster_key=cluster_key,
            en_cutoff=en_cutoff,
            p_thresh=p_thresh,
            add_to_existing=add_to_existing,
        )
        self._write_final_states(time=kwargs.get("time", None))

    @inject_docs(abs_prob=P.ABS_PROBS, diff_pot=P.DIFF_POT)
    def compute_absorption_probabilities(
        self,
        keys: Optional[Sequence[str]] = None,
        check_irred: bool = False,
        solver: Optional[str] = None,
        use_petsc: Optional[bool] = None,
        absorption_time_moments: Optional[str] = "first",
        n_jobs: Optional[int] = None,
        backend: str = "loky",
        tol: float = 1e-5,
    ) -> None:
        """
        Compute absorption probabilities of a Markov chain.

        For each cell, this computes the probability of it reaching any of the approximate recurrent classes.
        This also computes the entropy over absorption probabilities, which is a measure of cell plasticity, see
        [Setty19]_.

        Parameters
        ----------
        keys
            Comma separated sequence of keys defining the recurrent classes.
        check_irred
            Check whether the transition matrix is irreducible.
        solver
            Solver to use for the linear problem. Options are `['direct', 'gmres', 'lgmres', 'bicgstab', 'gcrotmk']`
            when :paramref:`use_petsc` `=False` or one of :class:`petsc4py.PETSc.KPS.Type` otherwise.

            Information on the :mod:`scipy` iterative solvers can be found in :func:`scipy.sparse.linalg` or for
            :mod:`petsc4py` solver found `here <https://www.mcs.anl.gov/petsc/documentation/linearsolvertable.html>`_.

            If is `None`, solver is chosen automatically, depending on the problem.
        use_petsc
            Whether to use solvers from :mod:`petsc4py` or :mod:`scipy`. Recommended for large problems.
            If `None`, it is determined automatically. If no installation is found, defaults
            to :func:`scipy.sparse.linalg.gmres`.
        absorption_time_moments
            Whether to compute mean absorption time and its variance. Valid options are `None`, `'first'`, `'second'`.
        n_jobs
            Number of parallel jobs to use when using an iterative solver.
            When :paramref:`use_petsc` `=True` or for quickly-solvable problems,
            we recommend higher number (>=4) of jobs in order to fully saturate the cores.
        backend
            Which backend to use for multiprocessing. See :class:`joblib.Parallel` for valid options.
        tol
            Convergence tolerance for the iterative solver. The default is fine for most cases, only consider
            decreasing this for severely ill-conditioned matrices.

        Returns
        -------
        None
            Nothing, but updates the following fields:

                - :paramref:`{abs_prob}`
                - :paramref:`{diff_pot}`
        """
        if self._get(P.FIN) is None:
            raise RuntimeError(
                "Compute final states first as `.compute_final_states()` or set them manually as "
                "`.set_final_states()`."
            )
        if absorption_time_moments not in (None, "first", "second"):
            raise ValueError(
                f"Expcted `absorption_moments` to be `None`, `'first'` or `'second'`, "
                f"found `{absorption_time_moments!r}`."
            )

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
        metastable_states_, colors_ = _process_series(
            series=self._get(P.FIN), keys=keys, colors=self._get(A.FIN_COLORS)
        )

        # define the dimensions of this problem
        n_cells = t.shape[0]
        n_macrostates = len(metastable_states_.cat.categories)

        #  create empty lineage object
        if self._get(P.ABS_PROBS) is not None:
            logg.debug(f"Overwriting `.{P.ABS_PROBS}`")

        self._set(
            A.ABS_RPOBS,
            Lineage(
                np.empty((1, len(colors_))),
                names=metastable_states_.cat.categories,
                colors=colors_,
            ),
        )

        # warn in case only one state is left
        keys = list(metastable_states_.cat.categories)
        if len(keys) == 1:
            logg.warning(
                "There is only one recurrent class, all cells will have probability 1 of going there"
            )

        # get indices corresponding to recurrent and transient states
        rec_indices, trans_indices, lookup_dict = _get_cat_and_null_indices(
            metastable_states_
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

        mat_x, abs_time_mean, abs_time_var = None, None, None

        if absorption_time_moments == "second":
            mat_x, abs_time_mean, abs_time_var = _calculate_absorption_time_moments(
                q,
                s,
                rec_indices,
                trans_indices,
                solver=solver,
                use_petsc=use_petsc,
                n_jobs=n_jobs,
                backend=backend,
                tol=tol,
                use_eye=False,
            )
        elif absorption_time_moments == "first":
            mean_time = _min_max_scale(
                _solve_lin_system(
                    q,
                    np.ones((q.shape[0], 1), dtype=np.float32),
                    solver=solver,
                    use_petsc=use_petsc,
                    n_jobs=1,
                    backend=backend,
                    tol=tol,
                    use_eye=True,
                )
            )
            abs_time_mean = np.ones((self.adata.n_obs,), dtype=np.float32)
            abs_time_mean[trans_indices] = 1 - mean_time.squeeze()
        if mat_x is None:
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
            )
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
            cl: np.where(metastable_states_ == cl)[0]
            for cl in metastable_states_.cat.categories
        }
        for col, cl_indices in enumerate(rec_classes_full.values()):
            abs_classes[trans_indices, col] = _abs_classes[:, col]
            abs_classes[cl_indices, col] = 1

        self._set(
            A.ABS_RPOBS,
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

        if abs_time_mean is not None:
            self._set(
                A.MEAN_ABS_TIME, pd.Series(abs_time_mean, index=self.adata.obs.index)
            )
        if abs_time_var is not None:
            self._set(
                A.VAR_ABS_TIME, pd.Series(abs_time_var, index=self.adata.obs.index)
            )

        self._write_absorption_probabilities(time=start)

    @d.get_sectionsf("lineage_drivers", sections=["Parameters", "Returns"])
    @d.get_full_descriptionf("lineage_drivers")
    def compute_lineage_drivers(
        self,
        lineages: Optional[Union[Sequence, str]] = None,
        cluster_key: Optional[str] = None,
        clusters: Optional[Union[Sequence, str]] = None,
        layer: str = "X",
        use_raw: bool = True,
        inplace: bool = True,
    ) -> Optional[pd.DataFrame]:
        """
        Compute driver genes per lineage.

        Correlates gene expression with lineage probabilities, for a given lineage and set of clusters.
        Often, it makes sense to restrict this to a set of clusters which are relevant
        for the lineage under consideration.

        Parameters
        ----------
        lineages
            Either a set of lineage names from :paramref:`absorption_probabilities` `.names` or `None`,
            in which case all lineages are considered.
        cluster_key
            Key from :paramref:`adata` `.obs` to obtain cluster annotations.
            These are considered for :paramref:`clusters`.
        clusters
            Restrict the correlations to these clusters.
        layer
            Key from :paramref:`adata` `.layers`.
        use_raw
            Whether or not to use :paramref:`adata` `.raw` to correlate gene expression.
            If using a layer other than `.X`, this must be set to `False`.
        inplace
            Whether to write to :paramref:`adata` or return a :class:`pandas.DataFrame` object.

        Returns
        -------
        None
            If :paramref:`inplace` `=False`, writes to :paramref:`adata` `.var` or :paramref:`adata` `.raw.var`,
            depending on the value of :paramref:`use_raw`.
            For each lineage specified, a key is added to `.var` and correlations are saved as
            `{direction} {lineage_names} corr`.
        :class:`pandas.DataFrame`
            If :paramref:`inplace` `=True`, a :class:`pandas.DataFrame` with the columns same as mentioned above.
        """

        # check that lineage probs have been computed
        abs_probs = self._get(P.ABS_PROBS)
        prefix = DirPrefix.BACKWARD if self.kernel.backward else DirPrefix.FORWARD

        if abs_probs is None:
            raise RuntimeError(
                "Compute absorption probabilities first as `.compute_absorption_probabilities()`."
            )

        # check all lin_keys exist in self.lin_names
        if isinstance(lineages, str):
            lineages = [lineages]
        if lineages is not None:
            _ = abs_probs[lineages]
        else:
            lineages = abs_probs.names

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

        # loop over lineages
        lin_corrs = {}
        for lineage in lineages:
            y = lin_probs[:, lineage].X.squeeze()
            correlations = _vec_mat_corr(data, y)

            if inplace:
                if use_raw:
                    self.adata.raw.var[f"{prefix} {lineage} corr"] = correlations
                else:
                    self.adata.var[f"{prefix} {lineage} corr"] = correlations
            else:
                lin_corrs[lineage] = correlations

        if not inplace:
            return pd.DataFrame(lin_corrs, index=var_names)

        field = "raw.var" if use_raw else "var"
        logg.info(
            f"Adding gene correlations to `.adata.{field}`\n    Finish", time=start
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

    def _write_final_states(self, time=None) -> None:
        self.adata.obs[self._fs_key] = self._get(P.FIN)
        self.adata.obs[_probs(self._fs_key)] = self._get(P.FIN_PROBS)

        self.adata.uns[_colors(self._fs_key)] = self._get(A.FIN_COLORS)

        logg.info(
            f"Adding `adata.obs[{_probs(self._fs_key)!r}]`\n"
            f"       `adata.obs[{self._fs_key!r}]`\n"
            f"       `.{P.FIN_PROBS}`\n"
            f"       `.{P.FIN}`",
            time=time,
        )

    @abstractmethod
    def _fit_final_states(self, *args, **kwargs):
        pass

    @inject_docs(fs=P.FIN, fsp=P.FIN_PROBS, ap=P.ABS_PROBS, dp=P.DIFF_POT)
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
        args
            Positional arguments.
        keys
            Names of final states for which to compute absorption probabilities.
        compute_absorption_probabilities
            Whether to compute absorption probabilities or just final states.
        kwargs
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
        self._fit_final_states(**kwargs)
        if compute_absorption_probabilities:
            self.compute_absorption_probabilities(keys=keys)

    def _write_absorption_probabilities(self, time: float) -> None:
        self.adata.obsm[self._abs_prob_key] = self._get(P.ABS_PROBS)

        abs_prob = self._get(P.ABS_PROBS)
        self.adata.obs[_dp(self._abs_prob_key)] = self._get(P.DIFF_POT)

        self.adata.uns[_lin_names(self._abs_prob_key)] = abs_prob.names
        self.adata.uns[_colors(self._abs_prob_key)] = abs_prob.colors

        logg.info(
            f"Adding `adata.obsm[{self._abs_prob_key!r}]`\n"
            f"       `adata.obs[{_dp(self._abs_prob_key)!r}]`\n"
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
        self, needle: str, haystack, attr: Optional[str] = None
    ) -> Optional[Any]:
        if needle in haystack:
            if attr is None:
                return haystack[needle]
            setattr(self, attr, haystack[needle])
        elif attr is not None:
            logg.debug(f"Unable to set attribute `.{attr}`. Skipping")

    def copy(self) -> "BaseEstimator":
        """Return a copy of self, including the underlying :paramref:`adata` object."""

        k = deepcopy(self.kernel)  # ensure we copy the adata object
        res = type(self)(k, read_from_adata=False)
        for k, v in self.__dict__.items():
            if k != "_kernel":
                res.__dict__[k] = copy(v)

        return res

    def __copy__(self) -> "BaseEstimator":
        return self.copy()

    def write(self, fname: Union[str, Path]) -> None:
        """
        Serialize self to a file.

        Parameters
        ----------
        fname
            Filename where to save the object.

        Returns
        -------
        None
            Nothing, just pickles itself to a file.
        """

        fname = str(fname)
        if not fname.endswith(".pickle"):
            fname += ".pickle"

        with open(fname, "wb") as fout:
            pickle.dump(self, fout)

    @staticmethod
    def read(fname: Union[str, Path]) -> "BaseEstimator":
        """
        Deserialize self from a file.

        Parameters
        ----------
        fname
            Filename from which to read the object.

        Returns
        -------
        :class:`cellrank.tl.estimators._base_estimator.BaseEstimator`
            An estimator.
        """

        with open(fname, "rb") as fin:
            return pickle.load(fin)
