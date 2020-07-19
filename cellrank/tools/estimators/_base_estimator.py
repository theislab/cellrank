# -*- coding: utf-8 -*-
"""Abstract base class for all kernel-holding estimators."""

from abc import ABC, abstractmethod
from typing import Any, Dict, Union, Optional, Sequence
from pathlib import Path
from functools import partial

import numpy as np
import pandas as pd
from pandas import Series
from scipy.stats import ranksums
from scipy.sparse import spmatrix
from pandas.api.types import infer_dtype, is_categorical_dtype

from anndata import AnnData

from cellrank import logging as logg
from cellrank.tools import Lineage
from cellrank.tools._utils import (
    _pairwise,
    _vec_mat_corr,
    _process_series,
    _solve_lin_system,
    _get_cat_and_null_indices,
    _merge_categorical_series,
    _convert_to_categorical_series,
)
from cellrank.tools._colors import (
    _map_names_and_colors,
    _convert_to_hex_colors,
    _create_categorical_colors,
)
from cellrank.tools._constants import (
    Direction,
    AbsProbKey,
    PrettyEnum,
    FinalStatesKey,
    _dp,
    _probs,
    _colors,
    _lin_names,
)
from cellrank.tools.kernels._kernel import KernelExpression
from cellrank.tools.estimators._property import Partitioner, LineageEstimatorMixin
from cellrank.tools.estimators._constants import A, F, P


class BaseEstimator(LineageEstimatorMixin, Partitioner, ABC):
    """
    Base class for all estimators.

    Params
    ------
    kernel
        Kernel object that stores a transition matrix.
    adata : :class:`anndata.AnnData`
        Optional annotated data object. If given, pre-computed lineages can be read in from this.
        Otherwise, read the object from the specified :paramref:`kernel`.
    inplace
        Whether to modify :paramref:`adata` object inplace or make a copy.
    read_from_adata
        Whether to read available attributes in :paramref:`adata`, if present.
    g2m_key
        Key from :paramref:`adata` `.obs`. Can be used to detect cell-cycle driven start- or endpoints.
    s_key
        Key from :paramref:`adata` `.obs`. Can be used to detect cell-cycle driven start- or endpoints.
    key_added
        Key in :paramref:`adata` where to store the final transition matrix.
    """

    def __init__(
        self,
        kernel: Union[KernelExpression, AnnData, spmatrix, np.ndarray],
        inplace: bool = True,
        read_from_adata: bool = False,
        obsp_key: Optional[str] = None,
        g2m_key: Optional[str] = "G2M_score",
        s_key: Optional[str] = "S_score",
        key_added: Optional[str] = None,
        *args,
        **kwargs,
    ):
        super().__init__(kernel, obsp_key=obsp_key, key_added=key_added)

        if isinstance(kernel, KernelExpression) and not inplace:
            self.kernel._adata = self.adata.copy()

        if self._direction == Direction.BACKWARD:
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

        self._set = partial(
            lambda o, n, v: setattr(o, n.s if isinstance(n, PrettyEnum) else n, v), self
        )
        self._get = partial(
            lambda o, n: getattr(o, n.s if isinstance(n, PrettyEnum) else n), self
        )

        if read_from_adata:
            self._read_from_adata(*args, **kwargs)

    # TODO: refactor...
    def _read_from_adata(
        self, g2m_key: Optional[str] = None, s_key: Optional[str] = None, **kwargs
    ) -> None:
        raise NotImplementedError()

        if f"eig_{self._direction}" in self._adata.uns.keys():
            self._eig = self._adata.uns[f"eig_{self._direction}"]
        else:
            logg.debug(f"`eig_{self._direction}` not found. Setting `.eig` to `None`")

        if self._rc_key in self._adata.obs.keys():
            self._meta_states = self._adata.obs[self._rc_key]
        else:
            logg.debug(
                f"`{self._rc_key}` not found in `adata.obs`. Setting `.metastable_states` to `None`"
            )

        if _colors(self._rc_key) in self._adata.uns.keys():
            self._meta_states_colors = self._adata.uns[_colors(self._rc_key)]
        else:
            logg.debug(
                f"`{_colors(self._rc_key)}` not found in `adata.uns`. "
                f"Setting `.metastable_states_colors`to `None`"
            )

        if self._lin_key in self._adata.obsm.keys():
            lineages = range(self._adata.obsm[self._lin_key].shape[1])
            colors = _create_categorical_colors(len(lineages))
            self._lin_probs = Lineage(
                self._adata.obsm[self._lin_key],
                names=[f"Lineage {i + 1}" for i in lineages],
                colors=colors,
            )
            self._adata.obsm[self._lin_key] = self._lin_probs
        else:
            logg.debug(
                f"`{self._lin_key}` not found in `adata.obsm`. Setting `.lin_probs` to `None`"
            )

        if _dp(self._lin_key) in self._adata.obs.keys():
            self._dp = self._adata.obs[_dp(self._lin_key)]
        else:
            logg.debug(
                f"`{_dp(self._lin_key)}` not found in `adata.obs`. Setting `.diff_potential` to `None`"
            )

        if g2m_key and g2m_key in self._adata.obs.keys():
            self._G2M_score = self._adata.obs[g2m_key]
        else:
            logg.debug(
                f"`{g2m_key}` not found in `adata.obs`. Setting `.G2M_score` to `None`"
            )

        if s_key and s_key in self._adata.obs.keys():
            self._S_score = self._adata.obs[s_key]
        else:
            logg.debug(
                f"`{s_key}` not found in `adata.obs`. Setting `.S_score` to `None`"
            )

        if _probs(self._rc_key) in self._adata.obs.keys():
            self._meta_states_probs = self._adata.obs[_probs(self._rc_key)]
        else:
            logg.debug(
                f"`{_probs(self._rc_key)}` not found in `adata.obs`. "
                f"Setting `.metastable_states_probs` to `None`"
            )

        if self._lin_probs is not None:
            if _lin_names(self._lin_key) in self._adata.uns.keys():
                self._lin_probs = Lineage(
                    np.array(self._lin_probs),
                    names=self._adata.uns[_lin_names(self._lin_key)],
                    colors=self._lin_probs.colors,
                )
                self._adata.obsm[self._lin_key] = self._lin_probs
            else:
                logg.debug(
                    f"`{_lin_names(self._lin_key)}` not found in `adata.uns`. "
                    f"Using default names"
                )

            if _colors(self._lin_key) in self._adata.uns.keys():
                self._lin_probs = Lineage(
                    np.array(self._lin_probs),
                    names=self._lin_probs.names,
                    colors=self._adata.uns[_colors(self._lin_key)],
                )
                self._adata.obsm[self._lin_key] = self._lin_probs
            else:
                logg.debug(
                    f"`{_colors(self._lin_key)}` not found in `adata.uns`. "
                    f"Using default colors"
                )

    def set_final_states(
        self,
        labels: Union[Series, Dict[Any, Any]],
        cluster_key: Optional[str] = None,
        en_cutoff: Optional[float] = None,
        p_thresh: Optional[float] = None,
        add_to_existing: bool = False,
        **kwargs,
    ):
        """
        Set the approximate recurrent classes, if they are known a priori.

        Params
        ------
        categories
            Either a categorical :class:`pandas.Series` with index as cell names, where `NaN` marks marks a cell
            belonging to a transient state or a :class:`dict`, where each key is the name of the recurrent class and
            values are list of cell names.
        cluster_key
            If a key to cluster labels is given, :paramref"`metastable_states` will ge associated
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

                - :paramref:`final_states`.
        """

        self._set_categorical_labels(
            attr_key=A.FIN.v,
            color_key=A.FIN_COLORS.v,
            pretty_attr_key=P.FIN.v,
            add_to_existing_error_msg=f"Compute final states first as `.{F.COMPUTE.fmt(P.FIN)}()`.",
            categories=labels,
            cluster_key=cluster_key,
            en_cutoff=en_cutoff,
            p_thresh=p_thresh,
            add_to_existing=add_to_existing,
        )
        self._write_final_states(time=kwargs.get("time", None))

    def compute_absorption_probabilities(
        self,
        keys: Optional[Sequence[str]] = None,
        check_irred: bool = False,
        solver: Optional[str] = None,
        use_petsc: Optional[bool] = None,
        preconditioner: Optional[str] = None,
        n_jobs: Optional[int] = None,
        backend: str = "loky",
        tol: float = 1e-5,
    ) -> None:
        """
        Compute absorption probabilities for a Markov chain.

        For each cell, this computes the probability of it reaching any of the approximate recurrent classes.
        This also computes the entropy over absorption probabilities, which is a measure of cell plasticity, see
        [Setty19]_.

        Params
        ------
        keys
            Comma separated sequence of keys defining the recurrent classes.
        n_cells
            TODO
        check_irred
            Check whether the transition matrix is irreducible.
        solver
            Solver to use for the linear problem. Options are `['direct', 'gmres', 'lgmres', 'bicgstab', 'gcrotmk']`
            when :paramref:`use_petsc` `=False` or one of :class:`petsc4py.PETSc.KPS.Type` otherwise.

            Information on the :mod:`scipy` iterative solvers can be found in :func:`scipy.sparse.linalg` or for
            :mod:`petsc4py` solver found `here <https://www.mcs.anl.gov/petsc/documentation/linearsolvertable.html/>`_.

            If is `None`, solver is chosen automatically, depending on the problem.
        use_petsc
            Whether to use solvers from :mod:`petsc4py` or :mod:`scipy`. Recommended for large problems.
            If `None`, it is determined automatically. If no installation is found, defaults
            to :func:`scipy.sparse.linalg.gmres`.
        preconditioner
            Preconditioner to use when :paramref:`use_petsc` `=True`.
            For available preconditioner types, see `petsc4py.PETSc.PC.Type`.
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

                - :paramref:`absorption_probabilities`
                - :paramref:`diff_potential`
        """  # TODO: SSoT - docs edition
        if self._get(P.FIN) is None:
            raise RuntimeError(
                f"Compute final states first as `.{F.COMPUTE.fmt(P.FIN)}()` or set them manually as "
                f"`.set_final_states()`."
            )
        if keys is not None:
            keys = sorted(set(keys))

        # Note: There are three relevant data structures here
        # - self.metastable_states: pd.Series which contains annotations for approx rcs. Associated colors in
        #   self.metastable_states_colors
        # - self.lin_probs: Linage object which contains the lineage probabilities with associated names and colors
        # -_metastable_states: pd.Series, temporary copy of self.approx rcs used in the context of this function.
        #   In this copy, some metastable_states may be removed or combined with others
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
        if self._get(P.ABS_RPOBS) is not None:
            logg.debug(f"Overwriting `.{P.ABS_RPOBS}`")

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
            use_petsc = solver != "direct" and n_cells >= 3e5

        logg.debug(f"Found `{n_cells}` cells and `{s.shape[1]}` absorbing states")

        # solve the linear system of equations
        mat_x = _solve_lin_system(
            q,
            s,
            solver=solver,
            use_petsc=use_petsc,
            preconditioner=preconditioner,
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
            cl: np.where(metastable_states_ == cl)
            for cl in metastable_states_.cat.categories
        }
        for col, cl_indices in enumerate(rec_classes_full.values()):
            abs_classes[trans_indices, col] = _abs_classes[:, col]
            abs_classes[cl_indices, col] = 1

        self._set(
            A.ABS_RPOBS,
            Lineage(
                abs_classes,
                names=self._get(P.ABS_RPOBS).names,
                colors=self._get(P.ABS_RPOBS).colors,
            ),
        )

        self._set(
            A.DIFF_POT,
            pd.Series(
                self._get(P.ABS_RPOBS).entropy(axis=1).X.squeeze(axis=1),
                index=self.adata.obs.index,
            ),
        )

        self._write_absorption_probabilities(time=start)

    def compute_lineage_drivers(
        self,
        lin_names: Optional[Union[Sequence, str]] = None,
        cluster_key: Optional[str] = None,
        clusters: Optional[Union[Sequence, str]] = None,
        layer: str = "X",
        use_raw: bool = True,
        inplace: bool = True,
    ):
        """
        Compute driver genes per lineage.

        Correlates gene expression with lineage probabilities, for a given lineage and set of clusters.
        Often, it makes sense to restrict this to a set of clusters which are relevant
        for the lineage under consideration.

        Params
        ------
        lin_names
            Either a set of lineage names from :paramref:`absorption_probabilities` `.names` or None,
            in which case all lineages are considered.
        cluster_key
            Key from :paramref:`adata` `.obs` to obtain cluster annotations.
            These are considered for :paramref:`clusters`. Default is `"clusters"` if a list of `clusters` is given.
        clusters
            Restrict the correlations to these clusters.
        layer
            Key from :paramref:`adata` `.layers`.
        use_raw
            Whether or not to use :paramref:`adata` `.raw` to correlate gene expression.
            If using a layer other than `.X`, this must be set to `False`.

        Returns
        -------
        :class:`pandas.DataFrame`, :class:`NoneType`
            Writes to :paramref:`adata` `.var` or :paramref:`adata` `.raw.var`,
            depending on the value of :paramref:`use_raw`.
            For each lineage specified, a key is added to `.var` and correlations are saved there.

            Returns `None` if :paramref:`inplace` `=True`, otherwise a :class:`pandas.DataFrame`.
        """

        # check that lineage probs have been computed
        abs_probs = self._get(P.ABS_RPOBS)
        prefix = "from" if self.kernel.backward else "from"

        if abs_probs is None:
            raise RuntimeError(
                f"Compute absorption probabilities first as `.{F.COMPUTE.fmt(P.ABS_RPOBS)}()`."
            )

        # check all lin_keys exist in self.lin_names
        if isinstance(lin_names, str):
            lin_names = [lin_names]
        if lin_names is not None:
            _ = abs_probs[lin_names]
        else:
            lin_names = abs_probs.names

        # use `cluster_key` and clusters to subset the data
        if clusters is not None:
            if cluster_key is None:
                cluster_key = "clusters"
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
            adata_comp = self.adata[subset_mask].copy()
            lin_probs = abs_probs[subset_mask, :]
        else:
            adata_comp = self.adata.copy()
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
            f"Computing correlations for lineages `{lin_names}` restricted to clusters `{clusters}` in "
            f"layer `{layer}` with `use_raw={use_raw}`"
        )

        # loop over lineages
        lin_corrs = {}
        for lineage in lin_names:
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

        Params
        ------
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
    ):
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

    def _write_final_states(self, time=None):
        self.adata.obs[self._fs_key] = self._get(P.FIN)
        self.adata.obs[_probs(self._fs_key)] = self._get(P.FIN_PROBS)

        self.adata.uns[_colors(self._fs_key)] = self._get(A.FIN_COLORS)

        logg.info(
            f"Adding `adata.obs[{_probs(self._fs_key)!r}]`\n"
            f"       `adata.obs[{self._fs_key!r}]`\n"
            f"       `.{P.FIN_PROBS}`\n"
            f"       `.{P.FIN}`\n",
            time=time,
        )

    def _write_absorption_probabilities(self, time):
        self.adata.obsm[self._abs_prob_key] = self._get(P.DIFF_POT)

        abs_prob = self._get(P.ABS_RPOBS)
        self.adata.obs[_dp(self._abs_prob_key)] = abs_prob

        self.adata.uns[_lin_names(self._abs_prob_key)] = abs_prob.names
        self.adata.uns[_colors(self._abs_prob_key)] = abs_prob.colors

        logg.info(
            f"Adding `adata.obs[{self._abs_prob_key!r}]`\n"
            f"       `adata.obs[{_probs(self._abs_prob_key)!r}]`\n"
            f"       `.{P.ABS_RPOBS}`\n"
            f"       `.{P.DIFF_POT}`\n"
            "    Finish",
            time=time,
        )

    @abstractmethod
    def copy(self) -> "BaseEstimator":
        """Return a copy of self."""
        pass

    def __copy__(self) -> "BaseEstimator":
        return self.copy()

    @abstractmethod
    def __getstate__(self) -> Dict[str, Any]:
        pass

    @abstractmethod
    def __setstate__(self, state: Dict[str, Any]) -> None:
        pass

    def write(ielf, fname: Union[str, Path]) -> None:
        """
        Serialize self to a file.

        Params
        ------
        fname
            Filename where to save the object.

        Returns
        -------
        None
            Nothing, just saves itself to as file.
        """

        raise NotImplementedError()

    @classmethod
    def read(cls, fname: Union[str, Path]) -> "BaseEstimator":
        """
        Deserialize self to a file.

        Params
        ------
        fname
            Filename from which to read the object.

        Returns
        -------
        An estimator.
        """

        raise NotImplementedError()
