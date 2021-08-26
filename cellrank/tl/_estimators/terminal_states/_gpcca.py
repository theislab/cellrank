from typing import Any, Tuple, Union, Optional, Sequence

from datetime import datetime

from anndata import AnnData
from cellrank import logging as logg
from cellrank.tl import Lineage
from cellrank.ul._docs import d
from cellrank.tl._utils import _fuzzy_to_discrete, _series_from_one_hot_matrix
from cellrank.tl._estimators.mixins import EigenMixin, SchurMixin, LinDriversMixin
from cellrank.tl.kernels._base_kernel import KernelExpression
from cellrank.tl._estimators.mixins._utils import register_plotter
from cellrank.tl._estimators.terminal_states._term_states_estimator import (
    TermStatesEstimator,
)

import numpy as np
import pandas as pd
from scipy.sparse import spmatrix


class GPCCA(TermStatesEstimator, LinDriversMixin, SchurMixin, EigenMixin):
    def __init__(
        self,
        obj: Union[AnnData, np.ndarray, spmatrix, KernelExpression],
        obsp_key: Optional[str] = None,
        **kwargs: Any,
    ):
        super().__init__(obj=obj, obsp_key=obsp_key, **kwargs)

        self._coarse_init_dist: Optional[pd.Series] = None
        self._coarse_stat_dist: Optional[pd.Series] = None
        self._coarse_tmat: Optional[pd.DataFrame] = None

        self._macrostates: Optional[pd.Series] = None
        self._macrostates_memberships: Optional[Lineage] = None
        self._macrostates_colors: Optional[np.ndarray] = None

        self._term_states_cont: Optional[Lineage] = None

    @d.dedent
    def compute_macrostates(
        self,
        n_states: Optional[Union[int, Sequence[int]]] = None,
        n_cells: Optional[int] = 30,
        cluster_key: Optional[str] = None,
        en_cutoff: Optional[float] = 0.7,
        p_thresh: float = 1e-15,
    ):
        """
        Compute the macrostates.

        Parameters
        ----------
        n_states
            Number of macrostates. If `None`, use the `eigengap` heuristic.
        %(n_cells)s
        cluster_key
            If a key to cluster labels is given, names and colors of the states will be associated with the clusters.
        %(en_cutoff_p_thresh)s

        Returns
        -------
        None
            Nothing, but updates the following fields:

                - :attr:`macrostates` - TODO.
                - :attr:`macrostates_memberships` - TODO.
                - :attr:`coarse_T` - TODO.
                - :attr:`coarse_initial_distribution` - TODO.
                - :attr:`coarse_stationary_distribution` - TODO.
        """

        n_states = self._n_states(n_states)
        if n_states == 1:
            self._compute_one_macrostate(
                n_cells=n_cells,
                cluster_key=cluster_key,
                p_thresh=p_thresh,
                en_cutoff=en_cutoff,
            )
            return

        if self._gpcca is None:
            self.compute_schur(n_states)
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
            # TODO: parse the error
            # this is the following case - we have 4 Schur vectors, user requests 5 states, but it splits the conj. ev.
            # in the try block, Schur decomposition with 5 vectors is computed, but it fails (no way of knowing)
            # so in this case, we increase it by 1
            n_states += 1
            logg.warning(f"{e}\nIncreasing `n_states={n_states}`")
            self._gpcca = self._gpcca.optimize(m=n_states)

        self._set_macrostates(
            memberships=self._gpcca.memberships,
            n_cells=n_cells,
            cluster_key=cluster_key,
            p_thresh=p_thresh,
            en_cutoff=en_cutoff,
            time=start,
        )

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
                f"split the conjugate eigenvalues. Increasing `n_states={n_states + 1}`"
            )
            n_states += 1  # cannot force recomputation of the Schur decomposition
            assert n_states not in self._invalid_n_states, "Sanity check failed."

        return n_states

    def compute_terminal_states(self, *args: Any, **kwargs: Any) -> None:
        pass

    plot_macrostates = register_plotter(
        discrete="macrostates", continuous="macrostates_memberships"
    )
    plot_terminal_states = register_plotter(
        discrete="terminal_states", continuous="_term_states_cont"
    )

    def fit(self, *args: Any, **kwargs: Any) -> None:
        # TOOO: call super + optionally abs prob?
        return NotImplemented

    def _n_states(self, n_states: Optional[Union[int, Sequence[int]]]) -> int:
        if n_states is None:
            if self.eigendecomposition is None:
                raise RuntimeError("TODO.")
            return self.eigendecomposition["eigengap"] + 1

        if isinstance(n_states, int):
            if n_states <= 0:
                raise ValueError("TODO")
            return n_states

        if self._gpcca is None:
            raise RuntimeError(
                "Compute Schur decomposition first as `.compute_schur()` when `use_min_chi=True`."
            )

        if not isinstance(n_states, Sequence):
            raise TypeError(
                f"Expected `n_states` to be a `Sequence`, found `{type(n_states).__name__!r}`."
            )
        if len(n_states) != 2:
            raise ValueError(
                f"Expected `n_states` to be of size `2`, found `{len(n_states)}`."
            )

        minn, maxx = sorted(n_states)
        if minn <= 1:
            raise ValueError(f"Minimum value must be > `1`, found `{minn}`.")
        elif minn == 2:
            logg.warning(
                "In most cases, 2 clusters will always be optimal. "
                "If you really expect 2 clusters, use `n_states=2` and `use_min_chi=False`. "
                "Setting the minimum to `3`"
            )
            minn = 3
        maxx = max(minn + 1, maxx)

        logg.info(f"Calculating minChi criterion in interval `[{minn}, {maxx}]`")

        return int(np.arange(minn, maxx + 1)[np.argmax(self._gpcca.minChi(minn, maxx))])

    def _compute_one_macrostate(
        self,
        n_cells: Optional[int],
        cluster_key: Optional[str],
        en_cutoff: Optional[float],
        p_thresh: float,
    ) -> None:
        start = logg.warning("For 1 macrostate, stationary distribution is computed")

        eig = self.eigendecomposition
        if (
            eig is not None
            and "stationary_dist" in eig
            and eig["params"]["which"] == "LR"
        ):
            stationary_dist = eig["stationary_dist"]
        else:
            self.compute_eigendecomposition(only_evals=False, which="LR")
            stationary_dist = self.eigendecomposition["stationary_dist"]

        self._set_macrostates(
            memberships=stationary_dist[:, None],
            n_cells=n_cells,
            cluster_key=cluster_key,
            check_row_sums=False,
            p_thresh=p_thresh,
            en_cutoff=en_cutoff,
            time=start,
        )

    @d.dedent
    def _set_macrostates(
        self,
        memberships: np.ndarray,
        n_cells: Optional[int] = 30,
        cluster_key: str = "clusters",
        en_cutoff: Optional[float] = 0.7,
        p_thresh: float = 1e-15,
        check_row_sums: bool = True,
        time: Optional[datetime] = None,
    ) -> None:
        """
        Map fuzzy clustering to pre-computed annotations to get names and colors.

        Given the fuzzy clustering, we would like to select the most likely cells from each state and use these to
        give each state a name and a color by comparing with pre-computed, categorical cluster annotations.

        Parameters
        ----------
        memberships
            Fuzzy clustering.
        %(n_cells)s
        cluster_key
            Key from :attr:`adata` ``.obs`` to get reference cluster annotations.
        en_cutoff
            Threshold to decide when we we want to warn the user about an uncertain name mapping. This happens when
            one fuzzy state overlaps with several reference clusters, and the most likely cells are distributed almost
            evenly across the reference clusters.
        p_thresh
            Only used to detect cell cycle stages. These have to be present in :attr:`adata` ``.obs`` as
            `'G2M_score'` and `'S_score'`.
        check_row_sums
            Check whether rows in `memberships` sum to `1`.

        Returns
        -------
        TODO.
        """

        if n_cells is None:
            # fmt: off
            logg.debug("Setting the macrostates using macrostate assignment")
            assignment = pd.Series(np.argmax(memberships, axis=1).astype(str), dtype="category")
            # sometimes, a category can be missing
            assignment = assignment.cat.reorder_categories([str(i) for i in range(memberships.shape[1])])
            not_enough_cells = []
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
        self._write_terminal_states(None, None, None, log=False)

        assignment, colors = self._set_categorical_labels(
            assignment, cluster_key=cluster_key, en_cutoff=en_cutoff, p_thresh=p_thresh
        )
        names = list(assignment.cat.categories)
        memberships = Lineage(memberships, names=names, colors=colors)

        self._set("_macrostates", value=assignment)
        self._set("_macrostates_colors", value=colors)
        self._set("_macrostates_memberships", value=memberships)

        # fmt: off
        if len(names) > 1:
            # not using stationary distribution
            g = self._gpcca
            self._set("_schur_vectors", value=g._p_X)
            self._set("_schur_matrix", value=g._p_R)
            tmat = pd.DataFrame(g.coarse_grained_transition_matrix, index=names, columns=names)
            self._set("_coarse_tmat", value=tmat)
            self._set("_coarse_init_dist", value=pd.Series(g.coarse_grained_input_distribution, index=names))
            if g.coarse_grained_stationary_probability is None:
                self._set("_coarse_stat_dist", value=None)
            else:
                self._set("_coarse_stat_dist", value=pd.Series(g.coarse_grained_stationary_probability, index=names))
        else:
            for attr in ["_schur_vectors", "_schur_matrix", "_coarse_tmat", "_course_init_dist", "_coarse_stat_dist"]:
                self._set(attr, value=None)
        # fmt: on
        # TODO: logg

        # _set_categorical_labels creates the names, we still need to remap the group names
        # TODO
        # orig_cats = macrostates.cat.categories
        # name_mapper = dict(zip(orig_cats, self._get(P.MACRO).cat.categories))
        # _print_insufficient_number_of_cells(
        #    [name_mapper.get(n, n) for n in not_enough_cells], n_cells
        # )

    @property
    def macrostates(self) -> Optional[pd.Series]:
        """TODO."""
        return self._macrostates

    @property
    def macrostates_memberships(self) -> Optional[Lineage]:
        """TODO."""
        return self._macrostates_memberships

    @property
    def coarse_initial_distribution(self) -> Optional[pd.Series]:
        """TODO."""
        return self._coarse_init_dist

    @property
    def coarse_stationary_distribution(self) -> Optional[pd.Series]:
        """TODO."""
        return self._coarse_stat_dist

    @property
    def coarse_T(self) -> Optional[pd.DataFrame]:
        """TODO."""
        return self._coarse_tmat
