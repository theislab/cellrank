from typing import Any, Tuple, Union, Mapping, Optional, Sequence
from typing_extensions import Literal

from datetime import datetime

from anndata import AnnData
from cellrank import logging as logg
from cellrank.tl import Lineage
from cellrank.ul._docs import d
from cellrank.tl._utils import (
    _eigengap,
    _fuzzy_to_discrete,
    _series_from_one_hot_matrix,
)
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

        self._term_states_memberships: Optional[Lineage] = None

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
            Number of macrostates. If `None`, use the *eigengap* heuristic.
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
                - :attr:`eigendecomposition` - TODO.
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
            # TODO: parse the error, optionally reraise
            # this is the following case - we have 4 Schur vectors, user requests 5 states, but it splits the conj. ev.
            # in the try block, Schur decomposition with 5 vectors is computed, but it fails (no way of knowing)
            # so in this case, we increase it by 1
            n_states += 1
            logg.warning(f"{e}\nUsing `n_states={n_states}`")
            self._gpcca = self._gpcca.optimize(m=n_states)

        self._set_macrostates(
            memberships=self._gpcca.memberships,
            n_cells=n_cells,
            cluster_key=cluster_key,
            p_thresh=p_thresh,
            en_cutoff=en_cutoff,
            time=start,
        )

    @d.dedent
    def compute_terminal_states(
        self,
        method: Literal[
            "stability", "top_n", "eigengap", "eigengap_coarse"
        ] = "stability",
        n_cells: int = 30,
        alpha: Optional[float] = 1,
        stability_threshold: float = 0.96,
        n_states: Optional[int] = None,
    ):
        """
        Automatically select terminal states from macrostates.

        Parameters
        ----------
        method
            How to select the terminal states. Valid option are:

                - `'eigengap'` - select the number of states based on the *eigengap* of :attr:`transition_matrix`.
                - `'eigengap_coarse'` - select the number of states based on the `*eigengap* of the diagonal
                  of :attr:`coarse_T`.
                - `'top_n'` - select top ``n_states`` based on the probability of the diagonal of :attr:`coarse_T`.
                - `'stability'` - select states which have a stability >= ``stability_threshold``.
                  The stability is given by the diagonal elements of :attr:`coarse_T`.
        %(n_cells)s
        alpha
            Weight given to the deviation of an eigenvalue from one.
            Only used when ``method = 'eigengap'`` or ``method = 'eigengap_coarse'``.
        stability_threshold
            Threshold used when ``method = 'stability'``.
        n_states
            Number of states used when ``method = 'top_n'``.

        Returns
        -------
        Nothing, just updates the following fields:

            - :attr:`terminal_states` - TODO.
            - :attr:`terminal_states_memberships` - TODO.
            - :attr:`terminal_states_probabilities` - TODO.
        """

        # TODO: ModeEnum
        if len(self._macrostates.cat.categories) == 1:
            logg.warning(
                "Found only one macrostate. Making it the single terminal state"
            )
            self.set_terminal_states_from_macrostates(None, n_cells=n_cells)
            return

        eig = self.eigendecomposition
        coarse_T = self.coarse_T

        # fmt: off
        if method == "eigengap":
            if eig is None:
                raise RuntimeError("Compute eigendecomposition first as `.compute_eigendecomposition()`.")
            n_states = _eigengap(eig["D"], alpha=alpha) + 1
        elif method == "eigengap_coarse":
            if coarse_T is None:
                raise RuntimeError("Compute macrostates first as `.compute_macrostates()`.")
            n_states = _eigengap(np.sort(np.diag(coarse_T)[::-1]), alpha=alpha)
        elif method == "top_n":
            if n_states is None:
                raise ValueError("Expected `n_states != None` for `method='top_n'`.")
            elif n_states <= 0:
                raise ValueError(f"Expected `n_states` to be positive, found `{n_states}`.")
        elif method == "stability":
            if stability_threshold is None:
                raise ValueError("Expected`stability_threshold != None` for `method='stability'`.")
            stability = pd.Series(np.diag(coarse_T), index=coarse_T.columns)
            names = stability[stability.values >= stability_threshold].index
            self.set_terminal_states_from_macrostates(names, n_cells=n_cells)
            return
        else:
            raise NotImplementedError(f"Method `{method}` is not yet implemented.")
        # fmt: on

        names = coarse_T.columns[np.argsort(np.diag(coarse_T))][-n_states:]
        self.set_terminal_states_from_macrostates(names, n_cells=n_cells)

    @d.dedent
    def set_terminal_states_from_macrostates(
        self,
        names: Optional[Union[str, Sequence[str], Mapping[str, str]]] = None,
        n_cells: int = 30,
    ):
        """
        Manually select terminal states from macrostates.

        Parameters
        ----------
        names
            Names of the macrostates to be marked as terminal. Multiple states can be combined using `','`,
            such as ``["Alpha, Beta", "Epsilon"]``.  If a :class:`dict`, keys correspond to the names
            of the macrostates and the values to the new names.  If `None`, select all macrostates.
        %(n_cells)s

        Returns
        -------
        # TODO: docrep
        Nothing, just updates the following fields:

            - :attr:`terminal_states` - TODO.
            - :attr:`terminal_states_probabilities` - TODO.
        """
        if n_cells <= 0:
            raise ValueError(f"Expected `n_cells` to be positive, found `{n_cells}`.")

        memberships = self.macrostates_memberships
        if memberships is None:
            raise RuntimeError("Compute macrostates first as `.compute_macrostates()`.")

        rename = True
        if names is None:
            names = memberships.names
            rename = False
        if isinstance(names, str):
            names = [names]
            rename = False
        if not isinstance(names, dict):
            names = {n: n for n in names}
            rename = False
        if not len(names):
            raise ValueError("No macrostates have been selected.")

        # we do this also here because if `rename_terminal_states` fails
        # invalid states would've been written to this object and nothing to adata
        names = {str(k): str(v) for k, v in names.items()}
        names_after_renaming = {names.get(n, n) for n in memberships.names}
        if len(names_after_renaming) != memberships.shape[1]:
            raise ValueError(
                f"After renaming, terminal state names will no longer be unique: `{names_after_renaming}`."
            )

        # this also checks that the names are correct before renaming
        is_singleton = memberships.shape[1] == 1
        memberships = memberships[list(names.keys())].copy()

        states = self._create_states(memberships, n_cells=n_cells, check_row_sums=False)
        if is_singleton:
            colors = self._macrostates_colors.copy()
            probs = memberships.X.squeeze() / memberships.X.max()
        else:
            colors = memberships[list(states.cat.categories)].colors
            probs = (memberships.X / memberships.X.max(0)).max(1)
        probs = pd.Series(probs, index=self.adata.obs_names)

        self._write_terminal_states(states, colors, probs, memberships)
        if rename:
            self.rename_terminal_states(names)

    def fit(self, *args: Any, **kwargs: Any) -> None:
        # TOOO: call super + optionally abs prob?
        return NotImplemented

    def _n_states(self, n_states: Optional[Union[int, Sequence[int]]]) -> int:
        if n_states is None:
            if self.eigendecomposition is None:
                raise RuntimeError(
                    "Compute eigendecomposition first as `.compute_eigendecomposition()` or "
                    "supply `n_states != None`."
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
            raise ValueError(f"Minimum value must be > `1`, found `{minn}`.")
        elif minn == 2:
            logg.warning(
                "In most cases, 2 clusters will always be optimal. "
                "If you really expect 2 clusters, use `n_states=2`. Setting the minimum to `3`"
            )
            minn = 3
        # fmt: on
        maxx = max(minn + 1, maxx)

        logg.info(f"Calculating minChi criterion in interval `[{minn}, {maxx}]`")

        return int(np.arange(minn, maxx + 1)[np.argmax(self._gpcca.minChi(minn, maxx))])

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
                f"split the conjugate eigenvalues. Using `n_states={n_states + 1}`"
            )
            n_states += 1  # cannot force recomputation of the Schur decomposition
            assert n_states not in self._invalid_n_states, "Sanity check failed."

        return n_states

    def _compute_one_macrostate(
        self,
        n_cells: Optional[int],
        cluster_key: Optional[str],
        en_cutoff: Optional[float],
        p_thresh: float,
    ) -> None:
        start = logg.info("For 1 macrostate, stationary distribution is computed")

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
        self._write_terminal_states(None, None, None, None, log=False)

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
            for attr in ["_schur_vectors", "_schur_matrix", "_coarse_tmat", "_coarse_init_dist", "_coarse_stat_dist"]:
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

    def _write_terminal_states(
        self,
        states: Optional[pd.Series],
        colors: Optional[np.ndarray],
        probs: Optional[pd.Series] = None,
        memberships: Optional[Lineage] = None,
        *,
        time: Optional[datetime] = None,
        log: bool = True,
    ) -> None:
        super()._write_terminal_states(states, colors, probs, time=time, log=log)

        self._set("_term_states_memberships", value=memberships)

    plot_macrostates = register_plotter(
        discrete="macrostates", continuous="macrostates_memberships"
    )
    plot_terminal_states = register_plotter(
        discrete="terminal_states", continuous="terminal_states_memberships"
    )

    @property
    def macrostates(self) -> Optional[pd.Series]:
        """TODO."""
        return self._macrostates

    @property
    def macrostates_memberships(self) -> Optional[Lineage]:
        """TODO."""
        return self._macrostates_memberships

    @property
    def terminal_states_memberships(self) -> Optional[Lineage]:
        """TODO."""
        return self._term_states_memberships

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
