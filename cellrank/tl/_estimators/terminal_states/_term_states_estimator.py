from typing import Any, Dict, Tuple, Union, Mapping, Optional, Sequence

from abc import ABC, abstractmethod
from datetime import datetime

from anndata import AnnData
from cellrank import logging as logg
from cellrank.ul._docs import d
from cellrank.tl._utils import _merge_categorical_series, _convert_to_categorical_series
from cellrank.tl._colors import (
    _map_names_and_colors,
    _convert_to_hex_colors,
    _create_categorical_colors,
)
from cellrank.tl._estimators import BaseEstimator
from cellrank.tl._estimators.mixins import CCDetectorMixin
from cellrank.tl.kernels._base_kernel import KernelExpression
from cellrank.tl._estimators.mixins._constants import Key

import numpy as np
import pandas as pd
from scipy.sparse import spmatrix
from pandas.api.types import infer_dtype, is_categorical_dtype


class TermStatesEstimator(CCDetectorMixin, BaseEstimator, ABC):
    def __init__(
        self,
        obj: Union[AnnData, np.ndarray, spmatrix, KernelExpression],
        obsp_key: Optional[str] = None,
        **kwargs: Any,
    ):
        super().__init__(obj=obj, obsp_key=obsp_key, **kwargs)
        self._term_states: Optional[pd.Series] = None
        self._term_states_probs: Optional[pd.Series] = None
        self._term_states_colors: Optional[np.ndarray] = None

    def to_adata(self) -> None:
        super().to_adata()

        key = Key.obs.term_states(self.backward)
        # TODO: set_or_debug
        if self.terminal_states is not None:
            self.adata.obs[key] = self.terminal_states
            self.adata.uns[Key.uns.colors(key)] = self._term_states_colors
        if self.terminal_states_probabilities is not None:
            self.adata.obs[Key.obs.probs(key)] = self.terminal_states_probabilities

    @d.dedent
    def set_terminal_states(
        self,
        labels: Union[pd.Series, Dict[str, Sequence[Any]]],
        cluster_key: Optional[str] = None,
        en_cutoff: Optional[float] = None,
        p_thresh: Optional[float] = None,
        add_to_existing: bool = False,
        **kwargs: Any,
    ) -> None:
        """
        Manually define terminal states.

        Parameters
        ----------
        labels
            Defines the terminal states. Valid options are:

                - categorical :class:`pandas.Series` where each category corresponds to a terminal state.
                  `NaN` entries denote cells that do not belong to any terminal state, i.e. these are either initial or
                  transient cells.
                - :class:`dict` where keys are terminal states and values are lists of cell barcodes corresponding to
                  annotations in :attr:`adata.AnnData.obs_names`.
                  If only 1 key is provided, values should correspond to terminal state clusters if a categorical
                  :class:`pandas.Series` can be found in :attr:`anndata.AnnData.obs`.

        cluster_key
            Key from :attr:`anndata.AnnData.obs` where categorical cluster labels are stored.
            These are used to associate names and colors with each terminal state.
            Each terminal state will be given the name and color corresponding to the cluster it mostly overlaps with.
        %(en_cutoff_p_thresh)s
        add_to_existing
            Whether the new terminal states should be added to pre-existing ones. Cells already assigned to a terminal
            state will be re-assigned to the new terminal state if there's a conflict between old and new annotations.
            This throws an error if no previous annotations corresponding to terminal states have been found.

        Returns
        -------
        Nothing, just updates the following fields:

            - :attr:`terminal_states` - TODO.
            - :attr:`terminal_states_probabilities` - TODO.
        """
        if add_to_existing:
            existing = self.terminal_states
            if existing is None:
                raise RuntimeError(
                    "Compute terminal states first as `.compute_terminal_states()` or "
                    "set them manually as `.set_terminal_states()`."
                )
        else:
            existing = None

        states, colors = self._set_categorical_labels(
            categories=labels,
            cluster_key=cluster_key,
            en_cutoff=en_cutoff,
            p_thresh=p_thresh,
            existing=existing,
        )
        self._write_terminal_states(states, colors, time=kwargs.get("time", None))

    def rename_terminal_states(self, new_names: Mapping[str, str]) -> None:
        """
        Rename the :attr:`terminal_states`.

        Parameters
        ----------
        new_names
            Mapping where keys corresponds to the old names and the values to the new names.
            The new names must be unique.

        Returns
        -------
        Nothing, just updates the names of :attr:`terminal_states`.
        """

        term_states = self.terminal_states
        if term_states is None:
            raise RuntimeError(
                "Compute terminal states first as `.compute_terminal_states()` or "
                "set them manually as `.set_terminal_states()`."
            )

        # fmt: off
        if not isinstance(new_names, Mapping):
            raise TypeError(f"Expected new names to be a `Mapping`, found `{type(new_names).__name__}`.")
        if not len(new_names):
            return

        new_names = {k: str(v) for k, v in new_names.items()}
        mask = np.isin(list(new_names.keys()), term_states.cat.categories)
        if not np.all(mask):
            invalid = list(np.array(list(new_names.keys()))[~mask])
            raise ValueError(f"Invalid old terminal states names: `{invalid}`.")

        names_after_renaming = [new_names.get(n, n) for n in term_states.cat.categories]
        if len(set(names_after_renaming)) != len(term_states.cat.categories):
            raise ValueError(f"After renaming, the names will not be unique: `{names_after_renaming}`.")

        self._term_states = term_states.cat.rename_categories(new_names)
        memberships = getattr(self, "_macrostates", None)
        if memberships is not None:  # GPCCA
            memberships.names = [new_names.get(n, n) for n in memberships.names]
        # fmt: on

        self.to_adata()

    def _set_categorical_labels(
        self,
        categories: Union[pd.Series, Dict[str, Any]],
        cluster_key: Optional[str] = None,
        en_cutoff: Optional[float] = None,
        p_thresh: Optional[float] = None,
        existing: Optional[pd.Series] = None,
    ) -> Tuple[pd.Series, np.ndarray]:
        # fmt: off
        if isinstance(categories, dict):
            key = next(iter(categories.keys()))
            if len(categories) == 1 and is_categorical_dtype(self.adata.obs.get(key, None)):
                vals = categories[key]
                if isinstance(vals, str) or not isinstance(vals, Sequence):
                    vals = (categories[key],)

                clusters = self.adata.obs[key]
                categories = {cat: self.adata[clusters == cat].obs_names for cat in vals}

            categories = _convert_to_categorical_series(categories, list(self.adata.obs_names))
        if not is_categorical_dtype(categories):
            raise TypeError(f"Object must be `categorical`, found `{infer_dtype(categories).__name__}`.")

        if existing is not None:
            categories = _merge_categorical_series(old=existing, new=categories)

        if cluster_key is not None:
            # check that we can load the reference series from adata
            if cluster_key not in self.adata.obs:
                raise KeyError(f"Unable to find clusters in `adata.obs[{cluster_key!r}]`.")
            series_query, series_reference = categories, self.adata.obs[cluster_key]

            # load the reference colors if they exist
            if Key.uns.colors(cluster_key) in self.adata.uns:
                colors_reference = _convert_to_hex_colors(self.adata.uns[Key.uns.colors(cluster_key)])
            else:
                colors_reference = _create_categorical_colors(len(series_reference.cat.categories))

            names, colors = _map_names_and_colors(
                series_reference=series_reference,
                series_query=series_query,
                colors_reference=colors_reference,
                en_cutoff=en_cutoff,
            )
            categories.cat.categories = names
        else:
            colors = _create_categorical_colors(len(categories.cat.categories))

        if p_thresh is not None:
            self._detect_cc_stages(categories, p_thresh=p_thresh)

        return categories, colors
        # fmt: on

    def _write_terminal_states(
        self,
        states: Optional[pd.Series],
        colors: Optional[np.ndarray],
        probs: Optional[pd.Series] = None,
        *,
        time: Optional[datetime],
        log: bool = True,
    ) -> None:
        key = Key.obs.term_states(self.backward)
        self._set("_term_states", self.adata.obs, key=key, value=states)
        self._set(
            "_term_states_probs",
            self.adata.obs,
            key=Key.obs.probs(key),
            value=probs,
        )
        self._set(
            "_term_states_colors", self.adata.uns, key=Key.uns.colors(key), value=colors
        )

        if log:
            logg.info(
                f"Adding `adata.obs[{key!r}]`\n"
                f"       `adata.obs[{Key.obs.probs(key)!r}]`\n"
                f"       `.terminal_states`\n"
                f"       `.terminal_states_probabilities`\n"
                "    Finish",
                time=time,
            )

    def fit(self, *args: Any, **kwargs: Any) -> None:
        # TODO: implement me
        return NotImplemented

    @abstractmethod
    def compute_terminal_states(self, *args: Any, **kwargs: Any) -> None:
        pass

    @property
    def terminal_states(self) -> Optional[pd.Series]:
        """TODO."""
        return self._term_states

    @property
    def terminal_states_probabilities(self) -> Optional[pd.Series]:
        """TODO."""
        return self._term_states_probs