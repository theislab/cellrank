from typing import Any, Dict, Tuple, Union, Mapping, Optional, Sequence

from abc import ABC
from types import MappingProxyType

from anndata import AnnData
from cellrank._key import Key
from cellrank.ul._docs import d
from cellrank.tl._utils import _merge_categorical_series, _convert_to_categorical_series
from cellrank.tl._colors import (
    _map_names_and_colors,
    _convert_to_hex_colors,
    _create_categorical_colors,
)
from cellrank.tl.estimators import BaseEstimator
from cellrank.tl.estimators._utils import SafeGetter
from cellrank.tl.kernels._base_kernel import KernelExpression
from cellrank.tl.estimators.mixins._utils import logger, shadow, register_plotter

import numpy as np
import pandas as pd
from scipy.sparse import spmatrix
from pandas.api.types import infer_dtype, is_categorical_dtype

from matplotlib.colors import to_hex


@d.dedent
class TermStatesEstimator(BaseEstimator, ABC):
    """
    Base class for all estimators predicting terminal states.

    Parameters
    ----------
    %(base_estimator.parameters)s
    """

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

    @property
    @d.get_summary(base="tse_term_states")
    def terminal_states(self) -> Optional[pd.Series]:
        """Categorical annotation of terminal states.

        By default, all cells in transient cells will be labelled as `NaN`.
        """
        return self._term_states

    @property
    @d.get_summary(base="tse_term_states_probs")
    def terminal_states_probabilities(self) -> Optional[pd.Series]:
        """Aggregated probability of cells to be in terminal states."""  # noqa: D401
        return self._term_states_probs

    @d.dedent
    def set_terminal_states(
        self,
        labels: Union[pd.Series, Dict[str, Sequence[Any]]],
        cluster_key: Optional[str] = None,
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
            Key in :attr:`anndata.AnnData.obs` in order to associate names and colors with :attr:`terminal_states`.
            Each terminal state will be given the name and color corresponding to the cluster it mostly overlaps with.
        add_to_existing
            Whether the new terminal states should be added to pre-existing ones. Cells already assigned to a terminal
            state will be re-assigned to the new terminal state if there's a conflict between old and new annotations.
            This throws an error if no previous annotations corresponding to terminal states have been found.

        Returns
        -------
        Nothing, just updates the following fields:

            - :attr:`terminal_states` - %(tse_term_states.summary)s
            - :attr:`terminal_states_probabilities` - %(tse_term_states_probs.summary)s
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
            existing=existing,
        )
        self._write_terminal_states(
            states,
            colors,
            **kwargs,
        )

    @d.get_sections(base="tse_rename_term_states", sections=["Parameters", "Returns"])
    @d.get_full_description(base="tse_rename_term_states")
    @d.dedent
    def rename_terminal_states(self, new_names: Mapping[str, str]) -> None:
        """
        Rename categories in :attr:`terminal_states`.

        Parameters
        ----------
        new_names
            Mapping where keys corresponds to the old names and the values to the new names.
            The new names must be unique.

        Returns
        -------
        Nothing, just updates the names of:

            - :attr:`terminal_states` - %(tse_term_states.summary)s
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

        old_names = term_states.cat.categories
        new_names = {str(k): str(v) for k, v in new_names.items()}
        mask = np.isin(list(new_names.keys()), old_names)
        if not np.all(mask):
            invalid = sorted(np.array(list(new_names.keys()))[~mask])
            raise ValueError(f"Invalid terminal states names: `{invalid}`. Valid names are: `{sorted(old_names)}`")

        names_after_renaming = [new_names.get(n, n) for n in old_names]
        if len(set(names_after_renaming)) != len(old_names):
            raise ValueError(f"After renaming, terminal states will no longer unique: `{names_after_renaming}`.")
        # fmt: on

        # GPCCA; alt. is to subclass
        self._term_states = term_states.cat.rename_categories(new_names)
        memberships = getattr(self, "terminal_states_memberships", None)
        if memberships is not None:
            memberships.names = [new_names.get(n, n) for n in memberships.names]

        self._write_terminal_states(
            self.terminal_states,
            self._term_states_colors,
            self.terminal_states_probabilities,
            log=False,
        )

    def _set_categorical_labels(
        self,
        categories: Union[pd.Series, Dict[str, Any]],
        cluster_key: Optional[str] = None,
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
            raise TypeError(f"Expected object to be `categorical`, found `{infer_dtype(categories)}`.")

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
            )
            categories.cat.categories = names
        else:
            colors = _create_categorical_colors(len(categories.cat.categories))

        return categories, colors
        # fmt: on

    @logger
    @shadow
    def _write_terminal_states(
        self,
        states: Optional[pd.Series],
        colors: Optional[np.ndarray],
        probs: Optional[pd.Series] = None,
        params: Dict[str, Any] = MappingProxyType({}),
    ) -> str:
        # fmt: off
        key = Key.obs.term_states(self.backward)
        self._set("_term_states", self.adata.obs, key=key, value=states)
        self._set("_term_states_probs", self.adata.obs, key=Key.obs.probs(key), value=probs)
        self._set("_term_states_colors", self.adata.uns, key=Key.uns.colors(key), value=colors)
        self.params[key] = dict(params)
        # fmt: on

        return (
            f"Adding `adata.obs[{key!r}]`\n"
            f"       `adata.obs[{Key.obs.probs(key)!r}]`\n"
            f"       `.terminal_states`\n"
            f"       `.terminal_states_probabilities`\n"
            "    Finish"
        )

    def _read_from_adata(self, adata: AnnData, **kwargs: Any) -> bool:
        ok = super()._read_from_adata(adata, **kwargs)
        if not ok:
            return False

        # fmt: off
        key = Key.obs.term_states(self.backward)
        with SafeGetter(self, allowed=KeyError) as sg:
            self._get("_term_states", self.adata.obs, key=key, where="obs", dtype=pd.Series)
            self._get("_term_states_probs", self.adata.obs, key=Key.obs.probs(key), where="obs", dtype=pd.Series)
            self._get("_term_states_colors", self.adata.uns, key=Key.uns.colors(key), where="uns",
                      dtype=(list, tuple, np.ndarray))
            self._term_states_colors = np.asarray([to_hex(c) for c in self._term_states_colors])
            self.params[key] = self._read_params(key)
        # fmt: on

        return sg.ok

    @d.dedent
    def compute_terminal_states(self, *args: Any, **kwargs: Any) -> None:
        """
        Compute terminal states of the process.

        This is an alias for :meth:`predict`.

        Parameters
        ----------
        args
            Positional arguments.
        kwargs
            Keyword arguments arguments.

        Return
        ------
        Nothing, just updates the following fields:

            - :attr:`terminal_states` - %(tse_term_states.summary)s
            - :attr:`terminal_states_probabilities` - %(tse_term_states_probs.summary)s
        """
        return self.predict(*args, **kwargs)

    plot_terminal_states = register_plotter(
        discrete="terminal_states", colors="_term_states_colors"
    )
