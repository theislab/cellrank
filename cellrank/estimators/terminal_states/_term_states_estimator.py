from typing import Any, Dict, Tuple, Union, Mapping, Optional, Sequence
from typing_extensions import Literal

from abc import ABC
from types import MappingProxyType

from anndata import AnnData
from cellrank._utils._key import Key
from cellrank._utils._docs import d
from cellrank._utils._utils import (
    _merge_categorical_series,
    _convert_to_categorical_series,
)
from cellrank._utils._colors import (
    _map_names_and_colors,
    _convert_to_hex_colors,
    _create_categorical_colors,
)
from cellrank.kernels._base_kernel import KernelExpression
from cellrank.estimators.mixins._utils import (
    SafeGetter,
    StatesHolder,
    logger,
    shadow,
    register_plotter,
)
from cellrank.estimators._base_estimator import BaseEstimator

import numpy as np
import pandas as pd
from scipy.sparse import spmatrix
from pandas.api.types import infer_dtype, is_categorical_dtype

from matplotlib.colors import to_hex

__all__ = ["TermStatesEstimator"]


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
        object: Union[AnnData, np.ndarray, spmatrix, KernelExpression],
        **kwargs: Any,
    ):
        super().__init__(object=object, **kwargs)
        self._init_states = StatesHolder()
        self._term_states = StatesHolder()

    @property
    @d.get_summary(base="tse_term_states")
    def terminal_states(self) -> Optional[pd.Series]:
        """
        Categorical annotation of terminal states.

        By default, all transient cells will be labeled as `NaN`.
        """
        return self._term_states.assignment

    @property
    @d.get_summary(base="tse_term_states_probs")
    def terminal_states_probabilities(self) -> Optional[pd.Series]:
        """Aggregated probability of cells to be in terminal states."""  # noqa: D401
        return self._term_states.probs

    @property
    @d.get_summary(base="tse_term_states")
    def initial_states(self) -> Optional[pd.Series]:
        """
        Categorical annotation of initial states.

        By default, all transient cells will be labeled as `NaN`.
        """
        return self._init_states.assignment

    @property
    @d.get_summary(base="tse_term_states_probs")
    def initial_states_probabilities(self) -> Optional[pd.Series]:
        """Aggregated probability of cells to be in initial states."""  # noqa: D401
        return self._init_states.probs

    @d.dedent
    def set_states(
        self,
        labels: Union[pd.Series, Dict[str, Sequence[Any]]],
        which: Literal["initial", "terminal"] = "terminal",
        cluster_key: Optional[str] = None,
        **kwargs: Any,
    ) -> "TermStatesEstimator":
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
                  annotations in :attr:`anndata.AnnData.obs_names`.
                  If only 1 key is provided, values should correspond to terminal state clusters if a categorical
                  :class:`pandas.Series` can be found in :attr:`anndata.AnnData.obs`.
        which: TODO(michalk8)

        cluster_key
            Key in :attr:`anndata.AnnData.obs` in order to associate names and colors with :attr:`terminal_states`.
            Each terminal state will be given the name and color corresponding to the cluster it mostly overlaps with.

        Returns
        -------
        Self and updates the following fields: TODO(michalk8)

            - :attr:`terminal_states` - %(tse_term_states.summary)s
            - :attr:`terminal_states_probabilities` - %(tse_term_states_probs.summary)s
        """
        states, colors = self._set_categorical_labels(
            categories=labels,
            cluster_key=cluster_key,
            existing=None,
        )
        self._write_states(
            which,
            states=states,
            colors=colors,
            **kwargs,
        )
        return self

    @d.get_sections(base="tse_rename_term_states", sections=["Parameters", "Returns"])
    @d.get_full_description(base="tse_rename_term_states")
    @d.dedent
    def rename_states(  # TODO(michalk8): override in GPCCA to handle macrostates
        self,
        new_names: Mapping[str, str],
        which: Literal["initial", "terminal"] = "terminal",
    ) -> "TermStatesEstimator":
        """
        Rename categories in :attr:`terminal_states` or :attr:`initial_states`.

        Parameters
        ----------
        new_names
            Mapping where keys corresponds to the old names and the values to the new names.
            The new names must be unique.
        which: TODO(michalk8)

        Returns
        -------
        Self and updates one the following field: TODO(michalk8)

            - :attr:`terminal_states` - %(tse_term_states.summary)s
        """
        backward = which == "initial"
        states = self.initial_states if backward else self.terminal_states
        if states is None:
            raise RuntimeError(
                f"Compute {which} states first as `.compute_states(..., which={which!r})` or "
                f"set them manually as `.set_states(..., which={which!r})`."
            )

        # fmt: off
        if not isinstance(new_names, Mapping):
            raise TypeError(f"Expected new names to be a `Mapping`, found `{type(new_names).__name__}`.")
        if not len(new_names):
            return self

        old_names = states.cat.categories
        new_names = {str(k): str(v) for k, v in new_names.items()}
        mask = np.isin(list(new_names.keys()), old_names)
        if not np.all(mask):
            invalid = sorted(np.array(list(new_names.keys()))[~mask])
            raise ValueError(f"Invalid {which} states names: `{invalid}`. Valid names are: `{sorted(old_names)}`")

        names_after_renaming = [new_names.get(n, n) for n in old_names]
        if len(set(names_after_renaming)) != len(old_names):
            raise ValueError(f"After renaming, {which} states will no longer unique: `{names_after_renaming}`.")
        # fmt: on

        if backward:
            self._init_states = self._init_states.set(
                assignment=states.cat.rename_categories(new_names)
            )
            memberships = self._init_states.memberships
            self._write_states(
                which,
                states=self.initial_states,
                colors=self._init_states.colors,
                probs=self.initial_states_probabilities,
                log=False,
            )
        else:
            self._term_states = self._term_states.set(
                assignment=states.cat.rename_categories(new_names)
            )
            memberships = self._term_states.memberships
            self._write_states(
                which,
                states=self.terminal_states,
                colors=self._term_states.colors,
                probs=self.terminal_states_probabilities,
                log=False,
            )
        if memberships is not None:
            memberships.names = [new_names.get(n, n) for n in memberships.names]

        return self

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
                clusters = clusters.cat.rename_categories({c: str(c) for c in clusters.cat.categories})
                vals = tuple(str(v) for v in vals)
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
            series_reference = series_reference.cat.rename_categories(
                {c: str(c) for c in series_reference.cat.categories}
            )

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
    def _write_states(
        self,
        which: Literal["initial", "terminal"],
        states: Optional[pd.Series],
        colors: Optional[np.ndarray],
        probs: Optional[pd.Series] = None,
        params: Dict[str, Any] = MappingProxyType({}),
    ) -> str:
        # fmt: off
        backward = which == "initial"
        key = Key.obs.term_states(backward)
        self._set(obj=self.adata.obs, key=key, value=states)
        self._set(obj=self.adata.obs, key=Key.obs.probs(key), value=probs)
        self._set(obj=self.adata.uns, key=Key.uns.colors(key), value=colors)
        if backward:
            self._init_states = self._init_states.set(assignment=states, probs=probs, colors=colors)
        else:
            self._term_states = self._term_states.set(assignment=states, probs=probs, colors=colors)
        self.params[key] = dict(params)
        # fmt: on

        return (
            f"Adding `adata.obs[{key!r}]`\n"
            f"       `adata.obs[{Key.obs.probs(key)!r}]`\n"
            f"       `.{which}_states`\n"
            f"       `.{which}_states_probabilities`\n"
            "    Finish"
        )

    def _read_from_adata(self, adata: AnnData, **kwargs: Any) -> bool:
        ok = super()._read_from_adata(adata, **kwargs)
        if not ok:
            return False

        # fmt: off
        for backward in [True, False]:
            key = Key.obs.term_states(backward)
            with SafeGetter(self, allowed=KeyError) as sg:
                assignment = self._get(obj=self.adata.obs, key=key, where="obs", dtype=pd.Series)
                probs = self._get(obj=self.adata.obs, key=Key.obs.probs(key), where="obs", dtype=pd.Series)
                colors = self._get(obj=self.adata.uns, key=Key.uns.colors(key), where="uns",
                                   dtype=(list, tuple, np.ndarray))
                colors = np.asarray([to_hex(c) for c in colors])
                if backward:
                    self._init_states = StatesHolder(assignment=assignment, probs=probs, colors=colors)
                else:
                    self._term_states = StatesHolder(assignment=assignment, probs=probs, colors=colors)
                self.params[key] = self._read_params(key)
        # fmt: on

        # status is based on `backward=False` by design
        return sg.ok

    @d.dedent
    def compute_states(self, *args: Any, **kwargs: Any) -> "TermStatesEstimator":
        """
        Compute initial or terminal states of the process.

        This is an alias for :meth:`predict`.

        Parameters
        ----------
        args
            Positional arguments.
        kwargs
            Keyword arguments arguments.

        Return
        ------
        Self and just updates the following fields: TODO

            - :attr:`terminal_states` - %(tse_term_states.summary)s
            - :attr:`terminal_states_probabilities` - %(tse_term_states_probs.summary)s
        """
        return self.predict(*args, **kwargs)

    plot_states = register_plotter(fwd_attr="_term_states", bwd_attr="_init_states")

    def _format_params(self) -> str:
        fmt = super()._format_params()
        ts = (
            None
            if self.terminal_states is None
            else sorted(self.terminal_states.cat.categories)
        )
        return fmt + f", terminal_states={ts}"
