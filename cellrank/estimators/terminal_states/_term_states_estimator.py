from typing import Any, Dict, Tuple, Union, Literal, Mapping, Optional, Sequence

from abc import ABC
from types import MappingProxyType

import scvelo as scv
from anndata import AnnData
from cellrank import logging as logg
from cellrank._utils._key import Key
from cellrank._utils._docs import d, inject_docs
from cellrank._utils._utils import (
    RandomKeys,
    _unique_order_preserving,
    _merge_categorical_series,
    _convert_to_categorical_series,
)
from cellrank._utils._colors import (
    _map_names_and_colors,
    _convert_to_hex_colors,
    _create_categorical_colors,
)
from cellrank._utils._lineage import Lineage
from cellrank.kernels._base_kernel import KernelExpression
from cellrank.estimators.mixins._utils import (
    PlotMode,
    SafeGetter,
    StatesHolder,
    logger,
    shadow,
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
        """Categorical annotation of terminal states.

        By default, all transient cells will be labeled as `NaN`.
        """
        return self._term_states.assignment

    @property
    @d.get_summary(base="tse_term_states_probs")
    def terminal_states_probabilities(self) -> Optional[pd.Series]:
        """Aggregated probability of cells to be in terminal states."""
        return self._term_states.probs

    @property
    @d.get_summary(base="tse_init_states")
    def initial_states(self) -> Optional[pd.Series]:
        """Categorical annotation of initial states.

        By default, all transient cells will be labeled as `NaN`.
        """
        return self._init_states.assignment

    @property
    @d.get_summary(base="tse_init_states_probs")
    def initial_states_probabilities(self) -> Optional[pd.Series]:
        """Aggregated probability of cells to be in initial states."""
        return self._init_states.probs

    @d.dedent
    def set_states(
        self,
        labels: Union[pd.Series, Dict[str, Sequence[Any]]],
        which: Literal["initial", "terminal"] = "terminal",
        cluster_key: Optional[str] = None,
        allow_overlap: bool = False,
        **kwargs: Any,
    ) -> "TermStatesEstimator":
        """
        Manually define initial or terminal states.

        Parameters
        ----------
        labels
            Defines the states. Valid options are:

                - categorical :class:`pandas.Series` where each category corresponds to an individual state.
                  `NaN` entries denote cells that do not belong to any state, i.e., transient cells.
                - :class:`dict` where keys are states and values are lists of cell barcodes corresponding to
                  annotations in :attr:`anndata.AnnData.obs_names`.
                  If only 1 key is provided, values should correspond to clusters if a categorical
                  :class:`pandas.Series` can be found in :attr:`anndata.AnnData.obs`.
        which
            Whether to set initial or terminal states.
        cluster_key
            Key in :attr:`anndata.AnnData.obs` in order to associate names and colors with :attr:`terminal_states` or
            :attr:`initial_states`. Each state will be given the name and color corresponding to the cluster it
            mostly overlaps with.
        %(allow_overlap)s

        Returns
        -------
        If ``which = 'terminal'``, returns self and updates the following fields:

            - :attr:`terminal_states` - %(tse_term_states.summary)s
            - :attr:`terminal_states_probabilities` - %(tse_term_states_probs.summary)s

        Otherwise, returns self and updates the following fields:

            - :attr:`initial_states` - %(tse_init_states.summary)s
            - :attr:`initial_states_probabilities` - %(tse_init_states_probs.summary)s
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
            allow_overlap=allow_overlap,
            **kwargs,
        )
        return self

    @d.get_sections(base="tse_rename_term_states", sections=["Parameters", "Returns"])
    @d.get_full_description(base="tse_rename_term_states")
    @d.dedent
    def rename_states(  # TODO(michalk8): override in GPCCA to handle macrostates
        self,
        old_new: Mapping[str, str],
        which: Literal["initial", "terminal"] = "terminal",
    ) -> "TermStatesEstimator":
        """
        Rename :attr:`terminal_states` or :attr:`initial_states`.

        Parameters
        ----------
        old_new
            Mapping where keys correspond to the old names and the values to the new names.
            The new names must be unique.
        %(which)s

        Returns
        -------
        If ``which = 'terminal'``, returns self and updates the following field:

            - :attr:`terminal_states` - %(tse_term_states.summary)s

        Otherwise, returns self and updates the following field:

            - :attr:`initial_states` - %(tse_init_states.summary)s
        """
        backward = which == "initial"
        states = self.initial_states if backward else self.terminal_states
        if states is None:
            raise RuntimeError(
                f"Compute {which} states first as `.predict_{which}_states(...)` or "
                f"set them manually as `.set_states(..., which={which!r})`."
            )

        # fmt: off
        if not isinstance(old_new, Mapping):
            raise TypeError(f"Expected new names to be a `Mapping`, found `{type(old_new)}`.")
        if not len(old_new):
            return self

        old_names = states.cat.categories
        old_new = {str(k): str(v) for k, v in old_new.items()}
        mask = np.isin(list(old_new.keys()), old_names)
        if not np.all(mask):
            invalid = sorted(np.array(list(old_new.keys()))[~mask])
            raise ValueError(f"Invalid {which} states names: `{invalid}`. Valid names are: `{sorted(old_names)}`.")

        names_after_renaming = [old_new.get(n, n) for n in old_names]
        if len(set(names_after_renaming)) != len(old_names):
            raise ValueError(f"After renaming, {which} states will no longer unique: `{names_after_renaming}`.")
        # fmt: on

        if backward:
            assignment = states.cat.rename_categories(old_new)
            memberships = self._init_states.memberships
            self._write_states(
                which,
                states=assignment,
                colors=self._init_states.colors,
                probs=self.initial_states_probabilities,
                log=False,
            )
            self._init_states = self._init_states.set(assignment=assignment)
        else:
            assignment = states.cat.rename_categories(old_new)
            memberships = self._term_states.memberships
            self._write_states(
                which,
                states=assignment,
                colors=self._term_states.colors,
                probs=self.terminal_states_probabilities,
                log=False,
            )
            self._term_states = self._term_states.set(assignment=assignment)
        if memberships is not None:
            memberships.names = [old_new.get(n, n) for n in memberships.names]

        return self

    @d.dedent
    @inject_docs(m=PlotMode)
    def plot_macrostates(
        self,
        which: Literal["macro", "initial", "terminal"] = "terminal",
        states: Optional[Union[str, Sequence[str]]] = None,
        color: Optional[str] = None,
        discrete: bool = True,
        mode: Literal["embedding", "time"] = PlotMode.EMBEDDING,
        time_key: str = "latent_time",
        same_plot: bool = True,
        title: Optional[Union[str, Sequence[str]]] = None,
        cmap: str = "viridis",
        **kwargs: Any,
    ) -> None:
        """Plot macrostates on an embedding or along pseudotime.

        Parameters
        ----------
        which
            Which type of macrostates to plot. Valid options are:

                - ``'macro'`` - plot the macrostates.
                - ``'initial'`` - plot the macrostates marked as initial states.
                - ``'terminal'`` - plot the macrostates marked as terminal terminal.
        states
            Subset of the macrostates to show. If ``None``, plot all macrostates.
        color
            Key in :attr:`anndata.AnnData.obs` or :attr:`anndata.AnnData.var` used to color the observations.
        discrete
            Whether to plot the data as continuous or discrete observations.
            If the data cannot be plotted as continuous observations, it will be plotted as discrete.
        time_key
            Key in :attr:`anndata.AnnData.obs` where pseudotime is stored. Only used when ``mode = {m.TIME!r}``.
        title
            Title of the plot.
        same_plot
            Whether to plot the data on the same plot or not. Only use when ``mode = {m.EMBEDDING!r}``.
            If `True` and ``discrete = False``, ``color`` is ignored.
        cmap
            Colormap for continuous annotations.
        kwargs
            Keyword arguments for :func:`scvelo.pl.scatter`.

        Returns
        -------
        %(just_plots)s
        """
        if which == "macro":
            obj: StatesHolder = self._macrostates
        elif which == "initial":
            obj = self._init_states
        elif which == "terminal":
            obj = self._term_states
        else:
            raise ValueError(
                f"Unable to plot `{which!r}` states. "
                f"Valid options are: `{['macro', 'initial', 'terminal']}`."
            )

        name = "macrostates" if which == "macro" else f"{which} states"
        if obj.assignment is None and obj.memberships is None:
            raise RuntimeError(f"Compute {name} first.")

        if not discrete and obj.memberships is None:
            logg.warning(f"Unable to plot {name} in continuous mode, using discrete")
            discrete = True

        data = obj.assignment if discrete else obj.memberships
        colors = obj.colors

        if discrete:
            return self._plot_discrete(
                _data=data,
                _colors=colors,
                _title=name,
                states=states,
                color=color,
                same_plot=same_plot,
                title=title,
                cmap=cmap,
                **kwargs,
            )
        return self._plot_continuous(
            _data=data,
            _colors=colors,
            _title=name,
            states=states,
            color=color,
            mode=mode,
            time_key=time_key,
            same_plot=same_plot,
            title=title,
            cmap=cmap,
            **kwargs,
        )

    def _plot_discrete(
        self,
        _data: pd.Series,
        _colors: Optional[np.ndarray] = None,
        _title: Optional[str] = None,
        states: Optional[Union[str, Sequence[str]]] = None,
        color: Optional[str] = None,
        title: Optional[Union[str, Sequence[str]]] = None,
        same_plot: bool = True,
        cmap: str = "viridis",
        **kwargs: Any,
    ) -> None:
        if not isinstance(_data, pd.Series):
            raise TypeError(
                f"Expected `data` to be of type `pandas.Series`, found `{type(_data)}`."
            )
        if not is_categorical_dtype(_data):
            raise TypeError(
                f"Expected `data` to be `categorical`, found `{infer_dtype(_data)}`."
            )

        names = list(_data.cat.categories)
        if _colors is None:
            _colors = _create_categorical_colors(len(names))
        if len(_colors) != len(names):
            raise ValueError(
                f"Expected `colors` to be of length `{len(names)}`, found `{len(_colors)}`."
            )
        color_mapper = dict(zip(names, _colors))

        states = _unique_order_preserving(states or names)
        if not len(states):
            raise ValueError("No states have been selected.")

        for name in states:
            if name not in names:
                raise ValueError(
                    f"Invalid name `{name!r}`. Valid options are: `{sorted(names)}`."
                )
        _data = _data.cat.set_categories(states)

        color = [] if color is None else (color,) if isinstance(color, str) else color
        color = _unique_order_preserving(color)

        same_plot = same_plot or len(names) == 1
        kwargs.setdefault("legend_loc", "on data")
        kwargs["color_map"] = cmap

        # fmt: off
        with RandomKeys(self.adata, n=1 if same_plot else len(states), where="obs") as keys:
            if same_plot:
                self.adata.obs[keys[0]] = _data
                self.adata.uns[f"{keys[0]}_colors"] = [color_mapper[name] for name in states]
                title = _title if title is None else title
            else:
                for key, cat in zip(keys, states):
                    self.adata.obs[key] = _data.cat.set_categories([cat])
                    self.adata.uns[f"{key}_colors"] = [color_mapper[cat]]
                title = [f"{_title} {name}" for name in states] if title is None else title

            if isinstance(title, str):
                title = [title]

            scv.pl.scatter(
                self.adata,
                color=color + keys,
                title=color + title,
                **kwargs,
            )
        # fmt: on

    def _plot_continuous(
        self,
        _data: Lineage,
        _colors: Optional[np.ndarray] = None,
        _title: Optional[str] = None,
        states: Optional[Union[str, Sequence[str]]] = None,
        color: Optional[str] = None,
        mode: Literal["embedding", "time"] = PlotMode.EMBEDDING,
        time_key: str = "latent_time",
        title: Optional[Union[str, Sequence[str]]] = None,
        same_plot: bool = True,
        cmap: str = "viridis",
        **kwargs: Any,
    ) -> None:
        mode = PlotMode(mode)
        if not isinstance(_data, Lineage):
            raise TypeError(
                f"Expected data to be of type `Lineage`, found `{type(_data)}`."
            )

        if states is None:
            states = _data.names
        if not len(states):
            raise ValueError("No lineages have been selected.")
        is_singleton = _data.shape[1] == 1
        _data = _data[states].copy()

        if mode == "time" and same_plot:
            logg.warning(
                "Invalid combination `mode='time'` and `same_plot=True`. Using `same_plot=False`"
            )
            same_plot = False

        _data_X = _data.X  # list(_data.T) behaves differently than a numpy.array
        if _data_X.shape[1] == 1:
            same_plot = False
            if np.allclose(_data_X, 1.0):
                # matplotlib shows even tiny perturbations in the colormap
                _data_X = np.ones_like(_data_X)

        for col in _data_X.T:
            mask = ~np.isclose(col, 1.0)
            # change the maximum value - the 1 is artificial and obscures the color scaling
            if np.any(mask):
                col[~mask] = np.nanmax(col[mask])

        # fmt: off
        color = [] if color is None else (color,) if isinstance(color, str) else color
        color = _unique_order_preserving(color)

        if mode == PlotMode.TIME:
            kwargs.setdefault("legend_loc", "best")
            if title is None:
                title = [f"{_title} {state}" for state in states]
            if time_key not in self.adata.obs:
                raise KeyError(f"Unable to find pseudotime in `adata.obs[{time_key!r}]`.")
            if len(color) and len(color) not in (1, _data_X.shape[1]):
                raise ValueError(f"Expected `color` to be of length `1` or `{_data_X.shape[1]}`, "
                                 f"found `{len(color)}`.")
            kwargs["x"] = self.adata.obs[time_key]
            kwargs["y"] = list(_data_X.T)
            kwargs["color"] = color if len(color) else None
            kwargs["xlabel"] = [time_key] * len(states)
            kwargs["ylabel"] = ["probability"] * len(states)
        elif mode == PlotMode.EMBEDDING:
            kwargs.setdefault("legend_loc", "on data")
            if same_plot:
                if color:
                    # https://github.com/theislab/scvelo/issues/673
                    logg.warning("Ignoring `color` when `mode='embedding'` and `same_plot=True`")
                title = [_title] if title is None else title
                kwargs["color_gradients"] = _data
            else:
                title = [f"{_title} {state}" for state in states] if title is None else title
                if isinstance(title, str):
                    title = [title]
                title = color + title
                kwargs["color"] = color + list(_data_X.T)
        else:
            raise NotImplementedError(f"Mode `{mode}` is not yet implemented.")
        # fmt: on

        # e.g. a stationary distribution
        if is_singleton and not np.allclose(_data_X, 1.0):
            kwargs.setdefault("perc", [0, 95])
            _ = kwargs.pop("color_gradients", None)

        scv.pl.scatter(
            self.adata,
            title=title,
            color_map=cmap,
            **kwargs,
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
        allow_overlap: bool = False,
    ) -> str:
        # fmt: off
        backward = which == "initial"
        if not allow_overlap:
            fwd_states, bwd_states = (self.terminal_states, states) if backward else (states, self.initial_states)
            if fwd_states is not None and bwd_states is not None:
                overlap = np.sum(~pd.isnull(fwd_states) & ~pd.isnull(bwd_states))
                if overlap:
                    raise ValueError(
                        f"Found `{overlap}` overlapping cells between initial and terminal states. "
                        f"If this is intended, use `allow_overlap=True`."
                    )

        key = Key.obs.term_states(self.backward, bwd=backward)
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
            key = Key.obs.term_states(self.backward, bwd=backward)
            with SafeGetter(self, allowed=KeyError) as sg:
                assignment = self._get(obj=adata.obs, key=key, shadow_attr="obs", dtype=pd.Series)
                probs = self._get(obj=adata.obs, key=Key.obs.probs(key), shadow_attr="obs", dtype=pd.Series)
                colors = self._get(obj=adata.uns, key=Key.uns.colors(key), shadow_attr="uns",
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

    def _format_params(self) -> str:
        fmt = super()._format_params()
        # fmt: off
        init_states = None if self.initial_states is None else sorted(self.initial_states.cat.categories)
        term_states = None if self.terminal_states is None else sorted(self.terminal_states.cat.categories)
        return fmt + f", initial_states={init_states}, terminal_states={term_states}"
        # fmt: on
