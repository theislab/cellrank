import abc
import types
from typing import Any, Dict, Literal, Optional, Sequence, Tuple, Union

import scvelo as scv

import numpy as np
import pandas as pd
import scipy.sparse as sp
from pandas.api.types import infer_dtype

from matplotlib.colors import to_hex

from anndata import AnnData

from cellrank import logging as logg
from cellrank._utils._colors import (
    _convert_to_hex_colors,
    _create_categorical_colors,
    _map_names_and_colors,
)
from cellrank._utils._docs import d, inject_docs
from cellrank._utils._key import Key
from cellrank._utils._lineage import Lineage
from cellrank._utils._utils import (
    RandomKeys,
    _convert_to_categorical_series,
    _merge_categorical_series,
    _unique_order_preserving,
)
from cellrank.estimators._base_estimator import BaseEstimator
from cellrank.estimators.mixins._utils import (
    PlotMode,
    SafeGetter,
    StatesHolder,
    logger,
    shadow,
)
from cellrank.kernels._base_kernel import KernelExpression

__all__ = ["TermStatesEstimator"]


@d.dedent
class TermStatesEstimator(BaseEstimator, abc.ABC):
    """Base class for all estimators predicting the initial and terminal states.

    Parameters
    ----------
    %(base_estimator.parameters)s
    """

    def __init__(
        self,
        object: Union[AnnData, np.ndarray, sp.spmatrix, KernelExpression],
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
        """Probability to be a terminal state."""
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
        """Probability to be an initial state."""
        return self._init_states.probs

    @d.dedent
    def set_terminal_states(
        self,
        states: Union[pd.Series, Dict[str, Sequence[Any]]],
        cluster_key: Optional[str] = None,
        allow_overlap: bool = False,
        **kwargs: Any,
    ) -> "TermStatesEstimator":
        """Set the :attr:`terminal_states`.

        Parameters
        ----------
        states
            States to select. Valid options are:

            - categorical :class:`~pandas.Series` where each category corresponds to an individual state.
              `NaN` entries denote cells that do not belong to any state, i.e., transient cells.
            - :class:`dict` where keys are states and values are lists of cell barcodes corresponding to
              annotations in :attr:`~anndata.AnnData.obs_names`.
              If only 1 key is provided, values should correspond to clusters if a categorical
              :class:`~pandas.Series` can be found in :attr:`~anndata.AnnData.obs`.
        cluster_key
            Key in :attr:`~anndata.AnnData.obs` to associate names and colors with :attr:`terminal_states`.
            Each state will be given the name and color corresponding to the cluster it mostly overlaps with.
        %(allow_overlap)s
        kwargs
            Additional keyword arguments.

        Returns
        -------
        Returns self and updates the following fields:

        - :attr:`terminal_states` - %(tse_term_states.summary)s
        - :attr:`terminal_states_probabilities` - %(tse_term_states_probs.summary)s
        """
        states, colors = self._set_categorical_labels(
            categories=states,
            cluster_key=cluster_key,
            existing=None,
        )
        self._write_states(
            "terminal",
            states=states,
            colors=colors,
            allow_overlap=allow_overlap,
            **kwargs,
        )
        return self

    @d.dedent
    def set_initial_states(
        self,
        states: Union[pd.Series, Dict[str, Sequence[Any]]],
        cluster_key: Optional[str] = None,
        allow_overlap: bool = False,
        **kwargs: Any,
    ) -> "TermStatesEstimator":
        """Set the :attr:`initial_states`.

        Parameters
        ----------
        states
            Which states to select. Valid options are:

            - categorical :class:`~pandas.Series` where each category corresponds to an individual state.
              `NaN` entries denote cells that do not belong to any state, i.e., transient cells.
            - :class:`dict` where keys are states and values are lists of cell barcodes corresponding to
              annotations in :attr:`~anndata.AnnData.obs_names`.
              If only 1 key is provided, values should correspond to clusters if a categorical
              :class:`~pandas.Series` can be found in :attr:`~anndata.AnnData.obs`.
        cluster_key
            Key in :attr:`~anndata.AnnData.obs` to associate names and colors :attr:`initial_states`.
            Each state will be given the name and color corresponding to the cluster it mostly overlaps with.
        %(allow_overlap)s
        kwargs
            Additional keyword arguments.

        Returns
        -------
        Returns self and updates the following fields:

        - :attr:`initial_states` - %(tse_init_states.summary)s
        - :attr:`initial_states_probabilities` - %(tse_init_states_probs.summary)s
        """
        states, colors = self._set_categorical_labels(
            categories=states,
            cluster_key=cluster_key,
            existing=None,
        )
        self._write_states(
            "initial",
            states=states,
            colors=colors,
            allow_overlap=allow_overlap,
            **kwargs,
        )
        return self

    @d.get_sections(base="tse_rename_term_states", sections=["Parameters", "Returns"])
    @d.get_full_description(base="tse_rename_term_states")
    @d.dedent
    def rename_terminal_states(self, old_new: Dict[str, str]) -> "TermStatesEstimator":
        """Rename the :attr:`terminal_states`.

        Parameters
        ----------
        old_new
            Dictionary that maps old names to unique new names.

        Returns
        -------
        Returns self and updates the following fields:

        - :attr:`terminal_states` - %(tse_term_states.summary)s
        """
        states = self.terminal_states
        if states is None:
            raise RuntimeError(
                "Compute terminal states first as `.predict_terminal_states()` or "
                "set them manually as `.set_terminal_states()`."
            )

        # fmt: off
        if not isinstance(old_new, dict):
            raise TypeError(f"Expected new names to be a `dict`, found `{type(old_new)}`.")
        if not len(old_new):
            return self

        old_names = states.cat.categories
        mask = np.isin(list(old_new.keys()), old_names)
        if not np.all(mask):
            invalid = sorted(np.array(list(old_new.keys()))[~mask])
            raise ValueError(f"Invalid terminal states names: `{invalid}`. Valid names are: `{sorted(old_names)}`.")

        names_after_renaming = [old_new.get(n, n) for n in old_names]
        if len(set(names_after_renaming)) != len(old_names):
            raise ValueError(f"After renaming, terminal states will no longer unique: `{names_after_renaming}`.")
        # fmt: on

        assignment = states.cat.rename_categories(old_new)
        memberships = self._term_states.memberships  # save before overwriting
        self._write_states(
            "terminal",
            states=assignment,
            colors=self._term_states.colors,
            probs=self.terminal_states_probabilities,
            log=False,
        )
        if memberships is not None:
            memberships.names = [old_new.get(n, n) for n in memberships.names]

        self._term_states = self._term_states.set(assignment=assignment, memberships=memberships)
        return self

    @d.dedent
    def rename_initial_states(self, old_new: Dict[str, str]) -> "TermStatesEstimator":
        """Rename the :attr:`initial_states`.

        Parameters
        ----------
        old_new
            Dictionary that maps old names to unique new names.

        Returns
        -------
        Returns self and updates the following fields:

        - :attr:`initial_states` - %(tse_init_states.summary)s
        """
        states = self.initial_states
        if states is None:
            raise RuntimeError(
                "Compute initial states first as `.predict_initial_states()` or "
                "set them manually as `.set_initial_states()`."
            )

        # fmt: off
        if not isinstance(old_new, dict):
            raise TypeError(f"Expected new names to be a `dict`, found `{type(old_new)}`.")
        if not len(old_new):
            return self

        old_names = states.cat.categories
        mask = np.isin(list(old_new.keys()), old_names)
        if not np.all(mask):
            invalid = sorted(np.array(list(old_new.keys()))[~mask])
            raise ValueError(f"Invalid terminal states names: `{invalid}`. Valid names are: `{sorted(old_names)}`.")

        names_after_renaming = [old_new.get(n, n) for n in old_names]
        if len(set(names_after_renaming)) != len(old_names):
            raise ValueError(f"After renaming, terminal states will no longer unique: `{names_after_renaming}`.")
        # fmt: on

        assignment = states.cat.rename_categories(old_new)
        memberships = self._init_states.memberships  # save overwriting
        self._write_states(
            "initial",
            states=assignment,
            colors=self._init_states.colors,
            probs=self.initial_states_probabilities,
            log=False,
        )
        if memberships is not None:
            memberships.names = [old_new.get(n, n) for n in memberships.names]

        self._init_states = self._init_states.set(assignment=assignment, memberships=memberships)
        return self

    @d.dedent
    @inject_docs(m=PlotMode)
    def plot_macrostates(
        self,
        which: Literal["all", "initial", "terminal"],
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
            Which macrostates to plot. Valid options are:

            - ``'all'`` - plot all macrostates.
            - ``'initial'`` - plot macrostates marked as :attr:`initial_states`.
            - ``'terminal'`` - plot macrostates marked as :attr:`terminal_states`.
        states
            Subset of the macrostates to show. If :obj:`None`, plot all macrostates.
        color
            Key in :attr:`~anndata.AnnData.obs` or :attr:`~anndata.AnnData.var` used to color the observations.
        discrete
            Whether to plot the data as continuous or discrete observations.
            If the data cannot be plotted as continuous observations, it will be plotted as discrete.
        mode
            Whether to plot the probabilities in an embedding or along the pseudotime.
        time_key
            Key in :attr:`~anndata.AnnData.obs` where pseudotime is stored. Only used when ``mode = {m.TIME!r}``.
        title
            Title of the plot.
        same_plot
            Whether to plot the data on the same plot or not. Only use when ``mode = {m.EMBEDDING!r}``.
            If `True` and ``discrete = False``, ``color`` is ignored.
        cmap
            Colormap for continuous annotations.
        kwargs
            Keyword arguments for :func:`~scvelo.pl.scatter`.

        Returns
        -------
        %(just_plots)s
        """
        if which == "all":
            obj: Optional[StatesHolder] = getattr(self, "_macrostates", None)
            if obj is None:
                raise RuntimeError(f"`{type(self).__name__}` cannot plot macrostates.")
        elif which == "initial":
            obj = self._init_states
        elif which == "terminal":
            obj = self._term_states
        else:
            raise ValueError(
                f"Unable to plot `{which!r}` states. " f"Valid options are: `{['all', 'initial', 'terminal']}`."
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
            raise TypeError(f"Expected `data` to be of type `pandas.Series`, found `{type(_data)}`.")
        if not isinstance(_data.dtype, pd.CategoricalDtype):
            raise TypeError(f"Expected `data` to be `categorical`, found `{infer_dtype(_data)}`.")

        names = list(_data.cat.categories)
        if _colors is None:
            _colors = _create_categorical_colors(len(names))
        if len(_colors) != len(names):
            raise ValueError(f"Expected `colors` to be of length `{len(names)}`, found `{len(_colors)}`.")
        color_mapper = dict(zip(names, _colors))

        states = _unique_order_preserving(states or names)
        if not len(states):
            raise ValueError("No states have been selected.")

        for name in states:
            if name not in names:
                raise ValueError(f"Invalid name `{name!r}`. Valid options are: `{sorted(names)}`.")
        _data = _data.cat.set_categories(states)

        color = [] if color is None else (color,) if isinstance(color, str) else color
        color = _unique_order_preserving(color)

        same_plot = same_plot or len(names) == 1
        kwargs.setdefault("legend_loc", "on data")
        kwargs["color_map"] = cmap

        # fmt: off
        with RandomKeys(self.adata, n=1 if same_plot else len(states), where="obs") as keys:
            if same_plot:
                outline = _data.cat.categories.to_list()
                _data = _data.cat.add_categories(["nan"]).fillna("nan")
                states.append("nan")
                color_mapper["nan"] = "#dedede"
                self.adata.obs[keys[0]] = _data
                self.adata.uns[f"{keys[0]}_colors"] = [color_mapper[name] for name in states]
                title = _title if title is None else title
            else:
                outline = None
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
                add_outline=outline,
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
            raise TypeError(f"Expected data to be of type `Lineage`, found `{type(_data)}`.")

        if states is None:
            states = _data.names
        if not len(states):
            raise ValueError("No lineages have been selected.")
        is_singleton = _data.shape[1] == 1
        _data = _data[states].copy()

        if mode == "time" and same_plot:
            logg.warning("Invalid combination `mode='time'` and `same_plot=True`. Using `same_plot=False`")
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
            if time_key is None:
                raise KeyError(
                    "The name of the column in `adata.obs` defining the pseudotime needs to be defined via the "
                    "`time_key` argument."
                )
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
            data = self.adata.obs.get(key, None)
            is_categorical = data is not None and isinstance(data.dtype, pd.CategoricalDtype)
            if len(categories) == 1 and is_categorical:
                vals = categories[key]
                if isinstance(vals, str) or not isinstance(vals, Sequence):
                    vals = (categories[key],)

                clusters = self.adata.obs[key]
                clusters = clusters.cat.rename_categories({c: str(c) for c in clusters.cat.categories})
                vals = tuple(str(v) for v in vals)
                categories = {cat: self.adata[clusters == cat].obs_names for cat in vals}

            categories = _convert_to_categorical_series(categories, list(self.adata.obs_names))
        if not isinstance(categories.dtype, pd.CategoricalDtype):
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
            cats = categories.cat.categories
            categories = categories.cat.rename_categories(dict(zip(cats, names)))
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
        params: Dict[str, Any] = types.MappingProxyType({}),
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
                        f"If this is intended, please use `allow_overlap=True`."
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
