from typing import Any, Dict, Tuple, Union, Mapping, Callable, Optional, Sequence
from typing_extensions import Literal, Protocol

from enum import auto
from wrapt import decorator

import scvelo as scv
from anndata import AnnData
from cellrank import logging as logg
from scanpy._utils import deprecated_arg_names
from cellrank.tl._enum import ModeEnum
from cellrank.ul._docs import d, inject_docs
from cellrank.tl._utils import RandomKeys, _unique_order_preserving
from cellrank.tl._colors import _create_categorical_colors
from cellrank.tl._lineage import Lineage

import numpy as np
import pandas as pd
from pandas.api.types import infer_dtype, is_categorical_dtype


class PlotMode(ModeEnum):  # noqa: D101
    EMBEDDING = auto()
    TIME = auto()


class BaseProtocol(Protocol):  # noqa: D101
    @property
    def adata(self) -> AnnData:  # noqa: D102
        ...

    @property
    def backward(self) -> bool:  # noqa: D102
        ...

    @property
    def params(self) -> Dict[str, Any]:  # noqa: D102
        ...

    def _set(
        self,
        attr: Optional[str] = None,
        obj: Optional[Union[pd.DataFrame, Mapping[str, Any]]] = None,
        key: Optional[str] = None,
        value: Optional[
            Union[np.ndarray, pd.Series, pd.DataFrame, Lineage, AnnData, Dict[str, Any]]
        ] = None,
        copy: bool = True,
        shadow_only: bool = False,
    ) -> None:
        ...

    def _get(
        self,
        attr: str,
        obj: Union[pd.DataFrame, Mapping[str, Any]],
        key: str,
        where: Optional[Literal["obs", "obsm", "var", "varm", "uns"]] = None,
        dtype: Optional[Union[type, Tuple[type, ...]]] = None,
        copy: bool = True,
        allow_missing: bool = False,
    ) -> None:
        ...

    def _create_params(
        self,
        locs: Optional[Mapping[str, Any]] = None,
        func: Optional[Callable] = None,
        remove: Sequence[str] = (),
    ) -> Dict[str, Any]:
        ...

    def _read_params(self, key: str) -> Dict[str, Any]:
        ...


class PlotterProtocol:  # noqa: D101
    @property
    def adata(self) -> AnnData:  # noqa: D102
        ...


@decorator()
def logger(
    wrapped: Callable[..., str], instance: Any, args: Any, kwargs: Dict[str, Any]
) -> str:
    """Handle logging for :class:`anndata.AnnData` writing functions of :class:`cellrank.estimators.BaseEstimator`."""
    log, time = kwargs.pop("log", True), kwargs.pop("time", None)
    msg = wrapped(*args, **kwargs)

    if log:
        logg.info(msg, time=time)

    return msg


@decorator()
def shadow(
    wrapped: Callable[..., str], instance: Any, args: Any, kwargs: Mapping[str, Any]
) -> str:
    """
    Duplicate function call with shadowed :class:`anndata.AnnData.object`.

    Used to easily create :class:`anndata.AnnData` serialization object in :class:`cellrank.estimators.BaseEstimator`.
    """
    res = wrapped(*args, **kwargs)
    with instance._shadow:
        try:
            # don't care what the "shadowed" function returns
            _ = wrapped(*args, **kwargs)
        except Exception as e:  # noqa: B902
            logg.error(
                f"Unable to duplicate function call using shadow `anndata.AnnData` object. Reason: `{e}`"
            )

    return res


def _plot_discrete(
    self: PlotterProtocol,
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
            f"Expected `data` to be of type `pandas.Series`, found `{type(_data).__name__}`."
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
    self: PlotterProtocol,
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
            f"Expected data to be of type `Lineage`, found `{type(_data).__name__}`."
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


# TODO(michalk8): in 2.0, remove %(time_mode)s and the deprecation
@d.dedent
@inject_docs(m=PlotMode)
@deprecated_arg_names({"cluster_key": "color"})
def _plot_dispatcher(
    self: PlotterProtocol,
    states: Optional[Union[str, Sequence[str]]] = None,
    color: Optional[str] = None,
    discrete: bool = False,
    mode: Literal["embedding", "time"] = PlotMode.EMBEDDING,
    time_key: str = "latent_time",
    same_plot: bool = True,
    title: Optional[Union[str, Sequence[str]]] = None,
    cmap: str = "viridis",
    **kwargs: Any,
) -> None:
    """
    Plot continuous or categorical observations in an embedding or along pseudotime.

    Parameters
    ----------
    states
        States to plot.
    color
        Key in :attr:`anndata.AnnData.obs`.
    discrete
        Whether to plot the data as continuous or discrete observations.
        If the data cannot be plotted as continuous observations, it will be plotted as discrete.
    %(time_mode)s
    time_key
        Key in :attr:`anndata.AnnData.obs` where pseudotime is stored. Only used when ``mode = {m.TIME!r}``.
    title
        Title of the plot(s).
    same_plot
        Whether to plot the data on the same plot or not. Only use when ``mode = {m.EMBEDDING!r}``.
        If `True` and ``discrete = False``, ``color`` is ignored.
    cmap
        Colormap for continuous data.
    kwargs
        Keyword arguments for :func:`scvelo.pl.scatter`.

    Returns
    -------
    %(just_plots)s
    """
    if discrete:
        return _plot_discrete(
            self,
            states=states,
            color=color,
            same_plot=same_plot,
            title=title,
            cmap=cmap,
            **kwargs,
        )

    return _plot_continuous(
        self,
        states=states,
        color=color,
        mode=mode,
        time_key=time_key,
        same_plot=same_plot,
        title=title,
        cmap=cmap,
        **kwargs,
    )


def register_plotter(
    *,
    discrete: Optional[str] = None,
    continuous: Optional[str] = None,
    colors: Optional[str] = None,
):
    """
    Register a plotting function.

    Used by the subclasses of :class`cellrank.estimators.BaseEstimator`.

    Parameters
    ----------
    discrete
        Attribute which stores categorical :class:`pandas.Series`.
    continuous
        Attribute which stores :class:`cellrank.tl.Lineage` containing continuous observations, such as
        macrostates memberships or absorption probabilities.
    colors
        Attribute where color is stored. Useful only when specifying ``discrete``,
        since :class:`cellrank.tl.Lineage` already contains color information.

    Returns
    -------
    The plotting function.
    """

    @decorator()
    def wrapper(
        wrapped: Callable[..., None],
        instance: Protocol,
        args: Any,
        kwargs: Dict[str, Any],
    ) -> None:
        # fmt: off
        disc: Optional[bool] = kwargs.pop("discrete", None)
        # prefer `continuous`
        disc = continuous is None if disc is None else disc

        if disc and discrete is None:
            logg.warning(f"Unable to plot `.{continuous}` in `discrete` mode. Using `continuous` mode")
            disc = False
        if not disc and continuous is None:
            logg.warning(f"Unable to plot `.{discrete}` in `continuous` mode. Using `discrete` mode")
            disc = True
        # fmt: on

        attr = discrete if disc else continuous
        data = getattr(instance, attr, None)
        if data is None:
            # prefer `continuous` for a prettier format
            attr = continuous if continuous is not None else discrete
            raise RuntimeError(f"Compute `.{attr}` first as `.compute_{attr}()`.")

        if colors is None:
            # extract colors from Lineage object, if present
            _colors = getattr(getattr(instance, continuous, None), "colors", None)
        else:
            _colors = getattr(instance, colors, None)

        return wrapped(
            *args,
            _data=data,
            _colors=_colors,
            _title=attr,
            discrete=disc,
            **kwargs,
        )

    if discrete is None and continuous is None:
        raise ValueError("At least 1 of `discrete` or `continuous` must be set.")

    return wrapper(_plot_dispatcher)
