from typing import Any, Union, Mapping, Callable, Optional, Sequence
from typing_extensions import Literal, Protocol

from wrapt import decorator

import scvelo as scv
from anndata import AnnData
from cellrank import logging as logg
from cellrank.tl import Lineage
from cellrank.tl._utils import RandomKeys, _unique_order_preserving
from cellrank.tl._colors import _create_categorical_colors
from cellrank.tl._estimators.mixins._constants import Key

import numpy as np
import pandas as pd
from pandas.api.types import infer_dtype, is_categorical_dtype


@decorator()
def logger(
    wrapped: Callable[..., str], instance: Any, args: Any, kwargs: Mapping[str, Any]
) -> str:
    log, time = kwargs.pop("log", True), kwargs.pop("time", None)
    msg = wrapped(*args, **kwargs)

    if log:
        logg.info(msg, time=time)

    return msg


@decorator()
def shadow(
    wrapped: Callable[..., str], instance: Any, args: Any, kwargs: Mapping[str, Any]
) -> None:
    res = wrapped(*args, **kwargs)
    with instance._shadow:
        _ = wrapped(*args, **kwargs)

    return res


class PlotterProtocol:
    @property
    def adata(self) -> AnnData:
        ...

    @property
    def backward(self) -> bool:
        ...


def _plot_discrete(
    self: PlotterProtocol,
    _data: pd.Series,
    _colors: Optional[np.ndarray] = None,
    lineages: Optional[Union[str, Sequence[str]]] = None,
    color: Optional[str] = None,
    title: Optional[Union[str, Sequence[str]]] = None,
    same_plot: bool = True,
    **kwargs: Any,
) -> None:
    if not isinstance(_data, pd.Series):
        raise TypeError(
            f"Expected `data` to be of type `pandas.Series`, found `{type(_data).__name__!r}`."
        )
    if not is_categorical_dtype(_data):
        raise TypeError(
            f"Expected `data` to be categorical, found `{infer_dtype(_data).__name__!r}`."
        )

    names = list(_data.cat.categories)
    if _colors is None:
        _colors = _create_categorical_colors(len(names))
    if len(_colors) != len(names):
        raise ValueError(
            f"Expected `colors` to be of length `{len(names)}`, found `{len(_colors)}`."
        )
    color_mapper = dict(zip(names, _colors))

    lineages = _unique_order_preserving(lineages or names)
    if not len(lineages):
        raise ValueError("No states have been selected.")

    for name in lineages:
        if name not in names:
            raise ValueError(
                f"Invalid name `{name!r}`. Valid options are: `{sorted(names)}`."
            )
    _data = _data.cat.set_categories(lineages)

    color = [] if color is None else (color,) if isinstance(color, str) else color
    color = _unique_order_preserving(color)

    same_plot = same_plot or len(names) == 1
    kwargs.setdefault("legend_loc", "on data")

    # fmt: off
    with RandomKeys(self.adata, n=1 if same_plot else len(lineages), where="obs") as keys:
        if same_plot:
            self.adata.obs[keys[0]] = _data
            self.adata.uns[f"{keys[0]}_colors"] = [color_mapper[name] for name in lineages]
            title = Key.obs.term_states(self.backward) if title is None else title
        else:
            for key, cat in zip(keys, lineages):
                self.adata.obs[key] = _data.cat.set_categories([cat])
                self.adata.uns[f"{key}_colors"] = [color_mapper[cat]]
            title = [f"{Key.initial(self.backward)} state {name}" for name in lineages] if title is None else title

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
    lineages: Optional[Union[str, Sequence[str]]] = None,
    color: Optional[str] = None,
    mode: Literal["embedding", "time"] = "embedding",
    time_key: str = "latent_time",
    title: Optional[Union[str, Sequence[str]]] = None,
    same_plot: bool = True,
    cmap: str = "viridis",
    **kwargs: Any,
) -> None:
    # TODO: cluster_key -> color
    if mode not in ("embedding", "time"):
        raise ValueError(
            f"Invalid mode `{mode!r}`. Valid options are `embedding` or `time`."
        )
    if not isinstance(_data, Lineage):
        raise TypeError(
            f"Expected the supplied data to be of type `Lineage`, "
            f"found `{type(_data).__name__!r}`."
        )

    if lineages is None:
        lineages = _data.names
    if not len(lineages):
        raise ValueError("No lineages have been selected.")
    is_singleton = _data.shape[1] == 1
    _data = _data[lineages]

    if mode == "time" and same_plot:
        logg.warning(
            "Invalid combination `mode='time'` and `same_plot=True`. Using `same_plot=False`"
        )
        same_plot = False

    _data_X = _data.X  # list(_data.T) behaves differently than a numpy.array
    if _data_X.shape[1] == 1:
        same_plot = False  # color grad looks empty
        if np.allclose(_data_X, 1.0):
            # matplotlib shows even tiny perturbations in the colormap
            _data_X = np.ones_like(_data_X)

    # TODO: reintroduce max scale? if yes, copy _data
    # fmt: off
    color = [] if color is None else (color,) if isinstance(color, str) else color
    color = _unique_order_preserving(color)

    if mode == "time":
        kwargs.setdefault("legend_loc", "best")
        if title is None:
            title = [f"{Key.where(self.backward)} {lin}" for lin in lineages]
        if time_key not in self.adata.obs:
            raise KeyError(f"Unable to find pseudotime in `adata.obs[{time_key!r}]`.")
        if len(color) and len(color) not in (1, _data_X.shape[1]):
            raise ValueError(f"Expected `color` to be of length `1` or `{_data_X.shape[1]}`, "
                             f"found `{len(color)}`.")
        kwargs["x"] = self.adata.obs[time_key]
        kwargs["y"] = list(_data_X.T)
        kwargs["color"] = color if len(color) else None
        kwargs["xlabel"] = [time_key] * len(lineages)
        kwargs["ylabel"] = ["probability"] * len(lineages)
    elif mode == "embedding":
        kwargs.setdefault("legend_loc", "on data")
        if same_plot:
            # TOOO: create an scvelo issue
            if color:
                logg.warning("Ignoring `cluster_key` when `mode='embedding'` and `same_plot=True`")
            title = [Key.initial(self.backward)] if title is None else title
            kwargs["color_gradients"] = _data
        else:
            title = [f"{Key.where(self.backward)} {lin}" for lin in lineages] if title is None else title
            title = color + title
            kwargs["color"] = color + list(_data_X.T)
    # fmt: on

    # TODO: can this even happen now?
    # e.g. a stationary distribution
    if is_singleton and not np.allclose(_data_X, 1.0):
        # TODO: maybe use always?
        kwargs.setdefault("perc", [0, 95])
        kwargs["color"] = _data_X
        _ = kwargs.pop("color_gradients", None)

    scv.pl.scatter(
        self.adata,
        title=title,
        color_map=cmap,
        **kwargs,
    )


def _plot_dispatcher(
    self: PlotterProtocol,
    *args: Any,
    discrete: Optional[bool] = None,
    **kwargs: Any,
) -> None:
    """TODO."""
    # TODO: expand kwargs for autocomplete
    if discrete:
        return _plot_discrete(self, *args, **kwargs)
    return _plot_continuous(self, *args, **kwargs)


def register_plotter(
    *,
    discrete: Optional[str] = None,
    continuous: Optional[str] = None,
    colors: Optional[str] = None,
):
    @decorator()
    def wrapper(
        wrapped: Callable[..., None],
        instance: Protocol,
        args: Any,
        kwargs: Mapping[str, Any],
    ) -> None:
        # fmt: off
        disc: Optional[bool] = kwargs.pop("discrete", None)
        disc = continuous is None if disc is None else disc

        if disc and discrete is None:
            # TODO: maybe lower the level
            logg.warning(f"Unable to plot `.{continuous}` in `discrete` mode. Using `continuous` mode")
            disc = False
        if not disc and continuous is None:
            logg.warning(f"Unable to plot `.{discrete}` in `continuous` mode. Using `discrete` mode")
            disc = True
        # fmt: on

        attr = discrete if disc else continuous
        data = getattr(instance, attr)
        if data is None:
            # prefer `continuous` for a prettier format
            attr = continuous if continuous is not None else discrete
            raise RuntimeError(f"Compute `.{attr}` first as `.compute_{attr}()`.")

        return wrapped(
            *args,
            _data=data,
            _colors=getattr(instance, colors, None)
            if isinstance(colors, str)
            else None,
            discrete=disc,
            **kwargs,
        )

    if discrete is None and continuous is None:
        raise ValueError("At least 1 of `discrete` or  `continuous` must be set.")

    return wrapper(_plot_dispatcher)
