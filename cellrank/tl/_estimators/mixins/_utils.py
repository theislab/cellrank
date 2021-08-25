from typing import Any, Union, Callable, Optional, Sequence
from typing_extensions import Protocol

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
from pandas.api.types import is_categorical_dtype


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
    cluster_key: Optional[str] = None,
    title: Optional[Union[str, Sequence[str]]] = None,
    same_plot: bool = True,
    **kwargs: Any,
) -> None:
    if not isinstance(_data, pd.Series):
        raise TypeError("TODO")
    if not is_categorical_dtype(_data):
        raise TypeError("TODO")

    names = list(_data.cat.categories)
    if _colors is None:
        _colors = _create_categorical_colors(len(names))
    if len(_colors) != len(names):
        raise ValueError("TODO")
    color_mapper = dict(zip(names, _colors))

    # these are states per-se, but I want to keep the arg names for dispatch the same
    if lineages is None:
        lineages = names
    lineages = _unique_order_preserving(lineages)
    if not len(lineages):
        raise ValueError("No states have been selected.")

    for lin in lineages:
        if lin not in names:
            raise ValueError(
                f"Invalid name `{lin!r}`. Valid options are: `{list(names)}`."
            )
    _data = _data.cat.set_categories(lineages)

    if cluster_key is None:
        cluster_key = []
    elif isinstance(cluster_key, str):
        cluster_key = [cluster_key]
    cluster_key = _unique_order_preserving(cluster_key)

    same_plot = same_plot or len(names) == 1
    kwargs["legend_loc"] = kwargs.get("legend_loc", "on data")

    with RandomKeys(
        self.adata, n=1 if same_plot else len(lineages), where="obs"
    ) as keys:
        if same_plot:
            self.adata.obs[keys[0]] = _data
            self.adata.uns[f"{keys[0]}_colors"] = [
                color_mapper[name] for name in lineages
            ]

            if title is None:
                title = Key.obs.term_states(self.backward)
            if isinstance(title, str):
                title = [title]
        else:
            for key, cat in zip(keys, lineages):
                self.adata.obs[key] = _data.cat.set_categories([cat])
                self.adata.uns[f"{key}_colors"] = [color_mapper[cat]]

            if title is None:
                title = [
                    f"{Key.initial(self.backward)} state {name}" for name in lineages
                ]

        scv.pl.scatter(
            self.adata,
            color=cluster_key + keys,
            title=cluster_key + title,
            **kwargs,
        )


def _plot_continuous(
    self: PlotterProtocol,
    _data: Lineage,
    _colors: Optional[np.ndarray] = None,
    lineages: Optional[Union[str, Sequence[str]]] = None,
    cluster_key: Optional[str] = None,
    title: Optional[Union[str, Sequence[str]]] = None,
    same_plot: bool = True,
    **kwargs: Any,
) -> None:
    raise NotImplementedError("TODO")


def _plot_dispatcher(
    self: PlotterProtocol,
    *args: Any,
    discrete: bool = False,
    **kwargs: Any,
) -> None:
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
        wrapped: Callable[..., None], instance: Protocol, args: Any, kwargs: Any
    ) -> None:
        disc = kwargs.pop("discrete", False)
        if disc and discrete is None:
            logg.warning(
                f"Unable to plot `.{continuous}` in `discrete` mode. Using `continuous` mode"
            )
            disc = False
        if not disc and continuous is None:
            logg.warning(
                f"Unable to plot `.{discrete}` in `continuous` mode. Using `discrete` mode"
            )
            disc = True

        attr = discrete if disc else continuous
        data = getattr(instance, attr)
        if data is None:
            raise RuntimeError(f"Compute `.{attr}` first as `.compute_{attr}()`.")

        return wrapped(
            *args,
            _data=data,
            _colors=getattr(instance, colors, None),
            discrete=discrete,
            **kwargs,
        )

    if discrete is None and continuous is None:
        raise ValueError("TODO.")

    return wrapper(_plot_dispatcher)
