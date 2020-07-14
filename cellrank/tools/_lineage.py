# -*- coding: utf-8 -*-
"""Lineage class module."""

from types import FunctionType
from typing import List, Tuple, Union, TypeVar, Callable, Iterable, Optional
from inspect import signature
from pathlib import Path
from functools import wraps
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import entropy

import matplotlib.colors as c
import matplotlib.pyplot as plt

from cellrank import logging as logg
from cellrank.tools._utils import (
    save_fig,
    _convert_lineage_name,
    _unique_order_preserving,
)
from cellrank.tools._colors import (
    _get_bg_fg_colors,
    _compute_mean_color,
    _create_categorical_colors,
)
from cellrank.tools._constants import Lin

ColorLike = TypeVar("ColorLike")
_ERROR_NOT_ITERABLE = "Expected `{}` to be iterable, found type `{}`."
_ERROR_WRONG_SIZE = "Expected `{}` to be of size `{{}}`, found `{{}}`."

_HT_CELLS = 10  # head and tail cells to show
_HTML_REPR_THRESH = 100
_DUMMY_CELL = "<td style='text-align: right;'>...</td>"
_ORDER = "C"


def _at_least_2d(array: np.ndarray, dim: int):
    return np.expand_dims(array, dim) if array.ndim < 2 else array


def wrap(numpy_func: Callable) -> Callable:
    """
    Wrap an numpy function.

    Modifies functionality of some function (e.g. ignoring `.squeeze`, retaining dimensions).

    Params
    ------
    numpy_func
        Function to be wrapped.

    Returns
    -------
    :class:`Callable`
        Wrapped function which takes a :class:`cellrank.tl.Lineage` and return :class:`cellrank.tl.Lineage`.

    """

    @wraps(numpy_func)
    def decorator(array, *args, **kwargs):
        if not isinstance(array, Lineage):
            raise TypeError(
                f"Expected array to be of type `Lineage`, found `{type(array).__name__}`."
            )
        if fname == "squeeze":
            return array
        if fname == "array_repr":
            return repr(array)

        if "axis" in kwargs:
            axis = kwargs["axis"]
        elif axis_ix < len(args):
            axis = args[axis_ix]
        else:
            axis = default_axis

        res = np.array(numpy_func(array.X, *args, **kwargs), copy=False)

        # handle expand_dim
        if res.ndim > 2:
            return array
        if array.ndim < res.ndim:
            names = (
                array.names
                if hasattr(array, "names")
                and array.names is not None
                and len(array.names) == res.shape[1]
                else [f"Lineage {i}" for i in range(res.shape[1])]
            )
            return Lineage(res, names=names)

        # handle reductions
        if not res.shape:
            return Lineage(np.array([[res]]), names=[fname], colors=["grey"])
        if res.shape == array.shape:
            return Lineage(res, names=array.names, colors=array.colors)

        res = np.expand_dims(res, axis)

        is_t = int(array._is_transposed)
        if is_t:
            res = res.T

        lin = None
        if res.shape[0] == array.shape[is_t]:
            lin = Lineage(
                res, names=[f"{fname} of {', '.join(array.names)}"], colors=["grey"],
            )

        if res.shape[1] == array.shape[1 - is_t]:
            lin = Lineage(res, names=[f"{fname} of {n}" for n in array.names])

        if lin is not None:
            return lin.T if is_t else lin

        raise RuntimeError(
            f"Unable to interpret result of function `{fname}` called "
            f"with args `{args}`, kwargs: `{kwargs}`."
        )

    params = signature(numpy_func).parameters
    if "axis" in params:
        axis_ix = list(params.keys()).index("axis") - 1
        default_axis = params["axis"].default
    else:
        axis_ix = 256
        default_axis = None
    assert (
        axis_ix >= 0
    ), f"Expected argument `'axis'` not to be first for function `{numpy_func.__name__}`."

    fname = numpy_func.__name__
    if fname == "amin":
        fname = "min"
    elif fname == "amax":
        fname = "max"

    return decorator


def _register_handled_functions():
    handled_fns = {}

    for attrname in dir(np):
        fn = getattr(np, attrname)
        if isinstance(fn, FunctionType):
            try:
                sig = signature(fn)
                if "axis" in sig.parameters.keys():
                    handled_fns[fn] = wrap(fn)
            except ValueError:
                pass

    handled_fns.pop(np.expand_dims, None)

    handled_fns[np.allclose] = wrap(np.allclose)
    handled_fns[np.array_repr] = wrap(np.array_repr)
    handled_fns[entropy] = wrap(entropy)  # qol change

    return handled_fns


_HANDLED_FUNCTIONS = _register_handled_functions()


class LineageMeta(type):
    """
    Metaclass for Lineaage.

    Registers functions which are handled by us and overloads common attibutes, such as `.sum` with these functions.
    """

    __overloaded_functions__ = dict(  # noqa
        sum=np.sum,
        mean=np.mean,
        min=np.min,
        argmin=np.argmin,
        max=np.max,
        argmax=np.argmax,
        std=np.std,
        squeeze=np.squeeze,
        entropy=entropy,
    )

    def __new__(cls, clsname, superclasses, attributedict):  # noqa
        res = type.__new__(cls, clsname, superclasses, attributedict)
        for attrname, fn in LineageMeta.__overloaded_functions__.items():
            wrapped_fn = _HANDLED_FUNCTIONS.get(fn, None)
            if wrapped_fn:
                setattr(res, attrname, wrapped_fn)

        return res


class Lineage(np.ndarray, metaclass=LineageMeta):
    """
    Lightweight :class:`numpy.ndarray` wrapper that adds names and colors.

    Params
    ------
    input_array
        Input array containing lineage probabilities, each lineage being stored in a column.
    names
        Names of the lineages.
    colors
        Colors of the lineages.
    """

    def __new__(
        cls,
        input_array: np.ndarray,
        *,
        names: Iterable[str],
        colors: Optional[Iterable[ColorLike]] = None,
    ) -> "Lineage":
        """Create and return a new object."""

        if not isinstance(input_array, np.ndarray):
            raise TypeError(
                f"Input array must be of type `numpy.ndarray`, found `{type(input_array).__name__!r}`."
            )

        if input_array.ndim == 1:
            input_array = np.expand_dims(input_array, -1)
        elif input_array.ndim > 2:
            raise ValueError(
                f"Input array must be 2-dimensional, found `{input_array.ndim}`."
            )

        if input_array.shape[0] == 0:
            raise ValueError("Expected number cells to be at least 1, found 0.")
        if input_array.shape[1] == 0:
            raise ValueError("Expected number of lineages to be at least 1, found 0.")

        obj = np.array(input_array, copy=True).view(cls)
        obj._n_lineages = obj.shape[1]
        obj._is_transposed = False
        obj.names = names  # these always create a copy, which is a good thing
        obj.colors = colors

        return obj

    def __array_finalize__(self, obj) -> None:
        if obj is None:
            return

        _names = getattr(obj, "_names", None)
        if _names is not None:
            self._names = _names
            self._n_lineages = len(_names)
            self._names_to_ixs = {n: i for i, n in enumerate(self.names)}
        else:
            self._names = None
            self._names_to_ixs = None
            self._n_lineages = getattr(
                obj, "_n_lineages", obj.shape[1] if obj.ndim == 2 else 0
            )

        self._colors = getattr(obj, "colors", None)
        self._is_transposed = getattr(obj, "_is_transposed", False)

    def _mixer(self, rows, mixtures):
        def update_entries(key):
            if key:
                res.append(self[rows, key].X.sum(1))
                # item = (key, rows) if self._is_transposed else (rows, key)
                # res.append(self[item].X.sum(int(not self._is_transposed)))
                names.append(" or ".join(self.names[key]))
                colors.append(_compute_mean_color(self.colors[key]))

        lin_kind = [_ for _ in mixtures if isinstance(_, Lin)]
        if len(lin_kind) > 1:
            raise ValueError(
                f"`Lin` enum is allowed only once in the expression, found `{lin_kind}`."
            )

        keys = [
            tuple(
                self._maybe_convert_names(
                    _convert_lineage_name(mixture), default=mixture
                )
            )
            if isinstance(mixture, str)
            else (mixture,)
            for mixture in mixtures
            if not (isinstance(mixture, Lin))
        ]
        keys = _unique_order_preserving(keys)

        # check the `keys` are unique
        overlap = [set(ks) for ks in keys]
        for c1, c2 in combinations(overlap, 2):
            overlap = c1 & c2
            if overlap:
                raise ValueError(
                    f"Found overlapping keys: `{self.names[list(overlap)]}`."
                )

        seen = set()
        names, colors, res = [], [], []
        for key in map(list, keys):
            seen.update(self.names[key])
            update_entries(key)

        if len(lin_kind) == 1:
            lin_kind = lin_kind[0]
            keys = [i for i, n in enumerate(self.names) if n not in seen]

            if lin_kind == Lin.OTHERS:
                for key in keys:
                    update_entries([key])
            elif lin_kind == Lin.REST:
                update_entries(keys)
                if keys:
                    names[-1] = str(lin_kind)
            else:
                raise ValueError(f"Invalid `Lin` enum `{lin_kind}`.")

        res = np.stack(res, axis=-1)
        return Lineage(res, names=names, colors=colors)

    def __getitem__(self, item) -> "Lineage":
        was_transposed = False
        if self._is_transposed:
            was_transposed = True
            self = self.T
            if isinstance(item, tuple):
                item = item[::-1]

        obj = self.__getitem(item)

        return obj.T if was_transposed else obj

    def __getitem(self, item):
        if isinstance(item, tuple):
            if len(item) > 2:
                raise ValueError(
                    f"Expected key to be of length `2`, found `{len(item)}`."
                )

            item = list(item)
            if item[0] is Ellipsis or item[0] is None:
                item[0] = range(self.shape[0])
            if len(item) == 2 and (item[1] is Ellipsis or item[1] is None):
                item[1] = range(self.shape[1])
            item = tuple(item)

        is_tuple_len_2 = (
            isinstance(item, tuple)
            and len(item) == 2
            and isinstance(
                item[0], (int, np.integer, range, slice, tuple, list, np.ndarray)
            )
        )
        if is_tuple_len_2:
            rows, col = item

            if isinstance(col, (int, np.integer, str)):
                col = [col]

            try:
                # slicing an array where row/col are like 2D indices
                if (
                    len(col) > 1
                    and len(rows) == self.shape[0]
                    and len(rows) == len(col)
                ):
                    col = self._maybe_convert_names(col, make_unique=False)
                    return Lineage(
                        # never remove this expand_dims - it's critical
                        np.expand_dims(self.X[rows, col], axis=-1),
                        names=["mixture"],
                        colors=["#000000"],
                    )
            except TypeError:  # because of range
                pass

            if isinstance(col, (list, tuple, np.ndarray)):
                if any(
                    map(
                        lambda i: isinstance(i, Lin)
                        or (isinstance(i, str) and ("," in i or "or" in i)),
                        col,
                    )
                ):
                    return self._mixer(rows, col)
                col = self._maybe_convert_names(col)
                item = rows, col
        else:
            if isinstance(item, (int, np.integer, str)):
                item = [item]

            col = range(len(self.names))
            if isinstance(item, (tuple, list, np.ndarray)):
                if any(
                    map(
                        lambda i: isinstance(i, Lin)
                        or (isinstance(i, str) and ("," in i or "or" in i)),
                        item,
                    )
                ):
                    return self._mixer(slice(None, None, None), item)
                elif any(map(lambda i: isinstance(i, str), item)):
                    item = (slice(None, None, None), self._maybe_convert_names(item))
                    col = item[1]

        shape, row_order, col_order = None, None, None
        if (
            is_tuple_len_2
            and not isinstance(item[0], slice)
            and not isinstance(item[1], slice)
        ):
            item_0 = (
                np.array(item[0]) if not isinstance(item[0], np.ndarray) else item[0]
            )
            item_1 = (
                np.array(item[1]) if not isinstance(item[1], np.ndarray) else item[1]
            )
            item_0 = _at_least_2d(item_0, -1)
            item_1 = _at_least_2d(item_1, 0)

            if item_1.dtype == np.bool:
                if item_0.dtype != np.bool:
                    if not issubclass(item_0.dtype.type, np.integer):
                        raise TypeError(f"Invalid type `{item_0.dtype.type}`.")
                    row_order = (
                        item_0[:, 0]
                        if item_0.shape[0] == self.shape[0]
                        else np.argsort(item_0[:, 0])
                    )
                    item_0 = _at_least_2d(np.isin(np.arange(self.shape[0]), item_0), -1)
                item = item_0 * item_1
                shape = np.max(np.sum(item, axis=0)), np.max(np.sum(item, axis=1))
                col = np.where(np.all(item_1, axis=0))[0]
            elif item_0.dtype == np.bool:
                if item_1.dtype != np.bool:
                    if not issubclass(item_1.dtype.type, np.integer):
                        raise TypeError(f"Invalid type `{item_1.dtype.type}`.")
                    col_order = (
                        item_1[0, :]
                        if item_1.shape[1] == self.shape[1]
                        else np.argsort(item_1[0, :])
                    )
                    item_1 = _at_least_2d(np.isin(np.arange(self.shape[1]), item_1), 0)

                item = item_0 * item_1
                shape = np.max(np.sum(item, axis=0)), np.max(np.sum(item, axis=1))
                col = np.where(np.all(item_1, axis=0))[0]
            else:
                item = (item_0, item_1)

        obj = super().__getitem__(item)

        if shape is not None:  # keep the resulting shape
            obj = obj.reshape(shape)
        if row_order is not None:
            obj = obj[row_order, :]
        if col_order is not None:
            obj = obj[:, col_order]

        if isinstance(obj, Lineage):
            obj._names = np.atleast_1d(self.names[col])
            obj._colors = np.atleast_1d(self.colors[col])
            if col_order is not None:
                obj._names = obj._names[col_order]
                obj._colors = obj._colors[col_order]
            obj._names_to_ixs = {name: i for i, name in enumerate(obj._names)}

        return obj

    def __array_function__(self, func, types, args, kwargs):
        if func not in _HANDLED_FUNCTIONS:
            return NotImplemented
        # Note: this allows subclasses that don't override
        # __array_function__ to handle MyArray objects
        if not all(issubclass(t, type(self)) for t in types):
            return NotImplemented

        return _HANDLED_FUNCTIONS[func](*args, **kwargs)

    @property
    def names(self) -> np.ndarray:
        """Lineage names. Must be unique."""
        return self._names

    @names.setter
    def names(self, value: Iterable[str]):
        if not isinstance(value, Iterable):
            raise TypeError(_ERROR_NOT_ITERABLE.format("names", type(value).__name__))
        for v in value:
            if not isinstance(v, str):
                raise TypeError(
                    f"Expected `names` to be strings, found type `{type(v).__name__}`."
                )

        value = self._check_shape(value, _ERROR_WRONG_SIZE.format("names"))
        if len(set(value)) != len(value):
            raise ValueError("Not all lineage names are unique.")

        self._names = self._prepare_annotation(value)
        self._names_to_ixs = {name: ix for ix, name in enumerate(self.names)}

    @property
    def colors(self) -> np.ndarray:
        """Lineage colors."""
        return self._colors

    @colors.setter
    def colors(self, value: Optional[Iterable[ColorLike]]):
        if value is None:
            value = _create_categorical_colors(self._n_lineages)
        elif not isinstance(value, Iterable):
            raise TypeError(_ERROR_NOT_ITERABLE.format("colors", type(value).__name__))

        value = self._check_shape(value, _ERROR_WRONG_SIZE.format("colors"))
        self._colors = self._prepare_annotation(
            value,
            checker=c.is_color_like,
            transformer=c.to_hex,
            checker_msg="Value `{}` is not a valid color.",
        )

    def _maybe_convert_names(
        self,
        names: Iterable[Union[int, str, bool]],
        is_singleton: bool = False,
        default: Optional[Union[int, str]] = None,
        make_unique: bool = True,
    ) -> Union[int, List[int], List[bool]]:
        if all(map(lambda n: isinstance(n, (bool, np.bool_)), names)):
            return list(names)
        res = []
        for name in names:
            if isinstance(name, str):
                if name in self._names_to_ixs:
                    name = self._names_to_ixs[name]
                elif default is not None:
                    if isinstance(default, str):
                        if default not in self._names_to_ixs:
                            raise KeyError(
                                f"Invalid lineage name: `{name}`. "
                                f"Valid names are: `{list(self.names)}`."
                            )
                        name = self._names_to_ixs[default]
                    else:
                        name = default
                else:
                    raise KeyError(
                        f"Invalid lineage name `{name}`. "
                        f"Valid names are: `{list(self.names)}`."
                    )
            res.append(name)

        if make_unique:
            res = _unique_order_preserving(res)

        return res[0] if is_singleton else res

    def _check_shape(
        self, array: Iterable[Union[str, ColorLike]], msg: str
    ) -> List[Union[str, ColorLike]]:
        array = list(array)
        if len(array) != self._n_lineages:
            raise ValueError(msg.format(self._n_lineages, len(array)))

        return array

    def _prepare_annotation(
        self,
        array: List[str],
        checker: Optional[Callable] = None,
        transformer: Optional[Callable] = None,
        checker_msg: Optional[str] = None,
    ) -> np.ndarray:
        if checker:
            assert checker_msg, "Please provide a message when `checker` is not `None`."
            for v in array:
                if not checker(v):
                    raise ValueError(checker_msg.format(v))

        if transformer is not None:
            array = np.array([transformer(v) for v in array])

        return np.array(array)

    @property
    def X(self) -> np.ndarray:
        """Convert self to numpy array, losing names and colors."""
        return np.array(self, copy=False)

    @property
    def T(self):
        """Transpose of self."""
        obj = self.transpose()
        obj._is_transposed = not self._is_transposed
        return obj

    def view(self, dtype=None, type=None) -> "LineageView":
        """Return a view of self."""
        return LineageView(self)

    def __repr__(self) -> str:
        return f'{super().__repr__()[:-1]},\n  names([{", ".join(self.names)}]))'

    def __str__(self):
        return f'{super().__str__()}\n names=[{", ".join(self.names)}]'

    def _repr_html_(self) -> str:
        def format_row(r):
            rng = (
                range(self.shape[1])
                if not self._is_transposed
                or (self._is_transposed and self.shape[1] <= _HTML_REPR_THRESH)
                else list(range(_HT_CELLS))
                + [...]
                + list(range(self.shape[1] - _HT_CELLS - 1, self.shape[1] - 1))
            )

            cells = "".join(
                f"<td style='text-align: right;'>" f"{self.X[r, c]:.06f}" f"</td>"
                if isinstance(c, int)
                else _DUMMY_CELL
                for c in rng
            )
            return f"<tr>{(names[r] if self._is_transposed else '') + cells}</tr>"

        def dummy_row():
            values = "".join(_DUMMY_CELL for _ in range(self.shape[1]))
            return f"<tr>{values}</tr>"

        if self.names is None or self.colors is None:
            raise RuntimeError(
                f"Name or colors are `None`. This can happen when running `array.view({type(self).__name__}`."
            )

        if self.ndim != 2:
            return repr(self)

        styles = [
            f"'background-color: {bg}; color: {fg}; text-align: center; word-wrap: break-word; max-width: 100px'"
            for bg, fg in map(_get_bg_fg_colors, self.colors)
        ]
        names = [f"<th style={style}>{n}</th>" for n, style in zip(self.names, styles)]
        header = f"<tr>{''.join(names)}</tr>"

        if self.shape[0] > _HTML_REPR_THRESH:
            body = "".join(format_row(i) for i in range(_HT_CELLS))
            body += dummy_row()
            body += "".join(
                format_row(i)
                for i in range(self.shape[0] - _HT_CELLS - 1, self.shape[0] - 1)
            )
        else:
            body = "".join(format_row(i) for i in range(self.shape[0]))

        cells = "cells" if self.shape[0] > 1 else "cell"
        lineages = "lineages" if self.shape[1] > 1 else "lineage"
        if self._is_transposed:
            cells, lineages = lineages, cells
        metadata = f"<p>{self.shape[0]} {cells} x {self.shape[1]} {lineages}</p>"

        if self._is_transposed:
            header = ""
        return f"<div style='scoped'><table>{header}{body}</table>{metadata}</div>"

    def __format__(self, format_spec):
        if self.shape == (1, 1):
            return format_spec.format(self.X[0, 0])
        if self.shape == (1,):
            return format_spec.format(self.X[0])
        return NotImplemented

    def __setstate__(self, state):
        *state, names, colors, is_t = state
        names = names[-1]
        colors = colors[-1]

        self._names = np.empty(names[1])
        self._colors = np.empty(colors[1])

        super().__setstate__(tuple(state))
        self._names.__setstate__(tuple(names))
        self._colors.__setstate__(tuple(colors))

        self._is_transposed = is_t
        self._n_lineages = len(self.names)
        self._names_to_ixs = {name: ix for ix, name in enumerate(self.names)}

    def __reduce__(self):
        res = list(super().__reduce__())

        names = self.names.__reduce__()
        colors = self.colors.__reduce__()
        res[-1] += (names, colors, self._is_transposed)

        return tuple(res)

    def copy(self, _="C") -> "Lineage":
        """Return a copy of itself."""
        obj = Lineage(
            self.T if self._is_transposed else self,
            names=np.array(self.names, copy=True, order=_ORDER),
            colors=np.array(self.colors, copy=True, order=_ORDER),
        )
        return obj.T if self._is_transposed else obj

    def __copy__(self):
        return self.copy()

    def plot_pie(
        self,
        reduction: Callable = np.mean,
        title: Optional[str] = None,
        figsize: Optional[Tuple[float, float]] = None,
        dpi: Optional[float] = None,
        save: Optional[Union[Path, str]] = None,
        **kwargs,
    ) -> None:
        """
        Plot a pie chart visualizing the aggregated lineage probabilities.

        Params
        ------
        reduction
            Function that will be applied per lineage.
        title
            Title of the figure.
        figsize
            Size of the figure.
        dpi
            Dots per inch.
        save
            Filename where to save the plots.
            If `None`, just shows the plot.

        Returns
        -------
        None
            Nothing, just plots the pie chart.
        """
        if not callable(reduction):
            raise TypeError(
                f"Expected `reduction` to be callable, found `{type(reduction).__name__}`."
            )

        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        title = reduction.__name__ if title is None else title

        reduction = reduction(self, axis=int(self._is_transposed))
        ax.pie(reduction.X.squeeze(), labels=self.names, colors=self.colors, **kwargs)
        ax.set_title(title)

        fig.show()

        if save is not None:
            save_fig(fig, save)

    def reduce(
        self,
        keys: Union[str, List[str], Tuple[str], np.ndarray],
        mode: str = "dist",
        dist_measure: str = "mutual_info",
        normalize_weights: str = "softmax",
        softmax_beta: float = 1,
        return_weights: bool = False,
    ) -> Union["Lineage", Tuple["Lineage", Optional[pd.DataFrame]]]:
        """
        Reduce metastable states to final/root states.

        Params
        ------
        keys
            List of keys that define the final/root states. The lineage will be reduced
            to these states by projecting the other states.
        mode
            Whether to use a distance measure to compute weights ('dist', or just rescale ('scale').
            Scaling is baseline for benchmarking.
        dist_measure
            Used to quantify similarity between query and reference states. Valid options are:

            - 'cosine_sim'
            - 'wasserstein_dist'
            - 'kl_div'
            - 'js_div'
            - 'mutual_inf'
            - 'equal'
        normalize_weights
            How to normalize the weights. Valid options are:

            - 'scale': divide by the sum (per row)
            - 'softmax': use a softmax with beta = 1
        softmax_beta
            Scaling factor in the softmax, used for normalizing the weights to sum to 1.
        return_weights
            If `True`, a :class:`pandas.DataFrame` of the weights used for the projection is returned.

        Returns
        -------
        :class:`cellrank.tl.Lineage`
            Lineage object, reduced to the final/root states.
        :class:`pandas.DataFrame`
            The weights used for the projection of shape `(n_query x n_reference)`.
        """
        if self._is_transposed:
            raise RuntimeError("This matrix seems to be transposed.")

        if isinstance(keys, str):
            keys = [keys]

        if not len(keys):
            raise ValueError("Unable to perform the reduction, no keys specified.")

        if set(keys) == set(self.names):
            logg.warning(
                "Unable to perform the reduction, `keys` specifies all lineages. Returning self"
            )
            return (self, None) if return_weights else self

        # check input parameters
        if return_weights and mode == "scale":
            logg.warning(
                "If `mode=='scale'`, no weights are computed. Returning `None`"
            )

        # check the lineage object
        if not np.allclose(np.sum(self, axis=1).X, 1.0):
            raise ValueError("Memberships do not sum to one row-wise.")

        # check the keys are all in L.names
        key_mask = np.array([key in self.names for key in keys])
        if not key_mask.all():
            raise ValueError(
                f"Invalid lineage names `{list(np.array(keys)[~key_mask])}`. Valid names are: `{list(self.names)}`."
            )

        # get query and reference
        mask = np.in1d(self.names, keys)
        reference = self[:, mask]
        query = self[:, ~mask]

        if mode == "dist":
            # compute a set of weights of shape (n_query x n_reference)
            if dist_measure == "cosine_sim":
                weights = _cosine_sim(reference.X, query.X)
            elif dist_measure == "wasserstein_dist":
                weights = _wasserstein_dist(reference.X, query.X)
            elif dist_measure == "kl_div":
                weights = _kl_div(reference.X, query.X)
            elif dist_measure == "js_div":
                weights = _js_div(reference.X, query.X)
            elif dist_measure == "mutual_info":
                weights = _mutual_info(reference.X, query.X)
            elif dist_measure == "equal":
                weights = _row_normalize(np.ones((query.shape[1], reference.shape[1])))
            else:
                raise ValueError(f"Invalid distance measure `{dist_measure!r}`.")

            # make some checks on the weights
            if weights.shape != (query.shape[1], reference.shape[1]):
                raise ValueError(
                    f"Expected weight matrix to be of shape `({query.shape[1]}, {reference.shape[1]})`, "
                    f"found `{weights.shape}`."
                )
            if not np.isfinite(weights).all():
                raise ValueError(
                    "Weights matrix contains elements that are not finite."
                )
            if (weights < 0).any():
                raise ValueError("Weights matrix contains negative elements.")

            if (weights == 0).any():
                logg.warning("Weights matrix contains exact zeros.")

            # normalize the weights to row-sum to one
            if normalize_weights == "scale":
                weights_n = _row_normalize(weights)
            elif normalize_weights == "softmax":
                weights_n = _softmax(_row_normalize(weights), softmax_beta)
            else:
                raise ValueError(
                    f"Normalization method `{normalize_weights!r}` not found. Valid options are: `'scale', 'softmax'`."
                )

            # check that the weights row-sum to one now
            if not np.allclose(weights_n.sum(1), 1.0):
                raise ValueError("Weights do not sum to 1 row-wise.")

            # use the weights to re-distribute probability mass form query to reference
            for i, w in enumerate(weights_n):
                reference += np.dot(query[:, i].X, w[None, :])

        elif mode == "scale":
            reference = _row_normalize(reference)
        else:
            raise ValueError(
                f"Invalid mode `{mode!r}`. Valid options are: `'dist', 'scale'`."
            )

        # check that the lineages row-sum to one now
        if not np.allclose(reference.sum(1), 1.0):
            raise ValueError("Reduced lineage rows do not sum to 1.")

        # potentially create a weights-df and return everything
        if return_weights:
            if mode == "dist":
                return (
                    reference,
                    pd.DataFrame(
                        data=weights_n, columns=reference.names, index=query.names
                    ),
                )
            return reference, None

        return reference


class LineageView(Lineage):
    """Simple view of :class:`cellrank.tools.Lineage`."""

    def __new__(cls, lineage: Lineage) -> "LineageView":
        """Create a LineageView."""
        if not isinstance(lineage, Lineage):
            raise TypeError(
                f"Cannot create a `{cls.__name__}` of `{type(lineage).__name__}`."
            )

        view = np.array(lineage, copy=False).view(cls)
        view._owner = lineage
        view._names = lineage.names
        view._n_lineages = len(view.names)
        view._names_to_ixs = lineage._names_to_ixs
        view._colors = lineage.colors
        view._is_transposed = lineage._is_transposed

        return view

    @property
    def names(self) -> np.ndarray:
        """Lineage names."""
        return super().names

    @property
    def owner(self) -> Lineage:
        """Return the lineage associated with this view."""
        return self._owner

    @names.setter
    def names(self, _):
        raise RuntimeError(f"Unable to set names of `{type(self).__name__}`.")

    @property
    def colors(self) -> np.ndarray:
        """Lineage colors."""
        return super().colors

    @colors.setter
    def colors(self, _):
        raise RuntimeError(f"Unable to set colors of `{type(self).__name__}`.")

    def view(self, dtype=None, type=None) -> "LineageView":
        """Return self."""
        return self

    def copy(self, _="C") -> Lineage:
        """Return a copy of itself."""
        was_trasposed = False
        if self._is_transposed:
            self = self.T
            was_trasposed = True

        obj = Lineage(
            self,
            names=np.array(self.names, copy=True, order=_ORDER),
            colors=np.array(self.colors, copy=True, order=_ORDER),
        )
        return obj.T if was_trasposed else obj


def _remove_zero_rows(a: np.ndarray, b: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    if a.shape[0] != b.shape[0]:
        raise ValueError("Lineage objects have unequal cell numbers")

    bool_a = (a == 0).any(axis=1)
    bool_b = (b == 0).any(axis=1)
    mask = ~np.logical_or(bool_a, bool_b)

    logg.warning(
        f"Removed {a.shape[0] - np.sum(mask)} rows because they contained zeros"
    )

    return a[mask, :], b[mask, :]


def _softmax(X, beta: float = 1):
    return np.exp(X * beta) / np.expand_dims(np.sum(np.exp(X * beta), axis=1), -1)


def _row_normalize(X):
    if isinstance(X, Lineage):
        return X / X.sum(1)  # Lineage is shape-preserving
    return X / np.expand_dims(X.sum(1), -1)


def _col_normalize(X, norm_ord=2):
    from numpy.linalg import norm

    return X / norm(X, ord=norm_ord, axis=0)


def _cosine_sim(reference, query):
    # the cosine similarity is symmetric

    # normalize these to have 2-norm 1
    reference_n, query_n = _col_normalize(reference, 2), _col_normalize(query, 2)

    return (reference_n.T @ query_n).T


def _point_wise_distance(reference, query, distance):
    # utility function for all point-wise distances/divergences
    # take care of rows that contain zeros
    reference_no_zero, query_no_zero = _remove_zero_rows(reference, query)

    # normalize these to be valid probability distributions (column-wise)
    reference_n, query_n = (
        _col_normalize(reference_no_zero, 1),
        _col_normalize(query_no_zero, 1),
    )

    # loop over query and reference columns and compute pairwise wasserstein distances
    weights = np.zeros((query.shape[1], reference.shape[1]))
    for i, q_d in enumerate(query_n.T):
        for j, r_d in enumerate(reference_n.T):
            weights[i, j] = 1.0 / distance(q_d, r_d)

    return weights


def _wasserstein_dist(reference, query):
    # the wasserstein distance is symmetric
    from scipy.stats import wasserstein_distance

    return _point_wise_distance(reference, query, wasserstein_distance)


def _kl_div(reference, query):
    # the KL divergence is not symmetric
    from scipy.stats import entropy

    return _point_wise_distance(reference, query, entropy)


def _js_div(reference, query):
    # the js divergence is symmetric
    from scipy.spatial.distance import jensenshannon

    return _point_wise_distance(reference, query, jensenshannon)


def _mutual_info(reference, query):
    # mutual information is not symmetric. We don't need to normalise the vectors, it's invariant under scaling.
    from sklearn.feature_selection import mutual_info_regression

    weights = np.zeros((query.shape[1], reference.shape[1]))
    for i, target in enumerate(query.T):
        weights[i, :] = mutual_info_regression(reference, target)

    return weights
