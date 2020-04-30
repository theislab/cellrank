# -*- coding: utf-8 -*-
from cellrank.tools._utils import (
    _create_categorical_colors,
    unique_order_preserving,
    _compute_mean_color,
)

from typing import Optional, Iterable, Callable, TypeVar, List, Union
from itertools import combinations

import matplotlib.colors as c
import numpy as np


ColorLike = TypeVar("ColorLike")
_ERROR_NOT_ITERABLE = "Expected `{}` to be iterable, found type `{}`."
_ERROR_WRONG_SIZE = "Expected `{}` to be of size `{{}}`, found `{{}}`."


def _at_least_2d(array: np.ndarray, dim: int):
    return np.expand_dims(array, dim) if array.ndim < 2 else array


class Lineage(np.ndarray):
    """
    Lightweight :class:`numpy.ndarray` wrapper that adds names and colors.

    Params
    ------
    input_array
        Input array containing lineage probabilities as columns.
    names
        Names of the lineages.
    colors
        Colors of the lineages
    """

    def __new__(
        cls,
        input_array: np.ndarray,
        *,
        names: Iterable[str],
        colors: Optional[Iterable[ColorLike]] = None,
    ) -> "Lineage":
        if not isinstance(input_array, np.ndarray):
            raise TypeError(
                f"Input array must be of type `numpy.ndarray`, found `{type(input_array).__name__!r}`"
            )

        obj = np.asarray(input_array).view(cls)
        if obj.ndim == 1:
            obj = np.expand_dims(obj, -1)
        elif obj.ndim > 2:
            raise ValueError(
                f"Input array must be 2-dimensional, found `{obj.ndim}` dimensions."
            )

        obj._n_lineages = obj.shape[1]
        obj.names = names
        obj.colors = colors

        return obj

    def view(self, dtype=None, type=None):
        return LineageView(self, names=self.names, colors=self.colors)

    def __array_finalize__(self, obj) -> None:
        if obj is None:
            return

        self._names = getattr(obj, "names", None)
        self._colors = getattr(obj, "colors", None)
        self._n_lineages = getattr(obj, "_n_lineages", 0)
        self._names_to_ixs = getattr(
            "obj",
            "_names_to_ixs",
            None
            if self.names is None
            else {name: ix for ix, name in enumerate(self.names)},
        )

    def _mixer(self, rows, mixtures):
        def update_entries(key):
            if key:
                res.append(self[rows, key].X.sum(1))
                names.append(" or ".join(self.names[key]))
                colors.append(_compute_mean_color(self.colors[key]))

        n_ellipses_or_none = len([_ for _ in mixtures if _ is Ellipsis or _ is None])
        if n_ellipses_or_none > 1:
            raise ValueError(f"`None` or `...` is allowed only once in the expression.")

        keys = [
            tuple(
                self._maybe_convert_names(
                    sorted({key.strip(" ") for key in mixture.strip(" ,").split(",")})
                )
            )
            if isinstance(mixture, str)
            else (mixture,)
            for mixture in mixtures
            if mixture is not Ellipsis and mixture is not None
        ]
        keys = unique_order_preserving(keys)

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

        if n_ellipses_or_none == 1:
            update_entries([i for i, n in enumerate(self.names) if n not in seen])
            if None in mixtures:
                names[-1] = "rest"

        res = np.stack(res, axis=-1)

        return Lineage(res, names=names, colors=colors)

    def __getitem__(self, item) -> "Lineage":
        is_tuple_len_2 = (
            isinstance(item, tuple)
            and len(item) == 2
            and isinstance(item[0], (int, range, slice, tuple, list, np.ndarray))
        )
        if is_tuple_len_2:
            rows, col = item

            if isinstance(col, (int, str)):
                col = [col]

            if isinstance(col, (list, tuple)):
                if any(
                    map(
                        lambda i: i is Ellipsis
                        or i is None
                        or (isinstance(i, str) and "," in i),
                        col,
                    )
                ):
                    return self._mixer(rows, col)
                col = self._maybe_convert_names(col)
                item = rows, col
        else:
            if isinstance(item, (int, str)):
                item = [item]

            col = range(len(self.names))
            if isinstance(item, (tuple, list)):
                if any(
                    map(
                        lambda i: i is Ellipsis
                        or i is None
                        or (isinstance(i, str) and "," in i),
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
                    row_order = np.argsort(item_0[:, 0])
                    item_0 = _at_least_2d(np.isin(np.arange(self.shape[0]), item_0), -1)
                item = item_0 * item_1
                shape = np.max(np.sum(item, axis=0)), np.max(np.sum(item, axis=1))
                col = np.where(np.all(item_1, axis=0))[0]
            elif item_0.dtype == np.bool:
                if item_1.dtype != np.bool:
                    if not issubclass(item_1.dtype.type, np.integer):
                        raise TypeError(f"Invalid type `{item_1.dtype.type}`.")
                    col_order = np.argsort(item_1[0, :])
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

        return obj

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
        self, names: Iterable[Union[int, str]], is_singleton: bool = False
    ) -> Union[int, List[int]]:
        res = []
        for name in names:
            if isinstance(name, str):
                if name not in self._names_to_ixs:
                    raise KeyError(
                        f"Invalid lineage name `{name}`. "
                        f"Valid names are: `{list(self._names_to_ixs.keys())}`."
                    )
                name = self._names_to_ixs[name]
            res.append(name)

        res = unique_order_preserving(res)

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
        return np.array(self)

    def __repr__(self) -> str:
        return f'{super().__repr__()[:-1]},\n  names([{", ".join(self.names)}]))'

    def __str__(self):
        return f'{super().__str__()}\n names=[{", ".join(self.names)}])'

    def __setstate__(self, state):
        *state, names, colors = state
        names = names[-1]
        colors = colors[-1]
        self._names = np.empty(names[1])
        self._colors = np.empty(colors[1])

        super().__setstate__(tuple(state))
        self._names.__setstate__(tuple(names))
        self._colors.__setstate__(tuple(colors))

        self._names_to_ixs = {name: ix for ix, name in enumerate(self.names)}

    def __reduce__(self):
        res = list(super().__reduce__())

        names = self.names.__reduce__()
        colors = self.colors.__reduce__()
        res[-1] += (names, colors)

        return tuple(res)

    def copy(self, order="C") -> "Lineage":
        return Lineage(
            np.array(self, copy=True, order=order),
            names=np.array(self.names, copy=True, order=order),
            colors=np.array(self.colors, copy=True, order=order),
        )


class LineageView(Lineage):
    def copy(self, order="C") -> Lineage:
        return Lineage(
            np.array(self, copy=True, order=order),
            names=np.array(self.names, copy=True, order=order),
            colors=np.array(self.colors, copy=True, order=order),
        )
