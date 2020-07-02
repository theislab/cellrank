# -*- coding: utf-8 -*-
"""Lineage class module."""

from typing import List, Tuple, Union, TypeVar, Callable, Iterable, Optional
from itertools import combinations

import matplotlib.colors as c

from scanpy import logging as logg

import numpy as np
import pandas as pd
from cellrank.tools._utils import (
    _compute_mean_color,
    _convert_lineage_name,
    _unique_order_preserving,
)
from cellrank.tools._colors import _get_bg_fg_colors, _create_categorical_colors
from cellrank.tools._constants import Lin

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

    def view(self, dtype=None, type=None) -> "LineageView":
        """Return a view of self."""
        return LineageView(self, names=self.names, colors=self.colors)

    def __array_finalize__(self, obj) -> None:
        if obj is None:
            return

        _names = getattr(obj, "_names", None)
        if _names is not None:
            self._n_lineages = len(_names)
            self.names = _names
        else:
            self._names = None
            self._names_to_ixs = None

        self._colors = getattr(obj, "colors", None)
        self._n_lineages = getattr(obj, "_n_lineages", 0)

    def _mixer(self, rows, mixtures):
        def update_entries(key):
            if key:
                res.append(self[rows, key].X.sum(1))
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
                        self.X[rows, col], names=["mixture"], colors=["#000000"]
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
                            # TODO: names_to_ixs is not being properly updated
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
        return np.array(self)

    def __repr__(self) -> str:
        return f'{super().__repr__()[:-1]},\n  names([{", ".join(self.names)}]))'

    def __str__(self):
        return f'{super().__str__()}\n names=[{", ".join(self.names)}])'

    def _repr_html_(self) -> str:
        def format_row(r):
            cells = "".join(
                f"<td style='text-align: right;'>"
                f"{super(type(self), self).__getitem__((r, c)):.06f}"
                f"</td>"
                for c in range(self.shape[1])
            )
            return f"<tr>{cells}</tr>"

        def dummy_row():
            values = "".join(
                "<td style='text-align: right;'>...</td>" for _ in range(self.shape[1])
            )
            return f"<tr>{values}</tr>"

        show_n_cells_head_tail = 10
        styles = [
            f"'background-color: {bg}; color: {fg}; text-align: center; word-wrap: break-word; max-width: 100px'"
            for bg, fg in map(_get_bg_fg_colors, self.colors)
        ]
        names = "".join(
            f"<th style={style}>{n}</th:w:w>" for n, style in zip(self.names, styles)
        )
        header = f"<tr>{names}</tr>"

        if self.shape[0] > 100:
            body = "".join(format_row(i) for i in range(show_n_cells_head_tail))
            body += dummy_row()
            body += "".join(
                format_row(i)
                for i in range(
                    self.shape[0] - show_n_cells_head_tail - 1, self.shape[0] - 1
                )
            )
        else:
            body = "".join(format_row(i) for i in range(self.shape[0]))

        cells = "cells" if self.shape[0] > 1 else "cell"
        lineages = "lineages" if self.shape[1] > 1 else "lineage"
        metadata = f"<p>{self.shape[0]} {cells} x {self.shape[1]} {lineages}</p>"

        return f"<div style='scoped'><table>{header}{body}</table>{metadata}</div>"

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
        """Return a copy of itself."""
        return Lineage(
            np.array(self, copy=True, order=order),
            names=np.array(self.names, copy=True, order=order),
            colors=np.array(self.colors, copy=True, order=order),
        )

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
        if not np.allclose(np.sum(self, axis=1), 1.0):
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
                weights = _cosine_sim(reference.copy(), query.copy())
            elif dist_measure == "wasserstein_dist":
                weights = _wasserstein_dist(reference.copy(), query.copy())
            elif dist_measure == "kl_div":
                weights = _kl_div(reference.copy(), query.copy())
            elif dist_measure == "js_div":
                weights = _js_div(reference.copy(), query.copy())
            elif dist_measure == "mutual_info":
                weights = _mutual_info(reference.copy(), query.copy())
            elif dist_measure == "equal":
                weights = np.ones((query.shape[1], reference.shape[1]))
                weights = _row_normalize(weights)
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

    def copy(self, order="C") -> Lineage:
        """Return a copy of itself."""
        return Lineage(
            np.array(self, copy=True, order=order),
            names=np.array(self.names, copy=True, order=order),
            colors=np.array(self.colors, copy=True, order=order),
        )


def _remove_zero_rows(a: Lineage, b: Lineage) -> Tuple[Lineage, Lineage]:
    if a.shape[0] != b.shape[0]:
        raise ValueError("Lineage objects have unequal cell numbers")

    bool_a = (a.X == 0).any(axis=1)
    bool_b = (b.X == 0).any(axis=1)
    mask = ~np.logical_or(bool_a, bool_b)

    logg.warning(
        f"Removed {a.shape[0] - np.sum(mask)} rows because they contained zeros"
    )

    return a[mask, :], b[mask, :]


def _softmax(X, beta: float = 1):
    return np.exp(X * beta) / np.sum(np.exp(X * beta), axis=1)[:, None]


def _row_normalize(X):
    return X / np.expand_dims(X.sum(1), -1)


def _col_normalize(X, norm_ord=2):
    from numpy.linalg import norm

    return X / norm(X, ord=norm_ord, axis=0)


def _cosine_sim(reference, query):
    # the cosine similarity is symmetric

    # normalize these to have 2-norm 1
    reference_n, query_n = _col_normalize(reference, 2), _col_normalize(query, 2)

    return (reference_n.X.T @ query_n.X).T


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
    for i, q_d in enumerate(query_n.T.X):
        for j, r_d in enumerate(reference_n.T.X):
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
    for i, target in enumerate(query.X.T):
        weights[i, :] = mutual_info_regression(reference.X, target)

    return weights
