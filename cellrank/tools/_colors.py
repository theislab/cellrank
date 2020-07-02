# -*- coding: utf-8 -*-
"""Color handling module."""

from typing import Any, List, Tuple, Union, Optional, Sequence

from matplotlib import cm as cm
from matplotlib import colors as mcolors

from scanpy import logging as logg

import numpy as np
from pandas import Series, DataFrame, to_numeric
from scipy.stats import entropy
from pandas._libs.lib import infer_dtype
from pandas.core.dtypes.common import is_categorical_dtype


def _create_colors(
    base_color: Union[str, Tuple[float, float, float]],
    n: int,
    hue_range: Optional[Tuple[float, float]] = (-0.1, 0.1),
    saturation_range: Optional[Tuple[float, float]] = (-0.3, 0.3),
    value_range: Optional[Tuple[float, float]] = (-0.3, 0.3),
    convert_to_rgb: bool = True,
    as_hex: bool = True,
) -> List[Any]:
    """
    Create variations of colors from base color.

    Params
    ------
    base_color
        Base color which serves as a starting point.
    n
        Number of colors to create.
    hue_range
        Minimum and maximum value to add to the base color's hue.
        If `None`, don't adjust the hue.
    saturation_range
        Minimum and maximum value to add to the base color's saturation.
        If `None`, don't adjust the saturation.
    value_range
        Minimum and maximum value to add to the base color's value.
        If `None`, don't adjust the value.
    convert_to_rgb
        Whether to convert colors from HSV to RGB.
    as_hex:
        Whether to return colors as hex string.

    Returns
    -------
    :class:`list`
        List of colors, either as a hex string or an RGB array.
    """

    if not mcolors.is_color_like(base_color):
        raise ValueError("Base color is not color-like.")
    if n <= 0:
        raise ValueError(f"Number of colors must be > 0, found `{n}`.")

    base_color = mcolors.rgb_to_hsv(mcolors.to_rgb(base_color))

    if n == 1:
        colors = [base_color]
    else:
        n *= 2  # sometimes the colors are too similar, we take every 2nd one
        colors = np.repeat(base_color[..., np.newaxis], n, axis=1).T

        for i, r in enumerate((hue_range, saturation_range, value_range)):
            if r is None:
                continue
            r_low, r_high = sorted(r)
            c = base_color[i]

            colors[:, i] = np.linspace(max(c + r_low, 0), min(c + r_high, 1), n)

    if convert_to_rgb:
        colors = map(mcolors.hsv_to_rgb, colors)
    if as_hex:
        colors = map(mcolors.to_hex, colors)

    return list(colors)[::2]  # we've created twice as much colors, select every other


def _convert_to_hex_colors(colors: Sequence[Any]) -> List[str]:
    if not all(mcolors.is_color_like(c) for c in colors):
        raise ValueError("Not all colors are color-like.")

    return [mcolors.to_hex(c) for c in colors]


def _create_categorical_colors(n_categories: Optional[int] = None):
    if n_categories is None:
        n_categories = 51
    if n_categories > 51:
        raise ValueError(f"Maximum number of colors (51) exceeded: `{n_categories}`.")

    colors = [cm.Set1(i) for i in range(cm.Set1.N)][:n_categories]
    colors += [cm.Set2(i) for i in range(cm.Set2.N)][: n_categories - len(colors)]
    colors += [cm.Set3(i) for i in range(cm.Set3.N)][: n_categories - len(colors)]
    colors += [cm.tab10(i) for i in range(cm.tab10.N)][: n_categories - len(colors)]
    colors += [cm.Paired(i) for i in range(cm.Paired.N)][: n_categories - len(colors)]

    return _convert_to_hex_colors(colors)


def _insert_categorical_colors(seen_colors: Union[np.ndarray, List], n_categories: int):
    seen_colors = set(_convert_to_hex_colors(seen_colors))
    candidates = list(
        filter(lambda c: c not in seen_colors, _create_categorical_colors())
    )[:n_categories]

    if len(candidates) != n_categories:
        raise RuntimeError(f"Unable to create `{n_categories}` categorical colors.")

    return candidates


def _contrasting_color(r: int, g: int, b: int) -> str:
    for val in [r, g, b]:
        assert 0 <= val <= 255

    return "#000000" if r * 0.299 + g * 0.587 + b * 0.114 > 186 else "#ffffff"


def _get_black_or_white(value: float, cmap) -> str:
    if not (0.0 <= value <= 1.0):
        raise ValueError(f"Value must be in range `[0, 1]`, found `{value}`.")

    r, g, b, *_ = [int(c * 255) for c in cmap(value)]
    return _contrasting_color(r, g, b)


def _get_bg_fg_colors(color, sat_scale: Optional[float] = None) -> Tuple[str, str]:
    if not mcolors.is_color_like(color):
        raise ValueError(f"Value `{color}` is not color-like.")

    color = np.squeeze(mcolors.to_rgba_array(color, alpha=1))[:3]
    if sat_scale is not None:
        h, s, v = mcolors.rgb_to_hsv(color)
        color = mcolors.hsv_to_rgb([h, s * sat_scale, v])

    return (
        mcolors.to_hex(color),
        _contrasting_color(*np.array(color * 255).astype(np.int)),
    )


def _map_names_and_colors(
    series_reference: Series,
    series_query: Series,
    colors_reference: Optional[np.array] = None,
    en_cutoff: Optional[float] = None,
) -> Union[Series, Tuple[Series, List[Any]]]:
    """
    Map annotations and colors from one series to another.

    Params
    ------
    series_reference
        Series object with categorical annotations.
    series_query
        Series for which we would like to query the category names.
    colors_reference
        If given, colors for the query categories are pulled from this color array.
    en_cutoff
        In case of a non-perfect overlap between categories of the two series,
        this decides when to label a category in the query as 'Unknown'.

    Returns
    -------
    :class:`pandas.Series`, :class:`list`
        Series with updated category names and a corresponding array of colors.
    """

    # checks: dtypes, matching indices, make sure colors match the categories
    if not is_categorical_dtype(series_reference):
        raise TypeError(
            f"Reference series must be `categorical`, found `{infer_dtype(series_reference)}`."
        )
    if not is_categorical_dtype(series_query):
        raise TypeError(
            f"Query series must be `categorical`, found `{infer_dtype(series_query)}`."
        )
    index_query, index_reference = series_query.index, series_reference.index
    if not np.all(index_reference == index_query):
        raise ValueError("Series indices do not match, cannot map names/colors.")

    process_colors = colors_reference is not None
    if process_colors:
        if len(series_reference.cat.categories) != len(colors_reference):
            raise ValueError(
                "Length of reference colors does not match length of reference series."
            )
        if not all(mcolors.is_color_like(c) for c in colors_reference):
            raise ValueError("Not all colors are color-like.")

    # create dataframe to store the associations between reference and query
    cats_query = series_query.cat.categories
    cats_reference = series_reference.cat.categories
    association_df = DataFrame(None, index=cats_query, columns=cats_reference)

    # populate the dataframe - compute the overlap
    for cl in cats_query:
        row = [
            np.sum(series_reference.loc[np.array(series_query == cl)] == key)
            for key in cats_reference
        ]
        association_df.loc[cl] = row
    association_df = association_df.apply(to_numeric)

    # find the mapping which maximizes overlap and compute entropy
    names_query = association_df.T.idxmax()
    association_df["entropy"] = entropy(association_df.T)
    association_df["name"] = names_query

    # assign query colors
    if process_colors:
        colors_query = []
        for name in names_query:
            mask = cats_reference == name
            color = np.array(colors_reference)[mask][0]
            colors_query.append(color)
        association_df["color"] = colors_query

    # next, we need to make sure that we have unique names and colors. In a first step, compute how many repetitions
    # we have
    names_query_series = Series(names_query, dtype="category")
    frequ = {
        key: np.sum(names_query == key) for key in names_query_series.cat.categories
    }

    names_query_new = np.array(names_query.copy())
    if process_colors:
        colors_query_new = np.array(colors_query.copy())

    # Create unique names by adding suffixes "..._1, ..._2" etc and unique colors by shifting the original color
    for key, value in frequ.items():
        if value == 1:
            continue  # already unique, skip

        # deal with non-unique names
        suffix = list(np.arange(1, value + 1).astype("str"))
        unique_names = [f"{key}_{rep}" for rep in suffix]
        names_query_new[names_query_series == key] = unique_names
        if process_colors:
            color = association_df[association_df["name"] == key]["color"].values[0]
            shifted_colors = _create_colors(color, value, saturation_range=None)
            colors_query_new[np.array(colors_query) == color] = shifted_colors

    association_df["name"] = names_query_new
    if process_colors:
        association_df["color"] = _convert_to_hex_colors(
            colors_query_new
        )  # original colors can be still there, convert to hex

    # issue a warning for mapping with high entropy
    if en_cutoff is not None:
        critical_cats = list(
            association_df.loc[association_df["entropy"] > en_cutoff, "name"].values
        )
        if len(critical_cats) > 0:
            logg.warning(
                f"The following states could not be mapped uniquely: `{', '.join(map(str, critical_cats))}`"
            )

    return (
        (association_df["name"], list(association_df["color"]))
        if process_colors
        else association_df["name"]
    )
