# -*- coding: utf-8 -*-
"""Color handling module."""

from typing import Any, List, Tuple, Union, Iterable, Optional, Sequence

import numpy as np
from pandas import Series, DataFrame, to_numeric
from scipy.stats import entropy
from pandas._libs.lib import infer_dtype
from pandas.core.dtypes.common import is_categorical_dtype

from matplotlib import cm as cm
from matplotlib import colors as mcolors

from cellrank import logging as logg


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

    Parameters
    ----------
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
        raise ValueError("Not all values are color-like.")

    return [mcolors.to_hex(c) for c in colors]


def _create_categorical_colors(n_categories: Optional[int] = None):
    from scanpy.plotting.palettes import vega_20_scanpy

    cmaps = [
        mcolors.ListedColormap(vega_20_scanpy),
        cm.Accent,
        mcolors.ListedColormap(np.array(cm.Dark2.colors)[[1, 2, 4, 5, 6]]),
        cm.Set1,
        cm.Set2,
        cm.Set3,
    ]
    max_cats = sum(c.N for c in cmaps)

    if n_categories is None:
        n_categories = max_cats
    if n_categories > max_cats:
        raise ValueError(
            f"Maximum number of colors ({max_cats}) exceeded: `{n_categories}`."
        )

    colors = []
    for cmap in cmaps:
        colors += [cmap(i) for i in range(cmap.N)][: n_categories - len(colors)]
        if len(colors) == n_categories:
            return _convert_to_hex_colors(colors)

    raise RuntimeError(f"Unable to create `{n_categories}` colors.")


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

    Parameters
    ----------
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
    if len(series_reference) != len(series_query):
        raise ValueError(
            f"Expected the reference and query to have same length,"
            f"found `{len(series_reference)}`, `{len(series_query)}`."
        )

    if en_cutoff is not None and en_cutoff < 0:
        raise ValueError(
            f"Expected entropy cutoff to be non-negative, found `{en_cutoff}`."
        )

    if not np.all(series_reference.index == series_query.index):
        raise ValueError("Series indices do not match, cannot map names and colors.")

    process_colors = colors_reference is not None

    if not len(series_query):
        res = Series([], dtype="category")
        return (res, []) if process_colors else res

    if process_colors:
        if len(colors_reference) < len(series_reference.cat.categories):
            raise ValueError(
                f"Length of reference colors `{len(colors_reference)}` is smaller than "
                f"length of reference categories `{len(series_reference.cat.categories)}`."
            )
        colors_reference = colors_reference[: len(series_reference.cat.categories)]
        if not all(mcolors.is_color_like(c) for c in colors_reference):
            raise ValueError("Not all values are valid colors.")
        if len(set(colors_reference)) != len(colors_reference):
            logg.warning("Color sequence contains non-unique elements")

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

    # find the mapping which maximizes overlap
    names_query = association_df.T.idxmax()
    if en_cutoff is not None:
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

    # warning: do NOT use np.array - if I pass for e.g. colors ['red'], the dtype will be '<U3'
    # but _create_colors convert them to hex, which will leave them trimmed to #ff or similar
    names_query_new = Series(names_query.copy())
    if process_colors:
        colors_query_new = Series(colors_query.copy())

    # Create unique names by adding suffixes "..._1, ..._2" etc and unique colors by shifting the original color
    for key, value in frequ.items():
        if value == 1:
            continue  # already unique, skip
        # deal with non-unique names
        unique_names = [f"{key}_{rep}" for rep in np.arange(1, value + 1)]
        # .value because of pandas 1.0.0
        names_query_new.iloc[(names_query_series == key).values] = unique_names
        if process_colors:
            color = association_df[association_df["name"] == key]["color"].values[0]
            shifted_colors = _create_colors(color, value, saturation_range=None)
            colors_query_new.iloc[(names_query_series == key).values] = shifted_colors

    # warnings: if it's categorical and assigning to `.cat.categories`, it will
    # take the categorical information, making the 2nd line below necessary
    names_query_new = names_query_new.astype("category")
    names_query_new.cat.reorder_categories(np.array(names_query_new), inplace=True)

    # issue a warning for mapping with high entropy
    if en_cutoff is not None:
        critical_cats = sorted(
            set(
                association_df.loc[association_df["entropy"] > en_cutoff, "name"].values
            )
        )
        if len(critical_cats) > 0:
            logg.warning(
                f"The following states could not be mapped uniquely: `{critical_cats}`"
            )

    return (
        (names_query_new, list(_convert_to_hex_colors(colors_query_new)))
        if process_colors
        else names_query_new
    )


def _compute_mean_color(color_list: List[str]) -> str:
    """Compute mean color."""

    if not all(map(lambda c: mcolors.is_color_like(c), color_list)):
        raise ValueError(f"Not all values are valid colors `{color_list}`.")

    color_list = np.array([mcolors.rgb_to_hsv(mcolors.to_rgb(c)) for c in color_list])

    return mcolors.to_hex(mcolors.hsv_to_rgb(np.mean(color_list, axis=0)))


def _colors_in_order(
    adata,
    clusters: Optional[Iterable[str]] = None,
    cluster_key: str = "clusters",
) -> List[Any]:
    """
    Get list of colors from AnnData in defined order.

    Extracts a list of colors from ``adata.uns[cluster_key]`` in the order defined by the ``clusters``.

    Parameters
    ----------
    %(adata)s
    clusters
        Subset of the clusters we want the color for. Must be a subset of ``adata.obs[cluster_key].cat.categories``.
    cluster_key
        Key from ``adata.obs``.

    Returns
    -------
    list
        List of colors in order defined by `clusters`.
    """

    assert (
        cluster_key in adata.obs.keys()
    ), f"Could not find {cluster_key} in `adata.obs`."

    if clusters is not None:
        assert np.all(
            np.in1d(clusters, adata.obs[cluster_key].cat.categories)
        ), "Not all `clusters` found."

    assert (
        f"{cluster_key}_colors" in adata.uns.keys()
    ), f"No colors associated to {cluster_key} in `adata.uns`."

    if clusters is None:
        clusters = adata.obs[cluster_key].cat.categories

    color_list = []
    all_clusters = adata.obs[cluster_key].cat.categories

    for cl in clusters:
        mask = np.in1d(all_clusters, cl)
        color_list.append(adata.uns[f"{cluster_key}_colors"][mask][0])

    return color_list
