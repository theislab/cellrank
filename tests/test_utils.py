# -*- coding: utf-8 -*-
from cellrank.tools._utils import (
    _merge_approx_rcs,
    _map_names_and_colors,
    _process_series,
)
from _helpers import assert_array_nan_equal

from pandas.api.types import is_categorical_dtype

import pytest
import pandas as pd
import numpy as np


class TestToolsUtils:
    def test_merge_rcs_not_categorical(self):
        x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")
        y = pd.Series(["b", np.nan, np.nan, "d", "a"])
        with pytest.raises(TypeError):
            _ = _merge_approx_rcs(x, y)

    def test_merge_rcs_different_index(self):
        x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")
        y = pd.Series(["b", np.nan, np.nan, "d", "a"], index=[5, 4, 3, 2, 1]).astype(
            "category"
        )
        with pytest.raises(ValueError):
            _ = _merge_approx_rcs(x, y)

    def test_merge_rcsn_ormal_run(self):
        x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")
        y = pd.Series(["b", np.nan, "a", "d", "a"]).astype("category")
        expected = pd.Series(["b", "b", "a", "d", "a"]).astype("category")

        res = _merge_approx_rcs(x, y, inplace=False)

        np.testing.assert_array_equal(res.values, expected.values)

    def test_merge_rcs_normal_run_inplace(self):
        x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")
        y = pd.Series(["b", np.nan, "a", "d", "a"]).astype("category")
        expected = pd.Series(["b", "b", "a", "d", "a"]).astype("category")

        _ = _merge_approx_rcs(x, y, inplace=True)

        assert _ is None
        np.testing.assert_array_equal(x.values, expected.values)

    def test_merge_rcs_normal_run_completely_different_categories(self):
        x = pd.Series(["a", "a", "a"]).astype("category")
        y = pd.Series(["b", "b", "b"]).astype("category")
        expected = pd.Series(["b", "b", "b"]).astype("category")

        res = _merge_approx_rcs(x, y, inplace=False)

        np.testing.assert_array_equal(res.values, expected.values)
        np.testing.assert_array_equal(res.cat.categories.values, ["b"])


class TestMapNamesAndColors:
    def test_simple_not_categorical(self):
        x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")
        y = pd.Series(["b", np.nan, np.nan, "d", "a"])

        with pytest.raises(TypeError):
            _ = _map_names_and_colors(x, y)

    def test_simple_wrong_index(self):
        x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")
        y = pd.Series(
            ["b", np.nan, np.nan, "d", "a"], index=["foo", "bar", "baz", "quux", "quas"]
        )

        with pytest.raises(TypeError):
            _ = _map_names_and_colors(x, y)

    def test_simple_not_color_like(self):
        x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")
        y = pd.Series(["b", np.nan, np.nan, "d", "a"]).astype("category")

        with pytest.raises(ValueError):
            _ = _map_names_and_colors(x, y, colors_reference=["foo", "bar"])

    def test_simple_not_invalid_color_length(self):
        x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")
        y = pd.Series(["b", np.nan, np.nan, "d", "a"]).astype("category")

        with pytest.raises(ValueError):
            _ = _map_names_and_colors(x, y, colors_reference=["red"])

    def test_simple_run(self):
        x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")
        y = pd.Series(["b", np.nan, np.nan, "d", "a"]).astype("category")
        expected = pd.Series(["a_1", "a_2", "b"])
        expected_index = pd.Index(["a", "b", "d"])

        res = _map_names_and_colors(x, y)

        assert isinstance(res, pd.Series)
        np.testing.assert_array_equal(res.values, expected.values)
        np.testing.assert_array_equal(res.index.values, expected_index.values)

    def test_return_colors(self):
        x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")
        y = pd.Series(["b", np.nan, np.nan, "d", "a"]).astype("category")

        res, colors = _map_names_and_colors(x, y)

        assert isinstance(res, pd.Series)
        assert isinstance(colors, list)
        np.testing.assert_array_equal(res, ["#bb2200", "#ee6655", "#008000"])


class TestProcessSeries:
    def test_not_categorical(self):
        x = pd.Series(["a", "b", np.nan, "b", np.nan])

        with pytest.raises(TypeError):
            _ = _process_series(x, ["foo"])

    def test_colors_wrong_number_of_colors(self):
        x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")

        with pytest.raises(ValueError):
            _ = _process_series(x, ["foo"], colors=["red"])

    def test_colors_not_colorlike(self):
        x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")

        with pytest.raises(ValueError):
            _ = _process_series(x, ["foo"], colors=["bar"])

    def test_keys_are_not_proper_categories(self):
        x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")

        with pytest.raises(ValueError):
            _ = _process_series(x, ["foo"])

    def test_keys_overlap(self):
        x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")

        with pytest.raises(ValueError):
            _ = _process_series(x, ["a", "b, a"])

    def test_normal_run(self):
        x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")
        expected = pd.Series(["a"] + [np.nan] * 4).astype("category")

        res = _process_series(x, keys=["a"])

        assert_array_nan_equal(expected, res)

    def test_repeat_key(self):
        x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")
        expected = pd.Series(["a"] + [np.nan] * 4).astype("category")

        res = _process_series(x, keys=["a, a, a"])

        assert_array_nan_equal(res, expected)

    def test_reoder_keys(self):
        x = pd.Series(["b", "c", "a", "d", "a"]).astype("category")
        expected = pd.Series(["a or b or d", np.nan] + ["a or b or d"] * 3).astype(
            "category"
        )

        res = _process_series(x, keys=["b, a, d"])

        assert_array_nan_equal(res, expected)

    def test_no_keys(self):
        x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")

        res = _process_series(x, keys=None)

        assert x is res

    def test_no_keys_colors(self):
        x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")
        colors = ["foo"]

        res, res_colors = _process_series(x, keys=None, colors=colors)

        assert x is res
        assert colors is res_colors

    def test_empty_keys(self):
        x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")

        res = _process_series(x, [])

        assert res.shape == x.shape
        assert np.all(pd.isnull(res))

    def test_return_colors(self):
        x = pd.Series(["b", "c", "a", "d", "a"]).astype("category")
        expected = pd.Series(["a or b", "c or d", "a or b", "c or d", "a or b"]).astype(
            "category"
        )

        res, colors = _process_series(
            x, keys=["b, a", "d, c"], colors=["red", "green", "blue", "white"]
        )

        assert isinstance(res, pd.Series)
        assert is_categorical_dtype(res)
        assert isinstance(colors, list)

        np.testing.assert_array_equal(res.values, expected.values)
        np.testing.assert_array_equal(colors, ["#804000", "#8080ff"])
