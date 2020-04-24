# -*- coding: utf-8 -*-
from cellrank.tools._utils import _merge_categorical_series

import pytest
import pandas as pd
import numpy as np


class TestToolsUtils:
    def test_merge_not_categorical(self):
        x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")
        y = pd.Series(["b", np.nan, np.nan, "d", "a"])
        with pytest.raises(TypeError):
            _ = _merge_categorical_series(x, y)

    def test_merge_different_index(self):
        x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")
        y = pd.Series(["b", np.nan, np.nan, "d", "a"], index=[5, 4, 3, 2, 1]).astype(
            "category"
        )
        with pytest.raises(ValueError):
            _ = _merge_categorical_series(x, y)

    def test_merge_normal_run(self):
        x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")
        y = pd.Series(["b", np.nan, "a", "d", "a"]).astype("category")
        expected = pd.Series(["b", "b", "a", "d", "a"]).astype("category")

        res = _merge_categorical_series(x, y, inplace=False)

        np.testing.assert_array_equal(res.values, expected.values)

    def test_merge_normal_run_inplace(self):
        x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")
        y = pd.Series(["b", np.nan, "a", "d", "a"]).astype("category")
        expected = pd.Series(["b", "b", "a", "d", "a"]).astype("category")

        _ = _merge_categorical_series(x, y, inplace=True)

        assert _ is None
        np.testing.assert_array_equal(x.values, expected.values)

    def test_merge_normal_run_completely_different_categories(self):
        x = pd.Series(["a", "a", "a"]).astype("category")
        y = pd.Series(["b", "b", "b"]).astype("category")
        expected = pd.Series(["b", "b", "b"]).astype("category")

        res = _merge_categorical_series(x, y, inplace=False)

        np.testing.assert_array_equal(res.values, expected.values)
        np.testing.assert_array_equal(res.cat.categories.values, ["b"])

    def test_merge_colors_not_colorlike(self):
        with pytest.raises(ValueError):
            x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")
            y = pd.Series(["b", np.nan, "a", "d", "a"]).astype("category")
            colors_x = ["red", "foo"]

            _ = _merge_categorical_series(x, y, colors_old=colors_x, inplace=True)

    def test_merge_colors_wrong_number_of_colors(self):
        with pytest.raises(ValueError):
            x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")
            y = pd.Series(["b", np.nan, "a", "d", "a"]).astype("category")
            colors_x = ["red"]

            _ = _merge_categorical_series(x, y, colors_old=colors_x, inplace=True)

    def test_merge_colors_wrong_dict(self):
        with pytest.raises(ValueError):
            x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")
            y = pd.Series(["b", np.nan, "a", "d", "a"]).astype("category")
            colors_x = {"a": "red", "foo": "blue"}

            _ = _merge_categorical_series(x, y, colors_old=colors_x, inplace=True)

    def test_merge_colors_simple_old(self):
        x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")
        y = pd.Series(["b", np.nan, "a", "d", "a"]).astype("category")
        colors_x = ["red", "blue"]

        colors_merged = _merge_categorical_series(
            x, y, colors_old=colors_x, inplace=True
        )

        np.testing.assert_array_equal(colors_merged, ["red", "blue", "#4daf4a"])

    def test_merge_colors_simple_old_no_inplace(self):
        x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")
        y = pd.Series(["b", np.nan, "a", "d", "a"]).astype("category")
        expected = pd.Series(["b", "b", "a", "d", "a"]).astype("category")
        colors_x = ["red", "blue"]

        merged, colors_merged = _merge_categorical_series(
            x, y, colors_old=colors_x, inplace=False
        )

        np.testing.assert_array_equal(merged.values, expected.values)
        np.testing.assert_array_equal(colors_merged, ["red", "blue", "#4daf4a"])

    def test_merge_colors_simple_new(self):
        x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")
        y = pd.Series(["b", np.nan, "a", "d", "a"]).astype("category")
        colors_y = ["red", "blue", "green"]

        colors_merged = _merge_categorical_series(
            x, y, colors_new=colors_y, inplace=True
        )

        np.testing.assert_array_equal(colors_merged, ["#e41a1c", "#377eb8", "green"])

    def test_merge_colors_both(self):
        x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")
        y = pd.Series(["b", np.nan, "a", "d", "a"]).astype("category")
        colors_x = ["red", "blue"]
        colors_y = ["green", "yellow", "black"]

        colors_merged = _merge_categorical_series(
            x, y, colors_old=colors_x, colors_new=colors_y, inplace=True
        )

        np.testing.assert_array_equal(colors_merged, ["red", "blue", "black"])

    def test_merge_colors_both_overwrite(self):
        x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")
        y = pd.Series(["b", np.nan, "a", "d", "a"]).astype("category")
        colors_x = ["red", "blue"]
        colors_y = ["green", "yellow", "black"]

        colors_merged = _merge_categorical_series(
            x,
            y,
            colors_old=colors_x,
            colors_new=colors_y,
            color_overwrite=True,
            inplace=True,
        )

        np.testing.assert_array_equal(colors_merged, ["green", "yellow", "black"])
