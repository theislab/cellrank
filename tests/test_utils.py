# -*- coding: utf-8 -*-
from cellrank.tools._utils import _merge_approx_rcs, _map_names_and_colors

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

        with pytest.raises(ValueError):
            _ = _map_names_and_colors(x, y)

    def test_simple_wrong_index(self):
        x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")
        y = pd.Series(
            ["b", np.nan, np.nan, "d", "a"], index=["foo", "bar", "baz", "quux", "quas"]
        )

        with pytest.raises(ValueError):
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
