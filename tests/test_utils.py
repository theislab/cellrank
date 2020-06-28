# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import pytest
from pandas.api.types import is_categorical_dtype

import scanpy as sc

from _helpers import assert_array_nan_equal
from cellrank.tools import Lineage
from cellrank.tools._utils import (
    _one_hot,
    _cluster_X,
    _process_series,
    _fuzzy_to_discrete,
    _merge_categorical_series,
    _series_from_one_hot_matrix,
)
from cellrank.tools._colors import _map_names_and_colors


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

        res, colors = _map_names_and_colors(x, y, colors_reference=["red", "green"])

        assert isinstance(res, pd.Series)
        assert isinstance(colors, list)
        np.testing.assert_array_equal(colors, ["#bb2200", "#ee6655", "#008000"])


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
        assert set(colors) == {"#804000", "#8080ff"}


class TestOneHot:
    def test_normal_run(self):
        _one_hot(n=10, cat=5)
        _one_hot(n=10, cat=None)

    def test_return_vector(self):
        a = _one_hot(n=10, cat=None)
        b = _one_hot(n=10, cat=5)

        b_check = np.zeros(10)
        b_check[5] = True

        assert a.dtype == "bool"
        assert b.dtype == "bool"
        assert (a == np.zeros(10)).all()
        assert (b == b_check).all()

    def test_index_error(self):
        with pytest.raises(IndexError):
            _ = _one_hot(10, 10)


class TestFuzzyToDiscrete:
    def test_normal_run(self):
        # create random data that sums to one row-wise
        a_fuzzy = np.random.standard_normal((100, 3))
        a_fuzzy = np.exp(a_fuzzy) / np.sum(np.exp(a_fuzzy), 1)[:, None]

        # check with both overlap handlings
        _fuzzy_to_discrete(a_fuzzy=a_fuzzy)
        _fuzzy_to_discrete(a_fuzzy=a_fuzzy, n_most_likely=30, remove_overlap=True)
        _fuzzy_to_discrete(a_fuzzy=a_fuzzy, n_most_likely=30, remove_overlap=False)

    def test_one_state(self):
        # create random data that sums to one row-wise
        a_fuzzy = np.random.standard_normal((100, 1))
        a_fuzzy = np.exp(a_fuzzy) / np.sum(np.exp(a_fuzzy), 1)[:, None]

        # check with both overlap handlings
        _fuzzy_to_discrete(a_fuzzy=a_fuzzy)

    def test_normalization(self):
        a_fuzzy = np.random.standard_normal((100, 3))
        with pytest.raises(ValueError):
            _fuzzy_to_discrete(a_fuzzy=a_fuzzy)

    def test_too_many_cells(self):
        a_fuzzy = np.random.standard_normal((100, 3))
        a_fuzzy = np.exp(a_fuzzy) / np.sum(np.exp(a_fuzzy), 1)[:, None]
        with pytest.raises(ValueError):
            _fuzzy_to_discrete(a_fuzzy=a_fuzzy, n_most_likely=50)

    def test_raise_threshold(self):
        a_fuzzy = np.repeat(np.array([0.9, 0.1])[None, :], 10, 0)
        with pytest.raises(ValueError):
            _fuzzy_to_discrete(a_fuzzy, n_most_likely=3, remove_overlap=True)
        with pytest.raises(ValueError):
            _fuzzy_to_discrete(a_fuzzy, n_most_likely=3, remove_overlap=False)

    def test_normal_output(self):
        a_fuzzy = np.array(
            [
                [0.3, 0.7, 0],
                [0.2, 0.5, 0.3],
                [0.1, 0.8, 0.1],
                [0.4, 0.4, 0.2],
                [0.5, 0.3, 0.2],
                [0.6, 0.3, 0.1],
                [0.3, 0.3, 0.4],
                [0.2, 0.2, 0.6],
            ]
        )

        # note: removing the overlap should have no effect in this case since there is none.
        # there should also be no critical clusters in this case
        a_actual_1, c_1 = _fuzzy_to_discrete(
            a_fuzzy, n_most_likely=2, remove_overlap=True
        )
        a_actual_2, c_2 = _fuzzy_to_discrete(
            a_fuzzy, n_most_likely=2, remove_overlap=False
        )
        a_expected = np.array(
            [
                [False, True, False],
                [False, False, False],
                [False, True, False],
                [False, False, False],
                [True, False, False],
                [True, False, False],
                [False, False, True],
                [False, False, True],
            ]
        )

        np.testing.assert_array_equal(a_actual_1, a_expected)
        np.testing.assert_array_equal(a_actual_2, a_expected)
        assert len(c_1) == 0
        assert len(c_2) == 0

    def test_critical_samples(self):
        a_fuzzy = np.array(
            [
                [0.3, 0.7, 0],
                [0.3, 0.6, 0.1],
                [0.0, 0.7, 0.3],
                [0.1, 0.9, 0],
                [0.4, 0.4, 0.2],
                [0.5, 0.3, 0.2],
                [0.6, 0.3, 0.1],
                [0.3, 0.3, 0.4],
                [0.2, 0.2, 0.6],
            ]
        )

        _, c_1 = _fuzzy_to_discrete(a_fuzzy, n_most_likely=3, remove_overlap=False)
        _, c_2 = _fuzzy_to_discrete(a_fuzzy, n_most_likely=3, remove_overlap=True)

        assert c_1 == np.array(2)
        np.testing.assert_array_equal(c_2, np.array([1, 2]))

    def test_passing_lineage_object(self):
        a_fuzzy = np.array(
            [
                [0.3, 0.7, 0],
                [0.2, 0.5, 0.3],
                [0.1, 0.8, 0.1],
                [0.4, 0.4, 0.2],
                [0.5, 0.3, 0.2],
                [0.6, 0.3, 0.1],
                [0.3, 0.3, 0.4],
                [0.2, 0.2, 0.6],
            ]
        )
        a_fuzzy_lin = Lineage(a_fuzzy, names=["0", "1", "2"])

        b_np, c_np = _fuzzy_to_discrete(a_fuzzy=a_fuzzy, n_most_likely=2)
        b_l, c_l = _fuzzy_to_discrete(a_fuzzy=a_fuzzy_lin, n_most_likely=2)

        np.testing.assert_array_equal(b_np, b_l)
        assert len(c_np) == 0
        assert len(c_l) == 0


class TestSeriesFromOneHotMatrix:
    def test_normal_run(self):
        a = np.array(
            [[0, 0, 1], [0, 0, 1], [0, 0, 0], [1, 0, 0], [1, 0, 0], [0, 1, 0]],
            dtype="bool",
        )
        res = _series_from_one_hot_matrix(a)

        assert_array_nan_equal(
            np.array(res).astype(np.float32),
            np.array([2, 2, np.nan, 0, 0, 1], dtype=np.float32),
        )
        np.testing.assert_array_equal(res.cat.categories, ["0", "1", "2"])

    def test_name_mismatch(self):
        a = np.array(
            [[0, 0, 1], [0, 0, 1], [0, 0, 0], [1, 0, 0], [1, 0, 0], [0, 1, 0]],
            dtype="bool",
        )
        names = ["0", "1"]

        with pytest.raises(ValueError):
            _series_from_one_hot_matrix(a, names=names)

    def test_dtype(self):
        a = np.array(
            [[0, 0, 1], [0, 0, 2], [0, 0, 0], [1, 0, 0], [1, 0, 0], [0, 1, 0]],
            dtype="int8",
        )

        with pytest.raises(TypeError):
            _series_from_one_hot_matrix(a)

    def test_not_one_hot(self):
        a = np.array(
            [[1, 0, 1], [0, 0, 1], [0, 0, 0], [1, 0, 0], [1, 0, 0], [0, 1, 0]],
            dtype="bool",
        )

        with pytest.raises(ValueError):
            _series_from_one_hot_matrix(a)

    def test_normal_return(self):
        a = np.array(
            [[0, 0, 1], [0, 0, 1], [0, 0, 0], [1, 0, 0], [1, 0, 0], [0, 1, 0]],
            dtype="bool",
        )
        actual_series = _series_from_one_hot_matrix(a)

        expected_series = pd.Series(index=range(6), dtype="category")
        expected_series.cat.add_categories(["0", "1", "2"], inplace=True)
        expected_series[0] = "2"
        expected_series[1] = "2"
        expected_series[3] = "0"
        expected_series[4] = "0"
        expected_series[5] = "1"

        assert actual_series.equals(expected_series)
        assert (actual_series.cat.categories == expected_series.cat.categories).all()


class TestClusterX:
    def test_normal_run(self):
        # create some data
        adata = sc.datasets.blobs(n_observations=100, n_variables=6)

        # kmeans, louvain, leiden
        labels_kmeans = _cluster_X(adata.X, n_clusters=5, method="kmeans")
        labels_louvain = _cluster_X(adata.X, n_clusters=5, method="louvain")

        assert len(labels_kmeans) == len(labels_louvain) == adata.n_obs

    def test_one_feature(self):
        # create some data
        adata = sc.datasets.blobs(n_observations=100, n_variables=1)

        # kmeans, louvain, leiden
        labels_kmeans = _cluster_X(adata.X, n_clusters=5, method="kmeans")
        labels_louvain = _cluster_X(adata.X, n_clusters=5, method="louvain")

        assert len(labels_kmeans) == len(labels_louvain) == adata.n_obs
