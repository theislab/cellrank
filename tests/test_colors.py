import numpy as np
import pandas as pd
import pytest
from matplotlib.colors import is_color_like

from cellrank._utils._colors import _create_categorical_colors, _map_names_and_colors


class TestColors:
    def test_create_categorical_colors_too_many_colors(self):
        with pytest.raises(ValueError, match=r".* exceeded the maximum"):
            _create_categorical_colors(1000)

    def test_create_categorical_colors_no_categories(self):
        c = _create_categorical_colors(0)

        assert c == []

    def test_create_categorical_colors_neg_categories(self):
        with pytest.raises(RuntimeError, match="Unable to create"):
            _create_categorical_colors(-1)

    def test_create_categorical_colors_normal_run(self):
        colors = _create_categorical_colors(62)

        assert len(colors) == 62
        assert all(isinstance(c, str) for c in colors), colors
        assert all(is_color_like(c) for c in colors), colors


class TestMappingColors:
    def test_mapping_colors_not_categorical(self):
        query = pd.Series(["foo", "bar", "baz"], dtype="str")
        reference = pd.Series(["foo", np.nan, "bar", "baz"], dtype="category")

        with pytest.raises(TypeError, match=r"Query series must be"):
            _map_names_and_colors(reference, query)

    def test_mapping_colors_invalid_size(self):
        query = pd.Series(["foo", "bar", "baz"], dtype="category")
        reference = pd.Series(["foo", np.nan, "bar", "baz"], dtype="category")

        with pytest.raises(ValueError, match=r".*to have the same length"):
            _map_names_and_colors(reference, query)

    def test_mapping_colors_different_index(self):
        query = pd.Series(["foo", "bar", "baz"], dtype="category", index=[2, 3, 4])
        reference = pd.Series(["foo", "bar", "baz"], dtype="category", index=[1, 2, 3])

        with pytest.raises(ValueError, match=r"Series indices do not match"):
            _map_names_and_colors(reference, query)

    def test_mapping_colors_invalid_colors(self):
        query = pd.Series(["foo", "bar", "baz"], dtype="category")
        reference = pd.Series(["foo", "bar", "baz"], dtype="category")

        with pytest.raises(ValueError, match=r"Not all values are valid colors"):
            _map_names_and_colors(reference, query, colors_reference=["red", "green", "foo"])

    def test_mapping_colors_too_few_colors(self):
        query = pd.Series(["foo", "bar", "baz"], dtype="category")
        reference = pd.Series(["foo", "bar", "baz"], dtype="category")

        with pytest.raises(ValueError, match=r"Length of reference colors"):
            _map_names_and_colors(reference, query, colors_reference=["red", "green"])

    def test_mapping_colors_simple_1(self):
        x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")
        y = pd.Series(["b", np.nan, np.nan, "d", "a"]).astype("category")
        expected = pd.Series(["a_1", "a_2", "b"])
        expected_index = pd.Index(["a", "b", "d"])

        res = _map_names_and_colors(x, y)

        assert isinstance(res, pd.Series)
        np.testing.assert_array_equal(res.values, expected.values)
        np.testing.assert_array_equal(res.index.values, expected_index.values)

    def test_mapping_colors_simple_2(self):
        query = pd.Series(["foo", "bar", "baz"], dtype="category")
        reference = pd.Series(["foo", "bar", "baz"], dtype="category")

        res = _map_names_and_colors(reference, query)

        assert isinstance(res, pd.Series)
        assert len(res) == 3
        assert isinstance(res.dtype, pd.CategoricalDtype)

    def test_mapping_colors_simple_colors(self):
        query = pd.Series(["foo", "bar", "baz"], dtype="category")
        reference = pd.Series(["foo", "bar", "baz"], dtype="category")

        res, c = _map_names_and_colors(reference, query, colors_reference=["red", "green", "blue"])

        assert isinstance(res, pd.Series)
        assert len(res) == 3
        assert isinstance(res.dtype, pd.CategoricalDtype)

        assert isinstance(c, list)
        assert c == ["#ff0000", "#008000", "#0000ff"]

    def test_mapping_colors_too_many_colors(self):
        query = pd.Series(["foo", "bar", "baz"], dtype="category")
        reference = pd.Series(["foo", "bar", "baz"], dtype="category")

        res, c = _map_names_and_colors(reference, query, colors_reference=["red", "green", "blue", "black"])

        assert isinstance(res, pd.Series)
        assert len(res) == 3
        assert isinstance(res.dtype, pd.CategoricalDtype)

        assert isinstance(c, list)
        assert c == ["#ff0000", "#008000", "#0000ff"]

    def test_mapping_colors_different_color_representation(self):
        query = pd.Series(["foo", "bar", "baz"], dtype="category")
        reference = pd.Series(["foo", "bar", "baz"], dtype="category")

        res, c = _map_names_and_colors(reference, query, colors_reference=[(1, 0, 0), "green", (0, 0, 1, 0)])

        assert isinstance(res, pd.Series)
        assert len(res) == 3
        assert isinstance(res.dtype, pd.CategoricalDtype)

        assert isinstance(c, list)
        assert c == ["#ff0000", "#008000", "#0000ff"]

    def test_mapping_colors_non_unique_colors(self):
        query = pd.Series(["foo", "bar", "baz"], dtype="category")
        reference = pd.Series(["foo", "bar", "baz"], dtype="category")

        res, c = _map_names_and_colors(reference, query, colors_reference=["red", "red", "red"])

        assert isinstance(res, pd.Series)
        assert len(res) == 3
        assert isinstance(res.dtype, pd.CategoricalDtype)

        assert isinstance(c, list)
        assert c == ["#ff0000", "#ff0000", "#ff0000"]

    def test_mapping_colors_same_reference(self):
        query = pd.Series(["foo", "bar", "baz"], dtype="category")
        reference = pd.Series(["foo", "foo", "foo"], dtype="category")

        r, c = _map_names_and_colors(reference, query, colors_reference=["red", "red", "red"])

        assert list(r.index) == ["bar", "baz", "foo"]
        assert list(r.values) == ["foo_1", "foo_2", "foo_3"]
        assert c == ["#b20000", "#d13200", "#f07300"]

    def test_mapping_colors_diff_query_reference(self):
        query = pd.Series(["bar", "bar", "bar"], dtype="category")
        reference = pd.Series(["foo", "foo", "foo"], dtype="category")

        r, c = _map_names_and_colors(reference, query, colors_reference=["red", "red", "red"])

        assert list(r.index) == ["bar"]
        assert list(r.values) == ["foo"]
        assert c == ["#ff0000"]

    def test_mapping_colors_empty(self):
        query = pd.Series([], dtype="category")
        reference = pd.Series([], dtype="category")

        r = _map_names_and_colors(reference, query)

        assert isinstance(r, pd.Series)
        assert isinstance(r.dtype, pd.CategoricalDtype)

    def test_mapping_colors_empty_with_color(self):
        query = pd.Series([], dtype="category")
        reference = pd.Series([], dtype="category")

        r, c = _map_names_and_colors(reference, query, colors_reference=[])

        assert isinstance(r, pd.Series)
        assert isinstance(r.dtype, pd.CategoricalDtype)
        assert isinstance(c, list)
        assert len(c) == 0

    def test_mapping_colors_negative_en_cutoff(self):
        query = pd.Series(["foo", "bar", "baz"], dtype="category")
        reference = pd.Series(["foo", "bar", "baz"], dtype="category")

        with pytest.raises(ValueError, match=".* entropy cutoff to be non-negative"):
            _map_names_and_colors(reference, query, en_cutoff=-1)

    def test_mapping_colors_0_en_cutoff(self):
        query = pd.Series(["bar", "bar", "bar"], dtype="category")
        reference = pd.Series(["bar", "bar", "bar"], dtype="category")

        # TODO: somehow extract the custom logger and check for logs
        r = _map_names_and_colors(reference, query, en_cutoff=0)

        assert isinstance(r, pd.Series)
        assert isinstance(r.dtype, pd.CategoricalDtype)
        assert list(r.index) == ["bar"]
        assert list(r.values) == ["bar"]

    def test_mapping_colors_merging(self):
        x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")
        y = pd.Series(["b", np.nan, np.nan, "d", "a"]).astype("category")

        res, colors = _map_names_and_colors(x, y, colors_reference=["red", "green"])

        assert isinstance(res, pd.Series)
        assert isinstance(colors, list)
        np.testing.assert_array_equal(colors, ["#b20000", "#e65c00", "#008000"])

    def test_mapping_colors_merging_more(self):
        x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")
        y = pd.Series(["b", np.nan, np.nan, "d", "a"]).astype("category")

        res, colors = _map_names_and_colors(x, y, colors_reference=["red", "green", "blue", "yellow"])

        assert isinstance(res, pd.Series)
        assert isinstance(colors, list)
        np.testing.assert_array_equal(colors, ["#b20000", "#e65c00", "#008000"])

    def test_mapping_colors_name_order_same_as_cat_order(self):
        x = pd.Series(["b", "a", np.nan, "a", np.nan]).astype("category")
        y = pd.Series(["a", np.nan, np.nan, "d", "b"]).astype("category")
        expected = pd.Series(["b", "a_1", "a_2"])
        expected_index = pd.Index(["a", "b", "d"])

        res = _map_names_and_colors(x, y)

        assert isinstance(res, pd.Series)
        assert isinstance(res.dtype, pd.CategoricalDtype)
        np.testing.assert_array_equal(res.values, expected.values)
        np.testing.assert_array_equal(res.index.values, expected_index.values)
        np.testing.assert_array_equal(res.cat.categories.values, res.values)
