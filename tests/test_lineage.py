# -*- coding: utf-8 -*-
import mock
import numpy as np
import pytest
from pandas import DataFrame

import matplotlib.colors as colors

import cellrank.tools._lineage as mocker
from cellrank.tools import Lineage
from cellrank.tools._utils import _compute_mean_color
from cellrank.tools._colors import _create_categorical_colors
from cellrank.tools._constants import Lin


class TestLineageCreation:
    def test_creation(self):
        x = np.random.random((10, 3))
        names = ["foo", "bar", "baz"]
        colors = ["#000000", "#ababab", "#ffffff"]
        l = Lineage(x, names=names, colors=colors)

        np.testing.assert_array_equal(l, x)
        np.testing.assert_array_equal(l.names, np.array(names))
        np.testing.assert_array_equal(l.colors, np.array(colors))

    def test_wrong_number_of_dimensions(self):
        with pytest.raises(ValueError):
            _ = Lineage(
                np.random.random((10, 3, 1)),
                names=["foo", "bar", "baz"],
                colors=[(0, 0, 0), "#ffffff", "#ff00FF"],
            )

    def test_names_length_mismatch(self):
        with pytest.raises(ValueError):
            _ = Lineage(
                np.random.random((10, 3)),
                names=["foo", "bar"],
                colors=[(0, 0, 0), (0.5, 0.5, 0.5), "foobar"],
            )

    def test_colors_length_mismatch(self):
        with pytest.raises(ValueError):
            _ = Lineage(
                np.random.random((10, 3)),
                names=["foo", "bar", "baz"],
                colors=[(0, 0, 0), (0.5, 0.5, 0.5)],
            )

    def test_wrong_colors(self):
        with pytest.raises(ValueError):
            _ = Lineage(
                np.random.random((10, 3)),
                names=["foo", "bar", "baz"],
                colors=[(0, 0, 0), (0.5, 0.5, 0.5), "foobar"],
            )

    def test_colors_setter(self):
        l = Lineage(
            np.random.random((10, 3)),
            names=["foo", "bar", "baz"],
            colors=[(0, 0, 0), (0.5, 0.5, 0.5), (1, 1, 1)],
        )

        colors = ["#ffffff", "#ffffff", "#ffffff"]
        l.colors = colors

        np.testing.assert_array_equal(l.colors, np.array(colors))

    def test_color_setter_wrong_colors(self):
        l = Lineage(
            np.random.random((10, 3)),
            names=["foo", "bar", "baz"],
            colors=[(0, 0, 0), (0.5, 0.5, 0.5), (1, 1, 1)],
        )

        with pytest.raises(ValueError):
            l.colors = ["#ffffff", "#ffffff", "foo"]

    def test_names_setter(self):
        l = Lineage(
            np.random.random((10, 3)),
            names=["foo", "bar", "baz"],
            colors=[(0, 0, 0), (0.5, 0.5, 0.5), (1, 1, 1)],
        )

        names = ["foo1", "bar1", "baz1"]
        l.names = names

        np.testing.assert_array_equal(l.names, np.array(names))

    def test_names_setter_wrong_type(self):
        l = Lineage(
            np.random.random((10, 3)),
            names=["foo", "bar", "baz"],
            colors=[(0, 0, 0), (0.5, 0.5, 0.5), (1, 1, 1)],
        )

        with pytest.raises(TypeError):
            l.names = ["foo1", "bar1", 3]

    def test_names_setter_wrong_size(self):
        l = Lineage(
            np.random.random((10, 3)),
            names=["foo", "bar", "baz"],
            colors=[(0, 0, 0), (0.5, 0.5, 0.5), (1, 1, 1)],
        )

        with pytest.raises(ValueError):
            l.names = ["foo1", "bar1"]

    def test_names_setter_non_unique(self):
        l = Lineage(
            np.random.random((10, 3)),
            names=["foo", "bar", "baz"],
            colors=[(0, 0, 0), (0.5, 0.5, 0.5), (1, 1, 1)],
        )

        with pytest.raises(ValueError):
            l.names = ["foo1", "bar1", "bar1"]

    def test_non_unique_names(self):
        with pytest.raises(ValueError):
            _ = Lineage(
                np.random.random((10, 3, 1)),
                names=["foo", "bar", "baz"],
                colors=[(0, 0, 0), "#ffffff", "#ff00FF"],
            )


class TestLineageAccessor:
    def test_subset_same_instance(self):
        x = np.random.random((10, 3))
        l = Lineage(
            x,
            names=["foo", "bar", "baz"],
            colors=[(0, 0, 0), (0.5, 0.5, 0.5), (1, 1, 1)],
        )

        y = l[0, 0]

        assert isinstance(y, Lineage)

    def test_singleton_column(self):
        x = np.random.random((10, 3))
        l = Lineage(
            x,
            names=["foo", "bar", "baz"],
            colors=[(0, 0, 0), (0.5, 0.5, 0.5), (1, 1, 1)],
        )

        y = l[:, 0]

        np.testing.assert_array_equal(x[:, 0], np.array(y)[:, 0])

    def test_singleton_column_name(self):
        x = np.random.random((10, 3))
        l = Lineage(
            x,
            names=["foo", "bar", "baz"],
            colors=[(0, 0, 0), (0.5, 0.5, 0.5), (1, 1, 1)],
        )

        y = l[:, "foo"]

        np.testing.assert_array_equal(x[:, 0], np.array(y)[:, 0])

    def test_singleton_column_first_index_assignment(self):
        x = np.random.random((10, 3))
        l = Lineage(
            x,
            names=["foo", "bar", "baz"],
            colors=[(0, 0, 0), (0.5, 0.5, 0.5), (1, 1, 1)],
        )

        y = l["baz"]

        np.testing.assert_array_equal(x[:, 2], np.array(y)[:, 0])
        np.testing.assert_array_equal(y.names, ["baz"])

    def test_singleton_row_and_column(self):
        x = np.random.random((10, 3))
        l = Lineage(
            x,
            names=["foo", "bar", "baz"],
            colors=[(0, 0, 0), (0.5, 0.5, 0.5), (1, 1, 1)],
        )

        y = l[0, "foo"]

        assert isinstance(l, Lineage)
        assert y.shape == (1, 1)
        assert x[0, 0] == y[0, 0]

    def test_mixed_columns(self):
        x = np.random.random((10, 3))
        l = Lineage(
            x,
            names=["foo", "bar", "baz"],
            colors=[(0, 0, 0), (0.5, 0.5, 0.5), (1, 1, 1)],
        )

        y = l[0, ["foo", 2, "bar"]]

        np.testing.assert_array_equal(x[[[0]], [0, 2, 1]], np.array(y))

    def test_remove_duplicates(self):
        x = np.random.random((10, 3))
        l = Lineage(
            x,
            names=["foo", "bar", "baz"],
            colors=[(0, 0, 0), (0.5, 0.5, 0.5), (1, 1, 1)],
        )

        y = l[0, ["foo", 2, "bar", 0, 0, "foo"]]

        np.testing.assert_array_equal(x[[[0]], [0, 2, 1]], np.array(y))

    def test_column_invalid_name(self):
        x = np.random.random((10, 3))
        l = Lineage(
            x,
            names=["foo", "bar", "baz"],
            colors=[(0, 0, 0), (0.5, 0.5, 0.5), (1, 1, 1)],
        )

        with pytest.raises(KeyError):
            y = l["quux"]

    def test_row_subset_with_ints(self):
        x = np.random.random((10, 3))
        l = Lineage(
            x,
            names=["foo", "bar", "baz"],
            colors=[(0, 0, 0), (0.5, 0.5, 0.5), (1, 1, 1)],
        )

        y = l[[1, 2, 3], :]

        np.testing.assert_array_equal(x[[1, 2, 3], :], np.array(y))

    def test_column_subset_boolean(self):
        x = np.random.random((10, 3))
        l = Lineage(
            x,
            names=["foo", "bar", "baz"],
            colors=[(0, 0, 0), (0.5, 0.5, 0.5), (1, 1, 1)],
        )

        y = l[:, [False, False, True]]

        np.testing.assert_array_equal(x[:, -1], y.X.squeeze())

    def test_column_subset_boolean_invalid_dim(self):
        x = np.random.random((10, 3))
        l = Lineage(
            x,
            names=["foo", "bar", "baz"],
            colors=[(0, 0, 0), (0.5, 0.5, 0.5), (1, 1, 1)],
        )

        with pytest.raises(IndexError):
            y = l[:, [True]]

    def test_row_subset_with_mask(self):
        x = np.random.random((10, 3))
        l = Lineage(
            x,
            names=["foo", "bar", "baz"],
            colors=[(0, 0, 0), (0.5, 0.5, 0.5), (1, 1, 1)],
        )

        mask = np.ones((x.shape[0]), dtype=np.bool)
        mask[:5] = False
        y = l[mask, :]

        np.testing.assert_array_equal(x[mask, :], np.array(y))

    def test_column_subset_with_ints(self):
        x = np.random.random((10, 3))
        l = Lineage(
            x,
            names=["foo", "bar", "baz"],
            colors=[(0, 0, 0), (0.5, 0.5, 0.5), (1, 1, 1)],
        )

        y = l[:, [2, 0]]

        np.testing.assert_array_equal(x[:, [2, 0]], np.array(y))

    def test_column_subset_with_mask(self):
        x = np.random.random((10, 3))
        l = Lineage(
            x,
            names=["foo", "bar", "baz"],
            colors=[(0, 0, 0), (0.5, 0.5, 0.5), (1, 1, 1)],
        )

        mask = np.ones((x.shape[1]), dtype=np.bool)
        mask[0] = False
        y = l[:, mask]

        np.testing.assert_array_equal(x[:, mask], np.array(y))

    def test_column_subset_with_names(self):
        x = np.random.random((10, 3))
        l = Lineage(
            x,
            names=["foo", "bar", "baz"],
            colors=[(0, 0, 0), (0.5, 0.5, 0.5), (1, 1, 1)],
        )

        y = l[:, ["foo", "bar"]]

        np.testing.assert_array_equal(x[:, [0, 1]], np.array(y))

    def test_comb_row_int_col_int(self):
        x = np.random.random((10, 3))
        l = Lineage(
            x,
            names=["foo", "bar", "baz"],
            colors=[(0, 0, 0), (0.5, 0.5, 0.5), (1, 1, 1)],
        )

        y = l[[0, 1], [1, 2]]

        np.testing.assert_array_equal(x[[0, 1], :][:, [1, 2]], np.array(y))

    def test_comb_row_int_col_mask(self):
        x = np.random.random((10, 3))
        l = Lineage(
            x,
            names=["foo", "bar", "baz"],
            colors=[(0, 0, 0), (0.5, 0.5, 0.5), (1, 1, 1)],
        )

        mask = np.ones((x.shape[1]), dtype=np.bool)
        mask[0] = False
        y = l[[0, 1], mask]

        np.testing.assert_array_equal(x[[0, 1], :][:, mask], np.array(y))

    def test_comb_row_int_col_names(self):
        x = np.random.random((10, 3))
        l = Lineage(
            x,
            names=["foo", "bar", "baz"],
            colors=[(0, 0, 0), (0.5, 0.5, 0.5), (1, 1, 1)],
        )

        y = l[[0, 1], "baz"]

        np.testing.assert_array_equal(x[[0, 1], :][:, [2]], np.array(y))

    def test_comb_row_mask_col_int(self):
        x = np.random.random((10, 3))
        l = Lineage(
            x,
            names=["foo", "bar", "baz"],
            colors=[(0, 0, 0), (0.5, 0.5, 0.5), (1, 1, 1)],
        )

        mask = np.ones((x.shape[0]), dtype=np.bool)
        mask[5:] = False
        y = l[mask, 0]

        np.testing.assert_array_equal(x[mask, :][:, [0]], np.array(y))

    def test_comb_row_mask_col_mask(self):
        x = np.random.random((10, 3))
        l = Lineage(
            x,
            names=["foo", "bar", "baz"],
            colors=[(0, 0, 0), (0.5, 0.5, 0.5), (1, 1, 1)],
        )

        row_mask = np.ones((x.shape[0]), dtype=np.bool)
        row_mask[5:] = False
        col_mask = np.ones((x.shape[1]), dtype=np.bool)
        y = l[row_mask, col_mask]

        np.testing.assert_array_equal(x[row_mask, :][:, col_mask], np.array(y))

    def test_comb_row_mask_col_names(self):
        x = np.random.random((10, 3))
        l = Lineage(
            x,
            names=["foo", "bar", "baz"],
            colors=[(0, 0, 0), (0.5, 0.5, 0.5), (1, 1, 1)],
        )

        mask = np.ones((x.shape[0]), dtype=np.bool)
        mask[5:] = False
        y = l[mask, ["baz", "bar"]]

        np.testing.assert_array_equal(x[mask, :][:, [2, 1]], np.array(y))

    def test_reordering(self):
        x = np.random.random((10, 3))
        l = Lineage(
            x, names=["foo", "bar", "baz"], colors=["#ff0000", "#00ff00", "#0000ff"]
        )

        y = l[["baz", "bar", "foo"]]

        np.testing.assert_array_equal(y.names, ["baz", "bar", "foo"])
        np.testing.assert_array_equal(y.colors, ["#0000ff", "#00ff00", "#ff0000"])

    def test_non_trivial_subset(self):
        x = np.random.random((10, 3))
        l = Lineage(
            x, names=["foo", "bar", "baz"], colors=["#ff0000", "#00ff00", "#0000ff"]
        )

        mask = np.ones((x.shape[0]), dtype=np.bool)
        mask[5:] = False
        y = l[mask, :][:, ["baz", "bar", "foo"]]

        np.testing.assert_array_equal(y, x[mask, :][:, ::-1])
        np.testing.assert_array_equal(y.names, ["baz", "bar", "foo"])
        np.testing.assert_array_equal(y.colors, ["#0000ff", "#00ff00", "#ff0000"])

    def test_non_trivial_subset_2(self):
        x = np.random.random((10, 3))
        l = Lineage(
            x, names=["foo", "bar", "baz"], colors=["#ff0000", "#00ff00", "#0000ff"]
        )

        mask = np.ones((x.shape[0]), dtype=np.bool)
        mask[5:] = False
        y = l[mask, ["baz", "bar", "foo"]]
        z = l[mask, :][:, ["baz", "bar", "foo"]]

        np.testing.assert_array_equal(y, z)
        np.testing.assert_array_equal(y, x[mask, :][:, [2, 1, 0]])
        np.testing.assert_array_equal(y.names, z.names)
        np.testing.assert_array_equal(y.colors, z.colors)

    def test_col_order(self):
        x = np.random.random((10, 5))
        l = Lineage(
            x,
            names=["foo", "bar", "baz", "quux", "wex"],
            colors=["#ff0000", "#00ff00", "#0000ff", "#aaaaaa", "#bbbbbb"],
        )

        y = l[["wex", "quux"]]

        np.testing.assert_array_equal(x[:, [4, 3]], y)
        np.testing.assert_array_equal(y.names, ["wex", "quux"])
        np.testing.assert_array_equal(y.colors, ["#bbbbbb", "#aaaaaa"])

    def test_automatic_color_assignment(self):
        x = np.random.random((10, 3))
        l = Lineage(x, names=["foo", "bar", "baz"])

        gt_colors = [colors.to_hex(c) for c in _create_categorical_colors(3)]

        np.testing.assert_array_equal(l.colors, gt_colors)

    def test_correct_names_to_ixs(self):
        x = np.random.random((10, 3))
        l = Lineage(x, names=["foo", "bar", "baz"])

        y = l[["baz", "bar"]]

        assert y._names_to_ixs == {"baz": 0, "bar": 1}

    def test_correct_order(self):
        x = np.random.random((10, 3))
        l = Lineage(x, names=["foo", "bar", "baz"])

        np.testing.assert_array_equal(l[["foo", "baz"]].X, l[["baz", "foo"]].X[:, ::-1])

    def test_mask_x_full_names_y(self):
        x = np.random.random((10, 3))
        l = Lineage(x, names=["Beta", "Epsilon", "Alpha"])
        cmapper = dict(zip(l.names, l.colors))
        mask = np.zeros(l.shape[0], dtype=np.bool)
        mask[0] = True
        mask[-1] = True

        y = l[mask, ["Epsilon", "Alpha", "Beta"]]

        assert y.shape == (2, 3)
        np.testing.assert_array_equal(y.names, ["Epsilon", "Alpha", "Beta"])
        np.testing.assert_array_equal(y.colors, [cmapper[n] for n in y.names])
        np.testing.assert_array_equal(y.X, x[[0, -1], :][:, [1, 2, 0]])


class TestLineageMixing:
    def test_overlap(self):
        x = Lineage(np.random.random((10, 4)), names=["foo", "bar", "baz", "quux"])

        with pytest.raises(ValueError):
            _ = x[["foo, bar", "foo"]]

    def test_overlap_mix(self):
        x = Lineage(np.random.random((10, 4)), names=["foo", "bar", "baz", "quux"])

        with pytest.raises(ValueError):
            _ = x[["foo, bar", 0]]

    def test_ellipsis_and_none(self):
        x = Lineage(np.random.random((10, 4)), names=["foo", "bar", "baz", "quux"])

        with pytest.raises(ValueError):
            _ = x[["foo, bar", Lin.REST, Lin.REST]]

    def test_rest(self):
        x = Lineage(np.random.random((10, 4)), names=["foo", "bar", "baz", "quux"])
        y = x[["foo, bar", Lin.REST]]

        expected = np.c_[np.sum(x.X[:, [0, 1]], axis=1), np.sum(x.X[:, [2, 3]], axis=1)]

        assert y.shape == (10, 2)
        np.testing.assert_array_equal(y.X, expected)
        np.testing.assert_array_equal(y.names, ["bar or foo", "rest"])
        np.testing.assert_array_equal(
            y.colors,
            [_compute_mean_color(x.colors[:2]), _compute_mean_color(x.colors[2:])],
        )

    def test_no_mixing(self):
        x = Lineage(np.random.random((10, 4)), names=["foo", "bar", "baz", "quux"])
        y = x[["foo", Lin.REST]]

        expected = np.c_[np.sum(x.X[:, [0]], axis=1), np.sum(x.X[:, [1, 2, 3]], axis=1)]

        assert y.shape == (10, 2)
        np.testing.assert_array_equal(y.X, expected)
        np.testing.assert_array_equal(y.names, ["foo", "rest"])
        np.testing.assert_array_equal(
            y.colors, [x.colors[0], _compute_mean_color(x.colors[1:])]
        )

    def test_no_rest_or_none(self):
        x = Lineage(np.random.random((10, 4)), names=["foo", "bar", "baz", "quux"])
        y = x[["foo, bar"]]

        expected = np.sum(x.X[:, [0, 1]], axis=1)[..., np.newaxis]

        assert y.shape == (10, 1)
        np.testing.assert_array_equal(y.X, expected)
        np.testing.assert_array_equal(y.names, ["bar or foo"])
        np.testing.assert_array_equal(y.colors, [_compute_mean_color(x.colors[:2])])

    def test_rest_all(self):
        x = Lineage(np.random.random((10, 4)), names=["foo", "bar", "baz", "quux"])
        y = x[[Lin.REST]]

        assert y.shape == (10, 1)
        np.testing.assert_array_equal(y.X[:, 0], np.sum(x.X, axis=1))
        np.testing.assert_array_equal(y.names, ["rest"])

    def test_row_subset(self):
        x = Lineage(np.random.random((10, 4)), names=["foo", "bar", "baz", "quux"])
        y = x[:5, ["foo, bar"]]

        expected = np.sum(x.X[:5, [0, 1]], axis=1)[..., np.newaxis]

        assert y.shape == (5, 1)
        np.testing.assert_array_equal(y.X, expected)
        np.testing.assert_array_equal(y.names, ["bar or foo"])
        np.testing.assert_array_equal(y.colors, [_compute_mean_color(x.colors[:2])])

    def test_rest_no_effect(self):
        names = ["foo", "bar", "baz", "quux"]
        x = Lineage(np.random.random((10, 4)), names=names)

        y = x[names + [Lin.REST]]

        np.testing.assert_array_equal(x.X, y.X)
        np.testing.assert_array_equal(x.names, y.names)
        np.testing.assert_array_equal(x.colors, y.colors)

    def test_others_all(self):
        names = ["foo", "bar", "baz", "quux"]
        x = Lineage(np.random.random((10, 4)), names=names)

        y = x[[Lin.OTHERS]]

        np.testing.assert_array_equal(x.X, y.X)
        np.testing.assert_array_equal(x.names, y.names)
        np.testing.assert_array_equal(x.colors, y.colors)

    def test_others(self):
        names = ["foo", "bar", "baz", "quux"]
        x = Lineage(np.random.random((10, 4)), names=names)

        y = x[["foo, baz"] + [Lin.OTHERS]]
        expected = x[["foo, baz", "bar", "quux"]]

        np.testing.assert_array_equal(y.X, expected.X)
        np.testing.assert_array_equal(y.names, expected.names)
        np.testing.assert_array_equal(y.colors, expected.colors)

    def test_others_no_effect(self):
        names = ["foo", "bar", "baz", "quux"]
        x = Lineage(np.random.random((10, 4)), names=names)

        y = x[names + [Lin.OTHERS]]

        np.testing.assert_array_equal(x.X, y.X)
        np.testing.assert_array_equal(x.names, y.names)
        np.testing.assert_array_equal(x.colors, y.colors)


class TestLineageNormalization:
    def test_empty_keys(self, lineage: Lineage):
        with pytest.raises(ValueError):
            lineage.reduce([])

    def test_not_summing_to_1(self, lineage: Lineage):
        lineage[0, 0] = 0
        with pytest.raises(ValueError):
            lineage.reduce(["foo"])

    def test_invalid_key(self, lineage: Lineage):
        with pytest.raises(ValueError):
            lineage.reduce(keys=["non_existent"])

    def test_all_names(self, lineage: Lineage):
        lin = lineage.reduce(lineage.names)
        np.testing.assert_array_equal(lin.X, lineage.X)

    def test_invalid_mode(self, lineage: Lineage):
        with pytest.raises(ValueError):
            lineage.reduce(["foo", "bar"], mode="foo")

    def test_invalid_dist_measure(self, lineage: Lineage):
        with pytest.raises(ValueError):
            lineage.reduce(["foo", "bar"], dist_measure="foo")

    def test_invalid_weight_normalize(self, lineage: Lineage):
        with pytest.raises(ValueError):
            lineage.reduce(["foo", "bar"], normalize_weights="foo")

    def test_return_weights_mode_scale(self, lineage: Lineage):
        lin, weights = lineage.reduce(["foo", "bar"], mode="scale", return_weights=True)

        assert isinstance(lin, Lineage)
        assert weights is None

    def test_return_weights_mode_dist(self, lineage: Lineage):
        lin, weights = lineage.reduce(["foo", "bar"], mode="dist", return_weights=True)

        assert isinstance(lin, Lineage)
        assert isinstance(weights, DataFrame)

    def test_normal_only_1(self, lineage: Lineage):
        lin = lineage.reduce("foo")

        assert lin.shape == (10, 1)
        np.testing.assert_allclose(np.sum(lin, axis=1), 1.0)
        np.testing.assert_array_equal(lin.names, ["foo"])
        np.testing.assert_array_equal(lin.colors, lineage[["foo"]].colors)

    def test_normal_run(self, lineage: Lineage):
        lin = lineage.reduce(["foo", "bar"])

        assert lin.shape == (10, 2)
        np.testing.assert_allclose(np.sum(lin, axis=1), 1.0)
        np.testing.assert_array_equal(lin.names, ["foo", "bar"])
        np.testing.assert_array_equal(lin.colors, lineage[["foo", "bar"]].colors)

    @mock.patch("cellrank.tools._lineage._cosine_sim")
    def test_cosine(self, mocker, lineage: Lineage):
        try:
            _ = lineage.reduce(["foo", "bar"], dist_measure="cosine_sim", mode="dist")
        except ValueError:
            pass
        finally:
            mocker.assert_called_once()

    @mock.patch("cellrank.tools._lineage._wasserstein_dist")
    def test_wasserstein(self, mocker, lineage: Lineage):
        try:
            _ = lineage.reduce(
                ["foo", "bar"], dist_measure="wasserstein_dist", mode="dist"
            )
        except ValueError:
            pass
        finally:
            mocker.assert_called_once()

    @mock.patch("cellrank.tools._lineage._kl_div")
    def test_kl_div(self, mocker, lineage: Lineage):
        try:
            _ = lineage.reduce(["foo", "bar"], dist_measure="kl_div", mode="dist")
        except ValueError:
            pass
        finally:
            mocker.assert_called_once()

    @mock.patch("cellrank.tools._lineage._js_div")
    def test_js_div(self, mocker, lineage: Lineage):
        try:
            _ = lineage.reduce(["foo", "bar"], dist_measure="js_div", mode="dist")
        except ValueError:
            pass
        finally:
            mocker.assert_called_once()

    @mock.patch("cellrank.tools._lineage._mutual_info")
    def test_mutual_info(self, mocker, lineage: Lineage):
        try:
            _ = lineage.reduce(["foo", "bar"], dist_measure="mutual_info", mode="dist")
        except ValueError:
            pass
        finally:
            mocker.assert_called_once()

    @mock.patch("cellrank.tools._lineage._row_normalize")
    def test_equal(self, mocker, lineage: Lineage):
        try:
            _ = lineage.reduce(["foo", "bar"], dist_measure="equal", mode="dist")
        except ValueError:
            pass
        finally:
            # should be twice, but we have extra check inside and we're mocking that does nothing
            mocker.assert_called_once()

    @mock.patch("cellrank.tools._lineage._row_normalize")
    def test_row_normalize(self, mocker, lineage: Lineage):
        try:
            _ = lineage.reduce(["foo", "bar"], mode="scale")
        except ValueError:
            pass
        finally:
            mocker.assert_called_once()

    @mock.patch("cellrank.tools._lineage._softmax")
    def test_softmax(self, mocker, lineage: Lineage):
        try:
            _ = lineage.reduce(["foo", "bar"], normalize_weights="softmax", mode="dist")
        except ValueError:
            pass
        finally:
            mocker.assert_called_once()


class TestLineageSameLengthIndexing:
    def test_same_names(self):
        x = Lineage(np.random.random((10, 4)), names=["foo", "bar", "baz", "quux"])
        y = x[np.arange(len(x)), ["foo"] * len(x)]

        expected = x["foo"]

        assert y.shape == (10, 1)
        np.testing.assert_array_equal(y.X, expected.X)
        np.testing.assert_array_equal(y.names, ["mixture"])
        np.testing.assert_array_equal(y.colors, ["#000000"])

    def test_same_indices(self):
        x = Lineage(np.random.random((10, 4)), names=["foo", "bar", "baz", "quux"])
        half = len(x) // 2
        y = x[[0] * len(x), ["foo"] * half + ["bar"] * half]

        expected = np.array([x[0, "foo"].X[0, 0]] * half + [x[0, "bar"].X[0, 0]] * half)

        assert y.shape == (10, 1)
        np.testing.assert_array_equal(y.X.squeeze(), expected)
        np.testing.assert_array_equal(y.names, ["mixture"])
        np.testing.assert_array_equal(y.colors, ["#000000"])
