# -*- coding: utf-8 -*-
import pytest

import numpy as np
import pandas as pd

from matplotlib.colors import is_color_like

from cellrank.tl._colors import _map_names_and_colors, _create_categorical_colors


class TestColors:
    def test_create_categorical_colors_too_many_colors(self):
        with pytest.raises(ValueError):
            _create_categorical_colors(1000)

    def test_create_categorical_colors_no_categories(self):
        c = _create_categorical_colors(0)

        assert c == []

    def test_create_categorical_colors_neg_categories(self):
        with pytest.raises(RuntimeError):
            _create_categorical_colors(-1)

    def test_create_categorical_colors_normal_run(self):
        colors = _create_categorical_colors(79)

        assert len(colors) == 79
        assert all(map(lambda c: isinstance(c, str), colors))
        assert all(map(lambda c: is_color_like(c), colors))

    def test_mapping_colors_not_categorical(self):
        query = pd.Series(["foo", "bar", "baz"], dtype="str")
        reference = pd.Series(["foo", np.nan, "bar", "baz"], dtype="category")

        with pytest.raises(TypeError):
            _map_names_and_colors(reference, query)

    def test_mapping_colors_invalid_size(self):
        query = pd.Series(["foo", "bar", "baz"], dtype="category")
        reference = pd.Series(["foo", np.nan, "bar", "baz"], dtype="category")

        with pytest.raises(ValueError):
            _map_names_and_colors(reference, query)

    def test_mapping_colors_different_index(self):
        query = pd.Series(["foo", "bar", "baz"], dtype="category", index=[2, 3, 4])
        reference = pd.Series(["foo", "bar", "baz"], dtype="category", index=[1, 2, 3])

        with pytest.raises(ValueError):
            _map_names_and_colors(reference, query)
