# -*- coding: utf-8 -*-
import pytest

from matplotlib.colors import is_color_like

from cellrank.tl._colors import _create_categorical_colors


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
