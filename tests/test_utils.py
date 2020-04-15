# -*- coding: utf-8 -*-
from cellrank.tools._utils import _merge_approx_rcs

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
