# -*- coding: utf-8 -*-
import pytest
from _helpers import create_model, assert_array_nan_equal, jax_not_installed_skip

import scanpy as sc
from anndata import AnnData

import numpy as np
import pandas as pd
from scipy.sparse import diags, random, csr_matrix
from pandas.api.types import is_categorical_dtype

from cellrank.tl import Lineage
from cellrank.pl._utils import _create_models, _create_callbacks
from cellrank.tl._utils import (
    _one_hot,
    _cluster_X,
    _process_series,
    _fuzzy_to_discrete,
    _merge_categorical_series,
    _series_from_one_hot_matrix,
)
from cellrank.ul.models import GAM
from cellrank.tl._colors import _compute_mean_color
from cellrank.tl.kernels._utils import (
    norm,
    np_max,
    np_sum,
    np_mean,
    _reconstruct_one,
    _calculate_starts,
    _get_probs_for_zero_vec,
    _predict_transition_probabilities_jax,
    _predict_transition_probabilities_numpy,
)


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

        np.testing.assert_array_equal(colors_merged, ["red", "blue", "#279e68"])

    def test_merge_colors_simple_old_no_inplace(self):
        x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")
        y = pd.Series(["b", np.nan, "a", "d", "a"]).astype("category")
        expected = pd.Series(["b", "b", "a", "d", "a"]).astype("category")
        colors_x = ["red", "blue"]

        merged, colors_merged = _merge_categorical_series(
            x, y, colors_old=colors_x, inplace=False
        )

        np.testing.assert_array_equal(merged.values, expected.values)
        np.testing.assert_array_equal(colors_merged, ["red", "blue", "#279e68"])

    def test_merge_colors_simple_new(self):
        x = pd.Series(["a", "b", np.nan, "b", np.nan]).astype("category")
        y = pd.Series(["b", np.nan, "a", "d", "a"]).astype("category")
        colors_y = ["red", "blue", "green"]

        colors_merged = _merge_categorical_series(
            x, y, colors_new=colors_y, inplace=True
        )

        np.testing.assert_array_equal(colors_merged, ["#1f77b4", "#ff7f0e", "green"])

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
        assert set(colors) == {
            _compute_mean_color(["red", "green"]),
            _compute_mean_color(["blue", "white"]),
        }


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


class TestPlottingUtils:
    def test_create_models_not_a_model_local(self):
        with pytest.raises(TypeError):
            _create_models({"foo": {"bar": 42}}, ["foo"], ["bar"])

    def test_create_models_not_a_model_gene_fallback(self, adata: AnnData):
        m = create_model(adata)
        with pytest.raises(TypeError):
            _create_models({"foo": {"baz": m}, "*": 42}, ["foo", "bar"], ["baz"])

    def test_create_models_not_a_model_lineage_fallback(self, adata: AnnData):
        m = create_model(adata)
        with pytest.raises(TypeError):
            _create_models({"foo": {"baz": m, "*": 42}}, ["foo"], ["bar", "baz"])

    def test_create_models_not_a_model_global(self):
        with pytest.raises(TypeError):
            _create_models(None, ["foo"], ["bar"])

    def test_create_models_no_models_gene(self):
        with pytest.raises(ValueError):
            _create_models({}, ["foo"], [])

    def test_create_models_no_models_lineage(self):
        with pytest.raises(ValueError):
            _create_models({"foo": {}}, ["foo"], ["bar"])

    def test_create_models_gene_incomplete(self, adata: AnnData):
        m = create_model(adata)
        with pytest.raises(ValueError):
            _create_models({"foo": {"baz": m}}, ["foo", "bar"], ["baz"])

    def test_create_models_lineage_incomplete(self, adata: AnnData):
        m = create_model(adata)
        with pytest.raises(ValueError):
            _create_models({"foo": {"baz": m}}, ["foo"], ["bar", "baz"])

    def test_create_model_no_genes(self, adata: AnnData):
        with pytest.raises(ValueError):
            m = create_model(adata)
            _create_models(m, [], ["foo"])

    def test_create_model_no_lineage(self, adata: AnnData):
        with pytest.raises(ValueError):
            m = create_model(adata)
            _create_models(m, ["foo"], [])

    def test_create_models_1_model(self, adata: AnnData):
        m = create_model(adata)
        models = _create_models(m, ["foo"], ["bar"])

        assert set(models.keys()) == {"foo"}
        assert set(models["foo"].keys()) == {"bar"}
        assert isinstance(models["foo"]["bar"], type(m))
        assert models["foo"]["bar"] is not m

    def test_create_models_gene_specific(self, adata: AnnData):
        m1 = create_model(adata)
        m2 = GAM(adata)

        models = _create_models({"foo": m1, "bar": m2}, ["foo", "bar"], ["baz"])
        assert set(models.keys()) == {"foo", "bar"}
        assert set(models["foo"].keys()) == {"baz"}
        assert set(models["bar"].keys()) == {"baz"}
        assert isinstance(models["foo"]["baz"], type(m1))
        assert models["foo"]["baz"] is not m1

        assert isinstance(models["bar"]["baz"], type(m2))
        assert models["bar"]["baz"] is not m2

    def test_create_models_gene_specific_fallback(self, adata: AnnData):
        m1 = create_model(adata)
        m2 = GAM(adata)

        models = _create_models(
            {"foo": m1, "*": m2}, ["foo", "bar", "baz", "quux"], ["quas", "wex"]
        )
        assert set(models.keys()) == {"foo", "bar", "baz", "quux"}
        for k, vs in models.items():
            assert set(vs.keys()) == {"quas", "wex"}

        for g in {"foo"}:
            for l in {"quas", "wex"}:
                assert isinstance(models[g][l], type(m1))
                assert models[g][l] is not m1

        for g in {"bar", "baz", "quux"}:
            for l in {"quas", "wex"}:
                assert isinstance(models[g][l], type(m2))
                assert models[g][l] is not m2

    def test_create_models_lineage_specific(self, adata: AnnData):
        m1 = create_model(adata)
        m2 = GAM(adata)

        models = _create_models(
            {"foo": {"bar": m1, "baz": m2}}, ["foo"], ["bar", "baz"]
        )
        assert set(models["foo"].keys()) == {"bar", "baz"}
        assert isinstance(models["foo"]["bar"], type(m1))
        assert models["foo"]["bar"] is not m1

        assert isinstance(models["foo"]["baz"], type(m2))
        assert models["foo"]["baz"] is not m2

    def test_create_models_lineage_specific_fallback(self, adata: AnnData):
        m1 = create_model(adata)
        m2 = GAM(adata)

        models = _create_models(
            {"foo": {"baz": m1, "*": m2}, "bar": {"quux": m2, "*": m1}},
            ["foo", "bar"],
            ["baz", "quux", "quas", "wex"],
        )
        assert set(models.keys()) == {"foo", "bar"}

        for k, vs in models.items():
            assert set(vs.keys()) == {"baz", "quux", "quas", "wex"}

        assert isinstance(models["foo"]["baz"], type(m1))
        assert models["foo"]["baz"] is not m1

        for l in {"quux", "quas", "wex"}:
            assert isinstance(models["foo"][l], type(m2))
            assert models["foo"][l] is not m2

        assert isinstance(models["bar"]["quux"], type(m2))
        assert models["bar"]["quux"] is not m2

        for l in {"baz", "quas", "wex"}:
            assert isinstance(models["bar"][l], type(m1))
            assert models["bar"][l] is not m1

    # TODO: test callback


class TestClusterX:
    def test_normal_run_louvain(self):
        # create some data
        adata = sc.datasets.blobs(n_observations=100, n_variables=6)

        # kmeans, louvain, leiden
        labels_kmeans = _cluster_X(adata.X, n_clusters=5, method="kmeans")
        labels_louvain = _cluster_X(adata.X, n_clusters=5, method="louvain")

        assert len(labels_kmeans) == len(labels_louvain) == adata.n_obs

    def test_normal_run_leiden(self):
        # create some data
        adata = sc.datasets.blobs(n_observations=100, n_variables=6)

        # kmeans, louvain, leiden
        labels_kmeans = _cluster_X(adata.X, n_clusters=5, method="kmeans")
        labels_louvain = _cluster_X(adata.X, n_clusters=5, method="leiden")

        assert len(labels_kmeans) == len(labels_louvain) == adata.n_obs

    def test_one_feature(self):
        # create some data
        adata = sc.datasets.blobs(n_observations=100, n_variables=1)

        # kmeans, louvain, leiden
        labels_kmeans = _cluster_X(adata.X, n_clusters=5, method="kmeans")
        labels_louvain = _cluster_X(adata.X, n_clusters=5, method="louvain")

        assert len(labels_kmeans) == len(labels_louvain) == adata.n_obs


class TestKernelUtils:
    @pytest.mark.parametrize("seed", range(10))
    def test_numba_mean_axis_0(self, seed):
        np.random.seed(seed)
        x = np.random.normal(size=(10, 10))

        np.testing.assert_allclose(np.mean(x, axis=0), np_mean(x, 0))

    @pytest.mark.parametrize("seed", range(10))
    def test_numba_mean_axis_1(self, seed):
        np.random.seed(seed)
        x = np.random.normal(size=(10, 10))

        np.testing.assert_allclose(np.mean(x, axis=1), np_mean(x, 1))

    @pytest.mark.parametrize("seed", range(10))
    def test_numba_max_axis_0(self, seed):
        np.random.seed(seed)
        x = np.random.normal(size=(10, 10))

        np.testing.assert_allclose(np.max(x, axis=0), np_max(x, 0))

    @pytest.mark.parametrize("seed", range(10))
    def test_numba_max_axis_1(self, seed):
        np.random.seed(seed)
        x = np.random.normal(size=(10, 10))

        np.testing.assert_allclose(np.max(x, axis=1), np_max(x, 1))

    @pytest.mark.parametrize("seed", range(10))
    def test_numba_sum_axis_0(self, seed):
        np.random.seed(seed)
        x = np.random.normal(size=(10, 10))

        np.testing.assert_allclose(np.sum(x, axis=0), np_sum(x, 0))

    @pytest.mark.parametrize("seed", range(10))
    def test_numba_sum_axis_1(self, seed):
        np.random.seed(seed)
        x = np.random.normal(size=(10, 10))

        np.testing.assert_allclose(np.sum(x, axis=1), np_sum(x, 1))

    @pytest.mark.parametrize("seed", range(10))
    def test_numba_norm_axis_0(self, seed):
        np.random.seed(seed)
        x = np.random.normal(size=(10, 10))

        np.testing.assert_allclose(np.linalg.norm(x, axis=0), norm(x, 0))

    @pytest.mark.parametrize("seed", range(10))
    def test_numba_norm_axis_1(self, seed):
        np.random.seed(seed)
        x = np.random.normal(size=(10, 10))

        np.testing.assert_allclose(np.linalg.norm(x, axis=1), norm(x, 1))

    def test_zero_unif_sum_to_1_vector(self):
        sum_to_1, zero = _get_probs_for_zero_vec(10)

        assert zero.shape == (10,)
        assert sum_to_1.shape == (10,)
        np.testing.assert_array_equal(zero, np.zeros_like(zero))
        np.testing.assert_allclose(sum_to_1, np.ones_like(sum_to_1) / sum_to_1.shape[0])
        assert np.isclose(sum_to_1.sum(), 1.0)

    def calculate_starts(self):
        starts = _calculate_starts(diags(np.ones(10)).tocsr().indptr, np.arange(10))

        np.testing.assert_array_equal(starts, np.arange(10))

    @pytest.mark.parametrize("seed, shuffle", zip(range(10), [False] * 5 + [True] * 5))
    def test_reconstruct_one(self, seed: int, shuffle: bool):
        m1 = random(100, 10, random_state=seed, density=0.5, format="csr")
        m1[:, 0] = 0.1
        m1 /= m1.sum(1)
        m1 = csr_matrix(m1)

        m2_data = np.random.normal(size=(m1.nnz))
        m2 = csr_matrix((m2_data, m1.indices, m1.indptr))

        if shuffle:
            ixs = np.arange(100)
            np.random.shuffle(ixs)
            data = np.c_[m1[ixs, :].data, m2[ixs, :].data].T
        else:
            ixs = None
            data = np.c_[m1.data, m2.data].T

        r1, r2 = _reconstruct_one(data, m1, ixs=ixs)

        np.testing.assert_array_equal(r1.A, m1.A)
        np.testing.assert_array_equal(r2.A, m2.A)

    @jax_not_installed_skip
    @pytest.mark.parametrize("seed", range(10))
    def test_numpy_and_jax(self, seed: int):
        np.random.seed(seed)
        x = np.random.normal(size=(100,))
        w = np.random.normal(size=(1, 100))

        np_res, _ = _predict_transition_probabilities_numpy(x[None, :], w, 1)
        jax_res = _predict_transition_probabilities_jax(x, w, 1)

        np.testing.assert_allclose(np_res, jax_res)
