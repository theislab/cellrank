# -*- coding: utf-8 -*-
import pytest
from _helpers import create_model

from anndata import AnnData

import numpy as np
from scipy.stats import rankdata
from sklearn.svm._classes import SVR

from cellrank.ul.models._utils import _rankdata, _extract_data, _get_knotlocs


class TestModel:
    def test_initialize(self, adata: AnnData):
        model = create_model(adata)

        assert isinstance(model.model, SVR)

    def test_prepare_invalid_gene(self, adata_cflare):
        model = create_model(adata_cflare)
        with pytest.raises(KeyError):
            model.prepare("foo", "0")

    def test_prepare_invalid_lineage(self, adata_cflare):
        model = create_model(adata_cflare)
        with pytest.raises(KeyError):
            model.prepare(adata_cflare.var_names[0], "foo")

    def test_prepare_invalid_data_key(self, adata_cflare):
        model = create_model(adata_cflare)
        with pytest.raises(KeyError):
            model.prepare(adata_cflare.var_names[0], "0", data_key="foo")

    def test_prepare_invalid_time_key(self, adata_cflare):
        model = create_model(adata_cflare)
        with pytest.raises(KeyError):
            model.prepare(adata_cflare.var_names[0], "0", time_key="foo")

    def test_prepare_invalid_time_range(self, adata_cflare):
        model = create_model(adata_cflare)
        with pytest.raises(ValueError):
            model.prepare(adata_cflare.var_names[0], "0", time_range=(0, 1, 2))

    def test_prepare_normal_run(self, adata_cflare):
        model = create_model(adata_cflare)
        model = model.prepare(adata_cflare.var_names[0], "0")

        assert isinstance(model.x, np.ndarray)
        assert isinstance(model.w, np.ndarray)
        assert isinstance(model.y, np.ndarray)
        assert isinstance(model.x_test, np.ndarray)
        assert len(model.x_test) == 200
        assert model.y_test is None
        assert model.conf_int is None

    def test_prepare_n_test_points(self, adata_cflare):
        model = create_model(adata_cflare)
        model = model.prepare(adata_cflare.var_names[0], "0", n_test_points=300)

        assert len(model.x_test) == 300

    def test_predict(self, adata_cflare):
        model = create_model(adata_cflare)
        model = model.prepare(adata_cflare.var_names[0], "0").fit()
        y_hat = model.predict()

        assert isinstance(model.y_test, np.ndarray)
        assert len(model.x_test) == len(model.y_test)
        assert y_hat is model.y_test

        assert model.conf_int is None

    def test_confidence_interval(self, adata_cflare):
        model = create_model(adata_cflare)
        model = model.prepare(adata_cflare.var_names[0], "0").fit()
        _ = model.predict()
        ci = model.confidence_interval()

        assert isinstance(model.conf_int, np.ndarray)
        assert len(model.y_test) == len(model.conf_int)
        assert ci is model.conf_int


class TestUtils:
    def test_extract_data_wrong_type(self):
        with pytest.raises(TypeError):
            _extract_data(None)

    def test_extract_data_raw_None(self, adata: AnnData):
        adata = AnnData(adata.X, raw=None)
        with pytest.raises(ValueError):
            _extract_data(adata, use_raw=True)

    def test_extract_data_invalid_layer(self, adata: AnnData):
        with pytest.raises(KeyError):
            _extract_data(adata, layer="foo", use_raw=False)

    def test_extract_data_normal_run(self, adata: AnnData):
        X = _extract_data(adata, use_raw=False)

        assert X is adata.X

    def test_extract_data_normal_run_layer(self, adata: AnnData):
        ms = _extract_data(adata, layer="Ms", use_raw=False)

        assert ms is adata.layers["Ms"]

    def test_extract_data_normal_run_raw(self, adata: AnnData):
        raw = _extract_data(adata, use_raw=True, layer="Ms")

        assert raw is adata.raw.X

    @pytest.mark.parametrize("method", ["average", "min", "max", "dense", "ordinal"])
    def test_rank_data(self, method):
        x = np.random.normal(size=(10,))

        np.testing.assert_array_equal(_rankdata(x), rankdata(x))

    def test_rank_data_invalid_method(self):
        with pytest.raises(AssertionError):
            _rankdata(np.random.normal(size=(10,)), method="foobar")

    def test_get_knots_invalid_n_knots(self):
        with pytest.raises(ValueError):
            _get_knotlocs([0, 1, 2], 0)

    def test_get_knots_non_finite_values(self):
        x = np.array([0, 1, 2, 3], dtype=np.float64)
        x[-1] = np.inf
        with pytest.raises(ValueError):
            _get_knotlocs(x, 1)

    def test_get_knots_wrong_shape(self):
        with pytest.raises(ValueError):
            _get_knotlocs(np.array([0, 1, 2, 3]).reshape((2, 2)), 1)

    def test_get_knots_only_same_value(self):
        with pytest.raises(ValueError):
            _get_knotlocs(np.array([42] * 10), 1)

    def test_get_knots_empty_pseudotime(self):
        with pytest.raises(ValueError):
            _get_knotlocs(np.array([]), 2)

    def test_get_knots_uniform(self):
        expected = np.linspace(0, 5, 3, endpoint=True)
        actual = _get_knotlocs(np.array([3, 5, 4, 0]), 3, uniform=True)

        np.testing.assert_array_equal(actual, expected)

    def test_get_knots_uniform_1_knot(self):
        actual = _get_knotlocs(np.array([3, 5, 4, 0]), 1, uniform=True)

        np.testing.assert_array_equal(actual, [5])

    def test_get_knots_1_knot(self):
        actual = _get_knotlocs(np.array([3, 5, 4, 0]), 1, uniform=False)

        np.testing.assert_array_equal(actual, [5])

    def test_get_knots_2d(self):
        expected = np.linspace(0, 5, 3, endpoint=True)
        actual = _get_knotlocs(np.array([3, 5, 4, 0]).reshape((-1, 1)), 3, uniform=True)

        assert actual.ndim == 1
        np.testing.assert_array_equal(actual, expected)

    @pytest.mark.parametrize("seed,n_knots", zip(range(10), range(2, 11)))
    def test_get_knots_unique(self, seed, n_knots):
        np.random.seed(seed)
        x = np.random.normal(size=(100,))
        actual = _get_knotlocs(x, n_knots=n_knots)

        assert actual.shape == (n_knots,)
        np.testing.assert_array_equal(actual, np.sort(actual))
        assert len(np.unique(actual)) == len(actual), actual
        assert (np.min(actual), np.max(actual)) == (np.min(x), np.max(x))

    def test_get_knots_heavy_tail(self):
        x = np.array([0] * 30 + list(np.linspace(0.1, 0.9, 30)) + [1] * 30)
        expected = np.array(
            [
                0.0,
                0.02222222,
                0.04444444,
                0.06666667,
                0.36360153,
                0.63639847,
                0.93333333,
                0.95555556,
                0.97777778,
                1.0,
            ]
        )

        actual = _get_knotlocs(x, 10, uniform=False)

        np.testing.assert_almost_equal(actual, expected)
