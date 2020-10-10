# -*- coding: utf-8 -*-
import pickle
from io import BytesIO
from copy import copy, deepcopy

import pytest
from _helpers import gamr_skip, create_model, assert_models_equal

from anndata import AnnData

import numpy as np
from scipy.stats import rankdata
from sklearn.svm import SVR

from cellrank.ul.models import GAM, GAMR, SKLearnModel
from cellrank.ul.models._utils import (
    _OFFSET_KEY,
    NormMode,
    _rankdata,
    _get_offset,
    _extract_data,
    _get_knotlocs,
)
from cellrank.ul.models._base_model import FailedModel


class TestModel:
    def test_wrong_type(self):
        with pytest.raises(TypeError):
            SKLearnModel(0, SVR())

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
    def test_rank_data(self, method: str):
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
    def test_get_knots_unique(self, seed: int, n_knots: int):
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

    @pytest.mark.parametrize(
        "method,seed", zip(list(NormMode), range(len(list(NormMode))))
    )
    def test_get_offset(self, method: str, seed: int):
        np.random.seed(seed)
        x = np.random.normal(size=(100, 50))

        offset = _get_offset(x, method=method, ref_ix=0)

        assert isinstance(offset, np.ndarray)
        assert offset.shape == (100,)
        assert np.all(np.isfinite(offset))

    def test_get_offset_degenerate_case(self):
        x = np.zeros((100, 2))

        offset = _get_offset(x, ref_ix=0)

        assert isinstance(offset, np.ndarray)
        np.testing.assert_array_equal(offset, np.ones((100,)))

    def test_get_offset_writing_to_adata(self, adata: AnnData):
        offset = _get_offset(adata, use_raw=False, ref_ix=0)

        assert _OFFSET_KEY in adata.obs
        np.testing.assert_array_equal(offset, adata.obs[_OFFSET_KEY].values)

    def test_get_offset_use_raw(self, adata: AnnData):
        offset = _get_offset(adata, use_raw=False, recompute=True, ref_ix=0)
        offset_raw = _get_offset(adata, use_raw=True, recompute=True, ref_ix=0)

        assert offset.shape == offset_raw.shape == (adata.n_obs,)
        assert not np.all(np.isclose(offset, offset_raw))

    def test_offset_automatic_ref_ix(self, adata: AnnData):
        offset = _get_offset(adata, ref_ix=None)

        assert offset.shape == (adata.n_obs,)
        assert np.all(np.isfinite(offset))


@gamr_skip
class TestGAMR:
    def test_invalid_n_knots(self, adata: AnnData):
        with pytest.raises(ValueError):
            _ = GAMR(adata, n_knots=0)

    def test_invalid_smoothing_penalty(self, adata: AnnData):
        with pytest.raises(ValueError):
            _ = GAMR(adata, smoothing_penalty=-0.001)

    def test_invalid_knotlocs(self, adata: AnnData):
        with pytest.raises(ValueError):
            _ = GAMR(adata, knotlocs="foobar")

    def test_normal_initialization(self, adata_cflare: AnnData):
        m = GAMR(adata_cflare)

        assert not m.prepared
        assert m._lineage is None
        assert m._gene is None

    def test_sharing_library(self, gamr_model: GAMR):
        actual = gamr_model.copy()

        assert actual._lib_name == gamr_model._lib_name
        assert actual._lib is gamr_model._lib

    def test_shallow_copy(self, gamr_model: GAMR):
        assert_models_equal(gamr_model, copy(gamr_model), deepcopy=False)

    def test_deep_copy(self, gamr_model: GAMR):
        assert_models_equal(gamr_model, deepcopy(gamr_model), deepcopy=True)

    def test_pickling(self, gamr_model: GAMR):
        fp = BytesIO()

        pickle.dump(gamr_model, fp)
        fp.flush()
        fp.seek(0)
        actual_model = pickle.load(fp)

        assert_models_equal(gamr_model, actual_model, pickled=True)


class TestSKLearnModel:
    def test_wrong_model_type(self, adata_cflare: AnnData):
        model = create_model(adata_cflare)
        with pytest.raises(TypeError):
            SKLearnModel(adata_cflare, model)

    def test_svr_correct_no_weights(self, adata_cflare: AnnData):
        model = (
            SKLearnModel(adata_cflare, SVR(), weight_name="")
            .prepare(adata_cflare.var_names[0], "0")
            .fit()
        )
        model_w = (
            SKLearnModel(adata_cflare, SVR())
            .prepare(adata_cflare.var_names[0], "0")
            .fit()
        )

        assert model._weight_name == ""
        assert model_w._weight_name == "sample_weight"

        assert not np.allclose(model.predict(), model_w.predict())

    def test_svr_invalid_weight_name(self, adata_cflare: AnnData):
        with pytest.raises(ValueError):
            SKLearnModel(adata_cflare, SVR(), weight_name="foobar")

    def test_svr_invalid_weight_name_no_raise_fit(self, adata_cflare: AnnData):
        model = SKLearnModel(
            adata_cflare, SVR(), weight_name="w", ignore_raise=True
        ).prepare(adata_cflare.var_names[0], "0")

        with pytest.raises(TypeError):
            model.fit()

    def test_svr_invalid_weight_name_no_raise(self, adata_cflare: AnnData):
        model = SKLearnModel(
            adata_cflare, SVR(), weight_name="foobar", ignore_raise=True
        )

        assert model._weight_name == "foobar"

    def test_svr_correct_weight_name(self, adata_cflare: AnnData):
        model = SKLearnModel(adata_cflare, SVR())

        assert model._weight_name == "sample_weight"


class TestFailedModel:
    def test_wrong_model_type(self):
        with pytest.raises(TypeError):
            FailedModel(SVR())

    def test_correct_gene_and_lineage(self, gamr_model):
        fm = FailedModel(gamr_model)

        assert fm.gene == gamr_model._gene
        assert fm.lineage == gamr_model._lineage

    def test_unable_to_do_anything(self, gamr_model):
        fm = FailedModel(gamr_model)

        with pytest.raises(RuntimeError):
            fm.prepare(fm.adata.var_names[0], "0")

        for fn in [
            "fit",
            "predict",
            "confidence_interval",
            "default_confidence_interval",
            "copy",
        ]:
            with pytest.raises(RuntimeError):
                getattr(fm, fn)()


class TestModelsIO:
    def test_shallow_copy_sklearn(self, sklearn_model: SKLearnModel):
        assert_models_equal(sklearn_model, copy(sklearn_model), deepcopy=False)

    def test_deep_copy_sklearn(self, sklearn_model: SKLearnModel):
        assert_models_equal(sklearn_model, deepcopy(sklearn_model), deepcopy=True)

    def test_pickling_sklearn(self, sklearn_model: SKLearnModel):
        fp = BytesIO()

        pickle.dump(sklearn_model, fp)
        fp.flush()
        fp.seek(0)
        actual_model = pickle.load(fp)

        assert_models_equal(sklearn_model, actual_model, pickled=True)

    def test_shallow_copy_pygam(self, pygam_model: GAM):
        assert_models_equal(pygam_model, copy(pygam_model), deepcopy=False)

    def test_deep_copy_pygam(self, pygam_model: GAM):
        assert_models_equal(pygam_model, deepcopy(pygam_model), deepcopy=True)

    def test_pickling_pygam(self, pygam_model: GAM):
        fp = BytesIO()

        pickle.dump(pygam_model, fp)
        fp.flush()
        fp.seek(0)
        actual_model = pickle.load(fp)

        assert_models_equal(pygam_model, actual_model, pickled=True)
