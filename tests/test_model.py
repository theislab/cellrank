# -*- coding: utf-8 -*-
import pytest
from _helpers import create_model

from anndata import AnnData

import numpy as np
from sklearn.svm._classes import SVR

from cellrank.ul.models._utils import _extract_data


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
