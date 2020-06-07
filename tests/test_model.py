# -*- coding: utf-8 -*-
import numpy as np
import pytest
from sklearn.svm._classes import SVR

from anndata import AnnData

from _helpers import create_model


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

    def test_prepare_invalid_start_cluster(self, adata_cflare):
        model = create_model(adata_cflare)
        with pytest.raises(KeyError):
            model.prepare(adata_cflare.var_names[0], "0", start_lineage="foo")

    def test_prepare_invalid_end_cluster(self, adata_cflare):
        model = create_model(adata_cflare)
        with pytest.raises(KeyError):
            model.prepare(adata_cflare.var_names[0], "0", end_lineage="foo")

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
