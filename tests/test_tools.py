# -*- coding: utf-8 -*-
import cellrank as cr
import pandas as pd
import pytest

from anndata import AnnData
from _helpers import create_model


class TestGeneImportance:
    def test_invalid_n_perms(self, adata: AnnData):
        model = create_model(adata)
        with pytest.raises(ValueError):
            _ = cr.tl.gene_importance(
                adata, model, adata.var_names[:10], "0", n_perms=-1
            )

    def test_no_perms(self, adata: AnnData):
        model = create_model(adata)
        res = cr.tl.gene_importance(
            adata, model, adata.var_names[:10], "0", n_perms=None
        )

        assert isinstance(res, pd.DataFrame)
        assert res.shape == (10, 1)
        assert set(res.index) == set(adata.var_names[:10])
        assert "imporance" in res.columns

    def test_perms_no_fdr_correction(self, adata: AnnData):
        model = create_model(adata)
        res = cr.tl.gene_importance(
            adata, model, adata.var_names, "0", n_perms=50, fdr_correction=None
        )

        assert "pval" in res.columns
        assert "qval" not in res.columns

    def test_perms_fdr_correction(self, adata: AnnData):
        model = create_model(adata)
        res = cr.tl.gene_importance(adata, model, adata.var_names, "0", n_perms=50)

        assert "pval" in res.columns
        assert "qval" in res.columns


class TestLineage:
    pass


class TestClusterFates:
    pass


class TestRootFinal:
    pass


class TestTransitionMatrix:
    pass
