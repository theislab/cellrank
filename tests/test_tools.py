# -*- coding: utf-8 -*-
from typing import Tuple

import cellrank as cr
import pandas as pd
import pytest

from anndata import AnnData
from scipy.spatial.distance import euclidean
from _helpers import create_model

from cellrank.tools import MarkovChain


class TestGeneImportance:
    def test_invalid_n_perms(self, adata_cr: AnnData):
        model = create_model(adata_cr)
        with pytest.raises(ValueError):
            _ = cr.tl.gene_importance(
                adata_cr, model, adata_cr.var_names[:10], "0", n_perms=-1
            )

    def test_no_perms(self, adata_cr: AnnData):
        model = create_model(adata_cr)
        res = cr.tl.gene_importance(
            adata_cr, model, adata_cr.var_names[:10], "0", n_perms=None
        )

        assert isinstance(res, pd.DataFrame)
        assert res.shape == (10, 1)
        assert set(res.index) == set(adata_cr.var_names[:10])
        assert "importance" in res.columns

    def test_perms_no_fdr_correction(self, adata_cr: AnnData):
        model = create_model(adata_cr)
        res = cr.tl.gene_importance(
            adata_cr, model, adata_cr.var_names, "0", n_perms=50, fdr_correction=None
        )

        assert "pval" in res.columns
        assert "qval" not in res.columns

    def test_perms_fdr_correction(self, adata_cr: AnnData):
        model = create_model(adata_cr)
        res = cr.tl.gene_importance(
            adata_cr, model, adata_cr.var_names, "0", n_perms=50
        )

        assert "pval" in res.columns
        assert "qval" in res.columns


class TestLineages:
    def test_no_root_cells(self, adata: AnnData):
        with pytest.raises(ValueError):
            cr.tl.lineages(adata)

    def test_normal_run(self, adata_cr: AnnData):
        cr.tl.lineages(adata_cr)

    def test_normal_run_copy(self, adata_cr: AnnData):
        adata_cr2 = cr.tl.lineages(adata_cr, copy=True)

        assert isinstance(adata_cr2, AnnData)
        assert adata_cr is not adata_cr2


class TextExcatMCTest:
    def test_invalid_cluster_obs(self, adata_cr: AnnData):
        with pytest.raises(ValueError):
            cr.tl.exact_mc_perm_test(adata_cr, "foo", "bar", "baz")

    def test_invalid_clusters(self, adata_cr: AnnData):
        with pytest.raises(ValueError):
            cr.tl.exact_mc_perm_test(adata_cr, "clusters", "bar", "baz")

    def test_invalid_n_perms(self, adata_cr: AnnData):
        with pytest.raises(ValueError):
            cr.tl.exact_mc_perm_test(
                adata_cr,
                "clusters",
                adata_cr.obs["clusters"].cat.categories[0],
                adata_cr.obs["clusters"].cat.categories[1],
                n_perms=-1,
            )

    def test_invalid_n_counts(self, adata_cr: AnnData):
        with pytest.raises(ValueError):
            cr.tl.exact_mc_perm_test(
                adata_cr,
                "clusters",
                adata_cr.obs["clusters"].cat.categories[0],
                adata_cr.obs["clusters"].cat.categories[1],
                n_counts=-1,
            )

    def test_normal_run(self, adata_cr: AnnData):
        dist, obs, pval = cr.tl.exact_mc_perm_test(
            adata_cr,
            "clusters",
            adata_cr.obs["clusters"].cat.categories[0],
            adata_cr.obs["clusters"].cat.categories[1],
            n_perms=50,
            n_counts=50,
        )

        assert isinstance(dist, list)
        assert isinstance(obs, float)
        assert isinstance(pval, float)


class TestRootFinal:
    pass


class TestTransitionMatrix:
    pass
