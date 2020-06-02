# -*- coding: utf-8 -*-
import pandas as pd
import pytest

from anndata import AnnData

import cellrank as cr
from _helpers import create_model
from cellrank.tools.kernels import Kernel


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
            adata_cr, model, adata_cr.var_names, "0", n_perms=10, fdr_correction=None
        )

        assert "pval" in res.columns
        assert "qval" not in res.columns

    def test_perms_fdr_correction(self, adata_cr: AnnData):
        model = create_model(adata_cr)
        res = cr.tl.gene_importance(
            adata_cr, model, adata_cr.var_names, "0", n_perms=10
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
        with pytest.raises(KeyError):
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
            n_perms=10,
            n_counts=10,
        )

        assert isinstance(dist, list)
        assert isinstance(obs, float)
        assert isinstance(pval, float)


class TestRootFinal:
    def test_find_root(self, adata: AnnData):
        cr.tl.find_root(adata)

    def test_find_final(self, adata: AnnData):
        cr.tl.find_final(adata)

    def test_invalid_cluster_key(self, adata: AnnData):
        with pytest.raises(KeyError):
            cr.tl.find_root(adata, cluster_key="foo")

    def test_invalid_weight(self, adata: AnnData):
        with pytest.raises(ValueError):
            cr.tl.find_root(adata, weight_connectivities=10)

    def test_invalid_percentile(self, adata: AnnData):
        with pytest.raises(ValueError):
            cr.tl.find_root(adata, percentile=110)


class TestTransitionMatrix:
    def test_invalid_velocity_key(self, adata: AnnData):
        with pytest.raises(KeyError):
            cr.tl.transition_matrix(adata, vkey="foo")

    def test_invalid_weight(self, adata: AnnData):
        with pytest.raises(ValueError):
            cr.tl.transition_matrix(adata, weight_connectivities=-1)

    def test_backward(self, adata: AnnData):
        kernel_add = cr.tl.transition_matrix(adata, backward=True)

        assert isinstance(kernel_add, Kernel)


class TestCytoTrace:
    def test_wrong_layer(self, adata: AnnData):
        with pytest.raises(KeyError):
            cr.tl.cyto_trace(adata, layer="foo")

    def test_normal_run(self, adata: AnnData):
        cr.tl.cyto_trace(adata)

        assert "gcs" in adata.obs.keys()
        assert "gene_corr" in adata.var.keys()
        assert "correlates" in adata.var.keys()

    def test_normal_run_copy(self, adata: AnnData):
        adata = cr.tl.cyto_trace(adata, copy=True)

        assert "gcs" in adata.obs.keys()
        assert "gene_corr" in adata.var.keys()
        assert "correlates" in adata.var.keys()
