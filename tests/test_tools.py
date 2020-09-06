# -*- coding: utf-8 -*-
import pytest

from anndata import AnnData

import pandas as pd

import cellrank as cr
import cellrank.tl._lineages
from cellrank.tl.kernels import Kernel
from cellrank.tl._constants import AbsProbKey, FinalStatesKey, _probs


class TestLineages:
    def test_no_root_cells(self, adata: AnnData):
        with pytest.raises(RuntimeError):
            cr.tl.lineages(adata)

    def test_normal_run(self, adata_cflare):
        cr.tl.lineages(adata_cflare)

        assert str(AbsProbKey.FORWARD) in adata_cflare.obsm.keys()

    def test_normal_run_copy(self, adata_cflare):
        adata_cr2 = cr.tl.lineages(adata_cflare, copy=True)

        assert isinstance(adata_cr2, AnnData)
        assert adata_cflare is not adata_cr2


class TestLineageDrivers:
    def test_no_abs_probs(self, adata: AnnData):
        with pytest.raises(RuntimeError):
            cellrank.tl._lineages.lineage_drivers(adata)

    def test_normal_run(self, adata_cflare):
        cr.tl.lineages(adata_cflare)
        cellrank.tl._lineages.lineage_drivers(adata_cflare, use_raw=False)

        for name in adata_cflare.obsm[AbsProbKey.FORWARD.s].names:
            assert f"to {name}" in adata_cflare.var

    def test_normal_run_raw(self, adata_cflare):
        adata_cflare.raw = adata_cflare.copy()
        cr.tl.lineages(adata_cflare)
        cellrank.tl._lineages.lineage_drivers(adata_cflare, use_raw=True)

        for name in adata_cflare.obsm[AbsProbKey.FORWARD.s].names:
            assert f"to {name}" in adata_cflare.raw.var

    def test_return_drivers(self, adata_cflare):
        cr.tl.lineages(adata_cflare)
        res = cellrank.tl._lineages.lineage_drivers(
            adata_cflare, return_drivers=True, use_raw=False
        )

        assert isinstance(res, pd.DataFrame)
        assert res.shape[0] == adata_cflare.n_vars
        assert set(res.columns) == set(adata_cflare.obsm[AbsProbKey.FORWARD.s].names)


class TextExcatMCTest:
    def test_invalid_cluster_obs(self, adata_cflare):
        with pytest.raises(KeyError):
            cr.tl._permutation_test(adata_cflare, "foo", "bar", "baz")

    def test_invalid_clusters(self, adata_cflare):
        with pytest.raises(ValueError):
            cr.tl._permutation_test(adata_cflare, "clusters", "bar", "baz")

    def test_invalid_n_perms(self, adata_cflare):
        with pytest.raises(ValueError):
            cr.tl._permutation_test(
                adata_cflare,
                "clusters",
                adata_cflare.obs["clusters"].cat.categories[0],
                adata_cflare.obs["clusters"].cat.categories[1],
                n_perms=-1,
            )

    def test_invalid_n_counts(self, adata_cflare):
        with pytest.raises(ValueError):
            cr.tl._permutation_test(
                adata_cflare,
                "clusters",
                adata_cflare.obs["clusters"].cat.categories[0],
                adata_cflare.obs["clusters"].cat.categories[1],
                n_counts=-1,
            )

    def test_normal_run(self, adata_cflare):
        dist, obs, pval = cr.tl._permutation_test(
            adata_cflare,
            "clusters",
            adata_cflare.obs["clusters"].cat.categories[0],
            adata_cflare.obs["clusters"].cat.categories[1],
            n_perms=10,
            n_counts=10,
        )

        assert isinstance(dist, list)
        assert isinstance(obs, float)
        assert isinstance(pval, float)


class TestRootFinal:
    def test_find_root(self, adata: AnnData):
        cr.tl.initial_states(adata)

        assert str(FinalStatesKey.BACKWARD) in adata.obs.keys()
        assert _probs(FinalStatesKey.BACKWARD) in adata.obs.keys()

    def test_find_final(self, adata: AnnData):
        cr.tl.terminal_states(adata, n_states=5, fit_kwargs=dict(n_cells=5))

        assert str(FinalStatesKey.FORWARD) in adata.obs.keys()
        assert _probs(FinalStatesKey.FORWARD) in adata.obs.keys()

    def test_invalid_cluster_key(self, adata: AnnData):
        with pytest.raises(KeyError):
            cr.tl.initial_states(adata, cluster_key="foo")

    def test_invalid_weight(self, adata: AnnData):
        with pytest.raises(ValueError):
            cr.tl.initial_states(adata, weight_connectivities=10)

    def test_invalid_invalid_mode(self, adata: AnnData):
        with pytest.raises(ValueError):
            cr.tl.initial_states(adata, mode="foobar")


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
