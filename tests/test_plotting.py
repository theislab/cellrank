# -*- coding: utf-8 -*-
from matplotlib.testing import setup
from matplotlib.testing.compare import compare_images
from pathlib import Path
from anndata import AnnData

import cellrank as cr

setup()

HERE: Path = Path(__file__).parent
ROOT = HERE / "expected_images"
FIGS = HERE / "figures"
DPI = 40
TOL = 50

cr.settings.figdir = FIGS


def compare(*, kind: str = "adata", use_model: bool = False, tol: int = TOL):
    def compare_fwd(
        func
    ):  # mustn't use functools.wraps - it think's `adata` is fixture
        def decorator(self, adata_mc_fwd, model):
            adata, mc = adata_mc_fwd
            data = adata if kind == "adata" else mc
            fpath = f"{func.__name__.replace('test_', '')}.png"
            if use_model:
                func(self, data, fpath)
            else:
                func(self, data, fpath, model)

            res = compare_images(ROOT / fpath, FIGS / fpath, tol=tol)
            assert res is None, res

        return decorator

    def compare_bwd(func):
        def decorator(self, adata_mc_bwd, model):
            adata, mc = adata_mc_bwd
            data = adata if kind == "adata" else mc
            fpath = f"{func.__name__.replace('test_', '')}.png"
            if use_model:
                func(self, data, fpath)
            else:
                func(self, data, fpath, model)

            res = compare_images(ROOT / fpath, FIGS / fpath, tol=tol)
            assert res is None, res

        return decorator

    if kind in ("adata", "mc_fwd"):
        return compare_fwd
    elif kind == "mc_bwd":
        return compare_bwd
    else:
        raise ValueError(
            f"Invalid kind `{kind}`. Valid options are `'adata`, `'mc_fwd'`, `'mc_bwd'`."
        )


class TestClusterFates:
    @compare()
    def test_bar(self, adata: AnnData, fpath: Path):
        cr.pl.cluster_fates(adata, "clusters", dpi=DPI, save=fpath)

    @compare()
    def test_bar_cluster_subset(self, adata: AnnData, fpath: Path):
        cr.pl.cluster_fates(
            adata, "clusters", clusters=["Astrocytes", "GABA"], dpi=DPI, save=fpath
        )

    @compare()
    def test_bar_lineage_subset(self, adata: AnnData, fpath: Path):
        cr.pl.cluster_fates(adata, "clusters", lineages=["0"], dpi=DPI, save=fpath)

    @compare()
    def test_paga_pie(self, adata: AnnData, fpath: Path):
        cr.pl.cluster_fates(adata, "clusters", mode="paga_pie", dpi=DPI, save=fpath)

    @compare()
    def test_paga_pie_embedding(self, adata: AnnData, fpath: Path):
        cr.pl.cluster_fates(
            adata, "clusters", mode="paga_pie", basis="umap", dpi=DPI, save=fpath
        )

    @compare(tol=100)
    def test_paga(self, adata: AnnData, fpath: Path):
        cr.pl.cluster_fates(adata, "clusters", mode="paga", dpi=DPI, save=fpath)

    @compare()
    def test_paga_lineage_subset(self, adata: AnnData, fpath: Path):
        cr.pl.cluster_fates(
            adata, "clusters", mode="paga", lineages=["0"], dpi=DPI, save=fpath
        )

    @compare()
    def test_violin(self, adata: AnnData, fpath: Path):
        cr.pl.cluster_fates(adata, "clusters", mode="violin", dpi=DPI, save=fpath)

    @compare()
    def test_violin_cluster_subset(self, adata: AnnData, fpath: Path):
        cr.pl.cluster_fates(adata, "clusters", mode="violin", dpi=DPI, save=fpath)

    @compare()
    def test_violin_lineage_subset(self, adata: AnnData, fpath: Path):
        cr.pl.cluster_fates(
            adata, "clusters", mode="violin", lineages=["1"], dpi=DPI, save=fpath
        )


class TestClusterLineages:
    pass
