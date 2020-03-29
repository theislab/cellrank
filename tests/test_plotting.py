# -*- coding: utf-8 -*-
from matplotlib.testing import setup
from matplotlib.testing.compare import compare_images
from pathlib import Path
from anndata import AnnData

from _helpers import create_model

import cellrank as cr
import matplotlib.cm as cm

setup()

HERE: Path = Path(__file__).parent
ROOT = HERE / "expected_images"
FIGS = HERE / "figures"
DPI = 40
TOL = 50

cr.settings.figdir = FIGS


def compare(*, kind: str = "adata", backward: bool = False, tol: int = TOL):
    def compare_fwd(
        func
    ):  # mustn't use functools.wraps - it think's `adata` is fixture
        def decorator(self, adata_mc_fwd):
            adata, mc = adata_mc_fwd
            fpath = f"{func.__name__.replace('test_', '')}.png"
            func(self, adata if kind == "adata" else mc, fpath)

            res = compare_images(ROOT / fpath, FIGS / fpath, tol=tol)
            assert res is None, res

        return decorator

    def compare_bwd(func):
        def decorator(self, adata_mc_bwd):
            adata, mc = adata_mc_bwd
            fpath = f"{func.__name__.replace('test_', '')}.png"
            func(self, adata if kind == "adata" else mc, fpath)

            res = compare_images(ROOT / fpath, FIGS / fpath, tol=tol)
            assert res is None, res

        return decorator

    if kind not in ("adata", "mc"):
        raise ValueError(
            f"Invalid kind `{kind!r}`. Valid options are `'adata'`, `'mc'`."
        )

    if backward:
        return compare_bwd
    return compare_fwd


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

    @compare()
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
    @compare()
    def test_cluster_lineage(self, adata: AnnData, fpath: Path):
        model = create_model(adata)
        cr.pl.cluster_lineage(
            adata, model, adata.var_names[:10], "0", time_key="latent_time", save=fpath
        )

    @compare()
    def test_cluster_lineage_no_norm(self, adata: AnnData, fpath: Path):
        model = create_model(adata)
        cr.pl.cluster_lineage(
            adata,
            model,
            adata.var_names[:10],
            "0",
            time_key="latent_time",
            norm=False,
            save=fpath,
        )

    @compare()
    def test_cluster_lineage_data_key(self, adata: AnnData, fpath: Path):
        model = create_model(adata)
        cr.pl.cluster_lineage(
            adata,
            model,
            adata.var_names[:10],
            "0",
            time_key="latent_time",
            norm=False,
            save=fpath,
            data_key="Ms",
        )


class TestGeneTrend:
    @compare()
    def gene_trend(self, adata: AnnData, fpath: Path):
        model = create_model(adata)
        cr.pl.gene_trends(
            adata,
            model,
            adata.var_names[:3],
            time_key="latent_time",
            data_key="Ms",
            save=fpath,
        )


class TestHeatmap:
    @compare()
    def test_heatmap_lineages(self, adata: AnnData, fpath: Path):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            adata.var_names[:10],
            kind="lineages",
            time_key="latent_time",
            save=fpath,
        )

    @compare()
    def test_heatmap_genes(self, adata: AnnData, fpath: Path):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            adata.var_names[:10],
            kind="genes",
            time_key="latent_time",
            save=fpath,
        )

    @compare(backward=True)
    def test_heatmap_backward(self, adata: AnnData, fpath: Path):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            adata.var_names[:10],
            kind="genes",
            time_key="latent_time",
            save=fpath,
        )

    @compare()
    def test_heatmap_cluster_genes(self, adata: AnnData, fpath: Path):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            adata.var_names[:10],
            kind="lineages",
            time_key="latent_time",
            cluster_genes=True,
            save=fpath,
        )

    @compare()
    def test_heatmap_lineage_height(self, adata: AnnData, fpath: Path):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            adata.var_names[:10],
            kind="lineages",
            time_key="latent_time",
            lineage_height=0.2,
            save=fpath,
        )

    @compare()
    def test_heatmap_start_end_clusters(self, adata: AnnData, fpath: Path):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            adata.var_names[:10],
            kind="lineages",
            time_key="latent_time",
            start_clusters="0",
            end_clusters="1",
            save=fpath,
        )

    @compare()
    def test_heatmap_cmap(self, adata: AnnData, fpath: Path):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            adata.var_names[:5],
            kind="genes",
            time_key="latent_time",
            cmap=cm.viridis,
            save=fpath,
        )


class TestGraph:
    @compare()
    def test_graph(self, adata: AnnData, fpath: Path):
        cr.pl.graph(adata, "T_fwd", ixs=range(10), save=fpath)

    @compare()
    def test_graph_layout(self, adata: AnnData, fpath: Path):
        cr.pl.graph(adata, "T_fwd", ixs=range(10), layout="umap", save=fpath)

    @compare()
    def test_graph_keys(self, adata: AnnData, fpath: Path):
        cr.pl.graph(
            adata, "T_fwd", ixs=range(10), keys=("outgoing", "self_loops"), save=fpath
        )

    @compare()
    def test_graph_edge_weight_scale(self, adata: AnnData, fpath: Path):
        cr.pl.graph(adata, "T_fwd", ixs=range(10), edge_weight_scale=100, save=fpath)

    @compare()
    def test_graph_show_arrows(self, adata: AnnData, fpath: Path):
        cr.pl.graph(
            adata,
            "T_fwd",
            ixs=range(15),
            show_arrows=False,
            edge_weight_scale=100,
            save=fpath,
        )

    @compare()
    def test_graph_curved_edges(self, adata: AnnData, fpath: Path):
        cr.pl.graph(adata, "T_fwd", ixs=range(10), edge_use_curved=False, save=fpath)

    @compare()
    def test_graph_labels(self, adata: AnnData, fpath: Path):
        cr.pl.graph(adata, "T_fwd", ixs=range(10), labels=range(10), save=fpath)

    @compare()
    def test_graph_cmap(self, adata: AnnData, fpath: Path):
        cr.pl.graph(adata, "T_fwd", ixs=range(10), cont_cmap=cm.inferno, save=fpath)

    @compare()
    def test_graph_top_n_edges_incoming(self, adata: AnnData, fpath: Path):
        cr.pl.graph(
            adata,
            "T_fwd",
            ixs=range(10),
            top_n_edges=(2, True, "incoming"),
            edge_weight_scale=100,
            save=fpath,
        )

    @compare()
    def test_graph_top_n_edges_outgoing(self, adata: AnnData, fpath: Path):
        cr.pl.graph(
            adata,
            "T_fwd",
            ixs=range(10),
            top_n_edges=(2, False, "outgoing"),
            edge_weight_scale=100,
            save=fpath,
        )

    @compare()
    def test_graph_edge_normalize(self, adata: AnnData, fpath: Path):
        cr.pl.graph(adata, "T_fwd", ixs=range(10), edge_normalize=True, save=fpath)

    @compare()
    def test_graph_categorical_key(self, adata: AnnData, fpath: Path):
        cr.pl.graph(
            adata,
            "T_fwd",
            ixs=range(10),
            keys=["clusters"],
            keylocs=["obs"],
            save=fpath,
        )
