# -*- coding: utf-8 -*-
import os
from typing import Union
from pathlib import Path

import pytest
from packaging import version

import matplotlib.cm as cm
from matplotlib.testing import setup
from matplotlib.testing.compare import compare_images

import scvelo as scv
from anndata import AnnData

import cellrank as cr
from _helpers import create_model, resize_images_to_same_sizes
from cellrank.tools import GPCCA, CFLARE

setup()

HERE: str = Path(__file__).parent
GT_FIGS = HERE / "_ground_truth_figures"
FIGS = HERE / "figures"
DPI = 40
TOL = 300

GENES = [
    "Tcea1",
    "Tmeff2",
    "Ndufb3",
    "Rpl37a",
    "Arpc2",
    "Ptma",
    "Cntnap5b",
    "Cntnap5a",
    "Mpc2",
    "2010300C02Rik",
]

cr.settings.figdir = FIGS
scv.settings.figdir = str(FIGS)


try:
    from importlib_metadata import version as get_version
except ImportError:
    from importlib.metadata import version as get_version

scvelo_paga_skip = pytest.mark.skipif(
    version.parse(get_version(scv.__name__)) < version.parse("0.1.26.dev189+gc441c72"),
    reason="scVelo < `0.1.26.dev189+gc441c72` supports new PAGA, including node colors and confidence",
)
del version, get_version


def compare(
    *, kind: str = "adata", dirname: Union[str, Path] = None, tol: int = TOL,
):
    def _compare_images(expected_path: Union[str, Path], actual_path: Union[str, Path]):
        resize_images_to_same_sizes(expected_path, actual_path)
        res = compare_images(expected_path, actual_path, tol=tol)
        assert res is None, res

    def compare_cflare_fwd(
        func,
    ):  # mustn't use functools.wraps - it think's the fact that `adata` is fixture
        def decorator(self, adata_cflare_fwd):
            adata, mc = adata_cflare_fwd
            fpath = f"{func.__name__.replace('test_', '')}.png"
            path = str(fpath[7:] if fpath.startswith("scvelo_") else fpath)
            func(self, adata if kind == "adata" else mc, path)

            if dirname is not None:
                for file in os.listdir(FIGS / dirname):
                    _compare_images(GT_FIGS / dirname / file, FIGS / dirname / file)
            else:
                _compare_images(GT_FIGS / fpath, FIGS / fpath)

        return decorator

    def compare_gpcca_fwd(
        func,
    ):  # mustn't use functools.wraps - it think's the fact that `adata` is fixture
        def decorator(self, adata_gpcca_fwd):
            _, gpcca = adata_gpcca_fwd
            fpath = f"{func.__name__.replace('test_', '')}.png"
            path = fpath[7:] if fpath.startswith("scvelo_") else fpath
            func(self, gpcca, path)

            if dirname is not None:
                for file in os.listdir(FIGS / dirname):
                    _compare_images(GT_FIGS / dirname / file, FIGS / dirname / file)
            else:
                _compare_images(GT_FIGS / fpath, FIGS / fpath)

        assert (
            kind == "gpcca"
        ), "Function `compare_gpcca_fwd` only supports `kind='gpcca'`."

        return decorator

    if kind not in ("adata", "cflare", "gpcca"):
        raise ValueError(
            f"Invalid kind `{kind!r}`. Valid options are `'adata'`, `'cflare'` and `'gpcca'`."
        )

    if kind == "gpcca":
        return compare_gpcca_fwd

    return compare_cflare_fwd  # here we hand `kind='adata'`


class TestClusterFates:
    @compare()
    def test_bar(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(adata, "clusters", dpi=DPI, save=fpath)

    @compare()
    def test_bar_cluster_subset(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata, "clusters", clusters=["Astrocytes", "GABA"], dpi=DPI, save=fpath
        )

    @compare()
    def test_bar_lineage_subset(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(adata, "clusters", lineages=["0"], dpi=DPI, save=fpath)

    @compare(tol=250)
    def test_paga_pie(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(adata, "clusters", mode="paga_pie", dpi=DPI, save=fpath)

    @compare(tol=250)
    def test_paga_pie_title(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata, "clusters", mode="paga_pie", title="foo bar baz", dpi=DPI, save=fpath
        )

    @scvelo_paga_skip
    @compare()
    def test_paga_pie_embedding(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata, "clusters", mode="paga_pie", basis="umap", dpi=DPI, save=fpath
        )

    @scvelo_paga_skip
    @compare()
    def test_paga(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(adata, "clusters", mode="paga", dpi=DPI, save=fpath)

    @scvelo_paga_skip
    @compare()
    def test_paga_lineage_subset(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata, "clusters", mode="paga", lineages=["0"], dpi=DPI, save=fpath
        )

    @compare()
    def test_violin(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(adata, "clusters", mode="violin", dpi=DPI, save=fpath)

    @compare()
    def test_violin_cluster_subset(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(adata, "clusters", mode="violin", dpi=DPI, save=fpath)

    @compare()
    def test_violin_lineage_subset(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata, "clusters", mode="violin", lineages=["1"], dpi=DPI, save=fpath
        )

    @compare()
    def test_violin_lineage_subset(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata, "clusters", mode="violin", lineages=["1"], dpi=DPI, save=fpath
        )

    @compare()
    def test_paga_pie_legend_simple(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata,
            cluster_key="clusters",
            mode="paga_pie",
            save=fpath,
            legend_kwargs=(dict(loc="top")),
        )

    @scvelo_paga_skip
    @compare()
    def test_paga_pie_legend_position(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata,
            cluster_key="clusters",
            mode="paga_pie",
            basis="umap",
            save=fpath,
            legend_kwargs=(dict(loc="lower")),
            legend_loc="upper",
        )

    @scvelo_paga_skip
    @compare()
    def test_paga_pie_no_legend(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata,
            cluster_key="clusters",
            mode="paga_pie",
            basis="umap",
            save=fpath,
            legend_kwargs=(dict(loc=None)),
            legend_loc=None,
        )

    @scvelo_paga_skip
    @compare()
    def test_paga_pie_only_abs_prob(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata,
            cluster_key="clusters",
            mode="paga_pie",
            basis="umap",
            save=fpath,
            legend_kwargs=(dict(loc="center")),
            legend_loc=None,
        )

    @scvelo_paga_skip
    @compare()
    def test_paga_pie_only_clusters(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata,
            cluster_key="clusters",
            mode="paga_pie",
            basis="umap",
            save=fpath,
            legend_kwargs=(dict(loc=None)),
            legend_loc="on data",
        )

    @scvelo_paga_skip
    @compare()
    def test_paga_pie_legend_position_out(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata,
            cluster_key="clusters",
            mode="paga_pie",
            basis="umap",
            save=fpath,
            legend_kwargs=(dict(loc="lower left out")),
            legend_loc="center right out",
        )

    def test_paga_pie_wrong_legend_kind_1(self, adata_cflare_fwd):
        adata, _ = adata_cflare_fwd
        with pytest.raises(ValueError):
            cr.pl.cluster_fates(
                adata,
                cluster_key="clusters",
                mode="paga_pie",
                legend_kwargs=(dict(loc="foo")),
            )

    def test_paga_pie_wrong_legend_kind_2(self, adata_cflare_fwd):
        adata, _ = adata_cflare_fwd
        with pytest.raises(ValueError):
            cr.pl.cluster_fates(
                adata,
                cluster_key="clusters",
                mode="paga_pie",
                legend_kwargs=(dict(loc="lower foo")),
            )

    def test_paga_pie_wrong_legend_kind_3(self, adata_cflare_fwd):
        adata, _ = adata_cflare_fwd
        with pytest.raises(ValueError):
            cr.pl.cluster_fates(
                adata,
                cluster_key="clusters",
                mode="paga_pie",
                legend_kwargs=(dict(loc="lower left bar")),
            )

    def test_paga_pie_wrong_legend_kind_4(self, adata_cflare_fwd):
        adata, _ = adata_cflare_fwd
        with pytest.raises(ValueError):
            cr.pl.cluster_fates(
                adata,
                cluster_key="clusters",
                mode="paga_pie",
                legend_kwargs=(dict(loc="lower left foo bar")),
            )

    @compare()
    def test_mode_heatmap(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(adata, "clusters", mode="heatmap", dpi=DPI, save=fpath)

    @compare()
    def test_mode_heatmap_title(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata, "clusters", mode="heatmap", title="foo", dpi=DPI, save=fpath
        )

    @compare()
    def test_mode_heatmap_cmap(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata, "clusters", mode="heatmap", cmap="inferno", dpi=DPI, save=fpath
        )

    @compare()
    def test_mode_heatmap_xticks_rotation(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata, "clusters", mode="heatmap", xticks_rotation=90, dpi=DPI, save=fpath
        )

    @compare()
    def test_mode_heatmap_clusters(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata,
            "clusters",
            mode="heatmap",
            clusters=["Astrocytes", "GABA"],
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_mode_heatmap_lineages(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata, "clusters", mode="heatmap", lineages=["0"], dpi=DPI, save=fpath
        )

    @compare()
    def test_mode_clustermap(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(adata, "clusters", mode="clustermap", dpi=DPI, save=fpath)


class TestClusterLineages:
    @compare()
    def test_cluster_lineage(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.cluster_lineage(
            adata, model, GENES[:10], "0", time_key="latent_time", dpi=DPI, save=fpath,
        )

    @compare()
    def test_cluster_lineage_no_norm(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.cluster_lineage(
            adata,
            model,
            GENES[:10],
            "0",
            time_key="latent_time",
            norm=False,
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_cluster_lineage_data_key(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.cluster_lineage(
            adata,
            model,
            GENES[:10],
            "0",
            time_key="latent_time",
            data_key="Ms",
            norm=False,
            dpi=DPI,
            save=fpath,
        )


class TestHeatmap:
    @compare()
    def test_heatmap_lineages(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            GENES[:10],
            kind="lineages",
            time_key="latent_time",
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_heatmap_genes(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            GENES[:10],
            kind="genes",
            time_key="latent_time",
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_heatmap_cluster_genes(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            GENES[:10],
            kind="lineages",
            time_key="latent_time",
            cluster_genes=True,
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_heatmap_lineage_height(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            GENES[:10],
            kind="lineages",
            time_key="latent_time",
            lineage_height=0.2,
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_heatmap_start_end_clusters(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            GENES[:10],
            kind="lineages",
            time_key="latent_time",
            start_lineage="0",
            end_lineage="1",
            dpi=DPI,
            save=fpath,
        )

    @compare(tol=250)
    def test_heatmap_cmap(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            GENES[:5],
            kind="genes",
            time_key="latent_time",
            cmap=cm.viridis,
            dpi=DPI,
            save=fpath,
        )


class TestGeneTrend:
    @compare(dirname="trends_simple")
    def test_trends(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.gene_trends(
            adata,
            model,
            GENES[:3],
            data_key="Ms",
            dirname="trends_simple",
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_trends_same_plot(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.gene_trends(
            adata, model, GENES[:3], data_key="Ms", same_plot=True, dpi=DPI, save=fpath,
        )

    @compare()
    def test_trends_hide_cells(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.gene_trends(
            adata,
            model,
            GENES[0],
            data_key="Ms",
            same_plot=True,
            hide_cells=True,
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_trends_conf_int(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.gene_trends(
            adata,
            model,
            GENES[0],
            data_key="Ms",
            same_plot=True,
            conf_int=False,
            dpi=DPI,
            save=fpath,
        )

    @compare(dirname="trends_sharey")
    def test_trends_sharey(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.gene_trends(
            adata,
            model,
            GENES[:3],
            data_key="Ms",
            sharey=False,
            dirname="trends_sharey",
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_trends_cbar(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.gene_trends(
            adata,
            model,
            GENES[0],
            data_key="Ms",
            same_plot=True,
            show_cbar=False,
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_trends_lineage_cmap(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.gene_trends(
            adata,
            model,
            GENES[0],
            data_key="Ms",
            same_plot=True,
            lineage_cmap=cm.Set2,
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_trends_lineage_cell_color(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.gene_trends(
            adata,
            model,
            GENES[0],
            data_key="Ms",
            same_plot=True,
            cell_color="red",
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_trend_lw(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.gene_trends(
            adata,
            model,
            GENES[0],
            data_key="Ms",
            same_plot=True,
            lw=10,
            dpi=DPI,
            save=fpath,
        )


class TestGraph:
    @compare()
    def test_graph(self, adata: AnnData, fpath: str):
        cr.pl.graph(adata, "T_fwd", ixs=range(10), dpi=DPI, save=fpath)

    @compare()
    def test_graph_layout(self, adata: AnnData, fpath: str):
        cr.pl.graph(adata, "T_fwd", ixs=range(10), layout="umap", dpi=DPI, save=fpath)

    @compare()
    def test_graph_keys(self, adata: AnnData, fpath: str):
        cr.pl.graph(
            adata,
            "T_fwd",
            ixs=range(10),
            keys=("outgoing", "self_loops"),
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_graph_edge_weight_scale(self, adata: AnnData, fpath: str):
        cr.pl.graph(
            adata, "T_fwd", ixs=range(10), edge_weight_scale=100, dpi=DPI, save=fpath
        )

    @compare()
    def test_graph_show_arrows(self, adata: AnnData, fpath: str):
        cr.pl.graph(
            adata,
            "T_fwd",
            ixs=range(15),
            show_arrows=False,
            edge_weight_scale=100,
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_graph_curved_edges(self, adata: AnnData, fpath: str):
        cr.pl.graph(
            adata, "T_fwd", ixs=range(10), edge_use_curved=False, dpi=DPI, save=fpath
        )

    @compare()
    def test_graph_labels(self, adata: AnnData, fpath: str):
        cr.pl.graph(
            adata, "T_fwd", ixs=range(10), labels=range(10), dpi=DPI, save=fpath
        )

    @compare()
    def test_graph_cmap(self, adata: AnnData, fpath: str):
        cr.pl.graph(
            adata, "T_fwd", ixs=range(10), cont_cmap=cm.inferno, dpi=DPI, save=fpath
        )

    @compare()
    def test_graph_top_n_edges_incoming(self, adata: AnnData, fpath: str):
        cr.pl.graph(
            adata,
            "T_fwd",
            ixs=range(10),
            top_n_edges=(2, True, "incoming"),
            edge_weight_scale=100,
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_graph_top_n_edges_outgoing(self, adata: AnnData, fpath: str):
        cr.pl.graph(
            adata,
            "T_fwd",
            ixs=range(10),
            top_n_edges=(2, False, "outgoing"),
            edge_weight_scale=100,
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_graph_edge_normalize(self, adata: AnnData, fpath: str):
        cr.pl.graph(
            adata, "T_fwd", ixs=range(10), edge_normalize=True, dpi=DPI, save=fpath
        )

    @compare()
    def test_graph_categorical_key(self, adata: AnnData, fpath: str):
        cr.pl.graph(
            adata,
            "T_fwd",
            ixs=range(10),
            keys=["clusters"],
            keylocs=["obs"],
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_trends_size(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.gene_trends(
            adata,
            model,
            GENES[0],
            data_key="Ms",
            same_plot=True,
            size=30,
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_trends_margins(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.gene_trends(
            adata,
            model,
            GENES[0],
            data_key="Ms",
            same_plot=True,
            margins=0.2,
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_trends_cell_alpha(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.gene_trends(
            adata,
            model,
            GENES[0],
            data_key="Ms",
            same_plot=True,
            cell_alpha=0,
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_trends_lineage_alpha(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.gene_trends(
            adata,
            model,
            GENES[0],
            data_key="Ms",
            same_plot=True,
            lineage_alpha=1,
            dpi=DPI,
            save=fpath,
        )


class TestCFLARE:
    @compare(kind="cflare")
    def test_mc_eig(self, mc: CFLARE, fpath: str):
        mc.plot_spectrum(dpi=DPI, save=fpath)

    @compare(kind="cflare")
    def test_mc_complex_spectrum(self, mc: CFLARE, fpath: str):
        mc.plot_spectrum(real_only=False, dpi=DPI, save=fpath)

    @compare(kind="cflare")
    def test_mc_real_spectrum(self, mc: CFLARE, fpath: str):
        mc.plot_spectrum(real_only=True, dpi=DPI, save=fpath)

    @compare(kind="cflare")
    def test_mc_spectrum_title(self, mc: CFLARE, fpath: str):
        mc.plot_spectrum(title="foobar", real_only=False, dpi=DPI, save=fpath)

    @compare(kind="cflare")
    def test_scvelo_eig_embedding_clusters(self, mc: CFLARE, fpath: str):
        mc.plot_eig_embedding(cluster_key="clusters", dpi=DPI, save=fpath)

    @compare(kind="cflare")
    def test_scvelo_eig_embedding_left(self, mc: CFLARE, fpath: str):
        mc.plot_eig_embedding(dpi=DPI, save=fpath)

    @compare(kind="cflare")
    def test_scvelo_eig_embedding_right(self, mc: CFLARE, fpath: str):
        mc.plot_eig_embedding(left=False, dpi=DPI, save=fpath)

    @compare(kind="cflare")
    def test_scvelo_eig_embedding_use_2(self, mc: CFLARE, fpath: str):
        mc.plot_eig_embedding(use=2, dpi=DPI, save=fpath)

    @compare(kind="cflare")
    def test_scvelo_meta_states(self, mc: CFLARE, fpath: str):
        mc.plot_metastable_states(dpi=DPI, save=fpath)

    @compare(kind="cflare")
    def test_scvelo_meta_states(self, mc: CFLARE, fpath: str):
        mc.plot_metastable_states(cluster_key="clusters", dpi=DPI, save=fpath)

    @compare(kind="cflare")
    def test_scvelo_lin_probs(self, mc: CFLARE, fpath: str):
        mc.plot_lin_probs(dpi=DPI, save=fpath)

    @compare(kind="cflare")
    def test_scvelo_lin_probs_clusters(self, mc: CFLARE, fpath: str):
        mc.plot_lin_probs(cluster_key="clusters", dpi=DPI, save=fpath)

    @compare(kind="cflare")
    def test_scvelo_lin_probs_cmap(self, mc: CFLARE, fpath: str):
        mc.plot_lin_probs(cmap=cm.inferno, dpi=DPI, save=fpath)

    @compare(kind="cflare")
    def test_scvelo_lin_probs_lineages(self, mc: CFLARE, fpath: str):
        mc.plot_lin_probs(lineages=["0"], dpi=DPI, save=fpath)

    @compare(kind="cflare")
    def test_scvelo_lin_probs_time(self, mc: CFLARE, fpath: str):
        mc.plot_lin_probs(mode="time", dpi=DPI, save=fpath)


class TestGPCCA:
    @compare(kind="gpcca")
    def test_gpcca_complex_spectrum(self, mc: GPCCA, fpath: str):
        mc.plot_spectrum(real_only=False, dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_gpcca_real_spectrum(self, mc: GPCCA, fpath: str):
        mc.plot_spectrum(real_only=True, dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_gpcca_spectrum_title(self, mc: GPCCA, fpath: str):
        mc.plot_spectrum(title="foobar", real_only=True, dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_gpcca_schur_matrix(self, mc: GPCCA, fpath: str):
        mc.plot_schur_matrix(dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_gpcca_schur_matrix_title(self, mc: GPCCA, fpath: str):
        mc.plot_schur_matrix(title="foobar", dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_gpcca_schur_matrix_cmap(self, mc: GPCCA, fpath: str):
        mc.plot_schur_matrix(cmap=cm.inferno, dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_schur_emb(self, mc: GPCCA, fpath: str):
        mc.plot_schur_embedding(dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_schur_emb_use_2(self, mc: GPCCA, fpath: str):
        mc.plot_schur_embedding(use=1, dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_schur_emb_abs(self, mc: GPCCA, fpath: str):
        mc.plot_schur_embedding(abs_value=True, dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_schur_cluster_key(self, mc: GPCCA, fpath: str):
        mc.plot_schur_embedding(cluster_key="clusters", dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_gpcca_coarse_T(self, mc: GPCCA, fpath: str):
        mc.plot_coarse_T(
            show_initial_dist=False, show_stationary_dist=False, dpi=DPI, save=fpath
        )

    @compare(kind="gpcca")
    def test_gpcca_coarse_T_stat_dist(self, mc: GPCCA, fpath: str):
        mc.plot_coarse_T(
            show_initial_dist=False, show_stationary_dist=True, dpi=DPI, save=fpath
        )

    @compare(kind="gpcca")
    def test_gpcca_coarse_T_init_dist(self, mc: GPCCA, fpath: str):
        mc.plot_coarse_T(
            show_initial_dist=True, show_stationary_dist=False, dpi=DPI, save=fpath
        )

    @compare(kind="gpcca")
    def test_gpcca_coarse_T_no_cbar(self, mc: GPCCA, fpath: str):
        mc.plot_coarse_T(show_cbar=False, dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_gpcca_coarse_T_no_annot(self, mc: GPCCA, fpath: str):
        mc.plot_coarse_T(annotate=False, dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_gpcca_coarse_T_cmap(self, mc: GPCCA, fpath: str):
        mc.plot_coarse_T(cmap=cm.inferno, dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_gpcca_coarse_T_xtick_rot(self, mc: GPCCA, fpath: str):
        mc.plot_coarse_T(xtick_rotation=0, dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_meta_states(self, mc: GPCCA, fpath: str):
        mc.plot_metastable_states(dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_meta_states_lineages(self, mc: GPCCA, fpath: str):
        mc.plot_metastable_states(lineages=["0"], dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_meta_states_discrete(self, mc: GPCCA, fpath: str):
        mc.plot_metastable_states(discrete=True, dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_meta_states_cluster_key(self, mc: GPCCA, fpath: str):
        mc.plot_metastable_states(cluster_key="clusters", dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_meta_states_no_same_plot(self, mc: GPCCA, fpath: str):
        mc.plot_metastable_states(same_plot=False, dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_meta_states_cmap(self, mc: GPCCA, fpath: str):
        mc.plot_metastable_states(cmap=cm.inferno, same_plot=False, dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_meta_states_title(self, mc: GPCCA, fpath: str):
        mc.plot_metastable_states(title="foobar", dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_meta_states_time(self, mc: GPCCA, fpath: str):
        mc.plot_metastable_states(mode="time", dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_main_states(self, mc: GPCCA, fpath: str):
        mc.plot_main_states(dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_main_states_lineages(self, mc: GPCCA, fpath: str):
        mc.plot_main_states(lineages=["0"], dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_main_states_discrete(self, mc: GPCCA, fpath: str):
        mc.plot_main_states(discrete=True, dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_main_states_cluster_key(self, mc: GPCCA, fpath: str):
        mc.plot_main_states(cluster_key="clusters", dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_main_states_no_same_plot(self, mc: GPCCA, fpath: str):
        mc.plot_main_states(same_plot=False, dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_main_states_cmap(self, mc: GPCCA, fpath: str):
        mc.plot_main_states(cmap=cm.inferno, same_plot=False, dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_main_states_title(self, mc: GPCCA, fpath: str):
        mc.plot_main_states(title="foobar", dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_main_states_time(self, mc: GPCCA, fpath: str):
        mc.plot_main_states(mode="time", dpi=DPI, save=fpath)


class TestLineages:
    @compare()
    def test_scvelo_lineages(self, adata: AnnData, fpath: str):
        cr.pl.lineages(adata, dpi=DPI, save=fpath)

    @compare()
    def test_scvelo_lineages_subset(self, adata: AnnData, fpath: str):
        cr.pl.lineages(adata, lineages=["1"], dpi=DPI, save=fpath)

    @compare()
    def test_scvelo_lineages_time(self, adata: AnnData, fpath: str):
        cr.pl.lineages(adata, mode="time", dpi=DPI, save=fpath)

    @compare()
    def test_scvelo_lineages_cmap(self, adata: AnnData, fpath: str):
        cr.pl.lineages(adata, cmap=cm.inferno, dpi=DPI, save=fpath)

    @compare()
    def test_scvelo_lineages_subset(self, adata: AnnData, fpath: str):
        cr.pl.lineages(adata, cluster_key="clusters", dpi=DPI, save=fpath)


class TestSimilarityPlot:
    @compare()
    def test_similarity(self, adata: AnnData, fpath: str):
        cr.pl.similarity_plot(adata, "clusters", n_samples=10, dpi=DPI, save=fpath)

    @compare()
    def test_similarity_clusters(self, adata: AnnData, fpath: str):
        cr.pl.similarity_plot(
            adata,
            "clusters",
            n_samples=10,
            clusters=adata.obs["clusters"].cat.categories,
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_similarity_cmap(self, adata: AnnData, fpath: str):
        cr.pl.similarity_plot(
            adata, "clusters", n_samples=10, cmap=cm.inferno, dpi=DPI, save=fpath
        )

    @compare()
    def test_similarity_fontsize(self, adata: AnnData, fpath: str):
        cr.pl.similarity_plot(
            adata, "clusters", n_samples=10, fontsize=30, dpi=DPI, save=fpath
        )

    @compare()
    def test_similarity_rotation(self, adata: AnnData, fpath: str):
        cr.pl.similarity_plot(
            adata, "clusters", n_samples=10, rotation=90, dpi=DPI, save=fpath
        )


class TestComposition:
    @compare()
    def test_composition(self, adata: AnnData, fpath: str):
        cr.pl.composition(adata, "clusters", dpi=DPI, save=fpath)
