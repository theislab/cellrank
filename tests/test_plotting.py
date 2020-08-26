# -*- coding: utf-8 -*-
import os
from typing import Tuple, Union, Callable
from pathlib import Path

import pytest
from _helpers import create_model, resize_images_to_same_sizes
from packaging import version

import scvelo as scv
from anndata import AnnData

import numpy as np
import pandas as pd

import matplotlib.cm as cm
from matplotlib.testing import setup
from matplotlib.testing.compare import compare_images

import cellrank as cr
import cellrank.pl._lineages
from cellrank.tl._constants import AbsProbKey
from cellrank.tl.estimators import GPCCA, CFLARE
from cellrank.pl._cluster_fates import similarity_plot

setup()

HERE: str = Path(__file__).parent
GT_FIGS = HERE / "_ground_truth_figures"
FIGS = HERE / "figures"
DPI = 40
TOL = 150

# both are for `50` adata
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
RAW_GENES = [
    "Synpr",
    "Rps24",
    "Erc2",
    "Mbnl2",
    "Thoc7",
    "Itm2b",
    "Pcdh9",
    "Fgf14",
    "Rpl37",
    "Cdh9",
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
) -> Callable:
    def _compare_images(
        expected_path: Union[str, Path], actual_path: Union[str, Path]
    ) -> None:
        resize_images_to_same_sizes(expected_path, actual_path)
        res = compare_images(expected_path, actual_path, tol=tol)
        assert res is None, res

    def _prepare_fname(func: Callable) -> Tuple[str, str]:
        fpath = f"{func.__name__.replace('test_', '')}"
        # scvelo saves figures as pdf
        return fpath, str(fpath[7:] + ".png" if fpath.startswith("scvelo_") else fpath)

    def _assert_equal(fpath: str) -> None:
        # TODO: not an elegant solution, consider passing dirname to the functions
        if not fpath.endswith(".png"):
            fpath += ".png"
        if dirname is not None:
            for file in os.listdir(FIGS / dirname):
                _compare_images(GT_FIGS / dirname / file, FIGS / dirname / file)
        else:
            _compare_images(GT_FIGS / fpath, FIGS / fpath)

    def compare_cflare_fwd(
        func: Callable,
    ) -> Callable:  # mustn't use functools.wraps - it think's the fact that `adata` is fixture
        def decorator(self, adata_cflare_fwd) -> None:
            adata, mc = adata_cflare_fwd
            fpath, path = _prepare_fname(func)

            func(self, adata if kind == "adata" else mc, path)

            _assert_equal(fpath)

        return decorator

    def compare_gpcca_fwd(func: Callable) -> Callable:
        def decorator(self, adata_gpcca_fwd) -> None:
            _, gpcca = adata_gpcca_fwd
            fpath, path = _prepare_fname(func)

            func(self, gpcca, path)

            _assert_equal(fpath)

        assert (
            kind == "gpcca"
        ), "Function `compare_gpcca_fwd` only supports `kind='gpcca'`."

        return decorator

    def compare_lineage(func: Callable):
        def decorator(self, lineage):
            path, fpath = _prepare_fname(func)

            func(self, lineage, path)

            _assert_equal(fpath)

        assert (
            kind == "lineage"
        ), "Function `compare_lineage` only supports `kind='lineage'`."

        return decorator

    if kind not in ("adata", "cflare", "gpcca", "lineage"):
        raise ValueError(
            f"Invalid kind `{kind!r}`. Valid options are `'adata'`, `'cflare'` and `'gpcca'`."
        )

    if kind == "gpcca":
        return compare_gpcca_fwd
    if kind == "lineage":
        return compare_lineage

    return compare_cflare_fwd  # here we hand `kind='adata'`


class TestClusterFates:
    @compare()
    def test_bar(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata, cluster_key="clusters", mode="bar", dpi=DPI, save=fpath
        )

    @compare()
    def test_bar_cluster_subset(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata,
            cluster_key="clusters",
            mode="bar",
            clusters=["Astrocytes", "GABA"],
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_bar_lineage_subset(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata,
            cluster_key="clusters",
            mode="bar",
            lineages=["0"],
            dpi=DPI,
            save=fpath,
        )

    @compare(tol=250)
    def test_paga_pie(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata, cluster_key="clusters", mode="paga_pie", dpi=DPI, save=fpath
        )

    @compare(tol=250)
    def test_paga_pie_title(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata,
            cluster_key="clusters",
            mode="paga_pie",
            title="foo bar baz",
            dpi=DPI,
            save=fpath,
        )

    @scvelo_paga_skip
    @compare()
    def test_paga_pie_embedding(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata,
            cluster_key="clusters",
            mode="paga_pie",
            basis="umap",
            dpi=DPI,
            save=fpath,
        )

    @scvelo_paga_skip
    @compare()
    def test_paga(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata, cluster_key="clusters", mode="paga", dpi=DPI, save=fpath
        )

    @scvelo_paga_skip
    @compare()
    def test_paga_lineage_subset(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata,
            cluster_key="clusters",
            mode="paga",
            lineages=["0"],
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_violin(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata, cluster_key="clusters", mode="violin", dpi=DPI, save=fpath
        )

    @compare()
    def test_violin_cluster_subset(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata, cluster_key="clusters", mode="violin", dpi=DPI, save=fpath
        )

    @compare()
    def test_violin_lineage_subset(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata,
            cluster_key="clusters",
            mode="violin",
            lineages=["1"],
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_violin_lineage_subset(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata,
            cluster_key="clusters",
            mode="violin",
            lineages=["1"],
            dpi=DPI,
            save=fpath,
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

    def test_invalid_mode(self, adata_cflare_fwd):
        adata, _ = adata_cflare_fwd
        with pytest.raises(ValueError):
            cr.pl.cluster_fates(
                adata, cluster_key="clusters", mode="foobar",
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
        cr.pl.cluster_fates(
            adata, cluster_key="clusters", mode="heatmap", dpi=DPI, save=fpath
        )

    @compare()
    def test_mode_heatmap_title(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata,
            cluster_key="clusters",
            mode="heatmap",
            title="foo",
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_mode_heatmap_cmap(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata,
            cluster_key="clusters",
            mode="heatmap",
            cmap="inferno",
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_mode_heatmap_xticks_rotation(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata,
            cluster_key="clusters",
            mode="heatmap",
            xticks_rotation=90,
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_mode_heatmap_clusters(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata,
            cluster_key="clusters",
            mode="heatmap",
            clusters=["Astrocytes", "GABA"],
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_mode_heatmap_lineages(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata,
            cluster_key="clusters",
            mode="heatmap",
            lineages=["0"],
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_mode_clustermap(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata, cluster_key="clusters", mode="clustermap", dpi=DPI, save=fpath
        )


class TestClusterLineages:
    @compare()
    def test_cluster_lineage(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.cluster_lineage(
            adata, model, GENES[:10], "0", time_key="latent_time", dpi=DPI, save=fpath,
        )

    @compare()
    def test_cluster_lineage_raw(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.cluster_lineage(
            adata,
            model,
            RAW_GENES[:5],
            "0",
            time_key="latent_time",
            dpi=DPI,
            save=fpath,
            use_raw=True,
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
    @compare(dirname="heatmap_lineages")
    def test_heatmap_lineages(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            GENES[:10],
            mode="lineages",
            time_key="latent_time",
            dpi=DPI,
            save=fpath,
        )

    @compare(dirname="heatmap_lineages_raw")
    def test_heatmap_lineages_raw(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            RAW_GENES[:5],
            mode="lineages",
            time_key="latent_time",
            use_raw=True,
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
            mode="genes",
            time_key="latent_time",
            dpi=DPI,
            save=fpath,
        )

    @compare(dirname="heatmap_no_cluster_genes")
    def test_heatmap_no_cluster_genes(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            GENES[:10],
            cluster_genes=False,
            mode="lineages",
            time_key="latent_time",
            dpi=DPI,
            save=fpath,
        )

    @compare(dirname="heatmap_cluster_genes")
    def test_heatmap_cluster_genes(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            GENES[:10],
            mode="lineages",
            time_key="latent_time",
            cluster_genes=True,
            dpi=DPI,
            save=fpath,
        )

    @compare(dirname="heatmap_lineage_height")
    def test_heatmap_lineage_height(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            GENES[:10],
            mode="lineages",
            time_key="latent_time",
            lineage_height=0.2,
            dpi=DPI,
            save=fpath,
        )

    @compare(dirname="heatmap_time_range")
    def test_heatmap_time_range(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            GENES[:10],
            mode="lineages",
            time_key="latent_time",
            time_range=(0.2, 0.5),
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
            mode="genes",
            time_key="latent_time",
            cmap=cm.viridis,
            dpi=DPI,
            save=fpath,
        )

    @compare(dirname="heatmap_no_cbar_lineages")
    def test_heatmap_no_cbar_lineages(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            GENES[:5],
            mode="lineages",
            time_key="latent_time",
            show_cbar=False,
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_heatmap_no_cbar_genes(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            GENES[:5],
            mode="genes",
            time_key="latent_time",
            show_cbar=False,
            dpi=DPI,
            save=fpath,
        )

    @compare(dirname="heatmap_abs_probs_lineages")
    def test_heatmap_abs_probs_lineages(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            GENES[:5],
            mode="lineages",
            time_key="latent_time",
            show_absorption_probabilities=True,
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_heatmap_abs_probs_genes(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            GENES[:5],
            mode="genes",
            time_key="latent_time",
            show_absorption_probabilities=True,
            dpi=DPI,
            save=fpath,
        )

    @compare(dirname="heatmap_no_convolve")
    def test_heatmap_no_convolve(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            GENES[:5],
            mode="lineages",
            time_key="latent_time",
            n_convolve=None,
            dpi=DPI,
            save=fpath,
        )

    @compare(dirname="heatmap_no_scale_lineages")
    def test_heatmap_no_scale_lineages(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            GENES[:5],
            mode="lineages",
            time_key="latent_time",
            scale=False,
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_heatmap_no_scale_genes(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            GENES[:5],
            mode="genes",
            time_key="latent_time",
            scale=False,
            dpi=DPI,
            save=fpath,
        )

    @compare(dirname="heatmap_cluster_no_scale")
    def test_heatmap_cluster_no_scale(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            GENES[:5],
            mode="lineages",
            time_key="latent_time",
            scale=False,
            cluster_genes=True,
            dpi=DPI,
            save=fpath,
        )

    @compare(dirname="heatmap_no_cluster")
    def test_heatmap_no_cluster(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            GENES[:5],
            mode="lineages",
            time_key="latent_time",
            cluster_genes=False,
            dpi=DPI,
            save=fpath,
        )

    @compare(dirname="heatmap_cluster_key_no_abs_probs")
    def test_heatmap_cluster_key_no_abs_probs(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            GENES[:5],
            cluster_key="clusters",
            show_absorption_probabilities=False,
            mode="lineages",
            time_key="latent_time",
            dpi=DPI,
            save=fpath,
        )

    @compare(dirname="heatmap_cluster_key_abs_probs")
    def test_heatmap_cluster_key_abs_probs(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            GENES[:5],
            cluster_key="clusters",
            show_absorption_probabilities=True,
            mode="lineages",
            time_key="latent_time",
            dpi=DPI,
            save=fpath,
        )

    @compare(dirname="heatmap_multiple_cluster_keys")
    def test_heatmap_multiple_cluster_keys(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            GENES[:5],
            cluster_key=["clusters", "clusters_enlarged", "clusters"],
            show_absorption_probabilities=True,
            mode="lineages",
            time_key="latent_time",
            dpi=DPI,
            save=fpath,
        )

    @compare(dirname="heatmap_multiple_cluster_keys_show_all_genes")
    def test_heatmap_multiple_cluster_keys_show_all_genes(
        self, adata: AnnData, fpath: str
    ):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            GENES[:10],
            mode="lineages",
            time_key="latent_time",
            dpi=DPI,
            save=fpath,
        )

    @compare(dirname="heatmap_n_jobs")
    def test_heatmap_n_jobs(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            GENES[:5],
            n_jobs=2,
            backend="threading",
            cluster_key=["clusters", "clusters_enlarged", "clusters"],
            show_absorption_probabilities=True,
            mode="lineages",
            time_key="latent_time",
            dpi=DPI,
            save=fpath,
        )

    @compare(dirname="heatmap_n_jobs_multiprocessing")
    def test_heatmap_n_jobs_multiprocessing(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            GENES[:5],
            n_jobs=2,
            backend="loky",  # uses pickling of objects, such as Lineage
            cluster_key=["clusters", "clusters_enlarged", "clusters"],
            show_absorption_probabilities=True,
            mode="lineages",
            time_key="latent_time",
            dpi=DPI,
            save=fpath,
        )

    @compare(dirname="heatmap_keep_gene_order")
    def test_heatmap_keep_gene_order(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            GENES[:10],
            mode="lineages",
            time_key="latent_time",
            keep_gene_order=True,
            dpi=DPI,
            save=fpath,
        )

    @compare(dirname="heatmap_show_dendrogram")
    def test_heatmap_show_dendrogram(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            GENES[:10],
            mode="lineages",
            time_key="latent_time",
            cluster_genes=True,
            show_dendrogram=True,
            dpi=DPI,
            save=fpath,
        )


class TestHeatMapReturns:
    def test_heatmap_lineages_return_genes(self, adata_cflare: AnnData):
        model = create_model(adata_cflare)
        df = cr.pl.heatmap(
            adata_cflare,
            model,
            GENES[:10],
            mode="lineages",
            time_key="latent_time",
            return_genes=True,
            dpi=DPI,
        )

        assert isinstance(df, pd.DataFrame)
        np.testing.assert_array_equal(
            df.columns, adata_cflare.obsm[AbsProbKey.FORWARD.s].names
        )
        assert len(df) == 10
        assert set(df.iloc[:, 0].values) == set(GENES[:10])

    def test_heatmap_lineages_return_genes_large_number(self, adata_cflare: AnnData):
        model = create_model(adata_cflare)
        genes = adata_cflare.var_names[:100]
        df = cr.pl.heatmap(
            adata_cflare,
            model,
            genes,
            mode="lineages",
            time_key="latent_time",
            return_genes=True,
            dpi=DPI,
        )

        assert isinstance(df, pd.DataFrame)
        np.testing.assert_array_equal(
            df.columns, adata_cflare.obsm[AbsProbKey.FORWARD.s].names
        )
        assert len(df) == len(genes)
        assert set(df.iloc[:, 0].values) == set(genes)

    def test_heatmap_lineages_return_genes_same_order(self, adata_cflare: AnnData):
        model = create_model(adata_cflare)
        df = cr.pl.heatmap(
            adata_cflare,
            model,
            GENES[:10],
            keep_gene_order=True,
            mode="lineages",
            time_key="latent_time",
            return_genes=True,
            dpi=DPI,
        )

        assert isinstance(df, pd.DataFrame)
        np.testing.assert_array_equal(
            df.columns, adata_cflare.obsm[AbsProbKey.FORWARD.s].names
        )
        assert len(df) == 10
        assert set(df.iloc[:, 0].values) == set(GENES[:10])

        ref = df.iloc[:, 0].values
        for i in range(1, len(df.columns)):
            np.testing.assert_array_equal(df.iloc[:, i].values, ref)

    def test_heatmap_genes_return_genes(self, adata_cflare: AnnData):
        model = create_model(adata_cflare)
        df = cr.pl.heatmap(
            adata_cflare,
            model,
            GENES[:10],
            mode="genes",
            time_key="latent_time",
            cluster_genes=True,
            show_dendrogram=True,
            return_genes=True,
            dpi=DPI,
        )

        assert df is None


class TestGeneTrend:
    @compare()
    def test_trends(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.gene_trends(
            adata, model, GENES[:3], data_key="Ms", dpi=DPI, save=fpath,
        )

    @compare()
    def test_trends_raw(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.gene_trends(
            adata,
            model,
            RAW_GENES[:5],
            data_key="X",
            use_raw=True,
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

    @compare()
    def test_trends_sharey(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.gene_trends(
            adata, model, GENES[:3], data_key="Ms", sharey="row", dpi=DPI, save=fpath,
        )

    @compare()
    def test_trends_sharex(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.gene_trends(
            adata,
            model,
            GENES[:10],
            ncols=3,
            data_key="Ms",
            sharex="all",
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_trends_gene_as_title(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.gene_trends(
            adata,
            model,
            GENES[:10],
            gene_as_title=False,
            same_plot=True,
            data_key="Ms",
            sharex="all",
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_trends_gene_no_legend(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.gene_trends(
            adata,
            model,
            GENES[:10],
            legend_loc=None,
            data_key="Ms",
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_trends_no_cbar(self, adata: AnnData, fpath: str):
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
    def test_trends_abs_prob_cmap(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.gene_trends(
            adata,
            model,
            GENES[0],
            data_key="Ms",
            same_plot=False,
            hide_cells=False,
            abs_prob_cmap=cm.inferno,
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

    @compare()
    def test_trend_suptitle(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.gene_trends(
            adata,
            model,
            GENES[:10],
            suptitle="FOOBAR",
            data_key="Ms",
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

    @compare()
    def test_trend_time_range(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.gene_trends(
            adata,
            model,
            GENES[:10],
            data_key="Ms",
            same_plot=False,
            time_range=(0, 0.5),
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_trend_perc(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.gene_trends(
            adata,
            model,
            GENES[:10],
            data_key="Ms",
            same_plot=False,
            perc=(0, 50),
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


class TestCFLARE:
    @compare(kind="cflare")
    def test_mc_spectrum(self, mc: CFLARE, fpath: str):
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
    def test_mc_spectrum_evals(self, mc: CFLARE, fpath: str):
        mc.plot_spectrum(2, real_only=True, dpi=DPI, save=fpath)

    @compare(kind="cflare")
    def test_mc_spectrum_evals_complex(self, mc: CFLARE, fpath: str):
        mc.plot_spectrum(2, real_only=False, dpi=DPI, save=fpath)

    @compare(kind="cflare")
    def test_scvelo_eigendecomposition_clusters(self, mc: CFLARE, fpath: str):
        mc.plot_eigendecomposition(cluster_key="clusters", dpi=DPI, save=fpath)

    @compare(kind="cflare")
    def test_scvelo_eigendecomposition_left(self, mc: CFLARE, fpath: str):
        mc.plot_eigendecomposition(dpi=DPI, save=fpath)

    @compare(kind="cflare")
    def test_scvelo_eigendecomposition_right(self, mc: CFLARE, fpath: str):
        mc.plot_eigendecomposition(left=False, dpi=DPI, save=fpath)

    @compare(kind="cflare")
    def test_scvelo_eigendecomposition_use_2(self, mc: CFLARE, fpath: str):
        mc.plot_eigendecomposition(use=2, dpi=DPI, save=fpath)

    @compare(kind="cflare")
    def test_scvelo_final_states(self, mc: CFLARE, fpath: str):
        mc.plot_final_states(dpi=DPI, save=fpath)

    @compare(kind="cflare")
    def test_scvelo_final_states_clusters(self, mc: CFLARE, fpath: str):
        mc.plot_final_states(cluster_key="clusters", dpi=DPI, save=fpath)

    @compare(kind="cflare")
    def test_scvelo_lin_probs(self, mc: CFLARE, fpath: str):
        mc.plot_absorption_probabilities(dpi=DPI, save=fpath)

    @compare(kind="cflare")
    def test_scvelo_lin_probs_clusters(self, mc: CFLARE, fpath: str):
        mc.plot_absorption_probabilities(cluster_key="clusters", dpi=DPI, save=fpath)

    @compare(kind="cflare")
    def test_scvelo_lin_probs_cmap(self, mc: CFLARE, fpath: str):
        mc.plot_absorption_probabilities(cmap=cm.inferno, dpi=DPI, save=fpath)

    @compare(kind="cflare")
    def test_scvelo_lin_probs_lineages(self, mc: CFLARE, fpath: str):
        mc.plot_absorption_probabilities(lineages=["0"], dpi=DPI, save=fpath)

    @compare(kind="cflare")
    def test_scvelo_lin_probs_time(self, mc: CFLARE, fpath: str):
        mc.plot_absorption_probabilities(mode="time", dpi=DPI, save=fpath)


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
    def test_gpcca_spectrum_evals(self, mc: CFLARE, fpath: str):
        mc.plot_spectrum(2, real_only=True, dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_gpcca_spectrum_evals_complex(self, mc: CFLARE, fpath: str):
        mc.plot_spectrum(2, real_only=False, dpi=DPI, save=fpath)

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
        mc.plot_schur(dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_schur_emb_use_2(self, mc: GPCCA, fpath: str):
        mc.plot_schur(use=1, dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_schur_emb_abs(self, mc: GPCCA, fpath: str):
        mc.plot_schur(abs_value=True, dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_schur_cluster_key(self, mc: GPCCA, fpath: str):
        mc.plot_schur(cluster_key="clusters", dpi=DPI, save=fpath)

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
    def test_scvelo_gpcca_final_states(self, mc: GPCCA, fpath: str):
        mc.plot_final_states(dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_final_states_lineages(self, mc: GPCCA, fpath: str):
        mc.plot_final_states(lineages=["0"], dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_final_states_discrete(self, mc: GPCCA, fpath: str):
        mc.plot_final_states(discrete=True, dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_final_states_cluster_key(self, mc: GPCCA, fpath: str):
        mc.plot_final_states(cluster_key="clusters", dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_final_states_no_same_plot(self, mc: GPCCA, fpath: str):
        mc.plot_final_states(same_plot=False, dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_final_states_cmap(self, mc: GPCCA, fpath: str):
        mc.plot_final_states(cmap=cm.inferno, same_plot=False, dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_final_states_title(self, mc: GPCCA, fpath: str):
        mc.plot_final_states(title="foobar", dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_final_states_time(self, mc: GPCCA, fpath: str):
        mc.plot_final_states(mode="time", dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_abs_probs_disc_same(self, mc: GPCCA, fpath: str):
        mc.plot_absorption_probabilities(
            cluster_key="clusters", discrete=True, same_plot=True, dpi=DPI, save=fpath
        )

    @compare(kind="gpcca")
    def test_scvelo_gpcca_abs_probs_disc_not_same(self, mc: GPCCA, fpath: str):
        mc.plot_absorption_probabilities(
            cluster_key="clusters", discrete=True, same_plot=False, dpi=DPI, save=fpath
        )

    @compare(kind="gpcca")
    def test_scvelo_gpcca_abs_probs_cont_same_no_clusters(self, mc: GPCCA, fpath: str):
        mc.plot_absorption_probabilities(
            discrete=False, same_plot=True, dpi=DPI, save=fpath
        )

    @compare(kind="gpcca")
    def test_scvelo_gpcca_abs_probs_cont_same_clusters(self, mc: GPCCA, fpath: str):
        mc.plot_absorption_probabilities(
            cluster_key="clusters", discrete=False, same_plot=True, dpi=DPI, save=fpath
        )

    @compare(kind="gpcca")
    def test_scvelo_gpcca_abs_probs_cont_not_same(self, mc: GPCCA, fpath: str):
        mc.plot_absorption_probabilities(
            cluster_key="clusters", discrete=False, same_plot=False, dpi=DPI, save=fpath
        )


class TestLineages:
    @compare()
    def test_scvelo_lineages(self, adata: AnnData, fpath: str):
        cellrank.pl._lineages.lineages(adata, dpi=DPI, save=fpath)

    @compare()
    def test_scvelo_lineages_subset(self, adata: AnnData, fpath: str):
        cellrank.pl._lineages.lineages(adata, lineages=["1"], dpi=DPI, save=fpath)

    @compare()
    def test_scvelo_lineages_time(self, adata: AnnData, fpath: str):
        cellrank.pl._lineages.lineages(adata, mode="time", dpi=DPI, save=fpath)

    @compare()
    def test_scvelo_lineages_cmap(self, adata: AnnData, fpath: str):
        cellrank.pl._lineages.lineages(adata, cmap=cm.inferno, dpi=DPI, save=fpath)

    @compare()
    def test_scvelo_lineages_subset(self, adata: AnnData, fpath: str):
        cellrank.pl._lineages.lineages(
            adata, cluster_key="clusters", dpi=DPI, save=fpath
        )


class TestSimilarityPlot:
    @compare()
    def test_similarity(self, adata: AnnData, fpath: str):
        similarity_plot(adata, "clusters", n_samples=10, dpi=DPI, save=fpath)

    @compare()
    def test_similarity_clusters(self, adata: AnnData, fpath: str):
        similarity_plot(
            adata,
            "clusters",
            n_samples=10,
            clusters=adata.obs["clusters"].cat.categories,
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_similarity_cmap(self, adata: AnnData, fpath: str):
        similarity_plot(
            adata, "clusters", n_samples=10, cmap=cm.inferno, dpi=DPI, save=fpath
        )

    @compare()
    def test_similarity_fontsize(self, adata: AnnData, fpath: str):
        similarity_plot(
            adata, "clusters", n_samples=10, fontsize=30, dpi=DPI, save=fpath
        )

    @compare()
    def test_similarity_rotation(self, adata: AnnData, fpath: str):
        similarity_plot(
            adata, "clusters", n_samples=10, rotation=90, dpi=DPI, save=fpath
        )


class TestComposition:
    @compare()
    def test_composition(self, adata: AnnData, fpath: str):
        cr.pl.composition(adata, "clusters", dpi=DPI, save=fpath)


class TestLineage:
    @compare(kind="lineage")
    def test_pie(self, lineage: cr.tl.Lineage, fpath: str):
        lineage.plot_pie(dpi=DPI, save=fpath)

    @compare(kind="lineage")
    def test_pie_reduction(self, lineage: cr.tl.Lineage, fpath: str):
        lineage.plot_pie(reduction=np.var, dpi=DPI, save=fpath)

    @compare(kind="lineage")
    def test_pie_title(self, lineage: cr.tl.Lineage, fpath: str):
        lineage.plot_pie(title="FOOBAR", dpi=DPI, save=fpath)

    @compare(kind="lineage")
    def test_pie_t(self, lineage: cr.tl.Lineage, fpath: str):
        lineage.T.plot_pie(dpi=DPI, save=fpath)


# TODO: more model tests + cr.pl.lineage_drivers
class TestModel:
    @compare()
    def test_model_default(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        model.prepare(adata.var_names[0], "0")
        model.fit()
        model.predict()
        model.confidence_interval()
        model.plot(save=fpath)
