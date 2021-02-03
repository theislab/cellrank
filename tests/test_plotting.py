import os
from typing import Tuple, Union, Callable
from pathlib import Path

import pytest
from _helpers import (
    gamr_skip,
    create_model,
    create_failed_model,
    resize_images_to_same_sizes,
)
from packaging import version

import scvelo as scv
from anndata import AnnData

import numpy as np
import pandas as pd
from scipy.sparse import issparse
from pandas.api.types import is_categorical_dtype

import matplotlib.cm as cm
from matplotlib.testing import setup
from matplotlib.testing.compare import compare_images

import cellrank as cr
from cellrank.ul.models import GAMR
from cellrank.tl._constants import AbsProbKey
from cellrank.tl.estimators import GPCCA, CFLARE

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
scv.set_figure_params(transparent=True)


try:
    from importlib_metadata import version as get_version
except ImportError:
    # >=Python3.8
    from importlib.metadata import version as get_version

scvelo_paga_skip = pytest.mark.skipif(
    version.parse(get_version(scv.__name__)) < version.parse("0.1.26.dev189+gc441c72"),
    reason="scVelo < `0.1.26.dev189+gc441c72` supports new PAGA, including node colors and confidence",
)

del version, get_version


def compare(
    *,
    kind: str = "adata",
    dirname: Union[str, Path] = None,
    tol: int = TOL,
) -> Callable:
    def _compare_images(
        expected_path: Union[str, Path], actual_path: Union[str, Path]
    ) -> None:
        resize_images_to_same_sizes(expected_path, actual_path)
        res = compare_images(expected_path, actual_path, tol=tol)
        assert res is None, res

    # TODO: refactor (we can remove the prefix from scvelo
    def _prepare_fname(func: Callable) -> Tuple[str, str]:
        fpath = f"{func.__name__.replace('test_', '')}"
        # scvelo saves figures as pdf
        return fpath, str(fpath[7:] + ".png" if fpath.startswith("scvelo_") else fpath)

    def _assert_equal(fpath: str) -> None:
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
            adata, gpcca = adata_gpcca_fwd
            fpath, path = _prepare_fname(func)

            func(self, adata if kind == "adata" else gpcca, path)

            _assert_equal(fpath)

        return decorator

    def compare_gpcca_bwd(func: Callable) -> Callable:
        def decorator(self, adata_gpcca_bwd) -> None:
            adata, gpcca = adata_gpcca_bwd
            fpath, path = _prepare_fname(func)

            func(self, adata, path)

            _assert_equal(fpath)

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

    def compare_gamr(func: Callable):
        def decorator(self, gamr_model: GAMR):
            path, fpath = _prepare_fname(func)

            func(self, gamr_model, path)

            _assert_equal(fpath)

        assert kind == "gamr", "Function `compare_gamr` only supports `kind='gamr'`."

        return decorator

    if kind not in ("adata", "cflare", "gpcca", "lineage", "bwd", "gamr"):
        raise ValueError(
            f"Invalid kind `{kind!r}`. Valid options are: `['adata', 'cflare', 'gpcca', 'lineage', 'bwd', 'gamr']`."
        )

    if kind == "adata":
        # `kind='adata'` - don't changes this, otherwise some tests in `TestHighLvlStates` are meaningless
        return compare_gpcca_fwd
    if kind == "cflare":
        return compare_cflare_fwd
    if kind == "gpcca":
        return compare_gpcca_fwd
    if kind == "lineage":
        return compare_lineage
    if kind == "bwd":
        return compare_gpcca_bwd
    if kind == "gamr":
        return compare_gamr

    raise NotImplementedError(f"Invalid kind `{kind!r}`.")


class TestClusterFates:
    @compare()
    def test_bar(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata, cluster_key="clusters", mode="bar", dpi=DPI, save=fpath
        )

    @compare(kind="bwd")
    def test_bar_bwd(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata,
            cluster_key="clusters",
            backward=True,
            mode="bar",
            dpi=DPI,
            save=fpath,
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
    def test_violin_no_cluster_key(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(adata, mode="violin", cluster_key=None, dpi=DPI, save=fpath)

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
                adata,
                cluster_key="clusters",
                mode="foobar",
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
    def test_mode_heatmap_format(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata,
            cluster_key="clusters",
            mode="heatmap",
            fmt=".10f",
            dpi=DPI,
            save=fpath,
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
            xrot=45,
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

    @compare()
    def test_mode_clustermap_format(self, adata: AnnData, fpath: str):
        cr.pl.cluster_fates(
            adata,
            cluster_key="clusters",
            mode="clustermap",
            fmt=".10f",
            dpi=DPI,
            save=fpath,
        )


class TestClusterLineage:
    @compare()
    def test_cluster_lineage(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.cluster_lineage(
            adata,
            model,
            GENES[:10],
            "1",
            random_state=0,
            time_key="latent_time",
            dpi=DPI,
            save=fpath,
        )

    @compare(kind="bwd")
    def test_cluster_lineage_bwd(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.cluster_lineage(
            adata,
            model,
            GENES[:10],
            "0",
            random_state=0,
            backward=True,
            time_key="latent_time",
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_cluster_lineage_raw(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.cluster_lineage(
            adata,
            model,
            RAW_GENES[:5],
            "1",
            random_state=0,
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
            "1",
            random_state=0,
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
            "1",
            random_state=0,
            time_key="latent_time",
            data_key="Ms",
            norm=False,
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_cluster_lineage_random_state(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.cluster_lineage(
            adata,
            model,
            GENES[:10],
            "1",
            time_key="latent_time",
            random_state=42,
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_cluster_lineage_leiden(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.cluster_lineage(
            adata,
            model,
            GENES[:10],
            "1",
            random_state=0,
            time_key="latent_time",
            use_leiden=True,
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_cluster_lineage_2_failed_genes(self, adata: AnnData, fpath: str):
        fm = create_failed_model(adata)
        cr.pl.cluster_lineage(
            adata,
            {GENES[0]: fm, GENES[5]: fm, "*": fm.model},
            GENES[:10],
            "1",
            random_state=0,
            time_key="latent_time",
            key="foobar",
            dpi=DPI,
            save=fpath,
        )

        assert isinstance(adata.uns["foobar"], AnnData)
        assert adata.uns["foobar"].shape == (8, 200)

    def test_cluster_lineage_returns_fitted_models(self, adata_cflare: AnnData):
        fm = create_failed_model(adata_cflare)
        models = cr.pl.cluster_lineage(
            adata_cflare,
            {GENES[0]: fm, "*": fm.model},
            GENES[:10],
            "1",
            random_state=0,
            time_key="latent_time",
            return_models=True,
        )

        models = pd.DataFrame(models).T
        np.testing.assert_array_equal(models.index, GENES[:10])
        np.testing.assert_array_equal(models.columns, ["1"])
        assert isinstance(models.loc[GENES[0], "1"], cr.ul.models.FailedModel)

        mask = models.astype(bool)
        assert not mask.loc[GENES[0], "1"]
        mask.loc[GENES[0], "1"] = True

        assert np.all(mask)

    def test_cluster_lineage_random_state_same_pca(self, adata_cflare: AnnData):
        model = create_model(adata_cflare)
        cr.pl.cluster_lineage(
            adata_cflare,
            model,
            GENES[:10],
            "1",
            time_key="latent_time",
            random_state=42,
            key="foo",
        )

        cr.pl.cluster_lineage(
            adata_cflare,
            model,
            GENES[:10],
            "1",
            time_key="latent_time",
            random_state=42,
            key="bar",
        )

        np.allclose(
            adata_cflare.uns["foo"].obsm["X_pca"], adata_cflare.uns["bar"].obsm["X_pca"]
        )

    def test_cluster_lineage_writes(self, adata_cflare: AnnData):
        model = create_model(adata_cflare)
        cr.pl.cluster_lineage(adata_cflare, model, GENES[:10], "0", n_test_points=200)

        assert isinstance(adata_cflare.uns["lineage_0_trend"], AnnData)
        assert adata_cflare.uns["lineage_0_trend"].shape == (10, 200)
        assert is_categorical_dtype(adata_cflare.uns["lineage_0_trend"].obs["clusters"])

    def test_cluster_lineage_key(self, adata_cflare: AnnData):
        model = create_model(adata_cflare)
        cr.pl.cluster_lineage(
            adata_cflare, model, GENES[:10], "0", n_test_points=200, key="foobar"
        )

        assert isinstance(adata_cflare.uns["foobar"], AnnData)
        assert adata_cflare.uns["foobar"].shape == (10, 200)
        assert is_categorical_dtype(adata_cflare.uns["foobar"].obs["clusters"])


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

    @compare(kind="bwd", dirname="heatmap_lineages_bwd")
    def test_heatmap_lineages_bwd(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            GENES[:10],
            backward=True,
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

    @compare()
    def test_heatmap_cluster_genes(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            GENES[:10],
            lineages="1",
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
            cbar=False,
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
            cbar=False,
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

    @compare()
    def test_heatmap_cluster_no_scale(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            GENES[:5],
            lineages="1",
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

    @pytest.mark.skip("Hangs using pytest-xdist")
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

    @pytest.mark.skip("Hangs using pytest-xdist")
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

    @compare()
    def test_heatmap_show_dendrogram(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.heatmap(
            adata,
            model,
            GENES[:10],
            mode="lineages",
            lineages="1",
            time_key="latent_time",
            cluster_genes=True,
            dendrogram=True,
            dpi=DPI,
            save=fpath,
        )

    @compare(dirname="heatmap_lineages_1_lineage_failed")
    def test_heatmap_lineages_1_lineage_failed(self, adata: AnnData, fpath: str):
        fm = create_failed_model(adata)
        cr.pl.heatmap(
            adata,
            {g: {"0": fm, "*": fm.model} for g in GENES[:10]},
            GENES[:10],
            mode="lineages",
            time_key="latent_time",
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_heatmap_genes_1_gene_failed(self, adata: AnnData, fpath: str):
        fm = create_failed_model(adata)
        cr.pl.heatmap(
            adata,
            {GENES[0]: fm, "*": fm.model},
            GENES[:10],
            mode="genes",
            time_key="latent_time",
            dpi=DPI,
            save=fpath,
        )


class TestHeatmapReturns:
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

    def test_heatmap_lineages_return_models(self, adata_cflare: AnnData):
        model = create_model(adata_cflare)
        models = cr.pl.heatmap(
            adata_cflare,
            model,
            GENES[:10],
            mode="lineages",
            time_key="latent_time",
            return_models=True,
            dpi=DPI,
        )

        models = pd.DataFrame(models).T
        np.testing.assert_array_equal(models.index, GENES[:10])
        np.testing.assert_array_equal(
            models.columns, adata_cflare.obsm[AbsProbKey.FORWARD.s].names
        )
        assert np.all(models.astype(bool))

    def test_heatmap_lineages_return_models_and_genes(self, adata_cflare: AnnData):
        model = create_model(adata_cflare)
        models, df = cr.pl.heatmap(
            adata_cflare,
            model,
            GENES[:10],
            mode="lineages",
            time_key="latent_time",
            return_models=True,
            return_genes=True,
            dpi=DPI,
        )

        lnames = adata_cflare.obsm[AbsProbKey.FORWARD.s].names

        models = pd.DataFrame(models).T
        np.testing.assert_array_equal(models.index, GENES[:10])
        np.testing.assert_array_equal(models.columns, lnames)
        assert np.all(models.astype(bool))

        assert isinstance(df, pd.DataFrame)
        np.testing.assert_array_equal(df.columns, lnames)
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

    def test_heatmap_genes_return_no_genes(self, adata_cflare: AnnData):
        model = create_model(adata_cflare)
        df = cr.pl.heatmap(
            adata_cflare,
            model,
            GENES[:10],
            mode="genes",
            time_key="latent_time",
            cluster_genes=True,
            dendrogram=True,
            return_genes=True,
            dpi=DPI,
        )

        assert df is None


class TestGeneTrend:
    @compare()
    def test_trends(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.gene_trends(
            adata,
            model,
            GENES[:3],
            data_key="Ms",
            dpi=DPI,
            save=fpath,
        )

    @compare(kind="bwd")
    def test_trends_bwd(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.gene_trends(
            adata,
            model,
            GENES[:3],
            backward=True,
            data_key="Ms",
            dpi=DPI,
            save=fpath,
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
            adata,
            model,
            GENES[:3],
            data_key="Ms",
            same_plot=True,
            dpi=DPI,
            save=fpath,
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
            adata,
            model,
            GENES[:3],
            data_key="Ms",
            sharey="row",
            dpi=DPI,
            save=fpath,
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
            cbar=False,
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
    def test_trends_lw(self, adata: AnnData, fpath: str):
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
    def test_trends_suptitle(self, adata: AnnData, fpath: str):
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
    def test_trends_time_range(self, adata: AnnData, fpath: str):
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
    def test_trends_perc(self, adata: AnnData, fpath: str):
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

    @compare()
    def test_trends_perc_per_lineage(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.gene_trends(
            adata,
            model,
            GENES[:3],
            figsize=(5, 5),
            data_key="Ms",
            same_plot=False,
            perc=[(0, 50), (5, 95), (50, 100)],
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_trends_time_key(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.gene_trends(
            adata,
            model,
            GENES[:10],
            data_key="Ms",
            same_plot=False,
            time_key="dpt_pseudotime",
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_trends_show_lineage_ignores_no_transpose(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.gene_trends(
            adata,
            model,
            GENES[:5],
            transpose=False,
            data_key="Ms",
            same_plot=True,
            plot_kwargs=dict(lineage_probability=True),
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_trends_show_lineage_same_plot(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.gene_trends(
            adata,
            model,
            GENES[:5],
            transpose=True,
            data_key="Ms",
            same_plot=True,
            plot_kwargs=dict(lineage_probability=True),
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_trends_show_lineage_diff_plot(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.gene_trends(
            adata,
            model,
            GENES[0],
            data_key="Ms",
            same_plot=False,
            transpose=True,
            plot_kwargs=dict(lineage_probability=True),
            figsize=(5, 5),
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_trends_show_lineage_ci(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.gene_trends(
            adata,
            model,
            GENES[0],
            data_key="Ms",
            same_plot=True,
            transpose=True,
            plot_kwargs=dict(
                lineage_probability=True, lineage_probability_conf_int=True
            ),
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_trends_time_key_del_latent_time(self, adata: AnnData, fpath: str):
        # this ensures that the callback passes the correct values
        del adata.obs["latent_time"]
        assert "latent_time" not in adata.obs

        model = create_model(adata)
        cr.pl.gene_trends(
            adata,
            model,
            GENES[:10],
            data_key="Ms",
            same_plot=False,
            time_key="dpt_pseudotime",
            dpi=DPI,
            save=fpath,
        )

    def test_invalid_time_key(self, adata_cflare: AnnData):
        model = create_model(adata_cflare)
        with pytest.raises(KeyError):
            cr.pl.gene_trends(
                adata_cflare,
                model,
                GENES[:10],
                data_key="Ms",
                same_plot=False,
                time_key="foobar",
            )

    def test_all_models_failed(self, adata_cflare: AnnData):
        fm = create_failed_model(adata_cflare)
        with pytest.raises(RuntimeError):
            cr.pl.gene_trends(
                adata_cflare,
                fm,
                GENES[:10],
                data_key="Ms",
                mode="lineages",
                time_key="latent_time",
                dpi=DPI,
            )

    def test_return_models_no_failures(self, adata_cflare: AnnData):
        model = create_model(adata_cflare)
        models = cr.pl.gene_trends(
            adata_cflare,
            model,
            GENES[:10],
            data_key="Ms",
            lineages=["0", "1"],
            time_key="latent_time",
            dpi=DPI,
            return_models=True,
        )

        models = pd.DataFrame(models).T
        np.testing.assert_array_equal(models.index, GENES[:10])
        np.testing.assert_array_equal(models.columns, [str(i) for i in range(2)])
        assert np.all(models.astype(bool))

    def test_return_models_with_failures(self, adata_cflare: AnnData):
        fm = create_failed_model(adata_cflare)
        models = cr.pl.gene_trends(
            adata_cflare,
            {GENES[0]: {"0": fm, "*": fm.model}, "*": fm.model},
            GENES[:10],
            lineages=["0", "1"],
            time_key="latent_time",
            dpi=DPI,
            return_models=True,
        )

        models = pd.DataFrame(models).T
        np.testing.assert_array_equal(models.index, GENES[:10])
        np.testing.assert_array_equal(models.columns, [str(i) for i in range(2)])
        assert isinstance(models.loc[GENES[0], "0"], cr.ul.models.FailedModel)

        mask = models.astype(bool)
        assert not mask.loc[GENES[0], "0"]
        mask.loc[GENES[0], "0"] = True

        assert np.all(mask)

    @compare()
    def test_all_models_for_1_gene_failed(self, adata: AnnData, fpath: str):
        fm = create_failed_model(adata)
        cr.pl.gene_trends(
            adata,
            {GENES[0]: fm, "*": fm.model},
            GENES[:3],
            figsize=(5, 5),
            data_key="Ms",
            time_key="latent_time",
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_all_models_for_1_lineage_failed(self, adata: AnnData, fpath: str):
        fm = create_failed_model(adata)
        cr.pl.gene_trends(
            adata,
            {g: {"0": fm, "*": fm.model} for g in GENES[:3]},
            GENES[:3],
            figsize=(5, 5),
            data_key="Ms",
            time_key="latent_time",
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_all_models_for_1_gene_failed_same_plot(self, adata: AnnData, fpath: str):
        fm = create_failed_model(adata)
        cr.pl.gene_trends(
            adata,
            {GENES[0]: fm, "*": fm.model},
            GENES[:10],
            data_key="Ms",
            time_key="latent_time",
            same_plot=True,
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_failed_only_main_diagonal(self, adata: AnnData, fpath: str):
        fm = create_failed_model(adata)
        cr.pl.gene_trends(
            adata,
            {g: {str(ln): fm.model, "*": fm} for ln, g in enumerate(GENES[:3])},
            GENES[:3],
            lineages=["0", "1", "2"],
            data_key="Ms",
            time_key="latent_time",
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_failed_only_off_diagonal(self, adata: AnnData, fpath: str):
        fm = create_failed_model(adata)
        cr.pl.gene_trends(
            adata,
            {
                g: {str(ln): fm.model, "*": fm}
                for ln, g in zip(range(3)[::-1], GENES[:3])
            },
            GENES[:3],
            data_key="Ms",
            time_key="latent_time",
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_transpose(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.gene_trends(
            adata,
            model,
            GENES[:4],
            transpose=True,
            data_key="Ms",
            time_key="latent_time",
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_transpose_same_plot(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        cr.pl.gene_trends(
            adata,
            model,
            GENES[:3],
            transpose=True,
            same_plot=True,
            data_key="Ms",
            time_key="latent_time",
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_transpose_all_models_for_1_gene_failed(self, adata: AnnData, fpath: str):
        fm = create_failed_model(adata)
        cr.pl.gene_trends(
            adata,
            {GENES[0]: fm, "*": fm.model},
            GENES[:10],
            transpose=True,
            time_key="latent_time",
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_transpose_all_models_for_1_lineage_failed(
        self, adata: AnnData, fpath: str
    ):
        fm = create_failed_model(adata)
        cr.pl.gene_trends(
            adata,
            {g: {"0": fm, "*": fm.model} for g in GENES[:10]},
            GENES[:10],
            transpose=True,
            data_key="Ms",
            time_key="latent_time",
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_transpose_failed_only_off_diagonal(self, adata: AnnData, fpath: str):
        fm = create_failed_model(adata)
        cr.pl.gene_trends(
            adata,
            {
                g: {str(ln): fm.model, "*": fm}
                for ln, g in zip(range(3)[::-1], GENES[:3])
            },
            GENES[:3],
            transpose=True,
            data_key="Ms",
            time_key="latent_time",
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_transpose_all_models_for_1_lineage_failed_same_plot(
        self, adata: AnnData, fpath: str
    ):
        fm = create_failed_model(adata)
        cr.pl.gene_trends(
            adata,
            {g: {"0": fm, "*": fm.model} for g in GENES[:10]},
            GENES[:10],
            transpose=True,
            same_plot=True,
            data_key="Ms",
            time_key="latent_time",
            dpi=DPI,
            save=fpath,
        )


class TestGraph:
    @compare()
    def test_graph(self, adata: AnnData, fpath: str):
        cr.pl.graph(
            adata, "T_fwd", ixs=range(10), edge_use_curved=False, dpi=DPI, save=fpath
        )

    @compare(kind="bwd")
    def test_graph_bwd(self, adata: AnnData, fpath: str):
        cr.pl.graph(
            adata, "T_bwd", ixs=range(10), edge_use_curved=False, dpi=DPI, save=fpath
        )

    @compare()
    def test_graph_layout(self, adata: AnnData, fpath: str):
        cr.pl.graph(
            adata,
            "T_fwd",
            ixs=range(10),
            edge_use_curved=False,
            layout="umap",
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_graph_title(self, adata: AnnData, fpath: str):
        cr.pl.graph(
            adata,
            "T_fwd",
            ixs=range(10),
            title="foo bar baz quux",
            edge_use_curved=False,
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_graph_titles(self, adata: AnnData, fpath: str):
        cr.pl.graph(
            adata,
            "T_fwd",
            ixs=range(10),
            edge_use_curved=False,
            keys=["incoming", "self_loops"],
            title=["foo", "bar"],
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_graph_keys(self, adata: AnnData, fpath: str):
        cr.pl.graph(
            adata,
            "T_fwd",
            ixs=range(10),
            edge_use_curved=False,
            keys=("outgoing", "self_loops"),
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_graph_edge_weight_scale(self, adata: AnnData, fpath: str):
        cr.pl.graph(
            adata,
            "T_fwd",
            ixs=range(10),
            edge_use_curved=False,
            edge_weight_scale=100,
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_graph_show_arrows(self, adata: AnnData, fpath: str):
        cr.pl.graph(
            adata,
            "T_fwd",
            ixs=range(15),
            edge_use_curved=False,
            arrows=False,
            edge_weight_scale=100,
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_graph_curved_edges(self, adata: AnnData, fpath: str):
        cr.pl.graph(
            adata, "T_fwd", ixs=range(10), edge_use_curved=True, dpi=DPI, save=fpath
        )

    @compare()
    def test_graph_labels(self, adata: AnnData, fpath: str):
        cr.pl.graph(
            adata,
            "T_fwd",
            ixs=range(10),
            edge_use_curved=False,
            labels=range(10),
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_graph_cmap(self, adata: AnnData, fpath: str):
        cr.pl.graph(
            adata,
            "T_fwd",
            ixs=range(10),
            edge_use_curved=False,
            cont_cmap=cm.inferno,
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_graph_top_n_edges_incoming(self, adata: AnnData, fpath: str):
        cr.pl.graph(
            adata,
            "T_fwd",
            ixs=range(10),
            edge_use_curved=False,
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
            edge_use_curved=False,
            top_n_edges=(2, False, "outgoing"),
            edge_weight_scale=100,
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_graph_edge_normalize(self, adata: AnnData, fpath: str):
        cr.pl.graph(
            adata,
            "T_fwd",
            ixs=range(10),
            edge_use_curved=False,
            edge_normalize=True,
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_graph_edge_reductions(self, adata: AnnData, fpath: str):
        cr.pl.graph(
            adata,
            "T_fwd",
            ixs=range(10),
            edge_use_curved=False,
            edge_reductions=np.max,
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_graph_edge_reductions_restriction_incoming(
        self, adata: AnnData, fpath: str
    ):
        cr.pl.graph(
            adata,
            "T_fwd",
            ixs=range(10),
            keys="incoming",
            edge_use_curved=False,
            edge_reductions_restrict_to_ixs=range(20, 40),
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_graph_edge_reductions_restriction_outgoing(
        self, adata: AnnData, fpath: str
    ):
        cr.pl.graph(
            adata,
            "T_fwd",
            ixs=range(10),
            keys="outgoing",
            edge_use_curved=False,
            edge_reductions_restrict_to_ixs=range(20, 40),
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_graph_categorical_key(self, adata: AnnData, fpath: str):
        cr.pl.graph(
            adata,
            "T_fwd",
            ixs=range(10),
            edge_use_curved=False,
            keys="clusters",
            keylocs="obs",
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_graph_filter_edges(self, adata: AnnData, fpath: str):
        cr.pl.graph(
            adata,
            "T_fwd",
            ixs=range(10),
            filter_edges=(0.25, 0.5),
            edge_use_curved=False,
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_graph_dict_layout(self, adata: AnnData, fpath: str):
        cr.pl.graph(
            adata,
            "T_fwd",
            ixs=range(10),
            layout={i: (i, i) for i in range(10)},
            edge_use_curved=False,
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_graph_networkx_layout(self, adata: AnnData, fpath: str):
        import networkx as nx

        cr.pl.graph(
            adata,
            "T_fwd",
            ixs=range(10),
            layout=nx.layout.kamada_kawai_layout,
            edge_use_curved=False,
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_graph_precomputed_layour_pca(self, adata: AnnData, fpath: str):
        cr.pl.graph(
            adata,
            "T_fwd",
            ixs=range(10),
            layout="pca",
            edge_use_curved=False,
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
    def test_mc_real_spectrum_hide_xticks(self, mc: CFLARE, fpath: str):
        mc.plot_spectrum(real_only=True, show_all_xticks=False, dpi=DPI, save=fpath)

    @compare(kind="cflare")
    def test_mc_real_spectrum_hide_eigengap(self, mc: CFLARE, fpath: str):
        mc.plot_spectrum(real_only=True, show_eigengap=False, dpi=DPI, save=fpath)

    @compare(kind="cflare")
    def test_mc_spectrum_title(self, mc: CFLARE, fpath: str):
        mc.plot_spectrum(title="foobar", real_only=False, dpi=DPI, save=fpath)

    @compare(kind="cflare")
    def test_mc_marker(self, mc: CFLARE, fpath: str):
        mc.plot_spectrum(dpi=DPI, marker="X", save=fpath)

    @compare(kind="cflare")
    def test_mc_kwargs_linewidths(self, mc: CFLARE, fpath: str):
        mc.plot_spectrum(dpi=DPI, linewidths=20, save=fpath)

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
        mc.plot_terminal_states(dpi=DPI, save=fpath)

    @compare(kind="cflare")
    def test_scvelo_final_states_clusters(self, mc: CFLARE, fpath: str):
        mc.plot_terminal_states(cluster_key="clusters", dpi=DPI, save=fpath)

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
    def test_gpcca_real_spectrum_hide_eigengap(self, mc: GPCCA, fpath: str):
        mc.plot_spectrum(real_only=True, show_eigengap=False, dpi=DPI, save=fpath)

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
    def test_gpcca_coarse_T_stat_init_dist(self, mc: GPCCA, fpath: str):
        mc.plot_coarse_T(
            show_initial_dist=True, show_stationary_dist=True, dpi=DPI, save=fpath
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
        mc.plot_macrostates(dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_meta_states_lineages(self, mc: GPCCA, fpath: str):
        mc.plot_macrostates(lineages=["0"], dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_meta_states_discrete(self, mc: GPCCA, fpath: str):
        mc.plot_macrostates(discrete=True, dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_meta_states_cluster_key(self, mc: GPCCA, fpath: str):
        mc.plot_macrostates(cluster_key="clusters", dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_meta_states_no_same_plot(self, mc: GPCCA, fpath: str):
        mc.plot_macrostates(same_plot=False, dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_meta_states_cmap(self, mc: GPCCA, fpath: str):
        mc.plot_macrostates(cmap=cm.inferno, same_plot=False, dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_meta_states_title(self, mc: GPCCA, fpath: str):
        mc.plot_macrostates(title="foobar", dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_meta_states_time(self, mc: GPCCA, fpath: str):
        mc.plot_macrostates(mode="time", dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_final_states(self, mc: GPCCA, fpath: str):
        mc.plot_terminal_states(dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_final_states_lineages(self, mc: GPCCA, fpath: str):
        mc.plot_terminal_states(lineages=["0"], dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_final_states_discrete(self, mc: GPCCA, fpath: str):
        mc.plot_terminal_states(discrete=True, dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_final_states_cluster_key(self, mc: GPCCA, fpath: str):
        mc.plot_terminal_states(cluster_key="clusters", dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_final_states_no_same_plot(self, mc: GPCCA, fpath: str):
        mc.plot_terminal_states(same_plot=False, dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_final_states_cmap(self, mc: GPCCA, fpath: str):
        mc.plot_terminal_states(cmap=cm.inferno, same_plot=False, dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_final_states_title(self, mc: GPCCA, fpath: str):
        mc.plot_terminal_states(title="foobar", dpi=DPI, save=fpath)

    @compare(kind="gpcca")
    def test_scvelo_gpcca_final_states_time(self, mc: GPCCA, fpath: str):
        mc.plot_terminal_states(mode="time", dpi=DPI, save=fpath)

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


class TestHighLvlStates:
    @compare()
    def test_scvelo_terminal_states_disc(self, adata: AnnData, fpath: str):
        cr.pl.terminal_states(adata, discrete=True, dpi=DPI, save=fpath)

    @compare(kind="bwd")
    def test_scvelo_initial_states_disc(self, adata: AnnData, fpath: str):
        cr.pl.initial_states(adata, discrete=True, dpi=DPI, save=fpath)

    # only matters when kind='adata' was computed using GPCCA
    @compare()
    def test_scvelo_terminal_states_cont(self, adata: AnnData, fpath: str):
        cr.pl.terminal_states(adata, discrete=False, dpi=DPI, save=fpath)

    @compare()
    def test_scvelo_terminal_disc_same_subset(self, adata: AnnData, fpath: str):
        cr.pl.terminal_states(
            adata, discrete=True, same_plot=True, states="0", dpi=DPI, save=fpath
        )

    @compare()
    def test_scvelo_terminal_disc_not_same_subset(self, adata: AnnData, fpath: str):
        cr.pl.terminal_states(
            adata, discrete=True, same_plot=False, states="0", dpi=DPI, save=fpath
        )

    @compare()
    def test_scvelo_terminal_cont_same_subset(self, adata: AnnData, fpath: str):
        cr.pl.terminal_states(
            adata, discrete=False, same_plot=True, states="0", dpi=DPI, save=fpath
        )

    @compare()
    def test_scvelo_terminal_cont_not_same_subset(self, adata: AnnData, fpath: str):
        cr.pl.terminal_states(
            adata, discrete=False, same_plot=False, states="0", dpi=DPI, save=fpath
        )

    @compare()
    def test_scvelo_terminal_diff_plot(self, adata: AnnData, fpath: str):
        cr.pl.terminal_states(adata, same_plot=False, dpi=DPI, save=fpath)

    @compare()
    def test_scvelo_terminal_diff_plot_titles(self, adata: AnnData, fpath: str):
        cr.pl.terminal_states(
            adata, same_plot=False, title=["foo", "bar"] * 10, dpi=DPI, save=fpath
        )

    @compare()
    def test_scvelo_terminal_cluster_key_discrete(self, adata: AnnData, fpath: str):
        cr.pl.terminal_states(
            adata, discrete=True, cluster_key="clusters", dpi=DPI, save=fpath
        )

    @compare()
    def test_scvelo_terminal_time_mode(self, adata: AnnData, fpath: str):
        # only works in continuous mode
        cr.pl.terminal_states(
            adata,
            discrete=False,
            mode="time",
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_scvelo_terminal_time_mode_subset(self, adata: AnnData, fpath: str):
        # only works in continuous mode
        cr.pl.terminal_states(
            adata,
            states="0",
            discrete=False,
            mode="time",
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_scvelo_terminal_time_mode_clusters(self, adata: AnnData, fpath: str):
        # only works in continuous mode
        cr.pl.terminal_states(
            adata,
            discrete=False,
            cluster_key="clusters",
            mode="time",
            dpi=DPI,
            save=fpath,
        )


class TestLineage:
    def test_pie(self, lineage: cr.tl.Lineage):
        with pytest.raises(ValueError):
            lineage[:, 0].plot_pie(dpi=DPI)

    @compare(kind="lineage")
    def test_pie(self, lineage: cr.tl.Lineage, fpath: str):
        lineage.plot_pie(np.mean, dpi=DPI, save=fpath)

    @compare(kind="lineage")
    def test_pie_reduction(self, lineage: cr.tl.Lineage, fpath: str):
        lineage.plot_pie(np.var, dpi=DPI, save=fpath)

    @compare(kind="lineage")
    def test_pie_title(self, lineage: cr.tl.Lineage, fpath: str):
        lineage.plot_pie(np.mean, title="FOOBAR", dpi=DPI, save=fpath)

    @compare(kind="lineage")
    def test_pie_t(self, lineage: cr.tl.Lineage, fpath: str):
        lineage.T.plot_pie(np.mean, dpi=DPI, save=fpath)

    @compare(kind="lineage")
    def test_pie_autopct_none(self, lineage: cr.tl.Lineage, fpath: str):
        lineage.T.plot_pie(np.mean, dpi=DPI, save=fpath, autopct=None)

    @compare(kind="lineage")
    def test_pie_legend_loc(self, lineage: cr.tl.Lineage, fpath: str):
        lineage.plot_pie(np.mean, dpi=DPI, save=fpath, legend_loc="best")

    @compare(kind="lineage")
    def test_pie_legend_loc_one(self, lineage: cr.tl.Lineage, fpath: str):
        lineage.plot_pie(np.mean, dpi=DPI, save=fpath, legend_loc=None)

    @compare(kind="lineage")
    def test_pie_legend_kwargs(self, lineage: cr.tl.Lineage, fpath: str):
        lineage.plot_pie(
            np.mean,
            dpi=DPI,
            save=fpath,
            legend_loc="best",
            legend_kwargs={"fontsize": 20},
        )


class TestLineageDrivers:
    @compare()
    def test_scvelo_drivers_n_genes(self, adata: AnnData, fpath: str):
        cr.pl.lineage_drivers(adata, "0", n_genes=5, dpi=DPI, save=fpath)

    @compare(kind="bwd")
    def test_scvelo_drivers_n_genes(self, adata: AnnData, fpath: str):
        cr.pl.lineage_drivers(adata, "0", backward=True, n_genes=5, dpi=DPI, save=fpath)

    @compare()
    def test_scvelo_drivers_cmap(self, adata: AnnData, fpath: str):
        cr.pl.lineage_drivers(adata, "0", cmap="inferno", dpi=DPI, save=fpath)


class TestModel:
    @compare()
    def test_model_default(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        model.prepare(adata.var_names[0], "1")
        model.fit().predict()
        model.confidence_interval()
        model.plot(save=fpath, dpi=DPI)

    @compare(kind="bwd")
    def test_model_default_bwd(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        model.prepare(adata.var_names[0], "0", backward=True)
        model.fit().predict()
        model.confidence_interval()
        model.plot(save=fpath, dpi=DPI)

    @compare()
    def test_model_obs_data_key(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        gene = adata.X[:, 0]
        adata.obs["foo"] = gene.A if issparse(gene) else gene

        model.prepare("foo", "1", data_key="obs")
        model.fit().predict()
        model.confidence_interval()
        model.plot(save=fpath, dpi=DPI)

    @compare()
    def test_model_no_lineage(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        model.prepare(adata.var_names[0], None)
        model.fit().predict()
        model.confidence_interval()
        model.plot(save=fpath, dpi=DPI)

    @compare()
    def test_model_no_lineage_show_lin_probs(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        model.prepare(adata.var_names[0], None)
        model.fit().predict()
        model.plot(save=fpath, dpi=DPI, lineage_probability=True)

    @compare()
    def test_model_no_legend(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        model.prepare(adata.var_names[0], "1")
        model.fit().predict()
        model.confidence_interval()
        model.plot(save=fpath, dpi=DPI, loc=None)

    # TODO: parametrize (hide cells, ci)
    @compare()
    def test_model_show_lin_prob_cells_ci(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        model.prepare(adata.var_names[0], "1")
        model.fit().predict()
        model.confidence_interval()
        model.plot(
            save=fpath,
            dpi=DPI,
            hide_cells=False,
            conf_int=True,
            lineage_probability=True,
        )

    @compare()
    def test_model_show_lin_prob_cells_lineage_ci(self, adata: AnnData, fpath: str):
        model = create_model(adata)
        model.prepare(adata.var_names[0], "1")
        model.fit().predict()
        model.confidence_interval()
        model.plot(
            save=fpath,
            dpi=DPI,
            hide_cells=True,
            conf_int=True,
            lineage_probability=True,
            lineage_probability_conf_int=True,
        )


@gamr_skip
class TestGAMR:
    @compare(kind="gamr")
    def test_gamr_default(self, model: GAMR, fpath: str):
        model.prepare(model.adata.var_names[0], "1")
        model.fit().predict()
        model.plot(
            save=fpath,
            dpi=DPI,
        )

    @compare(kind="gamr")
    def test_gamr_ci_50(self, model: GAMR, fpath: str):
        model.prepare(model.adata.var_names[0], "1")
        model.fit().predict(level=0.5)
        model.plot(
            conf_int=True,
            save=fpath,
            dpi=DPI,
        )

    @compare(kind="gamr")
    def test_gamr_no_ci(self, model: GAMR, fpath: str):
        model.prepare(model.adata.var_names[0], "1")
        model.fit().predict(level=None)
        model.plot(
            conf_int=False,
            save=fpath,
            dpi=DPI,
        )

    @compare(kind="gamr")
    def test_gamr_no_cbar(self, model: GAMR, fpath: str):
        model.prepare(model.adata.var_names[0], "1")
        model.fit().predict(level=0.95)
        model.plot(
            cbar=False,
            save=fpath,
            dpi=DPI,
        )

    @compare(kind="gamr")
    def test_gamr_lineage_prob(self, model: GAMR, fpath: str):
        model.prepare(model.adata.var_names[0], "1")
        model.fit().predict(level=0.95)
        model.plot(
            lineage_probability=True,
            lineage_probability_conf_int=True,
            save=fpath,
            dpi=DPI,
        )

    @compare(kind="gamr")
    def test_trends_gam_ci_100(self, model: GAMR, fpath: str):
        cr.pl.gene_trends(
            model.adata,
            model,
            GENES[:3],
            conf_int=1,
            backward=False,
            data_key="Ms",
            dpi=DPI,
            save=fpath,
        )

    @compare(kind="gamr")
    def test_trends_gam_ci_20(self, model: GAMR, fpath: str):
        cr.pl.gene_trends(
            model.adata,
            model,
            GENES[:3],
            conf_int=0.2,
            backward=False,
            data_key="Ms",
            dpi=DPI,
            save=fpath,
        )


class TestComposition:
    @compare()
    def test_composition(self, adata: AnnData, fpath: str):
        cr.pl._utils.composition(adata, "clusters", dpi=DPI, save=fpath)

    @compare()
    def test_composition_kwargs_autopct(self, adata: AnnData, fpath: str):
        cr.pl._utils.composition(
            adata, "clusters", dpi=DPI, save=fpath, autopct="%1.0f%%"
        )


class TestFittedModel:
    @compare()
    def test_fitted_empty_model(self, _adata: AnnData, fpath: str):
        np.random.seed(42)
        fm = cr.ul.models.FittedModel(np.arange(100), np.random.normal(size=100))
        fm.plot(dpi=DPI, save=fpath)

    @compare()
    def test_fitted_model_conf_int(self, _adata: AnnData, fpath: str):
        np.random.seed(43)
        y_test = np.random.normal(size=100)

        fm = cr.ul.models.FittedModel(
            np.arange(100), y_test, conf_int=np.c_[y_test - 1, y_test + 1]
        )
        fm.plot(conf_int=True, dpi=DPI, save=fpath)

    @compare()
    def test_fitted_model_conf_int_no_conf_int_computed(
        self, _adata: AnnData, fpath: str
    ):
        np.random.seed(44)

        fm = cr.ul.models.FittedModel(
            np.arange(100),
            np.random.normal(size=100),
        )
        fm.plot(conf_int=True, dpi=DPI, save=fpath)

    @compare()
    def test_fitted_model_cells_with_weights(self, _adata: AnnData, fpath: str):
        np.random.seed(45)

        fm = cr.ul.models.FittedModel(
            np.arange(100),
            np.random.normal(size=100),
            x_all=np.random.normal(size=200),
            y_all=np.random.normal(size=200),
        )

        fm.plot(hide_cells=False, dpi=DPI, save=fpath)

    @compare()
    def test_fitted_model_weights(self, _adata: AnnData, fpath: str):
        np.random.seed(46)

        fm = cr.ul.models.FittedModel(
            np.arange(100),
            np.random.normal(size=100),
            x_all=np.random.normal(size=200),
            y_all=np.random.normal(size=200),
            w_all=np.random.normal(size=200),
        )

        fm.plot(hide_cells=False, dpi=DPI, save=fpath)

    @compare()
    def test_fitted_ignore_plot_smoothed_lineage(self, _adata: AnnData, fpath: str):
        np.random.seed(47)

        fm = cr.ul.models.FittedModel(
            np.arange(100),
            np.random.normal(size=100),
            x_all=np.random.normal(size=200),
            y_all=np.random.normal(size=200),
            w_all=np.random.normal(size=200),
        )

        fm.plot(
            lineage_probability=True,
            lineage_probability_conf_int=True,
            dpi=DPI,
            save=fpath,
        )

    @compare()
    def test_fitted_gene_trends(self, adata: AnnData, fpath: str):
        np.random.seed(48)

        fm1 = cr.ul.models.FittedModel(
            np.arange(100),
            np.random.normal(size=100),
            x_all=np.random.normal(size=200),
            y_all=np.random.normal(size=200),
            w_all=np.random.normal(size=200),
        )
        fm2 = cr.ul.models.FittedModel(
            np.arange(100),
            np.random.normal(size=100),
            x_all=np.random.normal(size=200),
            y_all=np.random.normal(size=200),
            w_all=np.random.normal(size=200),
        )
        cr.pl.gene_trends(
            adata,
            {GENES[0]: fm1, GENES[1]: fm2},
            GENES[:2],
            data_key="Ms",
            dpi=DPI,
            save=fpath,
        )

    @compare(tol=250)
    def test_fitted_cluster_fates(self, adata: AnnData, fpath: str):
        np.random.seed(49)

        model = cr.ul.models.FittedModel(
            np.arange(100),
            np.random.normal(size=100),
        )
        cr.pl.cluster_lineage(
            adata,
            model,
            GENES[:10],
            "1",
            n_points=100,
            time_key="latent_time",
            random_state=49,
            dpi=DPI,
            save=fpath,
        )

    @compare(dirname="fitted_heatmap")
    def test_fitted_heatmap(self, adata: AnnData, fpath: str):
        np.random.seed(49)

        fm = cr.ul.models.FittedModel(
            np.arange(100),
            np.random.normal(size=100),
        )
        cr.pl.heatmap(
            adata,
            fm,
            GENES[:10],
            mode="lineages",
            time_key="latent_time",
            dpi=DPI,
            save=fpath,
        )


class TestCircularProjection:
    def test_proj_too_few_lineages(self, adata_gpcca_fwd):
        adata, _ = adata_gpcca_fwd
        lineages = adata.obsm[AbsProbKey.FORWARD.s].names[:2]

        with pytest.raises(ValueError, match=r"Expected at least `3` lineages"):
            cr.pl.circular_projection(
                adata, keys=["clusters", "clusters"], lineages=lineages
            )

    @compare()
    def test_proj_duplicate_keys(self, adata: AnnData, fpath: str):
        cr.pl.circular_projection(
            adata, keys=["clusters", "clusters"], dpi=DPI, save=fpath
        )

        key = "X_fate_simplex_fwd"
        assert key in adata.obsm
        assert isinstance(adata.obsm[key], np.ndarray)
        assert adata.obsm[key].shape[1] == 2

    @compare()
    def test_proj_key_added(self, adata: AnnData, fpath: str):
        key = "foo"
        cr.pl.circular_projection(
            adata, keys=adata.var_names[0], key_added=key, dpi=DPI, save=fpath
        )

        assert key in adata.obsm
        assert isinstance(adata.obsm[key], np.ndarray)
        assert adata.obsm[key].shape[1] == 2

    @compare()
    def test_proj_hide_edges(self, adata: AnnData, fpath: str):
        cr.pl.circular_projection(
            adata, keys="dpt_pseudotime", show_edges=False, dpi=DPI, save=fpath
        )

    @compare()
    def test_proj_dont_normalize_by_mean(self, adata: AnnData, fpath: str):
        cr.pl.circular_projection(
            adata, keys="clusters", normalize_by_mean=False, dpi=DPI, save=fpath
        )

    @compare()
    def test_proj_use_raw(self, adata: AnnData, fpath: str):
        cr.pl.circular_projection(
            adata, keys=adata.raw.var_names[0], use_raw=True, dpi=DPI, save=fpath
        )

    @compare()
    def test_proj_ncols(self, adata: AnnData, fpath: str):
        cr.pl.circular_projection(
            adata, keys=adata.var_names[:2], ncols=1, dpi=DPI, save=fpath
        )

    @compare()
    def test_proj_labelrot(self, adata: AnnData, fpath: str):
        cr.pl.circular_projection(
            adata, keys="clusters", labelrot="default", dpi=DPI, save=fpath
        )

    @compare()
    def test_proj_labeldistance(self, adata: AnnData, fpath: str):
        cr.pl.circular_projection(
            adata, keys="clusters", labeldistance=1.5, dpi=DPI, save=fpath
        )

    @compare()
    def test_proj_text_kwargs(self, adata: AnnData, fpath: str):
        cr.pl.circular_projection(
            adata, keys="clusters", text_kwargs={"size": 20}, dpi=DPI, save=fpath
        )

    @compare()
    def test_proj_default_ordering(self, adata: AnnData, fpath: str):
        cr.pl.circular_projection(
            adata, keys="clusters", lineage_order="default", dpi=DPI, save=fpath
        )

    @compare()
    def test_proj_extra_keys(self, adata: AnnData, fpath: str):
        cr.pl.circular_projection(
            adata, keys=["priming_direction", "priming_degree"], dpi=DPI, save=fpath
        )

        assert "priming_direction_fwd" in adata.obs
        assert is_categorical_dtype(adata.obs["priming_direction_fwd"])
        assert "priming_degree_fwd" in adata.obs

    @compare()
    def test_proj_scvelo_kwargs(self, adata: AnnData, fpath: str):
        cr.pl.circular_projection(
            adata, keys="clusters", legend_loc="upper right", dpi=DPI, save=fpath
        )

    @compare()
    def test_proj_no_cbar(self, adata: AnnData, fpath: str):
        cr.pl.circular_projection(
            adata, keys=adata.var_names[0], colorbar=False, dpi=DPI, save=fpath
        )
