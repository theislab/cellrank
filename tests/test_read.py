# -*- coding: utf-8 -*-
from typing import Tuple, Callable
from pathlib import Path

import pytest

import scanpy as sc
from anndata import AnnData

import numpy as np

import cellrank as cr
from cellrank.tl._colors import _create_categorical_colors
from cellrank.tl._lineage import Lineage
from cellrank.tl._constants import AbsProbKey, _colors, _lin_names
from cellrank.tl.estimators._cflare import CFLARE


def test_fwd():
    def wrapper(func: Callable) -> Callable:
        def decorator(self, adata_cflare_fwd: Tuple[AnnData, CFLARE], tmpdir):
            adata, _ = adata_cflare_fwd
            adata = adata.copy()

            dirname = func.__name__
            path = tmpdir.mkdir(dirname).join("tmp.h5ad")

            return func(
                self,
                adata,
                path,
                str(AbsProbKey.FORWARD),
                adata.obsm[str(AbsProbKey.FORWARD)].shape[1],
            )

        return decorator

    return wrapper


class TestRead:
    @test_fwd()
    def test_no_lineage(self, adata: AnnData, path: Path, lin_key: str, _: int):
        del adata.obsm[lin_key]

        sc.write(path, adata)
        adata_new = cr.read(path)

        assert adata_new is not adata  # sanity check
        assert lin_key not in adata_new.obsm.keys()

    @test_fwd()
    def test_no_names(self, adata: AnnData, path: Path, lin_key: str, n_lins: int):
        names_key = _lin_names(lin_key)
        del adata.uns[names_key]

        sc.write(path, adata)
        adata_new = cr.read(path)
        lins = adata_new.obsm[lin_key]

        assert isinstance(lins, Lineage)
        np.testing.assert_array_equal(
            lins.names, [f"Lineage {i}" for i in range(n_lins)]
        )
        np.testing.assert_array_equal(lins.names, adata_new.uns[names_key])

    @test_fwd()
    def test_no_colors(self, adata: AnnData, path: Path, lin_key: str, n_lins: int):
        colors_key = _colors(lin_key)
        del adata.uns[colors_key]

        sc.write(path, adata)
        adata_new = cr.read(path)
        lins = adata_new.obsm[lin_key]

        assert isinstance(lins, Lineage)
        np.testing.assert_array_equal(lins.colors, _create_categorical_colors(n_lins))
        np.testing.assert_array_equal(lins.colors, adata_new.uns[colors_key])

    @test_fwd()
    def test_wrong_names_length(
        self, adata: AnnData, path: Path, lin_key: str, n_lins: int
    ):
        names_key = _lin_names(lin_key)
        adata.uns[names_key] = list(adata.uns[names_key])
        adata.uns[names_key] += ["foo", "bar", "baz"]

        sc.write(path, adata)
        adata_new = cr.read(path)
        lins = adata_new.obsm[lin_key]

        assert isinstance(lins, Lineage)
        np.testing.assert_array_equal(
            lins.names, [f"Lineage {i}" for i in range(n_lins)]
        )
        np.testing.assert_array_equal(lins.names, adata_new.uns[names_key])

    @test_fwd()
    def test_non_unique_names(self, adata: AnnData, path: Path, lin_key: str, _: int):
        names_key = _lin_names(lin_key)
        adata.uns[names_key][0] = adata.uns[names_key][1]

        sc.write(path, adata)
        with pytest.raises(ValueError):
            _ = cr.read(path)

    @test_fwd()
    def test_wrong_colors_length(
        self, adata: AnnData, path: Path, lin_key: str, n_lins: int
    ):
        colors_key = _colors(lin_key)
        adata.uns[colors_key] = list(adata.uns[colors_key])
        adata.uns[colors_key] += [adata.uns[colors_key][0]]

        sc.write(path, adata)
        adata_new = cr.read(path)
        lins = adata_new.obsm[lin_key]

        assert isinstance(lins, Lineage)
        np.testing.assert_array_equal(lins.colors, _create_categorical_colors(n_lins))
        np.testing.assert_array_equal(lins.colors, adata_new.uns[colors_key])

    @test_fwd()
    def test_colors_not_colorlike(
        self, adata: AnnData, path: Path, lin_key: str, n_lins: int
    ):
        colors_key = _colors(lin_key)
        adata.uns[colors_key][0] = "foo"

        sc.write(path, adata)
        adata_new = cr.read(path)
        lins = adata_new.obsm[lin_key]

        assert isinstance(lins, Lineage)
        np.testing.assert_array_equal(lins.colors, _create_categorical_colors(n_lins))
        np.testing.assert_array_equal(lins.colors, adata_new.uns[colors_key])

    @test_fwd()
    def test_normal_run(self, adata: AnnData, path: Path, lin_key: str, n_lins: int):
        colors = _create_categorical_colors(10)[-n_lins:]
        names = [f"foo {i}" for i in range(n_lins)]

        adata.uns[_colors(lin_key)] = colors
        adata.uns[_lin_names(lin_key)] = names

        sc.write(path, adata)
        adata_new = cr.read(path)
        lins_new = adata_new.obsm[lin_key]

        np.testing.assert_array_equal(lins_new.colors, colors)
        np.testing.assert_array_equal(lins_new.names, names)
