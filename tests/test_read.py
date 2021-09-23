from typing import Tuple, Callable

from pathlib import Path

import scanpy as sc
import cellrank as cr
from anndata import AnnData
from cellrank._key import Key
from cellrank.tl._colors import _create_categorical_colors
from cellrank.tl._lineage import Lineage
from cellrank.tl.estimators import CFLARE

import numpy as np


def test_fwd():
    def wrapper(func: Callable) -> Callable:
        def decorator(self, adata_cflare_fwd: Tuple[AnnData, CFLARE], tmpdir):
            adata, _ = adata_cflare_fwd
            adata = adata.copy()

            dirname = func.__name__
            path = tmpdir.mkdir(dirname).join("tmp.h5ad")
            key = Key.obs.term_states(False)

            return func(
                self,
                adata,
                path,
                Key.obsm.abs_probs(False),
                len(adata.obs[key].cat.categories),
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
    def test_no_colors(self, adata: AnnData, path: Path, lin_key: str, n_lins: int):
        colors_key = Key.uns.colors(Key.obs.term_states(False))
        del adata.uns[colors_key]

        sc.write(path, adata)
        adata_new = cr.read(path)
        lins = adata_new.obsm[lin_key]

        assert isinstance(lins, Lineage)
        np.testing.assert_array_equal(lins.colors, _create_categorical_colors(n_lins))
        np.testing.assert_array_equal(lins.colors, adata_new.uns[colors_key])

    @test_fwd()
    def test_wrong_colors_length(
        self, adata: AnnData, path: Path, lin_key: str, n_lins: int
    ):
        colors_key = Key.uns.colors(Key.obs.term_states(False))
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
        colors_key = Key.uns.colors(Key.obs.term_states(False))
        adata.uns[colors_key][0] = "foo"

        sc.write(path, adata)
        adata_new = cr.read(path)
        lins = adata_new.obsm[lin_key]

        assert isinstance(lins, Lineage)
        np.testing.assert_array_equal(lins.colors, _create_categorical_colors(n_lins))
        np.testing.assert_array_equal(lins.colors, adata_new.uns[colors_key])

    @test_fwd()
    def test_normal_run(self, adata: AnnData, path: Path, lin_key: str, n_lins: int):
        key = Key.obs.term_states(False)
        colors = _create_categorical_colors(10)[-n_lins:]
        names = adata.obs[key].cat.categories
        adata.uns[Key.uns.colors(key)] = colors

        sc.write(path, adata)
        adata_new = cr.read(path)
        lins_new = adata_new.obsm[lin_key]

        np.testing.assert_array_equal(lins_new.names, names)
        np.testing.assert_array_equal(lins_new.colors, colors)
