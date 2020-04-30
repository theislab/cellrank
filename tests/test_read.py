# -*- coding: utf-8 -*-
from cellrank.tools._markov_chain import MarkovChain
from cellrank.tools._constants import _lin_names, _colors, LinKey

import cellrank as cr
import scanpy as sc

from typing import Tuple
from anndata import AnnData


class TestRead:
    def test_no_lineage(self, adata_mc_fwd: Tuple[AnnData, MarkovChain]):
        adata, _ = adata_mc_fwd
        lin_key = str(LinKey.FORWARD)

        del adata.obsm[lin_key]
        sc.write("no_lineage.h5ad", adata)
        adata_new = cr.read("no_lineage.h5ad")

        assert adata_new is not adata  # sanity check
        assert lin_key not in adata_new.obsm.keys()

    def test_no_names(self, adata_mc_fwd: Tuple[AnnData, MarkovChain]):
        adata, _ = adata_mc_fwd
        lin_key = str(LinKey.FORWARD)
        names_key = _lin_names(lin_key)

    def test_no_colors(self, adata_mc_fwd: Tuple[AnnData, MarkovChain]):
        adata, _ = adata_mc_fwd
        lin_key = str(LinKey.FORWARD)
        colors_key = _colors(lin_key)

    def test_wrong_names_length(self, adata_mc_fwd: Tuple[AnnData, MarkovChain]):
        adata, _ = adata_mc_fwd
        lin_key = str(LinKey.FORWARD)
        names_key = _lin_names(lin_key)

    def test_not_unique_names(self, adata_mc_fwd: Tuple[AnnData, MarkovChain]):
        adata, _ = adata_mc_fwd
        lin_key = str(LinKey.FORWARD)
        names_key = _lin_names(lin_key)

    def test_wrong_colors_length(self, adata_mc_fwd: Tuple[AnnData, MarkovChain]):
        adata, _ = adata_mc_fwd
        lin_key = str(LinKey.FORWARD)
        colors_key = _colors(lin_key)

    def test_colors_not_colorlike(self, adata_mc_fwd: Tuple[AnnData, MarkovChain]):
        adata, _ = adata_mc_fwd
        lin_key = str(LinKey.FORWARD)
        colors_key = _colors(lin_key)
