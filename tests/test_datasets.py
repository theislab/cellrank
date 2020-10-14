# -*- coding: utf-8 -*-
import sys

import pytest

from anndata import AnnData

import cellrank as cr


@pytest.mark.skipif(
    sys.version_info[:2] != (3, 8) or sys.platform != "linux",
    reason="Reduce the number of downloads.",
)
class TestDataSet:
    # don't test the lung since it's 140MiB
    def test_pancreas(self, tmpdir_factory):
        adata = cr.datasets.pancreas(
            tmpdir_factory.mktemp("pancreas").join("data.h5ad")
        )

        assert isinstance(adata, AnnData)
        assert adata.shape == (2531, 27998)

    def test_pancreas_preprocessed(self, tmpdir_factory):
        adata = cr.datasets.pancreas_preprocessed(
            tmpdir_factory.mktemp("pancreas_preprocessed").join("data.h5ad")
        )

        assert isinstance(adata, AnnData)
        assert adata.shape == (2531, 2000)
