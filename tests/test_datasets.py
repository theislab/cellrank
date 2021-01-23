import sys

import pytest

from anndata import AnnData

import cellrank as cr


@pytest.mark.skipif(
    sys.version_info[:2] != (3, 8) or sys.platform != "linux",
    reason="Reduce the number of downloads.",
)
class TestDataSet:
    def test_pancreas(self, tmpdir_factory):
        adata = cr.datasets.pancreas(
            tmpdir_factory.mktemp("pancreas").join("adata.h5ad")
        )

        assert isinstance(adata, AnnData)
        assert adata.shape == (2531, 27998)

    def test_pancreas_preprocessed(self, tmpdir_factory):
        adata = cr.datasets.pancreas_preprocessed(
            tmpdir_factory.mktemp("pancreas_preprocessed").join("adata.h5ad")
        )

        assert isinstance(adata, AnnData)
        assert adata.shape == (2531, 2000)

    def test_lung(self, tmpdir_factory):
        adata = cr.datasets.lung(tmpdir_factory.mktemp("lung").join("adata.h5ad"))

        assert isinstance(adata, AnnData)
        assert adata.shape == (24882, 24051)
