import pytest

from anndata import AnnData

import numpy as np
from scipy.sparse import spmatrix

import cellrank.external as cre


class TestOTKernel:
    def test_method_not_implemented(self, adata_large: AnnData):
        source = np.zeros((adata_large.n_obs,), dtype=np.bool_)
        sink = np.zeros((adata_large.n_obs,), dtype=np.bool_)
        ok = cre.kernels.OTKernel(
            adata_large,
            source_idx=source,
            sink_idx=sink,
            g=np.ones((adata_large.n_obs,), dtype=np.float64),
        )

        with pytest.raises(
            NotImplementedError, match="Method `'unbal'` is not yet implemented."
        ):
            ok.compute_transition_matrix(1, 0.001, method="unbal")

    def test_some_both_source_sink(self, adata_large: AnnData):
        source = np.zeros((adata_large.n_obs,), dtype=np.bool_)
        source[0] = True
        sink = np.zeros((adata_large.n_obs,), dtype=np.bool_)
        sink[0] = True

        with pytest.raises(
            ValueError, match="Some cells are marked as both source and sink."
        ):
            cre.kernels.OTKernel(
                adata_large,
                source_idx=source,
                sink_idx=sink,
                g=np.ones((adata_large.n_obs,), dtype=np.float64),
            )

    def test_normal_run(self, adata_large: AnnData):
        source = (adata_large.obs["clusters"] == "Granule immature").values
        sink = (adata_large.obs["clusters"] == "Granule mature").values

        ok = cre.kernels.OTKernel(
            adata_large,
            source_idx=source,
            sink_idx=sink,
            g=np.ones((adata_large.n_obs,), dtype=np.float64),
        )
        ok = ok.compute_transition_matrix(1, 0.001)

        assert isinstance(ok, cre.kernels.OTKernel)
        assert isinstance(ok._transition_matrix, (np.ndarray, spmatrix))
        np.testing.assert_allclose(ok.transition_matrix.sum(1), 1.0)
        assert isinstance(ok.params, dict)
