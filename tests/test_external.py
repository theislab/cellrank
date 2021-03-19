import pytest

from anndata import AnnData

import numpy as np
import pandas as pd
from scipy.sparse import spmatrix

import cellrank.external as cre


class TestOTKernel:
    def test_method_not_implemented(self, adata_large: AnnData):
        terminal_states = np.full((adata_large.n_obs,), fill_value=np.nan, dtype=object)
        ixs = np.where(adata_large.obs["clusters"] == "Granule immature")[0]
        terminal_states[ixs] = "GI"

        ok = cre.kernels.OTKernel(
            adata_large,
            terminal_states=pd.Series(terminal_states).astype("category"),
            g=np.ones((adata_large.n_obs,), dtype=np.float64),
        )

        with pytest.raises(
            NotImplementedError, match="Method `'unbal'` is not yet implemented."
        ):
            ok.compute_transition_matrix(1, 0.001, method="unbal")

    def test_no_terminal_states(self, adata_large: AnnData):
        with pytest.raises(RuntimeError, match="Unable to initialize the kernel."):
            cre.kernels.OTKernel(
                adata_large,
                g=np.ones((adata_large.n_obs,), dtype=np.float64),
            )

    def test_normal_run(self, adata_large: AnnData):
        terminal_states = np.full((adata_large.n_obs,), fill_value=np.nan, dtype=object)
        ixs = np.where(adata_large.obs["clusters"] == "Granule immature")[0]
        terminal_states[ixs] = "GI"

        ok = cre.kernels.OTKernel(
            adata_large,
            terminal_states=pd.Series(terminal_states).astype("category"),
            g=np.ones((adata_large.n_obs,), dtype=np.float64),
        )
        ok = ok.compute_transition_matrix(1, 0.001)

        assert isinstance(ok, cre.kernels.OTKernel)
        assert isinstance(ok._transition_matrix, (np.ndarray, spmatrix))
        np.testing.assert_allclose(ok.transition_matrix.sum(1), 1.0)
        assert isinstance(ok.params, dict)
