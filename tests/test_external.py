from typing import Optional

import pytest

import cellrank.external as cre
from anndata import AnnData
from cellrank.kernels import ConnectivityKernel
from cellrank.external.kernels._utils import MarkerGenes

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix


def _wot_not_installed() -> bool:
    try:
        import wot

        return False
    except ImportError:
        return True


def _statot_not_installed() -> bool:
    try:
        import statot

        return False
    except ImportError:
        return True


wot_not_installed_skip = pytest.mark.skipif(
    _wot_not_installed(), reason="WOT is not installed."
)
statot_not_installed_skip = pytest.mark.skipif(
    _statot_not_installed(), reason="statOT is not installed."
)


@statot_not_installed_skip
class TestOTKernel:
    def test_method_not_implemented(self, adata_large: AnnData):
        terminal_states = np.full((adata_large.n_obs,), fill_value=np.nan, dtype=object)
        ixs = np.where(adata_large.obs["clusters"] == "Granule immature")[0]
        terminal_states[ixs] = "GI"

        ok = cre.kernels.StationaryOTKernel(
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
            cre.kernels.StationaryOTKernel(
                adata_large,
                g=np.ones((adata_large.n_obs,), dtype=np.float64),
            )

    def test_normal_run(self, adata_large: AnnData):
        terminal_states = np.full((adata_large.n_obs,), fill_value=np.nan, dtype=object)
        ixs = np.where(adata_large.obs["clusters"] == "Granule immature")[0]
        terminal_states[ixs] = "GI"

        ok = cre.kernels.StationaryOTKernel(
            adata_large,
            terminal_states=pd.Series(terminal_states).astype("category"),
            g=np.ones((adata_large.n_obs,), dtype=np.float64),
        )
        ok = ok.compute_transition_matrix(1, 0.001)

        assert isinstance(ok, cre.kernels.StationaryOTKernel)
        assert isinstance(ok._transition_matrix, csr_matrix)
        np.testing.assert_allclose(ok.transition_matrix.sum(1), 1.0)
        assert isinstance(ok.params, dict)

    @pytest.mark.parametrize("connectivity_kernel", (None, ConnectivityKernel))
    def test_compute_projection(
        self, adata_large: AnnData, connectivity_kernel: Optional[ConnectivityKernel]
    ):
        terminal_states = np.full((adata_large.n_obs,), fill_value=np.nan, dtype=object)
        ixs = np.where(adata_large.obs["clusters"] == "Granule immature")[0]
        terminal_states[ixs] = "GI"

        ok = cre.kernels.StationaryOTKernel(
            adata_large,
            terminal_states=pd.Series(terminal_states).astype("category"),
            g=np.ones((adata_large.n_obs,), dtype=np.float64),
        )
        ok = ok.compute_transition_matrix(1, 0.001)

        if connectivity_kernel is not None:
            ck = connectivity_kernel(adata_large).compute_transition_matrix()
            combined_kernel = 0.9 * ok + 0.1 * ck
            combined_kernel.compute_transition_matrix()
            combined_kernel.plot_projection()
        else:
            with pytest.raises(RuntimeError, match=r"Unable to find connectivities"):
                ok.plot_projection()


class TestGetMarkers:
    @pytest.mark.parametrize("kind", ["proliferation", "apoptosis"])
    @pytest.mark.parametrize("organism", ["human", "mouse", "foo"])
    def test_get_markers(self, organism: str, kind: str):
        if organism == "foo":
            with pytest.raises(NotImplementedError, match=r""):
                getattr(MarkerGenes, f"{kind}_markers")(organism)
        else:
            markers = getattr(MarkerGenes, f"{kind}_markers")(organism)
            assert isinstance(markers, tuple)
            assert np.all([isinstance(marker, str) for marker in markers])
            if organism == "human":
                np.testing.assert_array_equal(
                    markers, [marker.upper() for marker in markers]
                )
