from typing import Optional

import pytest

import scanpy as sc
import cellrank.external as cre
from anndata import AnnData
from cellrank.tl.kernels import ConnectivityKernel
from cellrank.external.kernels._utils import MarkerGenes
from cellrank.external.kernels._wot_kernel import LastTimePoint

import numpy as np
import pandas as pd
from scipy.sparse import spmatrix
from pandas.core.dtypes.common import is_categorical_dtype


class TestOTKernel:
    def test_no_connectivities(self, adata_large: AnnData):
        del adata_large.obsp["connectivities"]
        terminal_states = np.full((adata_large.n_obs,), fill_value=np.nan, dtype=object)
        ixs = np.where(adata_large.obs["clusters"] == "Granule immature")[0]
        terminal_states[ixs] = "GI"

        ok = cre.kernels.StationaryOTKernel(
            adata_large,
            terminal_states=pd.Series(terminal_states).astype("category"),
            g=np.ones((adata_large.n_obs,), dtype=np.float64),
        )
        ok = ok.compute_transition_matrix(1, 0.001)

        assert ok._conn is None
        np.testing.assert_allclose(ok.transition_matrix.sum(1), 1.0)

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
        assert isinstance(ok._transition_matrix, (np.ndarray, spmatrix))
        np.testing.assert_allclose(ok.transition_matrix.sum(1), 1.0)
        assert isinstance(ok.params, dict)


@pytest.mark.skip("wot on PyPI doesn't support passing cost matrices")
class TestWOTKernel:
    def test_no_connectivities(self, adata_large: AnnData):
        del adata_large.obsp["connectivities"]
        ok = cre.kernels.WOTKernel(
            adata_large, time_key="age(days)"
        ).compute_transition_matrix()

        assert ok._conn is None
        np.testing.assert_allclose(ok.transition_matrix.sum(1), 1.0)

    def test_inversion_updates_adata(self, adata_large: AnnData):
        key = "age(days)"
        ok = cre.kernels.WOTKernel(adata_large, time_key=key)
        assert is_categorical_dtype(adata_large.obs[key])
        assert adata_large.obs[key].cat.ordered
        np.testing.assert_array_equal(ok.experimental_time, adata_large.obs[key])
        orig_time = ok.experimental_time

        ok = ~ok

        inverted_time = ok.experimental_time
        assert is_categorical_dtype(adata_large.obs[key])
        assert adata_large.obs[key].cat.ordered
        np.testing.assert_array_equal(ok.experimental_time, adata_large.obs[key])
        np.testing.assert_array_equal(
            orig_time.cat.categories, inverted_time.cat.categories
        )
        np.testing.assert_array_equal(orig_time.index, inverted_time.index)

        with pytest.raises(AssertionError):
            np.testing.assert_array_equal(orig_time, inverted_time)

    @pytest.mark.parametrize("cmat", [None, "Ms", "X_pca", "good_shape", "bad_shape"])
    def test_cost_matrices(self, adata_large: AnnData, cmat: str):
        ok = cre.kernels.WOTKernel(adata_large, time_key="age(days)")
        if isinstance(cmat, str) and "shape" in cmat:
            cost_matrices = {
                (12.0, 35.0): np.random.normal(size=(73 + ("bad" in cmat), 127))
            }
        else:
            cost_matrices = cmat

        if cmat == "bad_shape":
            with pytest.raises(ValueError, match=r"Expected cost matrix for time pair"):
                ok.compute_transition_matrix(cost_matrices=cost_matrices)
        else:
            ok = ok.compute_transition_matrix(cost_matrices=cost_matrices)
            param = ok.params["cost_matrices"]
            if cmat == "Ms":
                assert param == "layer:Ms"
            elif cmat == "X_pca":
                assert param == "obsm:X_pca"
            elif cmat == "good_shape":
                # careful, param is `nstr`, which does not equal anything
                assert str(param) == "precomputed"
            else:
                assert param == "default"

    @pytest.mark.parametrize("n_iters", [3, 5])
    def test_growth_rates(self, adata_large: AnnData, n_iters: int):
        ok = cre.kernels.WOTKernel(adata_large, time_key="age(days)")
        ok = ok.compute_transition_matrix(growth_iters=n_iters)

        assert isinstance(ok.growth_rates, pd.DataFrame)
        np.testing.assert_array_equal(adata_large.obs_names, ok.growth_rates.index)
        np.testing.assert_array_equal(
            ok.growth_rates.columns, [f"g{i}" for i in range(n_iters + 1)]
        )
        np.testing.assert_array_equal(
            adata_large.obs["estimated_growth_rates"], ok.growth_rates[f"g{n_iters}"]
        )
        assert ok.params["growth_iters"] == n_iters

    @pytest.mark.parametrize("key_added", [None, "gr"])
    def test_birth_death_process(self, adata_large: AnnData, key_added: Optional[str]):
        np.random.seed(42)
        adata_large.obs["foo"] = np.random.normal(size=(adata_large.n_obs,))
        adata_large.obs["bar"] = np.random.normal(size=(adata_large.n_obs,))

        ok = cre.kernels.WOTKernel(adata_large, time_key="age(days)")
        gr = ok.compute_initial_growth_rates("foo", "bar", key_added=key_added)

        if key_added is None:
            assert isinstance(gr, pd.Series)
            np.testing.assert_array_equal(gr.index, adata_large.obs_names)
        else:
            assert gr is None
            assert "gr" in adata_large.obs

    @pytest.mark.parametrize("ltp", list(LastTimePoint))
    def test_last_time_point(self, adata_large: AnnData, ltp: LastTimePoint):
        key = "age(days)"
        ok = cre.kernels.WOTKernel(adata_large, time_key=key).compute_transition_matrix(
            last_time_point=ltp.s, conn_kwargs={"n_neighbors": 11}
        )
        ixs = np.where(adata_large.obs[key] == 35.0)[0]

        T = ok.transition_matrix[ixs, :][:, ixs].A
        if ltp == LastTimePoint.UNIFORM:
            np.testing.assert_allclose(T, np.ones_like(T) / float(len(ixs)))
        elif ltp == LastTimePoint.DIAGONAL:
            np.testing.assert_allclose(T, np.eye(len(ixs)))
        elif ltp == LastTimePoint.CONNECTIVITIES:
            adata_subset = adata_large[adata_large.obs[key] == 35.0]
            sc.pp.neighbors(adata_subset, n_neighbors=11)
            T_actual = (
                ConnectivityKernel(adata_subset)
                .compute_transition_matrix()
                .transition_matrix.A
            )
            np.testing.assert_allclose(T, T_actual)

    @pytest.mark.parametrize("organism", ["human", "mouse"])
    def test_compute_scores_default(self, adata_large: AnnData, organism: str):
        pk, ak = "p_score", "a_score"
        if organism == "human":
            adata_large.var_names = adata_large.var_names.str.upper()

        ok = cre.kernels.WOTKernel(adata_large, time_key="age(days)")
        assert pk not in ok.adata.obs
        assert ak not in ok.adata.obs

        ok.compute_initial_growth_rates(
            organism=organism, proliferation_key=pk, apoptosis_key=ak, use_raw=False
        )

        assert pk in ok.adata.obs
        assert ak in ok.adata.obs

    def test_normal_run(self, adata_large: AnnData):
        ok = cre.kernels.WOTKernel(adata_large, time_key="age(days)")
        ok = ok.compute_transition_matrix()

        assert isinstance(ok, cre.kernels.WOTKernel)
        assert isinstance(ok._transition_matrix, (np.ndarray, spmatrix))
        np.testing.assert_allclose(ok.transition_matrix.sum(1), 1.0)
        assert isinstance(ok.params, dict)
        assert isinstance(ok.growth_rates, pd.DataFrame)
        assert isinstance(ok.transport_maps, dict)
        np.testing.assert_array_equal(adata_large.obs_names, ok.growth_rates.index)
        np.testing.assert_array_equal(ok.growth_rates.columns, ["g0", "g1"])
        assert isinstance(ok.transport_maps[12.0, 35.0], AnnData)

    def test_copy(self, adata_large: AnnData):
        ok = cre.kernels.WOTKernel(adata_large, time_key="age(days)")
        ok = ok.compute_transition_matrix()

        ok2 = ok.copy()

        assert isinstance(ok2, cre.kernels.WOTKernel)
        assert ok is not ok2
        np.testing.assert_array_equal(ok.transition_matrix.A, ok2.transition_matrix.A)


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
