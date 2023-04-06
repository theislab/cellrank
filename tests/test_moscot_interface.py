import pytest

import scanpy as sc
import cellrank.external as cre
from anndata import AnnData
from cellrank.kernels import MoscotKernel, ConnectivityKernel
from cellrank.external.kernels._utils import MarkerGenes
from cellrank.kernels._transport_map_kernel import SelfTransitions

import numpy as np
import pandas as pd
from scipy.sparse import issparse, csr_matrix
from pandas.core.dtypes.common import is_categorical_dtype

import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import to_hex


def _moscot_not_installed() -> bool:
    try:
        import moscot

        return False
    except ImportError:
        return True


moscot_not_installed_skip = pytest.mark.skipif(
    _moscot_not_installed(), reason="moscot is not installed."
)


@moscot_not_installed_skip
class TestMoscotKernel:
    from moscot.problems import LineageProblem, TemporalProblem, SpatioTemporalProblem

    def test_init_from_adata(self, adata_moscot: AnnData):
        pass

    @pytest.mark.parametrize("backward", [True, False])
    def test_init_from_TemporalProblem(
        self, moscot_tp: TemporalProblem, backward: bool
    ):
        from moscot.problems import TemporalProblem

        mk = MoscotKernel.load(moscot_tp, backward=backward)
        assert mk.transport_maps is None
        assert mk.obs is None
        assert isinstance(mk.problem, TemporalProblem)

    @pytest.mark.parametrize("backward", [True, False])
    def test_init_from_SpatioTemporalProblem(
        self, moscot_stp: SpatioTemporalProblem, backward: bool
    ):
        from moscot.problems import SpatioTemporalProblem

        mk = MoscotKernel.load(moscot_stp, backward=backward)
        assert mk.transport_maps is None
        assert mk.obs is None
        assert isinstance(mk.problem, SpatioTemporalProblem)

    @pytest.mark.parametrize("backward", [True, False])
    def test_init_from_LineageProblem(self, moscot_lp: LineageProblem, backward: bool):
        from moscot.problems import LineageProblem

        mk = MoscotKernel.load(moscot_lp, backward=backward)
        assert mk.transport_maps is None
        assert mk.obs is None
        assert isinstance(mk.problem, LineageProblem)

    @pytest.mark.parametrize("backward", [True, False])
    def test_compute_transition_matrix(
        self, moscot_tp: TemporalProblem, backward: bool
    ):
        mk = MoscotKernel.load(moscot_tp, backward=backward)
        assert mk.transport_maps is None
        assert mk.obs is None

        mk = mk.compute_transition_matrix()
        assert isinstance(mk.transport_maps, dict)
        assert len(mk.transport_maps) == 2
        assert mk.transition_matrix is not None

    @pytest.mark.parametrize("backward", [True, False])
    def test_compute_transition_matrix_from_output(
        self, moscot_tp: TemporalProblem, backward: bool
    ):
        mk = MoscotKernel.load(moscot_tp, backward=backward)
        assert mk.transport_maps is None
        assert mk.obs is None

        tmaps = {p: moscot_tp[p].solution for p in moscot_tp.problems}
        mk = mk.compute_transition_matrix(tmaps=tmaps)
        assert isinstance(mk.transport_maps, dict)
        assert len(mk.transport_maps) == 2
        assert mk.transition_matrix is not None

    def test_inversion_updates_adata(self, moscot_tp: TemporalProblem):
        key = "day"
        mk = MoscotKernel.load(moscot_tp, time_key=key)
        assert is_categorical_dtype(moscot_tp.adata.obs[key])
        assert moscot_tp.adata.obs[key].cat.ordered
        np.testing.assert_array_equal(mk.experimental_time, moscot_tp.adata.obs[key])
        orig_time = mk.experimental_time

        mk = ~mk

        assert mk.transition_matrix is None
        assert mk.transport_maps is None

        inverted_time = mk.experimental_time
        assert is_categorical_dtype(moscot_tp.adata.obs[key])
        assert moscot_tp.adata.obs[key].cat.ordered
        np.testing.assert_array_equal(mk.experimental_time, moscot_tp.adata.obs[key])
        np.testing.assert_array_equal(
            orig_time.cat.categories, inverted_time.cat.categories
        )
        np.testing.assert_array_equal(orig_time.index, inverted_time.index)

        with pytest.raises(AssertionError):
            np.testing.assert_array_equal(orig_time, inverted_time)

    @pytest.mark.parametrize("cmap", ["inferno", "viridis"])
    def test_update_colors(self, moscot_tp: TemporalProblem, cmap: str):
        ckey = "day_colors"
        _ = MoscotKernel.load(moscot_tp, cmap=cmap)

        colors = moscot_tp.adata.uns[ckey]
        cmap = get_cmap(cmap)

        assert isinstance(colors, np.ndarray)
        assert colors.shape == (3,)
        np.testing.assert_array_equal(colors, [to_hex(cmap(0)), to_hex(cmap(cmap.N))])

    @pytest.mark.parametrize("st", list(SelfTransitions))
    def test_self_transitions(self, moscot_tp: TemporalProblem, st: SelfTransitions):
        et, lt = 0, 1
        key = "day"
        conn_kwargs = {"n_neighbors": 11}

        mk = MoscotKernel.load(moscot_tp).compute_transition_matrix(
            self_transitions=st,
            conn_weight=0.2,
            conn_kwargs=conn_kwargs,
            threshold=None,
        )

        ixs = np.where(moscot_tp.adata.obs[key] == lt)[0]
        T = mk.transition_matrix[ixs, :][:, ixs].A

        assert mk.params["self_transitions"] == str(st)
        if st == SelfTransitions.UNIFORM:
            np.testing.assert_allclose(T, np.ones_like(T) / float(len(ixs)))
        elif st == SelfTransitions.DIAGONAL:
            np.testing.assert_allclose(T, np.eye(len(ixs)))
        elif st in (SelfTransitions.CONNECTIVITIES, SelfTransitions.ALL):
            adata_subset = moscot_tp.adata[moscot_tp.adata.obs[key] == lt]
            sc.pp.neighbors(adata_subset, **conn_kwargs)
            T_actual = (
                ConnectivityKernel(adata_subset)
                .compute_transition_matrix()
                .transition_matrix.A
            )
            np.testing.assert_allclose(T, T_actual)
            if st == SelfTransitions.ALL:
                ixs = np.where(moscot_tp.adata.obs[key] == et)[0]
                T = mk.transition_matrix[ixs, :][:, ixs].A
                adata_subset = moscot_tp[moscot_tp.adata.obs[key] == et]
                sc.pp.neighbors(adata_subset, **conn_kwargs)
                T_actual = (
                    ConnectivityKernel(adata_subset)
                    .compute_transition_matrix()
                    .transition_matrix.A
                )
                np.testing.assert_allclose(T, 0.2 * T_actual)
        else:
            raise NotImplementedError(st)
