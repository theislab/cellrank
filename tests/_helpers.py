from typing import Tuple, Union, Optional

import os
import pytest
from PIL import Image
from pathlib import Path

import scanpy as sc
import scvelo as scv
import cellrank as cr
from anndata import AnnData
from cellrank.tl.kernels import VelocityKernel, PrecomputedKernel, ConnectivityKernel

import numpy as np
import pandas as pd
from sklearn.svm import SVR
from scipy.sparse import spdiags, issparse, csr_matrix
from pandas.testing import assert_frame_equal, assert_series_equal


def _jax_not_installed() -> bool:
    try:
        import jax
        import jaxlib

        return False
    except ImportError:
        return True


def _rpy2_mgcv_not_installed() -> bool:
    try:
        import rpy2
        from packaging import version
        from rpy2.robjects.packages import PackageNotInstalledError, importr

        try:
            from importlib_metadata import version as get_version
        except ImportError:
            # >=Python3.8
            from importlib.metadata import version as get_version

        try:
            assert version.parse(get_version(rpy2.__name__)) >= version.parse("3.3.0")
            _ = importr("mgcv")
            return False
        except (PackageNotInstalledError, AssertionError):
            pass

    except ImportError:
        pass

    return True


def bias_knn(
    conn: csr_matrix,
    pseudotime: np.ndarray,
    n_neighbors: int,
    k: int = 3,
    frac_to_keep: Optional[float] = None,
) -> csr_matrix:
    # frac_to_keep=None mimics original impl. (which mimics Palantir)
    k_thresh = max(0, min(int(np.floor(n_neighbors / k)) - 1, 30))
    conn_biased = conn.copy()

    # check whether the original graph was connected
    assert _is_connected(conn), "The underlying KNN graph is disconnected."

    for i in range(conn.shape[0]):
        # get indices, values and current pseudo t
        row_data = conn[i, :].data
        row_ixs = conn[i, :].indices
        current_t = pseudotime[i]

        if frac_to_keep is not None:
            k_thresh = max(0, min(30, int(np.floor(len(row_data) * frac_to_keep))))

        # get the 'candidates' - ixs of nodes not in the k_thresh closest neighbors
        p = np.flip(np.argsort(row_data))
        sorted_ixs = row_ixs[p]
        cand_ixs = sorted_ixs[k_thresh:]

        # compare pseudotimes and set indices to zero
        cand_t = pseudotime[cand_ixs]
        rem_ixs = cand_ixs[cand_t < current_t]
        conn_biased[i, rem_ixs] = 0

    conn_biased.eliminate_zeros()

    # check whether the biased graph is still connected
    assert _is_connected(conn_biased), "The biased KNN graph has become disconnected."

    return conn_biased


def density_normalization(velo_graph, trans_graph):
    # function copied from scanpy
    q = np.asarray(trans_graph.sum(axis=0))
    if not issparse(trans_graph):
        Q = np.diag(1.0 / q)
    else:
        Q = spdiags(1.0 / q, 0, trans_graph.shape[0], trans_graph.shape[0])
    velo_graph = Q @ velo_graph @ Q

    return velo_graph


def _is_connected(c) -> bool:
    import networkx as nx
    from scipy.sparse import issparse

    G = nx.from_scipy_sparse_matrix(c) if issparse(c) else nx.from_numpy_array(c)

    return nx.is_connected(G)


def create_kernels(
    adata: AnnData,
    velocity_variances: Optional[str] = None,
    connectivity_variances: Optional[str] = None,
) -> Tuple[VelocityKernel, ConnectivityKernel]:
    vk = VelocityKernel(adata)
    vk._mat_scaler = adata.obsp.get(
        velocity_variances, np.random.normal(size=(adata.n_obs, adata.n_obs))
    )

    ck = ConnectivityKernel(adata)
    ck._mat_scaler = adata.obsp.get(
        connectivity_variances, np.random.normal(size=(adata.n_obs, adata.n_obs))
    )

    vk._transition_matrix = csr_matrix(np.eye(adata.n_obs))
    ck._transition_matrix = np.eye(adata.n_obs, k=1) / 2 + np.eye(adata.n_obs) / 2
    ck._transition_matrix[-1, -1] = 1
    ck._transition_matrix = csr_matrix(ck._transition_matrix)

    np.testing.assert_allclose(
        np.sum(ck._transition_matrix.A, axis=1), 1
    )  # sanity check

    return vk, ck


# TODO: make it a fixture
def create_model(adata: AnnData) -> cr.ul.models.SKLearnModel:
    return cr.ul.models.SKLearnModel(adata, SVR(kernel="rbf"))


# TODO: make it a fixture
def create_failed_model(adata: AnnData) -> cr.ul.models.FailedModel:
    return cr.ul.models.FailedModel(create_model(adata), exc="foobar")


def resize_images_to_same_sizes(
    expected_image_path: Union[str, Path],
    actual_image_path: Union[str, Path],
    kind: str = "actual_to_expected",
) -> None:
    if not os.path.isfile(actual_image_path):
        raise OSError(f"Actual image path `{actual_image_path!r}` does not exist.")
    if not os.path.isfile(expected_image_path):
        raise OSError(f"Expected image path `{expected_image_path!r}` does not exist.")
    expected_image = Image.open(expected_image_path)
    actual_image = Image.open(actual_image_path)
    if expected_image.size != actual_image.size:
        if kind == "actual_to_expected":
            actual_image.resize(expected_image.size).save(actual_image_path)
        elif kind == "expected_to_actual":
            expected_image.resize(actual_image.size).save(expected_image)
        else:
            raise ValueError(
                f"Invalid kind of conversion `{kind!r}`."
                f"Valid options are `'actual_to_expected'`, `'expected_to_actual'`."
            )


def assert_array_nan_equal(
    actual: Union[np.ndarray, pd.Series], expected: Union[np.ndarray, pd.Series]
) -> None:
    """
    Test is 2 arrays or :class:`pandas.Series` are equal.

    Params
    ------
    actual
        The actual data.
    expected
        The expected result.

    Returns
    -------
    None
        Nothing, but raises an exception if arrays are not equal, including the locations of NaN values.
    """

    mask1 = ~(pd.isnull(actual) if isinstance(actual, pd.Series) else np.isnan(actual))
    mask2 = ~(
        pd.isnull(expected) if isinstance(expected, pd.Series) else np.isnan(expected)
    )
    np.testing.assert_array_equal(np.where(mask1), np.where(mask2))
    np.testing.assert_array_equal(actual[mask1], expected[mask2])


def assert_models_equal(
    expected: cr.ul.models.BaseModel,
    actual: cr.ul.models.BaseModel,
    pickled: bool = False,
    deepcopy: bool = True,
) -> None:
    assert actual is not expected
    if pickled:
        assert actual.adata is not expected.adata
    else:
        assert actual.adata is expected.adata

    assert actual.shape == expected.shape
    assert actual.adata.shape == expected.adata.shape

    assert expected.__dict__.keys() == actual.__dict__.keys()

    for attr in expected.__dict__.keys():
        val2, val1 = getattr(actual, attr), getattr(expected, attr)
        if attr == "_prepared":
            # we expect the expected model to be prepared only if deepcopied
            if deepcopy:
                assert val2 == val1
            else:
                assert not val2
                assert val1
        elif isinstance(val1, cr.tl.Lineage):
            if deepcopy or pickled:
                assert val2 is not val1
                assert_array_nan_equal(val2.X, val1.X)
            else:
                assert val2 is val1, (val2, val1, attr)
        elif isinstance(val1, (np.ndarray, pd.Series, pd.DataFrame)):
            if deepcopy or pickled:
                try:
                    assert val2 is not val1, attr
                    # can be array of strings, can't get NaN
                    assert_array_nan_equal(val2, val1)
                except:
                    np.testing.assert_array_equal(val2, val1)
            # e.g. for GAMR, we point to the offset and design matrix
            # however, the `x`, and so pointers are not modified
            elif val2 is not None:
                assert val2 is val1, attr
        # we don't expect any dictionaries as in estimators
        elif attr == "_model":
            assert val2 is not val1  # model is always deepcopied
        elif not isinstance(val2, AnnData) and not callable(val2):
            # callable because SKLearnModel has default conf int function
            assert val2 == val1, (val2, val1, attr)
        else:
            assert isinstance(val2, type(val1)), (val2, val1, attr)


def assert_estimators_equal(
    expected: cr.tl.estimators.BaseEstimator,
    actual: cr.tl.estimators.BaseEstimator,
    copy: bool = False,
    deep: bool = False,
    from_adata: bool = False,
) -> None:
    def check_arrays(x, y):
        if isinstance(x, (np.ndarray, list, tuple)):
            try:
                np.testing.assert_array_compare(np.array_equal, x, y, equal_nan=True)
            except:
                np.testing.assert_array_compare(np.allclose, x, y, equal_nan=True)
        elif isinstance(x, pd.Series):
            assert_series_equal(x, y, check_names=False)
        elif isinstance(x, pd.DataFrame):
            assert_frame_equal(x, y, check_dtype=False)

    assert actual is not expected
    if copy:
        if deep:
            assert actual.adata is not expected.adata
        else:
            assert actual.adata is expected.adata
    else:
        assert actual.adata is not expected.adata
    assert actual.kernel is not expected.kernel
    if from_adata:
        assert isinstance(actual.kernel, PrecomputedKernel)
    else:
        assert isinstance(actual.kernel, type(expected.kernel))

    assert actual.adata.shape == expected.adata.shape
    assert actual.adata is actual.kernel.adata
    assert actual.kernel.backward == expected.kernel.backward

    np.testing.assert_array_equal(
        actual.transition_matrix.A, expected.transition_matrix.A
    )

    k1 = sorted(expected.__dict__.keys())
    k2 = sorted(actual.__dict__.keys())
    np.testing.assert_array_equal(k1, k2)

    for attr in expected.__dict__.keys():
        if attr == "_invalid_n_states" and from_adata:
            continue
        actual_val, expected_val = getattr(actual, attr), getattr(expected, attr)
        if isinstance(actual_val, cr.tl.Lineage):
            assert actual_val is not expected_val, attr
            assert_array_nan_equal(actual_val.X, expected_val.X)
        elif isinstance(actual_val, (np.ndarray, pd.Series, pd.DataFrame, list, tuple)):
            assert actual_val is not expected_val, attr
            check_arrays(actual_val, expected_val)
        elif isinstance(actual_val, dict):
            if from_adata:
                # _params can sometimes contain extra empty dict if initialized from `adata`
                for k in set(actual_val.keys()) | set(expected_val.keys()):
                    v2, v1 = actual_val.get(k, {}), expected_val.get(k, {})
                    if isinstance(v1, (np.ndarray, pd.Series, pd.DataFrame)):
                        check_arrays(v2, v1)
                    else:
                        assert v2 == v1, (v2, v1, attr, k)
            else:
                np.testing.assert_array_equal(
                    sorted(actual_val.keys()), sorted(expected_val.keys())
                )
                for k in sorted(actual_val.keys()):
                    v2, v1 = actual_val[k], expected_val[k]
                    if isinstance(v1, (np.ndarray, pd.Series, pd.DataFrame)):
                        check_arrays(v2, v1)
                    else:
                        assert v2 == v1, (v2, v1, attr, k)
        elif attr not in ("_kernel", "_gpcca", "_adata", "_shadow_adata"):
            assert actual_val == expected_val, (actual_val, expected_val, attr)
        else:
            try:
                assert isinstance(actual_val, type(expected_val)), (
                    actual_val,
                    expected_val,
                    attr,
                )
            except AssertionError:
                # objects initialized from `adata` don't have `_gpcca`
                if attr != "_gpcca" and not from_adata:
                    raise


def random_transition_matrix(n: int) -> np.ndarray:
    """
    Create a random transition matrix.

    Parameters
    ----------
    n
        Number of states.

    Returns
    -------
    :class:`numpy.ndarray`
        Row-normalized transition matrix.
    """

    x = np.abs(np.random.normal(size=(n, n)))
    rsum = x.sum(axis=1)
    return x / rsum[:, np.newaxis]


def _create_dummy_adata(n_obs: int) -> AnnData:
    """
    Create a testing :class:`anndata.AnnData` object.

    Call this function to regenerate the ground truth objects.

    Parameters
    ----------
    n_obs
        Number of cells.

    Returns
    -------
    :class:`anndata.AnnData`
        The created adata object.
    """

    np.random.seed(42)
    adata = scv.datasets.toy_data(n_obs=n_obs)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=1000)
    adata.var["symbol"] = adata.var_names.str.cat(["gs"] * adata.n_vars, sep=":")

    raw = adata[:, 42 : 42 + 50].copy()
    raw.var["symbol"] = raw.var_names.str.cat(["gs:raw"] * raw.n_vars, sep=":")
    adata.raw = raw

    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    scv.tl.recover_dynamics(adata)
    scv.tl.velocity(adata, mode="dynamical")
    scv.tl.velocity_graph(adata, mode_neighbors="connectivities")
    scv.tl.latent_time(adata)

    adata.uns["iroot"] = 0
    sc.tl.dpt(adata)

    if "velocity_graph" in adata.uns:
        adata.obsp["velocity_graph"] = adata.uns.pop("velocity_graph")
        adata.obsp["velocity_graph_neg"] = adata.uns.pop("velocity_graph_neg")
    adata.obsp["connectivity_variances"] = np.ones((n_obs, n_obs), dtype=np.float64)
    adata.obsp["velocity_variances"] = np.ones((n_obs, n_obs), dtype=np.float64)

    sc.write(f"tests/_ground_truth_adatas/adata_{n_obs}.h5ad", adata)

    return adata


jax_not_installed_skip = pytest.mark.skipif(
    _jax_not_installed(), reason="JAX is not installed."
)
gamr_skip = pytest.mark.skipif(
    _rpy2_mgcv_not_installed(), reason="Cannot import `rpy2` or R's `mgcv` package."
)

if __name__ == "__main__":
    for size in [50, 100, 200]:
        _ = _create_dummy_adata(size)
