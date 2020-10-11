# -*- coding: utf-8 -*-
import os
from sys import version_info
from typing import Tuple, Union, Optional
from pathlib import Path

import pytest
from PIL import Image

import scanpy as sc
import scvelo as scv
from anndata import AnnData

import numpy as np
import pandas as pd
from sklearn.svm import SVR
from scipy.sparse import spdiags, issparse, csr_matrix

import cellrank as cr
from cellrank.tl.kernels import VelocityKernel, ConnectivityKernel


def _jax_not_installed() -> bool:
    try:
        import jax
        import jaxlib

        return False
    except ImportError:
        return True


def bias_knn(conn, pseudotime, n_neighbors, k=3):
    k_thresh = np.min([int(np.floor(n_neighbors / k)) - 1, 30])
    conn_biased = conn.copy()

    # check whether the original graph was connected
    assert _is_connected(conn), "The underlying KNN graph is disconnected."

    for i in range(conn.shape[0]):
        # get indices, values and current pseudo t
        row_data = conn[i, :].data
        row_ixs = conn[i, :].indices
        current_t = pseudotime[i]

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
    vk._mat_scaler = adata.uns.get(
        velocity_variances, np.random.normal(size=(adata.n_obs, adata.n_obs))
    )

    ck = ConnectivityKernel(adata)
    ck._mat_scaler = adata.uns.get(
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
                f"Invalid kind of conversion `{kind!r}`. Valid options are `'actual_to_expected'`, `'expected_to_actual'`."
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
    if not pickled:
        assert actual.adata is expected.adata
    else:
        assert actual.adata is not expected.adata

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
) -> None:
    assert actual is not expected
    assert actual.adata is not expected.adata
    assert actual.kernel is not expected.kernel
    if copy or version_info[:2] > (3, 6):
        # pickling of Enums doesn't work in Python3.6
        assert isinstance(actual.kernel, type(expected.kernel)), (
            type(actual.kernel),
            type(expected.kernel),
        )
    else:
        assert isinstance(actual.kernel, cr.tl.kernels.PrecomputedKernel)

    assert actual.adata.shape == expected.adata.shape
    assert actual.adata is actual.kernel.adata
    assert actual.kernel.backward == expected.kernel.backward

    np.testing.assert_array_equal(
        actual.transition_matrix.A, expected.transition_matrix.A
    )

    assert expected.__dict__.keys() == actual.__dict__.keys()

    for attr in expected.__dict__.keys():
        val2, val1 = getattr(actual, attr), getattr(expected, attr)
        if isinstance(val1, cr.tl.Lineage):
            assert val2 is not val1, attr
            assert_array_nan_equal(val2.X, val1.X)
        elif isinstance(val1, (np.ndarray, pd.Series, pd.DataFrame)):
            assert val2 is not val1, attr
            try:
                # can be array of strings, can't get NaN
                assert_array_nan_equal(val2, val1)
            except:
                np.testing.assert_array_equal(val2, val1)
        elif isinstance(val1, dict):
            assert val2.keys() == val1.keys()
            for v2, v1 in zip(val2.values(), val1.values()):
                if isinstance(v2, np.ndarray):
                    assert v2 is not v1, attr
                    np.testing.assert_array_equal(v2, v1)
                else:
                    assert v2 == v1, (v2, v1, attr)
        elif attr not in ("_kernel", "_gpcca"):
            assert val2 == val1, (val2, val1, attr)
        elif copy or version_info[:2] > (3, 6):
            # we can compare the kernel types, but for 3.6, it's saved as Precomputed
            assert isinstance(val2, type(val1)), (val2, val1, attr)


def random_transition_matrix(n: int) -> np.ndarray:
    """
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

    Call this function to regenerate the above objects.

    Params
    ------
    n_obs
        Number of cells.

    Returns
    -------
    :class:`anndata.AnnData`
        The created adata object.
    """

    np.random.seed(42)
    adata = scv.datasets.toy_data(n_obs=n_obs)
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=1000)
    adata.raw = adata[:, 42 : 42 + 50].copy()
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    scv.tl.recover_dynamics(adata)
    scv.tl.velocity(adata, mode="dynamical")
    scv.tl.velocity_graph(adata, mode_neighbors="connectivities")
    scv.tl.latent_time(adata)

    adata.uns["iroot"] = 0
    sc.tl.dpt(adata)

    adata.uns["connectivity_variances"] = np.ones((n_obs, n_obs), dtype=np.float64)
    adata.uns["velocity_variances"] = np.ones((n_obs, n_obs), dtype=np.float64)

    sc.write(f"tests/_ground_truth_adatas/adata_{n_obs}.h5ad", adata)

    return adata


jax_not_installed_skip = pytest.mark.skipif(
    _jax_not_installed(), reason="JAX is not installed."
)


if __name__ == "__main__":
    for size in [50, 100, 200]:
        _ = _create_dummy_adata(size)


def _import_rpy2_mgcv() -> bool:
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


gamr_skip = pytest.mark.skipif(
    _import_rpy2_mgcv(), reason="Cannot import `rpy2` or R's `mgcv` package."
)
