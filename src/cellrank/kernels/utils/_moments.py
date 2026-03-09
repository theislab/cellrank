"""k-NN moment computation formerly imported from scVelo."""

import warnings

import numpy as np
import scipy.sparse as sp

__all__: list[str] = []


def _knn_moments(
    connectivities: sp.spmatrix,
    data: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    r"""Compute k-NN smoothed mean and variance.

    Parameters
    ----------
    connectivities
        Row-normalized k-NN connectivities matrix of shape ``(n_cells, n_cells)``.
    data
        Dense data matrix of shape ``(n_cells, n_features)``.

    Returns
    -------
    Tuple of ``(mean, variance)``, each of shape ``(n_cells, n_features)``.
    The mean is the k-NN weighted average, the variance is the centered
    second moment: :math:`\\text{Var} = W X^2 - (W X)^2`.
    """
    mean = connectivities @ data
    var = connectivities @ (data**2) - mean**2
    return mean, var


def _row_normalize_connectivities(conn: sp.spmatrix) -> sp.csr_matrix:
    """Build a row-normalized binary connectivities matrix.

    The transformation mirrors scVelo's ``get_connectivities``:
    binarize the graph, set the diagonal to 1 (self-loop), then
    row-normalize so each row sums to 1.

    Parameters
    ----------
    conn
        Sparse connectivities matrix of shape ``(n_cells, n_cells)``.

    Returns
    -------
    Row-normalized sparse matrix of shape ``(n_cells, n_cells)``.
    """
    C = conn.copy()

    # binarize
    C = (C > 0).astype(np.float32)

    # add self-loops
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        C.setdiag(1)

    # row-normalize
    row_sums = np.asarray(C.sum(axis=1)).ravel()
    row_sums[row_sums == 0] = 1
    C = sp.diags(1.0 / row_sums) @ C

    return C.tocsr().astype(np.float32)
