# -*- coding: utf-8 -*-
from scipy.sparse.linalg import eigs

from cellrank.tools._utils import _eigengap
from cellrank.tools.kernels._kernel import KernelExpression
from typing import Optional, Dict

from abc import ABC, abstractmethod
import numpy as np

from anndata import AnnData
from scanpy import logging as logg
from scipy.sparse import issparse

from cellrank.tools._constants import Direction, RcKey, LinKey, Prefix


class BaseEstimator(ABC):
    def __init__(
        self,
        kernel: KernelExpression,
        adata: Optional[AnnData] = None,
        inplace: bool = True,
        read_from_adata: bool = True,
        g2m_key: Optional[str] = "G2M_score",
        s_key: Optional[str] = "S_score",
        key_added: Optional[str] = None,
    ):
        if adata is not None:
            self._adata = adata if inplace else adata.copy()
        else:
            logg.debug("DEBUG: Loading `adata` object from kernel.")
            self._adata = kernel.adata if inplace else kernel.adata.copy()

        if kernel.backward:
            self._direction: Direction = Direction.BACKWARD
            self._rc_key: str = str(RcKey.BACKWARD)
            self._lin_key: str = str(LinKey.BACKWARD)
            self._prefix: str = str(Prefix.BACKWARD)
        else:
            self._direction: Direction = Direction.FORWARD
            self._rc_key: str = str(RcKey.FORWARD)
            self._lin_key: str = str(LinKey.FORWARD)
            self._prefix: str = str(Prefix.FORWARD)

        # import transition matrix and parameters
        if kernel.transition_matrix is None:
            logg.debug("DEBUG: Computing transition matrix using default parameters.")
            kernel.compute_transition_matrix()
        kernel.write_to_adata(key_added=key_added)

        self._kernel = kernel
        self._T = kernel.transition_matrix
        self._is_sparse = issparse(self._T)
        self._n_states = self._T.shape[0]
        if self._n_states != self._adata.n_obs:
            raise ValueError(
                f"Expected `{self._n_states}` (based on transition matrix), "
                f"found `{self._adata.n_obs}` (based on `adata` object)."
            )

        # for copy
        self._g2m_key = g2m_key
        self._s_key = s_key
        self._key_added = key_added

        self._eig = None  # stores eigendecomposition
        self._dp = None  # stores differentiation potential

        if read_from_adata:
            logg.debug(
                "DEBUG: Reading `eig`, `approx_rcs` and `lin_probs` from `adata` object"
            )
            self._read_from_adata(g2m_key, s_key)

    def _compute_eig(
        self, k: int = 20, which: str = "LR", alpha: float = 1, only_evals: bool = False
    ) -> None:
        """
        Compute eigendecomposition of transition matrix.

        Uses a sparse implementation, if possible, and only computes the top k eigenvectors
        to speed up the computation. Computes both left and right eigenvectors.

        Params
        ------
        k
            Number of eigenvalues/vectors to compute.
        which
            Eigenvalues are in general complex. `'LR'` - largest real part, `'LM'` - largest magnitude.
        alpha
            Used to compute the `eigengap`. :paramref:`alpha` is the weight given
            to the deviation of an eigenvalue from one.

        Returns
        -------
        None
            Nothing, but updates the following fields: :paramref:`eigendecomposition`.
        """

        def get_top_k_evals():
            return D[np.flip(np.argsort(D.real))][
                :k
            ]  # TODO: @Marius - I think for non-sparse matrices [:k] is necessary

        def write_result(eig):
            # write to class and AnnData object
            if self._eig is not None:
                logg.debug("DEBUG: Overwriting `.eig`")
            else:
                logg.debug(
                    f"DEBUG: Adding `.eig` and `adata.uns['eig_{self._direction}']`"
                )

            self._eig = eig
            self._adata.uns[f"eig_{self._direction}"] = eig

        logg.info("Computing eigendecomposition of transition matrix")
        if self._is_sparse:
            logg.debug(f"DEBUG: Computing top `{k}` eigenvalues for sparse matrix")
            D, V_l = eigs(self._T.T, k=k, which=which)
            if only_evals:
                write_result(
                    {
                        "D": get_top_k_evals(),
                        "eigengap": _eigengap(D.real, alpha),
                        "params": {"which": which, "k": k, "alpha": alpha},
                    }
                )
                return
            _, V_r = eigs(self._T, k=k, which=which)
        else:
            logg.warning(
                "This transition matrix is not sparse, computing full eigendecomposition"
            )
            D, V_l = np.linalg.eig(self._T.T)
            if only_evals:
                write_result(
                    {
                        "D": get_top_k_evals(),
                        "eigengap": _eigengap(D.real, alpha),
                        "params": {"which": which, "k": k, "alpha": alpha},
                    }
                )
                return
            _, V_r = np.linalg.eig(self._T)

        # Sort the eigenvalues and eigenvectors and take the real part
        logg.debug("DEBUG: Sorting eigenvalues by their real part")
        p = np.flip(np.argsort(D.real))
        D, V_l, V_r = D[p], V_l[:, p], V_r[:, p]
        e_gap = _eigengap(D.real, alpha)

        write_result(
            {
                "D": D,
                "V_l": V_l,
                "V_r": V_r,
                "eigengap": e_gap,
                "params": {"which": which, "k": k, "alpha": alpha},
            }
        )

    @abstractmethod
    def compute_eig(self, k: int = 20, which: str = "LR", alpha: float = 1) -> None:
        pass

    def plot_spectrum(self):
        pass

    def plot_real_spectrum(self):
        pass

    @abstractmethod
    def _read_from_adata(
        self, g2m_key: Optional[str] = None, s_key: Optional[str] = None, **kwargs
    ) -> None:
        pass

    @abstractmethod
    def copy(self) -> "BaseEstimator":
        pass

    @property
    def eigendecomposition(self) -> Optional[Dict[str, np.ndarray]]:
        """
        A dictionary with following fields:

        - `'D'` eigenvalues of left eigenvectors
        - `'V_l'` left eigenvectors
        - `'V_r'` right eigenvectors
        """
        return self._eig

    @property
    def diff_potential(self) -> np.ndarray:
        """
        Differentiation potential for each lineage.
        """
        return self._dp

    @property
    def adata(self) -> AnnData:
        """
        The underlying annotated data object.
        """
        return self._adata

    @property
    def kernel(self) -> KernelExpression:
        """
        The underlying kernel expression.
        """
        return self._kernel

    def __copy__(self) -> "MarkovChain":
        return self.copy()

    def __len__(self) -> int:
        return self._n_states
