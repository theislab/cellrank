# -*- coding: utf-8 -*-
from typing import Optional, List
from anndata import AnnData

from cellrank.tools._estimators._base_estimator import BaseEstimator
from cellrank.tools.kernels._kernel import KernelExpression


class GPCCA(BaseEstimator):
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
        super().__init__(
            kernel,
            adata,
            inplace=inplace,
            read_from_adata=read_from_adata,
            g2m_key=g2m_key,
            s_key=s_key,
            key_added=key_added,
        )
        # TODO: add extra info here

    def compute_eig(self, k: int = 20, which: str = "LR", alpha: float = 1) -> None:
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
        self._compute_eig(k=k, which=which, alpha=alpha, only_evals=True)

    def _read_from_adata(
        self, g2m_key: Optional[str] = None, s_key: Optional[str] = None, **kwargs
    ) -> None:
        raise NotImplementedError()

    def plot_schur_embedding(self, **kwargs):
        raise NotImplementedError()

    def metastable_states(
        self, m: Optional[int] = None, cluster_key: Optional[str] = "louvain", **kwargs
    ):
        raise NotImplementedError()

    def plot_metastable_states(
        self, n_cells: Optional[int] = None, same_plot: bool = True, **kwargs
    ):
        raise NotImplementedError()

    def plot_coarse_T(
        self, stationary_distribution: bool = False, clustermap: bool = False, **kwargs
    ):
        raise NotImplementedError()

    def set_main_states(self, keys: Optional[List[str]] = None, **kwargs):
        raise NotImplementedError()

    def copy(self) -> "GPCCA":
        raise NotImplementedError()
