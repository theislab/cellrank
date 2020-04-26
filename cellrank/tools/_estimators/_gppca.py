# -*- coding: utf-8 -*-
from typing import Optional, List, Tuple, Dict, Union
from anndata import AnnData
from msmtools.analysis.dense.gpcca import GPCCA as _GPPCA
from scanpy import logging as logg

from cellrank.tools._estimators._base_estimator import BaseEstimator
from cellrank.tools.kernels._kernel import KernelExpression

import numpy as np
import scvelo as scv


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
        self._gpcca = None
        self._schur_vectors = None
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
            Used to compute the `eigengap`. gref:`alpha` is the weight given
            to the deviation of an eigenvalue from one.

        Returns
        -------
        None
            Nothing, but updates the following fields: gref:`eigendecomposition`.
        """
        self._compute_eig(k=k, which=which, alpha=alpha, only_evals=True)

    def _read_from_adata(
        self, g2m_key: Optional[str] = None, s_key: Optional[str] = None, **kwargs
    ) -> None:
        raise NotImplementedError()

    def plot_schur_embedding(
        self,
        use: Optional[Union[int, tuple, list]] = None,
        abs_value: bool = False,
        cluster_key: Optional[str] = None,
        **kwargs,
    ) -> None:
        """
        Plt Schur vectors in an embedding.

        Params
        ------
        use
            Which or how many Schur vectors to be plotted. If `None`, all will be chosen.
        abs_value
            Whether to take the absolute value before plotting.
        cluster_key
            Key from :paramref:`adata` `.obs` to plot cluster annotations.

        Returns
        -------
        None
            Nothing, just plots the Schur vectors.
        """

        if self.schur_vectors is None:
            raise RuntimeError("Compute Schur vectors as `.metastable_state()` first.")

        self._plot_vectors(
            self.schur_vectors,
            "schur",
            abs_value=abs_value,
            use=use,
            cluster_key=cluster_key,
            **kwargs,
        )

    def metastable_states(
        self,
        n_states: Union[int, Tuple[int, int], List[int], Dict[str, int]],
        initial_distribution: Optional[np.ndarray] = None,
        use_min_chi: bool = False,
        method: str = "krylov",
        which: str = "LM",
        cluster_key: Optional[str] = "louvain",
    ):
        """
        Params
        ------
        n_states:
        initial_distribution
        use_min_chi
        method
        which
            Eigenvalues are in general complex. `'LR'` - largest real part, `'LM'` - largest magnitude.
        cluster_key

        Returns
        -------
        None
            Nothings, but updates the following fields:

            -

        """
        if (
            use_min_chi
        ):  # TODO: @Marius - this is the cleanest option I could thought of
            if not isinstance(n_states, (dict, tuple, list)):
                raise TypeError(
                    f"Expected `n_states` to be either `dict`, `tuple` or a `list`, found `{type(n_states).__name__}`."
                )
            if len(n_states) != 2:
                raise ValueError(
                    f"Expected `n_states` to be of size `2`, found `{len(n_states)}`."
                )

            minn, maxx = (
                n_states["n_min"],
                n_states["n_max"] if isinstance(n_states, dict) else n_states,
            )
            logg.debug(f"DEBUG: Calculating min Chi in interval [{minn}, {maxx}]")
            n_states = np.arange(minn, maxx)[np.argmax(self._gpcca.minChi(minn, maxx))]

        start = logg.info("Computing metastable states")

        self._gpcca = _GPPCA(self._T, eta=initial_distribution, z=which, method=method)
        self._gpcca.optimize(m=n_states)
        self._schur_vectors = self._gpcca.schur_vectors

        logg.info("Adding `...`\n" "    Finish", start=start)

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

    @property
    def schur_vectors(self):
        return self._schur_vectors
