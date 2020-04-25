# -*- coding: utf-8 -*-
from scipy.sparse.linalg import eigs

from pathlib import Path
from cellrank.tools._utils import _eigengap, save_fig
from cellrank.tools.kernels._kernel import KernelExpression
from typing import Optional, Dict, Union, Tuple
from abc import ABC, abstractmethod

import numpy as np
import matplotlib.pyplot as plt

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

    def plot_spectrum(
        self,
        dpi: int = 100,
        figsize: Optional[Tuple[float, float]] = (5, 5),
        legend_loc: Optional[str] = None,
        save: Optional[Union[str, Path]] = None,
    ) -> None:
        """
        Plot the top eigenvalues in complex plane.

        Params
        ------
        dpi
            Dots per inch.
        figsize
            Size of the figure.
        save
            Filename where to save the plots. If `None`, just shows the plot.

        Returns
        -------
        None
            Nothing, just plots the spectrum in complex plane.
        """

        if self._eig is None:
            logg.warning(
                "No eigendecomposition found, computing with default parameters"
            )
            self.compute_eig()
        D = self._eig["D"]
        params = self._eig["params"]

        # create fiture and axes
        fig, ax = plt.subplots(nrows=1, ncols=1, dpi=dpi, figsize=figsize)

        # get the original data ranges
        lam_x, lam_y = D.real, D.imag
        x_min, x_max = np.min(lam_x), np.max(lam_x)
        y_min, y_max = np.min(lam_y), np.max(lam_y)
        x_range, y_range = x_max - x_min, y_max - y_min
        final_range = np.max([x_range, y_range]) + 0.05

        # define a function to make the data limits rectangular
        adapt_range = lambda min_, max_, range_: (
            min_ + (max_ - min_) / 2 - range_ / 2,
            min_ + (max_ - min_) / 2 + range_ / 2,
        )
        x_min_, x_max_ = adapt_range(x_min, x_max, final_range)
        y_min_, y_max_ = adapt_range(y_min, y_max, final_range)

        # plot the data and the unit circle
        ax.scatter(D.real, D.imag, marker=".", label="Eigenvalue")
        t = np.linspace(0, 2 * np.pi, 500)
        x_circle, y_circle = np.sin(t), np.cos(t)
        ax.plot(x_circle, y_circle, "k-", label="Unit circle")

        # set labels, ranges and legend
        ax.set_xlabel("Im($\lambda$)")
        ax.set_xlim(x_min_, x_max_)

        ax.set_ylabel("Re($\lambda$)")
        ax.set_ylim(y_min_, y_max_)

        key = "real part" if params["which"] == "LR" else "magnitude"

        ax.set_title(f"Top {params['k']} eigenvalues according to their {key}")
        ax.legend(loc=legend_loc)

        if save is not None:
            save_fig(fig, save)

        fig.show()

    def plot_real_spectrum(
        self,
        dpi: int = 100,
        figsize: Optional[Tuple[float, float]] = None,
        legend_loc: Optional[str] = None,
        save: Optional[Union[str, Path]] = None,
    ) -> None:
        """
        Plot the real part of the top eigenvalues.

        Params
        ------
        dpi
            Dots per inch.
        figsize
            Size of the figure.
        save
            Filename where to save the plots. If `None`, just shows the plot.

        Returns
        -------
        None
            Nothing, just plots the spectrum.
        """

        if self._eig is None:
            logg.warning(
                "No eigendecomposition found, computing with default parameters"
            )
            self.compute_eig()

        # Obtain the eigendecomposition, create the color code
        D, params = self._eig["D"], self._eig["params"]
        D_real, D_imag = D.real, D.imag
        ixs = np.arange(len(D))
        mask = D_imag == 0

        # plot the top eigenvalues
        fig, ax = plt.subplots(nrows=1, ncols=1, dpi=dpi, figsize=figsize)
        ax.scatter(ixs[mask], D_real[mask], marker="o", label="Real eigenvalue")
        ax.scatter(ixs[~mask], D_real[~mask], marker="o", label="Complex eigenvalue")

        # add dashed line for the eigengap, ticks, labels, title and legend
        ax.axvline(self._eig["eigengap"], label="Eigengap", ls="--")

        ax.set_xlabel("index")
        ax.set_xticks(range(len(D)))

        ax.set_ylabel("Re($\lambda_i$)")
        key = "real part" if params["which"] == "LR" else "magnitude"
        ax.set_title(
            f"Real part of top {params['k']} eigenvalues according to their {key}"
        )

        ax.legend(loc=legend_loc)

        if save is not None:
            save_fig(fig, save)

        fig.show()

    @abstractmethod
    def _read_from_adata(
        self, g2m_key: Optional[str] = None, s_key: Optional[str] = None, **kwargs
    ) -> None:
        pass

    @abstractmethod
    def copy(self) -> "BaseEstimator":
        # TODO: add copy of eig etc.
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
