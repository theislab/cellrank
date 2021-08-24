from typing import Any, Dict, Tuple, Union, Optional
from typing_extensions import Literal, Protocol

from pathlib import Path
from datetime import datetime

from anndata import AnnData
from cellrank import logging as logg
from cellrank.ul._docs import d
from cellrank.tl._utils import save_fig, _eigengap
from cellrank.tl._estimators.mixins._constants import Key
from cellrank.tl._estimators.mixins.decomposition._plot import VectorPlottable

import numpy as np
from scipy.sparse import issparse, csr_matrix
from scipy.sparse.linalg import eigs

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

EPS = np.finfo(np.float64).eps


class EigenProtocol(Protocol):
    @property
    def adata(self) -> AnnData:
        ...

    @property
    def backward(self) -> bool:
        ...

    @property
    def transition_matrix(self) -> Union[csr_matrix, np.ndarray]:
        ...

    def _write_eigendecomposition(self, decomp: Dict[str, Any], time: datetime) -> None:
        ...


class EigenMixin(VectorPlottable):
    """TODO."""

    def __init__(self):
        super().__init__()
        self._eigendecomposition = None

    @d.dedent
    def compute_eigendecomposition(
        self: EigenProtocol,
        k: int = 20,
        which: Literal["LR", "LM"] = "LR",
        alpha: float = 1.0,
        only_evals: bool = False,
        ncv: Optional[int] = None,
    ) -> None:
        """
        Compute eigendecomposition of transition matrix.

        Uses a sparse implementation, if possible, and only computes the top :math:`k` eigenvectors
        to speed up the computation. Computes both left and right eigenvectors.

        Parameters
        ----------
        k
            Number of eigenvalues/vectors to compute.
        %(eigen)s
        only_evals
            Compute only eigenvalues.
        ncv
            Number of Lanczos vectors generated.

        Returns
        -------
        Nothing, but updates the following field:

            - :attr:`eigendecomposition` ``['D']`` - the eigenvalues.
            - :attr:`eigendecomposition` ``['eigengap']`` - the eigengap.
            - :attr:`eigendecomposition` ``['params']`` - parameters used for the computation.

        If ``only_evals = False``, also updates:

            - :attr:`eigendecomposition` ``['V_l']`` - left eigenvectors.
            - :attr:`eigendecomposition` ``['V_r']`` - right eigenvectors.
            - :attr:`eigendecomposition` ``['stationary_dist']`` - stationary distribution of :attr:`transition_matrix`.
        """

        def get_top_k_evals() -> np.ndarray:
            return D[np.flip(np.argsort(D.real))][:k]

        start = logg.info("Computing eigendecomposition of the transition matrix")

        if issparse(self.transition_matrix):
            logg.debug(f"Computing top `{k}` eigenvalues of a sparse matrix")
            D, V_l = eigs(self.transition_matrix.T, k=k, which=which, ncv=ncv)
            if only_evals:
                self._write_eigendecomposition(
                    {
                        "D": get_top_k_evals(),
                        "eigengap": _eigengap(get_top_k_evals().real, alpha),
                        "params": {"which": which, "k": k, "alpha": alpha},
                    },
                    time=start,
                )
                return
            _, V_r = eigs(self.transition_matrix, k=k, which=which, ncv=ncv)
        else:
            logg.warning(
                "This transition matrix is not sparse, computing full eigendecomposition"
            )
            D, V_l = np.linalg.eig(self.transition_matrix.T)
            if only_evals:
                self._write_eigendecomposition(
                    {
                        "D": get_top_k_evals(),
                        "eigengap": _eigengap(D.real, alpha),
                        "params": {"which": which, "k": k, "alpha": alpha},
                    },
                    time=start,
                )
                return
            _, V_r = np.linalg.eig(self.transition_matrix)

        # Sort the eigenvalues and eigenvectors and take the real part
        logg.debug("Sorting eigenvalues by their real part")
        p = np.flip(np.argsort(D.real))
        D, V_l, V_r = D[p], V_l[:, p], V_r[:, p]
        e_gap = _eigengap(D.real, alpha)

        pi = np.abs(V_l[:, 0].real)
        pi /= np.sum(pi)

        self._write_eigendecomposition(
            {
                "D": D,
                "stationary_dist": pi,
                "V_l": V_l,
                "V_r": V_r,
                "eigengap": e_gap,
                "params": {"which": which, "k": k, "alpha": alpha},
            },
            time=start,
        )

    @d.dedent
    def plot_eigendecomposition(
        self, *args: Any, left: bool = False, **kwargs: Any
    ) -> None:
        """
        Plot eigenvectors in an embedding.

        Parameters
        ----------
        %(plot_vectors.parameters)s
        left
            Whether to plot left or right eigenvectors.

        Returns
        -------
        %(plot_vectors.returns)s
        """

        eig = self.eigendecomposition
        if eig is None:
            raise RuntimeError(
                "Compute `.eigendecomposition` first as `.compute_eigendecomposition()`."
            )

        side = "left" if left else "right"
        D, V = (
            eig["D"],
            eig.get(f"V_{side[0]}", None),
        )
        if V is None:
            raise RuntimeError(
                "Compute `.eigendecomposition` first as `.compute_eigendecomposition(..., only_evals=False)`."
            )

        # if irreducible, first right e-vec should be const.
        if not left:
            # quick check for irreducibility:
            if np.sum(np.isclose(D, 1, rtol=1e2 * EPS, atol=1e2 * EPS)) == 1:
                V[:, 0] = 1.0

        self._plot_vectors(
            "eigen",
            V,
            *args,
            D=D,
            **kwargs,
        )

    @d.dedent
    def plot_spectrum(
        self,
        n: Optional[int] = None,
        real_only: bool = False,
        show_eigengap: bool = True,
        show_all_xticks: bool = True,
        legend_loc: Optional[str] = None,
        title: Optional[str] = None,
        marker: str = ".",
        figsize: Optional[Tuple[float, float]] = (5, 5),
        dpi: int = 100,
        save: Optional[Union[str, Path]] = None,
        **kwargs: Any,
    ) -> None:
        """
        Plot the top eigenvalues in real or complex plane.

        Parameters
        ----------
        n
            Number of eigenvalues to show. If `None`, show all that have been computed.
        real_only
            Whether to plot only the real part of the spectrum.
        show_eigengap
            When ``real_only = True``, this determines whether to show the inferred eigengap as
            a dotted line.
        show_all_xticks
            When ``real_only = True``, this determines whether to show the indices of all eigenvalues on the x-axis.
        legend_loc
            Location parameter for the legend.
        title
            Title of the figure.
        marker
            Marker symbol used, valid options can be found in :mod:`matplotlib.markers`.
        %(plotting)s
        kwargs
            Keyword arguments for :func:`matplotlib.pyplot.scatter`.

        Returns
        -------
        %(just_plots)s
        """

        eig = self.eigendecomposition
        if eig is None:
            raise RuntimeError(
                "Compute `.eigendecomposition` first as `.compute_eigendecomposition()`."
            )

        if n is None:
            n = len(eig["D"])
        if n <= 0:
            raise ValueError(f"Expected `n` to be > 0, found `{n}`.")

        if real_only:
            fig = self._plot_real_spectrum(
                n,
                show_eigengap=show_eigengap,
                show_all_xticks=show_all_xticks,
                dpi=dpi,
                figsize=figsize,
                legend_loc=legend_loc,
                title=title,
                marker=marker,
                **kwargs,
            )
        else:
            fig = self._plot_complex_spectrum(
                n,
                dpi=dpi,
                figsize=figsize,
                legend_loc=legend_loc,
                title=title,
                marker=marker,
                **kwargs,
            )

        if save:
            save_fig(fig, save)

    def _plot_complex_spectrum(
        self,
        n: int,
        dpi: int = 100,
        figsize: Optional[Tuple[float, float]] = (None, None),
        legend_loc: Optional[str] = None,
        title: Optional[str] = None,
        marker: str = ".",
        **kwargs: Any,
    ):
        # define a function to make the data limits rectangular
        def adapt_range(min_, max_, range_):
            return (
                min_ + (max_ - min_) / 2 - range_ / 2,
                min_ + (max_ - min_) / 2 + range_ / 2,
            )

        eig = self.eigendecomposition
        D, params = eig["D"][:n], eig["params"]

        # create fiture and axes
        fig, ax = plt.subplots(nrows=1, ncols=1, dpi=dpi, figsize=figsize)

        # get the original data ranges
        lam_x, lam_y = D.real, D.imag
        x_min, x_max = np.min(lam_x), np.max(lam_x)
        y_min, y_max = np.min(lam_y), np.max(lam_y)
        x_range, y_range = x_max - x_min, y_max - y_min
        final_range = np.max([x_range, y_range]) + 0.05

        x_min_, x_max_ = adapt_range(x_min, x_max, final_range)
        y_min_, y_max_ = adapt_range(y_min, y_max, final_range)

        # plot the data and the unit circle
        ax.scatter(D.real, D.imag, marker=marker, label="eigenvalue", **kwargs)
        t = np.linspace(0, 2 * np.pi, 500)
        x_circle, y_circle = np.sin(t), np.cos(t)
        ax.plot(x_circle, y_circle, "k-", label="unit circle")

        # set labels, ranges and legend
        ax.set_xlabel(r"Re($\lambda$)")
        ax.set_xlim(x_min_, x_max_)

        ax.set_ylabel(r"Im($\lambda$)")
        ax.set_ylim(y_min_, y_max_)

        key = "real part" if params["which"] == "LR" else "magnitude"
        if title is None:
            title = f"top {n} eigenvalues according to their {key}"

        ax.set_title(title)
        if legend_loc not in (None, "none"):
            ax.legend(loc=legend_loc)

        return fig

    def _plot_real_spectrum(
        self,
        n: int,
        show_eigengap: bool = True,
        show_all_xticks: bool = True,
        dpi: int = 100,
        figsize: Optional[Tuple[float, float]] = None,
        legend_loc: Optional[str] = None,
        title: Optional[str] = None,
        marker: str = ".",
        **kwargs,
    ):
        eig = self.eigendecomposition
        D, params = eig["D"][:n], eig["params"]

        D_real, D_imag = D.real, D.imag
        ixs = np.arange(len(D))
        mask = D_imag == 0

        # plot the top eigenvalues
        fig, ax = plt.subplots(nrows=1, ncols=1, dpi=dpi, figsize=figsize)
        if np.any(mask):
            ax.scatter(
                ixs[mask],
                D_real[mask],
                marker=marker,
                label="real eigenvalue",
                **kwargs,
            )
        if np.any(~mask):
            ax.scatter(
                ixs[~mask],
                D_real[~mask],
                marker=marker,
                label="complex eigenvalue",
                **kwargs,
            )

        # add dashed line for the eigengap, ticks, labels, title and legend
        if show_eigengap and eig["eigengap"] < n:
            ax.axvline(eig["eigengap"], label="eigengap", ls="--", lw=2, c="k")

        ax.set_xlabel("index")
        if show_all_xticks:
            ax.set_xticks(np.arange(len(D)))
        else:
            ax.xaxis.set_major_locator(MultipleLocator(2.0))
            ax.xaxis.set_major_formatter(FormatStrFormatter("%d"))

        ax.set_ylabel(r"Re($\lambda_i$)")

        key = "real part" if params["which"] == "LR" else "magnitude"
        if title is None:
            title = f"real part of top {n} eigenvalues according to their {key}"

        ax.set_title(title)
        if legend_loc not in (None, "none"):
            ax.legend(loc=legend_loc)

        return fig

    def _write_eigendecomposition(
        self: EigenProtocol, decomp: Dict[str, Any], time: datetime
    ) -> None:
        self._eigendecomposition = decomp

        key = Key.uns.eigen(self.backward)
        self.adata.uns[key] = decomp

        logg.info(
            f"Adding `adata.uns[{key!r}]`\n"
            f"       `.eigendecomposition`\n"
            "    Finish",
            time=time,
        )

    @property
    def eigendecomposition(self) -> Optional[Dict[str, Any]]:
        """TODO. docrep"""
        return self._eigendecomposition
