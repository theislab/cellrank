from typing import Any, Dict, Tuple, Union, Optional
from typing_extensions import Literal, Protocol

from pathlib import Path
from pygpcca import GPCCA
from datetime import datetime
from pygpcca._sorted_schur import _check_conj_split

from anndata import AnnData
from cellrank import logging as logg
from cellrank.ul._docs import d
from cellrank.tl._utils import save_fig, _eigengap
from cellrank.tl._estimators.mixins._constants import Key
from cellrank.tl._estimators.mixins.decomposition._plot import VectorPlottable

import numpy as np
from scipy.sparse import issparse, csr_matrix

import matplotlib.pyplot as plt
from seaborn import heatmap
from mpl_toolkits.axes_grid1 import make_axes_locatable

EPS = np.finfo(np.float64).eps


class SchurProtocol(Protocol):
    @property
    def adata(self) -> AnnData:
        ...

    @property
    def backward(self) -> bool:
        ...

    @property
    def transition_matrix(self) -> Union[csr_matrix, np.ndarray]:
        ...

    def _write_schur_decomposition(
        self,
        decomp: Dict[str, Any],
        vectors: np.ndarray,
        matrix: np.ndarray,
        time: datetime,
    ) -> None:
        ...


class SchurMixin(VectorPlottable):
    """TODO."""

    def __init__(self, **kwargs: Any):
        super().__init__(**kwargs)

        self._gpcca: Optional[GPCCA] = None
        self._invalid_n_states: Optional[np.ndarray] = None

        self._schur_vectors: Optional[np.ndarray] = None
        self._schur_matrix: Optional[np.ndarray] = None
        self._eigendecomposition: Optional[Dict[str, Any]] = None

    @d.dedent
    def compute_schur(
        self: SchurProtocol,
        n_components: int = 10,
        initial_distribution: Optional[np.ndarray] = None,
        method: Literal["krylov", "brandts"] = "krylov",
        which: Literal["LR", "LM"] = "LR",
        alpha: float = 1.0,
    ):
        """
        Compute the Schur decomposition.

        Parameters
        ----------
        n_components
            Number of vectors to compute.
        initial_distribution
            Input probability distribution over all cells. If `None`, uniform is used.
        method
            Method for calculating the Schur vectors. Valid options are:

                - `'krylov'` - an iterative procedure that computes a partial, sorted Schur decomposition for
                  large, sparse matrices.
                - `'brandts'` - full sorted Schur decomposition of a dense matrix.

            For benefits of each method, see :class:`pygpcca.GPCCA`.
        %(eigen)s

        Returns
        -------
        None
            Nothing, but updates the following fields:

                - :attr:`schur_vectors` - TODO.
                - :attr:`schur_matrix` -  TODO.
                - :attr:`eigendecomposition` - TODO.
        """
        if n_components < 2:
            raise ValueError(
                f"Number of components must be `>=2`, found `{n_components}`."
            )

        if method not in ("brandts", "krylov"):
            raise ValueError(
                f"Invalid method `{method!r}`. Valid options are `'brandts'` or `'krylov'`."
            )

        try:
            import petsc4py
            import slepc4py
        except ImportError:
            method = "brandts"
            logg.warning(
                f"Unable to import `petsc4py` or `slepc4py`. Using `method={method!r}`"
            )

        tmat = self.transition_matrix
        if method == "brandts" and issparse(self.transition_matrix):
            logg.warning("For `method='brandts'`, dense matrix is required. Densifying")
            tmat = tmat.A

        self._gpcca = GPCCA(tmat, eta=initial_distribution, z=which, method=method)
        start = logg.info("Computing Schur decomposition")

        try:
            self._gpcca._do_schur_helper(n_components)
        except ValueError:
            logg.warning(
                f"Using `{n_components}` components would split a block of complex conjugates. "
                f"Increasing `n_components` to `{n_components + 1}`"
            )
            self._gpcca._do_schur_helper(n_components + 1)

        self._invalid_n_states = np.array(
            [
                i
                for i in range(2, len(self._gpcca._p_eigenvalues))
                if _check_conj_split(self._gpcca._p_eigenvalues[:i])
            ]
        )
        if len(self._invalid_n_states):
            logg.info(
                f"When computing macrostates, choose a number of states NOT in `{list(self._invalid_n_states)}`"
            )

        self._write_schur_decomposition(
            {
                "D": self._gpcca._p_eigenvalues,
                "eigengap": _eigengap(self._gpcca._p_eigenvalues, alpha),
                "params": {
                    "which": which,
                    "k": len(self._gpcca._p_eigenvalues),
                    "alpha": alpha,
                },
            },
            vectors=self._gpcca._p_X,
            matrix=self._gpcca._p_R,
            time=start,
        )

    @d.dedent
    def plot_schur_matrix(
        self,
        title: Optional[str] = "schur matrix",
        cmap: str = "viridis",
        figsize: Optional[Tuple[float, float]] = None,
        dpi: Optional[float] = 80,
        save: Optional[Union[str, Path]] = None,
        **kwargs: Any,
    ):
        """
        Plot the Schur matrix.

        Parameters
        ----------
        title
            Title of the figure.
        cmap
            Colormap to use.
        %(plotting)s
        kwargs
            Keyword arguments for :func:`seaborn.heatmap`.

        Returns
        -------
        %(just_plots)s
        """
        schur_matrix = self.schur_matrix
        if schur_matrix is None:
            raise RuntimeError("Compute `.schur_matrix` first as `.compute_schur()`.")

        fig, ax = plt.subplots(
            figsize=schur_matrix.shape if figsize is None else figsize, dpi=dpi
        )

        divider = make_axes_locatable(ax)  # square=True make the colorbar a bit bigger
        cbar_ax = divider.append_axes("right", size="2%", pad=0.1)

        mask = np.zeros_like(schur_matrix, dtype=np.bool_)
        mask[np.tril_indices_from(mask, k=-1)] = True
        mask[~np.isclose(schur_matrix, 0.0)] = False

        vmin, vmax = (
            np.min(schur_matrix[~mask]),
            np.max(schur_matrix[~mask]),
        )

        kwargs["fmt"] = kwargs.get("fmt", "0.2f")
        heatmap(
            data=schur_matrix,
            cmap=cmap,
            square=True,
            annot=True,
            vmin=vmin,
            vmax=vmax,
            cbar_ax=cbar_ax,
            cbar_kws={"ticks": np.linspace(vmin, vmax, 10)},
            mask=mask,
            xticklabels=[],
            yticklabels=[],
            ax=ax,
            **kwargs,
        )
        ax.set_title(title)

        if save is not None:
            save_fig(fig, path=save)

    def _write_schur_decomposition(
        self: SchurProtocol,
        decomp: Dict[str, Any],
        vectors: np.ndarray,
        matrix: np.ndarray,
        time: datetime,
    ) -> None:
        self._eigendecomposition = decomp
        self._schur_vectors = vectors
        self._schur_matrix = matrix

        key = Key.uns.eigen(self.backward)
        self.adata.uns[key] = decomp

        logg.info(
            f"Adding `adata.uns[{key!r}]`\n"
            f"       `.schur_vectors`\n"
            f"       `.schur_matrix`\n"
            f"       `.eigendecomposition`\n"
            "    Finish",
            time=time,
        )

    @d.dedent
    def plot_schur(self, *args: Any, **kwargs: Any) -> None:
        """
        Plot Schur vectors in an embedding.

        Parameters
        ----------
        %(plot_vectors.parameters)s

        Returns
        -------
        %(plot_vectors.returns)s
        """

        vectors = self.schur_vectors
        if vectors is None:
            raise RuntimeError("Compute `.schur_vectors` first as `.compute_schur()`.")

        self._plot_vectors(
            "schur",
            vectors,
            *args,
            **kwargs,
        )

    # TODO: use docrep
    @property
    def schur_vectors(self) -> Optional[np.ndarray]:
        """TODO."""
        return self._schur_vectors

    @property
    def schur_matrix(self) -> Optional[np.ndarray]:
        """TODO."""
        return self._schur_matrix

    @property
    def eigendecomposition(self) -> Optional[Dict[str, Any]]:
        """TODO."""
        return self._eigendecomposition
