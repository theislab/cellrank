import contextlib
import io
import pathlib
import types
from typing import Any, Dict, Literal, Mapping, Optional, Tuple, Union

from pygpcca import GPCCA
from pygpcca._sorted_schur import _check_conj_split

import numpy as np
import scipy.sparse as sp

import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable

from anndata import AnnData

from cellrank import logging as logg
from cellrank._utils._docs import d
from cellrank._utils._key import Key
from cellrank._utils._utils import _eigengap, save_fig
from cellrank.estimators.mixins._utils import BaseProtocol, SafeGetter, logger, shadow

__all__ = ["SchurMixin"]

EPS = np.finfo(np.float64).eps


class SchurProtocol(BaseProtocol):
    _gpcca: Optional[GPCCA] = None
    _invalid_n_states: Optional[np.ndarray] = None

    schur_vectors: Optional[np.ndarray] = None
    _schur_matrix: Optional[np.ndarray] = None
    _eigendecomposition: Optional[Dict[str, Any]] = None

    def transition_matrix(self) -> Union[np.ndarray, sp.spmatrix]:
        ...

    def _write_schur_decomposition(
        self,
        decomp: Dict[str, Any],
        vectors: np.ndarray,
        matrix: np.ndarray,
        params: Mapping[str, Any] = types.MappingProxyType({}),
    ) -> str:
        ...


class SchurMixin:
    """Mixin that computes Schur decomposition."""

    def __init__(self, **kwargs: Any):
        super().__init__(**kwargs)
        self._gpcca: Optional[GPCCA] = None
        self._invalid_n_states: Optional[np.ndarray] = None

        self._schur_vectors: Optional[np.ndarray] = None
        self._schur_matrix: Optional[np.ndarray] = None
        self._eigendecomposition: Optional[Dict[str, Any]] = None

    @property
    @d.get_summary(base="schur_vectors")
    @d.get_extended_summary(base="schur_vectors")
    def schur_vectors(self) -> Optional[np.ndarray]:
        """Real Schur vectors of the transition matrix.

        The real Schur decomposition is a generalization of the Eigendecomposition and can be computed for any
        real-valued, square matrix :math:`A`. It is given by :math:`A = Q R Q^T`, where :math:`Q` contains the
        real Schur vectors and :math:`R` is the Schur matrix. :math:`Q` is orthogonal and :math:`R` is quasi-upper
        triangular with 1x1 and 2x2 blocks on the diagonal.

        If PETSc and SLEPc are installed, only the leading Schur vectors are computed.
        """
        return self._schur_vectors

    @property
    @d.get_summary(base="schur_matrix")
    @d.dedent
    def schur_matrix(self) -> Optional[np.ndarray]:
        """Schur matrix.

        %(schur_vectors.summary_ext)s
        """
        return self._schur_matrix

    @property
    @d.get_summary(base="eigen")
    def eigendecomposition(self) -> Optional[Dict[str, Any]]:
        """Eigendecomposition of the :attr:`transition_matrix`.

        For non-symmetric real matrices, left and right eigenvectors will in general be different and complex.
        We compute both left and right eigenvectors.

        Returns
        -------
        A dictionary with the following keys:

        - ``'D'`` - the eigenvalues.
        - ``'eigengap'`` - the `eigengap <https://en.wikipedia.org/wiki/Eigengap>`__.
        - ``'params'`` - parameters used for the computation.
        - ``'V_l'`` - left eigenvectors (optional).
        - ``'V_r'`` - right eigenvectors (optional).
        - ``'stationary_dist'`` - stationary distribution of the :attr:`transition_matrix`, if present.
        """
        return self._eigendecomposition

    @d.dedent
    def compute_schur(
        self: SchurProtocol,
        n_components: int = 20,
        initial_distribution: Optional[np.ndarray] = None,
        method: Literal["krylov", "brandts"] = "krylov",
        which: Literal["LR", "LM"] = "LR",
        alpha: float = 1.0,
        verbose: Optional[bool] = None,
    ) -> "SchurMixin":
        """Compute the Schur decomposition.

        Parameters
        ----------
        n_components
            Number of Schur vectors to compute.
        initial_distribution
            Input distribution over all cells. If :obj:`None`, uniform distribution is used.
        method
            Method for calculating the Schur vectors. Valid options are:

            - ``'krylov'`` - an iterative procedure that computes a partial, sorted Schur decomposition for
              large, sparse matrices.
            - ``'brandts'`` - full sorted Schur decomposition of a dense matrix.

            For benefits of each method, see :class:`~pygpcca.GPCCA`.
        %(eigen)s
        verbose
            Whether to print extra information when computing the Schur decomposition.
            If :obj:`None`, it's disabled when ``method = 'krylov'``.

        Returns
        -------
        Self and just updates the following fields:

        - :attr:`schur_vectors` - %(schur_vectors.summary)s
        - :attr:`schur_matrix` -  %(schur_matrix.summary)s
        - :attr:`eigendecomposition` - %(eigen.summary)s
        """
        from cellrank._utils._linear_solver import _is_petsc_slepc_available

        if n_components < 2:
            logg.warning(
                f"Number of Schur vectors `>=2`, but only `{n_components}` " f"were requested. Using `n_components=2`"
            )
            n_components = 2

        if method not in ("brandts", "krylov"):
            raise ValueError(f"Invalid method `{method!r}`. Valid options are:`'brandts'` or `'krylov'`.")

        if not _is_petsc_slepc_available():
            method = "brandts"
            logg.warning(f"Unable to import `petsc4py` or `slepc4py`. Using `method={method!r}`")
        if verbose is None:
            verbose = method == "brandts"

        tmat = self.transition_matrix
        if method == "brandts" and sp.issparse(self.transition_matrix):
            logg.warning("For `method='brandts'`, dense matrix is required. Densifying")
            tmat = tmat.A

        self._gpcca = GPCCA(tmat, eta=initial_distribution, z=which, method=method)
        start = logg.info("Computing Schur decomposition")

        with (contextlib.nullcontext if verbose else contextlib.redirect_stdout)(io.StringIO()):
            try:
                self._gpcca._do_schur_helper(n_components)
            except ValueError as e:
                if "will split complex conjugate eigenvalues" not in str(e):
                    raise
                logg.warning(
                    f"Using `{n_components}` components would split a block of complex conjugate eigenvalues. "
                    f"Using `n_components={n_components + 1}`"
                )
                self._gpcca._do_schur_helper(n_components + 1)

        self._invalid_n_states = np.array(
            [i for i in range(2, len(self._gpcca._p_eigenvalues)) if _check_conj_split(self._gpcca._p_eigenvalues[:i])]
        )
        if len(self._invalid_n_states):
            logg.info(f"When computing macrostates, choose a number of states NOT in `{list(self._invalid_n_states)}`")

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
            params=self._create_params(),
            time=start,
        )
        return self

    @d.dedent
    def plot_schur_matrix(
        self,
        title: Optional[str] = "schur matrix",
        cmap: str = "viridis",
        figsize: Optional[Tuple[float, float]] = None,
        dpi: Optional[float] = 80,
        save: Optional[Union[str, pathlib.Path]] = None,
        **kwargs: Any,
    ) -> None:
        """Plot the Schur matrix.

        Parameters
        ----------
        title
            Title of the figure.
        cmap
            Colormap to use.
        %(plotting)s
        kwargs
            Keyword arguments for :func:`~seaborn.heatmap`.

        Returns
        -------
        %(just_plots)s
        """
        schur_matrix = self.schur_matrix
        if schur_matrix is None:
            raise RuntimeError("Compute Schur matrix first as `.compute_schur()`.")

        fig, ax = plt.subplots(figsize=schur_matrix.shape if figsize is None else figsize, dpi=dpi)

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
        sns.heatmap(
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

    @logger
    @shadow
    def _write_schur_decomposition(
        self: SchurProtocol,
        decomp: Dict[str, Any],
        vectors: np.ndarray,
        matrix: np.ndarray,
        params: Mapping[str, Any] = types.MappingProxyType({}),
    ) -> str:
        key = Key.uns.schur_matrix(self.backward)
        self._set("_schur_matrix", self.adata.uns, key=key, value=matrix)
        key = Key.obsm.schur_vectors(self.backward)
        self._set("_schur_vectors", self.adata.obsm, key=key, value=vectors)
        key = Key.uns.eigen(self.backward)
        self._set("_eigendecomposition", self.adata.uns, key=key, value=decomp)
        self.params[f"schur_decomposition_{Key.backward(self.backward)}"] = dict(params)

        return (
            f"Adding `adata.uns[{key!r}]`\n"
            f"       `.schur_vectors`\n"
            f"       `.schur_matrix`\n"
            f"       `.eigendecomposition`\n"
            "    Finish"
        )

    def _read_schur_decomposition(self: SchurProtocol, adata: AnnData, allow_missing: bool = True) -> bool:
        key = Key.uns.eigen(self.backward)
        with SafeGetter(self, allowed=KeyError) as sg:
            self._eigendecomposition = self._get(
                obj=adata.uns,
                key=key,
                shadow_attr="uns",
                dtype=dict,
                allow_missing=allow_missing,
            )
            key = Key.obsm.schur_vectors(self.backward)
            self._schur_vectors = self._get(
                obj=self.adata.obsm,
                key=key,
                shadow_attr="obsm",
                dtype=np.ndarray,
                allow_missing=allow_missing,
            )
            key = Key.uns.schur_matrix(self.backward)
            self._schur_matrix = self._get(
                obj=self.adata.uns,
                key=key,
                shadow_attr="uns",
                dtype=np.ndarray,
                allow_missing=allow_missing,
            )
            key = f"schur_decomposition_{Key.backward(self.backward)}"
            self.params[key] = self._read_params(key)

        return sg.ok
