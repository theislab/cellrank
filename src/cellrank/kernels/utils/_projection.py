import warnings
from typing import Any, Optional

import scvelo as scv
from scvelo.tools.velocity_embedding import quiver_autoscale

import numpy as np
import scipy.sparse as sp

from cellrank import logging as logg
from cellrank._utils._docs import d
from cellrank._utils._key import Key
from cellrank.kernels._utils import _get_basis

__all__ = ["TmatProjection"]


class TmatProjection:
    """Project transition matrix onto a low-dimensional embedding.

    Should be used for visualization purposes.

    Parameters
    ----------
    kexpr
        Kernel that contains a transition matrix.
    basis
        Key in :attr:`~anndata.AnnData.obsm` where the basis is stored.
    """

    def __init__(self, kexpr: "ExpernelExpression", basis: str = "umap"):  # noqa: F821
        from cellrank.kernels.mixins import ConnectivityMixin

        for kernel in kexpr.kernels:
            if not isinstance(kernel, ConnectivityMixin):
                logg.warning(
                    f"{kernel!r} is not a kNN based kernel. "
                    f"The embedding projection works best for kNN-based kernels"
                )
                break
        self._kexpr = kexpr
        self._basis = basis[2:] if basis.startswith("X_") else basis
        self._key: Optional[str] = None

    def project(
        self,
        key_added: Optional[str] = None,
        recompute: bool = False,
        connectivities: Optional[sp.spmatrix] = None,
    ) -> None:
        """Project transition matrix onto an embedding.

        This function is adapted from :func:`~scvelo.tl.velocity_embedding`.

        Parameters
        ----------
        key_added
            Key in :attr:`~anndata.AnnData.obsm` where to store the projection.
        recompute
            Whether to recompute an existing projection.
        connectivities
            Connectivity matrix to use for projection. If :obj:`None`, use ones from the underlying kernel, is possible.

        Returns
        -------
        Nothing, just updates :attr:`adata` with the projection and the parameters used for computation.
        """
        from cellrank.kernels.mixins import ConnectivityMixin

        self._key = Key.uns.kernel(self._kexpr.backward, key=key_added)
        ukey = f"{self._key}_params"
        key = f"{self._key}_{self._basis}"

        if not recompute and key in self._kexpr.adata.obsm:
            logg.info(f"Using precomputed projection `adata.obsm[{key!r}]`")
            return

        start = logg.info(f"Projecting transition matrix onto `{self._basis}`")
        if connectivities is None:
            try:
                connectivities, *_ = (c.connectivities for c in self._kexpr.kernels if isinstance(c, ConnectivityMixin))
            except ValueError:
                raise RuntimeError(
                    "Unable to find connectivities in the kernel. "
                    "Please supply them explicitly as `connectivities=...`."
                ) from None
        if not sp.isspmatrix_csr(connectivities):
            connectivities = connectivities.tocsr()

        emb = _get_basis(self._kexpr.adata, self._basis)
        T_emb = np.empty_like(emb)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for row_id, row in enumerate(self._kexpr.transition_matrix):
                conn_idxs = connectivities[row_id, :].indices
                dX = emb[conn_idxs] - emb[row_id, None]

                if np.any(np.isnan(dX)):
                    T_emb[row_id, :] = np.nan
                else:
                    probs = row[:, conn_idxs].A.squeeze() if sp.issparse(row) else row[conn_idxs]
                    dX /= np.linalg.norm(dX, axis=1)[:, None]
                    dX = np.nan_to_num(dX)
                    T_emb[row_id, :] = probs.dot(dX) - dX.sum(0) / dX.shape[0]

        T_emb /= 3 * quiver_autoscale(np.nan_to_num(emb), T_emb)

        embs = self._kexpr.adata.uns.get(ukey, {}).get("embeddings", [])
        if self._basis not in embs:
            embs = list(embs) + [self._basis]
            self._kexpr.adata.uns[ukey] = self._kexpr.adata.uns.get(ukey, {})
            self._kexpr.adata.uns[ukey]["embeddings"] = embs

        logg.info(
            f"Adding `adata.obsm[{key!r}]`\n    Finish",
            time=start,
        )
        self._kexpr.adata.obsm[key] = T_emb

    @d.dedent
    def plot(self, *args: Any, stream: bool = True, **kwargs: Any) -> None:
        """Plot projected transition matrix in a embedding.

        Parameters
        ----------
        args
            Positional argument for the plotting function.
        stream
            If :obj:`True``, use :func:`~scvelo.pl.velocity_embedding_stream`.
            Otherwise, use :func:`~scvelo.pl.velocity_embedding_grid`.
        kwargs
            Keyword argument for the above-mentioned plotting function.

        Returns
        -------
        %(just_plots)s
        """
        if stream:
            return scv.pl.velocity_embedding_stream(
                self._kexpr.adata, *args, basis=self._basis, vkey=self._key, **kwargs
            )
        return scv.pl.velocity_embedding_grid(self._kexpr.adata, *args, basis=self._basis, vkey=self._key, **kwargs)
