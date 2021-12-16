from typing import Any, Optional

import warnings

import scvelo as scv
from cellrank import logging as logg
from cellrank._key import Key
from cellrank.tl.kernels._utils import _get_basis
from scvelo.tools.velocity_embedding import quiver_autoscale

import numpy as np
from scipy.sparse import issparse


class Projector:
    def __init__(self, kexpr: "ExpernelExpression"):
        from cellrank.tl.kernels._mixins import ConnectivityMixin

        for kernel in kexpr.kernels:
            if not isinstance(kernel, ConnectivityMixin):
                raise TypeError(
                    f"{kernel!r} is not a kNN based kernel. The embedding projection "
                    "only works for kNN based kernels."
                )
        self._kexpr = kexpr
        self._key: Optional[str] = None
        self._basis: Optional[str] = None

    def project(
        self,
        basis: str = "umap",
        key_added: Optional[str] = None,
        recompute: bool = False,
    ) -> None:
        """
        Project transition matrix onto an embedding.

        Parameters
        ----------
        basis
            TODO
        key_added
            TODO
        recompute
            Whether to recompute the projection if it already exists.

        Returns
        -------
        Nothing, just updates TODO.
        """
        # modified from: https://github.com/theislab/scvelo/blob/master/scvelo/tools/velocity_embedding.py

        self._key = Key.uns.kernel(self._kexpr.backward, key=key_added)
        self._basis = basis
        ukey = f"{self._key}_params"
        key = self._key + "_" + basis

        if not recompute and key in self._kexpr.adata.obsm:
            logg.debug(f"Using precomputed projection `adata.obsm[{key!r}]`")
            return

        start = logg.info(f"Projecting transition matrix onto `{basis}`")
        emb = _get_basis(self._kexpr.adata, basis)
        T_emb = np.empty_like(emb)

        conn = self._kexpr.kernels[0]._conn

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for row_id, row in enumerate(self._kexpr.transition_matrix):
                conn_idxs = conn[row_id, :].indices

                dX = emb[conn_idxs] - emb[row_id, None]

                if np.any(np.isnan(dX)):
                    T_emb[row_id, :] = np.nan
                else:
                    probs = row[:, conn_idxs]
                    if issparse(probs):
                        probs = probs.A.squeeze()

                    dX /= np.linalg.norm(dX, axis=1)[:, None]
                    dX = np.nan_to_num(dX)
                    T_emb[row_id, :] = probs.dot(dX) - dX.sum(0) / dX.shape[0]

        T_emb /= 3 * quiver_autoscale(np.nan_to_num(emb), T_emb)

        embs = self._kexpr.adata.uns.get(ukey, {}).get("embeddings", [])
        if basis not in embs:
            embs = list(embs) + [basis]
            self._kexpr.adata.uns[ukey] = self._kexpr.adata.uns.get(ukey, {})
            self._kexpr.adata.uns[ukey]["embeddings"] = embs

        logg.info(
            f"Adding `adata.obsm[{key!r}]`\n    Finish",
            time=start,
        )
        self._kexpr.adata.obsm[key] = T_emb

    def plot(self, *args: Any, stream: bool = True, **kwargs: Any) -> None:
        """
        Plot projected transition matrix in a embedding.

        Parameters
        ----------
        args
            Positional argument for the plotting function.
        stream
            If ``True``, use :func:`scvelo.pl.velocity_emebedding_stream`.
            Otherwise, use :func:`scvelo.pl.velocity_embedding_grid`.
        kwargs
            Keyword argument for the plotting function.
        """
        if stream:
            return scv.pl.velocity_embedding_stream(
                self._kexpr.adata, *args, basis=self._basis, vkey=self._key, **kwargs
            )
        return scv.pl.velocity_embedding_grid(
            self._kexpr.adata, *args, basis=self._basis, vkey=self._key, **kwargs
        )
