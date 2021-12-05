from typing import Optional

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
                raise AttributeError(
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
        force_recompute: bool = False,
    ) -> None:
        # modified from: https://github.com/theislab/scvelo/blob/master/scvelo/tools/velocity_embedding.py

        self._key = Key.uns.kernel(self._kexpr.backward, key=key_added)
        self._basis = basis
        ukey = f"{self._key}_params"
        key = self._key + "_" + basis

        if not force_recompute and key in self._kexpr.adata.obsm:
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

    def plot(self, **kwargs) -> None:
        # TODO: catch errors
        scv.pl.velocity_embedding_stream(
            self._kexpr.adata, basis=self._basis, vkey=self._key, **kwargs
        )
