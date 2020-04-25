# -*- coding: utf-8 -*-
from cellrank.tools.kernels._kernel import KernelExpression
from typing import Optional

from abc import ABC
import numpy as np

from anndata import AnnData
from scanpy import logging as logg
from scipy.sparse import issparse


from cellrank.tools._lineage import Lineage
from cellrank.tools._constants import (
    Direction,
    RcKey,
    LinKey,
    Prefix,
    _probs,
    _colors,
    _lin_names,
)
from cellrank.tools._utils import _create_categorical_colors


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

        self._is_irreducible = None
        self._rec_classes = None
        self._trans_classes = None

        # for copy
        self._g2m_key = g2m_key
        self._s_key = s_key
        self._key_added = key_added

        # read eig, approx_rcs and lin_probs from adata if present
        self._eig, self._approx_rcs, self._approx_rcs_colors, self._lin_probs, self._dp, self._G2M_score, self._S_score, self._approx_rcs_probs = (
            [None] * 8
        )

        if read_from_adata:
            logg.debug(
                "DEBUG: Reading `eig`, `approx_rcs` and `lin_probs` from `adata` object"
            )
            self._read_from_adata(g2m_key, s_key)

    def _read_from_adata(
        self, g2m_key: Optional[str] = None, s_key: Optional[str] = None
    ) -> None:
        if f"eig_{self._direction}" in self._adata.uns.keys():
            self._eig = self._adata.uns[f"eig_{self._direction}"]
        else:
            logg.debug(
                f"DEBUG: `eig_{self._direction}` not found. Setting `.eig` to `None`"
            )

        if self._rc_key in self._adata.obs.keys():
            self._approx_rcs = self._adata.obs[self._rc_key]
        else:
            logg.debug(
                f"DEBUG: `{self._rc_key}` not found in `adata.obs`. Setting `.approx_rcs` to `None`"
            )

        if _colors(self._rc_key) in self._adata.uns.keys():
            self._approx_rcs_colors = self._adata.uns[_colors(self._rc_key)]
        else:
            logg.debug(
                f"DEBUG: `{_colors(self._rc_key)}` not found in `adata.uns`. "
                f"Setting `.approx_rcs_colors`to `None`"
            )

        if self._lin_key in self._adata.obsm.keys():
            lineages = range(self._adata.obsm[self._lin_key].shape[1])
            colors = _create_categorical_colors(len(lineages))
            self._lin_probs = Lineage(
                self._adata.obsm[self._lin_key],
                names=[f"Lineage {i + 1}" for i in lineages],
                colors=colors,
            )
            self._adata.obsm[self._lin_key] = self._lin_probs
        else:
            logg.debug(
                f"DEBUG: `{self._lin_key}` not found in `adata.obsm`. Setting `.lin_probs` to `None`"
            )

        if f"{self._lin_key}_dp" in self._adata.obs.keys():
            self._dp = self._adata.obs[f"{self._lin_key}_dp"]
        else:
            logg.debug(
                f"DEBUG: `{self._lin_key}_dp` not found in `adata.obs`. Setting `.dp` to `None`"
            )

        if g2m_key and g2m_key in self._adata.obs.keys():
            self._G2M_score = self._adata.obs[g2m_key]
        else:
            logg.debug(
                f"DEBUG: `{g2m_key}` not found in `adata.obs`. Setting `.G2M_score` to `None`"
            )

        if s_key and s_key in self._adata.obs.keys():
            self._S_score = self._adata.obs[s_key]
        else:
            logg.debug(
                f"DEBUG: `{s_key}` not found in `adata.obs`. Setting `.S_score` to `None`"
            )

        if _probs(self._rc_key) in self._adata.obs.keys():
            self._approx_rcs_probs = self._adata.obs[_probs(self._rc_key)]
        else:
            logg.debug(
                f"DEBUG: `{_probs(self._rc_key)}` not found in `adata.obs`. "
                f"Setting `.approx_rcs_probs` to `None`"
            )

        if self._lin_probs is not None:
            if _lin_names(self._lin_key) in self._adata.uns.keys():
                self._lin_probs = Lineage(
                    np.array(self._lin_probs),
                    names=self._adata.uns[_lin_names(self._lin_key)],
                    colors=self._lin_probs.colors,
                )
                self._adata.obsm[self._lin_key] = self._lin_probs
            else:
                logg.debug(
                    f"DEBUG: `{_lin_names(self._lin_key)}` not found in `adata.uns`. "
                    f"Using default names"
                )

            if _colors(self._lin_key) in self._adata.uns.keys():
                self._lin_probs = Lineage(
                    np.array(self._lin_probs),
                    names=self._lin_probs.names,
                    colors=self._adata.uns[_colors(self._lin_key)],
                )
                self._adata.obsm[self._lin_key] = self._lin_probs
            else:
                logg.debug(
                    f"DEBUG: `{_colors(self._lin_key)}` not found in `adata.uns`. "
                    f"Using default colors"
                )
