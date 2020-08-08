# -*- coding: utf-8 -*-
"""Precomputed kernel module."""
from copy import copy
from typing import Union, Optional

import numpy as np
from scipy.sparse import spmatrix, csr_matrix

from cellrank import logging as logg
from cellrank.utils._docs import d
from cellrank.utils._utils import _read_graph_data
from cellrank.tools.kernels import Kernel
from cellrank.tools._constants import Direction
from cellrank.tools.kernels._base_kernel import AnnData


@d.dedent
class PrecomputedKernel(Kernel):
    """
    Kernel which contains precomputed transition matrix.

    Parameters
    ----------
    transition_matrix
        Row-normalized transition matrix or a key in :paramref:`adata` `.obsp`.
    %(adata)s
    %(backward)s
    """

    def __init__(
        self,
        transition_matrix: Union[np.ndarray, spmatrix, str],
        adata: Optional[AnnData] = None,
        backward: bool = False,
        compute_cond_num: bool = False,
    ):
        from anndata import AnnData as _AnnData

        if isinstance(transition_matrix, str):
            if adata is None:
                raise ValueError(
                    "When `transition_matrix` specifies a key to `adata.obsp`, `adata` cannot be None."
                )
            transition_matrix = _read_graph_data(adata, transition_matrix)

        if not isinstance(transition_matrix, (np.ndarray, spmatrix)):
            raise TypeError(
                f"Expected transition matrix to be of type `numpy.ndarray` or `scipy.sparse.spmatrix`, "
                f"found `{type(transition_matrix).__name__!r}`."
            )

        if transition_matrix.shape[0] != transition_matrix.shape[1]:
            raise ValueError(
                f"Expected transition matrix to be square, found `{transition_matrix.shape}`."
            )

        if not np.allclose(np.sum(transition_matrix, axis=1), 1.0):
            raise ValueError("Not a valid transition matrix: not all rows sum to 1.")

        if adata is None:
            logg.debug("Creating empty dummy AnnData object")
            adata = _AnnData(
                csr_matrix((transition_matrix.shape[0], 1), dtype=np.float32)
            )

        super().__init__(adata, backward=backward, compute_cond_num=compute_cond_num)
        self._transition_matrix = csr_matrix(transition_matrix)
        self._maybe_compute_cond_num()

    def _read_from_adata(self, **kwargs):
        pass

    def copy(self) -> "PrecomputedKernel":
        """Return a copy of self."""
        pk = PrecomputedKernel(
            copy(self.transition_matrix), adata=self.adata, backward=self.backward
        )
        pk._cond_num = self.condition_number

        return pk

    def compute_transition_matrix(self, *args, **kwargs) -> "PrecomputedKernel":
        """Return self."""
        return self

    def __invert__(self) -> "PrecomputedKernel":
        # do not call parent's invert, since it removes the transition matrix
        self._direction = Direction.FORWARD if self.backward else Direction.BACKWARD
        return self

    def __repr__(self):
        return f"{'~' if self.backward and self._parent is None else ''}<Precomputed>"

    def __str__(self):
        return repr(self)
