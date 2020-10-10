# -*- coding: utf-8 -*-
"""Precomputed kernel module."""
from copy import copy
from typing import Union, Optional

import numpy as np
from scipy.sparse import eye as speye
from scipy.sparse import spmatrix, csr_matrix

from cellrank import logging as logg
from cellrank.ul._docs import d
from cellrank.ul._utils import _read_graph_data
from cellrank.tl.kernels import Kernel
from cellrank.tl._constants import Direction, _transition
from cellrank.tl.kernels._base_kernel import _RTOL, AnnData, KernelExpression


@d.dedent
class PrecomputedKernel(Kernel):
    """
    Kernel which contains a precomputed transition matrix.

    Parameters
    ----------
    transition_matrix
        Row-normalized transition matrix or a key in :paramref:`adata` ``.obsp``
        or a :class:`cellrank.tl.kernels.KernelExpression` with the computed transition matrix.
        If `None`, try to determine the key based on ``backward``.
    %(adata)s
    %(backward)s
    """

    def __init__(
        self,
        transition_matrix: Optional[
            Union[np.ndarray, spmatrix, KernelExpression, str]
        ] = None,
        adata: Optional[AnnData] = None,
        backward: bool = False,
        compute_cond_num: bool = False,
    ):
        from anndata import AnnData as _AnnData

        self._origin = "'array'"
        params = {}

        if transition_matrix is None:
            transition_matrix = _transition(
                Direction.BACKWARD if backward else Direction.FORWARD
            )
            logg.debug(f"Setting transition matrix key to `{transition_matrix!r}`")

        if isinstance(transition_matrix, str):
            if adata is None:
                raise ValueError(
                    "When `transition_matrix` specifies a key to `adata.obsp`, `adata` cannot be None."
                )
            self._origin = f"adata.obsp[{transition_matrix!r}]"
            transition_matrix = _read_graph_data(adata, transition_matrix)

        elif isinstance(transition_matrix, KernelExpression):
            if transition_matrix._transition_matrix is None:
                raise ValueError(
                    "Compute transition matrix first as `.compute_transition_matrix()`."
                )
            if adata is not None and adata is not transition_matrix.adata:
                logg.warning(
                    "Ignoring supplied `adata` object because it differs from the kernel's `adata` object."
                )

            # use `str` because it captures the params
            self._origin = str(transition_matrix).strip("~<>")
            params = transition_matrix.params.copy()
            backward = transition_matrix.backward
            adata = transition_matrix.adata
            transition_matrix = transition_matrix.transition_matrix

        if not isinstance(transition_matrix, (np.ndarray, spmatrix)):
            raise TypeError(
                f"Expected transition matrix to be of type `numpy.ndarray` or `scipy.sparse.spmatrix`, "
                f"found `{type(transition_matrix).__name__!r}`."
            )

        if transition_matrix.shape[0] != transition_matrix.shape[1]:
            raise ValueError(
                f"Expected transition matrix to be square, found `{transition_matrix.shape}`."
            )

        if not np.allclose(np.sum(transition_matrix, axis=1), 1.0, rtol=_RTOL):
            raise ValueError("Not a valid transition matrix, not all rows sum to 1")

        if adata is None:
            logg.warning("Creating empty `AnnData` object")
            adata = _AnnData(
                csr_matrix((transition_matrix.shape[0], 1), dtype=np.float32)
            )

        super().__init__(adata, backward=backward, compute_cond_num=compute_cond_num)

        self._params = params
        self._transition_matrix = csr_matrix(transition_matrix)
        self._maybe_compute_cond_num()

    def _read_from_adata(self, **kwargs):
        pass

    @d.dedent
    def copy(self) -> "PrecomputedKernel":
        """%(copy)s"""  # noqa
        pk = PrecomputedKernel(
            copy(self.transition_matrix),
            adata=self.adata,
            backward=self.backward,
            compute_cond_num=False,
        )
        pk._cond_num = self.condition_number
        pk._origin = self._origin
        pk._params = self._params.copy()

        return pk

    def compute_transition_matrix(self, *args, **kwargs) -> "PrecomputedKernel":
        """Return self."""
        return self

    def __invert__(self) -> "PrecomputedKernel":
        # do not call parent's invert, since it removes the transition matrix
        self._direction = Direction.FORWARD if self.backward else Direction.BACKWARD
        return self

    def __repr__(self):
        return f"{'~' if self.backward and self._parent is None else ''}<Precomputed[origin={self._origin}]>"

    def __str__(self):
        return repr(self)


@d.dedent
class DummyKernel(PrecomputedKernel):
    """
    Kernel with 1s on the diagonal.

    Parameters
    ----------
    %(adata)s
    %(backward)s
    """

    def __init__(
        self,
        adata: AnnData,
        backward: bool = False,
    ):
        super().__init__(
            speye(adata.n_obs, format="csr"),
            adata,
            backward=backward,
            compute_cond_num=False,
        )
