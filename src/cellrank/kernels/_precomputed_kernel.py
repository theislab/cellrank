from typing import Any, Optional, Union

import numpy as np
import scipy.sparse as sp

from anndata import AnnData

from cellrank import logging as logg
from cellrank._utils._docs import d
from cellrank._utils._key import Key
from cellrank._utils._utils import _read_graph_data
from cellrank.kernels._base_kernel import KernelExpression, UnidirectionalKernel

__all__ = ["PrecomputedKernel"]


@d.dedent
class PrecomputedKernel(UnidirectionalKernel):
    """Kernel which is initialized based on a precomputed transition matrix.

    This kernel serves as CellRank's interface with other methods that compute cell-cell transition matrices; you
    can use this kernel to input you own custom transition matrix and continue to use all CellRank functionality.
    In particular, you can use a precomputed kernel, just like any other kernel, to initialize an estimator and
    compute initial and terminal states, fate probabilities, and driver genes.

    Parameters
    ----------
    object
        Can be one of the following types:

        - :class:`~anndata.AnnData` - annotated data object.
        - :class:`~numpy.ndarray`, :class:`~scipy.sparse.spmatrix` - row-normalized transition matrix.
        - :class:`~cellrank.kernels.Kernel` - kernel.
        - :class:`str` - key in :attr:`~anndata.AnnData.obsp` where the transition matrix is stored.
          ``adata`` must be provided in this case.
        - :class:`bool` - directionality of the transition matrix that will be used to infer its storage location.
          If :obj:`None`, the directionality will be determined automatically. ``adata`` must be provided in this case.
    %(adata)s
        Must be provided when ``object`` is :class:`str` or :class:`bool`.
    obsp_key
        Key in :attr:`~anndata.AnnData.obsp` where the transition matrix is stored.
        If :obj:`None`, it will be determined automatically. Only used when ``object`` is :class:`~anndata.AnnData`.
    copy
        Whether to copy the stored transition matrix.
    backward
        Hint whether this is a forward, backward or a unidirectional kernel.
        Only used when ``object`` is :class:`~anndata.AnnData`.

    Notes
    -----
    If ``object`` is :class:`~anndata.AnnData` and neither ``obsp_key`` nor ``backward`` is specified,
    default forward and backward are tried and first one is used.
    """

    _SENTINEL = object()

    def __init__(
        self,
        object: Union[str, bool, np.ndarray, sp.spmatrix, AnnData, KernelExpression],
        adata: Optional[AnnData] = None,
        obsp_key: Optional[str] = None,
        **kwargs: Any,
    ):
        if isinstance(object, AnnData):
            self._from_adata(object, obsp_key=obsp_key, **kwargs)
        elif isinstance(object, KernelExpression):
            self._from_kernel(object, copy=kwargs.get("copy", False))
        elif isinstance(object, (np.ndarray, sp.spmatrix)):
            self._from_matrix(object, adata=adata, **kwargs)
        elif isinstance(adata, AnnData):
            if isinstance(object, str):
                self._from_adata(adata, obsp_key=object, **kwargs)
            elif object is None or isinstance(object, bool):
                kwargs["backward"] = object
                self._from_adata(adata, **kwargs)
            else:
                raise ValueError("Unable to interpret the data.")
        else:
            raise TypeError(
                f"Expected `object` to be either `str`, `bool`, `numpy.ndarray`, "
                f"`scipy.sparse.spmatrix`, `AnnData` or `KernelExpression`, "
                f"found `{type(object).__name__}`."
            )

    def _from_adata(
        self,
        adata: AnnData,
        obsp_key: Optional[str] = None,
        backward: Optional[bool] = _SENTINEL,
        copy: bool = False,
    ) -> None:
        if obsp_key is None:
            if backward is PrecomputedKernel._SENTINEL:
                for backward in [False, True]:  # prefer `forward`
                    obsp_key = Key.uns.kernel(backward)
                    if obsp_key in adata.obsp:
                        break
                else:
                    raise ValueError(
                        "Unable to find transition matrix in `adata.obsp`. "
                        "Please specify `obsp_key=...` or `backward=...`."
                    )
            else:
                obsp_key = Key.uns.kernel(backward)
            logg.info(f"Using transition matrix from `adata.obsp[{obsp_key!r}]`")

        tmat = _read_graph_data(adata, obsp_key)

        if backward is PrecomputedKernel._SENTINEL:
            # not ideal, since None/False share the same key, prefer `forward`
            if obsp_key == Key.uns.kernel(bwd=False):
                backward = False
            elif obsp_key == Key.uns.kernel(bwd=True):
                backward = True
            else:
                backward = None
            logg.info(f"Setting directionality `backward={backward}`")

        self._from_matrix(tmat, adata=adata, backward=backward, copy=copy)
        self.params["origin"] = f"adata.obsp[{obsp_key!r}]"
        self._init_kwargs = {"obsp_key": obsp_key, "backward": backward}

    def _from_kernel(self, kernel: KernelExpression, copy: bool = False) -> None:
        if kernel.transition_matrix is None:
            raise RuntimeError("Compute transition matrix first as `.compute_transition_matrix()`.")
        self._from_matrix(
            kernel.transition_matrix,
            backward=kernel.backward,
            adata=kernel.adata,
            copy=copy,
        )
        self._params = kernel.params.copy()
        self.params["origin"] = repr(kernel)

    def _from_matrix(
        self,
        matrix: Union[np.ndarray, sp.spmatrix],
        adata: Optional[AnnData] = None,
        backward: Optional[bool] = None,
        copy: bool = False,
    ) -> None:
        # fmt: off
        if adata is None:
            logg.warning(f"Creating empty `AnnData` object of shape `{matrix.shape[0], 1}`")
            adata = AnnData(sp.csr_matrix((matrix.shape[0], 1)))
        super().__init__(adata)
        self._backward: Optional[bool] = backward
        self.transition_matrix = matrix.copy() if copy else matrix
        self.params["origin"] = "array"
        # fmt: on

    def compute_transition_matrix(self, *_: Any, **__: Any) -> "PrecomputedKernel":
        """Do nothing and return self."""
        return self

    @property
    def backward(self) -> Optional[bool]:
        """Direction of the process."""
        return self._backward
