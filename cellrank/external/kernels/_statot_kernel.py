from typing import Any, Union, Optional

from typing_extensions import Literal

from anndata import AnnData

import numpy as np

from cellrank.ul._docs import d

_error = None
try:
    from statot import OTKernel as OTKernel_
except ImportError as e:
    from cellrank.external.kernels._import_error_kernel import (
        ErroredKernel as OTKernel_,
    )

    _error = e


@d.dedent
class OTKernel(OTKernel_, error=_error):
    """
    Stationary optimal transport kernel from [Zhang21]_.

    This class requires the `statOT` package, which can be installed as `pip install statot POT`.

    Parameters
    ----------
    %(adata)s
    source_idx
        Key in :attr:`anndata.AnnData.obs` containing boolean mask marking the source indices or the mask itself.
    sink_idx
        Key in :attr:`anndata.AnnData.obs` containing boolean marking the sink indices or the mask itself.
    g
        Key in :attr:`anndata.AnnData.obs` containing relative growth rates for cells or the array itself.
    kwargs
        Additional keyword arguments.
    """

    __import_error_message__ = "Unable to import the kernel. Please install `statOT` first as `pip install statot POT`."

    def __init__(
        self,
        adata: AnnData,
        source_idx: Union[str, np.ndarray],
        sink_idx: Union[str, np.ndarray],
        g: Union[str, np.ndarray],
        **kwargs: Any,
    ):
        if isinstance(source_idx, str):
            source_idx = adata.obs[source_idx].values
        if isinstance(sink_idx, str):
            sink_idx = adata.obs[sink_idx].values
        source_idx, sink_idx = np.asarray(source_idx), np.asarray(sink_idx)
        if not np.issubdtype(source_idx.dtype, np.bool_):
            raise TypeError(
                f"Expected `source_idx` to be a boolean array, found `{source_idx.dtype}`."
            )
        if not np.issubdtype(sink_idx.dtype, np.bool_):
            raise TypeError(
                f"Expected `sink_idx` to be a boolean array, found `{sink_idx.dtype}`."
            )
        if np.any(source_idx & sink_idx):
            raise ValueError("Some cells are marked as both source and sink.")

        super().__init__(adata, source_idx=source_idx, sink_idx=sink_idx, g=g, **kwargs)

    def compute_transition_matrix(
        self,
        eps: float,
        dt: float,
        basis: str = "X_pca",
        cost_norm_method: Optional[str] = None,
        method: Literal["ent", "quad", "unbal"] = "ent",
        tol: float = 0.0,
        thresh: float = 0.0,
        maxiter: int = 5000,
        C: Optional[np.ndarray] = None,
        verbose: bool = False,
        **kwargs: Any,
    ) -> "OTKernel":
        """
        Compute transition matrix using stationary OT [Zhang21]_.

        Parameters
        ----------
        eps
            Regularization parameter.
        dt
            Choice of the time step over which to fit the model.
        basis
            Key in :attr:`anndata.AnnData.obsm` where the basis is stored.
        cost_norm_method
            Cost normalization method to use. Use "mean" to ensure `mean(C) = 1` or refer to
            :func:`ot.utils.cost_normalization` for more information.
        method
            Choice of regularization. Valid options are:

                - `"ent"` - entropy.
                - `"quad"` - L2-norm.
                - `"unbal"` - unbalanced transport (not yet implemented).

        tol
            Relative tolerance for OT solver convergence.
        thresh
            Threshold for output transition probabilities.
        maxiter
            Maximum number of iterations for OT solver.
        C
            Cost matrix for optimal transport problem.
        verbose
            Detailed output on convergence of OT solver.
        kwargs
            Additional keyword arguments.

        Returns
        -------
        :class:`cellrank.external.kernels.OTKernel`
            Makes :paramref:`transition_matrix` available.
        """
        if method not in ("ent", "quad", "unbal"):
            raise ValueError(f"Invalid method `{method!r}`.")
        if method == "unbal":
            raise NotImplementedError("Method `'unbal'` is not yet implemented.")

        super().compute_transition_matrix(
            eps,
            dt,
            expr_key=basis,
            cost_norm_method=cost_norm_method,
            method=method,
            tol=tol,
            thresh=thresh,
            maxiter=maxiter,
            C=C,
            verbose=verbose,
            **kwargs,
        )

        return self