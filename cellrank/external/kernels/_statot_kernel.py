from typing import Any, Union, Optional

from typing_extensions import Literal

from anndata import AnnData

import numpy as np

from cellrank.ul._docs import d

try:
    from statot import OTKernel as OTKernel_

    print("x")
except ImportError:
    from cellrank.external.kernels._import_error_kernel import (
        ErroredKernel as OTKernel_,
    )


@d.dedent
class OTKernel(OTKernel_):
    """
    Stationary optimal transport kernel from [Zhang21]_.

    This class requires the `statot` package, which can be installed as `pip install statot`.

    Parameters
    ----------
    %(adata)s
    source_idx
        Key in :attr:`anndata.AnnData.obs` containing boolean mask marking the source indices or the mask itself.
    sink_idx
        Key in :attr:`anndata.AnnData.obs` containing boolean marking the sink indices or the mask itself.
    g
        Key in :attr:`anndata.AnnData.obs` containing the relative growth rates for cells or the array itself.
    kwargs
        Additional keyword arguments.
    """

    __import_error_message__ = "Unable to import the kernel. Please install statOT first as `pip install statot`."

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
        super().__init__(adata, source_idx=source_idx, sink_idx=sink_idx, g=g, **kwargs)

    def compute_transition_matrix(
        self,
        eps: float,
        dt: float,
        basis: str = "X_pca",
        cost_norm_method: Optional[str] = None,
        method: Literal["ent", "quad", "unbal", "marginals"] = "ent",
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
            Key in :attr:`anndata.AnnData.obsm`.
        cost_norm_method
            Cost normalization method to use. Use "mean" to ensure `mean(C) = 1` or refer to
            :func:`ot.utils.cost_normalization` for more information.
        method
            Choice of regularization. Valid options are:

                - `"ent"` - entropy.
                - `"quad"` - L2-norm.
                - `"unbal"` - unbalanced transport (not yet implemented).
                - `"marginals"` - marginals (returns just `mu` and `nu`).

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
        return super().compute_transition_matrix(
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
