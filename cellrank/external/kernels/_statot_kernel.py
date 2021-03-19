from typing import Any, Union, Optional

from typing_extensions import Literal

from anndata import AnnData

import numpy as np
import pandas as pd

from cellrank.ul._docs import d
from cellrank.tl.estimators import GPCCA
from cellrank.tl.kernels._precomputed_kernel import DummyKernel

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
    terminal_states
        Categorical :class:`pandas.Series` where non-`NaN` values mark terminal states.
        If `None`, terminal states are assumed to be already present in :paramref:`adata` ``['terminal_states']``.
    g
        Key in :attr:`anndata.AnnData.obs` containing relative growth rates for cells or the array itself.
    cluster_key
        If a key to cluster labels is given, `terminal_states` will be associated with these for naming and colors.
    kwargs
        Additional keyword arguments.
    """

    __import_error_message__ = "Unable to import the kernel. Please install `statOT` first as `pip install statot POT`."

    def __init__(
        self,
        adata: AnnData,
        g: Union[str, np.ndarray],
        terminal_states: Optional[Union[str, pd.Series]] = None,
        cluster_key: Optional[str] = None,
        **kwargs: Any,
    ):
        if terminal_states is not None:
            dk = DummyKernel(adata, backward=False)
            estim = GPCCA(dk, write_to_adata=True)
            estim.set_terminal_states(terminal_states, cluster_key=cluster_key)

        try:
            super().__init__(adata, g=g, **kwargs)
        except Exception as e:  # noqa: B902
            raise RuntimeError("Unable to initialize the kernel.") from e

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

                - `'ent'`: entropy.
                - `'quad'`: L2-norm.
                - `'unbal'`: unbalanced transport (not yet implemented).

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
