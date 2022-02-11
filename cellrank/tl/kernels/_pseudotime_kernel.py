from typing import Any, Union, Callable, Optional
from typing_extensions import Literal

from enum import auto

from anndata import AnnData
from cellrank import logging as logg
from cellrank.tl._enum import _DEFAULT_BACKEND, ModeEnum, Backend_t
from cellrank.ul._docs import d
from cellrank.tl._utils import _connected, _irreducible
from cellrank.tl.kernels._mixins import ConnectivityMixin
from cellrank.tl.kernels._base_kernel import BidirectionalKernel
from cellrank.tl.kernels.utils._pseudotime_scheme import (
    ThresholdSchemeABC,
    HardThresholdScheme,
    SoftThresholdScheme,
    CustomThresholdScheme,
)

import numpy as np

__all__ = ("PseudotimeKernel",)


class ThresholdScheme(ModeEnum):
    SOFT = auto()
    HARD = auto()


@d.dedent
class PseudotimeKernel(ConnectivityMixin, BidirectionalKernel):
    """
    Kernel which computes directed transition probabilities based on a KNN graph and pseudotime.

    The KNN graph contains information about the (undirected) connectivities among cells, reflecting their similarity.
    Pseudotime can be used to either remove edges that point against the direction of increasing pseudotime
    :cite:`setty:19` or to downweight them :cite:`stassen:21`.

    Parameters
    ----------
    %(adata)s
    %(backward)s
    time_key
        Key in :attr:`anndata.AnnData.obs` where the pseudotime is stored.
    kwargs
        Keyword arguments for :class:`cellrank.kernels.Kernel`.
    """

    def __init__(
        self,
        adata: AnnData,
        backward: bool = False,
        time_key: str = "dpt_pseudotime",
        **kwargs: Any,
    ):
        super().__init__(
            adata,
            backward=backward,
            time_key=time_key,
            **kwargs,
        )

    def _read_from_adata(self, time_key: str = "dpt_pseudotime", **kwargs: Any) -> None:
        super()._read_from_adata(**kwargs)
        # fmt: off
        self._time_key = time_key
        if time_key not in self.adata.obs:
            raise KeyError(f"Unable to find pseudotime in `adata.obs[{time_key!r}]`.")

        self._pseudotime = np.array(self.adata.obs[time_key]).astype(np.float64, copy=True)
        if np.any(np.isnan(self.pseudotime)):
            raise ValueError("Encountered NaN values in pseudotime.")
        # fmt: on

    @d.dedent
    def compute_transition_matrix(
        self,
        threshold_scheme: Union[
            Literal["soft", "hard"],
            Callable[[float, np.ndarray, np.ndarray], np.ndarray],
        ] = "hard",
        frac_to_keep: float = 0.3,
        b: float = 10.0,
        nu: float = 0.5,
        check_irreducibility: bool = False,
        n_jobs: Optional[int] = None,
        backend: Backend_t = _DEFAULT_BACKEND,
        show_progress_bar: bool = True,
        **kwargs: Any,
    ) -> "PseudotimeKernel":
        """
        Compute transition matrix based on KNN graph and pseudotemporal ordering.

        Depending on the choice of the ``threshold_scheme``, it is based on ideas by either *Palantir*
        :cite:`setty:19` or *VIA* :cite:`stassen:21`.

        Parameters
        ----------
        threshold_scheme
            Which method to use when biasing the graph. Valid options are:

                - `'hard'` - based on *Palantir* :cite:`setty:19` which removes some edges that point against
                  the direction of increasing pseudotime. To avoid disconnecting the graph, it does not
                  remove all edges that point against the direction of increasing pseudotime, but keeps the ones
                  that point to cells inside a close radius. This radius is chosen according to the local cell density.
                - `'soft'` - based on *VIA* :cite:`stassen:21` which downweights edges that points against the direction
                  of increasing pseudotime. Essentially, the further "behind" a query cell is in pseudotime with respect
                  to the current reference cell, the more penalized will be its graph-connectivity.
                - :class:`callable` - any function conforming to the signature of
                  :func:`cellrank.kernels.utils.ThresholdSchemeABC.__call__`.
        frac_to_keep
            Fraction of the closest neighbors (according to graph connectivities) are kept, no matter whether they lie
            in the pseudotemporal past or future. This is done to ensure that the graph remains connected.
            Only used when ``threshold_scheme = 'hard'``. Needs to fall within the interval `[0, 1]`.
        %(soft_scheme_kernel)s
        check_irreducibility
            Optional check for irreducibility of the final transition matrix.
        %(parallel)s
        kwargs
            Keyword arguments for ``threshold_scheme``.

        Returns
        -------
        Self and updates :attr:`transition_matrix` and :attr:`params`.
        """
        if self.pseudotime is None:
            raise ValueError("Compute pseudotime first.")  # CytoTraceKernel

        start = logg.info("Computing transition matrix based on pseudotime`")
        if isinstance(threshold_scheme, str):
            threshold_scheme = ThresholdScheme(threshold_scheme)
            if threshold_scheme == ThresholdScheme.SOFT:
                scheme = SoftThresholdScheme()
                kwargs["b"] = b
                kwargs["nu"] = nu
            elif threshold_scheme == ThresholdScheme.HARD:
                scheme = HardThresholdScheme()
                kwargs["frac_to_keep"] = frac_to_keep
            else:
                raise NotImplementedError(
                    f"Threshold scheme `{threshold_scheme}` is not yet implemented."
                )
        elif isinstance(threshold_scheme, ThresholdSchemeABC):
            scheme = threshold_scheme
        elif callable(threshold_scheme):
            scheme = CustomThresholdScheme(threshold_scheme)
        else:
            raise TypeError(
                f"Expected `threshold_scheme` to be either a `str` or a `callable`, found `{type(threshold_scheme)}`."
            )

        # fmt: off
        if self._reuse_cache({"dnorm": False, "scheme": str(threshold_scheme), **kwargs}, time=start):
            return self
        # fmt: on

        biased_conn = scheme.bias_knn(
            self._conn,
            self.pseudotime,
            n_jobs=n_jobs,
            backend=backend,
            show_progress_bar=show_progress_bar,
            **kwargs,
        )

        # make sure the biased graph is still connected
        if not _connected(biased_conn):
            logg.warning("Biased KNN graph is disconnected")
        if check_irreducibility and not _irreducible(biased_conn):
            logg.warning("Biased KNN graph is not irreducible")

        self.transition_matrix = biased_conn
        logg.info("    Finish", time=start)

        return self

    @property
    def pseudotime(self) -> Optional[np.array]:
        """Pseudotemporal ordering of cells."""
        return self._pseudotime

    def __invert__(self) -> "PseudotimeKernel":
        pk = self._copy_ignore("_transition_matrix")
        pk._pseudotime = np.max(pk.pseudotime) - pk.pseudotime
        pk._backward = not self.backward
        pk._params = {}
        return pk
