"""Pseudotime kernel module."""
from copy import copy
from typing import Any, Union, Callable, Optional

from typing_extensions import Literal

from anndata import AnnData

import numpy as np

from cellrank import logging as logg
from cellrank.ul._docs import d
from cellrank.tl._utils import _connected
from cellrank.tl.kernels import Kernel
from cellrank.tl._constants import Direction, ThresholdScheme
from cellrank.tl.kernels._base_kernel import _dtype
from cellrank.tl.kernels._pseudotime_schemes import (
    ThresholdSchemeABC,
    HardThresholdScheme,
    SoftThresholdScheme,
    CustomThresholdScheme,
)


@d.dedent
class PseudotimeKernel(Kernel):
    """
    Kernel which computes directed transition probabilities based on a KNN graph and pseudotime.

    The KNN graph contains information about the (undirected) connectivities among cells, reflecting their similarity.
    Pseudotime can be used to either remove edges that point against the direction of increasing pseudotime (see
    [Setty19]_, or to downweight them (see [VIA21]_).

    %(density_correction)s

    Parameters
    ----------
    %(adata)s
    %(backward)s
    time_key
        Key in :paramref:`adata` ``.obs`` where the pseudotime is stored.
    compute_cond_num
        Whether to compute condition number of the transition matrix. Note that this might be costly,
        since it does not use sparse implementation.
    """

    def __init__(
        self,
        adata: AnnData,
        backward: bool = False,
        time_key: str = "dpt_pseudotime",
        compute_cond_num: bool = False,
        check_connectivity: bool = False,
        **kwargs: Any,
    ):
        super().__init__(
            adata,
            backward=backward,
            time_key=time_key,
            compute_cond_num=compute_cond_num,
            check_connectivity=check_connectivity,
            **kwargs,
        )
        self._time_key = time_key

    def _read_from_adata(self, time_key: str, **kwargs):
        super()._read_from_adata(**kwargs)

        if time_key not in self.adata.obs.keys():
            raise KeyError(f"Could not find time key in `adata.obs[{time_key!r}]`.")

        self._pseudotime = np.array(self.adata.obs[time_key]).astype(_dtype)

        if np.any(np.isnan(self._pseudotime)):
            raise ValueError("Encountered NaN values in pseudotime.")

    @d.dedent
    def compute_transition_matrix(
        self,
        threshold_scheme: Union[Literal["soft", "hard"], Callable] = "soft",
        frac_to_keep: float = 0.3,
        b: float = 20.0,
        nu: float = 1.0,
        percentile: Optional[int] = 95,
        density_normalize: bool = True,
        check_irreducibility: bool = False,
        **kwargs: Any,
    ) -> "PseudotimeKernel":
        """
        Compute transition matrix based on KNN graph and pseudotemporal ordering.

        Depending on the choice of the `thresholding_scheme`, this is based on ideas by either Palantir (see [Setty19]_)
        or VIA (see [VIA21]_).

        When using a `'hard'` thresholding scheme, this based on ideas by *Palantir* (see [Setty19]_) which removes some
        edges that point against the direction of increasing pseudotime. To avoid disconnecting the graph, it does not
        remove all edges that point against the direction of increasing pseudotime but keeps the ones that point to
        cells inside a close radius. This radius is chosen according to the local cell density.

        When using a `'soft'` thresholding scheme, this is based on ideas by *VIA* (see [VIA21]_) which downweights
        edges that points against the direction of increasing pseudotime. Essentially, the further "behind" a query
        cell is in pseudotime with respect to the current reference cell, the more penalized will be its
        graph-connectivity.

        Parameters
        ----------
        frac_to_keep
            The `fract_to_keep` * n_neighbors closest neighbors (according to graph connectivities) are kept, no matter
            whether they lie in the pseudotemporal past or future. This is done to ensure that the graph remains
            connected. Only used when `threshold_scheme='hard'`.
        %(soft_scheme_kernel)s
        density_normalize
            Whether or not to use the underlying KNN graph for density normalization.
        check_irreducibility
            Optional check for irreducibility of the final transition matrix.

        Returns
        -------
        :class:`cellrank.tl.kernels.PseudotimeKernel`
            Makes :paramref:`transition_matrix` available.
        """
        start = logg.info(f"Computing transition matrix based on `{self._time_key}`")

        # get the connectivities and number of neighbors
        n_neighbors = (
            self.adata.uns.get("neighbors", {})
            .get("params", {})
            .get("n_neighbors", None)
        )
        if n_neighbors is None:
            logg.warning(
                "Could not find 'n_neighbors' in `adata.uns['neighbors']['params']`. Using an estimate"
            )
            n_neighbors = np.min(self._conn.sum(1))

        if isinstance(threshold_scheme, str):
            threshold_scheme = ThresholdScheme(threshold_scheme)
            if threshold_scheme == ThresholdScheme.SOFT:
                scheme = SoftThresholdScheme()
                kwargs["b"], kwargs["nu"] = b, nu
            elif threshold_scheme == ThresholdScheme.HARD:
                scheme = HardThresholdScheme()
                kwargs["frac_to_keep"], kwargs["n_neighs"] = frac_to_keep, n_neighbors
            else:
                raise NotImplementedError(
                    f"Threshold scheme `{threshold_scheme}` is not yet implemented."
                )
        elif isinstance(threshold_scheme, ThresholdSchemeABC):
            pass
        elif callable(threshold_scheme):
            scheme = CustomThresholdScheme(threshold_scheme)
        else:
            raise TypeError(
                f"Expected `threshold_scheme` to be either a `str` or a `callable`, found `{type(threshold_scheme)}`."
            )

        # fmt: off
        if self._reuse_cache({"dnorm": density_normalize, "scheme": str(threshold_scheme), **kwargs}, time=start):
            return self

        # handle backward case and run biasing function
        pseudotime = (
            np.max(self.pseudotime) - self.pseudotime
            if self._direction == Direction.BACKWARD
            else self.pseudotime
        )

        biased_conn = scheme.bias_knn(self._conn, pseudotime, **kwargs).astype(_dtype)
        # fmt: on

        # make sure the biased graph is still connected
        if not _connected(biased_conn):
            logg.warning("Biased KNN graph is disconnected")

        self._compute_transition_matrix(
            matrix=biased_conn,
            density_normalize=density_normalize,
            check_irreducibility=check_irreducibility,
        )
        logg.info("    Finish", time=start)

        return self

    @property
    def pseudotime(self) -> np.array:
        """Pseudotemporal ordering of cells."""
        return self._pseudotime

    def copy(self) -> "PseudotimeKernel":
        """Return a copy of self."""
        pk = PseudotimeKernel(
            self.adata, backward=self.backward, time_key=self._time_key
        )
        pk._pseudotime = copy(self.pseudotime)
        pk._params = copy(self._params)
        pk._cond_num = self.condition_number
        pk._transition_matrix = copy(self._transition_matrix)

        return pk

    def __invert__(self) -> "PseudotimeKernel":
        super().__invert__()
        self._pseudotime = np.max(self.pseudotime) - self.pseudotime
        return self
