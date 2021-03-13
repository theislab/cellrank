"""Experimental kernel module."""
from copy import copy
from typing import Any

import pandas as pd
from pandas.core.dtypes.common import is_integer_dtype, is_categorical_dtype

from cellrank import logging as logg
from cellrank.ul._docs import d
from cellrank.tl._utils import _connected
from cellrank.tl.kernels import Kernel
from cellrank.tl._constants import Direction
from cellrank.tl.kernels._base_kernel import AnnData, _dtype
from cellrank.tl.kernels._pseudotime_schemes import HardThresholdScheme


class ExperimentalTimeKernel(Kernel):
    """
    Kernel which computes directed transition probabilities based on a KNN graph and experimental time.

    %(density_correction)s

    Parameters
    ----------
    %(adata)s
    %(backward)s
    time_key
        Key in :paramref:`adata` ``.obs`` where the experimental time is stored. TODO.
    compute_cond_num
        Whether to compute condition number of the transition matrix. Note that this might be costly,
        since it does not use sparse implementation.
    """

    def __init__(
        self,
        adata: AnnData,
        backward: bool = False,
        time_key: str = "exp_time",
        compute_cond_num: bool = False,
        check_connectivity: bool = False,
    ):
        super().__init__(
            adata,
            backward=backward,
            time_key=time_key,
            compute_cond_num=compute_cond_num,
            check_connectivity=check_connectivity,
        )
        self._time_key = time_key

    def _read_from_adata(self, **kwargs: Any):
        super()._read_from_adata(**kwargs)

        time_key = kwargs.pop("time_key", "dpt_pseudotime")
        if time_key not in self.adata.obs.keys():
            raise KeyError(f"Could not find time key in `adata.obs[{time_key!r}]`.")

        exp_time = self.adata.obs[time_key].copy()
        if is_integer_dtype(exp_time):
            exp_time = pd.Categorical(
                exp_time, categories=sorted(set(exp_time))[::-1], ordered=True
            )
        if not is_categorical_dtype(exp_time):
            raise TypeError()

        if not exp_time.ordered:
            exp_time = exp_time.as_ordered()

        self._exp_time = pd.Categorical(exp_time, ordered=True)

    @d.dedent
    def compute_transition_matrix(
        self, density_normalize: bool = True, check_irreducibility: bool = False
    ) -> "ExperimentalTimeKernel":
        """
        TODO.

        Parameters
        ----------
        %(dnorm_irred)s

        Returns
        -------
        :class:`cellrank.tl.kernels.ExperimentalTime`
        """
        start = logg.info("Computing transition matrix based on experimental time")

        if self._reuse_cache({"dnorm": density_normalize}, time=start):
            return self

        # handle backward case and run biasing function
        exp_time = (
            self.experimental_time.reorder_categories(
                self.experimental_time.categories[::-1], ordered=True
            )
            if self._direction == Direction.BACKWARD
            else self.experimental_time
        )

        biased_conn = HardThresholdScheme().bias_knn(
            self._conn.copy(), exp_time, k=1, n_neighs=1
        )
        biased_conn = biased_conn.astype(_dtype)

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
    def experimental_time(self) -> pd.Categorical:
        """Experimental time."""
        return self._exp_time

    def copy(self) -> "ExperimentalTimeKernel":
        """Return a copy of self."""
        pk = ExperimentalTimeKernel(
            self.adata, backward=self.backward, time_key=self._time_key
        )
        pk._exp_time = copy(self.experimental_time)
        pk._params = copy(self._params)
        pk._cond_num = self.condition_number
        pk._transition_matrix = copy(self._transition_matrix)

        return pk

    def __invert__(self) -> "ExperimentalTimeKernel":
        super().__invert__()
        self._exp_time = self.experimental_time.reorder_categories(
            self.experimental_time.categories[::-1], ordered=True
        )
        return self
