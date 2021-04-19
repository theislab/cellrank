"""Experimental time kernel module."""
from abc import ABC
from copy import copy
from typing import Any

import numpy as np
import pandas as pd
from pandas.api.types import infer_dtype
from pandas.core.dtypes.common import (
    is_object_dtype,
    is_numeric_dtype,
    is_categorical_dtype,
)

from cellrank import logging as logg
from cellrank.ul._docs import d
from cellrank.tl.kernels import Kernel
from cellrank.tl.kernels._base_kernel import AnnData


@d.dedent
class ExperimentalTimeKernel(Kernel, ABC):
    """
    Base kernel class which computes directed transition probabilities based on experimental time.

    %(density_correction)s

    Parameters
    ----------
    %(adata)s
    %(backward)s
    time_key
        Key in :paramref:`adata` ``.obs`` where the experimental time is stored.
        The experimental time can be of either a numeric or an ordered categorical type.
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

    def _read_from_adata(self, **kwargs: Any) -> None:
        super()._read_from_adata(**kwargs)

        time_key = kwargs.pop("time_key", "exp_time")
        if time_key not in self.adata.obs.keys():
            raise KeyError(f"Could not find time key in `adata.obs[{time_key!r}]`.")

        exp_time = self.adata.obs[time_key].copy()
        if not is_categorical_dtype(exp_time):
            exp_time = np.array(exp_time)
            if is_object_dtype(exp_time):
                try:
                    exp_time = exp_time.astype(float)
                except ValueError as e:
                    raise RuntimeError(
                        f"Unable to convert `adata.obs[{time_key!r}]` to `float` dtype."
                    ) from e
            if not is_numeric_dtype(exp_time):
                raise TypeError(
                    f"Expected experimental time to be `numeric` or `categorical`, found `{infer_dtype(exp_time)}`."
                )
            exp_time = pd.Series(
                pd.Categorical(
                    exp_time,
                    categories=sorted(set(exp_time[~np.isnan(exp_time)])),
                    ordered=True,
                )
            )

        if not exp_time.cat.ordered:
            logg.warning("Ordering categories")
            exp_time.cat = exp_time.cat.as_ordered()

        self._exp_time = pd.Categorical(exp_time, ordered=True)

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
        if len(self.experimental_time.categories) > 1:
            minn, maxx = self.experimental_time.min(), self.experimental_time.max()
            self._exp_time = pd.Categorical(
                maxx - np.array(self.experimental_time) + minn, ordered=True
            )
        return self
