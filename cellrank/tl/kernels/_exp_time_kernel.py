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
    Base class which computes directed transition probabilities based on experimental time.

    %(density_correction)s

    Parameters
    ----------
    %(adata)s
    %(backward)s
    time_key
        Key in :attr:`adata` ``.obs`` where experimental time is stored.
        The experimental time can be of either of a numeric or an ordered categorical type.
    %(cond_num)s
    """

    def __init__(
        self,
        adata: AnnData,
        backward: bool = False,
        time_key: str = "exp_time",
        compute_cond_num: bool = False,
    ):
        super().__init__(
            adata,
            backward=backward,
            time_key=time_key,
            compute_cond_num=compute_cond_num,
            check_connectivity=False,
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
            logg.warning("Time categories are not ordered. Using ascending order")
            exp_time.cat = exp_time.cat.as_ordered()

        self._exp_time = pd.Series(
            pd.Categorical(exp_time, ordered=True), index=self.adata.obs_names
        )
        if self.experimental_time.isnull().any():
            raise ValueError("Experimental time contains NaN values.")

    @property
    def experimental_time(self) -> pd.Series:
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
        if len(self.experimental_time.cat.categories) > 1:
            minn, maxx = self.experimental_time.min(), self.experimental_time.max()
            self._exp_time = pd.Series(
                pd.Categorical(
                    maxx - np.array(self.experimental_time) + minn, ordered=True
                ),
                index=self.experimental_time.index,
            )
        return self
