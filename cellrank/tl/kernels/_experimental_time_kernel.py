from typing import Any, Optional

from abc import ABC
from copy import copy

from anndata import AnnData
from cellrank.ul._docs import d
from cellrank.tl.kernels import Kernel
from cellrank.tl.kernels._utils import _ensure_numeric_ordered

import numpy as np
import pandas as pd

from matplotlib.colors import Normalize, to_hex
from matplotlib.pyplot import get_cmap

__all__ = ("ExperimentalTimeKernel",)


@d.dedent
class ExperimentalTimeKernel(Kernel, ABC):
    """
    Kernel base class which computes directed transition probabilities based on experimental time.

    %(density_correction)s

    Parameters
    ----------
    %(adata)s
    %(backward)s
    time_key
        Key in :attr:`anndata.AnnData.obs` where experimental time is stored.
        The experimental time can be of either of a numeric or an ordered categorical type.
    %(cond_num)s
    """

    def __init__(
        self,
        adata: AnnData,
        backward: bool = False,
        time_key: str = "exp_time",
        compute_cond_num: bool = False,
        **kwargs: Any,
    ):
        super().__init__(
            adata,
            backward=backward,
            time_key=time_key,
            compute_cond_num=compute_cond_num,
            check_connectivity=False,
            **kwargs,
        )
        self._time_key = time_key

    def _read_from_adata(self, **kwargs: Any) -> None:
        super()._read_from_adata(**kwargs)

        time_key = kwargs.pop("time_key", "exp_time")
        self._exp_time = _ensure_numeric_ordered(self.adata, time_key)
        self.adata.obs[time_key] = self.experimental_time.values

        # fmt: off
        cmap = get_cmap(kwargs.pop("cmap", "gnuplot"))
        cats = self.experimental_time.cat.categories
        norm = Normalize(vmin=cats.min(), vmax=cats.max())
        self.adata.uns[f"{time_key}_colors"] = np.array([to_hex(c) for c in cmap(norm(cats) * cmap.N)])
        # fmt: on

    @d.dedent
    def plot_single_flow(
        self,
        cluster: str,
        cluster_key: str,
        time_key: Optional[str] = None,
        *args: Any,
        **kwargs: Any,
    ) -> None:
        """
        %(plot_single_flow.full_desc)s

        Parameters
        ----------
        %(plot_single_flow.parameters)s

        Returns
        -------
        %(plot_single_flow.returns)s
        """  # noqa: D400
        if time_key is None:
            time_key = self._time_key
        return super().plot_single_flow(cluster, cluster_key, time_key, *args, **kwargs)

    @property
    def experimental_time(self) -> pd.Series:
        """Experimental time."""
        return self._exp_time

    def copy(self) -> "ExperimentalTimeKernel":
        """Return a copy of self."""
        pk = type(self)(self.adata, backward=self.backward, time_key=self._time_key)
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
