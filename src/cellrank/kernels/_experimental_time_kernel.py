import abc
from typing import Any, Optional

import numpy as np
import pandas as pd

from matplotlib.colors import Normalize, to_hex
from matplotlib.pyplot import get_cmap

from anndata import AnnData

from cellrank._utils._docs import d
from cellrank.kernels._base_kernel import BidirectionalKernel
from cellrank.kernels._utils import _ensure_numeric_ordered, require_tmat

__all__ = ["ExperimentalTimeKernel"]


@d.dedent
class ExperimentalTimeKernel(BidirectionalKernel, abc.ABC):
    """Kernel base class which computes transition probabilities based on experimental time.

    %(density_correction)s

    Parameters
    ----------
    %(adata)s
    %(backward)s
    time_key
        Key in :attr:`~anndata.AnnData.obs` where experimental time is stored.
        The experimental time can be of either of a numeric or an ordered categorical type.
    kwargs
        Keyword arguments for the parent class.
    """

    def __init__(
        self,
        adata: AnnData,
        time_key: str,
        backward: bool = False,
        **kwargs: Any,
    ):
        super().__init__(
            adata,
            backward=backward,
            time_key=time_key,
            **kwargs,
        )

    def _read_from_adata(self, time_key: str, cmap: str = "gnuplot", **kwargs: Any) -> None:
        super()._read_from_adata(**kwargs)

        self._time_key = time_key
        self._exp_time = _ensure_numeric_ordered(self.adata, time_key)
        self.adata.obs[time_key] = self.experimental_time.values

        # fmt: off
        cmap = get_cmap(cmap)
        cats = self.experimental_time.cat.categories
        norm = Normalize(vmin=cats.min(), vmax=cats.max())
        self.adata.uns[f"{time_key}_colors"] = np.array([to_hex(c) for c in cmap(norm(cats))])
        # fmt: on

    @d.dedent
    @require_tmat
    def plot_single_flow(
        self,
        cluster: str,
        cluster_key: str,
        time_key: Optional[str] = None,
        *args: Any,
        **kwargs: Any,
    ) -> None:
        """%(plot_single_flow.full_desc)s

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

    def __invert__(self, *args: str) -> "ExperimentalTimeKernel":
        etk = self._copy_ignore(*(("_transition_matrix",) + args))
        minn, maxx = etk.experimental_time.min(), etk.experimental_time.max()
        etk._exp_time = pd.Series(
            pd.Categorical(maxx - np.array(etk.experimental_time) + minn, ordered=True),
            index=etk.experimental_time.index,
        )
        etk._backward = not self.backward
        etk._params = {}
        return etk
