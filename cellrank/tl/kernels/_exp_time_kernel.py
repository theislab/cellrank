from typing import Any, Dict, Tuple, Union, Mapping, Optional
from typing_extensions import Literal

from abc import ABC, abstractmethod
from copy import copy
from types import MappingProxyType

import scanpy as sc
from anndata import AnnData
from cellrank import logging as logg
from cellrank.tl._enum import ModeEnum
from cellrank.ul._docs import d, inject_docs
from cellrank.tl._utils import _normalize
from cellrank.tl.kernels import Kernel
from cellrank.tl.kernels._utils import _ensure_numeric_ordered
from cellrank.tl.kernels._base_kernel import KernelExpression

import numpy as np
import pandas as pd
from scipy.sparse import bmat, spdiags

from matplotlib.colors import Normalize, to_hex
from matplotlib.pyplot import get_cmap


class LastTimePoint(ModeEnum):  # noqa
    UNIFORM = "uniform"
    DIAGONAL = "diagonal"
    CONNECTIVITIES = "connectivities"


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


class TransportMapKernel(ExperimentalTimeKernel, ABC):
    """Kernel base class which computes transition matrix based on transport maps for consecutive time pairs."""

    def __init__(self, *args: Any, **kwargs: Any):
        super().__init__(*args, **kwargs)
        self._tmaps: Optional[Dict[Tuple[Any, Any], AnnData]] = None

    def compute_transition_matrix(
        self,
        last_time_point: LastTimePoint = LastTimePoint.DIAGONAL,
        conn_kwargs: Mapping[str, Any] = MappingProxyType({}),
        **kwargs: Any,
    ) -> KernelExpression:
        """TODO - ltp/conn kwargs docrep."""
        timepoints = self.experimental_time.cat.categories
        timepoints = list(zip(timepoints[:-1], timepoints[1:]))

        tmap = self._restich_tmaps(
            {(t1, t2): self._compute_tmap(t1, t2, **kwargs) for t1, t2 in timepoints},
            normalize=True,
            last_time_point=last_time_point,
            conn_kwargs=conn_kwargs,
        )
        self._compute_transition_matrix(
            matrix=tmap.X,
            density_normalize=False,
            check_irreducibility=False,
        )

        return self

    @abstractmethod
    def _compute_tmap(self, t1: Any, t2: Any, **kwargs: Any) -> AnnData:
        """
        TODO.

        Parameters
        ----------
        t1
        t2
        kwargs

        Returns
        -------
        TODO.
        """

    @inject_docs(ltp=LastTimePoint)
    def _restich_tmaps(
        self,
        tmaps: Mapping[Tuple[Any, Any], AnnData],
        last_time_point: LastTimePoint = LastTimePoint.DIAGONAL,
        conn_kwargs: Mapping[str, Any] = MappingProxyType({}),
        normalize: bool = True,
    ) -> AnnData:
        """
        TODO - intro + ltp/conn kwargs docrep.

        Parameters
        ----------
        tmaps
            Sorted transport maps as ``{{(t1, t2): tmat_2, (t2, t3): tmat_2, ...}}``.
            :attr:`anndata.AnnData.var_names` at the current time-point must be the same as
            :attr:`anndata.AnnData.obs_names` at the next time-point.
        last_time_point
            How to define transitions within the last time point. Valid options are:

                - `{ltp.UNIFORM!r}` - row-normalized matrix of 1s for transitions within the last time point.
                - `{ltp.DIAGONAL!r}` - diagonal matrix with 1s on the diagonal.
                - `{ltp.CONNECTIVITIES!r}` - use transitions from :class:`cellrank.tl.kernels.ConnectivityKernel`
                  derived from the last time point subset of :attr:`adata`.
        conn_kwargs
            Keyword arguments for :func:`scanpy.pp.neighbors`, when using ``last_time_point = {ltp.CONNECTIVITIES!r}``.
            Can contain `'density_normalize'` for
            :meth:`cellrank.tl.kernels.ConnectivityKernel.compute_transition_matrix`.
        normalize
            Whether to normalize the transition matrix so that rows sum to `1`.

        Returns
        -------
        Concatenated transport maps into an :class:`anndata.AnnData` object.
        """
        from cellrank.tl.kernels import ConnectivityKernel

        tmaps = self._validate_tmaps(tmaps)

        conn_kwargs = dict(conn_kwargs)
        conn_kwargs["copy"] = False
        _ = conn_kwargs.pop("key_added", None)
        density_normalize = conn_kwargs.pop("density_normalize", True)

        blocks = [[None] * (len(tmaps) + 1) for _ in range(len(tmaps) + 1)]
        nrows, ncols = 0, 0
        obs_names, obs = [], []

        # TODO(michalk8): what happens in bwd case or if not sorted?
        for i, tmap in enumerate(tmaps.values()):
            blocks[i][i + 1] = _normalize(tmap.X) if normalize else tmap.X
            nrows += tmap.n_obs
            ncols += tmap.n_vars
            obs_names.extend(tmap.obs_names)
            obs.append(tmap.obs)
        obs_names.extend(tmap.var_names)

        n = self.adata.n_obs - nrows
        if last_time_point == LastTimePoint.DIAGONAL:
            blocks[-1][-1] = spdiags([1] * n, 0, n, n)
        elif last_time_point == LastTimePoint.UNIFORM:
            blocks[-1][-1] = np.ones((n, n)) / float(n)
        elif last_time_point == LastTimePoint.CONNECTIVITIES:
            adata_subset = self.adata[tmap.var_names].copy()
            sc.pp.neighbors(adata_subset, **conn_kwargs)
            blocks[-1][-1] = (
                ConnectivityKernel(adata_subset)
                .compute_transition_matrix(density_normalize)
                .transition_matrix
            )
        else:
            raise NotImplementedError(
                f"Last time point mode `{last_time_point}` is not yet implemented."
            )

        # prevent the last block from disappearing
        n = blocks[0][1].shape[0]
        blocks[0][0] = spdiags([], 0, n, n)

        tmp = AnnData(bmat(blocks, format="csr"))
        tmp.obs_names = obs_names
        tmp.var_names = obs_names
        tmp = tmp[self.adata.obs_names, :][:, self.adata.obs_names]

        tmp.obs = pd.merge(
            tmp.obs,
            pd.concat(obs),
            left_index=True,
            right_index=True,
            how="left",
        )

        return tmp

    def _threshold_transition_matrix(
        self, threshold: Union[float, Literal["auto"]]
    ) -> None:
        """
        Remove small non-zero values from :attr:`transition_matrix`.

        Parameters
        ----------
        threshold
            How to remove small non-zero values from the transition matrix. Valid options are:

                - `'auto'` - find the maximum threshold value which will not remove every non-zero value from any row.
                - :class:`float` - value in `[0, 100]` corresponding to a percentage of non-zeros to remove.
                  Rows where all values are removed will have uniform distribution.
                - `None` - do not threshold.

        Returns
        -------
        Nothing, just updates :attr:`transition_matrix`.
        """
        tmat = self.transition_matrix
        if threshold == "auto":
            threshold = min(np.max(tmat[i].data) for i in range(tmat.shape[0]))
            logg.info(f"Using `threshold={threshold}`")
            tmat.data[tmat.data < threshold] = 0.0
        else:
            if not (0 <= threshold <= 100):
                raise ValueError(
                    f"Expected `threshold to be in `[0, 100]`, found `{threshold}`.`"
                )
            threshold = np.percentile(tmat.data, threshold)
            logg.info(f"Using `threshold={threshold}`")
            tmat.data[tmat.data <= threshold] = 0.0

        tmat.eliminate_zeros()

        self._compute_transition_matrix(
            matrix=tmat,
            density_normalize=False,
            check_irreducibility=False,
        )

    def _validate_tmaps(
        self, tmaps: Mapping[Tuple[Any, Any], AnnData]
    ) -> Mapping[Tuple[Any, Any], AnnData]:
        """TODO."""
        return tmaps

    @property
    def transport_maps(self) -> Optional[Dict[Tuple[Any, Any], AnnData]]:
        """Transport maps for consecutive time pairs."""
        return self._tmaps
