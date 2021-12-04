from typing import Any, Dict, Tuple, Union, Mapping, Optional
from typing_extensions import Literal

from abc import ABC, abstractmethod
from types import MappingProxyType
from contextlib import contextmanager

import scanpy as sc
from anndata import AnnData
from cellrank import logging as logg
from cellrank.tl._enum import ModeEnum
from cellrank.ul._docs import d, inject_docs
from cellrank.tl.kernels._base_kernel import KernelExpression
from cellrank.tl.kernels._experimental_time_kernel import ExperimentalTimeKernel

import numpy as np
import pandas as pd
from scipy.sparse import bmat, spdiags, spmatrix

__all__ = ("TransportMapKernel",)


class LastTimePoint(ModeEnum):
    UNIFORM = "uniform"
    DIAGONAL = "diagonal"
    CONNECTIVITIES = "connectivities"


class TransportMapKernel(ExperimentalTimeKernel, ABC):
    """Kernel base class which computes transition matrix based on transport maps for consecutive time pairs."""

    def __init__(self, *args: Any, **kwargs: Any):
        super().__init__(*args, **kwargs)
        self._tmaps: Optional[Dict[Tuple[Any, Any], AnnData]] = None

    @abstractmethod
    def _compute_tmap(
        self, t1: Any, t2: Any, **kwargs: Any
    ) -> Union[np.ndarray, spmatrix, AnnData]:
        """
        Compute transport matrix for a time point pair.

        Parameters
        ----------
        t1
            Earlier time point in :attr:`experimental_time`.
        t2
            Later time point in :attr:`experimental_time`.
        kwargs
            Additional keyword arguments.

        Returns
        -------
        Array of shape ``(n_cells_early, n_cells_late)``. If :class:`anndata.AnnData`, :attr:`anndata.AnnData.obs_names`
        and :attr:`anndata.AnnData.var_names` correspond to subsets of observations from :attr:`adata`.
        """

    @d.get_sections(base="tmk_tmat", sections=["Parameters"])
    @d.dedent
    def compute_transition_matrix(
        self,
        threshold: Optional[Union[float, Literal["auto"]]] = "auto",
        last_time_point: Literal[
            "uniform", "diagonal", "connectivities"
        ] = LastTimePoint.DIAGONAL,
        conn_kwargs: Mapping[str, Any] = MappingProxyType({}),
        **kwargs: Any,
    ) -> KernelExpression:
        """
        Compute transition matrix using transport maps.

        Parameters
        ----------
        %(tmk_thresh.parameters)s
        last_time_point
            How to define transitions within the last time point. Valid options are:

                - `{ltp.UNIFORM!r}` - row-normalized matrix of 1s for transitions within the last time point.
                - `{ltp.DIAGONAL!r}` - diagonal matrix with 1s on the diagonal.
                - `{ltp.CONNECTIVITIES!r}` - use transitions from :class:`cellrank.tl.kernels.ConnectivityKernel`
                  derived from the last time point subset of :attr:`adata`.
        conn_kwargs
            Keyword arguments for :func:`scanpy.pp.neighbors` when using ``last_time_point = {ltp.CONNECTIVITIES!r}``.
            Can have `'density_normalize'` for :meth:`cellrank.kernels.ConnectivityKernel.compute_transition_matrix`.

        Returns
        -------
        Self and updates :attr:`transition_matrix` and :attr:`params`.
        """
        ckwargs = dict(kwargs)
        ckwargs["threshold"] = threshold
        ckwargs["last_time_point"] = last_time_point
        if self._reuse_cache(ckwargs):
            return self

        timepoints = self.experimental_time.cat.categories
        timepoints = list(zip(timepoints[:-1], timepoints[1:]))

        # fmt: off
        tmap = self._restich_tmaps(
            {(t1, t2): self._tmat_to_adata(t1, t2, self._compute_tmap(t1, t2, **kwargs)) for t1, t2 in timepoints},
            last_time_point=last_time_point,
            conn_kwargs=conn_kwargs,
        )
        # fmt: on
        self.transition_matrix = tmap.X
        if threshold:
            self._threshold_transition_matrix(threshold)

        return self

    @d.dedent
    @inject_docs(ltp=LastTimePoint)
    def _restich_tmaps(
        self,
        tmaps: Mapping[Tuple[Any, Any], AnnData],
        last_time_point: LastTimePoint = LastTimePoint.DIAGONAL,
        conn_kwargs: Mapping[str, Any] = MappingProxyType({}),
    ) -> AnnData:
        """
        Group individual transport maps into 1 matrix aligned with :attr:`adata`.

        Parameters
        ----------
        tmaps
            Sorted transport maps as ``{{(t1, t2): tmat_2, (t2, t3): tmat_2, ...}}``.
        %(tmk_tmat.parameters)s

        Returns
        -------
        Merged transport maps into 1 :class:`anndata.AnnData` object.
        """
        from cellrank.tl.kernels import ConnectivityKernel

        last_time_point = LastTimePoint(last_time_point)
        tmaps = self._validate_tmaps(tmaps)

        conn_kwargs = dict(conn_kwargs)
        conn_kwargs["copy"] = False
        _ = conn_kwargs.pop("key_added", None)
        # same default as in `ConnectivityKernel`
        density_normalize = conn_kwargs.pop("density_normalize", True)

        blocks = [[None] * (len(tmaps) + 1) for _ in range(len(tmaps) + 1)]
        nrows, ncols = 0, 0
        obs_names, obs = [], []

        for i, tmap in enumerate(tmaps.values()):
            # tmap.X can be a view, shouldn't matter
            blocks[i][i + 1] = tmap.X
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
                .compute_transition_matrix(density_normalize=density_normalize)
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

    def _validate_tmaps(
        self,
        tmaps: Mapping[Tuple[Any, Any], AnnData],
        allow_reorder: bool = True,
    ) -> Mapping[Tuple[Any, Any], AnnData]:
        """
        Validate that transport maps conform to various invariants.

        Parameters
        ----------
        tmaps
            Transport maps where :attr:`anndata.AnnData.var_names` in the earlier time point must correspond to
            :attr:`anndata.AnnData.obs_names` in the later time point.
        allow_reorder
            Whether the target cells in the earlier transport map can be used to reorder
            the source cells of the later transport map.

        Returns
        -------
        Possibly reordered transport maps.

        Raises
        ------
        ValueError
            If ``allow_subset = False`` and the target/source cell names on earlier/later transport map don't match.
        KeyError
            If ``allow_subset = True`` and the target cells in earlier transport map are not found in the later one or
            if the transport maps haven't been computed for all cells in :attr:`adata`.
        """
        if not len(tmaps):
            raise ValueError("No transport maps have been computed.")

        tps, tmap2 = list(tmaps.keys()), None
        for (t1, t2), (t3, t4) in zip(tps[:-1], tps[1:]):
            tmap1, tmap2 = tmaps[t1, t2], tmaps[t3, t4]
            try:
                np.testing.assert_array_equal(tmap1.var_names, tmap2.obs_names)
            except AssertionError:
                if not allow_reorder:
                    raise ValueError(
                        f"Target cells of coupling with shape `{tmap1.shape}` at `{(t1, t2)}` does not "
                        f"match the source cells of coupling with shape `{tmap2.shape}` at `{(t3, t4)}`"
                    ) from None
                try:
                    # if a subset happens, we still assume that the transport maps cover all cells from `adata`
                    # i.e. the transport maps define a superset
                    tmaps[t3, t4] = tmap2[tmap1.var_names]  # keep the view
                except KeyError as e:
                    raise KeyError(
                        f"Unable to reorder transport map at `{(t3, t4)}` "
                        f"using transport map at `{(t1, t2)}`."
                    ) from e

        seen_obs = []
        for tmap in tmaps.values():
            seen_obs.extend(tmap.obs_names)
        seen_obs.extend(tmap.var_names)

        try:
            pd.testing.assert_series_equal(
                pd.Series(sorted(seen_obs)),
                pd.Series(sorted(self.adata.obs_names)),
                check_names=False,
                check_flags=False,
                check_index=True,
            )
        except AssertionError as e:
            raise KeyError(
                "Observations from transport maps don't match"
                "the observations from the underlying `AnnData` object."
            ) from e

        return tmaps

    @d.get_sections(base="tmk_thresh", sections=["Parameters"])
    def _threshold_transition_matrix(
        self,
        threshold: Union[float, Literal["auto"]],
        normalize: bool = True,
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
        else:
            if not (0 <= threshold <= 100):
                raise ValueError(
                    f"Expected `threshold` to be in `[0, 100]`, found `{threshold}`.`"
                )
            threshold = np.percentile(tmat.data, threshold)
        logg.info(f"Using `threshold={threshold}`")
        tmat.data[tmat.data <= threshold] = 0.0
        tmat.eliminate_zeros()

        if normalize:
            self._compute_transition_matrix(
                matrix=tmat,
                density_normalize=False,
                check_irreducibility=False,
            )
        else:
            # e.g. when in `_tmap_as_tmat`
            self._transition_matrix = tmat

    # TODO(michalk8): change t1/t2 to numeric/ordered (everywhere)
    def _tmat_to_adata(
        self, t1: Any, t2: Any, tmat: Union[np.ndarray, spmatrix, AnnData]
    ) -> AnnData:
        """Convert transport map ``tmat`` to :class:`anndata.AnnData`."""
        if isinstance(tmat, AnnData):
            return tmat

        tmat = AnnData(X=tmat)
        tmat.obs_names = self.adata[self.adata.obs[self._time_key] == t1].obs_names
        tmat.var_names = self.adata[self.adata.obs[self._time_key] == t2].obs_names

        return tmat

    @contextmanager
    def _tmap_as_tmat(
        self, threshold: Optional[Union[float, Literal["auto"]]] = None, **kwargs: Any
    ) -> None:
        """Temporarily set :attr:`transport_maps` as :attr:`transition_matrix`."""
        if self.transport_maps is None:
            raise RuntimeError(
                "Compute transport maps first as `.compute_transition_matrix()`."
            )

        tmat = self._transition_matrix
        try:
            self._transition_matrix = self._restich_tmaps(
                self.transport_maps, normalize=False, **kwargs
            ).X
            if threshold is not None:
                self._threshold_transition_matrix(threshold, normalize=False)
            yield
        finally:
            self._transition_matrix = tmat

    @d.dedent
    def plot_single_flow(
        self,
        cluster: str,
        cluster_key: str,
        time_key: Optional[str] = None,
        use_transport_maps: bool = False,
        threshold: Optional[Union[float, Literal["auto"]]] = None,
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
        if use_transport_maps:
            with self._tmap_as_tmat(threshold):
                return super().plot_single_flow(
                    cluster, cluster_key, time_key=time_key, **kwargs
                )
        return super().plot_single_flow(
            cluster, cluster_key, time_key=time_key, **kwargs
        )

    @property
    def transport_maps(self) -> Optional[Dict[Tuple[Any, Any], AnnData]]:
        """Transport maps for consecutive time pairs."""
        return self._tmaps
