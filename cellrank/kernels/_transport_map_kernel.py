from typing import Any, Dict, Tuple, Union, Mapping, Iterable, Optional, Sequence
from typing_extensions import Literal

from abc import ABC, abstractmethod
from types import MappingProxyType
from contextlib import contextmanager

import scanpy as sc
from anndata import AnnData
from cellrank import logging as logg
from cellrank.settings import settings
from cellrank._utils._docs import d, inject_docs
from cellrank._utils._enum import ModeEnum
from cellrank._utils._utils import _normalize
from cellrank.kernels._experimental_time_kernel import ExperimentalTimeKernel

import numpy as np
import pandas as pd
from scipy.sparse import bmat, spdiags, spmatrix, csr_matrix

__all__ = ["TransportMapKernel"]


class SelfTransitions(ModeEnum):
    UNIFORM = "uniform"
    DIAGONAL = "diagonal"
    CONNECTIVITIES = "connectivities"
    ALL = "all"


Numeric_t = Union[int, float]
Pair_t = Tuple[Numeric_t, Numeric_t]
Threshold_t = Union[int, float, Literal["auto", "auto_local"]]


class TransportMapKernel(ExperimentalTimeKernel, ABC):
    """Kernel base class which computes transition matrix based on transport maps for consecutive time point pairs."""

    def __init__(self, *args: Any, **kwargs: Any):
        super().__init__(*args, **kwargs)
        self._tmaps: Optional[Dict[Pair_t, AnnData]] = None

    @abstractmethod
    def _compute_tmap(
        self, t1: Numeric_t, t2: Numeric_t, **kwargs: Any
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

    @d.get_sections(base="tmk_thresh", sections=["Parameters"])
    def _threshold_transport_maps(
        self,
        tmaps: Dict[Pair_t, AnnData],
        threshold: Threshold_t,
        copy: bool = False,
    ) -> Optional[Mapping[Pair_t, AnnData]]:
        """
        Remove small non-zero values from :attr:`transition_matrix`.

        Parameters
        ----------
        threshold
            How to remove small non-zero values from the transition matrix. Valid options are:

                - `'auto'` - find the maximum threshold value which will not remove every non-zero value from any row.
                - `'auto_local'` - same as above, but done for each transport separately.
                - :class:`float` - value in `[0, 100]` corresponding to a percentage of non-zeros to remove in each
                  transport map.

            Rows where all values are removed will have uniform distribution and a warning will be issued.
        copy
            Whether to return a copy of the ``tmaps`` or modify in-place.

        Returns
        -------
        If ``copy = True``, returns a thresholded transport maps. Otherwise, modifies ``tmaps`` in-place.
        """
        if threshold == "auto":
            thresh = min(
                adata.X[i].max() for adata in tmaps.values() for i in range(adata.n_obs)
            )
            logg.info(f"Using automatic `threshold={thresh}`")
        elif threshold == "auto_local":
            logg.info("Using automatic `threshold` for each time point pair separately")
        elif not (0 <= threshold <= 100):
            raise ValueError(
                f"Expected `threshold` to be in `[0, 100]`, found `{threshold}`.`"
            )

        for key, adata in tmaps.items():
            if copy:
                adata = adata.copy()
            tmat = adata.X

            if threshold == "auto_local":
                thresh = min(tmat[i].max() for i in range(tmat.shape[0]))
                logg.debug(f"Using `threshold={thresh}` at `{key}`")
            elif isinstance(threshold, (int, float)):
                thresh = np.percentile(tmat.data, threshold)
                logg.debug(f"Using `threshold={thresh}` at `{key}`")

            tmat = csr_matrix(tmat, dtype=tmat.dtype)
            tmat.data[tmat.data < thresh] = 0.0
            zeros_mask = np.where(np.asarray(tmat.sum(1)).squeeze() == 0)[0]
            if np.any(zeros_mask):
                logg.warning(
                    f"After thresholding, `{len(zeros_mask)}` row(s) of transport at `{key}` are forced to be uniform"
                )
                for ix in zeros_mask:
                    start, end = tmat.indptr[ix], tmat.indptr[ix + 1]
                    size = end - start
                    tmat.data[start:end] = np.ones(size, dtype=tmat.dtype) / size

            # after `zeros_mask` has been handled, to have access to removed row indices
            tmat.eliminate_zeros()
            tmaps[key] = AnnData(tmat, obs=adata.obs, var=adata.var, dtype=tmat.dtype)

        return tmaps if copy else None

    @d.get_sections(base="tmk_tmat", sections=["Parameters"])
    @d.dedent
    @inject_docs(st=SelfTransitions)
    def compute_transition_matrix(
        self,
        threshold: Optional[Threshold_t] = "auto",
        self_transitions: Union[
            Literal["uniform", "diagonal", "connectivities", "all"],
            Sequence[Numeric_t],
        ] = SelfTransitions.CONNECTIVITIES,
        conn_weight: Optional[float] = None,
        conn_kwargs: Mapping[str, Any] = MappingProxyType({}),
        **kwargs: Any,
    ) -> "TransportMapKernel":
        """
        Compute transition matrix using transport maps.

        Parameters
        ----------
        %(tmk_thresh.parameters)s
        self_transitions
            How to define transitions within the diagonal blocks that correspond to transitions within the same
            time point. Valid options are:

                - `{st.UNIFORM!r}` - row-normalized matrix of 1s for transitions. Only applied to the last time point.
                - `{st.DIAGONAL!r}` - diagonal matrix with 1s on the diagonal. Only applied to the last time point.
                - `{st.CONNECTIVITIES!r}` - use transition matrix from :class:`cellrank.kernels.ConnectivityKernel`.
                  Only applied to the last time point.
                - :class:`typing.Sequence` - sequence of source time points defining which blocks should be weighted
                  by connectivities. Always applied to the last time point.
                - `{st.ALL!r}` - same as above, but for all time points.
        conn_weight
            Weight of connectivities self transitions. Only used when ``self_transitions = {st.ALL!r}`` or a sequence
            of source time points is passed.
        conn_kwargs
            Keyword arguments for :func:`scanpy.pp.neighbors` when using ``self_transitions`` use
            :class:`cellrank.kernels.ConnectivityKernel`. Can contain `'density_normalize'` for
            :meth:`cellrank.kernels.ConnectivityKernel.compute_transition_matrix`.

        Returns
        -------
        Self and updates the following attributes:

            - :attr:`transition_matrix` - transition matrix.
            - :attr:`transport_maps` - transport maps between consecutive time points.
        """
        cache_params = dict(kwargs)
        cache_params["threshold"] = threshold
        cache_params["self_transitions"] = self_transitions
        if SelfTransitions(self_transitions) == SelfTransitions.ALL:
            conn_kwargs["conn_weight"] = conn_weight
        if self._reuse_cache(cache_params):
            return self

        timepoints = self.experimental_time.cat.categories
        timepoints = list(zip(timepoints[:-1], timepoints[1:]))

        self._tmaps = {
            (t1, t2): self._tmat_to_adata(t1, t2, self._compute_tmap(t1, t2, **kwargs))
            for t1, t2 in timepoints
        }
        if threshold is not None:
            self._threshold_transport_maps(self._tmaps, threshold=threshold, copy=False)

        self.transition_matrix = self._restich_tmaps(
            self._tmaps,
            self_transitions=self_transitions,
            conn_weight=conn_weight,
            conn_kwargs=conn_kwargs,
        ).X

        return self

    @d.dedent
    def _restich_tmaps(
        self,
        tmaps: Mapping[Pair_t, AnnData],
        self_transitions: Union[
            str, SelfTransitions, Sequence[Numeric_t]
        ] = SelfTransitions.DIAGONAL,
        conn_weight: Optional[float] = None,
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
        Merged transport maps into one :class:`anndata.AnnData` object.
        """
        if isinstance(self_transitions, (str, SelfTransitions)):
            self_transitions = SelfTransitions(self_transitions)
            if self_transitions == SelfTransitions.ALL:
                self_transitions = tuple(t1 for (t1, _) in tmaps.keys())
        elif isinstance(self_transitions, Iterable):
            self_transitions = tuple(self_transitions)
        else:
            raise TypeError(
                f"Expected `self_transitions` to be a `str` or a `Sequence`, "
                f"found `{type(self_transitions).__name__}`."
            )
        tmaps = self._validate_tmaps(tmaps)

        conn_kwargs = dict(conn_kwargs)
        conn_kwargs["copy"] = False
        _ = conn_kwargs.pop("key_added", None)

        blocks = [[None] * (len(tmaps) + 1) for _ in range(len(tmaps) + 1)]
        nrows, ncols = 0, 0
        obs_names, obs = [], []

        for i, tmap in enumerate(tmaps.values()):
            blocks[i][i + 1] = tmap.X
            nrows += tmap.n_obs
            ncols += tmap.n_vars
            obs_names.extend(tmap.obs_names)
            obs.append(tmap.obs)
        obs_names.extend(tmap.var_names)

        n = self.adata.n_obs - nrows
        if self_transitions == SelfTransitions.DIAGONAL:
            blocks[-1][-1] = spdiags([1] * n, 0, n, n)
        elif self_transitions == SelfTransitions.UNIFORM:
            blocks[-1][-1] = np.ones((n, n)) / float(n)
        elif self_transitions == SelfTransitions.CONNECTIVITIES:
            blocks[-1][-1] = self._compute_connectivity_tmat(
                self.adata[tmap.var_names], **conn_kwargs
            )
        elif isinstance(self_transitions, tuple):
            verbosity = settings.verbosity
            try:  # ignore overly verbose logging
                settings.verbosity = 0
                if conn_weight is None or not (0 < conn_weight < 1):
                    raise ValueError(
                        "Please specify `conn_weight` in interval `(0, 1)`."
                    )
                for i, ((t1, _), tmap) in enumerate(tmaps.items()):
                    if t1 not in self_transitions:
                        continue
                    blocks[i][i] = conn_weight * self._compute_connectivity_tmat(
                        self.adata[tmap.obs_names], **conn_kwargs
                    )
                    blocks[i][i + 1] = (1 - conn_weight) * _normalize(blocks[i][i + 1])
                blocks[-1][-1] = self._compute_connectivity_tmat(
                    self.adata[tmap.var_names], **conn_kwargs
                )
            finally:
                settings.verbosity = verbosity
        else:
            raise NotImplementedError(
                f"Self transitions' mode `{self_transitions}` is not yet implemented."
            )

        if not isinstance(self_transitions, tuple):
            # prevent the last block from disappearing
            n = blocks[0][1].shape[0]
            blocks[0][0] = spdiags([], 0, n, n)

        tmp = AnnData(bmat(blocks, format="csr"), dtype="float64")
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
        tmaps: Dict[Pair_t, AnnData],
        allow_reorder: bool = True,
    ) -> Mapping[Pair_t, AnnData]:
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
                "Observations from transport maps don't match "
                "the observations from the underlying `AnnData` object."
            ) from e

        return tmaps

    def _tmat_to_adata(
        self, t1: Numeric_t, t2: Numeric_t, tmat: Union[np.ndarray, spmatrix, AnnData]
    ) -> AnnData:
        """Convert transport map ``tmat`` to :class:`anndata.AnnData`. Do nothing if ``tmat`` is of the correct type."""
        if isinstance(tmat, AnnData):
            return tmat

        tmat = AnnData(X=tmat, dtype=tmat.dtype)
        tmat.obs_names = self.adata[self.adata.obs[self._time_key] == t1].obs_names
        tmat.var_names = self.adata[self.adata.obs[self._time_key] == t2].obs_names

        return tmat

    @contextmanager
    def _tmap_as_tmat(self, **kwargs: Any) -> None:
        """Temporarily set :attr:`transport_maps` as :attr:`transition_matrix`."""
        if self.transport_maps is None:
            raise RuntimeError(
                "Compute transport maps first as `.compute_transition_matrix()`."
            )

        tmat = self._transition_matrix
        try:
            # fmt: off
            self._transition_matrix = self._restich_tmaps(self.transport_maps, **kwargs).X
            # fmt: on
            yield
        finally:
            self._transition_matrix = tmat

    @staticmethod
    def _compute_connectivity_tmat(
        adata: AnnData, **kwargs: Any
    ) -> Union[np.ndarray, spmatrix]:
        from cellrank.kernels import ConnectivityKernel

        # same default as in `ConnectivityKernel`
        density_normalize = kwargs.pop("density_normalize", True)
        sc.pp.neighbors(adata, **kwargs)

        return (
            ConnectivityKernel(adata)
            .compute_transition_matrix(density_normalize=density_normalize)
            .transition_matrix
        )

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
    def transport_maps(self) -> Optional[Dict[Pair_t, AnnData]]:
        """Transport maps for consecutive time pairs."""
        return self._tmaps

    def __invert__(self) -> "TransportMapKernel":
        tk = super().__invert__("_tmaps")  # don't copy transport maps
        tk._tmaps = None
        return tk
