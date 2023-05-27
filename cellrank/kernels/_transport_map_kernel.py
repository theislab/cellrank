from typing import (
    Any,
    Dict,
    Tuple,
    Union,
    Literal,
    Mapping,
    Iterable,
    Optional,
    Sequence,
)

from types import MappingProxyType
from tqdm.auto import tqdm
from contextlib import contextmanager

import scanpy as sc
from anndata import AnnData
from cellrank import logging as logg
from cellrank.settings import settings
from cellrank._utils._docs import d, inject_docs
from cellrank._utils._enum import ModeEnum
from cellrank._utils._utils import _normalize
from cellrank.kernels._base_kernel import UnidirectionalKernel

import numpy as np
import pandas as pd
import scipy.sparse as sp

__all__ = ["TransportMapKernel"]


class SelfTransitions(ModeEnum):
    UNIFORM = "uniform"
    DIAGONAL = "diagonal"
    CONNECTIVITIES = "connectivities"
    ALL = "all"


Numeric_t = Union[int, float]
Pair_t = Tuple[Numeric_t, Numeric_t]
Threshold_t = Union[int, float, Literal["auto", "auto_local"]]


class TransportMapKernel(UnidirectionalKernel):
    """Kernel base class which computes transition matrix based on transport maps for consecutive time point pairs."""

    def __init__(
        self,
        adata: AnnData,
        batch_key: str,
        couplings: Optional[
            Mapping[Pair_t, Union[np.ndarray, sp.spmatrix, AnnData]]
        ] = None,
        policy: Literal["sequential", "star"] = "sequential",
        reference: Optional[str] = None,
        **kwargs: AnnData,
    ):
        super().__init__(
            adata, batch_key=batch_key, reference=reference, backward=False, **kwargs
        )
        self._couplings: Optional[Dict[Pair_t, AnnData]] = None
        self._precomputed_couplings = couplings
        self._obs: Optional[pd.DataFrame] = None

        if policy == "sequential":
            batch = sorted(self.batch.cat.categories)
            self._keys = list(zip(batch[:-1], batch[1:]))
            self._reference = batch[-1]
        elif policy == "star":
            if reference not in self._batch_to_ix:
                raise ValueError("TODO")
            self._keys = [(c, reference) for c in self._batch_to_ix if c != reference]
            self._reference = reference
        else:
            raise NotImplementedError("TODO")

    def _read_from_adata(
        self,
        batch_key: str,
        reference: Optional[str] = None,
        **kwargs: Any,
    ) -> None:
        super()._read_from_adata(**kwargs)
        self._batch = self.adata.obs[batch_key].astype("category")
        cats = self.batch.cat.categories
        self._batch_to_ix = dict(zip(cats, range(len(cats))))

    def compute_coupling(
        self, src: Any, tgt: Any, **kwargs: Any
    ) -> Union[np.ndarray, sp.spmatrix, AnnData]:
        """Compute transport matrix for a time point pair.

        See :meth:`from_moscot` or :meth:`from_wot`.

        Parameters
        ----------
        src
            Source key in :attr:`batch`.
        tgt
            Target key in :attr:`batch`.
        kwargs
            Additional keyword arguments.

        Returns
        -------
        Array of shape ``(n_source_cells, n_target_cells)``. If :class:`anndata.AnnData`,
        :attr:`anndata.AnnData.obs_names` and :attr:`anndata.AnnData.var_names` correspond to subsets of observations
        from :attr:`adata`.
        """
        del kwargs
        if self._precomputed_couplings is None:
            raise NotImplementedError("TODO")
        return self._precomputed_couplings[src, tgt]

    @d.get_sections(base="tmk_thresh", sections=["Parameters"])
    def _sparsify_couplings(
        self,
        couplings: Dict[Pair_t, AnnData],
        threshold: Threshold_t,
        copy: bool = False,
    ) -> Optional[Mapping[Pair_t, AnnData]]:
        """Remove small non-zero values from :attr:`transition_matrix`.

        .. warning::
            Only dense :attr:`anndata.AnnData.X` will be sparsified.

        Parameters
        ----------
        couplings
            Optimal transport couplings to sparsify.
        threshold
            How to remove small non-zero values from the transition matrix. Valid options are:

                - `'auto'` - find the maximum threshold value which will not remove every non-zero value from any row.
                - `'auto_local'` - same as above, but done for each transport separately.
                - :class:`float` - value in :math:`[0, 100]` corresponding to a percentage of non-zeros to remove in
                  the couplings.

            Rows where all values are removed will have a uniform distribution and a warning will be issued.
        copy
            Whether to return a copy of the ``couplings`` or modify in-place.

        Returns
        -------
        If ``copy = True``, returns sparsified couplings. Otherwise, modifies ``couplings`` in-place.
        """
        if threshold == "auto":
            thresh = min(
                adata.X[i].max()
                for adata in couplings.values()
                for i in range(adata.n_obs)
            )
            logg.info(f"Using automatic `threshold={thresh}`")
        elif threshold == "auto_local":
            logg.info("Using automatic `threshold` for each time point pair separately")
        elif not (0 <= threshold <= 100):
            raise ValueError(
                f"Expected `threshold` to be in `[0, 100]`, found `{threshold}`.`"
            )

        for key, adata in couplings.items():
            if sp.issparse(adata.X):
                continue
            if copy:
                adata = adata.copy()
            tmat = adata.X

            if threshold == "auto_local":
                thresh = min(tmat[i].max() for i in range(tmat.shape[0]))
                logg.debug(f"Using `threshold={thresh}` at `{key}`")
            elif isinstance(threshold, (int, float)):
                thresh = np.percentile(tmat.data, threshold)
                logg.debug(f"Using `threshold={thresh}` at `{key}`")

            tmat = sp.csr_matrix(tmat, dtype=tmat.dtype)
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
            couplings[key] = AnnData(
                tmat, obs=adata.obs, var=adata.var, dtype=tmat.dtype
            )

        return couplings if copy else None

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
        """Compute transition matrix from optimal transport couplings.

        Parameters
        ----------
        threshold
            How to remove small non-zero values from the transition matrix. Valid options are:

                - `'auto'` - find the maximum threshold value which will not remove every non-zero value from any row.
                - `'auto_local'` - same as above, but done for each transport separately.
                - :class:`float` - value in :math:`[0, 100]` corresponding to a percentage of non-zeros to remove in
                  the couplings.

            Rows where all values are removed will have a uniform distribution and a warning will be issued.
        self_transitions
            How to define transitions within the diagonal blocks that correspond to transitions within the same
            time point. Valid options are:

                - `{st.UNIFORM!r}` - row-normalized matrix of 1s for transitions. Only applied to the last time point.
                - `{st.DIAGONAL!r}` - diagonal matrix with 1s on the diagonal. Only applied to the last time point.
                - `{st.CONNECTIVITIES!r}` - use transition matrix from :class:`cellrank.kernels.ConnectivityKernel`.
                  Only applied to the last time point.
                - :class:`typing.Sequence` - sequence of source time points defining which blocks should be weighted
                  by connectivities. Always applied to the last time point.
                - `{st.ALL!r}` - same as above, but for all keys.
        conn_weight
            Weight of connectivities self transitions. Only used when ``self_transitions = {st.ALL!r}`` or a sequence
            of source time points is passed.
        conn_kwargs
            Keyword arguments for :func:`scanpy.pp.neighbors` when using ``self_transitions`` use
            :class:`cellrank.kernels.ConnectivityKernel`. Can contain `'density_normalize'` for
            :meth:`cellrank.kernels.ConnectivityKernel.compute_transition_matrix`.
        kwargs
            Additional keyword arguments.

        Returns
        -------
        Self and updates the following attributes:

            - :attr:`transition_matrix` - transition matrix.
            - :attr:`transport_maps` - transport maps between consecutive time points.
        """
        cache_params = dict(kwargs)
        cache_params["threshold"] = threshold
        cache_params["self_transitions"] = str(self_transitions)
        if self._reuse_cache(cache_params):
            return self

        if (
            self_transitions == "all" or isinstance(self_transitions, (tuple, list))
        ) and (conn_weight is None or not (0 < conn_weight < 1)):
            raise ValueError(
                f"Expected `conn_weight` to be in interval `(0, 1)`, found `{conn_weight}`."
            )

        self._couplings = {}
        for src, tgt in tqdm(self._keys, unit="batch"):
            coupling = self.compute_coupling(src, tgt, **kwargs)
            self.couplings[src, tgt] = self._coupling_to_adata(src, tgt, coupling)

        if threshold is not None:
            self._sparsify_couplings(self.couplings, threshold=threshold, copy=False)

        tmap = self._restich_couplings(
            self.couplings,
            self_transitions=self_transitions,
            conn_weight=conn_weight,
            **conn_kwargs,
        )

        self.transition_matrix = tmap.X
        self._obs = tmap.obs

        return self

    @d.dedent
    def _restich_couplings(
        self,
        couplings: Mapping[Pair_t, AnnData],
        self_transitions: Union[
            str, SelfTransitions, Sequence[Numeric_t]
        ] = SelfTransitions.DIAGONAL,
        conn_weight: Optional[float] = None,
        **kwargs: Any,
    ) -> AnnData:
        """
        Group individual transport maps into 1 matrix aligned with :attr:`adata`.

        Parameters
        ----------
        couplings
            Optimal transport couplings.
        %(tmk_tmat.parameters)s

        Returns
        -------
        Merged transport maps into one :class:`anndata.AnnData` object.
        """
        if isinstance(self_transitions, (str, SelfTransitions)):
            self_transitions = SelfTransitions(self_transitions)
            if self_transitions == SelfTransitions.ALL:
                self_transitions = tuple(src for (src, _) in couplings.keys())
        elif isinstance(self_transitions, Iterable):
            self_transitions = tuple(self_transitions)
        else:
            raise TypeError(
                f"Expected `self_transitions` to be a `str` or a `Sequence`, "
                f"found `{type(self_transitions)}`."
            )
        self._validate_couplings(couplings)

        conn_kwargs = dict(kwargs)
        conn_kwargs["copy"] = False
        _ = conn_kwargs.pop("key_added", None)

        obs_names, obs = {}, {}
        blocks = [
            [None] * (len(self._batch_to_ix)) for _ in range(len(self._batch_to_ix))
        ]
        for (src, tgt), coupling in couplings.items():
            src_ix = self._batch_to_ix[src]
            tgt_ix = self._batch_to_ix[tgt]
            blocks[src_ix][tgt_ix] = coupling.X
            obs_names[src_ix] = coupling.obs_names
            obs[src_ix] = coupling.obs
            # to prevent blocks from disappearing
            n = np.sum(self.batch == src)
            blocks[src_ix][src_ix] = sp.spdiags([0] * n, 0, n, n)
        obs_names[tgt_ix] = coupling.var_names
        obs_names = np.ravel([obs_names[ix] for ix in range(len(blocks))])

        # if policy='sequential', reference is the last key
        # if policy='star', reference is supplied by the user
        ref_ix = self._batch_to_ix[self._reference]
        ref_mask = self.batch == self._reference
        n = int(np.sum(ref_mask))

        if self_transitions == SelfTransitions.DIAGONAL:
            blocks[ref_ix][ref_ix] = sp.spdiags([1] * n, 0, n, n)
        elif self_transitions == SelfTransitions.UNIFORM:
            blocks[ref_ix][ref_ix] = np.ones((n, n)) / float(n)
        elif self_transitions == SelfTransitions.CONNECTIVITIES:
            blocks[ref_ix][ref_ix] = self._compute_connectivity_tmat(
                self.adata[ref_mask], **conn_kwargs
            )
        elif isinstance(self_transitions, tuple):
            verbosity = settings.verbosity
            try:  # ignore overly verbose logging
                settings.verbosity = 0
                for (src, tgt), coupling in couplings.items():
                    # fmt: off
                    if src not in self_transitions:
                        continue
                    src_ix, tgt_ix = self._batch_to_ix[src], self._batch_to_ix[tgt]
                    blocks[src_ix][src_ix] = conn_weight * self._compute_connectivity_tmat(
                        self.adata[coupling.obs_names], **conn_kwargs
                    )
                    blocks[src_ix][tgt_ix] = (1 - conn_weight) * _normalize(blocks[src_ix][tgt_ix])
                    # fmt: on
                blocks[ref_ix][ref_ix] = self._compute_connectivity_tmat(
                    self.adata[ref_mask], **conn_kwargs
                )
            finally:
                settings.verbosity = verbosity
        else:
            raise NotImplementedError(
                f"Self transitions' mode `{self_transitions}` is not yet implemented."
            )

        tmp = AnnData(sp.bmat(blocks, format="csr"), dtype="float64")
        tmp.obs_names = obs_names
        tmp.var_names = obs_names
        tmp = tmp[self.adata.obs_names, :][:, self.adata.obs_names]

        tmp.obs = pd.merge(
            tmp.obs,
            pd.concat([obs.get(ix, None) for ix in range(len(blocks))]),
            left_index=True,
            right_index=True,
            how="left",
        )

        return tmp

    def _validate_couplings(
        self,
        couplings: Mapping[Pair_t, AnnData],
    ) -> None:
        """
        Validate that transport maps conform to various invariants.

        Parameters
        ----------
        couplings
            Transport maps where :attr:`anndata.AnnData.var_names` in the earlier time point must correspond to
            :attr:`anndata.AnnData.obs_names` in the later time point.

        Returns
        -------
        Possibly reordered transport maps.
        """
        seen_obs = []
        for coupling in couplings.values():
            seen_obs.extend(coupling.obs_names)
        seen_obs.extend(coupling.var_names)

        # this invariant holds for both policies
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

    def _coupling_to_adata(
        self, src: Any, tgt: Any, coupling: Union[np.ndarray, sp.spmatrix, AnnData]
    ) -> AnnData:
        """Convert the coupling to :class:`~anndata.AnnData`."""
        if not isinstance(coupling, AnnData):
            coupling = AnnData(X=coupling, dtype=coupling.dtype)
            coupling.obs_names = self.adata[self.batch == src].obs_names
            coupling.var_names = self.adata[self.batch == tgt].obs_names

        if sp.issparse(coupling.X) and not sp.isspmatrix_csr(coupling.X):
            coupling.X = coupling.X.tocsr()

        return coupling

    @contextmanager
    def _tmap_as_tmat(self, **kwargs: Any) -> None:
        """Temporarily set :attr:`transport_mats` as :attr:`transition_matrix`.

        Parameters
        ----------
        kwargs
            Keyword arguments for :meth:`_restich_tmaps`.

        Returns
        -------
        Nothing.
        """
        if self.couplings is None:
            raise RuntimeError(
                "Compute transport matrices first as `.compute_transition_matrix()`."
            )

        tmat = self._transition_matrix
        try:
            # fmt: off
            self._transition_matrix = self._restich_couplings(self.couplings, **kwargs).X
            # fmt: on
            yield
        finally:
            self._transition_matrix = tmat

    @staticmethod
    def _compute_connectivity_tmat(
        adata: AnnData, **kwargs: Any
    ) -> Union[np.ndarray, sp.spmatrix]:
        from cellrank.kernels import ConnectivityKernel

        # same default as in `ConnectivityKernel`
        density_normalize = kwargs.pop("density_normalize", True)
        sc.pp.neighbors(adata, **kwargs)

        return (
            ConnectivityKernel(adata)
            .compute_transition_matrix(density_normalize=density_normalize)
            .transition_matrix
        )

    @property
    def couplings(self) -> Optional[Dict[Pair_t, AnnData]]:
        """Optimal transport couplings."""
        return self._couplings

    @property
    def batch(self) -> pd.Series:
        """TODO."""
        return self._batch

    @property
    def obs(self) -> Optional[pd.DataFrame]:
        """Cell-level metadata."""
        return self._obs
