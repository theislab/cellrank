import itertools
import os
import pathlib
import types
from typing import (
    TYPE_CHECKING,
    Any,
    Dict,
    Iterable,
    Literal,
    Mapping,
    Optional,
    Sequence,
    Tuple,
    Union,
)

from tqdm.auto import tqdm

import numpy as np
import pandas as pd
import scipy.sparse as sp
from pandas.api.types import infer_dtype

import scanpy as sc
from anndata import AnnData

from cellrank import logging as logg
from cellrank._utils._docs import d, inject_docs
from cellrank._utils._enum import ModeEnum
from cellrank._utils._utils import _normalize
from cellrank.kernels._base_kernel import UnidirectionalKernel
from cellrank.settings import settings

__all__ = ["RealTimeKernel"]

if TYPE_CHECKING:
    from moscot.problems.spatiotemporal import SpatioTemporalProblem
    from moscot.problems.time import LineageProblem, TemporalProblem


class SelfTransitions(ModeEnum):
    UNIFORM = "uniform"
    DIAGONAL = "diagonal"
    CONNECTIVITIES = "connectivities"
    ALL = "all"


Key_t = Tuple[Any, Any]
Threshold_t = Union[int, float, Literal["auto", "auto_local"]]
Coupling_t = Union[np.ndarray, sp.spmatrix, AnnData]


# TODO(michalk8): subclass the `ExperimentalTimeKernel`
@d.dedent
class RealTimeKernel(UnidirectionalKernel):
    """Kernel which computes transition matrix using optimal transport couplings.

    .. seealso::
        - See :doc:`../../../notebooks/tutorials/kernels/500_real_time` on how to
          compute the :attr:`~cellrank.kernels.RealTimeKernel.transition_matrix` between experimental time points
          using `optimal transport <https://en.wikipedia.org/wiki/Transportation_theory_(mathematics)>`_.

    This class should be constructed using either:

    1. the :meth:`from_moscot` or :meth:`from_wot` method,
    2. explicitly passing the pre-computed couplings, or
    3. overriding the :meth:`compute_coupling` method.

    Parameters
    ----------
    adata
        Annotated data object.
    time_key
        Key in :attr:`~anndata.AnnData.obs` containing the experimental time.
    couplings
        Pre-computed transport couplings. The keys should correspond to a :class:`tuple` of categories from
        the :attr:`time`. If :obj:`None`, the keys will be constructed using the ``policy``
        and the :meth:`compute_coupling` method must be overriden.
    policy
        How to construct the keys from the :attr:`time` when ``couplings = None``:

        - ``policy = 'sequential'`` - the keys will be set to ``[(t1, t2), (t2, t3), ...]``.
        - ``policy = 'triu'`` - the keys will be set to ``[(t1, t2), (t1, t3), ..., (t2, t3), ...]``.
    kwargs
        Keyword arguments for the :class:`~cellrank.kernels.Kernel`.
    """

    def __init__(
        self,
        adata: AnnData,
        time_key: str,
        couplings: Optional[Mapping[Key_t, Optional[Coupling_t]]] = None,
        policy: Literal["sequential", "triu"] = "sequential",
        **kwargs: Any,
    ):
        super().__init__(adata, time_key=time_key, **kwargs)
        self._obs: Optional[pd.DataFrame] = None
        self._couplings, self._reference = self._get_default_coupling(couplings, policy)

    def _read_from_adata(
        self,
        time_key: str,
        **kwargs: Any,
    ) -> None:
        super()._read_from_adata(**kwargs)
        self._time = self.adata.obs[time_key].copy()
        if not isinstance(self._time.dtype, pd.CategoricalDtype):
            raise TypeError(f"Expected `adata.obs[{time_key!r}]` to be categorical, found `{infer_dtype(self._time)}`.")
        self._time = self._time.cat.remove_unused_categories()
        cats = self._time.cat.categories
        self._time_to_ix = dict(zip(cats, range(len(cats))))

    def compute_coupling(self, src: Any, tgt: Any, **kwargs: Any) -> Coupling_t:
        """Compute coupling for a given time pair.

        .. note::
            This implementation only looks up the values in :attr:`couplings`,
            see :meth:`from_moscot` or :meth:`from_wot` on how to initialize them.

        Parameters
        ----------
        src
            Source key in :attr:`time`.
        tgt
            Target key in :attr:`time`.
        kwargs
            Additional keyword arguments.

        Returns
        -------
        Array of shape ``(n_source_cells, n_target_cells)``.
        """
        del kwargs
        if (src, tgt) not in self.couplings:
            raise KeyError(f"Key `{src, tgt}` not found in the `couplings`.")
        if self.couplings[src, tgt] is None:
            raise ValueError(
                f"Coupling for `{src, tgt}` is not computed. " f"Consider overriding the `compute_coupling` method."
            )
        return self.couplings[src, tgt]

    @inject_docs(st=SelfTransitions)
    def compute_transition_matrix(
        self,
        threshold: Optional[Threshold_t] = "auto",
        self_transitions: Union[
            Literal["uniform", "diagonal", "connectivities", "all"],
            Sequence[Any],
        ] = SelfTransitions.CONNECTIVITIES,
        conn_weight: Optional[float] = None,
        conn_kwargs: Mapping[str, Any] = types.MappingProxyType({}),
        **kwargs: Any,
    ) -> "RealTimeKernel":
        """Compute transition matrix from optimal transport couplings.

        Parameters
        ----------
        threshold
            How to remove small non-zero values from the transition matrix. Valid options are:

            - ``'auto'`` - find the maximum threshold value which will not remove every non-zero value from any row.
            - ``'auto_local'`` - same as above, but done for each transport separately.
            - :class:`float` - value in :math:`[0, 100]` corresponding to a percentage of non-zeros to remove in
              the couplings.

            Rows where all values are removed will have a uniform distribution and a warning will be issued.
        self_transitions
            How to define transitions within the blocks that correspond to transitions within the same key.
            Valid options are:

            - ``{st.UNIFORM!r}`` - row-normalized matrix of :math:`1`.
            - ``{st.DIAGONAL!r}`` - identity matrix.
            - ``{st.CONNECTIVITIES!r}`` - transition matrix from the :class:`~cellrank.kernels.ConnectivityKernel`.
            - :class:`~typing.Sequence`` - sequence of source keys defining which blocks should be weighted
              by the connectivities.
            - ``{st.ALL!r}`` - same as above, but for all keys.

            The first 3 options are applied to the block specified by the :attr:`reference`.
        conn_weight
            Weight of connectivities' self transitions. Only used when ``self_transitions = {st.ALL!r}`` or
            a sequence of source keys is passed.
        conn_kwargs
            Keyword arguments for :func:`~scanpy.pp.neighbors` or
            :meth:`~cellrank.kernels.ConnectivityKernel.compute_transition_matrix` when using
            ``self_transitions = 'connectivities'``.
        kwargs
            Keyword arguments for :meth:`compute_coupling`.

        Returns
        -------
        Returns self and updates :attr:`transition_matrix`, :attr:`couplings` and :attr:`params`.
        """
        cache_params = dict(kwargs)
        cache_params["threshold"] = threshold
        cache_params["self_transitions"] = str(self_transitions)
        if self._reuse_cache(cache_params):
            return self

        if (self_transitions == "all" or isinstance(self_transitions, (tuple, list))) and (
            conn_weight is None or not (0 < conn_weight < 1)
        ):
            raise ValueError(f"Expected `conn_weight` to be in interval `(0, 1)`, found `{conn_weight}`.")

        for src, tgt in tqdm(self._couplings, unit="time pair"):
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

    @classmethod
    def from_moscot(
        cls,
        problem: "Union[TemporalProblem, LineageProblem, SpatioTemporalProblem]",
        sparse_mode: Optional[Literal["threshold", "percentile", "min_row"]] = None,
        sparsify_kwargs: Mapping[str, Any] = types.MappingProxyType({}),
        copy: bool = False,
        **kwargs: Any,
    ) -> "RealTimeKernel":
        """Construct the kernel from :mod:`moscot` :cite:`klein:23`.

        Parameters
        ----------
        problem
            :mod:`moscot` problem.
        sparse_mode
            Sparsification mode for :meth:`~moscot.base.output.BaseSolverOutput.sparsify`.
            If :obj:`None`, do not sparsify the outputs. Note that :meth:`compute_transition_matrix`
            can also sparsify the final :attr:`transition_matrix`.
        sparsify_kwargs
            Keyword arguments for the sparsification.
        copy
            Whether to copy the underlying arrays. Note that :class:`jax arrays <jax.Array>` are always copied.
        kwargs
            Keyword arguments for :class:`~cellrank.kernels.RealTimeKernel`.

        Returns
        -------
        The kernel.

        Examples
        --------
        .. code-block:: python

            import moscot as mt
            import cellrank as cr

            adata = mt.datasets.hspc()
            adata.obs["day"] = adata.obs["day"].astype("category")

            problem = mt.problems.TemporalProblem(adata)
            problem = problem.prepare(time_key="day").solve()

            rtk = cr.kernels.RealTimeKernel.from_moscot(problem)
            rtk = rtk.compute_transition_matrix()
        """
        from moscot.utils.subset_policy import SequentialPolicy, TriangularPolicy

        if not problem.solutions:
            raise RuntimeError("Problem contains no solutions. Please run `problem.solve(...)` first.")

        policy = problem._policy
        if isinstance(policy, SequentialPolicy):
            policy = "sequential"
        elif isinstance(policy, TriangularPolicy):
            policy = "triu"
        else:
            raise NotImplementedError(f"Handling `{type(policy)}` policy is not yet implemented.")

        couplings = {}
        for (t1, t2), solution in problem.solutions.items():
            adata_src = problem[t1, t2].adata_src
            adata_tgt = problem[t1, t2].adata_tgt
            if sparse_mode is not None:
                solution = solution.sparsify(mode=sparse_mode, **sparsify_kwargs)
            coupling = solution.transport_matrix
            if not (isinstance(coupling, np.ndarray) or sp.issparse(coupling)):
                # convert from, e.g., `jax`
                # we always copy because otherwise the array is read-only and we may sparsify
                coupling = np.array(coupling, copy=True)
            elif copy:
                coupling = coupling.copy()
            couplings[t1, t2] = AnnData(coupling, obs=adata_src.obs, var=adata_tgt.obs)

        return cls(
            problem.adata,
            couplings=couplings,
            time_key=problem.temporal_key,
            policy=policy,
            **kwargs,
        )

    @classmethod
    def from_wot(
        cls,
        adata: AnnData,
        path: Union[str, pathlib.Path],
        time_key: str,
        **kwargs: Any,
    ) -> "RealTimeKernel":
        """Construct the kernel from Waddington-OT :cite:`schiebinger:19`.

        Parameters
        ----------
        adata
            Annotated data object.
        path
            Directory where the couplings are stored.
        time_key
            Key in :attr:`~anndata.AnnData.obs` containing the experimental time.
        kwargs
            Keyword arguments for :class:`~cellrank.kernels.RealTimeKernel`.

        Returns
        -------
        The kernel.

        Examples
        --------
        .. code-block:: python

            import wot
            import cellrank as cr

            adata = cr.datasets.reprogramming_schiebinger(subset_to_serum=True)
            adata.obs["day"] = adata.obs["day"].astype(float).astype("category")

            ot_model = wot.ot.OTModel(adata, day_field="day")
            ot_model.compute_all_transport_maps(tmap_out="tmaps/")

            rtk = cr.kernels.RealTimeKernel.from_wot(adata, path="tmaps/", time_key="day")
            rtk = rtk.compute_transition_matrix()
        """
        path = pathlib.Path(path)
        dtype = type(adata.obs[time_key].iloc[0])

        couplings = {}
        for fname in path.glob("*h5ad"):
            name, _ = os.path.splitext(fname)
            *_, src, tgt = name.split("_")
            couplings[dtype(src), dtype(tgt)] = sc.read(fname)

        return cls(adata, couplings=couplings, time_key=time_key, policy="sequential", **kwargs)

    @d.dedent
    def _restich_couplings(
        self,
        couplings: Mapping[Key_t, AnnData],
        self_transitions: Union[str, SelfTransitions, Sequence[Any]] = SelfTransitions.DIAGONAL,
        conn_weight: Optional[float] = None,
        **kwargs: Any,
    ) -> AnnData:
        """Group individual transport maps into 1 matrix aligned with :attr:`adata`.

        Parameters
        ----------
        couplings
            Optimal transport couplings.
        self_transitions
            How to define transitions within the blocks that correspond to transitions within the same key.
        conn_weight
            Weight of connectivities' self transitions. Only used when ``self_transitions = {st.ALL!r}`` or
            a sequence of source keys is passed.
        kwargs
            Keyword arguments for :func:`~scanpy.pp.neighbors` or
            :meth:`~cellrank.kernels.ConnectivityKernel.compute_transition_matrix` when using
            ``self_transitions = 'connectivities'``.

        Returns
        -------
        Merged transport maps into one :class:`~anndata.AnnData` object.
        """
        if isinstance(self_transitions, (str, SelfTransitions)):
            self_transitions = SelfTransitions(self_transitions)
            if self_transitions == SelfTransitions.ALL:
                self_transitions = tuple(src for (src, _) in couplings)
        elif isinstance(self_transitions, Iterable):
            self_transitions = tuple(self_transitions)
        else:
            raise TypeError(
                f"Expected `self_transitions` to be a `str` or a `Sequence`, " f"found `{type(self_transitions)}`."
            )
        self._validate_couplings(couplings)

        conn_kwargs = dict(kwargs)
        conn_kwargs["copy"] = False
        _ = conn_kwargs.pop("key_added", None)

        obs_names, obs = {}, {}
        blocks = [[None] * (len(self._time_to_ix)) for _ in range(len(self._time_to_ix))]
        for (src, tgt), coupling in couplings.items():
            src_ix = self._time_to_ix[src]
            tgt_ix = self._time_to_ix[tgt]
            blocks[src_ix][tgt_ix] = coupling.X
            obs_names[src_ix] = coupling.obs_names
            obs[src_ix] = coupling.obs
            # to prevent blocks from disappearing
            n = np.sum(self._time == src)
            blocks[src_ix][src_ix] = sp.spdiags([0] * n, 0, n, n)
            if tgt == self._reference:
                obs_names[tgt_ix] = coupling.var_names

        # if policy='sequential', reference is the last key
        ref_ix = self._time_to_ix[self._reference]
        ref_mask = self._time == self._reference
        n = int(np.sum(ref_mask))

        if self_transitions == SelfTransitions.DIAGONAL:
            blocks[ref_ix][ref_ix] = sp.spdiags([1] * n, 0, n, n)
        elif self_transitions == SelfTransitions.UNIFORM:
            blocks[ref_ix][ref_ix] = np.ones((n, n)) / float(n)
        elif self_transitions == SelfTransitions.CONNECTIVITIES:
            blocks[ref_ix][ref_ix] = _compute_connectivity_tmat(self.adata[ref_mask], **conn_kwargs)
        elif isinstance(self_transitions, tuple):
            verbosity = settings.verbosity
            try:  # ignore overly verbose logging
                settings.verbosity = 0
                for (src, tgt), coupling in couplings.items():
                    # fmt: off
                    if src not in self_transitions:
                        continue
                    src_ix, tgt_ix = self._time_to_ix[src], self._time_to_ix[tgt]
                    blocks[src_ix][src_ix] = conn_weight * _compute_connectivity_tmat(
                        self.adata[coupling.obs_names], **conn_kwargs
                    )
                    blocks[src_ix][tgt_ix] = (1 - conn_weight) * _normalize(blocks[src_ix][tgt_ix])
                    # fmt: on
                blocks[ref_ix][ref_ix] = _compute_connectivity_tmat(self.adata[ref_mask], **conn_kwargs)
            finally:
                settings.verbosity = verbosity
        else:
            raise NotImplementedError(f"Self transitions' mode `{self_transitions}` is not yet implemented.")

        index = []
        for ix in range(len(blocks)):
            index.extend(obs_names[ix])

        tmp = AnnData(sp.bmat(blocks, format="csr"))
        tmp.obs_names = index
        tmp.var_names = index
        tmp = tmp[self.adata.obs_names, :][:, self.adata.obs_names]

        obs = pd.concat([obs.get(ix, None) for ix in range(len(blocks))])
        tmp.obs = pd.merge(
            tmp.obs,
            obs,
            left_index=True,
            right_index=True,
            how="left",
        )
        return tmp

    @d.get_sections(base="tmk_thresh", sections=["Parameters"])
    def _sparsify_couplings(
        self,
        couplings: Dict[Key_t, AnnData],
        threshold: Threshold_t,
        copy: bool = False,
    ) -> Optional[Mapping[Key_t, AnnData]]:
        """Remove small non-zero values from :attr:`transition_matrix`.

        .. warning::
            Only dense :attr:`~anndata.AnnData.X` will be sparsified.

        Parameters
        ----------
        couplings
            Optimal transport couplings to sparsify.
        threshold
            How to remove small non-zero values from the transition matrix. Valid options are:

            - ``'auto'`` - find the maximum threshold value which will not remove every non-zero value from any row.
            - ``'auto_local'`` - same as above, but done for each transport separately.
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
            thresh = min(adata.X[i].max() for adata in couplings.values() for i in range(adata.n_obs))
            logg.info(f"Using automatic `threshold={thresh}`")
        elif threshold == "auto_local":
            logg.info("Using automatic `threshold` for each coupling separately")
        elif not (0 <= threshold <= 100):
            raise ValueError(f"Expected `threshold` to be in `[0, 100]`, found `{threshold}`.`")

        for key, adata in couplings.items():
            if sp.issparse(adata.X):
                continue
            if copy:
                adata = adata.copy()
            tmat = adata.X

            if threshold == "auto_local":
                thresh = min(tmat[i].max() for i in range(tmat.shape[0]))
                logg.debug(f"Using `threshold={thresh}` at `{key}`")
            elif isinstance(threshold, (int, float, np.number)):
                thresh = np.percentile(tmat.data, threshold)
                logg.debug(f"Using `threshold={thresh}` at `{key}`")

            tmat = sp.csr_matrix(tmat, dtype=tmat.dtype)
            tmat.data[tmat.data < thresh] = 0.0
            tmat.eliminate_zeros()
            couplings[key] = AnnData(tmat, obs=adata.obs, var=adata.var)

        return couplings if copy else None

    def _validate_couplings(
        self,
        couplings: Mapping[Key_t, AnnData],
    ) -> None:
        """Validate that transport maps conform to various invariants.

        Parameters
        ----------
        couplings
            Optimal transport couplings.

        Returns
        -------
        Possibly reordered transport maps.
        """

        def assert_same(expected: Sequence[Any], actual: Sequence[Any], msg: Optional[str] = None) -> None:
            try:
                pd.testing.assert_series_equal(
                    pd.Series(sorted(expected)),
                    pd.Series(sorted(actual)),
                    check_names=False,
                    check_flags=False,
                    check_index=True,
                )
            except AssertionError as e:
                raise IndexError(msg) from e

        for (src, tgt), coupling in couplings.items():
            src_obs = self.adata.obs_names[self.time == src]
            tgt_obs = self.adata.obs_names[self.time == tgt]
            assert_same(
                src_obs,
                coupling.obs_names,
                msg=f"Source observations for `{src, tgt}` don't match with `adata.obs_names`.",
            )
            assert_same(
                tgt_obs,
                coupling.var_names,
                msg=f"Source observations for `{src, tgt}` don't match with `adata.var_names`.",
            )

        for key, coupling in couplings.items():
            tmat = coupling.X
            zeros_mask = np.where(np.asarray(tmat.sum(1)).squeeze() == 0)[0]
            if np.any(zeros_mask):
                logg.warning(
                    f"Coupling at `{key}` contains `{len(zeros_mask)}` empty row(s), e.g., due to "
                    f"unbalancedness or 0s in the source marginals. Using uniform distribution"
                )
                tmat[zeros_mask] = np.ones((len(zeros_mask), tmat.shape[1]), dtype=tmat.dtype) / tmat.shape[1]

    def _coupling_to_adata(self, src: Any, tgt: Any, coupling: Coupling_t) -> AnnData:
        """Convert the coupling to :class:`~anndata.AnnData`."""
        if not isinstance(coupling, AnnData):
            coupling = AnnData(X=coupling)
            coupling.obs_names = np.asarray(self.adata.obs_names)[self.time == src]
            coupling.var_names = np.asarray(self.adata.obs_names)[self.time == tgt]

        if sp.issparse(coupling.X) and not sp.isspmatrix_csr(coupling.X):
            coupling.X = coupling.X.tocsr()

        return coupling

    def _get_default_coupling(
        self,
        couplings: Optional[Mapping[Key_t, Optional[Coupling_t]]],
        policy: Literal["sequential", "triu"],
    ) -> Tuple[Dict[Key_t, Optional[Coupling_t]], Any]:
        cats = self._time.cat.categories
        if policy == "sequential":
            if couplings is None:
                couplings = {(src, tgt): None for src, tgt in zip(cats[:-1], cats[1:])}
            reference = sorted(k for ks in couplings for k in ks)[-1]
        elif policy == "triu":
            if couplings is None:
                couplings = {(src, tgt): None for src, tgt in itertools.product(cats, cats) if src < tgt}
            reference = sorted(k for ks in couplings for k in ks)[-1]
        else:
            raise NotImplementedError(f"Handling `{policy}` is not yet implemented.")

        for keys in couplings:
            for key in keys:
                if key not in cats:
                    raise ValueError(f"Key `{key}` is not in `{sorted(cats)}`.")

        return dict(couplings), reference

    @property
    def time(self) -> pd.Series:
        """Experimental time."""
        return self._time

    @property
    def couplings(self) -> Optional[Dict[Key_t, AnnData]]:
        """Optimal transport couplings."""
        return self._couplings

    @property
    def obs(self) -> Optional[pd.DataFrame]:
        """Cell-level metadata."""
        return self._obs


def _compute_connectivity_tmat(
    adata: AnnData, density_normalize: bool = True, **kwargs: Any
) -> Union[np.ndarray, sp.spmatrix]:
    from cellrank.kernels import ConnectivityKernel

    sc.pp.neighbors(adata, **kwargs)
    return ConnectivityKernel(adata).compute_transition_matrix(density_normalize=density_normalize).transition_matrix
