import abc
import copy
import pathlib
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import scipy.sparse as sp

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from anndata import AnnData

from cellrank import logging as logg
from cellrank._utils._docs import d, inject_docs
from cellrank._utils._utils import _normalize, _read_graph_data, save_fig
from cellrank.kernels._utils import require_tmat
from cellrank.kernels.mixins import BidirectionalMixin, IOMixin, UnidirectionalMixin
from cellrank.kernels.utils import FlowPlotter, RandomWalk, TmatProjection
from cellrank.kernels.utils._random_walk import Indices_t

__all__ = ["Kernel", "UnidirectionalKernel", "BidirectionalKernel"]

Tmat_t = Union[np.ndarray, sp.spmatrix]


class KernelExpression(IOMixin, abc.ABC):
    def __init__(
        self,
        parent: Optional["KernelExpression"] = None,
        **kwargs: Any,
    ):
        super().__init__()
        self._parent = parent
        self._normalize = parent is None
        self._transition_matrix = None
        self._params: Dict[str, Any] = {}
        self._init_kwargs = kwargs  # for `_read_from_adata`

    def __init_subclass__(cls, **_: Any) -> None:
        super().__init_subclass__()

    @abc.abstractmethod
    def compute_transition_matrix(self, *args: Any, **kwargs: Any) -> "KernelExpression":
        """Compute transition matrix.

        Parameters
        ----------
        args
            Positional arguments.
        kwargs
            Keyword arguments.

        Returns
        -------
        Modifies :attr:`transition_matrix` and returns self.
        """

    @abc.abstractmethod
    def copy(self, *, deep: bool = False) -> "KernelExpression":
        """Return a copy of self. The :attr:`adata` object is not copied."""

    @property
    @abc.abstractmethod
    def adata(self) -> AnnData:
        """Annotated data object."""

    @adata.setter
    @abc.abstractmethod
    def adata(self, value: Optional[AnnData]) -> None:
        pass

    @property
    @abc.abstractmethod
    def kernels(self) -> Tuple["KernelExpression", ...]:
        """Underlying base kernels."""

    @property
    @abc.abstractmethod
    def shape(self) -> Tuple[int, int]:
        """``(n_cells, n_cells)``."""

    @property
    @abc.abstractmethod
    def backward(self) -> Optional[bool]:
        """Direction of the process."""

    @abc.abstractmethod
    def __getitem__(self, ix: int) -> "KernelExpression":
        pass

    @abc.abstractmethod
    def __len__(self) -> int:
        pass

    @d.get_full_description(base="plot_single_flow")
    @d.get_sections(base="plot_single_flow", sections=["Parameters", "Returns"])
    @d.dedent
    @require_tmat
    def plot_single_flow(
        self,
        cluster: str,
        cluster_key: str,
        time_key: str,
        clusters: Optional[Sequence[Any]] = None,
        time_points: Optional[Sequence[Union[int, float]]] = None,
        min_flow: float = 0,
        remove_empty_clusters: bool = True,
        ascending: Optional[bool] = False,
        legend_loc: Optional[str] = "upper right out",
        alpha: Optional[float] = 0.8,
        xticks_step_size: Optional[int] = 1,
        figsize: Optional[Tuple[float, float]] = None,
        dpi: Optional[int] = None,
        save: Optional[Union[str, pathlib.Path]] = None,
        show: bool = True,
    ) -> Optional[plt.Axes]:
        """Visualize outgoing flow from a cluster of cells :cite:`mittnenzweig:21`.

        Parameters
        ----------
        cluster
            Cluster for which to visualize outgoing flow.
        cluster_key
            Key in :attr:`~anndata.AnnData.obs` where clustering is stored.
        time_key
            Key in :attr:`~anndata.AnnData.obs` where experimental time is stored.
        clusters
            Visualize flow only for these clusters. If :obj:`None`, use all clusters.
        time_points
            Visualize flow only for these time points. If :obj:`None`, use all time points.
        %(flow.parameters)s
        %(plotting)s
        show
            If :obj:`False`, return :class:`~matplotlib.axes.Axes`.

        Returns
        -------
        The axes object, if ``show = False``.
        %(just_plots)s

        Notes
        -----
        This function is a Python re-implementation of the following
        `original R function <https://github.com/tanaylab/embflow/blob/main/scripts/generate_paper_figures/plot_vein.r>`_
        with some minor stylistic differences.
        This function will not recreate the results from :cite:`mittnenzweig:21`, because there, the *Metacell* model
        :cite:`baran:19` was used to compute the flow, whereas here the transition matrix is used.
        """  # noqa: E501
        fp = FlowPlotter(self.adata, self.transition_matrix, cluster_key, time_key)
        fp = fp.prepare(cluster, clusters, time_points)

        ax = fp.plot(
            min_flow=min_flow,
            remove_empty_clusters=remove_empty_clusters,
            ascending=ascending,
            alpha=alpha,
            xticks_step_size=xticks_step_size,
            legend_loc=legend_loc,
            figsize=figsize,
            dpi=dpi,
        )

        if save is not None:
            save_fig(ax.figure, save)

        if not show:
            return ax

    @d.dedent
    @require_tmat
    def plot_random_walks(
        self,
        n_sims: int = 100,
        max_iter: Union[int, float] = 0.25,
        seed: Optional[int] = None,
        successive_hits: int = 0,
        start_ixs: Indices_t = None,
        stop_ixs: Indices_t = None,
        basis: str = "umap",
        cmap: Union[str, LinearSegmentedColormap] = "gnuplot",
        linewidth: float = 1.0,
        linealpha: float = 0.3,
        ixs_legend_loc: Optional[str] = None,
        n_jobs: Optional[int] = None,
        backend: str = "loky",
        show_progress_bar: bool = True,
        figsize: Optional[Tuple[float, float]] = None,
        dpi: Optional[int] = None,
        save: Optional[Union[str, pathlib.Path]] = None,
        **kwargs: Any,
    ) -> None:
        """Plot random walks in an embedding.

        This method simulates random walks on the Markov chain defined though the corresponding transition matrix. The
        method is intended to give qualitative rather than quantitative insights into the transition matrix. Random
        walks are simulated by iteratively choosing the next cell based on the current cell's transition probabilities.

        Parameters
        ----------
        n_sims
            Number of random walks to simulate.
        %(rw_sim.parameters)s
        start_ixs
            Cells from which to sample the starting points. If :obj:`None`, use all cells.
            %(rw_ixs)s
            For example ``{'dpt_pseudotime': [0, 0.1]}`` means that starting points for random walks
            will be sampled uniformly from cells whose pseudotime is in :math:`[0, 0.1]`.
        stop_ixs
            Cells which when hit, the random walk is terminated. If :obj:`None`, terminate after ``max_iters``.
            %(rw_ixs)s
            For example ``{'clusters': ['Alpha', 'Beta']}`` and ``successive_hits = 3`` means that the random walk will
            stop prematurely after cells in the above specified clusters have been visited successively 3 times in
            a row.
        basis
            Basis in :attr:`~anndata.AnnData.obsm` to use as an embedding.
        cmap
            Colormap for the random walk lines.
        linewidth
            Width of the random walk lines.
        linealpha
            Alpha value of the random walk lines.
        ixs_legend_loc
            Legend location for the start/top indices.
        %(parallel)s
        %(plotting)s
        kwargs
            Keyword arguments for :func:`~scvelo.pl.scatter`.

        Returns
        -------
        %(just_plots)s
        For each random walk, the first/last cell is marked by the start/end colors of ``cmap``.
        """
        rw = RandomWalk(self.adata, self.transition_matrix, start_ixs=start_ixs, stop_ixs=stop_ixs)
        sims = rw.simulate_many(
            n_sims=n_sims,
            max_iter=max_iter,
            seed=seed,
            n_jobs=n_jobs,
            backend=backend,
            successive_hits=successive_hits,
            show_progress_bar=show_progress_bar,
        )

        rw.plot(
            sims,
            basis=basis,
            cmap=cmap,
            linewidth=linewidth,
            linealpha=linealpha,
            ixs_legend_loc=ixs_legend_loc,
            figsize=figsize,
            dpi=dpi,
            save=save,
            **kwargs,
        )

    @require_tmat
    def plot_projection(
        self,
        basis: str = "umap",
        key_added: Optional[str] = None,
        recompute: bool = False,
        stream: bool = True,
        connectivities: Optional[sp.spmatrix] = None,
        **kwargs: Any,
    ) -> None:
        """Plot :attr:`transition_matrix` as a stream or a grid plot.

        Parameters
        ----------
        basis
            Key in :attr:`~anndata.AnnData.obsm` containing the basis.
        key_added
            If not :obj:`None`, save the result to :attr:`adata.obsm['{key_added}'] <anndata.AnnData.obsm>`.
            Otherwise, save the result to ``'T_fwd_{basis}'`` or ``'T_bwd_{basis}'``, depending on the direction.
        recompute
            Whether to recompute the projection if it already exists.
        stream
            If :obj:`True`, use :func:`~scvelo.pl.velocity_embedding_stream`.
            Otherwise, use :func:`~scvelo.pl.velocity_embedding_grid`.
        connectivities
            Connectivity matrix to use for projection. If :obj:`None`, use ones from the underlying kernel, is possible.
        kwargs
            Keyword argument for the above-mentioned plotting function.

        Returns
        -------
        Nothing, just plots and modifies :attr:`~anndata.AnnData.obsm` with a key based on the ``key_added``.
        """
        proj = TmatProjection(self, basis=basis)
        proj.project(key_added=key_added, recompute=recompute, connectivities=connectivities)
        proj.plot(stream=stream, **kwargs)

    def __add__(self, other: "KernelExpression") -> "KernelExpression":
        return self.__radd__(other)

    def __radd__(self, other: "KernelExpression") -> "KernelExpression":
        def same_level_add(k1: "KernelExpression", k2: "KernelExpression") -> bool:
            if not (isinstance(k1, KernelAdd) and isinstance(k2, KernelMul)):
                return False

            for kexpr in k1:
                if not isinstance(kexpr, KernelMul):
                    return False
                if not kexpr._bin_consts:
                    return False
            return True

        if not isinstance(other, KernelExpression):
            return NotImplemented

        s = self * 1.0 if isinstance(self, Kernel) else self
        o = other * 1.0 if isinstance(other, Kernel) else other

        # (c1 * x + c2 * y + ...) + (c3 * z) => (c1 * x + c2 * y + c3 + ... + * z)
        if same_level_add(s, o):
            return KernelAdd(*tuple(s) + (o,))
        if same_level_add(o, s):
            return KernelAdd(*tuple(o) + (s,))

        # add virtual constant
        if not isinstance(s, KernelMul):
            s = s * 1.0
        if not isinstance(o, KernelMul):
            o = o * 1.0

        return KernelAdd(s, o)

    def __mul__(self, other: Union[float, int, "KernelExpression"]) -> "KernelExpression":
        return self.__rmul__(other)

    def __rmul__(self, other: Union[int, float, "KernelExpression"]) -> "KernelExpression":
        def same_level_mul(k1: "KernelExpression", k2: "KernelExpression") -> bool:
            return isinstance(k1, KernelMul) and isinstance(k2, Constant) and k1._bin_consts

        if isinstance(other, (int, float, np.integer, np.floating)):
            other = Constant(self.adata, other)

        if not isinstance(other, KernelExpression):
            return NotImplemented

        # fmt: off
        s = self if isinstance(self, (KernelMul, Constant)) else KernelMul(Constant(self.adata, 1.0), self)
        o = other if isinstance(other, (KernelMul, Constant)) else KernelMul(Constant(other.adata, 1.0), other)
        # fmt: on

        # at this point, only KernelMul and Constant is possible (two constants are not possible)
        # (c1 * k) * c2 => (c1 * c2) * k
        if same_level_mul(s, o):
            c, expr = s._split_const
            return KernelMul(c * o, expr)
        if same_level_mul(o, s):
            c, expr = o._split_const
            return KernelMul(c * s, expr)

        return KernelMul(s, o)

    @d.get_sections(base="write_to_adata", sections=["Parameters"])
    @inject_docs()  # gets rid of {{}} in %(write_to_adata)s
    @d.dedent
    @require_tmat
    def write_to_adata(self, key: Optional[str] = None, copy: bool = False) -> None:
        """Write the transition matrix and parameters used for computation to the underlying :attr:`adata` object.

        Parameters
        ----------
        key
            Key used when writing transition matrix to :attr:`adata`.
            If :obj:`None`, the key will be determined automatically.
        copy
            Whether to copy the :attr:`transition_matrix`.

        Returns
        -------
        %(write_to_adata)s
        """
        from cellrank._utils._key import Key

        if self.adata is None:
            raise ValueError("Underlying annotated data object is not set.")

        key = Key.uns.kernel(self.backward, key=key)
        # retain the embedding info
        self.adata.uns[f"{key}_params"] = {
            **self.adata.uns.get(f"{key}_params", {}),
            **{"params": self.params},
            **{"init": self._init_kwargs},
        }
        self.adata.obsp[key] = self.transition_matrix.copy() if copy else self.transition_matrix

    @property
    def transition_matrix(self) -> Union[np.ndarray, sp.csr_matrix]:
        """Row-normalized transition matrix."""
        if self._parent is None and self._transition_matrix is None:
            self.compute_transition_matrix()
        return self._transition_matrix

    @transition_matrix.setter
    def transition_matrix(self, matrix: Tmat_t) -> None:
        """Set the transition matrix.

        Parameters
        ----------
        matrix
            Transition matrix. The matrix is row-normalized if necessary.

        Returns
        -------
        Nothing, just updates the :attr:`transition_matrix` and optionally normalizes it.
        """
        # fmt: off
        if matrix.shape != self.shape:
            raise ValueError(
                f"Expected matrix to be of shape `{self.shape}`, found `{matrix.shape}`."
            )

        def should_norm(mat: Tmat_t) -> bool:
            return not np.all(np.isclose(np.asarray(mat.sum(1)).squeeze(), 1.0, rtol=1e-12))

        if sp.issparse(matrix) and not sp.isspmatrix_csr(matrix):
            matrix = sp.csr_matrix(matrix)
        matrix = matrix.astype(np.float64, copy=False)

        force_normalize = (self._parent is None or self._normalize) and should_norm(matrix)
        if force_normalize:
            if np.any((matrix.data if sp.issparse(matrix) else matrix) < 0):
                raise ValueError("Unable to normalize matrix with negative values.")
            matrix = _normalize(matrix)
            if should_norm(matrix):  # some rows are all 0s/contain invalid values
                n_inv = np.sum(~np.isclose(np.asarray(matrix.sum(1)).squeeze(), 1.0, rtol=1e-12))
                raise ValueError(f"Transition matrix is not row stochastic, {n_inv} rows do not sum to 1.")
        # fmt: on

        self._transition_matrix = matrix

    @property
    def params(self) -> Dict[str, Any]:
        """Parameters which are used to compute the transition matrix."""
        if len(self.kernels) == 1:
            return self._params
        return {f"{k!r}:{i}": k.params for i, k in enumerate(self.kernels)}

    def _reuse_cache(self, expected_params: Dict[str, Any], *, time: Optional[Any] = None) -> bool:
        # fmt: off
        try:
            if expected_params == self._params:
                assert self.transition_matrix is not None
                logg.debug("Using cached transition matrix")
                logg.info("    Finish", time=time)
                return True
            return False
        except AssertionError:
            logg.warning("Transition matrix does not exist for the given parameters. Recomputing")
            return False
        except Exception as e:  # noqa: BLE001
            # e.g. the dict is not comparable
            logg.warning(f"Expected and actually parameters are not comparable, reason `{e}`. Recomputing")
            expected_params = {}  # clear the params
            return False
        finally:
            self._params = expected_params
        # fmt: on

    def _get_boundary(self, source: str, target: str, cluster_key: str, graph_key: str = "distances") -> List[int]:
        """Identify source observations at boundary to target cluster.

        Parameters
        ----------
        source
            Name of source cluster.
        target
            Name of target cluster.
        cluster_key
            Key in :attr:`~anndata.AnnData.obs` to obtain cluster annotations.
        graph_key
            Name of graph representation to use from :attr:`~anndata.AnnData.obsp`.

        Returns
        -------
        List of observation IDs at boundary to target cluster.
        """
        source_obs_mask = self.adata.obs[cluster_key].isin([source] if isinstance(source, str) else source)
        target_obs_mask = self.adata.obs[cluster_key].isin([target] if isinstance(target, str) else target)

        source_ids = np.where(source_obs_mask)[0]
        boundary_ids = []

        graph = self.adata.obsp[graph_key]
        for source_id in source_ids:
            obs_mask = graph[source_id, :].toarray().squeeze().astype(bool)

            if (obs_mask & target_obs_mask).any():
                boundary_ids.append(source_id)

        return boundary_ids

    def _get_empirical_velocity_field(
        self, boundary_ids: List[int], target_obs_mask, rep: str, graph_key: str = "distances"
    ) -> np.ndarray:
        """Compute an emprical estimate of velocity field between two clusters.

        Parameters
        ----------
        boundary_ids
            List of observation IDs at boundary to target cluster.
        target_obs_mask
            Boolean indicator identifying relevant observations from target.
        graph_key
            Name of graph representation to use from :attr:`~anndata.AnnData.obsp`.

        Returns
        -------
        Empirical velocity estimate.
        """
        obs_ids = np.arange(0, self.adata.n_obs)
        graph = self.adata.obsp[graph_key]
        features = self.adata.obsm[rep]
        empirical_velo = np.empty(shape=(len(boundary_ids), features.shape[1]))

        for idx, boundary_id in enumerate(boundary_ids):
            row = graph[boundary_id, :].toarray().squeeze()
            obs_mask = row.astype(bool) & target_obs_mask
            neighbors = obs_ids[obs_mask]
            weights = row[obs_mask]

            empirical_velo[idx, :] = np.sum(
                weights.reshape(-1, 1) * (features[neighbors, :] - features[boundary_id, :]), axis=0
            )

        empirical_velo = np.array(empirical_velo)
        obs_mask = np.isnan(empirical_velo).any(axis=1)
        empirical_velo = empirical_velo[~obs_mask, :]

        return empirical_velo

    def _get_vector_field_estimate(self, rep: str) -> np.ndarray:
        """Compute estimate of vector field under one step of the transition matrix.

        Parameters
        ----------
        rep
            Key in :attr:`~anndata.AnnData.obsm` to use as data representation.

        Returns
        -------
        Vector field estimate based on kernel dynamics.
        """
        extrapolated_gex = self.transition_matrix @ self.adata.obsm[rep]
        return extrapolated_gex - self.adata.obsm[rep]

    # TODO: Add definition/reference to paper
    def cbc(self, source: str, target: str, cluster_key: str, rep: str, graph_key: str = "distances") -> np.ndarray:
        """Compute cross-boundary correctness score between source and target cluster.

        Parameters
        ----------
        source
            Name of the source cluster.
        target
            Name of the target cluster.
        cluster_key
            Key in :attr:`~anndata.AnnData.obs` to obtain cluster annotations.
        rep
            Key in :attr:`~anndata.AnnData.obsm` to use as data representation.
        graph_key
            Name of graph representation to use from :attr:`~anndata.AnnData.obsp`.

        Returns
        -------
        Cross-boundary correctness score for each observation.
        """

        def _pearsonr(x: np.ndarray, y: np.ndarray) -> np.ndarray:
            x_centered = x - np.mean(x, axis=1, keepdims=True)
            y_centered = y - np.mean(y, axis=1, keepdims=True)
            denom = np.linalg.norm(x_centered, axis=1) * np.linalg.norm(y_centered, axis=1)

            return np.sum(x_centered * y_centered, axis=1) / denom

        target_obs_mask = self.adata.obs[cluster_key].isin([target])
        boundary_ids = self._get_boundary(source=source, target=target, cluster_key=cluster_key, graph_key=graph_key)
        empirical_velo = self._get_empirical_velocity_field(
            boundary_ids=boundary_ids, target_obs_mask=target_obs_mask, rep=rep, graph_key=graph_key
        )
        estimated_velo = self._get_vector_field_estimate(rep=rep)[boundary_ids, :]

        return _pearsonr(x=estimated_velo, y=empirical_velo)


@d.dedent
class Kernel(KernelExpression, abc.ABC):
    """Base kernel class.

    Parameters
    ----------
    %(adata)s
    parent
        Parent kernel expression.
    kwargs
        Keyword arguments for the parent.
    """

    def __init__(self, adata: AnnData, parent: Optional[KernelExpression] = None, **kwargs: Any):
        super().__init__(parent=parent, **kwargs)
        self._adata = adata
        self._n_obs = adata.n_obs
        self._read_from_adata(**kwargs)

    def _read_from_adata(self, **kwargs: Any) -> None:
        pass

    @property
    def adata(self) -> AnnData:
        """Annotated data object."""
        return self._adata

    @adata.setter
    def adata(self, adata: Optional[AnnData]) -> None:
        if adata is None:
            self._adata = None
            return
        if not isinstance(adata, AnnData):
            raise TypeError(f"Expected `adata` to be of type `AnnData`, found `{type(adata).__name__}`.")
        shape = (adata.n_obs, adata.n_obs)
        if self.shape != shape:
            raise ValueError(
                f"Expected new `AnnData` object to have same shape as the previous `{self.shape}`, found `{shape}`."
            )
        self._adata = adata

    @classmethod
    @d.dedent
    def from_adata(
        cls,
        adata: AnnData,
        key: str,
        copy: bool = False,
    ) -> "Kernel":
        """Read the kernel saved using :meth:`write_to_adata`.

        Parameters
        ----------
        %(adata)s
        key
            Key in :attr:`~anndata.AnnData.obsp` where the transition matrix is stored.
            The parameters should be stored in :attr:`adata.uns['{key}_params'] <anndata.AnnData.uns>`.
        copy
            Whether to copy the transition matrix.

        Returns
        -------
        The kernel with explicitly initialized properties:

        - :attr:`transition_matrix` - the transition matrix.
        - :attr:`params` - parameters used for computation.
        """
        transition_matrix = _read_graph_data(adata, key=key)
        try:
            params = adata.uns[f"{key}_params"]["params"].copy()
            init_params = adata.uns[f"{key}_params"]["init"].copy()
        except KeyError as e:
            raise KeyError(f"Unable to kernel parameters, reason: `{e}`") from e

        if copy:
            transition_matrix = transition_matrix.copy()

        kernel = cls(adata, **init_params)
        kernel.transition_matrix = transition_matrix
        kernel._params = params

        return kernel

    def copy(self, *, deep: bool = False) -> "Kernel":
        """Return a copy of self.

        Parameters
        ----------
        deep
            Whether to use :func:`~copy.deepcopy`.

        Returns
        -------
        Copy of self.
        """
        with self._remove_adata:
            k = copy.deepcopy(self)
        k.adata = self.adata.copy() if deep else self.adata
        return k

    def _copy_ignore(self, *attrs: str) -> "Kernel":
        # prevent copying attributes that are not necessary, e.g. during inversion
        sentinel, attrs = object(), set(attrs)
        objects = [
            (attr, obj)
            for attr, obj in ((attr, getattr(self, attr, sentinel)) for attr in attrs)
            if obj is not sentinel
        ]
        try:
            for attr, _ in objects:
                setattr(self, attr, None)
            return self.copy(deep=False)
        finally:
            for attr, obj in objects:
                setattr(self, attr, obj)

    def _format_params(self) -> str:
        n, _ = self.shape
        params = ", ".join(f"{k}={round(v, 3) if isinstance(v, float) else v!r}" for k, v in self.params.items())
        return f"n={n}, {params}" if params else f"n={n}"

    @property
    def transition_matrix(self) -> Union[np.ndarray, sp.csr_matrix]:
        """Row-normalized transition matrix."""
        return self._transition_matrix

    @transition_matrix.setter
    def transition_matrix(self, matrix: Any) -> None:
        KernelExpression.transition_matrix.fset(self, matrix)

    @property
    def kernels(self) -> Tuple["KernelExpression", ...]:
        """Underlying base kernels."""
        return (self,)

    @property
    def shape(self) -> Tuple[int, int]:
        """``(n_cells, n_cells)``."""
        return self._n_obs, self._n_obs

    def __getitem__(self, ix: int) -> "Kernel":
        if ix != 0:
            raise IndexError(ix)
        return self

    def __len__(self) -> int:
        return 1

    def __repr__(self) -> str:
        params = self._format_params()
        prefix = "~" if self.backward and self._parent is None else ""
        return f"{prefix}{self.__class__.__name__}[{params}]"

    def __str__(self) -> str:
        prefix = "~" if self.backward and self._parent is None else ""
        return f"{prefix}{self.__class__.__name__}[n={self.shape[0]}]"


class UnidirectionalKernel(UnidirectionalMixin, Kernel, abc.ABC):
    pass


class BidirectionalKernel(BidirectionalMixin, Kernel, abc.ABC):
    pass


class Constant(UnidirectionalKernel):
    def __init__(self, adata: AnnData, value: Union[int, float]):
        super().__init__(adata)
        self.transition_matrix = value

    def compute_transition_matrix(self, value: Union[int, float]) -> "Constant":
        self.transition_matrix = value
        return self

    @property
    def transition_matrix(self) -> Union[int, float]:
        return self._transition_matrix

    @transition_matrix.setter
    def transition_matrix(self, value: Union[int, float]) -> None:
        if not isinstance(value, (int, float, np.integer, np.floating)):
            raise TypeError(f"Value must be a `float` or `int`, found `{type(value).__name__}`.")
        if value <= 0:
            raise ValueError(f"Expected the scalar to be positive, found `{value}`.")

        self._transition_matrix = value
        self._params = {"value": value}

    def copy(self, *, deep: bool = False) -> "Constant":
        return Constant(self.adata, self.transition_matrix)

    # fmt: off
    def __radd__(self, other: Union[int, float, "KernelExpression"]) -> "Constant":
        if isinstance(other, (int, float, np.integer, np.floating)):
            other = Constant(self.adata, other)
        if isinstance(other, Constant):
            if self.shape != other.shape:
                raise ValueError(f"Expected kernel shape to be `{self.shape}`, found `{other.shape}`.")
            return Constant(self.adata, self.transition_matrix + other.transition_matrix)

        return super().__radd__(other)

    def __rmul__(self, other: Union[int, float, "KernelExpression"]) -> "Constant":
        if isinstance(other, (int, float, np.integer, np.floating)):
            other = Constant(self.adata, other)
        if isinstance(other, Constant):
            if self.shape != other.shape:
                raise ValueError(f"Expected kernel shape to be `{self.shape}`, found `{other.shape}`.")
            return Constant(self.adata, self.transition_matrix * other.transition_matrix)

        return super().__rmul__(other)
    # fmt: on

    def __repr__(self) -> str:
        return repr(round(self.transition_matrix, 3))

    def __str__(self) -> str:
        return str(round(self.transition_matrix, 3))


class NaryKernelExpression(BidirectionalMixin, KernelExpression):
    def __init__(self, *kexprs: KernelExpression, parent: Optional[KernelExpression] = None):
        super().__init__(parent=parent)
        self._validate(kexprs)

    @abc.abstractmethod
    def _combine_transition_matrices(self, t1: Tmat_t, t2: Tmat_t) -> Tmat_t:
        pass

    @property
    @abc.abstractmethod
    def _initial_value(self) -> Union[int, float, Tmat_t]:
        pass

    @property
    @abc.abstractmethod
    def _combiner(self) -> str:
        pass

    def compute_transition_matrix(self) -> "NaryKernelExpression":
        for kexpr in self:
            if kexpr.transition_matrix is None:
                if isinstance(kexpr, Kernel):
                    raise RuntimeError(
                        f"`{kexpr}` is uninitialized. Compute its "
                        f"transition matrix first as `.compute_transition_matrix()`."
                    )
                kexpr.compute_transition_matrix()

        tmat = self._initial_value
        for kexpr in self:
            tmat = self._combine_transition_matrices(tmat, kexpr.transition_matrix)
        self.transition_matrix = tmat

        return self

    def _validate(self, kexprs: Sequence[KernelExpression]) -> None:
        if not len(kexprs):
            raise ValueError("No kernels to combined.")

        shapes = {kexpr.shape for kexpr in kexprs}
        if len(shapes) > 1:
            raise ValueError(f"Expected all kernels to have the same shapes, found `{sorted(shapes)}`.")

        directions = {kexpr.backward for kexpr in kexprs}
        if True in directions and False in directions:
            raise ValueError("Unable to combine both forward and backward kernels.")

        self._backward = True if True in directions else False if False in directions else None
        self._kexprs = kexprs
        for kexpr in self:
            kexpr._parent = self

    @property
    def adata(self) -> AnnData:
        """Annotated data object."""
        return self[0].adata

    @adata.setter
    def adata(self, adata: Optional[AnnData]) -> None:
        # allow resetting (use for temp. pickling without adata)
        for kexpr in self:
            kexpr.adata = adata

    @property
    def kernels(self) -> Tuple["KernelExpression", ...]:
        """Underlying unique basic kernels."""
        kernels = []
        for kexpr in self:
            if isinstance(kexpr, Kernel) and not isinstance(kexpr, Constant):
                kernels.append(kexpr)
            elif isinstance(kexpr, NaryKernelExpression):  # recurse
                kernels.extend(kexpr.kernels)

        # return only unique kernels
        return tuple(set(kernels))

    def copy(self, *, deep: bool = False) -> "KernelExpression":
        kexprs = (k.copy(deep=deep) for k in self)
        return type(self)(*kexprs, parent=self._parent)

    @property
    def shape(self) -> Tuple[int, int]:
        """``(n_cells, n_cells)``."""
        # all kernels have the same shape
        return self[0].shape

    def _format(self, formatter: Callable[[KernelExpression], str]) -> str:
        return (
            f"{'~' if self.backward and self._parent is None else ''}("
            + f" {self._combiner} ".join(formatter(kexpr) for kexpr in self)
            + ")"
        )

    def __getitem__(self, ix: int) -> "KernelExpression":
        return self._kexprs[ix]

    def __len__(self) -> int:
        return len(self._kexprs)

    def __invert__(self) -> "KernelExpression":
        if self.backward is None:
            return self.copy()
        kexprs = tuple(~k if isinstance(k, BidirectionalMixin) else k.copy() for k in self)
        kexpr = type(self)(*kexprs, parent=self._parent)
        kexpr._transition_matrix = None
        return kexpr

    def __repr__(self) -> str:
        return self._format(repr)

    def __str__(self) -> str:
        return self._format(str)


class KernelAdd(NaryKernelExpression):
    def compute_transition_matrix(self) -> "KernelAdd":
        self._maybe_recalculate_constants()
        return super().compute_transition_matrix()

    def _combine_transition_matrices(self, t1: Tmat_t, t2: Tmat_t) -> Tmat_t:
        return t1 + t2

    def _maybe_recalculate_constants(self) -> None:
        """Normalize constants to sum to 1."""
        constants = self._bin_consts
        if constants:
            total = sum((c.transition_matrix for c in constants), 0.0)
            for c in constants:
                c.transition_matrix = c.transition_matrix / total
            for kexpr in self:  # don't normalize  (c * x)
                kexpr._normalize = False

    @property
    def _bin_consts(self) -> List[KernelExpression]:
        """Return constant expressions for each binary multiplication children with at least 1 constant."""
        return [c for k in self if isinstance(k, KernelMul) for c in k._bin_consts]

    @property
    def _combiner(self) -> str:
        return "+"

    @property
    def _initial_value(self) -> Union[int, float, Tmat_t]:
        return 0.0


class KernelMul(NaryKernelExpression):
    def _combine_transition_matrices(self, t1: Tmat_t, t2: Tmat_t) -> Tmat_t:
        if sp.issparse(t1):
            return t1.multiply(t2)
        if sp.issparse(t2):
            return t2.multiply(t1)
        return t1 * t2

    def _format(self, formatter: Callable[[KernelExpression], str]) -> str:
        fmt = super()._format(formatter)
        if fmt[0] == "~" or not self._bin_consts:
            return fmt
        if fmt.startswith("("):
            fmt = fmt[1:]
        if fmt.endswith(")"):
            fmt = fmt[:-1]
        return fmt

    @property
    def _bin_consts(self) -> List[KernelExpression]:
        """Return all constants if this expression contains only 2 subexpressions."""
        if len(self) != 2:
            return []

        return [k for k in self if isinstance(k, Constant)]

    @property
    def _split_const(self) -> Tuple[Optional[Constant], Optional[KernelExpression]]:
        """Return a constant and the other expression, iff this expression is of length 2 and contains a constant."""
        if not self._bin_consts:
            return None, None
        k1, k2 = self
        if isinstance(k1, Constant):
            return k1, k2
        return k2, k1

    @property
    def _combiner(self) -> str:
        return "*"

    @property
    def _initial_value(self) -> Union[int, float, Tmat_t]:
        return 1.0
