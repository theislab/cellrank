"""Kernel module."""
from typing import (
    Any,
    Dict,
    List,
    Type,
    Tuple,
    Union,
    Callable,
    Iterable,
    Optional,
    Sequence,
)

import warnings
from abc import ABC, abstractmethod
from copy import copy
from pathlib import Path
from functools import reduce

import scvelo as scv
from anndata import AnnData
from cellrank import logging as logg
from cellrank.ul._docs import d, inject_docs
from cellrank.tl._utils import (
    save_fig,
    _connected,
    _normalize,
    _symmetric,
    _get_neighs,
    _irreducible,
)
from scvelo.plotting.utils import default_size, plot_outline
from cellrank.tl._mixins._io import IOMixin
from cellrank.tl.kernels._utils import _get_basis, _filter_kwargs
from cellrank.tl.kernels._tmat_flow import FlowPlotter
from cellrank.tl.kernels._random_walk import RandomWalk

import numpy as np
from scipy.sparse import spdiags, issparse, spmatrix, csr_matrix, isspmatrix_csr
from pandas.api.types import infer_dtype, is_numeric_dtype
from pandas.core.dtypes.common import is_categorical_dtype

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_hex
from matplotlib.collections import LineCollection

_ERROR_DIRECTION_MSG = "Can only combine kernels that have the same direction."
_ERROR_EMPTY_CACHE_MSG = (
    "Fatal error: tried to used cached values, but the cache was empty."
)
_ERROR_CONF_ADAPT = (
    "Confidence adaptive operator is only supported for kernels, found type `{!r}`."
)
_ERROR_VAR_NOT_FOUND = "Variances not found in kernel `{!r}`."

_LOG_USING_CACHE = "Using cached transition matrix"

_RTOL = 1e-12
_n_dec = 2
_dtype = np.float64
_cond_num_tolerance = 1e-15
Indices_t = Optional[
    Union[Sequence[str], Dict[str, Union[str, Sequence[str], Tuple[float, float]]]]
]


class KernelExpression(IOMixin, ABC):
    """Base class for all kernels and kernel expressions."""

    def __init__(
        self,
        op_name: Optional[str] = None,
        backward: bool = False,
        compute_cond_num: bool = False,
    ):
        self._op_name = op_name
        self._transition_matrix = None
        self._backward = backward
        self._compute_cond_num = compute_cond_num
        self._cond_num = None
        self._params = {}
        self._normalize = True
        self._parent = None

    def __init_subclass__(cls, **kwargs: Any):
        super().__init_subclass__()

    @property
    def condition_number(self) -> Optional[int]:
        """Condition number of the transition matrix."""
        return self._cond_num

    @property
    def transition_matrix(self) -> Union[np.ndarray, spmatrix]:
        """
        Return row-normalized transition matrix.

        If not present, it is computed iff all underlying kernels have been initialized.
        """

        if self._parent is None and self._transition_matrix is None:
            self.compute_transition_matrix()

        return self._transition_matrix

    @property
    def backward(self) -> bool:
        """Direction of the process."""
        return self._backward

    @property
    @abstractmethod
    def adata(self) -> AnnData:
        """Annotated data object."""

    @adata.setter
    @abstractmethod
    def adata(self, value: AnnData) -> None:
        pass

    @property
    @abstractmethod
    def shape(self) -> Tuple[int, int]:
        """`(n_cells, n_cells)`."""

    @property
    def params(self) -> Dict[str, Any]:
        """Parameters which are used to compute the transition matrix."""
        if len(self.kernels) == 1:
            return self._params
        # we need some identifier
        return {f"{repr(k)}:{i}": k.params for i, k in enumerate(self.kernels)}

    def _format_params(self):
        return ", ".join(
            f"{k}={round(v, _n_dec) if isinstance(v, float) else v}"
            for k, v in self.params.items()
        )

    @transition_matrix.setter
    def transition_matrix(self, value: Union[np.ndarray, spmatrix]) -> None:
        """
        Set a new value of the transition matrix.

        Parameters
        ----------
        value
            The new transition matrix. If the expression has no parent, the matrix is normalized, if needed.

        Returns
        -------
        None
            Nothing, just updates the :attr:`transition_matrix` and optionally normalizes it.
        """
        should_norm = ~np.isclose(value.sum(1), 1.0, rtol=_RTOL).all()

        if self._parent is None:
            self._transition_matrix = _normalize(value) if should_norm else value
        else:
            # it's AND, not OR, because of combinations
            self._transition_matrix = (
                _normalize(value) if self._normalize and should_norm else value
            )

    @abstractmethod
    def compute_transition_matrix(
        self, *args: Any, **kwargs: Any
    ) -> "KernelExpression":
        """
        Compute a transition matrix.

        Parameters
        ----------
        args
            Positional arguments.
        kwargs
            Keyword arguments.

        Returns
        -------
        :class:`cellrank.tl.kernels.KernelExpression`
            Self.
        """

    @d.get_sections(base="write_to_adata", sections=["Parameters"])
    @inject_docs()  # get rid of {{}}
    @d.dedent
    def write_to_adata(self, key: Optional[str] = None) -> None:
        """
        Write the transition matrix and parameters used for computation to the underlying :attr:`adata` object.

        Parameters
        ----------
        key
            Key used when writing transition matrix to :attr:`adata`. If `None`, determine the key automatically.

        Returns
        -------
        None
            %(write_to_adata)s
        """
        from cellrank._key import Key

        if self._transition_matrix is None:
            raise ValueError(
                "Compute transition matrix first as `.compute_transition_matrix()`."
            )

        key = Key.uns.kernel(self.backward, key=key)
        # retain the embedding info
        self.adata.uns[f"{key}_params"] = {
            **self.adata.uns.get(f"{key}_params", {}),
            **{"params": self.params},
        }
        self.adata.obsp[key] = self.transition_matrix

    @abstractmethod
    def copy(self) -> "KernelExpression":
        """Return a copy of itself. Note that the underlying :attr:`adata` object is not copied."""

    def _maybe_compute_cond_num(self) -> None:
        """Optionally compute condition number."""
        if self._compute_cond_num and self._cond_num is None:
            logg.debug("Computing condition number")
            self._cond_num = np.linalg.cond(
                self._transition_matrix.toarray()
                if issparse(self._transition_matrix)
                else self._transition_matrix
            )
            if self._cond_num > _cond_num_tolerance:
                logg.warning(
                    f"Transition matrix may be ill-conditioned, its condition number is `{self._cond_num:.2e}`"
                )
            else:
                logg.info(f"Condition number is `{self._cond_num:.2e}`")

    def _reuse_cache(
        self, expected_params: Dict[str, Any], *, time: Optional[Any] = None
    ) -> bool:
        if expected_params == self._params:
            assert self.transition_matrix is not None, _ERROR_EMPTY_CACHE_MSG
            logg.debug(_LOG_USING_CACHE)
            logg.info("    Finish", time=time)
            return True

        self._params = expected_params
        return False

    @abstractmethod
    def _get_kernels(self) -> Iterable["Kernel"]:
        pass

    @property
    def kernels(self) -> List["Kernel"]:
        """Get the kernels of the kernel expression, except for constants."""
        return list(set(self._get_kernels()))

    def compute_projection(
        self,
        basis: str = "umap",
        key_added: Optional[str] = None,
        copy: bool = False,
    ) -> Optional[np.ndarray]:
        """
        Compute a projection of the transition matrix in the embedding.

        Projections can only be calculated for kNN based kernels. The projected matrix
        can be then visualized as::

            scvelo.pl.velocity_embedding(adata, vkey='T_fwd', basis='umap')

        Parameters
        ----------
        basis
            Basis in :attr:`anndata.AnnData.obsm` for which to compute the projection.
        key_added
            If not `None` and ``copy = False``, save the result to :attr:`anndata.AnnData.obsm` ``['{key_added}']``.
            Otherwise, save the result to `'T_fwd_{basis}'` or `T_bwd_{basis}`, depending on the direction.
        copy
            Whether to return the projection or modify :attr:`adata` inplace.

        Returns
        -------
        If ``copy=True``, the projection array of shape `(n_cells, n_components)`.
        Otherwise, it modifies :attr:`anndata.AnnData.obsm` with a key based on ``key_added``.
        """
        # modified from: https://github.com/theislab/scvelo/blob/master/scvelo/tools/velocity_embedding.py
        from cellrank._key import Key
        from scvelo.tools.velocity_embedding import quiver_autoscale

        if self._transition_matrix is None:
            raise RuntimeError(
                "Compute transition matrix first as `.compute_transition_matrix()`."
            )

        for kernel in self.kernels:
            if kernel._conn is None:
                raise AttributeError(
                    f"{kernel!r} is not a kNN based kernel. The embedding projection "
                    "only works for kNN based kernels."
                )

        start = logg.info(f"Projecting transition matrix onto `{basis}`")
        emb = _get_basis(self.adata, basis)
        T_emb = np.empty_like(emb)

        conn = self.kernels[0]._conn

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for row_id, row in enumerate(self.transition_matrix):
                conn_idxs = conn[row_id, :].indices

                dX = emb[conn_idxs] - emb[row_id, None]

                if np.any(np.isnan(dX)):
                    T_emb[row_id, :] = np.nan
                else:
                    probs = row[:, conn_idxs]
                    if issparse(probs):
                        probs = probs.A.squeeze()

                    dX /= np.linalg.norm(dX, axis=1)[:, None]
                    dX = np.nan_to_num(dX)
                    T_emb[row_id, :] = probs.dot(dX) - dX.sum(0) / dX.shape[0]

        T_emb /= 3 * quiver_autoscale(np.nan_to_num(emb), T_emb)

        if copy:
            return T_emb

        key = Key.uns.kernel(self.backward, key=key_added)
        ukey = f"{key}_params"

        embs = self.adata.uns.get(ukey, {}).get("embeddings", [])
        if basis not in embs:
            embs = list(embs) + [basis]
            self.adata.uns[ukey] = self.adata.uns.get(ukey, {})
            self.adata.uns[ukey]["embeddings"] = embs

        key = key + "_" + basis
        logg.info(
            f"Adding `adata.obsm[{key!r}]`\n    Finish",
            time=start,
        )
        self.adata.obsm[key] = T_emb

    @d.dedent
    def plot_random_walks(
        self,
        n_sims: int,
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
        save: Optional[Union[str, Path]] = None,
        **kwargs: Any,
    ) -> None:
        """
        Plot random walks in an embedding.

        This method simulates random walks on the Markov chain defined though the corresponding transition matrix. The
        method is intended to give qualitative rather than quantitative insights into the transition matrix. Random
        walks are simulated by iteratively choosing the next cell based on the current cell's transition probabilities.

        Parameters
        ----------
        n_sims
            Number of random walks to simulate.
        %(rw_sim.parameters)s
        start_ixs
            Cells from which to sample the starting points. If `None`, use all cells.
            %(rw_ixs)s
            For example ``{'dpt_pseudotime': [0, 0.1]}`` means that starting points for random walks
            will be sampled uniformly from cells whose pseudotime is in `[0, 0.1]`.
        stop_ixs
            Cells which when hit, the random walk is terminated. If `None`, terminate after ``max_iters``.
            %(rw_ixs)s
            For example ``{'clusters': ['Alpha', 'Beta']}`` and ``successive_hits = 3`` means that the random walk will
            stop prematurely after cells in the above specified clusters have been visited successively 3 times in a
            row.
        basis
            Basis in :attr:`anndata.AnnData.obsm` to use as an embedding.
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
            Keyword arguments for :func:`scvelo.pl.scatter`.

        Returns
        -------
        %(just_plots)s
        For each random walk, the first/last cell is marked by the start/end colors of ``cmap``.
        """

        def create_ixs(ixs: Indices_t, *, kind: str) -> Optional[np.ndarray]:
            if ixs is None:
                return None
            if isinstance(ixs, dict):
                # fmt: off
                if len(ixs) != 1:
                    raise ValueError(f"Expected to find only 1 cluster key, found `{len(ixs)}`.")
                key = next(iter(ixs.keys()))
                if key not in self.adata.obs:
                    raise KeyError(f"Unable to find data in `adata.obs[{key!r}]`.")

                vals = self.adata.obs[key]
                if is_categorical_dtype(vals):
                    ixs = np.where(np.isin(vals, ixs[key]))[0]
                elif is_numeric_dtype(vals):
                    if len(ixs[key]) != 2:
                        raise ValueError(f"Expected range to be of length `2`, found `{len(ixs[key])}`")
                    minn, maxx = sorted(ixs[key])
                    ixs = np.where((vals >= minn) & (vals <= maxx))[0]
                else:
                    raise TypeError(f"Expected `adata.obs[{key!r}]` to be numeric or categorical, "
                                    f"found `{infer_dtype(vals)}`.")
                # fmt: on
            elif isinstance(ixs, str):
                ixs = np.where(self.adata.obs_names == ixs)[0]
            else:
                ixs = np.where(np.isin(self.adata.obs_names, ixs))[0]

            if not len(ixs):
                logg.warning(f"No {kind} indices have been selected, using `None`")
                return None

            return ixs

        if self._transition_matrix is None:
            raise RuntimeError(
                "Compute transition matrix first as `.compute_transition_matrix()`."
            )
        emb = _get_basis(self.adata, basis)

        if isinstance(cmap, str):
            cmap = plt.get_cmap(cmap)
        if not isinstance(cmap, LinearSegmentedColormap):
            if not hasattr(cmap, "colors"):
                raise AttributeError(
                    "Unable to create a colormap, `cmap` does not have attribute `colors`."
                )
            cmap = LinearSegmentedColormap.from_list(
                "random_walk", colors=cmap.colors, N=max_iter
            )

        start_ixs = create_ixs(start_ixs, kind="start")
        stop_ixs = create_ixs(stop_ixs, kind="stop")
        rw = RandomWalk(self.transition_matrix, start_ixs=start_ixs, stop_ixs=stop_ixs)
        sims = rw.simulate_many(
            n_sims=n_sims,
            max_iter=max_iter,
            seed=seed,
            n_jobs=n_jobs,
            backend=backend,
            successive_hits=successive_hits,
            show_progress_bar=show_progress_bar,
        )

        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        scv.pl.scatter(self.adata, basis=basis, show=False, ax=ax, **kwargs)

        logg.info("Plotting random walks")
        for sim in sims:
            x = emb[sim][:, 0]
            y = emb[sim][:, 1]
            points = np.array([x, y]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            n_seg = len(segments)

            lc = LineCollection(
                segments,
                linewidths=linewidth,
                colors=[cmap(float(i) / n_seg) for i in range(n_seg)],
                alpha=linealpha,
                zorder=2,
            )
            ax.add_collection(lc)

        for ix in [0, -1]:
            ixs = [sim[ix] for sim in sims]
            plot_outline(
                x=emb[ixs][:, 0],
                y=emb[ixs][:, 1],
                outline_color=("black", to_hex(cmap(float(abs(ix))))),
                kwargs={
                    "s": kwargs.get("s", default_size(self.adata)) * 1.1,
                    "alpha": 0.9,
                },
                ax=ax,
                zorder=4,
            )

        if ixs_legend_loc not in (None, "none"):
            from cellrank.pl._utils import _position_legend

            h1 = ax.scatter([], [], color=cmap(0.0), label="start")
            h2 = ax.scatter([], [], color=cmap(1.0), label="stop")
            legend = ax.get_legend()
            if legend is not None:
                ax.add_artist(legend)
            _position_legend(ax, legend_loc=ixs_legend_loc, handles=[h1, h2])

        if save is not None:
            save_fig(fig, save)

    @d.get_full_description(base="plot_single_flow")
    @d.get_sections(base="plot_single_flow", sections=["Parameters", "Returns"])
    @d.dedent
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
        save: Optional[Union[str, Path]] = None,
        show: bool = True,
    ) -> Optional[plt.Axes]:
        """
        Visualize outgoing flow from a cluster of cells :cite:`mittnenzweig:21`.

        Parameters
        ----------
        cluster
            Cluster for which to visualize outgoing flow.
        cluster_key
            Key in :attr:`anndata.AnnData.obs` where clustering is stored.
        time_key
            Key in :attr:`anndata.AnnData.obs` where experimental time is stored.
        clusters
            Visualize flow only for these clusters. If `None`, use all clusters.
        time_points
            Visualize flow only for these time points. If `None`, use all time points.
        %(flow.parameters)s
        %(plotting)s
        show
            If `False`, return :class:`matplotlib.pyplot.Axes`.

        Returns
        -------
        The axes object, if ``show = False``.
        %(just_plots)s

        Notes
        -----
        This function is a Python reimplementation of the following
        `original R function <https://github.com/tanaylab/embflow/blob/main/scripts/generate_paper_figures/plot_vein.r>`_
        with some minor stylistic differences.
        This function will not recreate the results from :cite:`mittnenzweig:21`, because there, the Metacell model
        :cite:`baran:19` was used to compute the flow, whereas here the transition matrix is used.
        """  # noqa: E501
        if self._transition_matrix is None:
            raise RuntimeError(
                "Compute transition matrix first as `.compute_transition_matrix()`."
            )

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

    def __xor__(self, other: "KernelExpression") -> "KernelExpression":
        return self.__rxor__(other)

    def __rxor__(self, other: "KernelExpression") -> "KernelExpression":
        def convert(obj):
            if _is_adaptive_type(obj):
                if obj._mat_scaler is None:
                    raise ValueError(_ERROR_VAR_NOT_FOUND.format(obj))
                return KernelMul(
                    [
                        ConstantMatrix(
                            obj.adata, 1, obj._mat_scaler, backward=obj.backward
                        ),
                        obj,
                    ]
                )
            if _is_adaptive_type(_is_bin_mult(obj, return_constant=False)):
                e, c = _get_expr_and_constant(obj)
                if e._mat_scaler is None:
                    raise ValueError(_ERROR_VAR_NOT_FOUND.format(e))
                return KernelMul(
                    [ConstantMatrix(e.adata, c, e._mat_scaler, backward=e.backward), e]
                )

            return obj

        if (
            not _is_adaptive_type(self)
            and not isinstance(self, KernelAdaptiveAdd)
            and not _is_adaptive_type(_is_bin_mult(self, return_constant=False))
        ):
            raise TypeError(_ERROR_CONF_ADAPT.format(self.__class__.__name__))
        if (
            not _is_adaptive_type(other)
            and not isinstance(other, KernelAdaptiveAdd)
            and not _is_adaptive_type(_is_bin_mult(other, return_constant=False))
        ):
            raise TypeError(_ERROR_CONF_ADAPT.format(other.__class__.__name__))

        s = convert(self)
        o = convert(other)

        if isinstance(s, KernelAdaptiveAdd):
            exprs = list(s) + (list(o) if isinstance(o, KernelAdaptiveAdd) else [o])
            # (c1 * x + c2 * y + ...) + (c3 * z) => (c1 * x + c2 * y + c3 + ... + * z)
            if all(_is_bin_mult(k, ConstantMatrix) for k in exprs):
                return KernelAdaptiveAdd(exprs)

        # same but reverse
        if isinstance(o, KernelAdaptiveAdd):
            exprs = (list(s) if isinstance(s, KernelAdaptiveAdd) else [s]) + list(o)
            if all(_is_bin_mult(k, ConstantMatrix) for k in exprs):
                return KernelAdaptiveAdd(exprs)

        return KernelAdaptiveAdd([s, o])

    def __add__(self, other: "KernelExpression") -> "KernelExpression":
        return self.__radd__(other)

    def __radd__(self, other: "KernelExpression") -> "KernelExpression":
        if not isinstance(other, KernelExpression):
            raise TypeError(
                f"Expected type `KernelExpression`, found `{other.__class__.__name__}`."
            )

        s = self * 1 if isinstance(self, Kernel) else self
        o = other * 1 if isinstance(other, Kernel) else other

        if isinstance(s, KernelSimpleAdd):
            exprs = list(s) + [o]
            # (c1 * x + c2 * y + ...) + (c3 * z) => (c1 * x + c2 * y + c3 + ... + * z)
            if all(_is_bin_mult(k) for k in exprs):
                return KernelSimpleAdd(exprs)

        # same but reverse
        if isinstance(o, KernelSimpleAdd):
            exprs = [s] + list(o)
            if all(_is_bin_mult(k) for k in exprs):
                return KernelSimpleAdd(exprs)

        # (c1 + c2) => c3
        if isinstance(s, Constant) and isinstance(o, Constant):
            assert s.backward == o.backward, _ERROR_DIRECTION_MSG
            return Constant(
                s.adata, s.transition_matrix + o.transition_matrix, backward=s.backward
            )

        ss = s * 1 if not isinstance(s, KernelMul) else s
        oo = o * 1 if not isinstance(o, KernelMul) else o

        return KernelSimpleAdd([ss, oo])

    def __rmul__(
        self, other: Union[int, float, "KernelExpression"]
    ) -> "KernelExpression":
        return self.__mul__(other)

    def __mul__(
        self, other: Union[float, int, "KernelExpression"]
    ) -> "KernelExpression":

        if isinstance(other, (int, float)):
            other = Constant(self.adata, other, backward=self.backward)

        if not isinstance(other, KernelExpression):
            raise TypeError(
                f"Expected type `KernelExpression`, found `{other.__class__.__name__}`."
            )

        # (c1 * c2) => c3
        if isinstance(self, Constant) and isinstance(
            other, Constant
        ):  # small optimization
            assert self.backward == other.backward, _ERROR_DIRECTION_MSG
            return Constant(
                self.adata,
                self.transition_matrix * other.transition_matrix,
                backward=self.backward,
            )

        s = (
            self
            if isinstance(self, (KernelMul, Constant))
            else KernelMul([Constant(self.adata, 1, backward=self.backward), self])
        )
        o = (
            other
            if isinstance(other, (KernelMul, Constant))
            else KernelMul([Constant(other.adata, 1, backward=other.backward), other])
        )

        cs, co = _is_bin_mult(s), _is_bin_mult(o)
        if cs and isinstance(o, Constant):
            cs._transition_matrix *= o.transition_matrix
            return s

        if co and isinstance(s, Constant):
            co._transition_matrix *= s.transition_matrix
            return o

        return KernelMul([s, o])

    def __invert__(self: "KernelExpression") -> "KernelExpression":
        # mustn't return a copy because transition matrix
        self._transition_matrix = None
        self._params = {}
        self._backward = not self.backward

        return self

    def __copy__(self) -> "KernelExpression":
        return self.copy()

    def __deepcopy__(self, memodict={}) -> "KernelExpression":  # noqa
        res = self.copy()
        if self._parent is None:
            res.adata = self.adata.copy()
        memodict[id(self)] = res

        return res


class UnaryKernelExpression(KernelExpression, ABC):
    """Base class for unary kernel expressions, such as kernels or constants."""

    def __init__(
        self,
        adata,
        backward: bool = False,
        op_name: Optional[str] = None,
        compute_cond_num: bool = False,
        **kwargs,
    ):
        super().__init__(op_name, backward=backward, compute_cond_num=compute_cond_num)
        assert (
            op_name is None
        ), "Unary kernel does not support any kind operation associated with it."
        self._adata = adata
        self._n_obs = adata.n_obs

    @property
    @d.dedent
    def adata(self) -> AnnData:
        """
        Annotated data object.

        Returns
        -------
        %(adata_ret)s
        """
        return self._adata

    @property
    def shape(self) -> Tuple[int, int]:
        """`(n_cells, n_cells)`."""
        return self._n_obs, self._n_obs

    @adata.setter
    def adata(self, adata: Optional[AnnData]) -> None:
        if adata is None:
            self._adata = None
            return
        if not isinstance(adata, AnnData):
            raise TypeError(
                f"Expected argument of type `anndata.AnnData`, found `{type(adata).__name__!r}`."
            )
        shape = (adata.n_obs, adata.n_obs)
        if self.shape != shape:
            raise ValueError(
                f"Expected the new object to have same shape as previous object `{self.shape}`, "
                f"found `{shape}`."
            )
        self._adata = adata

    def __repr__(self):
        return f"{'~' if self.backward and self._parent is None else ''}<{self.__class__.__name__}>"

    def __str__(self):
        params_fmt = self._format_params()
        if params_fmt:
            return (
                f"{'~' if self.backward and self._parent is None else ''}"
                f"<{self.__class__.__name__}[{params_fmt}]>"
            )
        return repr(self)


class NaryKernelExpression(KernelExpression, ABC):
    """Base class for n-ary kernel expressions."""

    def __init__(self, kexprs: List[KernelExpression], op_name: Optional[str] = None):
        assert len(kexprs), "No kernel expressions specified."

        backward = kexprs[0].backward
        assert all(k.backward == backward for k in kexprs), _ERROR_DIRECTION_MSG

        # use OR instead of AND
        super().__init__(
            op_name,
            backward=backward,
            compute_cond_num=any(k._compute_cond_num for k in kexprs),
        )

        # copies of constants are necessary because of the recalculation
        self._kexprs = [copy(k) if isinstance(k, Constant) else k for k in kexprs]

        for kexprs in self._kexprs:
            kexprs._parent = self

    @property
    def shape(self) -> Tuple[int, int]:
        """`(n_cells, n_cells)`."""
        return self.kernels[0].shape

    def _maybe_recalculate_constants(self, type_: Type):
        if type_ == Constant:
            accessor = "transition_matrix"
        elif type_ == ConstantMatrix:
            accessor = "_value"
        else:
            raise RuntimeError(
                f"Unable to determine accessor for type `{type(type_).__name__}`."
            )

        constants = [_is_bin_mult(k, type_) for k in self]
        if all(c is not None for c in constants):
            assert all(
                isinstance(getattr(c, accessor), (int, float)) for c in constants
            )
            total = sum(getattr(c, accessor) for c in constants) + 0.0
            for c in constants:
                c._recalculate(getattr(c, accessor) / total)

            for kexpr in self._kexprs:  # don't normalize  (c * x)
                kexpr._normalize = False

    def _get_kernels(self) -> Iterable["Kernel"]:
        for k in self:
            if isinstance(k, Kernel) and not isinstance(k, (Constant, ConstantMatrix)):
                yield k
            elif isinstance(k, NaryKernelExpression):
                yield from k._get_kernels()

    @property
    @d.dedent
    def adata(self) -> AnnData:
        """
        Annotated data object.

        Returns
        -------
        %(adata_ret)s
        """
        # we can do this because Constant requires adata as well
        return self._kexprs[0].adata

    @adata.setter
    def adata(self, adata: Optional[AnnData]) -> None:
        for kexpr in self._kexprs:
            kexpr.adata = adata

    def __invert__(self) -> "NaryKernelExpression":
        super().__invert__()
        self._kexprs = [~kexpr for kexpr in self]
        return self

    def __len__(self) -> int:
        return len(self._kexprs)

    def __getitem__(self, item) -> "KernelExpression":
        return self._kexprs[item]

    def __repr__(self) -> str:
        return (
            f"{'~' if self.backward and self._parent is None else ''}("
            + f" {self._op_name} ".join(repr(kexpr) for kexpr in self._kexprs)
            + ")"
        )

    def __str__(self) -> str:
        return (
            f"{'~' if self.backward and self._parent is None else ''}("
            + f" {self._op_name} ".join(str(kexpr) for kexpr in self._kexprs)
            + ")"
        )


@d.dedent
class Kernel(UnaryKernelExpression, ABC):
    """
    A base class from which all kernels are derived.

    These kernels read from a given AnnData object, usually the KNN graph and additional variables, to compute a
    weighted, directed graph. Every kernel object has a direction. The kernels defined in the derived classes are not
    strictly kernels in the mathematical sense because they often only take one input argument - however, they build
    on other functions which have computed a similarity based on two input arguments. The role of the kernels defined
    here is to add directionality to these symmetric similarity relations or to transform them.

    Parameters
    ----------
    %(adata)s
    %(backward)s
    %(cond_num)s
    check_connectivity
        Check whether the underlying KNN graph is connected.
    kwargs
        Keyword arguments which can specify key to be read from :attr:`adata` object.
    """

    def __init__(
        self,
        adata: AnnData,
        backward: bool = False,
        compute_cond_num: bool = False,
        check_connectivity: bool = False,
        **kwargs: Any,
    ):
        super().__init__(
            adata, backward, op_name=None, compute_cond_num=compute_cond_num, **kwargs
        )
        self._read_from_adata(check_connectivity=check_connectivity, **kwargs)

    # TODO: move to a mixin class
    def _read_from_adata(
        self,
        conn_key: Optional[str] = "connectivities",
        read_conn: bool = True,
        **kwargs: Any,
    ) -> None:
        """
        Import the base-KNN graph and optionally check for symmetry and connectivity.

        Parameters
        ----------
        conn_key
            Key in :attr:`anndata.AnnData.uns` where connectivities are stored.
        read_conn
            Whether to read connectivities or set them to `None`. Useful when not exposing density normalization or
            when KNN connectivities are not used to compute the transition matrix.
        kwargs
            Additional keyword arguments.

        Returns
        -------
        Nothing, just sets :attr:`_conn`.
        """
        if not read_conn:
            self._conn = None
            return

        self._conn = _get_neighs(
            self.adata, mode="connectivities", key=conn_key
        ).astype(_dtype)

        check_connectivity = kwargs.pop("check_connectivity", False)
        if check_connectivity:
            start = logg.debug("Checking the KNN graph for connectedness")
            if not _connected(self._conn):
                logg.warning("KNN graph is not connected", time=start)
            else:
                logg.debug("KNN graph is connected", time=start)

        start = logg.debug("Checking the KNN graph for symmetry")
        if not _symmetric(self._conn):
            logg.warning("KNN graph is not symmetric", time=start)
        else:
            logg.debug("KNN graph is symmetric", time=start)

    def _density_normalize(
        self, other: Union[np.ndarray, spmatrix]
    ) -> Union[np.ndarray, spmatrix]:
        """
        Density normalization by the underlying KNN graph.

        Parameters
        ----------
        other
            Matrix to normalize.

        Returns
        -------
        :class:`np.ndarray` or :class:`scipy.sparse.spmatrix`
            Density normalized transition matrix.

        Raises
        ------
        ValueError
            If KNN connectivities are not set.
        """
        if self._conn is None:
            raise ValueError(
                "Unable to density normalize the transition matrix "
                "because KNN connectivities are not set."
            )

        logg.debug("Density-normalizing the transition matrix")

        q = np.asarray(self._conn.sum(axis=0))
        Q = spdiags(1.0 / q, 0, other.shape[0], other.shape[0])

        return Q @ other @ Q

    def _get_kernels(self) -> Iterable["Kernel"]:
        yield self

    def _compute_transition_matrix(
        self,
        matrix: Union[np.ndarray, spmatrix],
        density_normalize: bool = True,
        check_irreducibility: bool = False,
    ):
        if matrix.shape[0] != matrix.shape[1]:
            raise ValueError(f"Expected a square matrix, found `{matrix.shape}`.")
        if matrix.shape[0] != self.adata.n_obs:
            raise ValueError(
                f"Expected matrix to be of shape `{(self.adata.n_obs, self.adata.n_obs)}`, "
                f"found `{matrix.shape}`."
            )

        matrix = matrix.astype(_dtype)
        if issparse(matrix) and not isspmatrix_csr(matrix):
            matrix = csr_matrix(matrix)

        # density correction based on node degrees in the KNN graph
        if density_normalize:
            matrix = self._density_normalize(matrix)

        # check for zero-rows
        problematic_indices = np.where(np.array(matrix.sum(1)).flatten() == 0)[0]
        if len(problematic_indices):
            logg.warning(
                f"Detected `{len(problematic_indices)}` absorbing states in the transition matrix. "
                f"This matrix won't be irreducible"
            )
            matrix[problematic_indices, problematic_indices] = 1.0

        if check_irreducibility:
            _irreducible(matrix)

        # setting this property automatically row-normalizes
        self.transition_matrix = matrix
        self._maybe_compute_cond_num()


@d.dedent
class Constant(Kernel):
    """
    Kernel representing a multiplication by a constant number.

    Parameters
    ----------
    %(adata)s
    value
        Constant value by which to multiply. Must be a positive number.
    %(backward)s
    """

    def __init__(
        self, adata: AnnData, value: Union[int, float], backward: bool = False
    ):
        super().__init__(adata, backward=backward)
        if not isinstance(value, (int, float)):
            raise TypeError(
                f"Value must be a `float` or `int`, found `{type(value).__name__}`."
            )
        if value <= 0:
            raise ValueError(f"Expected the constant to be positive, found `{value}`.")
        self._recalculate(value)

    def _recalculate(self, value) -> None:
        self._transition_matrix = value
        self._params = {"value": value}

    def _read_from_adata(self, **kwargs: Any) -> None:
        pass

    def compute_transition_matrix(self, *args, **kwargs) -> "Constant":
        """Return self."""
        return self

    @d.dedent
    def copy(self) -> "Constant":
        """%(copy)s"""  # noqa: D400, D401
        return Constant(self.adata, self.transition_matrix, self.backward)

    def __invert__(self) -> "Constant":
        # do not call parent's invert, since it removes the transition matrix
        self._backward = not self.backward
        return self

    def __repr__(self) -> str:
        return repr(round(self.transition_matrix, _n_dec))

    def __str__(self) -> str:
        return str(round(self.transition_matrix, _n_dec))


class ConstantMatrix(Kernel):
    """Kernel representing multiplication by a constant matrix."""

    def __init__(
        self,
        adata: AnnData,
        value: Union[int, float],
        variances: Union[np.ndarray, spmatrix],
        backward: bool = False,
    ):
        super().__init__(adata, backward=backward)
        conn_shape = self._conn.shape
        if variances.shape != conn_shape:
            raise ValueError(
                f"Expected variances of shape `{conn_shape}`, found `{variances.shape}`."
            )
        if not isinstance(value, (int, float)):
            raise TypeError(
                f"Value must be on `float` or `int`, found `{type(value).__name__!r}`."
            )
        if value <= 0:
            raise ValueError(f"Expected the constant to be positive, found `{value}`.")

        self._value = value
        self._mat_scaler = (
            csr_matrix(variances) if not issparse(variances) else variances
        )
        self._recalculate(value)

    @d.dedent
    def copy(self) -> "ConstantMatrix":
        """%(copy)s"""  # noqa: D400, D401
        return ConstantMatrix(
            self.adata, self._value, copy(self._mat_scaler), self.backward
        )

    def _recalculate(self, value) -> None:
        self._value = value
        self._params = {"value": value}
        self._transition_matrix = value * self._mat_scaler

    def compute_transition_matrix(self, *args, **kwargs) -> "ConstantMatrix":
        """Return self."""
        return self

    def __invert__(self) -> "ConstantMatrix":
        # do not call parent's invert, since it removes the transition matrix
        self._backward = not self.backward
        return self

    def __repr__(self) -> str:
        return repr(round(self._value, _n_dec))

    def __str__(self) -> str:
        return str(round(self._value, _n_dec))


class SimpleNaryExpression(NaryKernelExpression):
    """Base class for n-ary operations."""

    def __init__(self, kexprs: List[KernelExpression], op_name: str, fn: Callable):
        super().__init__(kexprs, op_name=op_name)
        self._fn = fn

    def compute_transition_matrix(self, *args, **kwargs) -> "SimpleNaryExpression":
        """Compute and combine the transition matrices."""
        # must be done before, because the underlying expression don't have to be normed
        if isinstance(self, KernelSimpleAdd):
            self._maybe_recalculate_constants(Constant)
        elif isinstance(self, KernelAdaptiveAdd):
            self._maybe_recalculate_constants(ConstantMatrix)

        for kexpr in self:
            if kexpr._transition_matrix is None:
                if isinstance(kexpr, Kernel):
                    raise RuntimeError(
                        f"Kernel `{kexpr}` is uninitialized. "
                        f"Compute its transition matrix first as `.compute_transition_matrix()`."
                    )
                kexpr.compute_transition_matrix()
            elif isinstance(kexpr, Kernel):
                logg.debug(_LOG_USING_CACHE)

        self.transition_matrix = csr_matrix(
            self._fn([kexpr.transition_matrix for kexpr in self])
        )

        # only the top level expression and kernels will have condition number computed
        if self._parent is None:
            self._maybe_compute_cond_num()

        return self

    @d.dedent
    def copy(self) -> "SimpleNaryExpression":
        """%(copy)s"""  # noqa: D400, D401
        constructor = type(self)
        kwargs = {"op_name": self._op_name, "fn": self._fn}

        # preserve the type information so that combination can properly work
        # we test for this
        sne = constructor(
            [copy(k) for k in self], **_filter_kwargs(constructor, **kwargs)
        )

        # TODO: copying's a bit buggy - the parent stays the same
        # we could disallow it for inner expressions
        sne._transition_matrix = copy(self._transition_matrix)
        sne._condition_number = self._cond_num
        sne._normalize = self._normalize

        return sne


class KernelAdd(SimpleNaryExpression):
    """Base class that represents the addition of :class:`KernelExpression`."""

    def __init__(self, kexprs: List[KernelExpression], op_name: str):
        super().__init__(kexprs, op_name=op_name, fn=Reductor(np.add, 0))


class KernelSimpleAdd(KernelAdd):
    """Addition of :class:`KernelExpression`."""

    def __init__(self, kexprs: List[KernelExpression]):
        super().__init__(kexprs, op_name="+")


class KernelAdaptiveAdd(KernelAdd):
    """Adaptive addition of :class:`KernelExpression`."""

    def __init__(self, kexprs: List[KernelExpression]):
        super().__init__(kexprs, op_name="^")


class KernelMul(SimpleNaryExpression):
    """Multiplication of :class:`KernelExpression`."""

    def __init__(self, kexprs: List[KernelExpression]):
        super().__init__(kexprs, op_name="*", fn=Reductor(np.multiply, 1))


class Reductor:
    """
    Class that reduces an iterable.

    We don't define a decorator because of pickling.

    Parameters
    ----------
    func
        Reduction function.
    initial
        Initial value for :func:`reduce`.
    """

    def __init__(self, func: Callable, initial: Union[int, float]):
        self._func = func
        self._initial = initial

    def __call__(self, seq: Iterable):  # noqa
        return reduce(self._func, seq, self._initial)


def _get_expr_and_constant(k: KernelMul) -> Tuple[KernelExpression, Union[int, float]]:
    """
    Get the value of a constant in binary multiplication.

    Parameters
    ----------
    k
        Binary multiplication involving a constant and a kernel.

    Returns
    -------
    :class:`KernelExpression`, int or float
        The expression which is being multiplied and the value of the constant.
    """

    if not isinstance(k, KernelMul):
        raise TypeError(
            f"Expected expression to be of type `KernelMul`, found `{type(k).__name__}`."
        )
    if len(k) != 2:
        raise ValueError(
            f"Expected expression to be binary, found `{len(k)}` subexpressions."
        )
    e1, e2 = k[0], k[1]

    if isinstance(e1, Constant):
        return e2, e1.transition_matrix
    elif isinstance(e2, Constant):
        return e1, e2.transition_matrix
    else:
        raise ValueError(
            "Expected one of the subexpressions to be `Constant`, found "
            f"`{type(e1).__name__}` and `{type(e2).__name__}`."
        )


def _is_bin_mult(
    k: KernelExpression,
    const_type: Union[Type, Tuple[Type]] = Constant,
    return_constant: bool = True,
) -> Optional[KernelExpression]:
    """
    Check if an expression is a binary multiplication.

    Parameters
    ----------
    k
        Kernel expression to check.
    const_type
        Type of the constant.
    return_constant
        Whether to return the constant in the expression or the multiplied expression.

    Returns
    -------
    None
        If the expression is not a binary multiplication.
    :class:`cellrank.tl.kernels.KernelExpression`
        Depending on ``return_constant``, it either returns the constant multiplier or the expression being multiplied.
    """

    if not isinstance(k, KernelMul):
        return None

    if len(k) != 2:
        return None

    lhs, rhs = k[0], k[1]

    if isinstance(lhs, const_type):
        return lhs if return_constant else rhs
    if isinstance(rhs, const_type):
        return rhs if return_constant else lhs

    return None


def _is_adaptive_type(k: KernelExpression) -> bool:
    return isinstance(k, Kernel) and not isinstance(k, (Constant, ConstantMatrix))
