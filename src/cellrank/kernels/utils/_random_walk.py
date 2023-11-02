import itertools
import pathlib
from typing import Any, List, Literal, Mapping, Optional, Sequence, Tuple, Union

import scvelo as scv

import numpy as np
import pandas as pd
import scipy.sparse as sp
from pandas.api.types import infer_dtype

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, to_hex

from anndata import AnnData

from cellrank import logging as logg
from cellrank._utils._docs import d
from cellrank._utils._parallelize import parallelize
from cellrank._utils._utils import save_fig
from cellrank.kernels._utils import _get_basis

__all__ = ["RandomWalk"]

Indices_t = Optional[Union[Sequence[str], Mapping[str, Union[str, Sequence[str], Tuple[float, float]]]]]


@d.dedent
class RandomWalk:
    """Class that simulates a random walk on a Markov chain.

    Parameters
    ----------
    %(adata)s
    transition_matrix
        Row-stochastic transition matrix.
    start_ixs
        Indices from which to uniformly sample the starting points. If :obj:`None`, use all points.
    stop_ixs
        Indices which when hit, the random walk is terminated.
    """

    def __init__(
        self,
        adata: AnnData,
        transition_matrix: Union[np.ndarray, sp.spmatrix],
        start_ixs: Optional[Sequence[int]] = None,
        stop_ixs: Optional[Sequence[int]] = None,
    ):
        if transition_matrix.ndim != 2 or (transition_matrix.shape[0] != transition_matrix.shape[1]):
            raise ValueError(f"Expected transition matrix to be a square matrix, found `{transition_matrix.ndim}`.")
        if transition_matrix.shape[0] != adata.n_obs:
            raise ValueError(
                f"Expected transition matrix to be of shape `{adata.n_obs, adata.n_obs}`,"
                f"found `{transition_matrix.shape}`."
            )
        if not np.allclose(transition_matrix.sum(1), 1.0):
            raise ValueError("Transition matrix is not row-stochastic.")

        self._adata = adata
        self._tmat = transition_matrix
        self._ixs = np.arange(self._tmat.shape[0])
        self._is_sparse = sp.issparse(self._tmat)

        start_ixs = self._normalize_ixs(start_ixs, kind="start")
        stop_ixs = self._normalize_ixs(stop_ixs, kind="stop")
        self._stop_ixs = set([] if stop_ixs is None or not len(stop_ixs) else stop_ixs)
        self._starting_dist = np.ones_like(self._ixs) if start_ixs is None else np.isin(self._ixs, start_ixs)
        _sum = np.sum(self._starting_dist)
        if _sum == 0:
            raise ValueError("No starting indices have been selected.")

        self._starting_dist = self._starting_dist.astype(transition_matrix.dtype) / _sum

    @d.get_sections(base="rw_sim", sections=["Parameters"])
    def simulate_one(
        self,
        max_iter: Union[int, float] = 0.25,
        seed: Optional[int] = None,
        successive_hits: int = 0,
    ) -> np.ndarray:
        """Simulate one random walk.

        Parameters
        ----------
        max_iter
            Maximum number of steps of a random walk. If a :class:`float`, it can be specified
            as a fraction of the number of cells.
        seed
            Random seed.
        successive_hits
            Number of successive hits in the ``stop_ixs`` required to stop prematurely.

        Returns
        -------
        Array of shape ``(max_iter + 1,)`` of states that have been visited.
        If ``stop_ixs`` was specified, the array may have a smaller shape.
        """
        max_iter = self._max_iter(max_iter)
        if successive_hits < 0:
            raise ValueError(f"Expected number of successive hits to be positive, found `{successive_hits}`.")

        rs = np.random.RandomState(seed)
        ix = rs.choice(self._ixs, p=self._starting_dist)
        sim, cnt = [ix], -1

        for _ in range(max_iter):
            ix = self._sample(ix, rs=rs)
            sim.append(ix)
            cnt = (cnt + 1) if self._should_stop(ix) else -1
            if cnt >= successive_hits:
                break

        return np.array(sim)

    def _simulate_many(
        self,
        sims: np.ndarray,
        max_iter: Union[int, float] = 0.25,
        seed: Optional[int] = None,
        successive_hits: int = 0,
        queue: Optional[Any] = None,
    ) -> List[np.ndarray]:
        res = []
        for s in sims:
            sim = self.simulate_one(
                max_iter=max_iter,
                seed=None if seed is None else seed + s,
                successive_hits=successive_hits,
            )
            res.append(sim)
            if queue is not None:
                queue.put(1)

        if queue is not None:
            queue.put(None)

        return res

    @d.dedent
    def simulate_many(
        self,
        n_sims: int,
        max_iter: Union[int, float] = 0.25,
        seed: Optional[int] = None,
        successive_hits: int = 0,
        n_jobs: Optional[int] = None,
        backend: str = "loky",
        show_progress_bar: bool = True,
    ) -> List[np.ndarray]:
        """Simulate many random walks.

        Parameters
        ----------
        n_sims
            Number of random walks to simulate.
        %(rw_sim.params)s
        %(parallel)s

        Returns
        -------
        List of arrays of shape ``(max_iter + 1,)`` of states that have been visited.
        If ``stop_ixs`` was specified, the arrays may have smaller shape.
        """
        if n_sims <= 0:
            raise ValueError(f"Expected number of simulations to be positive, found `{n_sims}`.")
        max_iter = self._max_iter(max_iter)
        start = logg.info(f"Simulating `{n_sims}` random walks of maximum length `{max_iter}`")

        simss = parallelize(
            self._simulate_many,
            collection=np.arange(n_sims),
            n_jobs=n_jobs,
            backend=backend,
            show_progress_bar=show_progress_bar,
            as_array=False,
            unit="sim",
        )(max_iter=max_iter, seed=seed, successive_hits=successive_hits)
        simss = list(itertools.chain.from_iterable(simss))

        logg.info("    Finish", time=start)

        return simss

    @d.dedent
    def plot(
        self,
        sims: List[np.ndarray],
        basis: str = "umap",
        cmap: Union[str, LinearSegmentedColormap] = "gnuplot",
        linewidth: float = 1.0,
        linealpha: float = 0.3,
        ixs_legend_loc: Optional[str] = None,
        figsize: Optional[Tuple[float, float]] = None,
        dpi: Optional[int] = None,
        save: Optional[Union[str, pathlib.Path]] = None,
        **kwargs: Any,
    ) -> None:
        """Plot simulated random walks.

        Parameters
        ----------
        sims
            The simulated random walks.
        basis
            Basis used for plotting.
        cmap
            Colormap for the random walks.
        linewidth
            Line width for the random walks.
        linealpha
            Line alpha.
        ixs_legend_loc
            Position of the legend describing start- and endpoints.
        %(plotting)s
        kwargs
            Keyword arguments for :func:`~scvelo.pl.scatter`.

        Returns
        -------
        %(just_plots)s
        """
        emb = _get_basis(self._adata, basis)
        if isinstance(cmap, str):
            cmap = plt.get_cmap(cmap)
        if not isinstance(cmap, LinearSegmentedColormap):
            if not hasattr(cmap, "colors"):
                raise AttributeError("Unable to create a colormap, `cmap` does not have attribute `colors`.")
            cmap = LinearSegmentedColormap.from_list(
                "random_walk",
                colors=cmap.colors,
                N=max(map(len, sims)),
            )

        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        scv.pl.scatter(self._adata, basis=basis, show=False, ax=ax, **kwargs)

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
            from scvelo.plotting.utils import default_size, plot_outline

            plot_outline(
                x=emb[ixs][:, 0],
                y=emb[ixs][:, 1],
                outline_color=("black", to_hex(cmap(float(abs(ix))))),
                kwargs={
                    "s": kwargs.get("size", default_size(self._adata)) * 1.1,
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

    def _normalize_ixs(self, ixs: Indices_t, *, kind: Literal["start", "stop"]) -> Optional[np.ndarray]:
        if ixs is None:
            return None

        if isinstance(ixs, dict):
            # fmt: off
            if len(ixs) != 1:
                raise ValueError(f"Expected to find only 1 cluster key, found `{len(ixs)}`.")
            key = next(iter(ixs.keys()))
            if key not in self._adata.obs:
                raise KeyError(f"Unable to find data in `adata.obs[{key!r}]`.")

            vals = self._adata.obs[key]
            if isinstance(vals.dtype, pd.CategoricalDtype):
                ixs = np.where(np.isin(vals, ixs[key]))[0]
            elif np.issubdtype(vals.dtype, np.number):
                if len(ixs[key]) != 2:
                    raise ValueError(f"Expected range to be of length `2`, found `{len(ixs[key])}`")
                minn, maxx = sorted(ixs[key])
                ixs = np.where((vals >= minn) & (vals <= maxx))[0]
            else:
                raise TypeError(f"Expected `adata.obs[{key!r}]` to be numeric or categorical, "
                                f"found `{infer_dtype(vals)}`.")
            # fmt: on
        elif isinstance(ixs, str):
            ixs = np.where(self._adata.obs_names == ixs)[0]
        elif isinstance(ixs[0], str):
            ixs = np.where(np.isin(self._adata.obs_names, ixs))[0]
        elif isinstance(ixs[0], bool):
            if len(ixs) != self._adata.n_obs:
                raise ValueError(
                    f"Expected `bool` {kind} indices of length" f"`{self._adata.n_obs}`, found `{len(ixs)}`."
                )
            ixs = np.where(ixs)[0]
        elif isinstance(ixs[0], int):
            ixs = list(set(ixs))
            if max(ixs) >= self._adata.n_obs:
                raise IndexError(max(ixs))
            if min(ixs) < -self._adata.n_obs:
                raise IndexError(min(ixs))
        else:
            raise TypeError(
                f"Expected {kind} indices to be either `dict` or a sequence of "
                f"`int`, `str`, `bool`, found `{type(ixs).__nam__}`."
            )

        if not len(ixs):
            raise ValueError(f"No {kind} indices have been selected.")

        return ixs

    def _should_stop(self, ix: int) -> bool:
        return ix in self._stop_ixs

    def _sample(self, ix: int, *, rs: np.random.RandomState) -> int:
        return rs.choice(
            self._ixs,
            p=self._tmat[ix].A.squeeze() if self._is_sparse else self._tmat[ix],
        )

    def _max_iter(self, max_iter: Union[int, float]) -> int:
        if isinstance(max_iter, float):
            max_iter = int(np.ceil(max_iter * len(self._ixs)))
        if max_iter <= 1:
            raise ValueError(f"Expected number of iterations to be > 1, found `{max_iter}`.")
        return max_iter
