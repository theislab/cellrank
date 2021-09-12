from typing import Any, List, Tuple, Union, Mapping, Optional, Sequence

from functools import lru_cache
from dataclasses import dataclass
from statsmodels.nonparametric.smoothers_lowess import lowess

from anndata import AnnData
from cellrank import logging as logg
from cellrank.ul._docs import d
from cellrank.tl._utils import _unique_order_preserving
from cellrank.tl._colors import _create_categorical_colors
from cellrank.tl.kernels._utils import _ensure_numeric_ordered

import numpy as np
import pandas as pd
from scipy.stats import logistic
from scipy.sparse import issparse, spmatrix
from pandas.api.types import infer_dtype
from scipy.interpolate import interp1d
from pandas.core.dtypes.common import is_categorical_dtype

import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.collections import PolyCollection

Numeric_t = Union[float, int]


@dataclass(frozen=True)
class Point:  # noqa: D101
    x: float
    xt: float


@d.dedent
class FlowPlotter:
    """
    Class that plots outgoing flow for a specific cluster :cite:`mittnenzweig:21`.

    It should be able to recreate (to a high degree) figures such as Fig. 4a in the above mentioned paper.

    Parameters
    ----------
    %(adata)s
    tmat
        Matrix of shape ``(adata.n_obs, adata.n_obs)``.
    cluster_key
        Key in :attr:`adata` ``.obs`` where clustering is stored.
    time_key
        Key in :attr:`adata` ``.obs`` where experimental time is stored.
    """

    TIME_KEY = "time"

    def __init__(
        self,
        adata: AnnData,
        tmat: Union[np.ndarray, spmatrix],
        cluster_key: str,
        time_key: str,
    ):
        self._adata = adata
        self._tmat = tmat
        self._ckey = cluster_key
        self._tkey = time_key

        self._cluster: Optional[str] = None
        self._clusters: Optional[Sequence[Any]] = None

        self._flow: Optional[pd.DataFrame] = None
        self._cmat: Optional[pd.DataFrame] = None

        if self._ckey not in self._adata.obs:
            raise KeyError(f"Unable to find clusters in `adata.obs[{self._ckey!r}]`.")
        if not is_categorical_dtype(self._adata.obs[self._ckey]):
            raise TypeError(
                f"Expected `adata.obs[{self._ckey!r}]` to be categorical, "
                f"found `{infer_dtype(self._adata.obs[self._ckey])}`."
            )
        self._adata.obs[self._tkey] = _ensure_numeric_ordered(self._adata, self._tkey)

    def prepare(
        self,
        cluster: str,
        clusters: Optional[Sequence[Any]] = None,
        time_points: Optional[Sequence[Numeric_t]] = None,
    ) -> "FlowPlotter":
        """
        Prepare itself for plotting by computing flow and contingency matrix.

        Parameters
        ----------
        cluster
            Source cluster for flow calculation.
        clusters
            Target clusters for flow calculation. If `None`, use all clusters.
        time_points
            Restrict flow calculation only to these time points. If `None`, use all time points.

        Returns
        -------
        Returns self and modifies internal internal attributes.
        """
        if clusters is None:
            self._clusters = self.clusters.cat.categories
        else:
            clusters = _unique_order_preserving([cluster] + list(clusters))
            mask = self.clusters.isin(clusters).values

            self._adata = self._adata[mask]
            if not self._adata.n_obs:
                raise ValueError("No valid clusters have been selected.")
            self._tmat = self._tmat[mask, :][:, mask]
            self._clusters = [c for c in clusters if c in self.clusters.cat.categories]

        if cluster not in self._clusters:
            raise ValueError(f"Invalid source cluster `{cluster!r}`.")

        if len(self._clusters) < 2:
            raise ValueError(
                f"Expected at least `2` clusters, found `{len(clusters)}`."
            )

        if time_points is not None:
            time_points = _unique_order_preserving(time_points)
            if len(time_points) < 2:
                raise ValueError(
                    f"Expected at least `2` time points, found `{len(time_points)}`."
                )

            mask = self.time.isin(time_points)

            self._adata = self._adata[mask]
            if not self._adata.n_obs:
                raise ValueError("No valid time points have been selected.")
            self._tmat = self._tmat[mask, :][:, mask]

        time_points = list(
            zip(self.time.cat.categories[:-1], self.time.cat.categories[1:])
        )

        logg.info(
            f"Computing flow from `{cluster}` into `{len(self._clusters) - 1}` cluster(s) "
            f"in `{len(time_points)}` time points"
        )
        self._cluster = cluster
        self._cmat = self.compute_contingency_matrix()
        self._flow = self.compute_flow(time_points, cluster)

        return self

    def compute_flow(
        self,
        time_points: Sequence[Tuple[Numeric_t, Numeric_t]],
        cluster: Optional[str] = None,
    ) -> pd.DataFrame:
        """
        Compute outgoing flow.

        Parameters
        ----------
        time_points
            Time point pair for which to calculate the flow.
        cluster
            Cluster for which to calculate the outgoing flow. If `None`, calculate the flow for all clusters.

        Returns
        -------
        Dataframe of shape ``(n_time_points, n_clusters)`` if ``cluster != None`` or
        a dataframe of shape ``(n_time_points * n_clusters, n_clusters)`` otherwise.
        The dataframe's index is a multi-index and the 1st level corresponds to time, the 2nd level to source clusters.
        """

        def default_helper(t1: Numeric_t, t2: Numeric_t) -> pd.DataFrame:
            subset, row_cls, col_cls = self._get_time_subset(t1, t2)

            df = pd.DataFrame(subset.A if issparse(subset) else subset)
            df = df.groupby(row_cls).sum().T.groupby(col_cls).sum().T

            res = pd.DataFrame(np.zeros((n, n)), index=categories, columns=categories)
            res.loc[df.index, df.columns] = df
            res.fillna(0, inplace=True)

            return res

        def cluster_helper(t1: Numeric_t, t2: Numeric_t) -> pd.DataFrame:
            subset, row_cls, col_cls = self._get_time_subset(t1, t2, cluster=cluster)

            df = pd.DataFrame(subset.A if issparse(subset) else subset).sum(0)
            df = df.groupby(col_cls).sum()
            df = pd.DataFrame([df], index=[cluster], columns=df.index)

            res = pd.DataFrame(np.zeros((1, n)), index=[cluster], columns=categories)
            res.loc[df.index, df.columns] = df
            res.fillna(0, inplace=True)

            return res

        categories = self.clusters.cat.categories
        n = len(categories)
        callback = cluster_helper if cluster is not None else default_helper
        flows, times = [], []

        for t1, t2 in time_points:
            flow = callback(t1, t2)
            times.extend([t1] * len(flow))
            flows.append(flow)

        flow = pd.concat(flows)
        flow.set_index([times, flow.index], inplace=True)
        flow /= flow.sum(1).values[:, None]
        flow.fillna(0, inplace=True)

        return flow

    def compute_contingency_matrix(self) -> pd.DataFrame:
        """Row-normalized contingency matrix of shape ``(n_clusters, n_time_points)``."""
        cmat = pd.crosstab(self.clusters, self.time)
        return (cmat / cmat.sum(0).values[None, :]).fillna(0)

    @d.get_sections(base="flow", sections=["Parameters"])
    def plot(
        self,
        min_flow: float = 0,
        remove_empty_clusters: bool = True,
        ascending: Optional[bool] = False,
        alpha: float = 0.8,
        xticks_step_size: Optional[int] = 1,
        legend_loc: Optional[str] = "upper right out",
        figsize: Optional[Tuple[float, float]] = None,
        dpi: Optional[int] = None,
    ) -> plt.Axes:
        """
        Plot outgoing flow.

        Parameters
        ----------
        min_flow
            Only show flow edges with flow greater than this value. Flow values are always in `[0, 1]`.
        remove_empty_clusters
            Whether to remove clusters with no incoming flow edges.
        ascending
            Whether to sort the cluster by ascending or descending incoming flow.
            If `None`, use the order as in defined by ``clusters``.
        alpha
            Alpha value for cell proportions.
        xticks_step_size
            Show only every n-th ticks on x-axis. If `None`, don't show any ticks.
        legend_loc
            Position of the legend. If `None`, do not show the legend.

        Returns
        -------
        The axes object.
        """
        if self._flow is None or self._cmat is None:
            raise RuntimeError(
                "Compute flow and contingency matrix first as `.prepare()`."
            )

        flow, cmat = self._flow, self._cmat
        try:
            if remove_empty_clusters:
                self._remove_min_clusters(min_flow)
            logg.info(
                f"Plotting flow from `{self._cluster}` into `{len(self._flow.columns) - 1}` cluster(s) "
                f"in `{len(self._cmat.columns) - 1}` time points"
            )
            return self._plot(
                self._rename_times(),
                ascending=ascending,
                min_flow=min_flow,
                alpha=alpha,
                xticks_step_size=xticks_step_size,
                legend_loc=legend_loc,
                figsize=figsize,
                dpi=dpi,
            )
        finally:
            self._flow = flow
            self._cmat = cmat

    def _get_time_subset(
        self, t1: Numeric_t, t2: Numeric_t, cluster: Optional[str] = None
    ) -> Tuple[Union[np.ndarray, spmatrix], pd.Series, pd.Series]:
        if cluster is None:
            row_ixs = np.where(self.time == t1)[0]
        else:
            row_ixs = np.where((self.time == t1) & (self.clusters == cluster))[0]

        col_ixs = np.where(self.time == t2)[0]
        row_cls = self.clusters.values[row_ixs]
        col_cls = self.clusters.values[col_ixs]

        return self._tmat[row_ixs, :][:, col_ixs], row_cls, col_cls

    def _remove_min_clusters(self, min_flow: float) -> None:
        logg.debug("Removing clusters with no incoming flow edges")
        columns = (self._flow.loc[(slice(None), self._cluster), :] > min_flow).any()
        columns = columns[columns].index
        if not len(columns):
            raise ValueError(
                "After removing clusters with no incoming flow edges, none remain."
            )
        self._flow = self._flow[columns]

    def _rename_times(self) -> Sequence[Numeric_t]:
        # make sure we have enough horizontal space to draw the flow (i.e. time points are at least 1 unit apart)
        old_times = self._cmat.columns
        tmp = np.array(old_times)
        tmp = (tmp - tmp.min()) / (tmp.max() - tmp.min())
        tmp /= np.min(tmp[1:] - tmp[:-1])
        time_mapper = dict(zip(old_times, tmp))
        self._flow.index = pd.MultiIndex.from_tuples(
            [(time_mapper[t], c) for t, c in self._flow.index]
        )
        self._cmat.columns = tmp
        return old_times

    def _order_clusters(
        self, cluster: str, ascending: Optional[bool] = False
    ) -> Tuple[List[Any], List[Any]]:
        if ascending is not None:
            tmp = [[], []]
            total_flow = (
                self._flow.loc[(slice(None), cluster), :]
                .sum()
                .sort_values(ascending=ascending)
            )
            for i, c in enumerate(c for c in total_flow.index if c != cluster):
                tmp[i % 2].append(c)
            return tmp[0][::-1], tmp[1]

        clusters = [c for c in self._clusters if c != cluster]
        return clusters[: len(clusters) // 2], clusters[len(clusters) // 2 :]

    def _calculate_y_offsets(
        self, clusters: Sequence[Any], delta: float = 0.2
    ) -> Mapping[Any, float]:
        offset = [0]
        for i in range(1, len(clusters)):
            offset.append(
                offset[-1]
                + delta
                + np.max(self._cmat.loc[clusters[i]] + self._cmat.loc[clusters[i - 1]])
            )
        return dict(zip(clusters, offset))

    def _plot_smoothed_proportion(
        self,
        ax: plt.Axes,
        clusters: Sequence[Any],
        y_offset: Mapping[Any, float],
        alpha: float = 0.8,
    ) -> Tuple[Mapping[Any, np.ndarray], Mapping[Any, PolyCollection]]:
        start_t, end_t = self._cmat.columns.min(), self._cmat.columns.max()
        x = np.array(self._cmat.columns)  # fitting
        # extrapolation
        e = np.linspace(start_t, end_t, int(1 + (end_t - start_t) * 100))

        smoothed_proportion, handles = {}, {}
        for clust in clusters:
            y = self._cmat.loc[clust]
            f = interp1d(x, y)
            fe = f(e)
            lo = lowess(fe, e, frac=0.3, is_sorted=True, return_sorted=False)
            smoothed_proportion[clust] = lo

            handles[clust] = ax.fill_between(
                e,
                y_offset[clust] + lo,
                y_offset[clust] - lo,
                color=self.cmap[clust],
                label=clust,
                alpha=alpha,
                edgecolor=None,
            )

        return smoothed_proportion, handles

    def _draw_flow_edge(
        self,
        ax,
        x1: Point,
        x2: Point,
        y1: Point,
        y2: Point,
        start_color: Tuple[float, float, float],
        end_color: Tuple[float, float, float],
        flow: float,
        alpha: float = 0.8,
    ) -> None:
        # transcribed from: https://github.com/tanaylab/embflow/blob/main/scripts/generate_paper_figures/plot_vein.r
        dx = x2.xt - x1.x
        dy = y2.xt - y1.x
        dxt = x2.x - x1.x
        dyt = y2.x - y1.xt

        start_color = np.asarray(to_rgb(start_color))
        end_color = np.asarray(to_rgb(end_color))
        delta = 0.05

        beta0 = _lcdf(0)
        beta_f = _lcdf(1) - _lcdf(0)

        rs = np.arange(0, 1, delta)
        beta = (_lcdf(rs) - beta0) / beta_f
        beta5 = (_lcdf(rs + delta) - beta0) / beta_f

        sx1 = x1.x + rs * dx
        sy1 = y1.x + beta * dy
        sx2 = x1.x + (rs + delta) * dx
        sy2 = y1.x + beta5 * dy

        sx1t = x1.x + flow + rs * dxt
        sy1t = y1.xt + beta * dyt
        sx2t = x1.x + flow + (rs + delta) * dxt
        sy2t = y1.xt + beta5 * dyt

        xs = np.c_[sx1, sx2, sx2t, sx1t]
        ys = np.c_[sy1, sy2, sy2t, sy1t]

        start_alpha, end_alpha = 0.2, alpha
        if start_alpha > end_alpha:
            start_alpha, end_alpha = end_alpha, start_alpha
        col = np.c_[
            (start_color * (1 - rs[:, None])) + (end_color * rs[:, None]),
            np.linspace(start_alpha, end_alpha, len(rs)),
        ]

        for x, y, c in zip(xs, ys, col):
            ax.fill(x, y, c=c, edgecolor=None)

    def _plot(
        self,
        old_times: Sequence[Numeric_t],
        ascending: Optional[bool],
        min_flow: float = 0,
        alpha: float = 0.8,
        xticks_step_size: Optional[int] = 1,
        legend_loc: Optional[str] = "upper right out",
        figsize: Optional[Tuple[float, float]] = None,
        dpi: Optional[int] = None,
    ) -> plt.Axes:
        from cellrank.pl._utils import _position_legend

        def r(num: float) -> int:
            return max(0, int(round(num, 2) * 100) - 1)

        def draw_edges(
            curr_t: Numeric_t,
            next_t: Numeric_t,
            clusters: Sequence[Any],
            *,
            bottom: bool,
        ):
            smooth_cluster = float(smoothed_proportions[self._cluster][r(curr_t)])
            flow = self._flow.loc[curr_t]
            for clust in clusters:
                fl = flow.loc[self._cluster, clust]
                if fl > min_flow:
                    fl = np.clip(fl, 0, 0.95)
                    smooth_cluster_fl = smoothed_proportions[self._cluster][
                        r(curr_t + fl)
                    ]

                    if bottom:
                        self._draw_flow_edge(
                            ax,
                            x1=Point(curr_t, 0),
                            x2=Point(next_t - fl, next_t - fl - 0.05),
                            y1=Point(
                                cluster_offset - smooth_cluster,
                                cluster_offset - smooth_cluster_fl,
                            ),
                            y2=Point(
                                y_offset[clust]
                                + smoothed_proportions[clust][r(next_t)],
                                y_offset[clust]
                                + smoothed_proportions[clust][r(next_t - fl - 0.05)],
                            ),
                            flow=fl,
                            start_color=self.cmap[self._cluster],
                            end_color=self.cmap[clust],
                            alpha=alpha,
                        )
                    else:
                        self._draw_flow_edge(
                            ax,
                            x1=Point(curr_t + fl, 0),
                            x2=Point(next_t - 0.05, next_t),
                            y1=Point(
                                cluster_offset + smooth_cluster_fl,
                                cluster_offset + smooth_cluster,
                            ),
                            y2=Point(
                                y_offset[clust]
                                - smoothed_proportions[clust][r(next_t - fl - 0.05)],
                                y_offset[clust]
                                - smoothed_proportions[clust][r(next_t)],
                            ),
                            flow=-fl,
                            start_color=self.cmap[self._cluster],
                            end_color=self.cmap[clust],
                            alpha=alpha,
                        )

        if xticks_step_size is not None:
            xticks_step_size = max(1, xticks_step_size)
        times = self._cmat.columns
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

        clusters_bottom, clusters_top = self._order_clusters(self._cluster, ascending)
        all_clusters = clusters_bottom + [self._cluster] + clusters_top

        y_offset = self._calculate_y_offsets(all_clusters)
        cluster_offset = y_offset[self._cluster]

        smoothed_proportions, handles = self._plot_smoothed_proportion(
            ax, all_clusters, y_offset, alpha=alpha
        )

        for curr_t, next_t in zip(times[:-1], times[1:]):
            draw_edges(curr_t, next_t, clusters_bottom, bottom=True)
            draw_edges(curr_t, next_t, clusters_top, bottom=False)

        ax.margins(0.025)
        ax.set_title(self._cluster)
        ax.set_xlabel(self._tkey)
        ax.set_ylabel(self._ckey)
        if xticks_step_size is None:
            ax.set_xticks([])
        else:
            ax.set_xticks(times[::xticks_step_size])
            ax.set_xticklabels(old_times[::xticks_step_size])
        ax.set_yticks([])

        if legend_loc not in (None, "none"):
            _position_legend(
                ax,
                legend_loc=legend_loc,
                handles=[handles[c] for c in all_clusters[::-1]],
            )

        return ax

    @property
    def clusters(self) -> pd.Series:
        """Clusters."""
        return self._adata.obs[self._ckey]

    @property
    def time(self) -> pd.Series:
        """Time points."""
        return self._adata.obs[self._tkey]

    @property
    @lru_cache(1)
    def cmap(self) -> Mapping[str, Any]:
        """Colormap for :attr:`clusters`."""
        return dict(
            zip(
                self.clusters.cat.categories,
                self._adata.uns.get(
                    f"{self._ckey}_colors",
                    _create_categorical_colors(len(self.clusters.cat.categories)),
                ),
            )
        )


def _lcdf(
    x: Union[int, float, np.ndarray], loc: float = 0.5, scale: float = 0.2
) -> float:
    return logistic.cdf(x, loc=loc, scale=scale)
