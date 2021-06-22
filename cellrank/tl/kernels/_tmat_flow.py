from typing import Any, List, Tuple, Union, Mapping, Optional, Sequence
from functools import lru_cache
from dataclasses import dataclass

from statsmodels.nonparametric.smoothers_lowess import lowess

from anndata import AnnData
from cellrank.tl._utils import _unique_order_preserving

import numpy as np
import pandas as pd
from scipy.stats import logistic
from scipy.sparse import issparse, spmatrix
from scipy.interpolate import interp1d

import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb


@dataclass(frozen=True)
class Point:
    x: float
    xt: float


class FlowPlotter:
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

        self._cluster = None
        self._clusters = None

        self._flow = None
        self._cmat = None

    def prepare(
        self,
        cluster: str,
        clusters: Optional[Sequence[Any]] = None,
        time_points: Optional[Sequence[Tuple[Any, Any]]] = None,
    ) -> "FlowPlotter":
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

        self._cluster = cluster
        self._cmat = self.compute_contigency_matrix()
        self._flow = self.compute_flow(time_points, cluster)

        return self

    def compute_flow(
        self, time_points: Sequence[Tuple[Any, Any]], cluster: Optional[str] = None
    ) -> pd.DataFrame:
        def default_helper(t1: Any, t2: Any) -> pd.DataFrame:
            subset, row_cls, col_cls = self._get_time_subset(t1, t2)

            df = pd.DataFrame(subset.A if issparse(subset) else subset)
            df = df.groupby(row_cls).sum().T.groupby(col_cls).sum().T

            res = pd.DataFrame(np.zeros((n, n)), index=categories, columns=categories)
            res.loc[df.index, df.columns] = df

            return res.fillna(0)

        def cluster_helper(t1: Any, t2: Any) -> pd.DataFrame:
            subset, row_cls, col_cls = self._get_time_subset(t1, t2, cluster=cluster)

            df = pd.DataFrame(subset.A if issparse(subset) else subset).sum(0)
            df = df.groupby(col_cls).sum()
            df = pd.DataFrame([df], index=[cluster], columns=df.index)

            res = pd.DataFrame(np.zeros((1, n)), index=[cluster], columns=categories)
            res.loc[df.index, df.columns] = df

            return res.fillna(0)

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

        # `[time point, 1] x clusters`
        flow /= flow.sum(1)[:, None]

        return flow.fillna(0)

    def compute_contigency_matrix(self) -> pd.DataFrame:
        cmat = pd.crosstab(self.clusters, self.time)
        return (cmat / cmat.sum(0)[None, :]).fillna(0)

    def _get_time_subset(
        self, t1: Any, t2: Any, cluster: Optional[str] = None
    ) -> Tuple[Union[np.ndarray, spmatrix], pd.Series, pd.Series]:
        if cluster is None:
            row_ixs = np.where(self.time == t1)[0]
        else:
            row_ixs = np.where((self.time == t1) & (self.clusters == cluster))[0]

        col_ixs = np.where(self.time == t2)[0]
        row_cls = self.clusters.values[row_ixs]
        col_cls = self.clusters.values[col_ixs]

        return self._tmat[row_ixs, :][:, col_ixs], row_cls, col_cls

    @property
    def adata(self) -> AnnData:
        return self._adata

    @property
    def clusters(self) -> pd.Series:
        return self._adata.obs[self._ckey]

    @property
    def time(self) -> pd.Series:
        return self._adata.obs[self._tkey]

    @property
    @lru_cache(1)
    def cmap(self) -> Mapping[str, Any]:
        return dict(
            zip(
                self.clusters.cat.categories,
                self._adata.uns[f"{self._ckey}_colors"],
            )
        )

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
    ) -> Mapping[Any, np.ndarray]:
        start_t, end_t = self._cmat.columns.min(), self._cmat.columns.max()
        x = np.array(self._cmat.columns)  # fitting
        # extrapolation
        e = np.linspace(start_t, end_t, int(1 + (end_t - start_t) * 100))

        smoothed_proportion = {}
        for clust in clusters:
            y = self._cmat.loc[clust]
            f = interp1d(x, y)
            fe = f(e)
            lo = lowess(fe, e, frac=0.3, is_sorted=True, return_sorted=False)
            smoothed_proportion[clust] = lo

            ax.fill_between(
                e,
                y_offset[clust] + lo,
                y_offset[clust] - lo,
                color=self.cmap[clust],
                label=clust,
                alpha=alpha,
                edgecolor=None,
            )

        return smoothed_proportion

    def plot(
        self,
        ascending: Optional[bool],
        min_flow: float = 0,
        alpha: float = 0.8,
        legend_loc: Optional[str] = "upper right out",
        figsize: Optional[Tuple[float, float]] = None,
        dpi: Optional[int] = None,
    ) -> Tuple[plt.Figure, plt.Axes]:
        from cellrank.pl._utils import _position_legend

        def r(num: float) -> int:
            return max(0, int(round(num, 2) * 100) - 1)

        def draw_edges(
            curr_t: Any, next_t: Any, clusters: Sequence[Any], *, bottom: bool
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

        old_times = times = self._cmat.columns
        tmp = np.array(times)
        tmp = (tmp - tmp.min()) / (tmp.max() - tmp.min())
        tmp /= np.min(tmp[1:] - tmp[:-1])
        time_mapper = dict(zip(times, tmp))
        self._flow.index = pd.MultiIndex.from_tuples(
            [(time_mapper[t], c) for t, c in self._flow.index]
        )
        self._cmat.columns = tmp
        times = self._cmat.columns

        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

        clusters_bottom, clusters_top = self._order_clusters(self._cluster, ascending)
        all_clusters = clusters_bottom + [self._cluster] + clusters_top

        y_offset = self._calculate_y_offsets(all_clusters)
        cluster_offset = y_offset[self._cluster]

        smoothed_proportions = self._plot_smoothed_proportion(
            ax, all_clusters, y_offset, alpha=alpha
        )

        for curr_t, next_t in zip(times[:-1], times[1:]):
            draw_edges(curr_t, next_t, clusters_bottom, bottom=True)
            draw_edges(curr_t, next_t, clusters_top, bottom=False)

        ax.margins(0.025)
        ax.set_title(self._cluster)
        ax.set_xlabel(self._tkey)
        ax.set_ylabel(self._ckey)
        ax.set_xticks(times)
        ax.set_xticklabels(old_times)
        ax.set_yticks([])
        if legend_loc not in (None, "none"):
            _position_legend(ax, legend_loc)

        return fig, ax


def _lcdf(
    x: Union[int, float, np.ndarray], loc: float = 0.5, scale: float = 0.2
) -> float:
    return logistic.cdf(x, loc=loc, scale=scale)
