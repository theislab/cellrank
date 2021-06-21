from typing import Any, Tuple, Union, Optional, Sequence

from statsmodels.nonparametric.smoothers_lowess import lowess

from anndata import AnnData

import numpy as np
import pandas as pd
from scipy.stats import logistic
from scipy.sparse import issparse, spmatrix
from scipy.interpolate import interp1d

import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb


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

    def prepare(self) -> "FlowPlotter":
        return self

    def flow(
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
            times.extend([str(t1)] * len(flow))
            flows.append(flow)

        flow = pd.concat(flows)
        flow.set_index([times, flow.index], inplace=True)

        # only normalization that makes sense
        flow /= flow.sum(1)[:, None]

        return flow

    def contingency_matrix(self) -> pd.DataFrame:
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
    def clusters(self) -> pd.Series:
        return self._adata.obs[self._ckey]

    @property
    def time(self) -> pd.Series:
        return self._adata.obs[self._tkey]


def _lcdf(
    x: Union[int, float, np.ndarray], loc: float = 0.5, scale: float = 0.2
) -> float:
    return logistic.cdf(x, loc=loc, scale=scale)


# TODO: bundle args, what are those?`
def _draw_sig_edge(
    ax,
    x1: float,
    x2: float,
    x2t: float,
    y1: float,
    y1t: float,
    y2: float,
    y2t: float,
    start_color: Tuple[float, float, float],
    end_color: Tuple[float, float, float],
    fl: float,
    alpha: float = 0.8,
) -> None:
    dx = x2t - x1
    dy = y2t - y1
    dxt = x2 - x1
    dyt = y2 - y1t

    start_color = np.asarray(to_rgb(start_color))
    end_color = np.asarray(to_rgb(end_color))
    delta = 0.05

    beta0 = _lcdf(0)
    beta_f = _lcdf(1) - _lcdf(0)

    rs = np.arange(0, 1, delta)
    beta = (_lcdf(rs) - beta0) / beta_f
    beta5 = (_lcdf(rs + delta) - beta0) / beta_f

    sx1 = x1 + rs * dx
    sy1 = y1 + beta * dy
    sx2 = x1 + (rs + delta) * dx
    sy2 = y1 + beta5 * dy

    sx1t = x1 + fl + rs * dxt
    sy1t = y1t + beta * dyt
    sx2t = x1 + fl + (rs + delta) * dxt
    sy2t = y1t + beta5 * dyt

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


# TODO:
# 1. flow normalization?
# 2. cluster filtering tweak + flow filtering
# 3. z-order (in case of potential overlap
def _plot_flow(
    cm,
    cluster: str,
    cluster_key: str,
    clusters: Optional[Sequence[str]],
    ascending: Optional[bool],
    time_key: str,
    type_agn: pd.DataFrame,
    type_flow: pd.DataFrame,
    min_flow: float = 0,
    alpha: float = 0.8,
    legend_loc: Optional[str] = "upper right out",
    figsize: Optional[Tuple[float, float]] = None,
    dpi: Optional[int] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    from cellrank.pl._utils import _position_legend

    # TODO: 1 function
    def plot_edges_bottom(j: int, t: Any, ixs: np.array):
        cum_y = float(smoo_y2[cluster][f"{float(t):.2f}"])
        flow = type_flow.loc[str(t)]
        for i in ixs:
            col_i = cols[i]
            fl = flow.loc[cluster, col_i]
            if fl > min_flow:
                fl = np.clip(fl, 0, 0.95)
                cum_yt = float(smoo_y2[cluster][f"{float(t)+fl:.2f}"])
                ix2 = f"{float((x[j + 1])):.2f}"
                ix3 = f"{float(x[j + 1] - fl - 0.05):.2f}"

                _draw_sig_edge(
                    ax,
                    x1=t,
                    x2=x[j + 1] - fl,
                    x2t=x[j + 1] - fl - 0.05,
                    fl=fl,
                    y1=base_foc - cum_y,
                    y1t=base_foc - cum_yt,
                    y2=base_y[col_i] + smoo_y2[col_i][ix2],
                    y2t=base_y[col_i] + smoo_y2[col_i][ix3],
                    start_color=cm[cluster],
                    end_color=cm[col_i],
                    alpha=alpha,
                )

    def plot_edges_top(j: int, t: Any, ixs: np.array):
        cum_y = float(smoo_y2[cluster][f"{float(t):.2f}"])
        flow = type_flow.loc[str(t)]
        for i in ixs:
            col_i = cols[i]
            fl = flow.loc[cluster, col_i]
            if fl > min_flow:
                fl = np.clip(fl, 0, 0.95)
                cum_yt = float(smoo_y2[cluster][f"{float(t)+fl:.2f}"])
                ix2 = f"{float((x[j + 1])):.2f}"
                ix3 = f"{float(x[j + 1] - fl - 0.05):.2f}"

                _draw_sig_edge(
                    ax,
                    x1=t + fl,
                    x2t=x[j + 1],
                    x2=x[j + 1] - 0.05,
                    y1=base_foc + cum_yt,
                    y1t=base_foc + cum_y,
                    y2t=base_y[col_i] - smoo_y2[col_i][ix2],
                    y2=base_y[col_i] - smoo_y2[col_i][ix3],
                    fl=-fl,
                    start_color=cm[cluster],
                    end_color=cm[col_i],
                    alpha=alpha,
                )

    # TODO: clean file, remove constants
    t1 = type_agn.columns[0]
    t2 = type_agn.columns[-1]

    # TODO: extract to function
    if ascending is not None:
        top_bottom = [[], []]
        agg = (
            type_flow.loc[(slice(None), cluster), :]
            .sum()
            .sort_values(ascending=ascending)
        )
        for i, c in enumerate(c for c in agg.index if c != cluster):
            top_bottom[i % 2].append(c)
        # TODO: document
        cols = top_bottom[0][::-1] + [cluster] + top_bottom[1]
    else:
        cols = [c for c in clusters if c != cluster]
        # TODO: sort left/right by total descending flow
        cols = cols[: len(cols) // 2] + [cluster] + cols[len(cols) // 2 :]
    cols = np.array(cols)

    foc_agn = type_agn.loc[list(cols), t1:t2]

    base_y = [0]
    for i in range(1, len(cols)):
        # TODO: cleaner impl.
        base_y.append(
            base_y[-1]
            + 0.2
            + np.max(foc_agn.loc[str(cols[i])] + foc_agn.loc[str(cols[i - 1])])
        )
    base_y = dict(zip(cols, base_y))

    x = np.array(type_agn.columns)
    length = int(1 + (t2 - t1) * 100)
    e = np.linspace(t1, t2, length)
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    smoo_y2 = {}

    base = [0]
    for i, c in enumerate(cols):
        y = foc_agn.loc[c, t1:t2]
        f = interp1d(x, y)
        fe = f(e)
        lo = lowess(fe, e, frac=0.3, is_sorted=True, return_sorted=False)
        # TODO: find cleaner way
        smoo_y2[c] = {f"{float(k):.2f}": v for k, v in zip(e, lo)}

        ax.fill_between(
            e,
            lo + base_y[c],
            -lo + base_y[c],
            color=cm[c],
            label=c,
            alpha=alpha,
            edgecolor=None,
        )
        base.append(base[i] + 1)

    foc_i = np.where(cols == cluster)[0][0]
    base_foc = base_y[cluster]

    for j, t in enumerate(x[:-1]):
        plot_edges_bottom(j, t, np.arange(foc_i)[::-1])
        plot_edges_top(j, t, np.arange(foc_i + 1, len(cols)))

    ax.margins(0.025)
    ax.set_title(cluster)
    ax.set_xlabel(time_key)
    ax.set_xticks(x)
    ax.set_yticks([])
    ax.set_ylabel(cluster_key)
    if legend_loc not in (None, "none"):
        _position_legend(ax, legend_loc)

    return fig, ax
