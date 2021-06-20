from typing import Any, Tuple, Union, Sequence

from statsmodels.nonparametric.smoothers_lowess import lowess

from anndata import AnnData

import numpy as np
import pandas as pd
from scipy.stats import logistic
from scipy.sparse import issparse, spmatrix
from scipy.interpolate import interp1d

import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb


def _compute_flow(
    T: Union[np.ndarray, spmatrix],
    clusters: pd.Series,
    time: pd.Series,
    time_points: Sequence[Tuple[Any, Any]],
) -> pd.DataFrame:
    def flow(t1: Any, t2: Any) -> pd.DataFrame:
        row_ixs = np.where(time == t1)[0]
        col_ixs = np.where(time == t2)[0]
        row_cls = clusters.values[row_ixs]
        col_cls = clusters.values[col_ixs]

        subset = T[row_ixs, :][:, col_ixs]
        if issparse(subset):
            df = pd.DataFrame(subset.A)
        else:
            df = pd.Datarame(subset)
        df = df.groupby(row_cls).sum().T.groupby(col_cls).sum()
        # TODO: normalize? if not, how to handle offset when plotting?
        vals = df.values / df.values.sum(1)[:, None]

        res = pd.DataFrame(np.zeros((n, n)), index=categories, columns=categories)
        res.loc[df.index, df.columns] = vals
        return res.fillna(0)

    categories = clusters.cat.categories
    n = len(categories)
    tfs = []

    for t1, t2 in time_points:
        type_flow = flow(t1, t2)
        # TODO: (multi)index? (or at lea
        type_flow["t1"] = str(t1)
        tfs.append(type_flow)

    return pd.concat(tfs)


# TODO: rename
def _norm_cont(adata: AnnData, cluster_key: str, time_key: str) -> pd.DataFrame:
    df = pd.crosstab(adata.obs[cluster_key], adata.obs[time_key])
    return (df / df.sum(0)[None, :]).fillna(0)


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
    col = np.c_[
        (start_color * (1 - rs[:, None])) + (end_color * rs[:, None]),
        np.linspace(0.15, 0.75, len(rs)),
    ]

    for x, y, c in zip(xs, ys, col):
        ax.fill(x, y, c=c, edgecolor=None)


# TODO:
# 1. flow normalization?
# 2. cluster filtering tweak + flow filtering
# 3. z-order (in case of potential overlap
def _plot_flow(
    adata: AnnData,
    cluster: str,
    cluster_key: str,
    time_key: str,
    type_agn: pd.DataFrame,
    type_flow: pd.DataFrame,
) -> None:
    def plot_edges_bottom(j: int, t: Any, ixs: np.array):
        cum_y = float(smoo_y2[cluster][f"{float(t):.2f}"])
        flow = type_flow[type_flow["t1"].astype(float) == float(t)]
        for i in ixs:
            col_i = cols[i]
            fl = flow.loc[cluster, col_i]
            # TODO: threshold (+careful normalization)
            if fl > 0:
                fl = np.clip(fl, 0, 0.95)
                ix2 = f"{float((x[j + 1])):.2f}"
                ix3 = f"{float(x[j + 1] - fl - 0.05):.2f}"

                _draw_sig_edge(
                    ax,
                    x1=t,
                    x2=x[j + 1] - fl,
                    x2t=x[j + 1] - fl - 0.05,
                    fl=fl,
                    y1=base_foc - cum_y,
                    y1t=base_foc - cum_y,
                    y2=base_y[col_i] + smoo_y2[col_i][ix2],
                    y2t=base_y[col_i] + smoo_y2[col_i][ix3],
                    start_color=cm[cluster],
                    end_color=cm[col_i],
                )

    def plot_edges_top(j: int, t: Any, ixs: np.array):
        cum_y = float(smoo_y2[cluster][f"{float(t):.2f}"])
        flow = type_flow[type_flow["t1"].astype(float) == float(t)]
        for i in ixs:
            col_i = cols[i]
            fl = flow.loc[cluster, col_i]
            # TODO: threshold (+careful normalization)
            if fl > 0:
                fl = np.clip(fl, 0, 0.95)
                ix2 = f"{float((x[j + 1])):.2f}"
                ix3 = f"{float(x[j + 1] - fl - 0.05):.2f}"

                _draw_sig_edge(
                    ax,
                    x1=t + fl,
                    x2t=x[j + 1],
                    x2=x[j + 1] - fl - 0.05,
                    y1=base_foc + cum_y,
                    y1t=base_foc + cum_y,
                    y2t=base_y[col_i] - smoo_y2[col_i][ix2],
                    y2=base_y[col_i] - smoo_y2[col_i][ix3],
                    fl=-fl,
                    start_color=cm[cluster],
                    end_color=cm[col_i],
                )

    # TODO: clean file, remove constant
    t1 = type_agn.columns[0]
    t2 = type_agn.columns[-1]
    T_minflow_for_type = 0.005
    agg = type_flow.loc[cluster].select_dtypes(exclude=["object"]).sum()
    cols = agg.index[agg > T_minflow_for_type]
    cols = [c for c in cols if c != cluster]
    # TODO: sort left/right by total descending flow
    cols = cols[: len(cols) // 2] + [cluster] + cols[len(cols) // 2 :]
    cols = np.array(cols)
    foc_agn = type_agn.select_dtypes(exclude=["object"]).loc[list(cols), t1:t2]

    base_y = [0]
    for i in range(1, len(cols)):
        # TODO: tweak
        base_y.append(
            base_y[-1]
            + 0.2  # TODO
            + np.max(foc_agn.loc[str(cols[i])] + foc_agn.loc[str(cols[i - 1])])
        )
    base_y = dict(zip(cols, base_y))

    x = np.array(type_agn.columns)
    length = int(1 + (t2 - t1) * 100)
    e = np.linspace(t1, t2, length)
    # TODO: expose DPI
    fig, ax = plt.subplots(dpi=180)

    smoo_y2 = {}
    # TODO: what if no colors present?
    cm = dict(
        zip(adata.obs[cluster_key].cat.categories, adata.uns[f"{cluster_key}_colors"])
    )

    base = [0]

    for i, c in enumerate(cols):
        y = foc_agn.loc[c, t1:t2]
        f = interp1d(x, y)
        fe = f(e)
        lo = lowess(fe, e, frac=0.3, is_sorted=True, return_sorted=False)
        # TODO: find cleaner way
        smoo_y2[c] = {f"{float(k):.2f}": v for k, v in zip(e, lo)}

        ax.fill_between(e, lo + base_y[c], -lo + base_y[c], color=cm[c], label=c)
        base.append(base[i] + 1)

    foc_i = np.where(cols == cluster)[0][0]
    base_foc = base_y[cluster]

    for j, t in enumerate(x[:-1]):
        plot_edges_bottom(j, t, np.arange(foc_i)[::-1])
        plot_edges_top(j, t, np.arange(foc_i + 1, len(cols)))

    ax.set_title(cluster)
    ax.set_xlabel(time_key)
    ax.set_xticks(x)
    ax.set_yticks([])
    ax.set_ylabel(cluster_key)
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))  # TODO: expose pos
