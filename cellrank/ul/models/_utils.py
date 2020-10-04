# -*- coding: utf-8 -*-
from typing import List, Union, TypeVar, Optional

import numpy as np
from scipy.stats import rankdata
from scipy.sparse import issparse, spmatrix
from sklearn.utils.sparsefuncs import csc_median_axis_0

import cellrank.logging as logg
from cellrank.ul._docs import d
from cellrank.ul._utils import valuedispatch
from cellrank.tl._constants import ModeEnum
from cellrank.ul._parallelize import parallelize

# sources:
# edgeR: https://github.com/jianjinxu/edgeR/blob/master/R/calcNormFactors.R

AnnData = TypeVar("AnnData")


class NormMode(ModeEnum):  # noqa
    TMM = "tmm"
    RLE = "rle"
    UPPER_QUANT = "upper_quant"
    NONE = "none"


@d.dedent
def _extract_data(
    data: AnnData, layer: str, use_raw: bool = True
) -> Union[np.ndarray, spmatrix]:
    """
    Extract expression data from an object.

    Parameters
    ----------
    data
        Annotated data object or an array.
    layer
        Key in ``adata.layers`` accessed when ``use_raw=False``. Only used when ``data`` is :class:`anndata.AnnData`.
    use_raw
        Whether to access ``adata.raw``. Only used when ``data`` is :class:`anndata.AnnData`.

    Returns
    -------
    :class:`numpy.ndarray` or :class`scipy.sparse.spmatrix`
        The extracted expression data.
    """

    from anndata import AnnData as _AnnData

    if isinstance(data, _AnnData):
        if use_raw:
            if not hasattr(data, "raw"):
                raise AttributeError("No `.raw` attribute found.")
            elif data.raw is None:
                raise ValueError("Attribute `.raw` is None.")
            x = data.raw.X
        elif layer is not None:
            if layer not in data.layers:
                raise KeyError(
                    f"Layer `{layer!s}` not found in `adata.layers`. "
                    f"Valid options are: `{list(data.layers.keys())}`."
                )
            x = data.layers[layer]
        else:
            x = data.X
    elif not isinstance(data, (np.ndarray, spmatrix)):
        raise TypeError(
            f"Expected parameter `data` to be either `anndata.AnnData`, `numpy.ndarray` "
            f"or `scipy.sparse.spmatrix`, found `{type(data).__name__!r}`."
        )
    else:
        x = data

    return x


def _calculate_norm_factors(
    data: Union[AnnData, np.ndarray, spmatrix],
    method: str = NormMode.TMM.s,
    layer: Optional[str] = None,
    use_raw: bool = True,
    library_size: Optional[np.ndarray] = None,
    ref_ix: Optional[int] = None,
    logratio_trim: float = 0.3,
    sum_trim: float = 0.05,
    weight: bool = True,
    a_cutoff: float = -1e10,
    p: float = 0.75,
) -> np.ndarray:
    """
    Calculate normalization factors according to edgeR TODO: link.

    Parameters
    ----------
    data
        Data of shape `(n_cells, n_genes)` containing e.g. the counts.
    method
    layer
        Layer in ``adata.layers`` or `None` for ``adata.X``. Only used when ``use_raw=False``.
    use_raw
        Whether to access ``adata.raw``.
    library_size
    ref_ix
    logratio_trim
    sum_trim
    weight
    a_cutoff
    p

    Returns
    -------
    :class:`numpy.ndarray`
        Array of shape `(data.shape[0],)` containing the factors.
    """
    method = NormMode(method)

    x = _extract_data(data, layer, use_raw)

    if library_size is None:
        library_size = np.array(x.sum(1)).squeeze()
    elif library_size.shape != (x.shape[0],):
        raise ValueError(
            f"Expected `library_size` to be of shape `({x.shape[0]},)`, found `{library_size.shape}`."
        )

    f = _dispatch_computation(
        method,
        x=x,
        library_size=library_size,
        ref_ix=ref_ix,
        logratio_trim=logratio_trim,
        sum_trim=sum_trim,
        weight=weight,
        a_cutoff=a_cutoff,
        p=p,
    )

    return f / np.exp(np.mean(np.log(f)))


@valuedispatch
def _dispatch_computation(mode, *_args, **_kwargs):
    raise NotImplementedError(mode)


@_dispatch_computation.register(NormMode.TMM)
def _(
    x: Union[np.ndarray, spmatrix], library_size, ref_ix=None, **kwargs
) -> np.ndarray:
    if ref_ix is None:
        assert library_size is not None
        f75 = _calc_factor_quant(x, library_size=library_size, p=0.75)
        ref_ix = np.argmin(np.abs(f75 - np.mean(f75)))

    return parallelize(
        _calc_factor_weighted,
        collection=np.arange(x.shape[0]),
        show_progress_bar=False,
        as_array=False,
        extractor=lambda res: np.array([r for rs in res for r in rs]),
        backend="threading",  # TODO: numba...
        n_jobs=4,
    )(
        obs_=x,
        obs_lib_size_=library_size,
        ref=x[ref_ix],
        ref_lib_size=library_size[ref_ix],
        **kwargs,
    )


@_dispatch_computation.register(NormMode.UPPER_QUANT)
def _calc_factor_quant(
    x: Union[np.ndarray, spmatrix], library_size: np.ndarray, p: float, **_
) -> np.ndarray:
    library_size = np.array(library_size).reshape((-1, 1))
    if not issparse(x):
        # this is same as below, which is R's default
        # return mquantiles(x / library_size, prob=q, alphap=1, betap=1, axis=1).data.squeeze()
        return np.quantile(x / library_size, p, axis=1)

    # not very efficient
    y = x.multiply(1.0 / library_size).tocsr()
    return np.array([np.quantile(y[i].A[0], p) for i in range(y.shape[0])])


@_dispatch_computation.register(NormMode.RLE)
def _(x: Union[np.ndarray, spmatrix], **_) -> np.ndarray:
    gm = np.array(np.exp(np.sum(np.log1p(x), axis=0))).squeeze()
    mask = (gm > 0) & np.isfinite(gm)
    if not issparse(x):
        return np.median(x[:, mask] / gm[mask], axis=1)

    return csc_median_axis_0(x.tocsr()[:, mask].multiply(1.0 / gm[mask]).tocsr().T)


@_dispatch_computation.register(NormMode.NONE)
def _(x: Union[np.ndarray, spmatrix], **_) -> np.ndarray:
    return np.ones((x.shape[0],), dtype=x.dtype)


def _calc_factor_weighted(
    ixs: np.ndarray,
    obs_: Union[np.ndarray, spmatrix],
    ref: Union[np.ndarray, spmatrix],
    obs_lib_size_: np.ndarray,
    ref_lib_size: Optional[float] = None,
    logratio_trim: float = 0.3,
    sum_trim: float = 0.05,
    weight: bool = True,
    a_cutoff: float = -1e10,
    queue=None,
    **_,
) -> List[float]:
    # TODO: vectorize, not parallelize
    if issparse(ref):
        ref = ref.A.squeeze(0)  # 1 x genes
    if ref_lib_size is None:
        ref_lib_size = np.sum(ref)

    r_scaled = ref / ref_lib_size
    log_ref = np.log2(r_scaled)

    res = []

    for ix in ixs:
        if queue is not None:
            queue.put(1)

        obs = obs_[ix]

        if issparse(obs):
            obs = obs.A.squeeze(0)  # 1 x genes

        obs_lib_size = obs_lib_size_[ix]

        if obs_lib_size is None:
            obs_lib_size = np.sum(obs)

        o_scaled = obs / obs_lib_size
        log_obs = np.log2(o_scaled)

        log_r = log_obs - log_ref
        abs_e = (log_obs + log_ref) / 2.0
        v = (obs_lib_size - obs) / o_scaled + (ref_lib_size - ref) / r_scaled

        mask = np.isfinite(log_r) & np.isfinite(abs_e) & (abs_e > a_cutoff)

        log_r = log_r[mask]
        abs_e = abs_e[mask]
        v = v[mask]

        if not len(log_r) or np.max(np.abs(log_r)) < 1e-6:
            res.append(1.0)
            continue

        n = len(log_r)
        lol = np.floor(n * logratio_trim) + 1.0
        hil = n + 1 - lol

        los = np.floor(n * sum_trim) + 1.0
        his = n + 1 - los

        rank_log_r, rank_abs_e = rankdata(log_r), rankdata(abs_e)
        mask = (
            (rank_log_r >= lol)
            & (rank_log_r <= hil)
            & (rank_abs_e >= los)
            & (rank_abs_e <= his)
        )

        f = (
            (np.nansum(log_r[mask] / v[mask]) / np.nansum(1.0 / v[mask]))
            if weight
            else np.nanmean(log_r[mask])
        )

        res.append(1.0 if not np.isfinite(f) else 2 ** f)

    if queue is not None:
        queue.put(None)

    return res


def _get_knotlocs(
    pseudotime: np.ndarray,
    n_knots: int,
) -> np.ndarray:
    """
    Find knot locations.

    Parameters
    ----------
    pseudotime
        Pseudotemporal ordering of cells.
    n_knots
        Number of knots.

    Returns
    -------
    :class:`numpy.ndarray`
        Array of shape `(n_knots,)` containing the locations of knots along the pseudotime.
    """

    if n_knots <= 0:
        raise ValueError(f"Expected number of knots to be positive, found `{n_knots}`.")

    if np.any(~np.isfinite(pseudotime)):
        raise ValueError("Not all pseudotime values are finite.")

    pseudotime = np.asarray(pseudotime)
    if pseudotime.ndim != 1:
        raise ValueError(
            f"Expected pseudotime to have `1` dimension, found `{pseudotime.ndim}`."
        )

    x = np.quantile(
        pseudotime, q=np.arange(n_knots, dtype=np.float64) / max(n_knots - 1, 1)
    )
    u, ix, c = np.unique(x, return_index=True, return_counts=True)

    if len(u) != len(x):
        knotlocs = []
        for start, end, size in zip(x[ix], x[ix[1:]], c):
            knotlocs.extend(np.linspace(start, end, size, endpoint=False))
        knotlocs.extend(
            np.linspace(knotlocs[-1], x[ix[-1]], c[-1] + 1, endpoint=True)[1:]
        )
        knotlocs = np.array(knotlocs)
    else:
        knotlocs = x

    knotlocs[0] = np.min(pseudotime)
    knotlocs[-1] = np.max(pseudotime)

    logg.debug(f"Setting knot locations to `{list(knotlocs)}`")

    return knotlocs


@d.dedent
def _get_offset(
    adata: AnnData, use_raw: bool = True, layer: Optional[str] = None, **kwargs
) -> np.ndarray:
    """
    Return an offset for GAM.

    Parameters
    ----------
    %(adata)s
    use_raw
        Whether to access ``adata.raw``.
    layer
        Layer in ``adata.layers`` or `None` for ``adata.X``. Only used when ``use_raw=False``.
    **kwargs
        Keyword arguments for :func:`cellrank.ul.models._utils._calc_norm_factors`.

    Returns
    -------
    :class:`numpy.ndarray`
        Array of shape `(adata.n_obs,)` containing the offset.
    """

    data = _extract_data(adata, layer=layer, use_raw=use_raw)
    try:
        nf = _calculate_norm_factors(adata, layer=layer, use_raw=use_raw, **kwargs)
    except Exception as e:
        logg.debug(f"Unable to calculate the normalization factors. Reason: `{e}`")
        nf = np.ones(adata.n_obs, dtype=np.float64)

    # TODO; handle zeros
    return np.log(nf * np.array(data.sum(1)).squeeze())
