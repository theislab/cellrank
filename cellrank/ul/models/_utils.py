# -*- coding: utf-8 -*-
from typing import Union, TypeVar, Optional, Sequence

import numpy as np
from numba import njit, prange
from scipy.sparse import issparse, spmatrix
from sklearn.utils.sparsefuncs import csc_median_axis_0

import cellrank.logging as logg
from cellrank.ul._docs import d, inject_docs
from cellrank.ul._utils import valuedispatch
from cellrank.tl._constants import ModeEnum
from cellrank.ul._parallelize import parallelize

AnnData = TypeVar("AnnData")
_OFFSET_KEY = "cellrank_offset"


class NormMode(ModeEnum):  # noqa
    TMM = "tmm"
    RLE = "rle"
    UPPER_QUANT = "upper_quant"
    NONE = "none"


@d.dedent
def _extract_data(
    data: AnnData, layer: Optional[str] = None, use_raw: bool = True
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


@njit(fastmath=True, cache=True)
def _rankdata(a: np.ndarray, method: str = "average") -> np.ndarray:
    assert method in (
        "average",
        "min",
        "max",
        "dense",
        "ordinal",
    ), "Invalid ranking method."

    arr = np.ravel(np.asarray(a))
    # algo = 'mergesort' if method == 'ordinal' else 'quicksort'
    sorter = np.argsort(arr)  # , kind=algo)

    inv = np.empty(sorter.size, dtype=np.intp)
    inv[sorter] = np.arange(sorter.size, dtype=np.intp)

    if method == "ordinal":
        return inv + 1.0

    arr = arr[sorter]
    obs = np.ones(len(arr), dtype=np.bool_)
    obs[1:] = arr[1:] != arr[:-1]
    dense = obs.cumsum()[inv]

    if method == "dense":
        return dense * 1.0

    # cumulative counts of each unique value
    count = np.append(np.nonzero(obs)[0], len(obs))

    if method == "max":
        return count[dense] * 1.0

    if method == "min":
        return count[dense - 1] + 1.0

    # average method
    return 0.5 * (count[dense] + count[dense - 1] + 1.0)


@inject_docs(m=NormMode)
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
    perc: float = 0.75,
) -> np.ndarray:
    """
    Calculate normalization factors according to
    `edgeR <https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/calcNormFactors>`__.

    Parameters
    ----------
    data
        Data of shape `(n_cells, n_genes)` containing e.g. the counts.
    method
        One of following:

            - method={m.TMM.s!r} - weighted trimmed mean of M-values from [Robinson10]_.
            - method={m.RLE.s!r} - relative log expression from [Anders10]_.
            - method={m.UPPER_QUANT.s!r} - upper-quartile normalization method from [Bullard10]_.
            - method={m.NONE.s!r} - all the factors are set to 1.
    layer
        Layer in ``adata.layers`` or `None` for ``adata.X``. Only used when ``use_raw=False``.
    use_raw
        Whether to access ``adata.raw``.
    library_size
        Library size. If `None`, it will be set as the sum of values for each cell.
    ref_ix
        Index of a reference cell. If `None`, it will be determined automatically.
    logratio_trim
        Amount of trim to use on log-ratios ('M' values) when ``method={m.TMM.s!r}``.
    sum_trim
        Amount of trim to use on the combined absolute levels ('A' values) when ``method={m.TMM.s!r}``.
    weight
        Whether to compute asymptotic binomial precision weights when ``method={m.TMM.s!r}``.
    a_cutoff
        Cutoff on 'A' values to use before trimming when ``method={m.TMM.s!r}``.
    perc
        Percentile of the counts that is aligned when ``method={m.UPPER_QUANT.s!r}``.

    Returns
    -------
    :class:`numpy.ndarray`
        Array of shape `(data.shape[0],)` containing the factors.

    References
    ----------
    .. [Anders10] Anders, S. *et al.* (2010),
        *Differential expression analysis for sequence count data*,
        `Genome Biology <https://doi.org/10.1186/gb-2010-11-10-r106>`__.
    .. [Bullard10] Bullard, J. H. *et al.* (2010),
        *Evaluation of statistical methods for normalization and differential expression in mRNA-Seq experiments*,
        `BMC Bioinformatics <https://doi.org/10.1186/1471-2105-11-94>`__.
    """  # noqa

    method = NormMode(method)

    if not (0 <= perc <= 1):
        raise ValueError(
            f"Expected the percentile to be within interval `[0, 1]`, found `{perc}`."
        )

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
        p=perc,
    )

    return f / np.exp(np.mean(np.log(f)))


@valuedispatch
def _dispatch_computation(mode, *_args, **_kwargs):
    raise NotImplementedError(mode)


@_dispatch_computation.register(NormMode.TMM)
def _(
    x: Union[np.ndarray, spmatrix],
    library_size: np.ndarray,
    ref_ix: Optional[int] = None,
    **kwargs,
) -> np.ndarray:
    if ref_ix is None:
        f75 = _calc_factor_quant(x, library_size=library_size, p=0.75)
        ref_ix = np.argmin(np.abs(f75 - np.mean(f75)))

    ref = x[ref_ix]
    if issparse(ref):
        ref = ref.A.squeeze(0)  # (genes,)

    return parallelize(
        _calc_factor_weighted,
        collection=np.arange(x.shape[0]),
        show_progress_bar=False,
        as_array=False,
        extractor=np.concatenate,
        backend="threading",
        n_jobs=4,
    )(
        obs_=x,
        obs_lib_size_=library_size,
        ref=ref,
        ref_lib_size=library_size[ref_ix],
        **kwargs,
    )


@_dispatch_computation.register(NormMode.UPPER_QUANT)
def _calc_factor_quant(
    x: Union[np.ndarray, spmatrix], library_size: np.ndarray, p: float, **_
) -> np.ndarray:
    # unused, because we fix the reference cell to 0

    library_size = np.array(library_size).reshape((-1, 1))
    if not issparse(x):
        # this is same as below, which is R's default
        # return mquantiles(x / library_size, prob=q, alphap=1, betap=1, axis=1).data.squeeze()
        return np.quantile(x / library_size, p, axis=1)

    y = x.multiply(1.0 / library_size).tocsr()
    return _calc_factor_quant_sparse_helper(y.data, y.indptr, p=p, n_cols=y.shape[1])


@njit(parallel=True, fastmath=True, cache=True)
def _calc_factor_quant_sparse_helper(
    data: np.ndarray,
    indptr: np.ndarray,
    p: float,
    n_cols: int,
) -> np.ndarray:
    n = len(indptr) - 1
    out = np.empty(n)

    # maybe disable parallel
    for ix in prange(n):
        d = data[indptr[ix] : indptr[ix + 1]]
        tmp = np.zeros(shape=(n_cols,))
        tmp[: len(d)] = d
        out[ix] = np.quantile(tmp, p)

    return out


@_dispatch_computation.register(NormMode.RLE)
def _(x: Union[np.ndarray, spmatrix], **_) -> np.ndarray:
    # unused

    # log 0 -> inf -> masked out anyway, so we select only genes in every cell
    # this is the exact replica an in edgeR
    mask = np.array((x > 0).sum(0)).squeeze() == x.shape[0]
    tmp = x[:, mask]
    if issparse(tmp):
        tmp = tmp.A

    gm = np.array(np.exp(np.mean(np.log(tmp), axis=0))).squeeze()
    gm_mask = (gm > 0) & np.isfinite(gm)

    if not issparse(x):
        return np.median(x[:, mask][:, gm_mask] / gm[gm_mask], axis=1)

    return csc_median_axis_0(
        x.tocsr()[:, mask][:, gm_mask].multiply(1.0 / gm[gm_mask]).tocsr().T
    )


@_dispatch_computation.register(NormMode.NONE)
def _(x: Union[np.ndarray, spmatrix], **_) -> np.ndarray:
    # unused
    return np.ones((x.shape[0],), dtype=x.dtype)


def _calc_factor_weighted(
    ixs: np.ndarray,
    obs_: Union[np.ndarray, spmatrix],
    obs_lib_size_: np.ndarray,
    ref: np.ndarray,
    ref_lib_size: float,
    logratio_trim: float = 0.3,
    sum_trim: float = 0.05,
    weight: bool = True,
    a_cutoff: float = -1e10,
    queue=None,
    **_,
) -> np.ndarray:
    log_ref = np.log2(ref / ref_lib_size)
    res = np.empty(shape=(len(ixs),))

    for i, ix in enumerate(ixs):
        if queue is not None:
            queue.put(1)

        obs = obs_[ix]
        if issparse(obs):
            obs = obs.A.squeeze(0)  # (genes,)
        obs_lib_size = obs_lib_size_[ix]

        res[i] = _calc_factor_weighted_helper(
            obs=obs,
            obs_lib_size=obs_lib_size,
            ref=ref,
            ref_lib_size=ref_lib_size,
            log_ref=log_ref,
            a_cutoff=a_cutoff,
            logratio_trim=logratio_trim,
            sum_trim=sum_trim,
            weight=weight,
        )

    if queue is not None:
        queue.put(None)

    return res


@njit(fastmath=True, cache=True)
def _calc_factor_weighted_helper(
    obs: np.ndarray,
    obs_lib_size: float,
    ref: np.ndarray,
    ref_lib_size: float,
    log_ref: np.ndarray,
    a_cutoff: float,
    logratio_trim: float,
    sum_trim: float,
    weight: bool,
) -> float:
    o_scaled = obs / obs_lib_size
    log_obs = np.log2(o_scaled)

    log_r = log_obs - log_ref
    abs_e = (log_obs + log_ref) / 2.0
    v = ((obs_lib_size - obs) / obs_lib_size / obs) + (
        ref_lib_size - ref
    ) / ref_lib_size / ref

    mask = np.isfinite(log_r) & np.isfinite(abs_e) & (abs_e > a_cutoff)

    log_r = log_r[mask]
    abs_e = abs_e[mask]
    v = v[mask]

    if not len(log_r) or np.max(np.abs(log_r)) < 1e-6:
        return 1.0

    n = len(log_r)
    lol = np.floor(n * logratio_trim) + 1.0
    hil = n + 1 - lol

    los = np.floor(n * sum_trim) + 1.0
    his = n + 1 - los

    rank_log_r, rank_abs_e = _rankdata(log_r), _rankdata(abs_e)
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

    return 2 ** f if np.isfinite(f) else 1.0


def _get_knotlocs(
    pseudotime: Union[np.ndarray, Sequence],
    n_knots: int,
    uniform: bool = False,
) -> np.ndarray:
    """
    Find knot locations.

    The first and last knots are always placed at the beginning and the of the ``pseudotime``.

    Parameters
    ----------
    pseudotime
        Pseudotemporal ordering of cells.
    n_knots
        Number of knots.
    uniform
        Whether to place the knots uniformly across the ``pseudotime`` or place them based on ``pseudotime``'s density.

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

    unique = np.unique(pseudotime)
    if len(unique) in (0, 1):
        raise ValueError(f"All pseudotime values are the same: `{unique}`.")

    if pseudotime.ndim == 2 and pseudotime.shape[1] == 1:
        pseudotime = pseudotime.squeeze(1)
    if pseudotime.ndim != 1:
        raise ValueError(
            f"Expected `pseudotime` to have `1` dimension, found `{pseudotime.ndim}`."
        )

    if uniform:
        # replicate the result from not uniform
        return (
            np.linspace(np.min(pseudotime), np.max(pseudotime), n_knots, endpoint=True)
            if n_knots > 1
            else np.array([np.max(pseudotime)])
        )

    x = np.quantile(
        pseudotime, q=np.arange(n_knots, dtype=np.float64) / max(n_knots - 1, 1)
    )
    x[0] = np.min(pseudotime)
    x[-1] = np.max(pseudotime)

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

    logg.debug(f"Setting knot locations to `{list(knotlocs)}`")

    return knotlocs


@inject_docs(key=_OFFSET_KEY)
@d.dedent
def _get_offset(
    adata: AnnData,
    layer: Optional[str] = None,
    use_raw: bool = True,
    ref_ix: Optional[int] = None,
    recompute: bool = False,
    **kwargs,
) -> np.ndarray:
    """
    Return an offset for GAM.

    Parameters
    ----------
    %(adata)s
    layer
        Layer in ``adata.layers`` or `None` for ``adata.X``. Only used when ``use_raw=False``.
    use_raw
        Whether to access ``adata.raw``.
    ref_ix
        Index of a reference cell. If `None`, it will be determined automatically.
    recompute
        Whether to recompute the offset if already found in ``adata.obs[{key!r}]``.
    **kwargs
        Keyword arguments for :func:`cellrank.ul.models._utils._calc_norm_factors`.

    Returns
    -------
    :class:`numpy.ndarray`
        Array of shape `(adata.n_obs,)` containing the offset.
    """

    from anndata import AnnData as _AnnData

    if not recompute and isinstance(adata, _AnnData) and _OFFSET_KEY in adata.obs:
        logg.debug(f"Fetching offset from `adata.obs[{_OFFSET_KEY!r}]`")
        return adata.obs[_OFFSET_KEY].values.copy()

    logg.debug(f"Calculating offset for `{adata.shape[0]}` cells")

    data = _extract_data(adata, layer=layer, use_raw=use_raw)
    try:
        nf = _calculate_norm_factors(
            adata, layer=layer, use_raw=use_raw, ref_ix=ref_ix, **kwargs
        )
    except Exception as e:
        logg.debug(
            f"Unable to calculate the normalization factors, setting them to `1`. Reason: `{e}`"
        )
        nf = np.ones(adata.n_obs, dtype=np.float64)

    offset = np.log(nf * np.array(data.sum(1)).squeeze())
    offset[offset == 0] = 1.0

    mask = ~np.isfinite(offset)
    if np.any(mask):
        logg.warning(f"`{np.sum(mask)}` elements are not finite. Setting them to `1`")
        offset[mask] = 1.0

    if isinstance(adata, _AnnData):
        adata.obs[_OFFSET_KEY] = offset

    return offset
