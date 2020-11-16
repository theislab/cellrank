# -*- coding: utf-8 -*-
"""Datasets module."""

import os
from typing import Tuple, Union, TypeVar
from pathlib import Path

from scanpy import read

from cellrank import logging as logg
from cellrank.ul._docs import d, inject_docs
from cellrank.tl._constants import ModeEnum

AnnData = TypeVar("AnnData")


class ReprogrammingSubset(ModeEnum):  # noqa: D101
    FULL = "full"
    K48 = "48k"
    K85 = "85k"


_datasets = {
    "pancreas": (
        "https://ndownloader.figshare.com/files/25060877",
        (2531, 27998),
    ),
    "pancreas_preprocessed": (
        "https://ndownloader.figshare.com/files/25030028",
        (2531, 2000),
    ),
    "lung": (
        "https://ndownloader.figshare.com/files/25038224",
        (24882, 24051),
    ),
    "reprogramming": (
        "https://ndownloader.figshare.com/files/25503773",
        (104679, 22630),
    ),
}


def _load_dataset_from_url(
    fpath: Union[os.PathLike, str], url: str, expected_shape: Tuple[int, int], **kwargs
) -> AnnData:

    fpath = str(fpath)
    if not fpath.endswith(".h5ad"):
        fpath += ".h5ad"

    if os.path.isfile(fpath):
        logg.debug(f"Loading dataset from `{fpath!r}`")
    else:
        logg.debug(f"Downloading dataset from `{url!r}` as `{fpath!r}`")

    dirname, _ = os.path.split(fpath)
    try:
        if not os.path.isdir(dirname):
            logg.debug(f"Creating directory `{dirname!r}`")
            os.makedirs(dirname, exist_ok=True)
    except OSError as e:
        logg.debug(f"Unable to create directory `{dirname!r}`. Reason `{e}`")

    kwargs.setdefault("sparse", True)
    kwargs.setdefault("cache", True)

    adata = read(fpath, backup_url=url, **kwargs)

    if adata.shape != expected_shape:
        raise ValueError(
            f"Expected `anndata.AnnData` object to have shape `{expected_shape}`, found `{adata.shape}`."
        )

    adata.var_names_make_unique()

    return adata


@d.get_sectionsf("dataset", sections=["Parameters"])
@d.dedent
def pancreas(
    path: Union[str, Path] = "datasets/endocrinogenesis_day15.5.h5ad",
    **kwargs,
) -> AnnData:
    """
    Development of the murine pancreas at E15.5 from [Panc19]_.

    scRNA-seq dataset comprising 2531 cells recorded using 10x Chromium in a single time point. Data was filtered
    to remove heavily cycling populations and to focus on the late stages of endocrinogenesis.

    Contains raw spliced and un-spliced count data, low-dimensional embedding coordinates as well as original
    cluster annotations.

    Parameters
    ----------
    path
        Path where to save the dataset.
    **kwargs
        Keyword arguments for :func:`scanpy.read`.

    Returns
    -------
    %(adata)s
    """

    return _load_dataset_from_url(path, *_datasets["pancreas"], **kwargs)


@d.dedent
def pancreas_preprocessed(
    path: Union[str, Path] = "datasets/endocrinogenesis_day15.5_preprocessed.h5ad",
    **kwargs,
) -> AnnData:
    """
    Development of the murine pancreas at E15.5 from [Panc19]_, preprocessed according to the \
    `basic tutorial <https://cellrank.readthedocs.io/en/latest/pancreas_basic.html>`__.

    Parameters
    ----------
    %(dataset.parameters)s

    Returns
    -------
    %(adata)s
    """

    return _load_dataset_from_url(path, *_datasets["pancreas_preprocessed"], **kwargs)


@d.dedent
def lung(
    path: Union[str, Path] = "datasets/lung_regeneration.h5ad",
    **kwargs,
) -> AnnData:
    """
    Regeneration of murine lung epithelial cells at 13 time points from [Lung20]_.

    scRNA-seq dataset comprising 24,051 cells recorded using Dropseq [Macosko15]_ at 13 time points spanning days
    2-15 past lung bleomycin injury. Data was filtered to remove control cells as well as later time points which  are
    more spaced out. We wanted to focus on the densely sampled days where RNA velocity [Manno18]_ [Bergen20]_ can be
    used to predict the future cellular state.

    Contains raw spliced and un-spliced count data, low-dimensional embedding coordinates as well as original
    cluster annotations.

    Parameters
    ----------
    %(dataset.parameters)s

    Returns
    -------
    %(adata)s
    """

    return _load_dataset_from_url(path, *_datasets["lung"], **kwargs)


@inject_docs(s=ReprogrammingSubset)
@d.dedent
def reprogramming(
    subset: str = ReprogrammingSubset.FULL.s,
    path: Union[str, Path] = "datasets/reprogramming.h5ad",
    **kwargs,
) -> AnnData:
    """
    Reprogramming of mouse embryonic fibroblasts to induced endoderm progenitors at 8 time points from [Morris18]_.

    scRNA-seq dataset comprising `104,887` cell recorded using 10X Chromium and Dropseq [Macosko15]_ at 8 time points
    spanning days 0-28 past reprogramming initiation.

    Contains raw spliced and un-spliced count data, low-dimensional embedding coordinates as well as clonal information
    from CellTagging [Morris18]_. Moreover, contains the following attr:`anndata.AnnData.obs`: annotations:

        - `'reprogramming_day'` - time-point information.
        - `'reprogramming'` - whether this clone is enriched for cells from successfully reprogrammed populations.
        - `'CellTagDN_XXk'` - CellTag from day `N` from the `XXk` cells ``subset``.

    Parameters
    ---------
    subset
        Whether to return the full object or just a subset. Can be one of:

            - `{s.FULL.s!r}` - return the complete dataset containing `104,887` cells.
            - `{s.K85.s!r}` - return the subset as described in [Morris18]_ Fig. 1, containing `85,010` cells.
            - `{s.K48.s!r}` - return the subset as described in [Morris18]_ Fig. 3, containing `48,515` cells.

    %(dataset.parameters)s

    Returns
    -------
    %(adata)s

    Notes
    -----
    The dataset has approximately 1.5GiB and the subsetting is performed only locally after the full download.
    """
    subset = ReprogrammingSubset(subset)
    adata = _load_dataset_from_url(path, *_datasets["reprogramming"], **kwargs)

    if subset == ReprogrammingSubset.FULL:
        return adata
    if subset == ReprogrammingSubset.K48:
        return adata[~adata.obs["cluster"].isnull()].copy()
    if subset == ReprogrammingSubset.K85:
        return adata[~adata.obs["timecourse"].isnull()].copy()

    raise NotImplementedError(
        f"Subsetting option `{subset.s!r}` is not yet implemented."
    )
