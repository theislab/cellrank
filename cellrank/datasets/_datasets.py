# -*- coding: utf-8 -*-
"""Datasets module."""

import os
from typing import Tuple, Union, TypeVar
from pathlib import Path

from scanpy import read

from cellrank import logging as logg
from cellrank.ul._docs import d

AnnData = TypeVar("AnnData")


_datasets = {
    "pancreas": (
        "https://ndownloader.figshare.com/files/25028786?private_link=bb68d6e4686df12a0e2a",
        (2531, 27998),
    ),
    "pancreas_preprocessed": (
        "https://ndownloader.figshare.com/files/25030028?private_link=c3a8004ca127b370cc15",
        (2531, 2000),
    ),
    "lung": (
        "https://ndownloader.figshare.com/files/25038224?private_link=89d53249778c18d45e9f",
        (24882, 24051),
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
def _lung(
    path: Union[str, Path] = "datasets/lung_regeneration.h5ad",
    **kwargs,
) -> AnnData:
    """
    Regeneration of murine lung epithelial cells at 13 time points from [Lung20]_.

    scRNA-seq dataset comprising of 24,051 cells recorded using Dropseq [Macosko15]_ at 13 time points spanning days
    2-15 past lung bleomycin injury. Data was filtered to remove control cells as well as later time points which  are
    more spaced out. We wanted to focus on the densely sampled days where RNA velocity [Manno18]_ [Bergen20]_ can be
    used to predict the future cellular state.

    Contains raw spliced and un-spliced count data, low-dimensional embedding coordinates as well as  original
    cluster annotations.

    Parameters
    ----------
    %(dataset.parameters)s

    Returns
    -------
    %(adata)s
    """

    return _load_dataset_from_url(path, *_datasets["lung"], **kwargs)
