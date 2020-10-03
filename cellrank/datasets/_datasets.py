# -*- coding: utf-8 -*-
"""Datasets module."""

import os
from typing import Union, TypeVar
from pathlib import Path

from scanpy import read

from cellrank import logging as logg
from cellrank.ul._docs import d

AnnData = TypeVar("AnnData")

_datasets = dict(  # noqa
    pancreas=(
        "datasets/endocrinogenesis_day15.5.h5ad",
        "https://github.com/theislab/cellrank_notebooks/raw/master/datasets/pancreas/endocrinogenesis_day15.5.h5ad",
    ),
    lung=(
        "datasets/lung_regeneration.h5ad",
        "https://github.com/theislab/cellrank_notebooks/raw/master/datasets/lung/regeneration.h5ad",
    ),
    pancreas_preprocessed=(
        None,
        "https://github.com/theislab/cellrank_notebooks/raw/master/"
        "datasets/pancreas/endocrinogenesis_day15.5_preprocessed.h5ad",
    ),
)


def _load_dataset_from_url(fpath: Union[os.PathLike, str], url: str) -> AnnData:
    if os.path.isfile(fpath):
        logg.debug(f"Loading dataset from `{fpath!r}`")
    else:
        logg.debug(f"Downloading dataset from `{url!r}` as `{fpath!r}`")

    dirname, _ = os.path.split(fpath)
    try:
        logg.debug(f"Creating directory `{dirname!r}`")
        os.makedirs(dirname, exist_ok=True)
    except OSError as e:
        logg.debug(f"Unable to create directory `{dirname!r}`. Reason `{e}`")

    adata = read(fpath, backup_url=url, sparse=True, cache=True)
    adata.var_names_make_unique()

    return adata


@d.dedent
def pancreas() -> AnnData:
    """
    Development of the murine pancreas at E15.5 from [Panc19]_.

    scRNA-seq dataset comprising 2531 cells recorded using 10x Chromium in a single time point. Data was filtered
    to remove heavily cycling populations and to focus on the late stages of endocrinogenesis. Contains raw spliced and
    un-spliced count data, low-dimensional embedding coordinates as well as original cluster annotations.

    Returns
    -------
    %(adata)s
    """

    return _load_dataset_from_url(*_datasets["pancreas"])


@d.dedent
def lung() -> AnnData:
    """
    Regeneration of murine lung epithelial cells at 13 time points from [Lung20]_.

    scRNA-seq dataset comprising 24,051 cells recorded using Dropseq CITE at 13 time points spanning days 2-15 past
    lung bleomycin injury. Data was filtered to remove control cells as well as later time points which are more spaced
    out. We wanted to focus on the densely sampled days where RNA velocity CITE can be used to predict the future
    cellular state. Contains raw spliced and un-spliced count data, low-dimensional embedding coordinates as well as
    original cluster annotations.


    Returns
    -------
    %(adata)s
    """

    return _load_dataset_from_url(*_datasets["lung"])


@d.dedent
def pancreas_preprocessed(
    path: Union[str, Path] = "datasets/endocrinogenesis_day15.5_preprocessed.h5ad"
) -> AnnData:
    """
    Development of the murine pancreas at E15.5 from [Panc19]_, preprocessed according to the \
    `basic tutorial <https://cellrank.readthedocs.io/en/latest/pancreas_basic.html>`__.

    Parameters
    ----------
    path
        Path where to save the dataset.

    Returns
    -------
    %(adata)s
    """

    _, url = _datasets["pancreas_preprocessed"]
    return _load_dataset_from_url(path, url)
