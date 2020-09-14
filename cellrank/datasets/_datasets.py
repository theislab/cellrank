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
    Pancreatic endocrinogenesis from [Panc19]_.

    Note that we subset the original data to focus on endocrine development downstream of the Ngn3 low EP cluster,
    i.e. we only consider cells that have high probability of becoming endocrine.

    Returns
    -------
    %(adata)s
    """

    return _load_dataset_from_url(*_datasets["pancreas"])


@d.dedent
def pancreas_preprocessed(
    path: Union[str, Path] = "datasets/endocrinogenesis_day15.5_preprocessed.h5ad"
) -> AnnData:
    """
    Pancreatic endocrinogenesis from [Panc19]_.

    Note that we subset the original data to focus on endocrine development downstream of the Ngn3 low EP cluster,
    i.e. we only consider cells that have high probability of becoming endocrine.

    This dataset has been preprocessed according to the basic tutorial.

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
