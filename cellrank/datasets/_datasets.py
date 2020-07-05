# -*- coding: utf-8 -*-
"""Datasets module."""

import os
from typing import Union, TypeVar

from scanpy import read

from cellrank import logging as logg

AnnData = TypeVar("AnnData")

_datasets = dict(  # noqa
    pancreas=(
        "datasets/endocrinogenesis_day15.5.h5ad",
        "https://github.com/theislab/cellrank_notebooks/raw/master/datasets/pancreas/endocrinogenesis_day15.5.h5ad",
    )
)


def _load_dataset_from_url(fpath: Union[os.PathLike, str], url: str) -> AnnData:
    if os.path.isfile(fpath):
        logg.debug(f"Loading dataset from `{fpath!r}`")
    else:
        logg.debug(f"Downloading dataset from `{url!r}`")

    adata = read(fpath, backup_url=url, sparse=True, cache=True)
    adata.var_names_make_unique()

    return adata


def pancreas() -> AnnData:
    """
    Pancreatic endocrinogenesis from [Panc19]_.

    Note that we subsetted the original data to focus on endocrine development downstream of the Ngn3 low EP cluster,
    i.e. we only consider cells that have high probability of becoming endocrine.

    Returns
    -------
    :class:`anndata.AnnData`
        The annotated data object.
    """

    return _load_dataset_from_url(*_datasets["pancreas"])
