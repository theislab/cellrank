# -*- coding: utf-8 -*-
import os
from typing import Union

from scanpy import read
from anndata import AnnData

_datasets = dict(
    pancreas=(
        "datasets/pancreas/endocrinogenesis_day15.5.h5ad",
        "https://github.com/theislab/cellrank_notebooks/raw/master/datasets/pancreas/endocrinogenesis_day15.5.h5ad",
    )
)


def _load_dataset_from_url(fpath: Union[os.PathLike, str], url: str) -> AnnData:
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
