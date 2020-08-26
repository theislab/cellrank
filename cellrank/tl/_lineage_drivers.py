# -*- coding: utf-8 -*-
"""Gene importance module."""

from types import MappingProxyType
from typing import Any, List, Tuple, Union, Mapping, TypeVar, Optional, Sequence

from statsmodels.stats.multitest import multipletests

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler

from cellrank import logging as logg
from cellrank.ul._docs import d
from cellrank.pl._utils import (
    _model_type,
    _get_backend,
    _callback_type,
    _create_models,
    _create_callbacks,
)
from cellrank.ul._utils import _get_n_cores, _check_collection
from cellrank.tl._constants import _DEFAULT_BACKEND, AbsProbKey
from cellrank.tl.estimators import GPCCA
from cellrank.ul._parallelize import parallelize
from cellrank.tl.estimators._constants import P
from cellrank.tl.kernels._precomputed_kernel import DummyKernel

AnnData = TypeVar("AnnData")
Queue = TypeVar("Queue")


def _gi_permute(
    ix: Optional[int],
    perms: Sequence[int],
    x: np.ndarray,
    y: np.ndarray,
    seed: Optional[int],
    queue: Optional[Queue],
    **kwargs,
) -> List[np.ndarray]:
    """
    Permutes the features and calculate importances using :class:`sklearn.ensemble.RandomForestRegressor`.

    Parameters
    ----------
    ix
        Offset for random ``seed``. Only used when ``seed`` is not None.
    perms
        Permutation indices. Only used to get number of permutations.
    x
        Features to permute.
    y
        Dependent variable, such as pseudotime.
    seed
        Random seed for reproducibility.
    **kwargs
        Keyword arguments for :class:`sklearn.ensemble.RandomForestRegressor`.

    Returns
    -------
    :class:`list`
        List of gene importances of each features for every permutation.
    """

    state = np.random.RandomState(None if seed is None else seed + ix)
    imps = []

    for _ in perms:
        imps.append(
            RandomForestRegressor(random_state=state.randint(2 ** 32), **kwargs)
            .fit(state.permutation(x), y)
            .feature_importances_
        )

        if queue is not None:
            queue.put(1)

    if queue is not None:
        queue.put(None)

    return imps


def _gi_process(
    genes: Sequence[str],
    models: _model_type,
    callbacks: _callback_type,
    lineage: str,
    queue: Optional[Queue],
    **kwargs,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate the smoothed gene expression which will be used to get the gene importances.

    Parameters
    ----------
    genes
        Genes for which to calculate the importances.
    models
        Gene and lineage specific models.
    callbacks
        Gene and lineage specific prepare callbacks.
    lineage
        Name of the lineage for which to calculate the gene importances.
    queue
        Signalling queue in the parent process/thread used to update the progress bar.
    **kwargs
        Keyword arguments for :meth:`cellrank.ul.models.BaseModel.prepare`.

    Returns
    -------
    :class:`numpy.ndarray`, :class:`numpy.ndarray`
        The independent variables and the optionally normalized smoothed expression.
    """

    res = []

    for gene in genes:
        cb = callbacks[gene][lineage]
        model = cb(models[gene][lineage], gene=gene, lineage=lineage, **kwargs).fit()
        res.append([model.x_test.squeeze(), model.predict()])

        if queue is not None:
            queue.put(1)

    if queue is not None:
        queue.put(None)

    res = np.array(res).swapaxes(1, 2)  # genes x points x 2
    x, y = res[..., 0], res[..., 1]

    return x, y


@d.dedent
def _gene_importance(
    adata: AnnData,
    model: _model_type,
    genes: Sequence[str],
    lineage: str,
    backward: bool = False,
    n_points: int = 200,
    time_key: str = "latent_time",
    norm: bool = True,
    n_perms: Optional[int] = None,
    fdr_correction: Optional[str] = "fdr_bh",
    alpha: float = 0.05,
    seed: Optional[int] = None,
    return_model: bool = False,
    show_progress_bar: bool = True,
    n_jobs: Optional[int] = 1,
    backend: str = _DEFAULT_BACKEND,
    callback: _callback_type = None,
    rf_kwargs: Mapping[str, Any] = MappingProxyType({"criterion": "mse"}),
    **kwargs,
) -> Union[pd.DataFrame, Tuple[RandomForestRegressor, pd.DataFrame]]:
    """
    Calculate gene importance within a lineage according to [SCORPIUS16]_.

    SCORPIUS detects potential lineage drivers by using `Random Forests` to predict the pseudotemporal ordering of cells
    given the gene expression data. We can asses how important each gene is for this prediction, which we define as the
    `'importance'` of this gene. P-values are computed via a permutation test and q-values are computed via an FDR
    correction.

    We adapted SCORPIUS to work with soft lineage assignments given by our lineage probabilities computed using
    :func:`cellrank.tl.lineages`.

    Parameters
    ----------
    %(adata)s
    genes
        Genes in ``adata.var_names``.
    lineage
        Name of the lineage for which to calculate gene importance.
    %(backward)s
    n_points
        Number of points used for prediction.
        If `None`, use original data points, not uniformly distributed across the pseudotime.
    time_key
        Key in ``adata.obs`` where the pseudotime is stored.
    norm
        Normalize each trend to `0` mean, `1` variance.
    n_perms
        Number of permutations to perform. If `None`, don't calculate the p-values.
    fdr_correction
        Method used to correct for false discovery rate.
        For available methods, see :func:`statsmodels.stats.multitest.multipletests`.
    alpha
        Family-wise error rate for FDR correction.
    seed
        Seed for :class:`sklearn.ensemble.RandomForestRegressor` and the permutations.
    return_model
        Whether to return the fitted model.
    %(parallel)s
    %(model_callback)s
    rf_kwargs
        Keyword arguments for :class:`sklearn.ensemble.RandomForestRegressor`.
    **kwargs
        Keyword arguments for :meth:`cellrank.ul.models.BaseModel.prepare`.

    Returns
    -------
    :class:`pandas.DataFrame`
        Dataframe with `'importance'` column which contains genes' importances:
            - If ``n_perm!=None``, it also contains `'pval'` column with the calculated p-values.
            - If ``fdr_correction!= None``, it also contains `'qval'` column with the corrected p-value`.
    :class:`sklearn.ensemble.RandomForestRegressor`, :class:`pandas.DataFrame`
        Same as above, but also returns the fitted model.
    """

    ln_key = str(AbsProbKey.BACKWARD if backward else AbsProbKey.FORWARD)
    if ln_key not in adata.obsm:
        raise KeyError(f"Lineage key `{ln_key!r}` not found in `adata.obsm`.")

    _ = adata.obsm[ln_key][lineage]

    if n_perms is not None:
        if n_perms < 0:
            raise ValueError(
                f"Number of permutations must be `>= 0`, found `{n_perms}`."
            )
    _check_collection(adata, genes, "var_names", use_raw=kwargs.get("use_raw", False))

    n_jobs = _get_n_cores(n_jobs, len(genes))

    kwargs["backward"] = backward
    kwargs["time_key"] = time_key
    kwargs["n_test_points"] = n_points

    models = _create_models(model, genes, [lineage])
    callbacks = _create_callbacks(adata, callback, genes, [lineage])
    backend = _get_backend(model, backend)

    start = logg.info(f"Calculating gene trends using `{n_jobs}` core(s)")
    data = parallelize(
        _gi_process,
        genes,
        n_jobs=n_jobs,
        unit="gene",
        as_array=False,
        extractor=np.hstack,
        backend=backend,
        show_progress_bar=show_progress_bar,
    )(models, callbacks, lineage, **kwargs).T
    logg.info("    Finish", time=start)

    x, pt = data[..., 1], data[:, 0, 0]
    if np.all(pt[..., np.newaxis] != data[..., 0]):
        raise RuntimeError("Sanity check failed: pseudotime differs for genes.")
    if np.all(np.sort(pt) != pt):
        raise RuntimeError("Sanity check failed: pseudotime is not sorted.")

    if norm:
        _ = StandardScaler(copy=False).fit_transform(x)

    rf_kwargs = dict(rf_kwargs)
    if rf_kwargs.get("n_jobs", None) is None:
        rf_kwargs["n_jobs"] = n_jobs

    logg.debug("Fitting random forest")
    model = RandomForestRegressor(random_state=seed, **rf_kwargs).fit(x, pt)

    importances = pd.DataFrame(
        model.feature_importances_.T, index=genes, columns=["importance"]
    )  # shape: n_genex x 1

    if n_perms is None or n_perms == 0:
        importances.sort_values("importance", ascending=False, inplace=True)
        return (importances, model) if return_model else importances

    # since we use most of the threads for permutations, don't want to use many of them
    rf_kwargs["n_jobs"] = 4

    # shape: n_perms x n_genes
    start = logg.info(f"Running permutation test using `{n_jobs}` core(s)")
    perms = parallelize(
        _gi_permute,
        list(range(n_perms)),
        unit="perm",
        n_jobs=n_jobs,
        backend=backend,
        extractor=np.vstack,
        use_ixs=True,
        show_progress_bar=show_progress_bar,
    )(x, pt, seed, **rf_kwargs)
    if perms.shape != (n_perms, len(genes)):
        raise RuntimeError("Sanity check failed: the number of permutations differ.")
    logg.info("    Finish", time=start)

    importances["pval"] = np.mean(importances.values < np.ravel(perms), axis=1)
    sort_by = "pval"

    if fdr_correction is not None:
        _, qvals, _, _ = multipletests(
            importances["pval"].values,
            alpha=alpha,
            method=fdr_correction,
            is_sorted=False,
            returnsorted=False,
        )
        importances["qval"] = qvals
        sort_by = "qval"

    importances.sort_values(
        [sort_by, "importance"], ascending=[True, False], inplace=True
    )

    return (importances, model) if return_model else importances


@d.dedent
def lineage_drivers(
    adata: AnnData,
    backward: bool = False,
    lineages: Optional[Union[Sequence, str]] = None,
    cluster_key: Optional[str] = None,
    clusters: Optional[Union[Sequence, str]] = None,
    layer: str = "X",
    use_raw: bool = True,
    return_drivers: bool = False,
) -> Optional[pd.DataFrame]:  # noqa
    """
    %(lineage_drivers.full_desc)s

    Parameters
    ----------
    %(adata)s
    %(backward)s
    %(lineage_drivers.parameters)s

    Returns
    -------
    %(lineage_drivers.returns)s
    """

    # create dummy kernel and estimator
    pk = DummyKernel(adata, backward=backward)
    g = GPCCA(pk, read_from_adata=True)
    if g._get(P.ABS_PROBS) is None:
        raise RuntimeError(
            "Compute absorption probabilities first as `cellrank.tl.lineages()`."
        )

    # call the underlying function to compute and store the lineage drivers
    return g.compute_lineage_drivers(
        lineages=lineages,
        cluster_key=cluster_key,
        clusters=clusters,
        layer=layer,
        use_raw=use_raw,
        return_drivers=return_drivers,
    )
