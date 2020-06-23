# -*- coding: utf-8 -*-
"""Gene importance module."""

from types import MappingProxyType
from typing import Any, Dict, List, Tuple, Union, Mapping, Optional, Sequence

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from statsmodels.stats.multitest import multipletests

from scanpy import logging as logg
from anndata import AnnData

from cellrank.tools import GPCCA
from cellrank.utils._utils import _get_n_cores, check_collection
from cellrank.utils.models import Model
from cellrank.tools.kernels import ConnectivityKernel
from cellrank.plotting._utils import _model_type, _create_models, _is_any_gam_mgcv
from cellrank.tools._constants import LinKey
from cellrank.utils._parallelize import parallelize


def _gi_permute(
    ix: Optional[int],
    perms: Sequence[int],
    x: np.ndarray,
    y: np.ndarray,
    seed: Optional[int],
    queue,
    **kwargs,
) -> List[np.ndarray]:
    """
    Permutes the features and calculate importances using :class:`sklearn.ensemble.RandomForestRegressor`.

    Params
    ------
    ix
        Offset for random :paramref:`seed`. Only used when :paramref:`seed` is not None.
    perms
        Permutation indices. Only used to get number of permutations.
    x
        Features to permute.
    y
        Dependent variable, such as pseudotime.
    seed
        Random seed for reproducibility.
    kwargs
        Keyword arguments for :class:`sklearn.ensemble.RandomForestRegressor`.

    Returns
    -------
    :class:`list`
        List of importances of each features for every permutation.
    """

    state = np.random.RandomState(None if seed is None else seed + ix)
    imps = []

    for _ in perms:
        imps.append(
            RandomForestRegressor(random_state=state.randint(2 ** 32), **kwargs)
            .fit(state.permutation(x), y)
            .feature_importances_
        )
        queue.put(1)
    queue.put(None)

    return imps


def _gi_process(
    genes: Sequence[str],
    models: Dict[str, Dict[str, Model]],
    lineage_name: str,
    norm: bool,
    queue,
    **kwargs,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate the smoothed gene expression which will be used to get the gene importances.

    Params
    ------
    genes
        Genes for which to calculate the importances.
    models
        Gene and lineage specific models.
    lineage_name
        Name of the lineage for which to calculate the gene importances.
    queue
        Signalling queue in the parent process/thread used to update the progress bar.
    kwargs
        Keyword arguments for :meth:`cellrank.ul.models.Model.prepare`.

    Returns
    -------
    :class:`numpy.ndarray`, :class:`numpy.ndarray`
        The independent variables and the optionally normalized smoothed expression.
    """

    res = []

    for gene in genes:
        model = models[gene][lineage_name].prepare(gene, lineage_name, **kwargs).fit()
        res.append([model.x_test.squeeze(), model.predict()])
        queue.put(1)
    queue.put(None)  # sentinel

    res = np.array(res).swapaxes(1, 2)  # genes x points x 2
    x, res = res[..., 0], res[..., 1]

    if not norm:
        return x, res

    mean = np.expand_dims(np.mean(res, 1), -1)
    sd = np.expand_dims(np.sqrt(np.var(res, 1)), -1)

    return x, (res - mean) / sd


def gene_importance(
    adata: AnnData,
    model: _model_type,
    genes: Sequence[str],
    lineage_name: str,
    n_points: int = 200,
    time_key: str = "latent_time",
    final: bool = True,
    norm: bool = True,
    n_perms: Optional[int] = None,
    fdr_correction: Optional[str] = "fdr_bh",
    alpha: float = 0.05,
    n_jobs: Optional[int] = 1,
    seed: Optional[int] = None,
    return_model: bool = False,
    backend: str = "multiprocessing",
    show_progress_bar: bool = True,
    rf_kwargs: Mapping[str, Any] = MappingProxyType({"criterion": "mse"}),
    **kwargs,
) -> Union[pd.DataFrame, Tuple[RandomForestRegressor, pd.DataFrame]]:
    """
    Calculate gene importance within a lineage according to [SCORPIUS16]_.

    SCORPIUS detects potential lineage drivers by using `Random Forests` to predict the pseudotemporal ordering of cells
    given the gene expression data. We can asses how important each gene is for this prediction, which we define as the
    'importance' of this gene. p-values are computed via a permutation test and q-values are computed via an FDR
    correction.

    We adapted SCORPIUS to work with soft lineage assignments given by our lineage probabilities computed using
    :func:`cellrank.tl.lineages`.

    Params
    ------
    adata : :class:`anndata.AnnData`
        Annotated data object.
    genes
        Genes in :paramref:`adata` `.var_names`.
    lineage_name
        Name of the lineage for which to calculate gene importance.
    n_points
        Number of points used for prediction.

        If `None`, use original data points, not uniformly distributed across the pseudotime.
    final
        Whether to consider cells going to final states or vice versa.
    time_key
        Key in :paramref:`adata` `.obs` where the pseudotime is stored.
    norm
        Normalize each trend to `0` mean, `1` variance.
    n_perms
        Number of permutations to perform.

        Use `None` if you don't want to calculate the p-values.
    fdr_correction
        Method used to correct for false discovery rate.

        For available methods, see :func:`statsmodels.stats.multitest.multipletests`.
    alpha
        Family-wise error rate for FDR correction.
    n_jobs
        Number of parallel jobs. If `-1`, use all available cores. If `None` or `1`, the execution is sequential.
    seed
        Seed for :class:`sklearn.ensemble.RandomForestRegressor` and the permutations.
    return_model
        Whether to also return the fitted model.
    backend
        Which backend to use for multiprocessing.
        See :class:`joblib.Parallel` for valid options.
    show_progress_bar
        Whether to show a progress bar tracking models fitted.
    rf_kwargs
        Keyword arguments for :class:`sklearn.ensemble.RandomForestRegressor`.
    **kwargs:
        Keyword arguments for :meth:`cellrank.ul.models.Model.prepare`.

    Returns
    -------
    :class:`pandas.DataFrame`
        Dataframe with `'importance'` column which contains genes' importances.

        - If :paramref:`n_perm` `!=None`, it also contains `'pval'` columns with calculated `p-values`.
        - If :paramref:`fdr_correction` `!= None`, it also contains `'qval'` column, containing the
          corrected `p-values`.
    (:class:`sklearn.ensemble.RandomForestRegressor`, :class:`pandas.DataFrame`)
        Same as above, but also returns the fitted model.
    """

    ln_key = str(LinKey.FORWARD if final else LinKey.BACKWARD)
    if ln_key not in adata.obsm:
        raise KeyError(f"Lineages key `{ln_key!r}` not found in `adata.obsm`.")

    _ = adata.obsm[ln_key][lineage_name]

    if n_perms is not None:
        if n_perms < 0:
            raise ValueError(
                f"Number of permutations must be `>= 0`, found `{n_perms}`."
            )
    check_collection(adata, genes, "var_names")

    n_jobs = _get_n_cores(n_jobs, len(genes))

    kwargs["final"] = final
    kwargs["time_key"] = time_key
    kwargs["n_test_points"] = n_points

    models = _create_models(model, genes, [lineage_name])
    if _is_any_gam_mgcv(models):
        logg.debug(
            "DEBUG: Setting backend to multiprocessing because model is `GamMGCV`"
        )
        backend = "multiprocessing"

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
    )(models, lineage_name, norm, **kwargs).T
    logg.info("    Finish", time=start)

    x, y = data[..., 1], data[:, 0, 0]
    if np.all(y[..., np.newaxis] != data[..., 0]):
        raise RuntimeError("Sanity check failed: pseudotime differs for genes.")
    if np.all(np.sort(y) != y):
        raise RuntimeError("Sanity check failed: pseudotime is not sorted.")

    rf_kwargs = dict(rf_kwargs)
    if rf_kwargs.get("n_jobs", None) is None:
        rf_kwargs["n_jobs"] = n_jobs

    logg.debug("DEBUG: Running random forest")
    model = RandomForestRegressor(random_state=seed, **rf_kwargs).fit(x, y)

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
    )(x, y, seed, **rf_kwargs)
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


def lineage_drivers(
    adata: AnnData,
    final: bool = True,
    lin_names: Optional[Union[Sequence, str]] = None,
    cluster_key: Optional[str] = None,
    clusters: Optional[Union[Sequence, str]] = None,
    layer: str = "X",
    use_raw: bool = True,
    inplace: bool = True,
):
    """
    Compute driver genes per lineage.

    Correlates gene expression with lineage probabilities, for a given lineage and set of clusters.
    Often, it makes sense to restrict this to a set of clusters which are relevant
    for the lineage under consideration.

    Params
    ------
    adata
        Annodated data matrix
    final
        If True, use forward process, else backward
    lin_names
        Either a set of lineage names from :paramref:`lineage_probabilities` `.names` or None,
        in which case all lineages are considered.
    cluster_key
        Key from :paramref:`adata` `.obs` to obtain cluster annotations.
        These are considered for :paramref:`clusters`. Default is `"clusters"` if a list of `clusters` is given.
    clusters
        Restrict the correlations to these clusters.
    layer
        Key from :paramref:`adata` `.layers`.
    use_raw
        Whether or not to use :paramref:`adata` `.raw` to correlate gene expression.
        If using a layer other than `.X`, this must be set to `False`.

    Returns
    -------
    :class:`pandas.DataFrame`, :class:`NoneType`
        Writes to :paramref:`adata` `.var` or :paramref:`adata` `.raw.var`,
        depending on the value of :paramref:`use_raw`.
        For each lineage specified, a key is added to `.var` and correlations are saved there.

        Returns `None` if :paramref:`inplace` `=True`, otherwise a :class:`pandas.DataFrame`.
    """

    # create dummy kernel and estimator
    kernel = ConnectivityKernel(adata, backward=not final)
    g = GPCCA(kernel)
    g._lin_probs = adata.obsm[g._lin_key]

    # call the underlying function to compute and store the lineage drivers
    g.compute_lineage_drivers(
        lin_names=lin_names,
        cluster_key=cluster_key,
        clusters=clusters,
        layer=layer,
        use_raw=use_raw,
        inplace=inplace,
    )
