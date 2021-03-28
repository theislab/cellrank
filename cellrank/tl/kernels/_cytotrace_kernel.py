from typing import Any

from typing_extensions import Literal

from anndata import AnnData

import numpy as np
from scipy.stats import gmean, hmean

from cellrank import logging as logg
from cellrank.ul._docs import d
from cellrank.tl._utils import _correlation_test_helper
from cellrank.tl._constants import ModeEnum
from cellrank.tl.kernels._pseudotime_kernel import PseudotimeKernel


def _ct(key: str) -> str:
    return f"ct_{key}"


class CytoTRACEAggregation(ModeEnum):  # noqa
    MEAN = "mean"
    MEDIAN = "median"
    GMEAN = "gmean"
    HMEAN = "hmean"


@d.dedent
class CytoTRACEKernel(PseudotimeKernel):
    """
    Kernel which computes directed transition probabilities based on a KNN graph and the CytoTRACE score [Cyto20]_.

    The KNN graph contains information about the (undirected) connectivities among cells, reflecting their similarity.
    CytoTRACE can be used to estimate cellular plasticity and in turn, a pseudotemporal ordering of cells from more
    plastic to less plastic states.
    This kernel internally uses the :class:`cellrank.tl.kernels.PseudotimeKernel` to direct the KNN graph
    on the basis of the CytoTRACE-derived pseudotime.

    %(density_correction)s

    Parameters
    ----------
    %(adata)s
    %(backward)s
    %(cytotrace.parameters)s

    compute_cond_num
        Whether to compute condition number of the transition matrix. Note that this might be costly,
        since it does not use sparse implementation.

    Examples
    --------
    Workflow::

        import scvelo as scv
        import cellrank as cr

        adata = cr.datasets.pancreas()

        sc.pp.filter_genes(adata, min_cells=10)
        adata.raw = adata.copy()
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata)

        if 'spliced' not in adata.layers or 'unspliced' not in adata.layers:
            # use the following trick to get scvelo's moments function working
            adata.layers['spliced'] = adata.X
            adata.layers['unspliced'] = adata.X

        scv.pp.moments(adata, n_pcs=None, n_neighbors=None)
    """

    def __init__(
        self,
        adata: AnnData,
        backward: bool = False,
        layer: str = "Ms",
        aggregation: Literal["mean", "median", "hmean", "gmean"] = "mean",
        use_raw: bool = False,
        compute_cond_num: bool = False,
        check_connectivity: bool = False,
    ):
        super().__init__(
            adata,
            backward=backward,
            time_key=_ct("pseudotime"),
            compute_cond_num=compute_cond_num,
            check_connectivity=check_connectivity,
            layer=layer,
            aggregation=aggregation,
            use_raw=use_raw,
        )
        self._time_key = _ct("pseudotime")  # quirk or PT kernel

    def _read_from_adata(
        self,
        time_key: str,
        layer: str = "Ms",
        aggregation: Literal["mean", "median", "hmean", "gmean"] = "mean",
        use_raw: bool = True,
        **kwargs: Any,
    ):
        self.compute_cytotrace(layer=layer, aggregation=aggregation, use_raw=use_raw)

        super()._read_from_adata(time_key=time_key, **kwargs)

    @d.get_sections(base="cytotrace", sections=["Parameters"])
    def compute_cytotrace(
        self,
        layer: str = "Ms",
        aggregation: Literal["mean", "median", "hmean", "gmean"] = "mean",
        use_raw: bool = False,
    ) -> None:
        """
        Re-implementation of the CytoTRACE algorithm [Cyto20]_ to estimate cellular plasticity.

        Computes the number of genes expressed per cell and ranks genes according to their correlation with this
        measure. Next, it selects to top-correlating genes and aggregates their (imputed) expression to obtain
        the CytoTRACE score. A high score stands for high differentiation potential (naive, plastic cells) and
        a low score stands for low differentiation potential (mature, differentiation cells).

        Note that this will not exactly reproduce the results of the original CytoTRACE algorithm [Cyto20]_ because we
        allow for any normalization and imputation techniques whereas CytoTRACE has build-in specific methods for that.

        Parameters
        ----------
        layer
            Key in :attr:`anndata.AnnData.layers` or `'X'` for :attr:`anndata.AnnData.X`
            from where to get the expression.
        aggregation
            How to aggregate expression of the top-correlating genes. Valid options are:

                - `'mean'`: arithmetic mean.
                - `'median'`: median.
                - `'gmean'`: geometric mean.
                - `'hmean'`: harmonic mean.

        use_raw
            Whether to use the :attr:`anndata.AnnData.raw` to compute the number of genes expressed per cell
            (#genes/cell) and the correlation of gene expression across cells with #genes/cell.

        Returns
        -------
        Nothing, just modifies :attr:`anndata.AnnData.obs` with the following keys:

            - `'ct_score'`: the normalized CytoTRACE score.
            - `'ct_pseudotime'`: associated pseudotime, essentially `1 - CytoTRACE score`.
            - `'ct_num_exp_genes'`: the number of genes expressed per cell, basis of the CytoTRACE score.

        It also modifies :attr:`anndata.AnnData.var` with the following keys:

            - `'ct_gene_corr'`: the correlation as specified above.
            - `'ct_correlates'`: indication of the genes used to compute the CytoTRACE score, i.e. the ones that
              correlated best with `'num_exp_genes'`.
        """
        # check use_raw
        aggregation = CytoTRACEAggregation(aggregation)
        if use_raw and self.adata.raw is None:
            logg.warning("`adata.raw` is `None`. Setting `use_raw=False`")
            use_raw = False
        if use_raw and self.adata.raw.n_vars != self.adata.n_vars:
            logg.warning(
                f"`adata.raw` has different number of genes ({self.adata.raw.n_vars}) "
                f"than `adata` ({self.adata.n_vars}). Setting `use_raw=False`"
            )
            use_raw = False

        adata_mraw = self.adata.raw if use_raw else self.adata
        if layer != "X" and layer not in self.adata.layers:
            raise KeyError(
                f"Unable to find `{layer!r}` in `adata.layers`. "
                f"Valid option are: `{sorted({'X'} | set(self.adata.layers.keys()))}`."
            )

        msg = f"Computing CytoTRACE score with `{self.adata.n_vars}` genes"
        if self.adata.n_vars < 10000:
            msg += ". Consider using more than `10000` genes"
        start = logg.info(msg)

        # compute number of expressed genes per cell
        logg.debug(
            f"Computing number of genes expressed per cell with `use_raw={use_raw}`"
        )
        num_exp_genes = np.array((adata_mraw.X > 0).sum(axis=1)).reshape(-1)
        self.adata.obs[_ct("num_exp_genes")] = num_exp_genes

        # fmt: off
        # compute correlation with all genes
        logg.debug("Correlating all genes with number of genes expressed per cell")
        gene_corr, _, _, _ = _correlation_test_helper(adata_mraw.X.T, num_exp_genes[:, None])

        # annotate the top 200 genes in terms of correlation
        logg.debug("Finding the top `200` most correlated genes")
        self.adata.var[_ct("gene_corr")] = gene_corr
        top_200 = self.adata.var.sort_values(by=_ct("gene_corr"), ascending=False).index[:200]
        self.adata.var[_ct("correlates")] = False
        self.adata.var.loc[top_200, _ct("correlates")] = True

        # compute mean/median over top 200 genes, aggregate over genes and shift to [0, 1] range
        logg.debug(f"Aggregating imputed gene expression using aggregation `{aggregation}` in layer `{layer}`")
        corr_mask = self.adata.var[_ct("correlates")]
        imputed_exp = self.adata[:, corr_mask].X if layer == "X" else self.adata[:, corr_mask].layers[layer]

        # aggregate across the top 200 genes
        if aggregation == CytoTRACEAggregation.MEAN:
            cytotrace_score = np.mean(imputed_exp, axis=1)
        elif aggregation == CytoTRACEAggregation.MEDIAN:
            cytotrace_score = np.median(imputed_exp, axis=1)
        elif aggregation == CytoTRACEAggregation.GMEAN:
            cytotrace_score = gmean(imputed_exp, axis=1)
        elif aggregation == CytoTRACEAggregation.HMEAN:
            cytotrace_score = hmean(imputed_exp, axis=1)
        else:
            raise NotImplementedError(f"Aggregation method `{aggregation}` is not yet implemented.")
        # fmt: on

        # scale to 0-1 range
        cytotrace_score -= np.min(cytotrace_score)
        cytotrace_score /= np.max(cytotrace_score)
        self.adata.obs[_ct("score")] = cytotrace_score
        self.adata.obs[_ct("pseudotime")] = 1 - cytotrace_score

        self.adata.uns[_ct("params")] = {
            "aggregation": aggregation.s,
            "layer": layer,
            "use_raw": use_raw,
        }

        logg.info(
            f"Adding `adata.obs[{_ct('score')!r}]`\n"
            f"       `adata.obs[{_ct('pseudotime')!r}]`\n"
            f"       `adata.obs[{_ct('num_exp_genes')!r}]`\n"
            f"       `adata.var[{_ct('gene_corr')!r}]`\n"
            f"       `adata.var[{_ct('correlates')!r}]`\n"
            f"       `adata.uns[{_ct('params')!r}]`\n"
            f"    Finish",
            time=start,
        )
