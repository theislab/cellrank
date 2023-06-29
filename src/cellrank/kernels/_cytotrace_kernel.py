import enum
from typing import Any, List, Literal, Optional, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.stats as st

from anndata import AnnData

from cellrank import logging as logg
from cellrank._utils._docs import d, inject_docs
from cellrank._utils._enum import ModeEnum
from cellrank._utils._key import Key
from cellrank._utils._utils import _correlation_test
from cellrank.kernels._pseudotime_kernel import PseudotimeKernel

__all__ = ["CytoTRACEKernel"]


class CytoTRACEAggregation(ModeEnum):
    MEAN = enum.auto()
    MEDIAN = enum.auto()
    GMEAN = enum.auto()
    HMEAN = enum.auto()


@d.dedent
class CytoTRACEKernel(PseudotimeKernel):
    """Kernel which computes directed transition probabilities using the CytoTRACE score :cite:`gulati:20`.

    .. seealso::
        - See :doc:`../../../notebooks/tutorials/kernels/400_cytotrace` on how to
          compute the :attr:`~cellrank.kernels.CytoTRACEKernel.transition_matrix` based on the CytoTRACE score.

    The k-NN graph contains information about the (undirected) connectivities among cells, reflecting their similarity.
    CytoTRACE can be used to estimate cellular plasticity and in turn, a pseudotemporal ordering of cells from more
    plastic to less plastic states. It relies on the assumption that differentiated cells express, on average,
    fewer genes than naive cells.

    This kernel internally uses the :class:`~cellrank.kernels.PseudotimeKernel` to direct the k-NN graph
    on the basis of the CytoTRACE pseudotime.

    %(density_correction)s

    Parameters
    ----------
    %(adata)s
    %(backward)s
    kwargs
        Keyword arguments for the :class:`~cellrank.kernels.PseudotimeKernel`.

    Examples
    --------
    .. code-block::

        import scvelo as scv
        import cellrank as cr

        adata = cr.datasets.pancreas()

        sc.pp.filter_genes(adata, min_cells=10)
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata)

        # CytoTRACE by default uses imputed data - a simple way to compute
        # k-NN imputed data is to use scVelo's moments function.
        # However, note that this function expects `spliced` counts because
        # it's designed for RNA velocity, so we're using a simple hack here:
        if 'spliced' not in adata.layers or 'unspliced' not in adata.layers:
            adata.layers['spliced'] = adata.X
            adata.layers['unspliced'] = adata.X
        scv.pp.moments(adata)

        ctk = cr.kernels.CytoTRACEKernel(adata)
        ckt = ctk.compute_cytotrace().compute_transition_matrix()
    """

    def __init__(
        self,
        adata: AnnData,
        backward: bool = False,
        **kwargs: Any,
    ):
        _ = kwargs.pop("time_key", None)
        super().__init__(
            adata,
            backward=backward,
            time_key=Key.cytotrace("pseudotime"),
            **kwargs,
        )

    def _read_from_adata(
        self,
        time_key: str,
        **kwargs: Any,
    ) -> None:
        try:
            super()._read_from_adata(time_key=time_key, **kwargs)
        except KeyError as e:
            self._pseudotime: Optional[np.ndarray] = None
            if "Unable to find pseudotime" not in str(e):
                raise
            super(PseudotimeKernel, self)._read_from_adata(**kwargs)

    @d.get_sections(base="cytotrace", sections=["Parameters"])
    @inject_docs(ct=CytoTRACEAggregation)
    def compute_cytotrace(
        self,
        layer: Optional[str] = "Ms",
        aggregation: Literal["mean", "median", "hmean", "gmean"] = CytoTRACEAggregation.MEAN,
        use_raw: bool = False,
        n_genes: int = 200,
    ) -> "CytoTRACEKernel":
        """Re-implementation of the CytoTRACE algorithm :cite:`gulati:20` to estimate cellular plasticity.

        Computes the number of genes expressed per cell and ranks genes according to their correlation with this
        measure. Next, it selects to top-correlating genes and aggregates their (imputed) expression to obtain
        the CytoTRACE score. A high score stands for high differentiation potential (naive, plastic cells) and
        a low score stands for low differentiation potential (mature, differentiation cells).

        Parameters
        ----------
        layer
            Key in :attr:`~anndata.AnnData.layers` or ``'X'`` for :attr:`~anndata.AnnData.X`
            from where to get the expression.
        aggregation
            How to aggregate expression of the top-correlating genes. Valid options are:

            - ``{ct.MEAN!r}`` - arithmetic mean.
            - ``{ct.MEDIAN!r}`` - median.
            - ``{ct.HMEAN!r}`` - harmonic mean.
            - ``{ct.GMEAN!r}`` - geometric mean.
        use_raw
            Whether to use the :attr:`~anndata.AnnData.raw` to compute the number of genes expressed per cell
            (#genes/cell) and the correlation of gene expression across cells with #genes/cell.
        n_genes
            Number of top positively correlated genes to compute the CytoTRACE score.

        Returns
        -------
        Nothing, just modifies :attr:`~anndata.AnnData.obs` with the following keys:

        - ``'ct_score'`` - the normalized CytoTRACE score.
        - ``'ct_pseudotime'`` - associated pseudotime, essentially ``1 - CytoTRACE score``.
        - ``'ct_num_exp_genes'`` - the number of genes expressed per cell, basis of the CytoTRACE score.

        It also modifies :attr:`~anndata.AnnData.var` with the following keys:

        - ``'ct_gene_corr'`` - the correlation as specified above.
        - ``'ct_correlates'`` - indication of the genes used to compute the CytoTRACE score, i.e. the ones that
          correlated positively with ``'ct_num_exp_genes'``.

        Notes
        -----
        This will not exactly reproduce the results of the original CytoTRACE algorithm :cite:`gulati:20` because we
        allow for any normalization and imputation techniques whereas CytoTRACE has built-in specific methods for that.
        """
        from cellrank._utils import Lineage

        aggregation = CytoTRACEAggregation(aggregation)

        if use_raw and self.adata.raw is None:
            logg.warning("`adata.raw` is `None`. Setting `use_raw=False`")
            use_raw = False
        if layer not in (None, "X") and layer not in self.adata.layers:
            raise KeyError(
                f"Unable to find `{layer!r}` in `adata.layers`. "
                f"Valid option are: `{sorted({'X'} | set(self.adata.layers.keys()))}`."
            )

        adata_mraw = self.adata.raw if use_raw else self.adata
        msg = f"Computing CytoTRACE score with `{adata_mraw.n_vars}` genes"
        if adata_mraw.n_vars < 10000:
            msg += ". Consider using more than `10000` genes"
        start = logg.info(msg)

        # compute number of expressed genes per cell
        num_exp_genes = np.asarray((adata_mraw.X > 0).sum(axis=1)).squeeze()

        logg.debug("Correlating all genes with number of genes expressed per cell")
        gene_corr = _correlation_test(
            adata_mraw.X,
            Lineage(num_exp_genes[:, None], names=["gene"]),
            gene_names=adata_mraw.var_names,
        )["gene_corr"]

        cytotrace_score, pos_top_genes = self._compute_score(
            gene_corr,
            layer=layer,
            aggregation=aggregation,
            ascending=False,
            n_genes=n_genes,
        )
        cytotrace_score -= cytotrace_score.min()
        cytotrace_score /= cytotrace_score.max()

        self.adata.obs[Key.cytotrace("score")] = cytotrace_score
        self.adata.obs[Key.cytotrace("pseudotime")] = self._pseudotime = 1 - cytotrace_score
        self.adata.obs[Key.cytotrace("num_exp_genes")] = num_exp_genes
        self.adata.var[Key.cytotrace("gene_corr")] = gene_corr
        self.adata.var[Key.cytotrace("correlates")] = False
        self.adata.var.loc[pos_top_genes, Key.cytotrace("correlates")] = True
        self.adata.uns[Key.cytotrace("params")] = {
            "aggregation": str(aggregation),
            "layer": layer,
            "use_raw": use_raw,
            "n_genes": len(pos_top_genes),
        }

        logg.info(
            f"Adding `adata.obs[{Key.cytotrace('score')!r}]`\n"
            f"       `adata.obs[{Key.cytotrace('pseudotime')!r}]`\n"
            f"       `adata.obs[{Key.cytotrace('num_exp_genes')!r}]`\n"
            f"       `adata.var[{Key.cytotrace('gene_corr')!r}]`\n"
            f"       `adata.var[{Key.cytotrace('correlates')!r}]`\n"
            f"       `adata.uns[{Key.cytotrace('params')!r}]`\n"
            f"    Finish",
            time=start,
        )

        return self

    def _compute_score(
        self,
        gene_corr: pd.Series,
        *,
        layer: Optional[str] = None,
        aggregation: CytoTRACEAggregation.MEAN,
        ascending: bool = False,
        n_genes: int = 200,
    ) -> Tuple[np.ndarray, List[str]]:
        from sklearn.utils.sparsefuncs import csc_median_axis_0

        # fmt: off
        modifier = "negatively" if ascending else "positively"
        if n_genes <= 0:
            raise ValueError(f"Expected number of {modifier} correlated genes to be positive, found `{n_genes}`.")
        top_genes = [g for g in gene_corr.sort_values(ascending=ascending).index if g in self.adata.var_names]
        top_genes = top_genes[:n_genes]
        invalid_genes = int(gene_corr.loc[top_genes].isnull().sum())
        # fmt: on

        if invalid_genes:
            raise ValueError(
                f"Top `{len(top_genes)}` {modifier} correlated genes contain `{invalid_genes}` NaN values."
            )
        if not len(top_genes):
            raise ValueError("No genes have been selected.")

        if len(top_genes) != n_genes:
            logg.warning(
                f"Unable to get requested top {modifier} correlated `{n_genes}`. " f"Using top `{len(top_genes)}` genes"
            )

        # fmt: off
        imputed_exp = self.adata[:, top_genes].X if layer == "X" else self.adata[:, top_genes].layers[layer]
        if sp.issparse(imputed_exp) and aggregation not in (CytoTRACEAggregation.MEAN, CytoTRACEAggregation.MEDIAN):
            imputed_exp = imputed_exp.A

        if aggregation == CytoTRACEAggregation.MEAN:
            cytotrace_score = np.asarray(imputed_exp.mean(axis=1)).reshape((-1,))
        elif aggregation == CytoTRACEAggregation.MEDIAN:
            if sp.issparse(imputed_exp):
                cytotrace_score = np.asarray(csc_median_axis_0(imputed_exp.T.tocsc())).reshape((-1,))
            else:
                cytotrace_score = np.median(imputed_exp, axis=1)
        elif aggregation == CytoTRACEAggregation.GMEAN:
            cytotrace_score = st.gmean(imputed_exp, axis=1)
        elif aggregation == CytoTRACEAggregation.HMEAN:
            cytotrace_score = st.hmean(imputed_exp, axis=1)
        else:
            raise NotImplementedError(f"Aggregation method `{aggregation}` is not yet implemented.")
        # fmt: on
        return cytotrace_score, top_genes
