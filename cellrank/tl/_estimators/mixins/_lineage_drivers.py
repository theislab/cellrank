from typing import Any, Dict, Tuple, Union, Mapping, Optional, Protocol, Sequence

from pathlib import Path
from datetime import datetime

import scanpy as sc
import scvelo as scv
from anndata import AnnData
from cellrank import logging as logg
from cellrank.tl import Lineage
from cellrank.ul._docs import d, inject_docs
from cellrank.tl._utils import TestMethod, save_fig, _correlation_test
from cellrank.tl._colors import _create_categorical_colors
from cellrank.tl._estimators.mixins._constants import Key
from cellrank.tl._estimators.mixins._absorption_probabilities import AbsProbsMixin

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import rc_context, patheffects
from matplotlib.axes import Axes
from matplotlib.patches import ArrowStyle


class LinDriversProtocol(Protocol):
    def _write_lineage_drivers(
        self, names: Sequence[str], use_raw: bool, time: datetime
    ) -> None:
        ...

    @property
    def adata(self) -> AnnData:
        ...

    @property
    def backward(self) -> bool:
        ...

    @property
    def eigendecomposition(self) -> Dict[str, Any]:
        ...

    @property
    def absorption_probabilities(self) -> Optional[Lineage]:
        ...

    @property
    def lineage_drivers(self) -> Optional[pd.DataFrame]:
        ...

    def _set(
        self,
        attr: str,
        obj: Optional[Union[pd.DataFrame, Mapping[str, Any]]] = None,
        key: Optional[str] = None,
        value: Optional[Union[np.ndarray, pd.Series, pd.DataFrame, Lineage]] = None,
    ):
        ...


class LinDriversMixin(AbsProbsMixin):
    def __init__(self, **kwargs: Any):
        super().__init__(**kwargs)
        self._lineage_drivers: Optional[pd.DataFrame] = None

    @d.dedent
    @inject_docs(tm=TestMethod)
    def compute_lineage_drivers(
        self: LinDriversProtocol,
        lineages: Optional[Union[str, Sequence]] = None,
        method: str = TestMethod.FISCHER.s,
        cluster_key: Optional[str] = None,
        clusters: Optional[Union[str, Sequence]] = None,
        layer: Optional[str] = None,
        use_raw: bool = False,
        confidence_level: float = 0.95,
        n_perms: int = 1000,
        seed: Optional[int] = None,
        return_drivers: bool = True,
        **kwargs: Any,
    ) -> Optional[pd.DataFrame]:
        """
        Compute driver genes per lineage.

        Correlates gene expression with lineage probabilities, for a given lineage and set of clusters.
        Often, it makes sense to restrict this to a set of clusters which are relevant for the specified lineages.

        Parameters
        ----------
        lineages
            Set of lineage names from :attr:`absorption_probabilities`. If `None`, use all lineages.
        method
            Mode to use when calculating p-values and confidence intervals. Valid options are:

                - `{tm.FISCHER.s!r}` - use Fischer transformation :cite:`fischer:21`.
                - `{tm.PERM_TEST.s!r}` - use permutation test.
        cluster_key
            Key from :attr:`anndata.AnnData.obs` to obtain cluster annotations. These are considered for ``clusters``.
        clusters
            Restrict the correlations to these clusters.
        layer
            Key from :attr:`anndata.AnnData.layers` from which to get the expression.
            If `None` or `'X'`, use :attr:`anndata.AnnData.X`.
        use_raw
            Whether or not to use :attr:`anndata.AnnData.raw` to correlate gene expression.
        confidence_level
            Confidence level for the confidence interval calculation. Must be in interval `[0, 1]`.
        n_perms
            Number of permutations to use when ``method = {tm.PERM_TEST.s!r}``.
        seed
            Random seed when ``method = {tm.PERM_TEST.s!r}``.
        return_drivers
            TODO: remove me
            Whether to return the drivers. This also contains the lower and upper ``confidence_level`` bounds.
        %(parallel)s

        Returns
        -------
        %(correlation_test.returns)s
        Only if ``return_drivers = True``.

        TODO: use varm
        Otherwise, updates :attr:`anndata.AnnData.varm` or :attr:`anndata.AnnData.raw.varm`, depending ``use_raw``:

            - ``['{{direction}}_{{lineage}}'] ['corr']`` - the potential lineage drivers.
            - ``['{{direction}}_{{lineage}}'] ['qval']`` - the corrected p-values.

        Also updates the following fields:

            - :attr:`lineage_drivers` - the same as described above.
        """

        # check that lineage probs have been computed
        method = TestMethod(method)
        abs_probs = self.absorption_probabilities
        if abs_probs is None:
            raise RuntimeError(
                "Compute `.absorption_probabilities` first as `.compute_absorption_probabilities()`."
            )

        if abs_probs.shape[1] == 1:
            logg.warning(
                "There is only 1 lineage present. Using stationary distribution instead"
            )
            stat_dist = self.eigendecomposition.get("stationary_dist", None)
            if stat_dist is None:
                raise RuntimeError(
                    "No stationary distribution found in `.eigendecomposition['stationary_dist']`."
                )
            abs_probs = Lineage(
                stat_dist,
                names=abs_probs.names,
                colors=abs_probs.colors,
            )

        # check all lin_keys exist in self.lin_names
        if isinstance(lineages, str):
            lineages = [lineages]
        if lineages is not None:
            _ = abs_probs[lineages]
        else:
            lineages = abs_probs.names

        if not len(lineages):
            raise ValueError("No lineages have been selected.")

        # use `cluster_key` and clusters to subset the data
        if clusters is not None:
            if cluster_key not in self.adata.obs.keys():
                raise KeyError(f"Key `{cluster_key!r}` not found in `adata.obs`.")
            if isinstance(clusters, str):
                clusters = [clusters]

            all_clusters = np.array(self.adata.obs[cluster_key].cat.categories)
            cluster_mask = np.array([name not in all_clusters for name in clusters])

            if any(cluster_mask):
                raise KeyError(
                    f"Clusters `{list(np.array(clusters)[cluster_mask])}` not found in "
                    f"`adata.obs[{cluster_key!r}]`."
                )
            subset_mask = np.in1d(self.adata.obs[cluster_key], clusters)
            adata_comp = self.adata[subset_mask]
            lin_probs = abs_probs[subset_mask, :]
        else:
            adata_comp = self.adata
            lin_probs = abs_probs

        # check that the layer exists, and that use raw is only used with layer X
        if layer in (None, "X"):
            if use_raw and self.adata.raw is None:
                logg.warning("No raw attribute set. Using `.X` instead")
                use_raw = False
            data = adata_comp.raw.X if use_raw else adata_comp.X
            var_names = adata_comp.raw.var_names if use_raw else adata_comp.var_names
        else:
            if layer not in self.adata.layers:
                raise KeyError(f"Layer `{layer!r}` not found in `adata.layers`.")
            if use_raw:
                logg.warning("If `use_raw=True`, layer must be `None`.")
                use_raw = False
            data = adata_comp.layers[layer]
            var_names = adata_comp.var_names

        start = logg.info(
            f"Computing correlations for lineages `{lineages}` restricted to clusters `{clusters}` in "
            f"layer `{'X' if layer is None else layer}` with `use_raw={use_raw}`"
        )

        lin_probs = lin_probs[lineages]
        drivers = _correlation_test(
            data,
            lin_probs,
            gene_names=var_names,
            method=method,
            n_perms=n_perms,
            seed=seed,
            confidence_level=confidence_level,
            **kwargs,
        )
        self._write_lineage_drivers(drivers, use_raw=use_raw, time=start)

        if return_drivers:
            return drivers

    def _write_lineage_drivers(
        self: LinDriversProtocol,
        drivers: pd.DataFrame,
        use_raw: bool,
        *,
        time: datetime,
    ) -> None:
        self._lineage_drivers = drivers

        key = Key.varm.lineage_drivers(self.backward)
        if use_raw:
            field = "raw.varm"
            self.adata.raw.varm[key] = drivers
        else:
            field = "varm"
            self.adata.varm[key] = drivers

        logg.info(
            f"Adding `adata.{field}[{key!r}]`\n"
            f"       `.lineage_drivers`\n"
            "    Finish",
            time=time,
        )

    @d.get_sections(base="plot_lineage_drivers", sections=["Parameters"])
    @d.dedent
    def plot_lineage_drivers(
        self: LinDriversProtocol,
        lineage: str,
        n_genes: int = 8,
        ncols: Optional[int] = None,
        use_raw: bool = False,
        title_fmt: str = "{gene} qval={qval:.4e}",
        figsize: Optional[Tuple[float, float]] = None,
        dpi: Optional[int] = None,
        save: Optional[Union[str, Path]] = None,
        **kwargs: Any,
    ) -> None:
        """
        Plot lineage drivers discovered by :meth:`compute_lineage_drivers`.

        Parameters
        ----------
        lineage
            Lineage for which to plot the driver genes.
        n_genes
            Top most correlated genes to plot.
        ncols
            Number of columns.
        use_raw
            Whether to look in :attr:`adata` ``.raw.var`` or :attr:`adata` ``.var``.
        title_fmt
            Title format. Possible keywords include `{gene}`, `{qval}`, `{corr}` for gene name,
            q-value and correlation, respectively.
        %(plotting)s
        kwargs
            Keyword arguments for :func:`scvelo.pl.scatter`.

        Returns
        -------
        %(just_plots)s
        """

        def prepare_format(
            gene: str,
            *,
            pval: Optional[float],
            qval: Optional[float],
            corr: Optional[float],
        ) -> Dict[str, Any]:
            kwargs = {}
            if "{gene" in title_fmt:
                kwargs["gene"] = gene
            if "{pval" in title_fmt:
                kwargs["pval"] = float(pval) if pval is not None else np.nan
            if "{qval" in title_fmt:
                kwargs["qval"] = float(qval) if qval is not None else np.nan
            if "{corr" in title_fmt:
                kwargs["corr"] = float(corr) if corr is not None else np.nan

            return kwargs

        lin_drivers = self.lineage_drivers
        if lin_drivers is None:
            raise RuntimeError(
                "Compute `.lineage_drivers` first as `.compute_lineage_drivers()`."
            )

        key = f"{lineage} corr"
        if key not in lin_drivers:
            raise KeyError(
                f"Lineage `{key!r}` not found in `{list(lin_drivers.columns)}`."
            )

        if n_genes <= 0:
            raise ValueError(f"Expected `n_genes` to be positive, found `{n_genes}`.")

        kwargs.pop("save", None)
        genes = lin_drivers.sort_values(by=key, ascending=False).head(n_genes)

        ncols = 4 if ncols is None else ncols
        nrows = int(np.ceil(len(genes) / ncols))

        fig, axes = plt.subplots(
            ncols=ncols,
            nrows=nrows,
            dpi=dpi,
            figsize=(ncols * 6, nrows * 4) if figsize is None else figsize,
        )
        axes = np.ravel([axes])

        _i = 0
        for _i, (gene, ax) in enumerate(zip(genes.index, axes)):
            data = genes.loc[gene]
            scv.pl.scatter(
                self.adata,
                color=gene,
                ncols=ncols,
                use_raw=use_raw,
                ax=ax,
                show=False,
                title=title_fmt.format(
                    **prepare_format(
                        gene,
                        pval=data.get(f"{lineage} pval", None),
                        qval=data.get(f"{lineage} qval", None),
                        corr=data.get(f"{lineage} corr", None),
                    )
                ),
                **kwargs,
            )

        for j in range(_i + 1, len(axes)):
            axes[j].remove()

        if save is not None:
            save_fig(fig, save)

    @d.dedent
    def plot_lineage_drivers_correlation(
        self: LinDriversProtocol,
        lineage_x: str,
        lineage_y: str,
        color: Optional[str] = None,
        gene_sets: Optional[Dict[str, Sequence[str]]] = None,
        gene_sets_colors: Optional[Sequence[str]] = None,
        use_raw: bool = False,
        cmap: str = "RdYlBu_r",
        fontsize: int = 12,
        adjust_text: bool = False,
        legend_loc: Optional[str] = "best",
        figsize: Optional[Tuple[float, float]] = (4, 4),
        dpi: Optional[int] = None,
        save: Optional[Union[str, Path]] = None,
        show: bool = True,
        **kwargs: Any,
    ) -> Optional[Axes]:
        """
        Show scatter plot of gene-correlations between two lineages.

        Optionally, you can pass a :class:`dict` of gene names that will be annotated in the plot.

        Parameters
        ----------
        lineage_x
            Name of the lineage on the x-axis.
        lineage_y
            Name of the lineage on the y-axis.
        color
            Key in :attr:`adata` ``.var``.
        gene_sets
            Gene sets annotations of the form `{'gene_set_name': ['gene_1', 'gene_2'], ...}`.
        gene_sets_colors
            List of colors where each entry corresponds to a gene set from ``genes_sets``.
            If `None` and keys in ``gene_sets`` correspond to lineage names, use the lineage colors.
            Otherwise, use default colors.
        use_raw
            Whether to access :attr:`adata` ``.raw.var`` or :attr:`adata` ``.var``.
        cmap
            Colormap to use.
        fontsize
            Size of the text when plotting ``gene_sets``.
        adjust_text
            Whether to automatically adjust text in order to reduce overlap.
        legend_loc
            Position of the legend. If `None`, don't show the legend.
            Only used when ``gene_sets!=None``.
        %(plotting)s
        show
            If `False`, return :class:`matplotlib.pyplot.Axes`.
        kwargs
            Keyword arguments for :func:`scanpy.pl.scatter`.

        Returns
        -------
        :class:`matplotlib.pyplot.Axes`
            The axis object if ``show=False``.
        %(just_plots)s

        Notes
        -----
        This plot is based on the following
        `notebook <https://github.com/theislab/gastrulation_analysis/blob/main/6_cellrank.ipynb>`_ by Maren Büttner.
        """
        from cellrank.pl._utils import _position_legend

        if use_raw and self.adata.raw is None:
            logg.warning("No raw attribute set. Setting `use_raw=False`")
            use_raw = False
        adata = self.adata.raw if use_raw else self.adata

        # silent assumption: `.compute_lineage_drivers()` always writes to AnnData
        key1, key2 = f"to {lineage_x} corr", f"to {lineage_y} corr"
        if key1 not in adata.var or key2 not in adata.var:
            haystack = "adata.raw.var" if use_raw else "adata.var"
            raise RuntimeError(
                f"Unable to find correlations in `{haystack}[{key1!r}]` or `{haystack}[{key2!r}]`."
                f"Compute `.lineage_drivers` first as "
                f"`.compute_lineage_drivers([{lineage_x!r}, {lineage_y!r}], use_raw={use_raw})`."
            )
        # produce the actual scatter plot
        ctx = {"figure.figsize": figsize, "figure.dpi": dpi}
        for key in list(ctx.keys()):
            if ctx[key] is None:
                del ctx[key]
        with rc_context(ctx):
            ax = sc.pl.scatter(
                # TODO(michalk8): why?
                adata.to_adata() if adata.is_view else adata,
                x=key1,
                y=key2,
                color_map=cmap,
                use_raw=False,
                show=False,
                color=color,
                **kwargs,
            )
        fig = ax.figure

        # add some lines to highlight the origin
        xmin, xmax = np.nanmin(adata.var[key1]), np.nanmax(adata.var[key1])
        ymin, ymax = np.nanmin(adata.var[key2]), np.nanmax(adata.var[key2])
        ax.hlines(0, xmin=xmin, xmax=xmax, color="grey", alpha=0.5, zorder=-1)
        ax.vlines(0, ymin=ymin, ymax=ymax, color="grey", alpha=0.5, zorder=-1)

        # annotate the passed set of genes
        if gene_sets is not None:
            if gene_sets_colors is None:
                try:
                    # fmt: off
                    sets = list(gene_sets.keys())
                    gene_sets_colors = self.absorption_probabilities[sets].colors
                    # fmt: on
                except KeyError:
                    logg.warning(
                        "Unable to determine gene sets colors from lineages. Using default colors"
                    )
                    gene_sets_colors = _create_categorical_colors(len(gene_sets))
            if len(gene_sets_colors) != len(gene_sets):
                raise ValueError(
                    f"Expected `gene_sets_colors` to be of length `{len(gene_sets)}`, "
                    f"found `{len(gene_sets_colors)}`."
                )

            path_effect = [
                patheffects.Stroke(linewidth=2, foreground="w", alpha=0.8),
                patheffects.Normal(),
            ]
            annots = []
            for (key, values), color in zip(gene_sets.items(), gene_sets_colors):
                arrowprops = (
                    {
                        "arrowstyle": ArrowStyle.CurveFilledB(
                            head_width=0.1, head_length=0.2
                        ),
                        "fc": color,
                        "ec": color,
                    }
                    if adjust_text
                    else None
                )
                values = set(values) & set(adata.var_names)
                for value in values:
                    x = adata.var.loc[value, key1]
                    y = adata.var.loc[value, key2]

                    annot = ax.annotate(
                        value,
                        xy=(x, y),
                        va="top",
                        ha="left",
                        path_effects=path_effect,
                        arrowprops=arrowprops,
                        size=fontsize,
                        c=color,
                    )
                    annots.append(annot)
                if values:
                    ax.scatter([], [], color=color, label=key)

            if adjust_text:
                try:
                    import adjustText

                    start = logg.info("Adjusting text position")
                    adjustText.adjust_text(
                        annots,
                        x=adata.var[key1].values,
                        y=adata.var[key2].values,
                        ax=ax,
                    )
                    logg.info("    Finish", time=start)
                except ImportError:
                    logg.error(
                        "Please install `adjustText` first as `pip install adjustText`"
                    )
            if len(annots) and legend_loc not in (None, "none"):
                _position_legend(ax, legend_loc=legend_loc)

        if save is not None:
            save_fig(fig, path=save)

        if not show:
            return ax

    def _write_terminal_states(
        self: LinDriversProtocol,
        states: Optional[pd.Series],
        colors: Optional[np.ndarray],
        probs: Optional[pd.Series] = None,
        *,
        time: Optional[datetime],
        log: bool = True,
    ) -> None:
        super()._write_terminal_states(states, colors, probs, time=time, log=log)

        key = Key.varm.lineage_drivers(self.backward)
        self._set("_lineage_drivers", self.adata.varm, key=key, value=None)

    def _write_absorption_probabilities(
        self: LinDriversProtocol,
        abs_probs: Lineage,
        abs_times: Optional[pd.DataFrame],
        *,
        time: datetime,
        log: bool = True,
    ) -> None:
        super()._write_absorption_probabilities(
            abs_probs, abs_times, time=time, log=log
        )

        key = Key.varm.lineage_drivers(self.backward)
        self._set("_lineage_drivers", self.adata.varm, key=key, value=None)

    @property
    def lineage_drivers(self) -> Optional[pd.DataFrame]:
        """TODO."""
        return self._lineage_drivers