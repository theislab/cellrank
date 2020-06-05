# -*- coding: utf-8 -*-
"""Estimator module."""

from abc import ABC, abstractmethod
from typing import Any, Dict, List, Tuple, Union, Iterable, Optional, Sequence
from pathlib import Path

import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt

import scvelo as scv
from scanpy import logging as logg
from anndata import AnnData

import numpy as np
import pandas as pd
from pandas import Series
from scipy.stats import ranksums
from scipy.sparse import issparse
from pandas.api.types import infer_dtype, is_categorical_dtype
from scipy.sparse.linalg import eigs
from cellrank.tools._utils import (
    save_fig,
    _eigengap,
    _make_cat,
    partition,
    _vec_mat_corr,
    _complex_warning,
    _merge_categorical_series,
    _convert_to_categorical_series,
)
from cellrank.tools._colors import (
    _map_names_and_colors,
    _convert_to_hex_colors,
    _create_categorical_colors,
)
from cellrank.tools._lineage import Lineage
from cellrank.tools._constants import LinKey, Prefix, StateKey, Direction, _colors
from cellrank.tools.kernels._kernel import KernelExpression


class BaseEstimator(ABC):
    """Base class for all lineage probabilities estimators."""

    def __init__(
        self,
        kernel: KernelExpression,
        adata: Optional[AnnData] = None,
        inplace: bool = True,
        read_from_adata: bool = True,
        g2m_key: Optional[str] = "G2M_score",
        s_key: Optional[str] = "S_score",
        key_added: Optional[str] = None,
    ):
        if adata is not None:
            self._adata = adata if inplace else adata.copy()
        else:
            logg.debug("DEBUG: Loading `adata` object from kernel.")
            self._adata = kernel.adata if inplace else kernel.adata.copy()

        self._is_irreducible = None
        self._rec_classes = None
        self._trans_classes = None

        if kernel.backward:
            self._direction: Direction = Direction.BACKWARD
            self._rc_key: str = str(StateKey.BACKWARD)
            self._lin_key: str = str(LinKey.BACKWARD)
            self._prefix: str = str(Prefix.BACKWARD)
        else:
            self._direction: Direction = Direction.FORWARD
            self._rc_key: str = str(StateKey.FORWARD)
            self._lin_key: str = str(LinKey.FORWARD)
            self._prefix: str = str(Prefix.FORWARD)

        # import transition matrix and parameters
        if kernel.transition_matrix is None:
            logg.debug("DEBUG: Computing transition matrix using default parameters.")
            kernel.compute_transition_matrix()
        kernel.write_to_adata(key_added=key_added)

        self._kernel = kernel
        self._T = kernel.transition_matrix
        self._is_sparse = issparse(self._T)

        self._n_states = self._T.shape[0]
        if self._n_states != self._adata.n_obs:
            raise ValueError(
                f"Expected `{self._n_states}` (based on transition matrix), "
                f"found `{self._adata.n_obs}` (based on `adata` object)."
            )

        self._lin_probs = None

        # for copy
        self._g2m_key = g2m_key
        self._s_key = s_key
        self._key_added = key_added

        self._eig = None  # stores eigendecomposition
        self._dp = None  # stores differentiation potential
        self._G2M_score = None  # cell-cycle
        self._S_score = None

        if read_from_adata:
            logg.debug("DEBUG: Reading data from `adata` object")
            self._read_from_adata(g2m_key, s_key)

    def _compute_eig(
        self,
        k: int = 20,
        which: str = "LR",
        alpha: float = 1,
        only_evals: bool = False,
        ncv: Optional[int] = None,
    ) -> None:
        """
        Compute eigendecomposition of transition matrix.

        Uses a sparse implementation, if possible, and only computes the top k eigenvectors
        to speed up the computation. Computes both left and right eigenvectors.

        Params
        ------
        k
            Number of eigenvalues/vectors to compute.
        which
            Eigenvalues are in general complex. `'LR'` - largest real part, `'LM'` - largest magnitude.
        alpha
            Used to compute the `eigengap`. :paramref:`alpha` is the weight given
            to the deviation of an eigenvalue from one.
        only_evals
            Return only eigenvalues
        ncv
            Number of Lanczos vectors generated

        Returns
        -------
        None
            Nothing, but updates the following fields: :paramref:`eigendecomposition`.
        """

        def get_top_k_evals():
            return D[np.flip(np.argsort(D.real))][:k]

        logg.info("Computing eigendecomposition of transition matrix")
        if self._is_sparse:
            logg.debug(f"DEBUG: Computing top `{k}` eigenvalues for sparse matrix")
            D, V_l = eigs(self._T.T, k=k, which=which, ncv=ncv)
            if only_evals:
                self._write_eig_to_adata(
                    {
                        "D": get_top_k_evals(),
                        "eigengap": _eigengap(get_top_k_evals().real, alpha),
                        "params": {"which": which, "k": k, "alpha": alpha},
                    }
                )
                return
            _, V_r = eigs(self._T, k=k, which=which, ncv=ncv)
        else:
            logg.warning(
                "This transition matrix is not sparse, computing full eigendecomposition"
            )
            D, V_l = np.linalg.eig(self._T.T)
            if only_evals:
                self._write_eig_to_adata(
                    {
                        "D": get_top_k_evals(),
                        "eigengap": _eigengap(D.real, alpha),
                        "params": {"which": which, "k": k, "alpha": alpha},
                    }
                )
                return
            _, V_r = np.linalg.eig(self._T)

        # Sort the eigenvalues and eigenvectors and take the real part
        logg.debug("DEBUG: Sorting eigenvalues by their real part")
        p = np.flip(np.argsort(D.real))
        D, V_l, V_r = D[p], V_l[:, p], V_r[:, p]
        e_gap = _eigengap(D.real, alpha)

        pi = np.abs(V_l[:, 0].real)
        pi /= np.sum(pi)

        self._write_eig_to_adata(
            {
                "D": D,
                "stationary_dist": pi,
                "V_l": V_l,
                "V_r": V_r,
                "eigengap": e_gap,
                "params": {"which": which, "k": k, "alpha": alpha},
            }
        )

    @abstractmethod
    def compute_eig(self, k: int = 20, which: str = "LR", alpha: float = 1) -> None:
        """
        Compute eigendecomposition of the transition matrix.

        Uses a sparse implementation, if possible, and only computes the top k eigenvalues.

        Params
        ------
        k
            Number of eigenvalues/vectors to compute.
        which
            Eigenvalues are in general complex. `'LR'` - largest real part, `'LM'` - largest magnitude.
        alpha
            Used to compute the `eigengap`. paramref:`alpha` is the weight given
            to the deviation of an eigenvalue from one.
        ncv
            Number of Lanczos vectors generated.

        Returns
        -------
        None
            Nothing, but updates the following fields: paramref:`eigendecomposition`.
        """
        pass

    def compute_partition(self) -> None:
        """
        Compute communication classes for the Markov chain.

        Returns
        -------
        None
            Nothing, but updates the following fields:
                - :paramref:`recurrent_classes`
                - :paramref:`transient_classes`
                - :paramref:`irreducible`
        """

        start = logg.info("Computing communication classes")

        rec_classes, trans_classes = partition(self._T)

        self._is_irreducible = len(rec_classes) == 1 and len(trans_classes) == 0

        if not self._is_irreducible:
            self._trans_classes = _make_cat(
                trans_classes, self._n_states, self._adata.obs_names
            )
            self._rec_classes = _make_cat(
                rec_classes, self._n_states, self._adata.obs_names
            )
            self._adata.obs[f"{self._rc_key}_rec_classes"] = self._rec_classes
            self._adata.obs[f"{self._rc_key}_trans_classes"] = self._trans_classes
            logg.info(
                f"Found `{(len(rec_classes))}` recurrent and `{len(trans_classes)}` transient classes\n"
                f"Adding `.recurrent_classes`\n"
                f"       `.transient_classes`\n"
                f"       `.irreducible`\n"
                f"    Finish",
                time=start,
            )
        else:
            logg.warning(
                "The transition matrix is irreducible - cannot further partition it\n    Finish",
                time=start,
            )

    def plot_spectrum(
        self,
        real_only: bool = False,
        dpi: int = 100,
        figsize: Optional[Tuple[float, float]] = (5, 5),
        legend_loc: Optional[str] = None,
        title: Optional[str] = None,
        save: Optional[Union[str, Path]] = None,
    ) -> None:
        """
        Plot the top eigenvalues in complex plane.

        Params
        ------
        real_only
            Whether to plot only the real part of the spectrum.
        dpi
            Dots per inch.
        figsize
            Size of the figure.
        legend_loc
            Location parameter for the legend
        title
            Figure title
        save
            Filename where to save the plots. If `None`, just shows the plot.

        Returns
        -------
        None
            Nothing, just plots the spectrum in complex plane.
        """

        # define a function to make the data limits rectangular
        def adapt_range(min_, max_, range_):
            return (
                min_ + (max_ - min_) / 2 - range_ / 2,
                min_ + (max_ - min_) / 2 + range_ / 2,
            )

        if self._eig is None:
            logg.warning(
                "No eigendecomposition found, computing with default parameters"
            )
            self.compute_eig()

        if real_only:
            self._plot_real_spectrum(
                dpi=dpi, figsize=figsize, legend_loc=legend_loc, save=save, title=title
            )
            return

        D = self._eig["D"]
        params = self._eig["params"]

        # create fiture and axes
        fig, ax = plt.subplots(nrows=1, ncols=1, dpi=dpi, figsize=figsize)

        # get the original data ranges
        lam_x, lam_y = D.real, D.imag
        x_min, x_max = np.min(lam_x), np.max(lam_x)
        y_min, y_max = np.min(lam_y), np.max(lam_y)
        x_range, y_range = x_max - x_min, y_max - y_min
        final_range = np.max([x_range, y_range]) + 0.05

        x_min_, x_max_ = adapt_range(x_min, x_max, final_range)
        y_min_, y_max_ = adapt_range(y_min, y_max, final_range)

        # plot the data and the unit circle
        ax.scatter(D.real, D.imag, marker=".", label="eigenvalue")
        t = np.linspace(0, 2 * np.pi, 500)
        x_circle, y_circle = np.sin(t), np.cos(t)
        ax.plot(x_circle, y_circle, "k-", label="unit circle")

        # set labels, ranges and legend
        ax.set_xlabel(r"Re($\lambda$)")
        ax.set_xlim(x_min_, x_max_)
        ax.set_ylabel(r"Im($\lambda$)")
        ax.set_ylim(y_min_, y_max_)
        key = "real part" if params["which"] == "LR" else "magnitude"

        # set the figure title
        if title is None:
            fig_title = f"top {params['k']} eigenvalues according to their {key}"
        else:
            fig_title = title
        ax.set_title(fig_title)

        # set legend location
        ax.legend(loc=legend_loc)

        if save is not None:
            save_fig(fig, save)

        fig.show()

    def _plot_real_spectrum(
        self,
        dpi: int = 100,
        figsize: Optional[Tuple[float, float]] = None,
        legend_loc: Optional[str] = None,
        title: Optional[str] = None,
        save: Optional[Union[str, Path]] = None,
    ) -> None:
        """
        Plot the real part of the top eigenvalues.

        Params
        ------
        dpi
            Dots per inch.
        figsize
            Size of the figure.
        legend_loc
            Location parameter for the legend.
        title
            Figure title.
        save
            Filename where to save the plots. If `None`, just shows the plot.

        Returns
        -------
        None
            Nothing, just plots the spectrum.
        """

        if self._eig is None:
            logg.warning(
                "No eigendecomposition found, computing with default parameters"
            )
            self.compute_eig()

        # Obtain the eigendecomposition, create the color code
        D, params = self._eig["D"], self._eig["params"]
        D_real, D_imag = D.real, D.imag
        ixs = np.arange(len(D))
        mask = D_imag == 0

        # plot the top eigenvalues
        fig, ax = plt.subplots(nrows=1, ncols=1, dpi=dpi, figsize=figsize)
        ax.scatter(ixs[mask], D_real[mask], marker="o", label="real eigenvalue")
        ax.scatter(ixs[~mask], D_real[~mask], marker="o", label="complex eigenvalue")

        # add dashed line for the eigengap, ticks, labels, title and legend
        ax.axvline(self._eig["eigengap"], label="eigengap", ls="--")

        ax.set_xlabel("index")
        ax.set_xticks(range(len(D)))

        ax.set_ylabel(r"Re($\lambda_i$)")
        key = "real part" if params["which"] == "LR" else "magnitude"

        # set the title
        if title is None:
            fig_title = (
                f"real part of top {params['k']} eigenvalues according to their {key}"
            )
        else:
            fig_title = title

        ax.set_title(fig_title)

        ax.legend(loc=legend_loc)

        if save is not None:
            save_fig(fig, save)

        fig.show()

    def _set_categorical_labels(
        self,
        attr_key: str,
        pretty_attr_key: str,
        cat_key: str,
        add_to_existing_error_msg: str,
        categories: Union[Series, Dict[Any, Any]],
        cluster_key: Optional[str] = None,
        en_cutoff: Optional[float] = None,
        p_thresh: Optional[float] = None,
        add_to_existing: bool = False,
    ):
        if isinstance(categories, dict):
            categories = _convert_to_categorical_series(
                categories, list(self.adata.obs_names)
            )
        if not is_categorical_dtype(categories):
            raise TypeError(
                f"Object must be `categorical`, found `{infer_dtype(categories)}`."
            )

        if add_to_existing:
            if getattr(self, attr_key) is None:
                raise RuntimeError(add_to_existing_error_msg)
            categories = _merge_categorical_series(
                getattr(self, attr_key), categories, inplace=False
            )

        if cluster_key is not None:
            logg.debug(f"DEBUG: Creating colors based on `{cluster_key}`")

            # check that we can load the reference series from adata
            if cluster_key not in self._adata.obs:
                raise KeyError(
                    f"Cluster key `{cluster_key!r}` not found in `adata.obs`."
                )
            series_query, series_reference = categories, self._adata.obs[cluster_key]

            # load the reference colors if they exist
            if _colors(cluster_key) in self._adata.uns.keys():
                colors_reference = _convert_to_hex_colors(
                    self._adata.uns[_colors(cluster_key)]
                )
            else:
                colors_reference = _create_categorical_colors(
                    len(series_reference.cat.categories)
                )

            approx_rcs_names, colors = _map_names_and_colors(
                series_reference=series_reference,
                series_query=series_query,
                colors_reference=colors_reference,
                en_cutoff=en_cutoff,
            )
            setattr(self, f"{attr_key}_colors", colors)
            categories.cat.categories = approx_rcs_names
        else:
            setattr(
                self,
                f"{attr_key}_colors",
                _create_categorical_colors(len(categories.cat.categories)),
            )

        if p_thresh is not None:
            self._detect_cc_stages(categories, p_thresh=p_thresh)

        # write to class and adata
        if getattr(self, attr_key) is not None:
            logg.debug(f"DEBUG: Overwriting `.{pretty_attr_key}`")

        setattr(self, attr_key, categories)
        self._adata.obs[cat_key] = categories.values
        self._adata.uns[_colors(cat_key)] = getattr(self, f"{attr_key}_colors")

    def _plot_vectors(
        self,
        vectors: np.ndarray,
        kind: str,
        abs_value: bool = False,
        cluster_key: Optional[str] = None,
        use: Optional[Union[int, Tuple[int], List[int]]] = None,
        **kwargs,
    ):
        if kind not in ("eigen", "schur"):
            raise ValueError(
                f"Invalid kind `{kind!r}`. Valid options are `'eigen', 'schur'`."
            )
        is_schur = kind == "schur"

        # check whether dimensions are consistent
        if self.adata.n_obs != vectors.shape[0]:
            raise ValueError(
                f"Number of cells ({self.adata.n_obs}) is inconsistent with the 1."
                f"dimensions of vectors ({vectors.shape[0]})."
            )

        if use is None:
            use = list(range(is_schur, vectors.shape[1] + is_schur - 1))
        elif isinstance(use, int):
            use = list(range(is_schur, use + is_schur))
        elif not isinstance(use, (tuple, list, range)):
            raise TypeError(
                f"Argument `use` must be either `int`, `tuple`, `list` or `range`,"
                f"found `{type(use).__name__}`."
            )
        else:
            if not all(map(lambda u: isinstance(u, int), use)):
                raise TypeError("Not all values in `use` argument are integers.")
        use = list(use)
        if not use:
            raise ValueError("No vectors to plot.")

        muse = max(use)
        if muse >= vectors.shape[1]:
            vec = "Schur " if is_schur else "eigen"
            raise ValueError(
                f"Maximum specified {vec}vector ({muse}) is larger "
                f"than the number of computed {vec}vectors ({vectors.shape[1]})."
            )
        V_ = vectors[:, use]
        if is_schur:
            title = [f"Schur vector {i}" for i in use]
        else:
            D = kwargs.pop("D")
            V_ = _complex_warning(V_, use, use_imag=kwargs.pop("use_imag", False))
            title = [fr"$\lambda_{i}$={d:.02f}" for i, d in zip(use, D[use])]

        if abs_value:
            V_ = np.abs(V_)

        color = list(V_.T)
        if cluster_key is not None:
            color = [cluster_key] + color

        # actual plotting with scvelo
        logg.debug(f"DEBUG: Showing `{use}` Schur vectors")
        scv.pl.scatter(self._adata, color=color, title=title, **kwargs)

    @abstractmethod
    def _read_from_adata(
        self, g2m_key: Optional[str] = None, s_key: Optional[str] = None, **kwargs
    ) -> None:
        pass

    def _detect_cc_stages(self, rc_labels: Series, p_thresh: float = 1e-15) -> None:
        """
        Detect cell-cycle driven start or endpoints.

        Params
        ------
        rc_labels
            Approximate recurrent classes.
        p_thresh
            P-value threshold for the rank-sum test for the group to be considered cell-cycle driven.
        Returns
        -------
        None
            Nothing, but warns if a group is cell-cycle driven.
        """

        # initialise the groups (start or end clusters) and scores
        groups = rc_labels.cat.categories
        scores = []
        if self._G2M_score is not None:
            scores.append(self._G2M_score)
        if self._S_score is not None:
            scores.append(self._S_score)

        # loop over groups and scores
        for group in groups:
            mask = rc_labels == group
            for score in scores:
                a, b = score[mask], score[~mask]
                result = ranksums(a, b)
                if result.statistic > 0 and result.pvalue < p_thresh:
                    logg.warning(f"Group `{group}` appears to be cell-cycle driven")
                    break

    def _plot_probabilities(
        self,
        attr: str,
        error_msg: str,
        lineages: Optional[Union[str, Iterable[str]]] = None,
        cluster_key: Optional[str] = None,
        mode: str = "embedding",
        time_key: str = "latent_time",
        show_dp: bool = True,
        title: Optional[str] = None,
        same_plot: bool = False,
        color_map: Union[str, mpl.colors.ListedColormap] = cm.viridis,
        **kwargs,
    ) -> None:
        probs: Optional[Lineage] = getattr(self, attr)

        if probs is None:
            raise RuntimeError(error_msg)
        if isinstance(lineages, str):
            lineages = [lineages]

        # retrieve the lineage data
        if lineages is None:
            lineages = probs.names
            A = probs
        else:
            for lineage in lineages:
                if lineage not in probs.names:
                    raise ValueError(
                        f"Invalid lineage name `{lineages!r}`. Valid options are `{list(probs.names)}`."
                    )
            A = probs[lineages]

        show_dp = show_dp and self._dp is not None
        # change the maximum value - the 1 is artificial and obscures the color scaling
        for col in A.T:
            mask = col != 1
            if np.sum(mask) > 0:
                max_not_one = np.max(col[mask])
                col[~mask] = max_not_one

        if mode == "time":
            if time_key not in self._adata.obs.keys():
                raise KeyError(f"Time key `{time_key!r}` not found in `adata.obs`.")
            t = self.adata.obs[time_key]
            cluster_key = None

        rc_titles = (
            [f"{self._prefix} {rc}" for rc in lineages]
            + (["Differentiation potential"] if show_dp else [])
            if title is None
            else title
        )

        if cluster_key is not None:
            color = [cluster_key] + list(A.T) + ([self._dp] if show_dp else [])
            titles = [cluster_key] + rc_titles
        else:
            color = list(A.T) + ([self._dp] if show_dp else [])
            titles = rc_titles

        if mode == "embedding":
            if same_plot:
                scv.pl.scatter(
                    self.adata,
                    title=(
                        "metastable states (backward)"
                        if self.kernel.backward
                        else "metastable states (forward)"
                    )
                    if title is None
                    else title,
                    color_gradients=A,
                    color_map=color_map,
                    **kwargs,
                )
            else:
                scv.pl.scatter(
                    self.adata,
                    color=color,
                    title=titles,
                    color_map=color_map,
                    **kwargs,
                )
        elif mode == "time":
            xlabel, ylabel = [time_key] * len(titles), ["probability"] * len(titles)
            scv.pl.scatter(
                self._adata,
                x=t,
                color_map=color_map,
                y=list(A.T) + [self._dp],
                title=titles,
                xlabel=xlabel,
                ylabel=ylabel,
                **kwargs,
            )
        else:
            raise ValueError(
                f"Invalid mode `{mode!r}`. Valid options are: `'embedding', 'time'`."
            )

    def compute_lineage_drivers(
        self,
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

        # check that lineage probs have been computed
        if self._lin_probs is None:
            raise RuntimeError(
                "Compute lineage probabilities first as `.compute_lin_probs()` or `.set_main_states`."
            )

        # check all lin_keys exist in self.lin_names
        if isinstance(lin_names, str):
            lin_names = [lin_names]
        if lin_names is not None:
            _ = self._lin_probs[lin_names]
        else:
            lin_names = self._lin_probs.names

        # use `cluster_key` and clusters to subset the data
        if clusters is not None:
            if cluster_key is None:
                cluster_key = "clusters"
            if cluster_key not in self._adata.obs.keys():
                raise KeyError(f"Key `{cluster_key!r}` not found in `adata.obs`.")
            if isinstance(clusters, str):
                clusters = [clusters]
            all_clusters = np.array(self._adata.obs[cluster_key].cat.categories)
            cluster_mask = np.array([name not in all_clusters for name in clusters])
            if any(cluster_mask):
                raise KeyError(
                    f"Clusters `{list(np.array(clusters)[cluster_mask])}` not found in "
                    f"`adata.obs[{cluster_key!r}]`."
                )
            subset_mask = np.in1d(self._adata.obs[cluster_key], clusters)
            adata_comp = self._adata[subset_mask].copy()
            lin_probs = self._lin_probs[subset_mask, :]
        else:
            adata_comp = self._adata.copy()
            lin_probs = self._lin_probs

        # check that the layer exists, and that use raw is only used with layer X
        if layer != "X":
            if layer not in self._adata.layers:
                raise KeyError(f"Layer `{layer!r}` not found in `adata.layers`.")
            if use_raw:
                raise ValueError("For `use_raw=True`, layer must be 'X'.")
            data = adata_comp.layers[layer]
            var_names = adata_comp.var_names
        else:
            if use_raw and self._adata.raw is None:
                logg.warning("No raw attribute set. Using `.X` instead")
                use_raw = False
            data = adata_comp.raw.X if use_raw else adata_comp.X
            var_names = adata_comp.raw.var_names if use_raw else adata_comp.var_names

        start = logg.info(
            f"Computing correlations for lineages `{lin_names}` restricted to clusters `{clusters}` in "
            f"layer `{layer}` with `use_raw={use_raw}`"
        )

        # loop over lineages
        lin_corrs = {}
        for lineage in lin_names:
            y = lin_probs[:, lineage].X.squeeze()
            correlations = _vec_mat_corr(data, y)

            if inplace:
                if use_raw:
                    self._adata.raw.var[f"{self._prefix} {lineage} corr"] = correlations
                else:
                    self._adata.var[f"{self._prefix} {lineage} corr"] = correlations
            else:
                lin_corrs[lineage] = correlations

        if not inplace:
            return pd.DataFrame(lin_corrs, index=var_names)

        field = "raw.var" if use_raw else "var"
        logg.info(
            f"Adding gene correlations to `.adata.{field}`\n    Finish", time=start
        )

    def _write_eig_to_adata(self, eig):
        # write to class and AnnData object
        if self._eig is not None:
            logg.debug("DEBUG: Overwriting `.eigendecomposition`")
        else:
            logg.debug(
                f"DEBUG: Adding `.eigendecomposition` and `adata.uns['eig_{self._direction}']`"
            )

        self._eig = eig
        self._adata.uns[f"eig_{self._direction}"] = eig

    @abstractmethod
    def copy(self) -> "BaseEstimator":
        """Return a copy of self."""
        pass

    @property
    def irreducible(self) -> Optional[bool]:
        """Whether the Markov chain is irreducible or not."""
        return self._is_irreducible

    @property
    def recurrent_classes(self) -> Optional[List[List[Any]]]:
        """The recurrent classes of the Markov chain."""  # noqa
        return self._rec_classes

    @property
    def transient_classes(self) -> Optional[List[List[Any]]]:
        """The recurrent classes of the Markov chain."""  # noqa
        return self._trans_classes

    @property
    def eigendecomposition(self) -> Optional[Dict[str, np.ndarray]]:
        """
        A dictionary with the following fields:

        - `'D'` eigenvalues of left eigenvectors
        - `'V_l'` left eigenvectors
        - `'V_r'` right eigenvectors
        """  # noqa
        return self._eig

    @property
    def diff_potential(self) -> np.ndarray:
        """Differentiation potential for each lineage."""  # noqa
        return self._dp

    @property
    def lineage_probabilities(self) -> Lineage:
        """An array with names and colors, where each column represents one lineage."""  # noqa
        return self._lin_probs

    @property
    def adata(self) -> AnnData:
        """The underlying annotated data object."""  # noqa
        return self._adata

    @property
    def kernel(self) -> KernelExpression:
        """The underlying kernel expression."""  # noqa
        return self._kernel

    def __copy__(self) -> "BaseEstimator":
        return self.copy()

    def __len__(self) -> int:
        return self._n_states

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}[n={len(self)}, kernel={repr(self._kernel)}]"

    def __str__(self) -> str:
        return f"{self.__class__.__name__}[n={len(self)}, kernel={str(self._kernel)}]"
