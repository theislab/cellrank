# -*- coding: utf-8 -*-
"""Models module."""

import re
from abc import ABC, abstractmethod
from copy import copy as _copy
from copy import deepcopy
from types import MappingProxyType
from typing import Any, Tuple, Union, TypeVar, Iterable, Optional
from inspect import signature
from collections import Mapping, defaultdict

import numpy as np
import pandas as pd
from pygam import GAM as pGAM
from pygam import (
    GammaGAM,
    LinearGAM,
    PoissonGAM,
    InvGaussGAM,
    LogisticGAM,
    ExpectileGAM,
)
from pygam.terms import s
from sklearn.base import BaseEstimator
from scipy.ndimage.filters import convolve

import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

import cellrank.logging as logg
from cellrank.utils._docs import d
from cellrank.tools._utils import save_fig, _densify_squeeze
from cellrank.utils._utils import _minmax
from cellrank.tools._lineage import Lineage
from cellrank.tools._constants import ModeEnum, AbsProbKey
from cellrank.tools.kernels._utils import _filter_kwargs

AnnData = TypeVar("AnnData")


_dup_spaces = re.compile(r" +")
_r_lib = None
_r_lib_name = None

_valid_distributions = ()

_valid_link_functions = ()


class GamDistribution(ModeEnum):  # noqa
    NORMAL = "normal"
    BINOMIAL = "binomial"
    POISSON = "poisson"
    GAMMA = "gamma"
    GAUSS = "gaussian"
    INV_GAUSS = "inv_gauss"


class GamLinkFunction(ModeEnum):  # noqa
    IDENTITY = "identity"
    LOGIT = "logit"
    INV = "inverse"
    LOG = "log"
    INV_SQUARED = "inverse-squared"


_gams = defaultdict(
    lambda: pGAM,
    {
        (GamDistribution.NORMAL, GamLinkFunction.IDENTITY): LinearGAM,
        (GamDistribution.BINOMIAL, GamLinkFunction.LOGIT): LogisticGAM,
        (GamDistribution.POISSON, GamLinkFunction.LOG): PoissonGAM,
        (GamDistribution.GAMMA, GamLinkFunction.LOG): GammaGAM,
        (GamDistribution.INV_GAUSS, GamLinkFunction.LOG): InvGaussGAM,
    },
)


@d.get_sectionsf("base_model", sections=["Parameters"])
@d.dedent
class BaseModel(ABC):
    """
    Base class for other model classes.

    Parameters
    ----------
    %(adata)s
    model
        Underlying model.
    """

    def __init__(
        self, adata: AnnData, model: Any,
    ):
        self._adata = adata
        self._model = model
        self._gene = None
        self._lineage = None

        self._x_all = None
        self._y_all = None
        self._w_all = None

        self._x = None
        self._y = None
        self._w = None

        self._x_test = None
        self._y_test = None
        self._x_hat = None
        self._y_hat = None

        self._conf_int = None

        self._dtype = np.float32

    @property
    @d.dedent
    def adata(self) -> AnnData:
        """
        Returns
        -------
        %(adata)s
        """  # noqa
        return self._adata

    @property
    def model(self) -> Any:
        """The underlying model."""  # noqa
        return self._model

    @property
    @d.get_summaryf("base_model_x_all")
    def x_all(self) -> np.ndarray:
        """Unfiltered independent variables of shape `(n_cells, 1)`."""
        return self._x_all

    @property
    @d.get_summaryf("base_model_y_all")
    def y_all(self) -> np.ndarray:
        """Unfiltered dependent variables of shape `(n_cells, 1)`."""
        return self._y_all

    @property
    @d.get_summaryf("base_model_w_all")
    def w_all(self) -> np.ndarray:
        """Unfiltered weights of shape `(n_cells,)`."""
        return self._w_all

    @property
    @d.get_summaryf("base_model_x")
    def x(self) -> np.ndarray:
        """Filtered independent variables of shape `(n_filtered_cells, 1)` used for fitting."""  # noqa
        return self._x

    @property
    @d.get_summaryf("base_model_y")
    def y(self) -> np.ndarray:
        """Filtered dependent variables of shape `(n_filtered_cells, 1)` used for fitting."""  # noqa
        return self._y

    @property
    @d.get_summaryf("base_model_w")
    def w(self) -> np.ndarray:
        """Filtered weights of shape `(n_filtered_cells,)` used for fitting."""  # noqa
        return self._w

    @property
    @d.get_summaryf("base_model_x_test")
    def x_test(self) -> np.ndarray:
        """Independent variables of shape `(n_samples, 1)` used for prediction."""
        return self._x_test

    @property
    @d.get_summaryf("base_model_y_test")
    def y_test(self) -> np.ndarray:
        """Prediction values of shape `(n_samples,)` for :paramref:`x_test`."""
        return self._y_test

    @property
    @d.get_summaryf("base_model_x_hat")
    def x_hat(self) -> np.ndarray:
        """Filtered independent variables used when calculating default confidence interval, usually same as :paramref:`x`."""  # noqa
        return self._x_hat

    @property
    @d.get_summaryf("base_model_y_hat")
    def y_hat(self) -> np.ndarray:
        """Filtered dependent variables used when calculating default confidence interval, usually same as :paramref:`y`."""  # noqa
        return self._y_hat

    @property
    @d.get_summaryf("base_model_conf_int")
    def conf_int(self) -> np.ndarray:
        """Array of shape `(n_samples, 2)` containing the lower and upper bounds of the confidence interval."""  # noqa
        return self._conf_int

    @d.dedent
    def prepare(
        self,
        gene: str,
        lineage: str,
        backward: bool = False,
        time_range: Optional[Union[float, Tuple[float, float]]] = None,
        data_key: str = "X",
        time_key: str = "latent_time",
        use_raw: bool = False,
        threshold: Optional[float] = None,
        weight_threshold: Union[float, Tuple[float, float]] = (0.01, 0.01),
        filter_dropouts: Optional[float] = None,
        n_test_points: int = 200,
    ) -> "BaseModel":
        """
        Prepare the model to be ready for fitting.

        Parameters
        ----------
        gene
            Gene in :paramref:`adata` `.var_names` or in :paramref:`adata` `.raw.var_names`.
        lineage
            Name of a lineage in :paramref:`adata` `.uns`:paramref:`lineage_key`.
        %(backward)s
        %s(time_range)s
        data_key
            Key in :attr:`paramref.adata` `.layers` or `'X'` for :paramref:`adata` `.X`.
        time_key
            Key in :paramref:`adata` `.obs` where the pseudotime is stored.
        use_raw
            Whether to access :paramref:`adata` `.raw` or not.
        threshold
            Consider only cells with :paramref:`weights` > :paramref:`threshold` when estimating the testing endpoint.
            If `None`, use median of :paramref:`w`.
        weight_threshold
            Set all weights below :paramref:`weight_threshold` to either 0, or if :class:`tuple`, to the second value.
        filter_dropouts
            Filter out all cells with expression lower than this.
        n_test_points
            Number of testing points. If `None`, use the original points based on :paramref:`threshold`.

        Returns
        -------
        None
            Nothing, but updates the following fields:

                - :paramref:`x` - %(base_model_x.summary)s
                - :paramref:`y` - %(base_model_y.summary)s
                - :paramref:`w` - %(base_model_w.summary)s

                - :paramref:`x_all` - %(base_model_x_all.summary)s
                - :paramref:`y_all` - %(base_model_y_all.summary)s
                - :paramref:`w_all` - %(base_model_w_all.summary)s

                - :paramref:`x_test` - %(base_model_x_test.summary)s
        """
        if use_raw and self.adata.raw is None:
            raise AttributeError("AnnData object has no attribute `.raw`.")

        if data_key not in ["X", "obs"] + list(self.adata.layers.keys()):
            raise KeyError(
                f"Data key must be a key of `adata.layers`: `{list(self.adata.layers.keys())}`, '`obs`' or `'X'`."
            )
        if time_key not in self.adata.obs:
            raise KeyError(f"Time key `{time_key!r}` not found in `adata.obs`.")

        if data_key != "obs":
            if use_raw and gene not in self.adata.raw.var_names:
                raise KeyError(f"Gene `{gene!r}` not found in `adata.raw.var_names`.")
            if not use_raw and gene not in self.adata.var_names:
                raise KeyError(f"Gene `{gene!r}` not found in `adata.var_names`.")
        else:
            if gene not in self.adata.obs:
                raise KeyError(f"Unable to find key `{gene!r}` in `adata.obs`.")

        lineage_key = str(AbsProbKey.BACKWARD if backward else AbsProbKey.FORWARD)
        if lineage_key not in self.adata.obsm:
            raise KeyError(f"Lineage key `{lineage_key!r}` not found in `adata.obsm`.")
        if not isinstance(self.adata.obsm[lineage_key], Lineage):
            raise TypeError(
                f"Expected `adata.obsm[{lineage_key!r}]` to be of type `cellrank.tl.Lineage`, "
                f"found `{type(self.adata.obsm[lineage_key]).__name__!r}`."
            )

        if not isinstance(time_range, (type(None), float, int, tuple)):
            raise TypeError(
                f"Expected time range to be either `None`, "
                f"`float` or `tuple`, found `{type(time_range).__name__!r}`."
            )
        if isinstance(time_range, tuple):
            if len(time_range) != 2:
                raise ValueError(
                    f"Expected time range to be a `tuple` of length `2`, found `{len(time_range)}`."
                )
            if not isinstance(
                time_range[0], (float, int, type(None))
            ) or not isinstance(time_range[1], (float, int, type(None))):
                raise TypeError(
                    f"Expected values in time ranges be of types `float` or `int` or `None`, "
                    f"got `{type(time_range[0]).__name__!r}` and `{type(time_range[1]).__name__!r}`."
                )
            val_start, val_end = time_range
        elif time_range is None:
            val_start, val_end = None, None
        else:
            val_start, val_end = None, time_range

        if lineage is not None:
            _ = self.adata.obsm[lineage_key][lineage]

        x = np.array(self.adata.obs[time_key]).astype(np.float64)

        adata = self.adata.raw.to_adata() if use_raw else self.adata
        gene_ix = np.where(adata.var_names == gene)[0]

        if data_key == "X":
            y = adata.X[:, gene_ix]
        elif data_key == "obs":
            y = adata.obs[gene].values
        elif data_key in adata.layers:
            y = adata.layers[data_key][:, gene_ix]
        else:
            raise NotImplementedError(f"Data key `{data_key!r}` is not implemented.")

        if lineage is not None:
            weight_threshold, val = (
                weight_threshold
                if isinstance(weight_threshold, (tuple, list))
                else (weight_threshold, 0)
            )
            w = _densify_squeeze(self.adata.obsm[lineage_key][lineage].X, self._dtype)
            w[w < weight_threshold] = val
        else:
            w = np.ones(len(x), dtype=self._dtype)

        if use_raw:
            correct_ixs = np.isin(self.adata.obs_names, adata.obs_names)
            x = x[correct_ixs]
            w = w[correct_ixs]

        del adata

        self._x_all, self._y_all, self._w_all = (
            _densify_squeeze(x, self._dtype)[:, np.newaxis],
            _densify_squeeze(y, self._dtype)[:, np.newaxis],
            w,
        )
        x, y, w = self.x_all[:], self.y_all[:], self.w_all[:]

        # sanity checks
        if self._x_all.shape[0] != self._y_all.shape[0]:
            raise ValueError(
                f"Independent variable's first dimension ({self._x_all.shape[0]}) "
                f"differs from dependent variable's first dimension ({self._y_all.shape[0]})."
            )
        if self._x_all.shape[0] != self._w_all.shape[0]:
            raise ValueError(
                f"Independent variable's first dimension ({self._x_all.shape[0]}) "
                f"differs from weights' first dimension ({self._w_all.shape[0]})."
            )

        x, ixs = np.unique(x, return_index=True)  # GAMR (mgcv) needs unique
        y = y[ixs]
        w = w[ixs]

        ixs = np.argsort(x)
        x, y, w = x[ixs], y[ixs], w[ixs]

        if val_start is None:
            val_start = np.nanmin(self.adata.obs[time_key])
        if val_end is None:
            if threshold is None:
                threshold = np.nanmedian(w)
            w_test = w[w > threshold]
            n_window = 10 if n_test_points is None else n_test_points // 20
            tmp = convolve(w_test, np.ones(n_window) / n_window, mode="nearest")
            val_end = x[w > threshold][np.nanargmax(tmp)]

        if val_start > val_end:
            val_start, val_end = val_end, val_start
        val_start, val_end = (
            max(val_start, np.min(self.adata.obs[time_key])),
            min(val_end, np.max(self.adata.obs[time_key])),
        )

        fil = (x >= val_start) & (x <= val_end)
        x_test = (
            np.linspace(val_start, val_end, n_test_points)
            if n_test_points is not None
            else x[fil]
        )
        x, y, w = x[fil], y[fil], w[fil]

        if filter_dropouts is not None:
            tmp = y.squeeze()
            fil = (tmp >= filter_dropouts) & (
                ~np.isclose(tmp, filter_dropouts).astype(np.bool)
            )
            x, y, w = x[fil], y[fil], w[fil]

        self._x, self._y, self._w = (
            self._reshape_and_retype(x),
            self._reshape_and_retype(y),
            self._reshape_and_retype(w).squeeze(-1),
        )
        self._x_test = self._reshape_and_retype(x_test)

        if self.x.shape[0] == 0:
            raise RuntimeError("Unable to proceed, no values to fit.")

        if len({self.x.shape[0], self.y.shape[0], self.w.shape[0]}) != 1:
            raise RuntimeError(
                f"Values have different shapes: `{self.x.shape[0]}`, `{self.y.shape[0]}`, `{self.w.shape[0]}`."
            )

        self._gene = gene
        self._lineage = lineage

        return self

    @abstractmethod
    @d.get_sectionsf("base_model_fit", sections=["Parameters"])
    @d.get_full_descriptionf("base_model_fit")
    def fit(
        self,
        x: Optional[np.ndarray] = None,
        y: Optional[np.ndarray] = None,
        w: Optional[np.ndarray] = None,
        **kwargs,
    ) -> "BaseModel":
        """
        Fit the model.

        Parameters
        ----------
        x
            Independent variables, array of shape `(n_samples, 1)`.
        y
            Dependent variables, array of shape `(n_samples, 1)`
        w
            Optional weights of :paramref:`x`, array of shape `(n_samples,)`
        **kwargs
            Keyword arguments for underlying :paramref:`model`'s fitting function.

        Returns
        -------
        :class:`cellrank.ul.models.BaseModel`
            Fits the model and returns self.
        """

        self._check("_x", x)
        self._check("_y", y)
        self._check("_w", w, ndim=1)

        if self._x.shape != self._y.shape:
            raise ValueError(
                f"Inputs and targets differ in shape: `{self._x.shape}` vs. `{self._y.shape}`."
            )
        if self._y.shape[0] != self._w.shape[0]:
            raise ValueError(
                f"Inputs and weights differ in shape: `{self._y.shape[0]}` vs. `{self._w.shape[0]}`."
            )

        return self

    @abstractmethod
    @d.get_sectionsf("base_model_predict", sections=["Parameters", "Returns"])
    @d.get_full_descriptionf("base_model_predict")
    @d.dedent
    def predict(
        self,
        x_test: Optional[np.ndarray] = None,
        key_added: Optional[str] = "_x_test",
        **kwargs,
    ) -> np.ndarray:
        """
        Run the prediction.

        Parameters
        ----------
        x_test
            Array of shape `(n_samples,)` used for prediction.
        key_added
            Attribute name where to save the :paramref:`x_test` for later use. If `None`, don't save it.
        **kwargs
            Keyword arguments for underlying :paramref:`model`'s prediction method.

        Returns
        -------
        :class:`numpy.ndarray`
            Updates and returns the following:

                - :paramref:`y_test` - %(base_model_y_test.summary)s
        """
        pass

    @abstractmethod
    @d.get_sectionsf("base_model_ci", sections=["Parameters", "Returns"])
    @d.get_summaryf("base_model_ci")
    @d.get_full_descriptionf("base_model_ci")
    @d.dedent
    def confidence_interval(
        self, x_test: Optional[np.ndarray] = None, **kwargs
    ) -> np.ndarray:
        """
        Calculate the confidence interval.

        Use :meth:`default_confidence_interval` function if underlying :paramref:`model` has not method
        for confidence interval calculation.

        Parameters
        ----------
        x_test
            Array of shape `(n_samples,)` used for confidence interval calculation.
        **kwargs
            Keyword arguments for underlying :paramref:`model`'s confidence method
            or for :meth:`default_confidence_interval`.

        Returns
        -------
        :class:`numpy.ndarray`
            Updates the following fields:

                - :paramref:`conf_int` - %(base_model_conf_int.summary)s
        """
        pass

    @d.dedent
    def default_confidence_interval(
        self, x_test: Optional[np.ndarray] = None, **kwargs,
    ) -> np.ndarray:
        """
        Calculate the confidence interval, if the underlying :paramref:`model` has no method for it.

        Parameters
        ----------
        %(base_model_ci.parameters)s

        Returns
        -------
        %(base_model_ci.returns)s
                - :paramref:`x_hat` - %(base_model_x_hat.summary)s
                - :paramref:`y_hat` - %(base_model_y_hat.summary)s
        """

        use_ixs = self.w > 0
        x_hat = self.x[use_ixs]
        if x_test is None:
            x_test = self.x_test

        self._y_hat = self.predict(x_hat, key_added="_x_hat", **kwargs)
        self._y_test = self.predict(x_test, key_added="_x_test", **kwargs)

        n = np.sum(use_ixs)
        sigma = np.sqrt(((self.y_hat - self.y[use_ixs].squeeze()) ** 2).sum() / (n - 2))

        stds = (
            np.sqrt(
                1
                + 1 / n
                + ((self.x_test - np.mean(self.x)) ** 2)
                / ((self.x - np.mean(self.x)) ** 2).sum()
            )
            * sigma
            / 2
        )
        stds = np.squeeze(stds)

        self._conf_int = np.c_[self._y_test - stds, self._y_test + stds]

        return self.conf_int

    @d.dedent
    def plot(
        self,
        figsize: Tuple[float, float] = (15, 10),
        same_plot: bool = False,
        hide_cells: bool = False,
        perc: Tuple[float, float] = None,
        abs_prob_cmap: mcolors.ListedColormap = cm.viridis,
        cell_color: str = "black",
        color: str = "black",
        alpha: float = 0.8,
        lineage_alpha: float = 0.2,
        title: Optional[str] = None,
        size: int = 15,
        lw: float = 2,
        show_cbar: bool = True,
        margins: float = 0.015,
        xlabel: str = "pseudotime",
        ylabel: str = "expression",
        show_conf_int: bool = True,
        dpi: int = None,
        fig: mpl.figure.Figure = None,
        ax: mpl.axes.Axes = None,
        return_fig: bool = False,
        save: Optional[str] = None,
    ) -> Optional[mpl.figure.Figure]:
        """
        Plot the smoothed gene expression.

        Parameters
        ----------
        figsize
            Size of the figure.
        same_plot
            Whether to plot all trends in the same plot.
        hide_cells
            Whether to hide the cells.
        perc
            Percentile by which to clip the absorption probabilities.
        abs_prob_cmap
            Colormap to use when coloring in the absorption probabilities.
        cell_color
            Color for the cells when not coloring absorption probabilities.
        color
            Color for the lineages.
        alpha
            Alpha channel for cells.
        lineage_alpha
            Alpha channel for lineage confidence intervals.
        title
            Title of the plot.
        size
            Size of the points.
        lw
            Line width for the smoothed values.
        show_cbar
            Whether to show colorbar.
        margins
            Margins around the plot.
        xlabel
            Label on the x-axis.
        ylabel
            Label on the y-axis.
        show_conf_int
            Whether to show the confidence interval.
        dpi
            Dots per inch.
        fig
            Figure to use, if `None`, create a new one.
        ax: :class:`matplotlib.axes.Axes`
            Ax to use, if `None`, create a new one.
        return_fig
            If `True`, return the figure object.
        save
            Filename where to save the plot. If `None`, just shows the plots.

        Returns
        -------
        %(just_plots)s
        """

        if fig is None or ax is None:
            fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

        if dpi is not None:
            fig.set_dpi(dpi)

        vmin, vmax = _minmax(self.w, perc)
        if not hide_cells:
            _ = ax.scatter(
                self.x_all.squeeze(),
                self.y_all.squeeze(),
                c=cell_color
                if same_plot or np.allclose(self.w_all, 1.0)
                else self.w_all.squeeze(),
                s=size,
                cmap=abs_prob_cmap,
                vmin=vmin,
                vmax=vmax,
                alpha=alpha,
            )

        if title is None:
            title = f"{self._gene} @ {self._lineage}"

        ax.plot(self.x_test, self.y_test, color=color, lw=lw, label=title)

        ax.set_title(title)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        ax.margins(margins)

        if show_conf_int and self.conf_int is not None:
            ax.fill_between(
                self.x_test.squeeze(),
                self.conf_int[:, 0],
                self.conf_int[:, 1],
                alpha=lineage_alpha,
                color=color,
                linestyle="--",
            )

        if show_cbar and not hide_cells and not same_plot:
            norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
            cax, _ = mpl.colorbar.make_axes(ax, aspect=200)
            _ = mpl.colorbar.ColorbarBase(
                cax, norm=norm, cmap=abs_prob_cmap, label="absorption probability"
            )

        if save is not None:
            save_fig(fig, save)

        if return_fig:
            return fig

    def _reshape_and_retype(self, arr: np.ndarray) -> np.ndarray:
        """
        Convert 1D or 2D array to 2D array of shape `(n, 1)` and set the data type.

        Parameters
        ----------
        arr
            Array to convert.

        Returns
        -------
        :class:`np.ndarray`
            Array of shape `(n, 1)` with dtype as :paramref:`_dtype`.
        """

        if arr.ndim not in (1, 2):
            raise ValueError(
                f"Expected array to be 1 or 2 dimensional, found `{arr.ndim}` dimension(s)."
            )
        elif arr.ndim == 2 and arr.shape[1] != 1:
            raise ValueError(
                f"Expected the 2nd dimension to be 1, found `{arr.shape[1]}.`"
            )

        return np.reshape(arr, (-1, 1)).astype(self._dtype)

    def _check(
        self, attr_name: Optional[str], arr: np.ndarray, ndim: int = 2
    ) -> Optional[np.ndarray]:
        """
        Check if the attribute exists with the correct dimension and optionally set it.

        Parameters
        ----------
        attr_name
            Attribute name to check.
        arr
            Value to set. If `None`, just perform the checking.
        ndim
            Expected number of dimensions of the :paramref:`value`.

        Returns
        -------
        :class:`numpy.ndarray`
            The attribute under :paramref:`attr_name`.
        """

        if attr_name is None:
            return
        if arr is None:  # already called prepare
            if not hasattr(self, attr_name):
                raise AttributeError(f"No attribute `{attr_name!r}` found.")
            if getattr(self, attr_name).ndim != ndim:
                raise ValueError(
                    f"Expected attribute `{attr_name!r}` to have `{ndim}` dimensions, "
                    f"found `{getattr(self, attr_name).ndim}` dimensions."
                )
            return getattr(self, attr_name)

        setattr(self, attr_name, self._reshape_and_retype(arr))
        if attr_name.startswith("_"):
            try:
                getattr(self, attr_name[1:])
            except AttributeError:
                setattr(
                    self,
                    attr_name[1:],
                    property(lambda self: getattr(self, attr_name)),
                )
        return getattr(self, attr_name)

    def _copy_attributes(self, dst: "BaseModel") -> None:
        for attr in [
            "_x_all",
            "_y_all",
            "_w_all",
            "_x",
            "_y",
            "_w",
            "_x_test",
            "_y_test",
            "_x_hat",
            "_y_hat",
            "_conf_int",
        ]:
            setattr(dst, attr, _copy(getattr(self, attr)))

    @abstractmethod
    @d.dedent
    def copy(self) -> "BaseModel":  # noqa
        """%(copy)s"""  # noqa
        pass

    def __copy__(self) -> "BaseModel":
        return self.copy()

    def __deepcopy__(self, memodict={}) -> "BaseModel":  # noqa
        res = self.copy()
        res._adata = res.adata.copy()
        self._copy_attributes(res)
        memodict[id(self)] = res
        return res

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        return "{}[{}]".format(
            self.__class__.__name__,
            None
            if self.model is None
            else _dup_spaces.sub(" ", str(self.model).replace("\n", " ")).strip(),
        )


@d.dedent
class SKLearnModel(BaseModel):
    """
    Wrapper around :mod:`sklearn` model.

    Parameters
    ----------
    %(adata)s
    model
        Instance of :mod:`sklearn` model.
    weight_name
        Name of the weight argument for :paramref:`model` `.fit`.
    ignore_raise
        Do not raise an exception if weight argument is not found in the fittng function of :paramref:`model`.
        This is useful in case when weight is passed in `**kwargs` and cannot be determined from signature.
    """

    _fit_names = ("fit", "__init__")
    _predict_names = ("predict", "__call__")
    _weight_names = ("w", "weights", "sample_weight", "sample_weights")
    _conf_int_names = ("conf_int", "confidence_intervals")

    def __init__(
        self,
        adata: AnnData,
        model: BaseEstimator,
        weight_name: Optional[str] = None,
        ignore_raise: bool = False,
    ):
        if not isinstance(model, BaseEstimator):
            raise TypeError(
                f"Expected model to be of type `BaseEstimator`, found `{type(model).__name__!r}`."
            )

        super().__init__(adata, model)

        fit_name = self._find_func(self._fit_names)
        predict_name = self._find_func(self._predict_names)
        ci_name = self._find_func(self._conf_int_names, use_default=True, default=None)
        self._weight_name = None

        if weight_name is None:
            self._weight_name = self._find_arg_name(fit_name, self._weight_names)
        else:
            params = signature(getattr(self.model, fit_name)).parameters
            if not ignore_raise and weight_name not in params:
                raise ValueError(
                    f"Unable to detect `{weight_name!r}` in the signature of `{fit_name!r}`."
                    f"If it's in `kwargs`, set `ignore_raise=True.`"
                )
            self._weight_name = weight_name

        if self._weight_name is None:
            raise RuntimeError(
                f"Unable to determine weights for function `{fit_name!r}`, searched `{self._weight_names}`. "
                f"Consider specifying it manually as `weight_name=...`."
            )

        self._fit_fn = getattr(self.model, fit_name)
        self._pred_fn = getattr(self.model, predict_name)
        self._ci_fn = None if ci_name is None else getattr(self.model, ci_name, None)

    @d.dedent
    def fit(
        self,
        x: Optional[np.ndarray] = None,
        y: Optional[np.ndarray] = None,
        w: Optional[np.ndarray] = None,
        **kwargs,
    ) -> "SKLearnModel":
        """
        %(base_model_fit.full_desc)s

        Parameters
        ----------
        %(base_model_fit.parameters)s

        Returns
        -------
        :class:`cellrank.ul.models.SKLearnModel`
            Fits the model and returns self.
        """  # noqa

        super().fit(x, y, w, **kwargs)

        if self._weight_name is not None:
            kwargs[self._weight_name] = self._w

        self._model = self._fit_fn(self.x, self.y, **kwargs)

        return self

    @d.dedent
    def predict(
        self, x_test: Optional[np.ndarray] = None, key_added: str = "_x_test", **kwargs
    ) -> np.ndarray:
        """
        %(base_model_predict.full_desc)s

        Parameters
        ----------
        %(base_model_predict.parameters)s

        Returns
        -------
        %(base_model_predict.returns)s
        """  # noqa

        x_test = self._check(key_added, x_test)

        self._y_test = self._pred_fn(x_test, **kwargs)
        self._y_test = np.squeeze(self._y_test)

        return self.y_test

    @d.dedent
    def confidence_interval(
        self, x_test: Optional[np.ndarray] = None, **kwargs
    ) -> np.ndarray:
        """
        %(base_model_ci.full_desc)s

        Parameters
        ----------
        %(base_model_ci.parameters)s

        Returns
        -------
        %(base_model_ci.returns)s
        """  # noqa

        if self._ci_fn is None:
            return self.default_confidence_interval(x_test=x_test, **kwargs)

        x_test = self._check("_x_test", x_test)
        self._conf_int = self._ci_fn(x_test, **kwargs)

        return self.conf_int

    def _find_func(
        self,
        func_names: Iterable[str],
        use_default: bool = False,
        default: Optional[str] = None,
    ) -> Optional[str]:
        """
        Find a function in :parmref:`model` from given names.

        If `None` is found, use :parmaref;`default` or raise a :class:`RuntimeError`.

        Parameters
        ----------
        func_names
            Function names to search. The first one found is returned.
        use_default
            Whether to return the default is it's `None` or raise :class:`RuntimeError`.
        default
            The default functional name to use if `None` was found.

        Returns
        -------
        str, None
            Name of the function or the default name.
        """

        for name in func_names:
            if hasattr(self.model, name) and callable(getattr(self.model, name)):
                return name
        if use_default:
            return default
        raise RuntimeError(
            f"Unable to find function and no default specified, searched for `{list(func_names)}`."
        )

    def _find_arg_name(
        self, func_name: Optional[str], param_names: Iterable[str]
    ) -> Optional[str]:
        """
        Find an argument in :paramref:`model`'s :paramref:`func_name`.

        Parameters
        ----------
        func_name
            Function name of :paramref:`model`.
        param_names
            Parameter names to search. The first one found is returned.

        Returns
        -------
        str, None
            The parameter name or `None`, if `None` was found or :paramref:`func_name` was `None`.
        """

        if func_name is None:
            return None

        for param in signature(getattr(self.model, func_name)).parameters:
            if param in param_names:
                return param

        return None

    @property
    def model(self) -> BaseEstimator:
        """The underlying :mod:`sklearn` model."""  # noqa
        return self._model

    @d.dedent
    def copy(self) -> "SKLearnModel":
        """%(copy)s"""  # noqa
        return SKLearnModel(self.adata, deepcopy(self._model))


@d.dedent
class GAM(BaseModel):
    """
    Fit Generalized Additive Models (GAMs) from package :mod:`pygam`.

    Parameters
    ----------
    %(adata)s
    n_splines
        Number of splines.
    spline_order
        Order of the splines.
    distribution
        Name of the distribution. Available distributions can be found
        `here <https://pygam.readthedocs.io/en/latest/notebooks/tour_of_pygam.html#Distribution:>`_.
    link
        Name of the link function. Available functions can be found
        `here <https://pygam.readthedocs.io/en/latest/notebooks/tour_of_pygam.html#Link-function:>`_.
    max_iter
        Maximum number of iterations for optimization.
    expectile
        Expectile for :class:`pygam.pygam.ExpectileGAM`. This forces the distribution to be `'normal'`
        and link function to `'identity'`. Must be in `(0, 1)`.
    use_default_conf_int
        Whether to use :meth:`default_confidence_interval` to calculate the confidence interval or
        use the :paramref:`model`'s method.
    grid
        Whether to perform a grid search. Keys correspond to a parameter names and values to range to be searched.
        If an empty :class:`dict`, don't perform a grid search. If `None`, uses a default grid.
    spline_kwargs
        Keyword arguments for :class:`pygam.terms.s`.
    **kwargs
        Keyword arguments for :class:`pygam.pygam.GAM`.
    """

    def __init__(
        self,
        adata: AnnData,
        n_splines: Optional[int] = 10,
        spline_order: int = 3,
        distribution: str = "gamma",
        link: str = "log",
        max_iter: int = 2000,
        expectile: Optional[float] = None,
        use_default_conf_int: bool = False,
        grid: Optional[Mapping] = MappingProxyType({}),
        spline_kwargs: Mapping = MappingProxyType({}),
        **kwargs,
    ):
        term = s(
            0,
            spline_order=spline_order,
            n_splines=n_splines,
            penalties=["derivative", "l2"],
            **_filter_kwargs(s, **{**{"lam": 3}, **spline_kwargs}),
        )
        link = GamLinkFunction(link)
        distribution = GamDistribution(distribution)
        if distribution == GamDistribution.GAUSS:
            distribution = GamDistribution.NORMAL

        if expectile is not None:
            if not (0 < expectile < 1):
                raise ValueError(
                    f"Expected `expectile` to be in `(0, 1)`, found `{expectile}`."
                )
            if distribution != "normal" or link != "identity":
                logg.warning(
                    f"Expectile GAM works only with `normal` distribution and `identity` link function,"
                    f"found `{distribution!r}` distribution and {link!r} link functions."
                )
            model = ExpectileGAM(
                term, expectile=expectile, max_iter=max_iter, verbose=False, **kwargs
            )
        else:
            gam = _gams[
                distribution, link
            ]  # doing it like this ensure that user can specify scale
            kwargs["link"] = link.s
            kwargs["distribution"] = distribution.s
            model = gam(
                term,
                max_iter=max_iter,
                verbose=False,
                **_filter_kwargs(gam.__init__, **kwargs),
            )
        super().__init__(adata, model=model)
        self._use_default_conf_int = use_default_conf_int
        self._grid = object()  # sentinel value, `None` performs a grid search

        if grid is None:
            self._grid = None
        elif isinstance(grid, (dict, MappingProxyType)):
            if len(grid):
                self._grid = dict(grid)
        else:
            raise TypeError(
                f"Expected `grid` to be `dict` or `None`, found `{type(grid).__name__!r}`."
            )

    @d.dedent
    def fit(
        self,
        x: Optional[np.ndarray] = None,
        y: Optional[np.ndarray] = None,
        w: Optional[np.ndarray] = None,
        **kwargs,
    ) -> "GAM":
        """
        %(base_model_fit.full_desc)s

        Parameters
        ----------
        %(base_model_fit.parameters)s

        Returns
        -------
        :class:`cellrank.ul.models.GAM`
            Fits the model and returns self.
        """  # noqa

        super().fit(x, y, w, **kwargs)

        if self._grid is not None:

            grid = {} if not isinstance(self._grid, dict) else self._grid
            try:
                # workaround for: https://github.com/dswah/pyGAM/issues/273
                self.model.fit(self.x, self.y, weights=self.w, **kwargs)
                self.model.gridsearch(
                    self.x,
                    self.y,
                    weights=self.w,
                    keep_best=True,
                    progress=False,
                    **grid,
                    **kwargs,
                )
                return self
            except Exception as e:
                logg.error(
                    f"Grid search failed, reason: `{e}`. Fitting with default values"
                )

        try:
            self.model.fit(self.x, self.y, weights=self.w, **kwargs)
            return self
        except Exception as e:
            raise RuntimeError(
                f"Unable to fit `{type(self).__name__}` for gene "
                f"`{self._gene!r}` in lineage `{self._lineage!r}`. Reason: `{e}`"
            ) from e

    @d.dedent
    def predict(
        self,
        x_test: Optional[np.ndarray] = None,
        key_added: Optional[str] = "_x_test",
        **kwargs,
    ) -> np.ndarray:
        """
        %(base_model_predict.full_desc)s

        Parameters
        ----------
        %(base_model_predict.parameters)s

        Returns
        -------
        %(base_model_predict.returns)s
        """  # noqa

        x_test = self._check(key_added, x_test)

        self._y_test = self.model.predict(x_test, **kwargs)
        self._y_test = np.squeeze(self._y_test)

        return self.y_test

    @d.dedent
    def confidence_interval(
        self, x_test: Optional[np.ndarray] = None, **kwargs
    ) -> np.ndarray:
        """
        %(base_model_ci.summary)s

        Parameters
        ----------
        %(base_model_ci.parameters)s

        Returns
        -------
        %(base_model_ci.returns)s
        """  # noqa

        x_test = self._check("_x_test", x_test)
        if self._use_default_conf_int:
            self._conf_int = self.default_confidence_interval(x_test=x_test, **kwargs)
        else:
            self._conf_int = self.model.confidence_intervals(x_test, **kwargs)

        return self.conf_int

    @d.dedent
    def copy(self) -> "BaseModel":
        """%(copy)s"""  # noqa
        res = GAM(self.adata)

        res._use_default_conf_int = self._use_default_conf_int
        res._grid = deepcopy(self._grid)
        res._model = deepcopy(self.model)

        return res


@d.dedent
class GAMR(BaseModel):
    """
    Wrapper around R's `mgcv <https://cran.r-project.org/web/packages/mgcv/>`_ or
    `gam <https://cran.r-project.org/web/packages/gam/>`_ package for fitting Generalized Additive Models (GAMs).

    Parameters
    ----------
    %(adata)s
    n_splines
        Number of splines for the GAM.
    smoothing_param
        Smoothing parameter. Increasing this increases the smootheness of splines.
    distribution
        Family in `rpy2.robjects.r`, such as `'gaussian'` or `'poisson'`.
    backend
        R library used to fit GAMs. Valid options are `'mgcv'` and `'gam'`.
        Option `'gam'` ignores the number of splines, as well as family and smoothing parameter.
    """  # noqa

    _fallback_backends = {
        "gam": "mgcv",
        "mgcv": "gam",
    }

    def __init__(
        self,
        adata: AnnData,
        n_splines: int = 5,
        smoothing_param: float = 2,
        distribution: str = "gaussian",
        backend: str = "mgcv",
        perform_import_check: bool = True,
    ):
        super().__init__(adata, model=None)
        self._n_splines = n_splines
        self._sp = smoothing_param
        self._lib = None
        self._lib_name = None
        self._family = distribution

        if backend not in self._fallback_backends.keys():
            raise ValueError(
                f"Invalid backend library `{backend!r}`. Valid options are `{list(self._fallback_backends.keys())}`."
            )

        if (
            perform_import_check
        ):  # it's a bit costly to import, copying just passes the reference
            self._lib, self._lib_name = _maybe_import_r_lib(
                backend
            ) or _maybe_import_r_lib(self._fallback_backends[backend], raise_exc=True)

    @d.dedent
    def fit(
        self,
        x: Optional[np.ndarray] = None,
        y: Optional[np.ndarray] = None,
        w: Optional[np.ndarray] = None,
        **kwargs,
    ) -> "GAMR":
        """
        %(base_model_fit.full_desc)s

        Parameters
        ----------
        %(base_model_fit.parameters)s

        Returns
        -------
        :class:`cellrank.ul.models.GAMR`
            Fits the model and returns self. Updates the following fields by filtering out `0` weights :paramref:`w`:

                - :paramref:`x` - %(base_model_x.summary)s
                - :paramref:`y` - %(base_model_y.summary)s
                - :paramref:`w` - %(base_model_w.summary)s
        """  # noqa

        from rpy2 import robjects
        from rpy2.robjects import pandas2ri, Formula

        super().fit(x, y, w, **kwargs)

        use_ixs = self.w > 0
        self._x = self.x[use_ixs]
        self._y = self.y[use_ixs]
        self._w = self.w[use_ixs]

        family = getattr(robjects.r, self._family, None)
        if family is None:
            family = robjects.r.gaussian

        pandas2ri.activate()
        df = pandas2ri.py2rpy(pd.DataFrame(np.c_[self.x, self.y], columns=["x", "y"]))

        if self._lib_name == "mgcv":
            self._model = self._lib.gam(
                Formula(f'y ~ s(x, k={self._n_splines}, bs="cs")'),
                data=df,
                sp=self._sp,
                family=family,
                weights=pd.Series(self.w),
            )
        elif self._lib_name == "gam":
            self._model = self._lib.gam(
                Formula("y ~ s(x)"),
                data=df,
                sp=self._sp,
                family=family,
                weights=pd.Series(self.w),
            )
        else:
            raise NotImplementedError(
                f"No fitting implemented for R library `{self._lib_name!r}`."
            )

        pandas2ri.deactivate()

        return self

    @d.dedent
    def predict(
        self, x_test: Optional[np.ndarray] = None, key_added: str = "_x_test", **kwargs
    ) -> np.ndarray:
        """
        %(base_model_predict.full_desc)s

        Parameters
        ----------
        %(base_model_predict.parameters)s

        Returns
        -------
        %(base_model_predict.returns)s
        """  # noqa

        from rpy2 import robjects
        from rpy2.robjects import pandas2ri

        if self.model is None:
            raise RuntimeError(
                "Trying to call an uninitialized model. To initialize it, run `.fit()` first."
            )
        if self._lib is None:
            raise RuntimeError(
                f"Unable to fit the model, R package `{self._lib_name}` is not imported."
            )

        x_test = self._check(key_added, x_test)

        pandas2ri.activate()
        self._y_test = (
            np.array(
                robjects.r.predict(
                    self.model,
                    newdata=pandas2ri.py2rpy(pd.DataFrame(x_test, columns=["x"])),
                )
            )
            .squeeze()
            .astype(self._dtype)
        )
        pandas2ri.deactivate()

        return self.y_test

    @d.dedent
    def confidence_interval(
        self, x_test: Optional[np.ndarray] = None, **kwargs
    ) -> np.ndarray:
        """
        %(base_model_ci.summary)s

        This method uses :meth:`default_confidence_interval`.

        Parameters
        ----------
        %(base_model_ci.parameters)s

        Returns
        -------
        %(base_model_ci.returns)s
        """  # noqa
        return self.default_confidence_interval(x_test=x_test, **kwargs)

    @d.dedent
    def copy(self) -> "GAMR":
        """%(copy)s"""  # noqa
        res = GAMR(
            self.adata,
            self._n_splines,
            self._sp,
            distribution=self._family,
            perform_import_check=False,
        )
        res._lib = self._lib
        res._lib_name = self._lib_name
        return res

    def __getstate__(self) -> dict:
        return {k: v for k, v in self.__dict__.items() if k != "_lib"}

    def __setstate__(self, state: dict):
        self.__dict__ = state
        self._lib, self._lib_name = _maybe_import_r_lib(self._lib_name, raise_exc=True)


def _maybe_import_r_lib(
    name: str, raise_exc: bool = False
) -> Optional[Tuple[Any, str]]:
    global _r_lib, _r_lib_name

    if name == _r_lib_name and _r_lib is not None:
        return _r_lib, _r_lib_name

    try:
        from logging import ERROR
        from rpy2.robjects import r
        from rpy2.rinterface_lib.callbacks import logger
        from rpy2.robjects.packages import PackageNotInstalledError, importr

        logger.setLevel(ERROR)
        r["options"](warn=-1)

    except ImportError as e:
        raise ImportError(
            "Unable to import `rpy2`, install it first as `pip install rpy2`."
        ) from e

    try:
        _r_lib = importr(name)
        _r_lib_name = name
        return _r_lib, _r_lib_name
    except PackageNotInstalledError as e:
        if not raise_exc:
            return
        raise RuntimeError(
            f"Install R library `{name!r}` first as `install.packages({name!r}).`"
        ) from e
