import abc
import collections
import copy
import enum
import re
import warnings
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple, Union

import wrapt

import numpy as np
import pandas as pd
import scipy.sparse as sp
from pandas.api.types import infer_dtype
from scipy.ndimage import convolve

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

from anndata import AnnData
from scanpy.plotting._utils import add_colors_for_categorical_sample_annotation

from cellrank import logging as logg
from cellrank._utils._docs import d
from cellrank._utils._enum import ModeEnum
from cellrank._utils._lineage import Lineage
from cellrank._utils._utils import _densify_squeeze, _minmax, save_fig, valuedispatch
from cellrank.kernels.mixins import IOMixin

__all__ = ["BaseModel"]

_dup_spaces = re.compile(r" +")  # used on repr for underlying model's repr
ArrayLike = Union[np.ndarray, sp.spmatrix, List, Tuple]


class UnknownModelError(RuntimeError):
    pass


class FailedReturnType(ModeEnum):
    PREPARE = enum.auto()
    FIT = enum.auto()
    PREDICT = enum.auto()
    CONFIDENCE_INTERVAL = enum.auto()
    DEFAULT_CONFIDENCE_INTERVAL = enum.auto()
    PLOT = enum.auto()


class ColorType(ModeEnum):
    CONT = enum.auto()
    CAT = enum.auto()
    STR = enum.auto()


def _handle_exception(return_type: FailedReturnType, func: Callable) -> Callable:
    def handle(*, exception_handler: Callable):
        @wrapt.decorator
        def wrapper(wrapped, instance, args, kwargs):
            try:
                if isinstance(instance, FailedModel):
                    instance.reraise()
                return wrapped(*args, **kwargs)
            except Exception as e:  # noqa: BLE001
                return exception_handler(instance, e, *args, **kwargs)

        return wrapper

    def array_output(instance: "BaseModel", exc: BaseException, *_args, **kwargs) -> np.ndarray:
        if not instance._is_bulk:
            raise exc from None

        if kwargs.get("x_test", None) is not None:
            n_obs = kwargs["x_test"].shape[0]
        elif instance.x_test is not None:
            n_obs = instance.x_test.shape[0]
        else:
            n_obs = 200

        return np.full((n_obs, array_shape), np.nan, dtype=instance._dtype)

    def model_output(instance: "BaseModel", exc: BaseException, *_args, **_kwargs) -> "FailedModel":
        if not isinstance(instance, FailedModel):
            instance = FailedModel(instance, exc=exc)
        if not instance._is_bulk:
            instance.reraise()

        return instance

    def no_output(instance: "BaseModel", exc: BaseException, *_args, **_kwargs) -> None:
        if not instance._is_bulk:
            raise exc from None

    @valuedispatch
    def wrapper(mode: FailedReturnType, *_args, **_kwargs):
        raise NotImplementedError(mode)

    @wrapper.register(FailedReturnType.PREPARE)
    @wrapper.register(FailedReturnType.FIT)
    def _(func: Callable) -> Callable:
        return handle(exception_handler=model_output)(func)

    @wrapper.register(FailedReturnType.PREDICT)
    @wrapper.register(FailedReturnType.CONFIDENCE_INTERVAL)
    @wrapper.register(FailedReturnType.DEFAULT_CONFIDENCE_INTERVAL)
    def _(func: Callable) -> Callable:
        return handle(exception_handler=array_output)(func)

    @wrapper.register(FailedReturnType.PLOT)
    def _(func: Callable):
        return handle(exception_handler=no_output)(func)

    return_type = FailedReturnType(return_type)
    array_shape = 1 + (
        return_type
        in (
            FailedReturnType.CONFIDENCE_INTERVAL,
            FailedReturnType.DEFAULT_CONFIDENCE_INTERVAL,
        )
    )

    return wrapper(return_type, func)


class BaseModelMeta(abc.ABCMeta):
    """Metaclass for all base models."""

    def __new__(cls, clsname: str, superclasses: Tuple[type, ...], attributedict: Dict[str, Any]):
        """Create a new instance.

        Parameters
        ----------
        clsname
            Name of class to be constructed.
        superclasses
            List of superclasses.
        attributedict
            Dictionary of attributes.
        """
        obj = super().__new__(cls, clsname, superclasses, attributedict)
        for fun_name in list(FailedReturnType):
            setattr(
                obj,
                fun_name,
                _handle_exception(FailedReturnType(fun_name), getattr(obj, fun_name)),
            )

        return obj


@d.get_sections(base="base_model", sections=["Parameters"])
@d.dedent
class BaseModel(IOMixin, abc.ABC, metaclass=BaseModelMeta):
    """Base class for all model classes.

    Parameters
    ----------
    %(adata)s
    model
        The underlying model that is used for fitting and prediction.
    """

    _dtype = np.float64

    def __init__(
        self,
        adata: Optional[AnnData],
        model: Any,
    ):
        if not isinstance(adata, AnnData) and not isinstance(self, FittedModel):
            # FittedModel doesn't need it
            raise TypeError(f"Expected `adata` to be of type `anndata.AnnData`, found `{type(adata).__name__}`.")
        super().__init__()

        self._adata = adata
        self._n_obs = 0 if adata is None else adata.n_obs
        self._model = model
        self._gene = None
        self._use_raw = False
        self._data_key = None
        self._lineage = None
        self._prepared = False

        self._obs_names = None
        self._is_bulk = False

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

    @d.get_summary(base="base_model_prepared")
    @property
    def prepared(self):
        """Whether the model is prepared for fitting."""
        return self._prepared

    @property
    @d.dedent
    def adata(self) -> AnnData:
        """Annotated data object."""
        return self._adata

    @adata.setter
    def adata(self, adata: Optional[AnnData]) -> None:
        self._adata = adata

    @property
    def shape(self) -> Tuple[int]:
        """Number of cells in :attr:`adata`."""
        return (self._n_obs,)

    @property
    def model(self) -> Any:
        """Underlying model."""
        return self._model

    @property
    @d.get_summary(base="base_model_x_all")
    def x_all(self) -> np.ndarray:
        """Unfiltered independent variables of shape ``(n_cells, 1)``."""
        return self._x_all

    @property
    @d.get_summary(base="base_model_y_all")
    def y_all(self) -> np.ndarray:
        """Unfiltered dependent variables of shape ``(n_cells, 1)``."""
        return self._y_all

    @property
    @d.get_summary(base="base_model_w_all")
    def w_all(self) -> np.ndarray:
        """Unfiltered weights of shape ``(n_cells,)``."""
        return self._w_all

    @property
    @d.get_summary(base="base_model_x")
    def x(self) -> np.ndarray:
        """Filtered independent variables of shape ``(n_filtered_cells, 1)`` used for fitting."""
        return self._x

    @property
    @d.get_summary(base="base_model_y")
    def y(self) -> np.ndarray:
        """Filtered dependent variables of shape ``(n_filtered_cells, 1)`` used for fitting."""
        return self._y

    @property
    @d.get_summary(base="base_model_w")
    def w(self) -> np.ndarray:
        """Filtered weights of shape ``(n_filtered_cells,)`` used for fitting."""
        return self._w

    @property
    @d.get_summary(base="base_model_x_test")
    def x_test(self) -> np.ndarray:
        """Independent variables of shape ``(n_samples, 1)`` used for prediction."""
        return self._x_test

    @property
    @d.get_summary(base="base_model_y_test")
    def y_test(self) -> np.ndarray:
        """Prediction values of shape ``(n_samples,)`` for :attr:`x_test`."""
        return self._y_test

    @property
    @d.get_summary(base="base_model_x_hat")
    def x_hat(self) -> np.ndarray:
        """Filtered independent variables used when calculating default confidence interval, usually same as :attr:`x`."""  # noqa: E501
        return self._x_hat

    @property
    @d.get_summary(base="base_model_y_hat")
    def y_hat(self) -> np.ndarray:
        """Filtered dependent variables used when calculating default confidence interval, usually same as :attr:`y`."""
        return self._y_hat

    @property
    @d.get_summary(base="base_model_conf_int")
    def conf_int(self) -> np.ndarray:
        """Array of shape ``(n_samples, 2)`` containing the lower and upper bound of the confidence interval."""
        return self._conf_int

    @d.get_sections(base="base_model_prepare", sections=["Parameters", "Returns"])
    @d.get_full_description(base="base_model_prepare")
    @d.dedent
    def prepare(
        self,
        gene: str,
        lineage: Optional[str],
        time_key: str,
        backward: bool = False,
        time_range: Optional[Union[float, Tuple[float, float]]] = None,
        data_key: Optional[str] = "X",
        use_raw: bool = False,
        threshold: Optional[float] = None,
        weight_threshold: Union[float, Tuple[float, float]] = (0.01, 0.01),
        filter_cells: Optional[float] = None,
        n_test_points: int = 200,
    ) -> "BaseModel":
        """Prepare the model to be ready for fitting.

        Parameters
        ----------
        gene
            Gene in :attr:`~anndata.AnnData.var_names`.
        lineage
            Name of the lineage. If :obj:`None`, all weights will be set to :math:`1`.
        time_key
            Key in :attr:`~anndata.AnnData.obs` where the pseudotime is stored.
        %(backward)s
        %(time_range)s
        data_key
            Key in :attr:`~anndata.AnnData.layers` or ``'X'`` for :attr:`~anndata.AnnData.X`.
            If ``use_raw = True``, it's always set to ``'X'``.
        use_raw
            Whether to access :attr:`~anndata.AnnData.raw`.
        threshold
            Consider only cells with weights > ``threshold`` when estimating the test endpoint.
            If :obj:`None`, use the median of the weights.
        weight_threshold
            Set all weights below ``weight_threshold`` to ``weight_threshold`` if a :class:`float`,
            or to the second value, if a :class:`tuple`.
        filter_cells
            Filter out all cells with expression values lower than this threshold.
        n_test_points
            Number of test points. If :obj:`None`, use the original points based on ``threshold``.

        Returns
        -------
        Nothing, just updates the following fields:

        - :attr:`x` - %(base_model_x.summary)s
        - :attr:`y` - %(base_model_y.summary)s
        - :attr:`w` - %(base_model_w.summary)s
        - :attr:`x_all` - %(base_model_x_all.summary)s
        - :attr:`y_all` - %(base_model_y_all.summary)s
        - :attr:`w_all` - %(base_model_w_all.summary)s
        - :attr:`x_test` - %(base_model_x_test.summary)s
        - :attr:`prepared` - %(base_model_prepared.summary)s
        """
        if use_raw:
            if self.adata.raw is None:
                raise AttributeError("AnnData object has no attribute `.raw`.")
            if data_key != "X":
                data_key = "X"
        self._use_raw = use_raw

        if data_key not in ["X", "obs", None] + list(self.adata.layers.keys()):
            raise KeyError(
                f"Data key must be a key of `adata.layers`: `{list(self.adata.layers.keys())}`, "
                f"`adata.X` or `adata.obs`."
            )

        probs = Lineage.from_adata(self.adata, backward=backward)
        if lineage is not None:
            probs = probs[lineage]

        if time_key not in self.adata.obs:
            raise KeyError(f"Time key `{time_key!r}` not found in `adata.obs`.")

        if data_key == "obs":
            if gene not in self.adata.obs:
                raise KeyError(f"Unable to find data in `adata.obs[{gene!r}]`.")
        else:
            if use_raw and gene not in self.adata.raw.var_names:
                raise KeyError(f"Gene `{gene!r}` not found in `adata.raw.var_names`.")
            if not use_raw and gene not in self.adata.var_names:
                raise KeyError(f"Gene `{gene!r}` not found in `adata.var_names`.")

        if not isinstance(time_range, (type(None), float, int, tuple)):
            raise TypeError(
                f"Expected time range to be either `None`, "
                f"`float` or `tuple`, found `{type(time_range).__name__!r}`."
            )
        if isinstance(time_range, tuple):
            if len(time_range) != 2:
                raise ValueError(f"Expected time range to be a `tuple` of length `2`, found `{len(time_range)}`.")
            if not isinstance(time_range[0], (float, int, type(None))) or not isinstance(
                time_range[1], (float, int, type(None))
            ):
                raise TypeError(
                    f"Expected values in time ranges be of types `float` or `int` or `None`, "
                    f"got `{type(time_range[0]).__name__!r}` and `{type(time_range[1]).__name__!r}`."
                )
            val_start, val_end = time_range
        elif time_range is None:
            val_start, val_end = None, None
        else:
            val_start, val_end = None, time_range

        if isinstance(weight_threshold, (int, float, np.number)):
            weight_threshold = (weight_threshold, weight_threshold)
        if len(weight_threshold) != 2:
            raise ValueError(f"Expected `weight_threshold` to be of size `2`, found `{len(weight_threshold)}`.")

        self._obs_names = self.adata.obs_names.values[:]

        x = np.array(self.adata.obs[time_key]).astype(self._dtype)

        adata = self.adata.raw.to_adata() if use_raw else self.adata
        gene_ix = np.where(adata.var_names == gene)[0]

        if data_key in ("X", None):
            y = adata.X[:, gene_ix]
            self._data_key = None
        elif data_key == "obs":
            y = adata.obs[gene].values
            # don't set data_key, it's just when `cell_color = ...`
        elif data_key in adata.layers:
            y = adata.layers[data_key][:, gene_ix]
            self._data_key = data_key
        else:
            raise NotImplementedError(f"Data key `{data_key!r}` is not implemented.")

        if lineage is not None:
            weight_threshold, val = weight_threshold
            w = _densify_squeeze(probs.X, self._dtype)
            w[w < weight_threshold] = val
        else:
            w = np.ones(len(x), dtype=self._dtype)

        if use_raw:
            correct_ixs = np.isin(self.adata.obs_names, adata.obs_names)
            x = x[correct_ixs]
            y = y[correct_ixs]
            w = w[correct_ixs]
            self._obs_names = self._obs_names[correct_ixs]

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
        self._obs_names = self._obs_names[ixs]

        if np.allclose(w, w[0]):
            # degenerate case
            val_start = w[0] - 1
            val_end = w[0] + 1

        if val_start is None:
            val_start = np.min(x)
        if val_end is None:
            if threshold is None:
                threshold = np.nanmedian(w)
            # use `>=` because weights can all be 1
            w_test = w[w >= threshold]
            n_window = n_test_points // 20
            tmp = convolve(w_test, np.ones(n_window) / n_window, mode="nearest")
            val_end = x[w >= threshold][-1 if lineage is None else np.nanargmax(tmp)]

        if val_start > val_end:
            val_start, val_end = val_end, val_start
        val_start, val_end = (max(val_start, np.min(x)), min(val_end, np.max(x)))

        fil = (x >= val_start) & (x <= val_end)
        x_test = np.linspace(val_start, val_end, n_test_points) if n_test_points is not None else x[fil]
        x, y, w = x[fil], y[fil], w[fil]
        self._obs_names = self._obs_names[fil]

        if filter_cells is not None:
            tmp = y.squeeze()
            fil = (tmp >= filter_cells) & (~np.isclose(tmp, filter_cells).astype(bool))
            x, y, w = x[fil], y[fil], w[fil]
            self._obs_names = self._obs_names[fil]

        self._x, self._y, self._w = (
            self._reshape_and_retype(x),
            self._reshape_and_retype(y),
            self._reshape_and_retype(w).squeeze(-1),
        )
        self._x_test = self._reshape_and_retype(x_test)
        self._y_test = None
        self._x_hat = None
        self._y_hat = None
        self._conf_int = None

        if self.x.shape[0] == 0:
            raise RuntimeError("Unable to proceed, no values to fit.")

        if len({self.x.shape[0], self.y.shape[0], self.w.shape[0]}) != 1:
            raise RuntimeError(
                f"Values have different shapes: `{self.x.shape[0]}`, `{self.y.shape[0]}`, `{self.w.shape[0]}`."
            )

        self._gene = gene
        self._lineage = lineage
        self._prepared = True

        return self

    @abc.abstractmethod
    @d.get_sections(base="base_model_fit", sections=["Parameters"])
    @d.get_full_description(base="base_model_fit")
    def fit(
        self,
        x: Optional[np.ndarray] = None,
        y: Optional[np.ndarray] = None,
        w: Optional[np.ndarray] = None,
        **kwargs: Any,
    ) -> "BaseModel":
        """Fit the model.

        Parameters
        ----------
        x
            Independent variables, array of shape ``(n_samples, 1)``. If :obj:`None`, use :attr:`x`.
        y
            Dependent variables, array of shape ``(n_samples, 1)``. If :obj:`None`, use :attr:`y`.
        w
            Optional weights of :attr:`x`, array of shape ``(n_samples,)``. If :obj:`None`, use :attr:`w`.
        kwargs
            Keyword arguments for underlying :attr:`model`'s fitting function.

        Returns
        -------
        Fits the :attr:`model` and returns self.
        """
        if not self.prepared:
            raise RuntimeError("The model has not been prepared yet, call `.prepare()` first.")

        self._check("_x", x)
        self._check("_y", y)
        self._check("_w", w, ndim=1)

        if self._x.shape != self._y.shape:
            raise ValueError(f"Inputs and targets differ in shape: `{self._x.shape}` vs. `{self._y.shape}`.")
        if self._y.shape[0] != self._w.shape[0]:
            raise ValueError(f"Inputs and weights differ in shape: `{self._y.shape[0]}` vs. `{self._w.shape[0]}`.")

        return self

    @abc.abstractmethod
    @d.get_sections(base="base_model_predict", sections=["Parameters", "Returns"])
    @d.get_full_description(base="base_model_predict")
    @d.dedent
    def predict(
        self,
        x_test: Optional[np.ndarray] = None,
        key_added: Optional[str] = "_x_test",
        **kwargs: Any,
    ) -> np.ndarray:
        """Run the prediction.

        Parameters
        ----------
        x_test
            Array of shape ``(n_samples,)`` used for prediction. If :obj:`None`, use :attr:`x_test`.
        key_added
            Attribute name where to save the :attr:`x_test` for later use. If :obj:`None`, don't save it.
        kwargs
            Keyword arguments for underlying :attr:`model`'s prediction method.

        Returns
        -------
        Returns and updates the following fields:

        - :attr:`y_test` - %(base_model_y_test.summary)s
        """

    @abc.abstractmethod
    @d.get_sections(base="base_model_ci", sections=["Parameters", "Returns"])
    @d.get_summary(base="base_model_ci")
    @d.get_full_description(base="base_model_ci")
    @d.dedent
    def confidence_interval(self, x_test: Optional[np.ndarray] = None, **kwargs) -> np.ndarray:
        """Calculate the confidence interval.

        Use :meth:`default_confidence_interval` function if underlying :attr:`model` has no method
        for confidence interval calculation.

        Parameters
        ----------
        x_test
            Array of shape ``(n_samples,)`` used for confidence interval calculation.
            If :obj:`None`, use :attr:`x_test`.
        kwargs
            Keyword arguments for underlying :attr:`model`'s confidence method
            or for :meth:`default_confidence_interval`.

        Returns
        -------
        Returns self and updates the following fields:

        - :attr:`conf_int` - %(base_model_conf_int.summary)s
        """

    @d.dedent
    def default_confidence_interval(
        self,
        x_test: Optional[np.ndarray] = None,
        **kwargs: Any,
    ) -> np.ndarray:
        """Calculate the confidence interval, if the underlying :attr:`model` has no method for it.

        This formula is taken from :cite:`desalvo:70`, eq. 5.

        Parameters
        ----------
        %(base_model_ci.parameters)s

        Returns
        -------
        %(base_model_ci.returns)s

        Also updates the following fields:

        - :attr:`x_hat` - %(base_model_x_hat.summary)s
        - :attr:`y_hat` - %(base_model_y_hat.summary)s
        """
        use_ixs = self.w > 0
        x_hat = self.x[use_ixs]
        if x_test is None:
            x_test = self.x_test

        self._y_hat = self.predict(x_hat, key_added="_x_hat", **kwargs)
        self._y_test = self.predict(x_test, key_added="_x_test", **kwargs)

        n = np.sum(use_ixs)
        sigma_hat = np.sqrt(((self.y_hat - self.y[use_ixs].squeeze()) ** 2).sum() / (n - 2))
        mean = np.mean(self.x)

        # fmt: off
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            stds = sigma_hat * np.sqrt(1 + 1 / n + ((self.x_test - mean) ** 2) / ((self.x - mean) ** 2).sum())
        # fmt: on
        stds = np.squeeze(stds)

        self._conf_int = np.c_[self._y_test - stds / 2.0, self._y_test + stds / 2.0]

        return self.conf_int

    @d.dedent
    def plot(
        self,
        figsize: Tuple[float, float] = (8, 5),
        same_plot: bool = False,
        hide_cells: bool = False,
        perc: Tuple[float, float] = None,
        fate_prob_cmap: colors.ListedColormap = cm.viridis,
        cell_color: Optional[str] = None,
        lineage_color: str = "black",
        alpha: float = 0.8,
        lineage_alpha: float = 0.2,
        title: Optional[str] = None,
        size: int = 15,
        lw: float = 2,
        cbar: bool = True,
        margins: float = 0.015,
        xlabel: str = "pseudotime",
        ylabel: str = "expression",
        conf_int: bool = True,
        lineage_probability: bool = False,
        lineage_probability_conf_int: Union[bool, float] = False,
        lineage_probability_color: Optional[str] = None,
        obs_legend_loc: Optional[str] = "best",
        dpi: int = None,
        fig: mpl.figure.Figure = None,
        ax: mpl.axes.Axes = None,
        return_fig: bool = False,
        save: Optional[str] = None,
        **kwargs: Any,
    ) -> Optional[mpl.figure.Figure]:
        """Plot the smoothed gene expression.

        Parameters
        ----------
        figsize
            Size of the figure.
        same_plot
            Whether to plot all trends in the same plot.
        hide_cells
            Whether to hide the cells.
        perc
            Percentile by which to clip the fate probabilities.
        fate_prob_cmap
            Colormap to use when coloring in the fate probabilities.
        cell_color
            Key in :attr:`~anndata.AnnData.obs` or :attr:`~anndata.AnnData.var_names` used for coloring the cells.
        lineage_color
            Color for the lineage.
        alpha
            Alpha value in :math:`[0, 1]` for the transparency of cells.
        lineage_alpha
            Alpha value in :math:`[0, 1]` for the transparency lineage confidence intervals.
        title
            Title of the plot.
        size
            Size of the points.
        lw
            Line width for the smoothed values.
        cbar
            Whether to show the colorbar.
        margins
            Margins around the plot.
        xlabel
            Label on the x-axis.
        ylabel
            Label on the y-axis.
        conf_int
            Whether to show the confidence interval.
        lineage_probability
            Whether to show smoothed lineage probability as a dashed line.
            Note that this will require 1 additional model fit.
        lineage_probability_conf_int
            Whether to compute and show smoothed lineage probability confidence interval.
        lineage_probability_color
            Color to use when plotting the smoothed ``lineage_probability``.
            If :obj:`None`, it's the same as ``lineage_color``. Only used when ``show_lineage_probability = True``.
        obs_legend_loc
            Location of the legend when ``cell_color`` corresponds to a categorical variable.
        dpi
            Dots per inch.
        fig
            Figure to use. If :obj:`None`, create a new one.
        ax
            Ax to use. If :obj:`None`, create a new one.
        return_fig
            If :obj:`True`, return the figure object.
        save
            Filename where to save the plot. If :obj:`None`, just shows the plots.
        kwargs
            Keyword arguments for :meth:`~matplotlib.axes.Axes.legend`.

        Returns
        -------
        %(just_plots)s
        """
        if self.y_test is None:
            raise RuntimeError("Run `.predict()` first.")

        if fig is None or ax is None:
            fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

        if dpi is not None:
            fig.set_dpi(dpi)

        conf_int = conf_int and self.conf_int is not None
        hide_cells = hide_cells or self.x_all is None or self.w_all is None or self.y_all is None

        lineage_probability_color = lineage_color if lineage_probability_color is None else lineage_probability_color

        scaler = kwargs.pop(
            "scaler",
            self._create_scaler(
                lineage_probability,
                show_conf_int=conf_int,
            ),
        )

        if lineage_probability and ylabel in ("expression", self._gene):
            ylabel = f"scaled {ylabel}"

        vmin, vmax = None, None
        key, color, typp, mapper = self._get_colors(cell_color, same_plot=same_plot)

        if not hide_cells:
            cbar = cbar and (typp == ColorType.CONT)
            if typp == ColorType.CONT:
                vmin, vmax = _minmax(color, perc)

            _ = ax.scatter(
                self.x_all.squeeze(),
                scaler(self.y_all.squeeze()),
                c=color,
                s=size,
                cmap=fate_prob_cmap,
                vmin=vmin,
                vmax=vmax,
                alpha=alpha,
            )

        if title is None:
            title = f"{self._gene} @ {self._lineage}" if self._lineage is not None else f"{self._gene}"

        ax.plot(self.x_test, scaler(self.y_test), color=lineage_color, lw=lw, label=title)

        if title is not None:
            ax.set_title(title)
        if ylabel is not None:
            ax.set_ylabel(ylabel)
        if xlabel is not None:
            ax.set_xlabel(xlabel)

        ax.margins(margins)

        if conf_int:
            ax.fill_between(
                self.x_test.squeeze(),
                scaler(self.conf_int[:, 0]),
                scaler(self.conf_int[:, 1]),
                alpha=lineage_alpha,
                color=lineage_color,
                linestyle="--",
            )

        if lineage_probability and not isinstance(self, FittedModel) and not np.allclose(self.w, 1.0):
            from cellrank.pl._utils import _is_any_gam_mgcv

            model = copy.deepcopy(self)
            model._y = self._reshape_and_retype(self.w).copy()
            model = model.fit()

            if not lineage_probability_conf_int:
                y = model.predict()
            elif _is_any_gam_mgcv(model):
                y = model.predict(
                    level=lineage_probability_conf_int if isinstance(lineage_probability_conf_int, float) else 0.95
                )
            else:
                y = model.predict()
                model.confidence_interval()

                ax.fill_between(
                    model.x_test.squeeze(),
                    model.conf_int[:, 0],
                    model.conf_int[:, 1],
                    alpha=lineage_alpha,
                    color=lineage_probability_color,
                    linestyle="--",
                )

            handle = ax.plot(
                model.x_test,
                y,
                color=lineage_probability_color,
                lw=lw,
                linestyle="--",
                zorder=-1,
                label="probability",
            )

            self._maybe_add_legend(fig, ax, mapper=handle, title=None, **kwargs)

        if cbar and not hide_cells and not same_plot and not np.allclose(self.w_all, 1.0):
            norm = colors.Normalize(vmin=vmin, vmax=vmax)
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="2%", pad=0.1)
            _ = mpl.colorbar.ColorbarBase(
                cax,
                norm=norm,
                cmap=fate_prob_cmap,
                ticks=np.linspace(norm.vmin, norm.vmax, 5),
            )
            cax.set_ylabel(key)
        elif typp == ColorType.CAT:
            self._maybe_add_legend(fig, ax, mapper, title=key, loc=obs_legend_loc, is_line=False)

        if save is not None:
            save_fig(fig, save)

        if return_fig:
            return fig

    def _maybe_add_legend(
        self,
        fig: mpl.figure.Figure,
        ax: mpl.axes.Axes,
        mapper: Union[Sequence[Any], Mapping[str, Any]],
        title: Optional[str] = None,
        loc: Optional[str] = "best",
        is_line: bool = True,
        **kwargs: Any,
    ) -> None:
        from cellrank.pl._utils import _position_legend

        if loc in ("none", None):
            return

        if isinstance(mapper, dict):
            handles = [
                ax.plot([], [], label=name, color=color)[0] if is_line else ax.scatter([], [], label=name, color=color)
                for name, color in mapper.items()
            ]
        else:
            handles = mapper

        legend = _position_legend(
            ax,
            legend_loc=loc,
            handles=handles,
            title=title,
            **kwargs,
        )
        fig.add_artist(legend)

    def _reshape_and_retype(self, arr: np.ndarray) -> np.ndarray:
        """Convert 1D or 2D array to 2D array of shape ``(n, 1)`` and set the data type.

        Parameters
        ----------
        arr
            Array to convert.

        Returns
        -------
        Array of shape ``(n, 1)`` with dtype set to :attr:`_dtype`.
        """
        if arr.ndim not in (1, 2):
            raise ValueError(f"Expected array to be 1 or 2 dimensional, found `{arr.ndim}` dimension(s).")
        if arr.ndim == 2 and arr.shape[1] != 1:
            raise ValueError(f"Expected the 2nd dimension to be 1, found `{arr.shape[1]}.`")

        return np.reshape(arr, (-1, 1)).astype(self._dtype)

    def _check(self, attr_name: Optional[str], arr: Optional[np.ndarray], ndim: int = 2) -> Optional[np.ndarray]:
        """Check if the attribute exists with the correct dimension and optionally set it.

        Parameters
        ----------
        attr_name
            Attribute name to check.
        arr
            Value to set. If :obj:`None`, just perform the checking.
        ndim
            Expected number of dimensions of the ``arr``.

        Returns
        -------
        The attribute under ``attr_name``.
        """
        if attr_name is None:
            return
        if arr is None:  # already called prepare
            if not hasattr(self, attr_name):
                raise AttributeError(f"No attribute `{attr_name!r}` found.")
            if not isinstance(getattr(self, attr_name, None), np.ndarray):
                raise AttributeError(
                    f"Expected `{attr_name!r}` to be `numpy.ndarray`, "
                    f"found `{type(getattr(self, attr_name, None)).__name__!r}`."
                )
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

    def _deepcopy_attributes(self, dst: "BaseModel") -> None:
        # __deepcopy__ will usually call `_shallowcopy_attributes` twice, since it calls `.copy()`,
        # which should always call it (it copies the lineage and gene names + deepcopies the model)

        self._shallowcopy_attributes(dst)
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
            "_prepared",
        ]:
            setattr(dst, attr, copy.copy(getattr(self, attr)))

    def _shallowcopy_attributes(self, dst: "BaseModel") -> None:
        for attr in ["_gene", "_lineage", "_use_raw", "_data_key", "_is_bulk"]:
            setattr(dst, attr, copy.copy(getattr(self, attr)))
        # user is not exposed to this
        dst._obs_names = self._obs_names
        # always deepcopy the model (we're manipulating it in multiple threads/processed)
        dst._model = copy.deepcopy(self.model)

    @abc.abstractmethod
    @d.dedent
    def copy(self) -> "BaseModel":
        """%(copy)s"""  # noqa

    def __copy__(self) -> "BaseModel":
        return self.copy()

    def __deepcopy__(self, memodict={}) -> "BaseModel":  # noqa
        # deepcopy expects that `.copy()` makes a really shallow copy (i.e. only references to the arrays)
        # it should also not copy the `.prepared` attribute, since copying is happening mostly during
        # parallelization and it serves as 1 extra sanity check (a precaution that's not necessary, per-se, but highly
        # desirable)
        res = self.copy()
        # we don't copy the adata object for 2 reasons:
        # 1. these objects are meant to be as lightweight as possible
        # 2. in `.plot`, we deepcopy the model when plotting smoothed probabilities
        self._deepcopy_attributes(res)
        memodict[id(self)] = res

        return res

    def __bool__(self):
        return True

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        return "<{}[gene={!r}, lineage={!r}, model={}]>".format(
            self.__class__.__name__,
            self._gene,
            self._lineage,
            None if self.model is None else _dup_spaces.sub(" ", str(self.model).replace("\n", " ")).strip(),
        )

    def _get_colors(
        self,
        key: Optional[str],
        *,
        same_plot: bool = False,
    ) -> Tuple[Optional[str], Optional[Union[str, np.ndarray]], ColorType, Optional[Dict[str, Any]],]:
        """
        Get color array.

        Parameters
        ----------
        key
            Key in :attr:`~anndata.AnnData.obs`, :attr:`~anndata.AnnData.var_names` or
            :attr:`~anndata.AnnData.raw.var_names`. Search first starts in `.obs`, then `.raw.var_names` and lastly
            `.layers`, using :meth:`~anndata.AnnData.obs_vector`.

        Returns
        -------
        Triple of the following:

        - name of the colorbar label.
        - array of values to map.
        - type of the color.
        - color mapper for categorical colors.
        """
        if colors.is_color_like(key):
            return None, key, ColorType.STR, None

        if key in self.adata.obs:
            if isinstance(self.adata.obs[key].dtype, pd.CategoricalDtype):
                add_colors_for_categorical_sample_annotation(
                    self.adata,
                    key=key,
                    force_update_colors=False,
                    palette=None,
                )
                col_dict = collections.defaultdict(
                    lambda: colors.to_rgb("grey"),
                    zip(
                        self.adata.obs[key].cat.categories,
                        [colors.to_rgb(i) for i in self.adata.uns[f"{key}_colors"]],
                    ),
                )
                return (
                    key,
                    np.array([col_dict[v] for v in self.adata.obs[key]]),
                    ColorType.CAT,
                    col_dict,
                )

            if np.issubdtype(self.adata.obs[key].dtype, np.number):
                return key, self.adata.obs[key].values, ColorType.CONT, None

            logg.debug(f"Unable to interpret cell color from type `{infer_dtype(self.adata.obs[key])}`")
            return None, "black", ColorType.STR, None

        if self._use_raw and key in self.adata.raw.var_names:
            return (
                key,
                _densify_squeeze(self.adata.raw[:, key], np.float64),
                ColorType.CONT,
                None,
            )

        try:
            # can in principle return data from `.obs`, in which case it's mostly invalid
            return (
                key,
                self.adata.obs_vector(key, layer=self._data_key),
                ColorType.CONT,
                None,
            )
        except KeyError:
            logg.debug(
                f"Key `{key!r}` not found in `adata.obs` or "
                f"`adata{'.raw' if self._use_raw else ''}.var_names`. Ignoring`"
            )

        if same_plot or np.allclose(self.w_all, 1.0):
            return None, "black", ColorType.STR, None

        return "fate probability", np.squeeze(self.w_all), ColorType.CONT, None

    def _create_scaler(self, show_lineage_probability: bool, show_conf_int: bool):
        if not show_lineage_probability:
            return lambda _: _

        minn, maxx = self._return_min_max(show_conf_int)
        return lambda x: (x - minn) / (maxx - minn)

    def _return_min_max(self, show_conf_int: bool):
        if self.y_test is None:
            raise RuntimeError("Run `.predict()` first.")

        vals = [self.y_test]
        # FittedModel Does not need to have these
        if self.y_all is not None:
            vals.append(self.y_all)
        if show_conf_int and self.conf_int is not None:
            vals.append(self.conf_int)

        minn = min(map(np.min, vals))
        maxx = max(map(np.max, vals))

        return minn, maxx


class FailedModel(BaseModel):
    """Model representing a failure of the original :attr:`model`.

    Parameters
    ----------
    model
        The original model which has failed.
    exc
        The exception that caused the :attr:`model` to fail or a :class:`str` containing the message.
        In the latter case, :meth:`~cellrank.models.FailedModel.reraise` a :class:`RuntimeError` with that message.
        If :obj:`None`, :class`UnknownModelError` will eventually be raised.
    """

    # in a functional programming language like Haskell essentially BaseModel would be a Maybe monad and this Nothing
    def __init__(self, model: BaseModel, exc: Optional[Union[BaseException, str]] = None):
        if not isinstance(model, BaseModel):
            raise TypeError(f"Expected `model` to be of type `BaseModel`, found `{type(model).__name__!r}`.")
        if exc is not None:
            if isinstance(exc, str):
                exc = RuntimeError(exc)
            if not isinstance(exc, BaseException):
                raise TypeError(
                    f"Expected `exc` to be either a string or a `BaseException`, found `{type(exc).__name__!r}`."
                )
        else:
            exc = UnknownModelError()

        super().__init__(model.adata, model)

        self._gene = model._gene
        self._lineage = model._lineage
        self._prepared = model._prepared
        self._is_bulk = model._is_bulk

        self._exc = exc

    def prepare(
        self,
        *_args,
        **_kwargs,
    ) -> "FailedModel":
        """Do nothing and return self."""
        return self

    def fit(
        self,
        x: Optional[np.ndarray] = None,
        y: Optional[np.ndarray] = None,
        w: Optional[np.ndarray] = None,
        **kwargs: Any,
    ) -> "FailedModel":
        """Do nothing and return self."""
        return self

    def predict(
        self,
        x_test: Optional[np.ndarray] = None,
        key_added: Optional[str] = "_x_test",
        **kwargs: Any,
    ) -> np.ndarray:
        """Do nothing."""

    def confidence_interval(self, x_test: Optional[np.ndarray] = None, **kwargs) -> np.ndarray:
        """Do nothing."""

    def default_confidence_interval(
        self,
        x_test: Optional[np.ndarray] = None,
        **kwargs: Any,
    ) -> np.ndarray:
        """Do nothing."""

    def reraise(self) -> None:
        """Raise the original exception with additional model information."""
        # retain the exception type and also the original exception
        raise type(self._exc)(f"Fatal model failure `{self}`.") from self._exc

    def _get_colors(
        self,
        key: Optional[str],
        *,
        same_plot: bool = False,
    ) -> Tuple[Optional[str], Optional[Union[str, np.ndarray]], ColorType, Optional[Dict[str, Any]],]:
        return None, "black", ColorType.STR, None

    def _return_min_max(self, show_conf_int: bool):
        return np.inf, -np.inf

    @d.dedent
    def copy(self) -> "FailedModel":
        """%(copy)s"""  # noqa
        return FailedModel(self.model.copy(), exc=self._exc)

    def __bool__(self):
        return False

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}[origin={repr(self.model).strip('<>')}]>"


@d.dedent
class FittedModel(BaseModel):
    """Class representing an already fitted model. Useful when the smoothed gene expression was computed externally.

    Parameters
    ----------
    x_test
        %(base_model_x_test.summary)s
    y_test
        %(base_model_y_test.summary)s
    conf_int
        %(base_model_conf_int.summary)s
        If :obj:`None`, always set ``conf_int=False`` in :meth:`~cellrank.models.BaseModel.plot`.
    x_all
        %(base_model_x_all.summary)s
        If :obj:`None`, always sets ``hide_cells=True`` in :meth:`~cellrank.models.BaseModel.plot`.
    y_all
        %(base_model_y_all.summary)s
        If :obj:`None`, always sets `hide_cells=True`` in :meth:`~cellrank.models.BaseModel.plot`.
    w_all
        %(base_model_w_all.summary)s
        If :obj:`None` and :attr:`x_all` and :attr:`y_all` are present, it will be set an array of :math:`1`.
    """

    def __init__(
        self,
        x_test: ArrayLike,
        y_test: ArrayLike,
        conf_int: Optional[ArrayLike] = None,
        x_all: Optional[ArrayLike] = None,
        y_all: Optional[ArrayLike] = None,
        w_all: Optional[ArrayLike] = None,
    ):
        super().__init__(None, None)

        self._x_test = _densify_squeeze(x_test, self._dtype)[:, np.newaxis]
        self._y_test = _densify_squeeze(y_test, self._dtype)

        if self.x_test.ndim != 2 or self.x_test.shape[1] != 1:
            raise ValueError(f"Expected `x_test` to be of shape `(..., 1)`, found `{self.x_test.shape}`.")

        if self.y_test.shape != (self.x_test.shape[0],):
            raise ValueError(
                f"Expected `y_test` to be of shape `({self.x_test.shape[0]},)`, " f"found `{self.y_test.shape}`."
            )

        if conf_int is not None:
            self._conf_int = _densify_squeeze(conf_int, self._dtype)
            if self.conf_int.shape != (self.x_test.shape[0], 2):
                raise ValueError(f"Expected `conf_int` to be of shape `({self.x_test.shape[0]}, 2)`.")
        else:
            logg.debug("No `conf_int` have been supplied, will be ignored during plotting")

        if x_all is not None and y_all is not None:
            self._x_all = _densify_squeeze(x_all, self._dtype)[:, np.newaxis]
            if self.x_all.ndim != 2 or self.x_all.shape[1] != 1:
                raise ValueError(f"Expected `x_all` to be of shape `(..., 1)`, found `{self.x_all.shape}`.")

            self._y_all = _densify_squeeze(y_all, self._dtype)[:, np.newaxis]
            if self.y_all.shape != self.x_all.shape:
                raise ValueError(
                    f"Expected `y_all` to be of shape `({self.x_all.shape[0]}, 1)`, found `{self.y_all.shape}`."
                )

            if w_all is not None:
                self._w_all = _densify_squeeze(w_all, self._dtype)
                if self.w_all.ndim != 1 or self.w_all.shape[0] != self.x_all.shape[0]:
                    raise ValueError(
                        f"Expected `w_all` to be of shape `({self.x_all.shape[0]},)`, " f"found `{self.w_all.shape}`."
                    )
            else:
                logg.debug("Setting `w_all` to an array of `1`")
                self._w_all = np.ones_like(self.x_all, dtype=np.float64).squeeze(1)
        else:
            logg.debug(
                "None or partially incomplete `x_all` and `y_all` have been supplied, "
                "will be ignored during plotting"
            )

        self._prepared = True

    def _get_colors(
        self,
        key: Optional[str],
        *,
        same_plot: bool = False,
    ) -> Tuple[Optional[str], Optional[Union[str, np.ndarray]], ColorType, Optional[Dict[str, Any]],]:
        # w_all does not need to be defined
        if same_plot or self.w_all is None or np.allclose(self.w_all, 1.0):
            return None, "black", ColorType.STR, None
        return "fate probability", np.squeeze(self.w_all), ColorType.CONT, None

    def prepare(self, *_args, **_kwargs) -> "FittedModel":
        """Do nothing and return self."""
        return self

    def fit(
        self,
        x: Optional[np.ndarray] = None,
        y: Optional[np.ndarray] = None,
        w: Optional[np.ndarray] = None,
        **kwargs: Any,
    ) -> "FittedModel":
        """Do nothing and return self."""
        return self

    @d.dedent
    def predict(
        self,
        x_test: Optional[np.ndarray] = None,
        key_added: Optional[str] = "_x_test",
        **kwargs: Any,
    ) -> np.ndarray:
        """%(base_model_y_test.summary)s."""
        return self._y_test

    @d.dedent
    def confidence_interval(self, x_test: Optional[np.ndarray] = None, **kwargs) -> np.ndarray:
        """%(base_model_conf_int.summary)s Raise a :class:`RuntimeError` if not present."""
        if self.conf_int is None:
            raise RuntimeError(
                "No confidence interval has been supplied. " "Use `conf_int=...` when instantiating this class."
            )
        return self.conf_int

    @d.dedent
    def default_confidence_interval(
        self,
        x_test: Optional[np.ndarray] = None,
        **kwargs: Any,
    ) -> np.ndarray:
        """%(base_model_conf_int.summary)s Raise a :class:`RuntimeError` if not present."""
        return self.confidence_interval()

    @d.dedent
    def copy(self) -> "FittedModel":
        """%(copy)s"""  # noqa
        # here we return a deepcopy since it doesn't make sense to make a shallow one
        return FittedModel.from_model(self)

    @staticmethod
    def from_model(model: BaseModel) -> "FittedModel":
        """Create a :class:`~cellrank.models.FittedModel` from a :class:`~cellrank.models.BaseModel`."""
        if not isinstance(model, BaseModel):
            raise TypeError(f"Expected `model` to be of type `BaseModel`, found `{type(model).__name__!r}`.")

        fm = FittedModel(
            model.x_test,
            model.y_test,
            conf_int=model.conf_int,
            x_all=model.x_all,
            y_all=model.y_all,
            w_all=model.w_all,
        )

        fm._gene = model._gene
        fm._lineage = model._lineage
        fm._is_bulk = model._is_bulk
        fm._prepared = True

        return fm
