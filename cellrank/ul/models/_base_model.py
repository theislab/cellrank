"""Base class for all models."""
import re
from abc import ABC, ABCMeta, abstractmethod
from copy import copy as _copy
from copy import deepcopy
from typing import Any, List, Tuple, Union, TypeVar, Callable, Optional

import wrapt

import numpy as np
from scipy.sparse import spmatrix
from scipy.ndimage import convolve

import matplotlib as mpl
from matplotlib import cm
from matplotlib import colors as mcolors
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from cellrank import logging as logg
from cellrank.tl import Lineage
from cellrank.ul._docs import d
from cellrank.tl._utils import save_fig
from cellrank.ul._utils import Pickleable, _minmax, valuedispatch, _densify_squeeze
from cellrank.tl._constants import ModeEnum, AbsProbKey

AnnData = TypeVar("AnnData")
_dup_spaces = re.compile(r" +")  # used on repr for underlying model's repr
ArrayLike = Union[np.ndarray, spmatrix, List, Tuple]


class UnknownModelError(RuntimeError):  # noqa
    pass


class FailedReturnType(ModeEnum):  # noqa
    PREPARE = "prepare"
    FIT = "fit"
    PREDICT = "predict"
    CONFIDENCE_INTERVAL = "confidence_interval"
    DEFAULT_CONFIDENCE_INTERVAL = "default_confidence_interval"
    PLOT = "plot"


def _handle_exception(return_type: FailedReturnType, func: Callable) -> Callable:
    def handle(*, exception_handler: Callable):
        @wrapt.decorator
        def wrapper(wrapped, instance, args, kwargs):
            try:
                if isinstance(instance, FailedModel):
                    instance.reraise()
                return wrapped(*args, **kwargs)
            except Exception as e:  # noqa: B902
                return exception_handler(instance, e, *args, **kwargs)

        return wrapper

    def array_output(
        instance: "BaseModel", exc: BaseException, *_args, **kwargs
    ) -> np.ndarray:
        if not instance._is_bulk:
            raise exc from None

        if kwargs.get("x_test", None) is not None:
            n_obs = kwargs["x_test"].shape[0]
        elif instance.x_test is not None:
            n_obs = instance.x_test.shape[0]
        else:
            n_obs = 200

        return np.full((n_obs, array_shape), np.nan, dtype=instance._dtype)

    def model_output(
        instance: "BaseModel", exc: BaseException, *_args, **_kwargs
    ) -> "FailedModel":
        if not isinstance(instance, FailedModel):
            instance = FailedModel(instance, exc=exc)
        if not instance._is_bulk:
            instance.reraise()

        return instance

    def no_output(instance: "BaseModel", exc: BaseException, *_args, **_kwargs) -> None:
        if not instance._is_bulk:
            raise exc from None
        return None

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


class BaseModelMeta(ABCMeta):
    """Metaclass for all base models."""

    def __new__(cls, clsname, superclasses, attributedict):
        """
        Create a new instance.

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
                fun_name.s,
                _handle_exception(FailedReturnType(fun_name), getattr(obj, fun_name.s)),
            )

        return obj


@d.get_sections(base="base_model", sections=["Parameters"])
@d.dedent
class BaseModel(Pickleable, ABC, metaclass=BaseModelMeta):
    """
    Base class for all model classes.

    Parameters
    ----------
    %(adata)s
    model
        The underlying model that is used for fitting and prediction.
    """

    def __init__(
        self,
        adata: AnnData,
        model: Any,
    ):
        from anndata import AnnData as _AnnData

        if not isinstance(adata, _AnnData) and not isinstance(self, FittedModel):
            # FittedModel doesn't need it
            raise TypeError(
                f"Expected `adata` to be of type `anndata.AnnData`, found `{type(adata).__name__!r}`."
            )

        self._adata = adata
        self._model = model
        self._gene = None
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

        self._dtype = np.float32

    @d.get_summary(base="base_model_prepared")
    @property
    def prepared(self):
        """Whether the model is prepared for fitting."""
        return self._prepared

    @property
    @d.dedent
    def adata(self) -> AnnData:
        """
        Annotated data object.

        Returns
        -------
        %(adata)s
        """
        return self._adata

    @property
    def model(self) -> Any:
        """The underlying model."""  # noqa
        return self._model

    @property
    @d.get_summary(base="base_model_x_all")
    def x_all(self) -> np.ndarray:
        """Unfiltered independent variables of shape `(n_cells, 1)`."""
        return self._x_all

    @property
    @d.get_summary(base="base_model_y_all")
    def y_all(self) -> np.ndarray:
        """Unfiltered dependent variables of shape `(n_cells, 1)`."""
        return self._y_all

    @property
    @d.get_summary(base="base_model_w_all")
    def w_all(self) -> np.ndarray:
        """Unfiltered weights of shape `(n_cells,)`."""
        return self._w_all

    @property
    @d.get_summary(base="base_model_x")
    def x(self) -> np.ndarray:
        """Filtered independent variables of shape `(n_filtered_cells, 1)` used for fitting."""  # noqa
        return self._x

    @property
    @d.get_summary(base="base_model_y")
    def y(self) -> np.ndarray:
        """Filtered dependent variables of shape `(n_filtered_cells, 1)` used for fitting."""  # noqa
        return self._y

    @property
    @d.get_summary(base="base_model_w")
    def w(self) -> np.ndarray:
        """Filtered weights of shape `(n_filtered_cells,)` used for fitting."""  # noqa
        return self._w

    @property
    @d.get_summary(base="base_model_x_test")
    def x_test(self) -> np.ndarray:
        """Independent variables of shape `(n_samples, 1)` used for prediction."""
        return self._x_test

    @property
    @d.get_summary(base="base_model_y_test")
    def y_test(self) -> np.ndarray:
        """Prediction values of shape `(n_samples,)` for :paramref:`x_test`."""
        return self._y_test

    @property
    @d.get_summary(base="base_model_x_hat")
    def x_hat(self) -> np.ndarray:
        """Filtered independent variables used when calculating default confidence interval, usually same as :paramref:`x`."""  # noqa
        return self._x_hat

    @property
    @d.get_summary(base="base_model_y_hat")
    def y_hat(self) -> np.ndarray:
        """Filtered dependent variables used when calculating default confidence interval, usually same as :paramref:`y`."""  # noqa
        return self._y_hat

    @property
    @d.get_summary(base="base_model_conf_int")
    def conf_int(self) -> np.ndarray:
        """Array of shape `(n_samples, 2)` containing the lower and upper bounds of the confidence interval."""
        return self._conf_int

    @d.get_sections(base="base_model_prepare", sections=["Parameters", "Returns"])
    @d.get_full_description(base="base_model_prepare")
    @d.dedent
    def prepare(
        self,
        gene: str,
        lineage: Optional[str],
        backward: bool = False,
        time_range: Optional[Union[float, Tuple[float, float]]] = None,
        data_key: str = "X",
        time_key: str = "latent_time",
        use_raw: bool = False,
        threshold: Optional[float] = None,
        weight_threshold: Union[float, Tuple[float, float]] = (0.01, 0.01),
        filter_cells: Optional[float] = None,
        n_test_points: int = 200,
    ) -> "BaseModel":
        """
        Prepare the model to be ready for fitting.

        Parameters
        ----------
        gene
            Gene in :paramref:`adata` ``.var_names`` or in :paramref:`adata` ``.raw.var_names``.
        lineage
            Name of a lineage in :paramref:`adata` ``.obsm[lineage_key]``. If `None`, all weights will be set to `1`.
        %(backward)s
        %(time_range)s
        data_key
            Key in :paramref:`adata` ``.layers`` or `'X'` for :paramref:`adata` ``.X``.
            If ``use_raw=True``, it's always set to `'X'`.
        time_key
            Key in :paramref:`adata` ``.obs`` where the pseudotime is stored.
        use_raw
            Whether to access :paramref:`adata` ``.raw`` or not.
        threshold
            Consider only cells with weights > ``threshold`` when estimating the test endpoint.
            If `None`, use the median of the weights.
        weight_threshold
            Set all weights below ``weight_threshold`` to ``weight_threshold`` if a :class:`float`,
            or to the second value, if a :class:`tuple`.
        filter_cells
            Filter out all cells with expression values lower than this threshold.
        n_test_points
            Number of test points. If `None`, use the original points based on ``threshold``.

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

                - :paramref:`prepared` - %(base_model_prepared.summary)s
        """

        if use_raw:
            if self.adata.raw is None:
                raise AttributeError("AnnData object has no attribute `.raw`.")
            if data_key != "X":
                data_key = "X"

        if data_key not in ["X", "obs"] + list(self.adata.layers.keys()):
            raise KeyError(
                f"Data key must be a key of `adata.layers`: `{list(self.adata.layers.keys())}`, "
                f"`adata.X` or `adata.obs`."
            )

        lineage_key = str(AbsProbKey.BACKWARD if backward else AbsProbKey.FORWARD)
        if lineage_key not in self.adata.obsm:
            raise KeyError(f"Lineage key `{lineage_key!r}` not found in `adata.obsm`.")
        if not isinstance(self.adata.obsm[lineage_key], Lineage):
            raise TypeError(
                f"Expected `adata.obsm[{lineage_key!r}]` to be of type `cellrank.tl.Lineage`, "
                f"found `{type(self.adata.obsm[lineage_key]).__name__!r}`."
            )

        if lineage is not None:
            _ = self.adata.obsm[lineage_key][lineage]

        if time_key not in self.adata.obs:
            raise KeyError(f"Time key `{time_key!r}` not found in `adata.obs`.")

        if data_key == "obs":
            if gene not in self.adata.obs:
                raise KeyError(f"Unable to find key `{gene!r}` in `adata.obs`.")
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

        if isinstance(weight_threshold, (int, float)):
            weight_threshold = (weight_threshold, weight_threshold)
        if len(weight_threshold) != 2:
            raise ValueError(
                f"Expected `weight_threshold` to be of size `2`, found `{len(weight_threshold)}`."
            )

        if lineage is not None:
            _ = self.adata.obsm[lineage_key][lineage]

        self._obs_names = self.adata.obs_names.values[:]

        x = np.array(self.adata.obs[time_key]).astype(self._dtype)

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
            weight_threshold, val = weight_threshold
            w = _densify_squeeze(self.adata.obsm[lineage_key][lineage].X, self._dtype)
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
        x_test = (
            np.linspace(val_start, val_end, n_test_points)
            if n_test_points is not None
            else x[fil]
        )
        x, y, w = x[fil], y[fil], w[fil]
        self._obs_names = self._obs_names[fil]

        if filter_cells is not None:
            tmp = y.squeeze()
            fil = (tmp >= filter_cells) & (
                ~np.isclose(tmp, filter_cells).astype(np.bool)
            )
            x, y, w = x[fil], y[fil], w[fil]
            self._obs_names = self._obs_names[fil]

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
        self._prepared = True

        return self

    @abstractmethod
    @d.get_sections(base="base_model_fit", sections=["Parameters"])
    @d.get_full_description(base="base_model_fit")
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
            Independent variables, array of shape `(n_samples, 1)`. If `None`, use :paramref:`x`.
        y
            Dependent variables, array of shape `(n_samples, 1)`. If `None`, use :paramref:`y`.
        w
            Optional weights of :paramref:`x`, array of shape `(n_samples,)`. If `None`, use :paramref:`w`.
        kwargs
            Keyword arguments for underlying :paramref:`model`'s fitting function.

        Returns
        -------
        :class:`cellrank.ul.models.BaseModel`
            Fits the model and returns self.
        """
        if not self.prepared:
            raise RuntimeError(
                "The model has not been prepared yet, call `.prepare()` first."
            )

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
    @d.get_sections(base="base_model_predict", sections=["Parameters", "Returns"])
    @d.get_full_description(base="base_model_predict")
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
            Array of shape `(n_samples,)` used for prediction. If `None`, use :paramref:`x_test`.
        key_added
            Attribute name where to save the :paramref:`x_test` for later use. If `None`, don't save it.
        kwargs
            Keyword arguments for underlying :paramref:`model`'s prediction method.

        Returns
        -------
        :class:`numpy.ndarray`
            Updates and returns the following:

                - :paramref:`y_test` - %(base_model_y_test.summary)s
        """

    @abstractmethod
    @d.get_sections(base="base_model_ci", sections=["Parameters", "Returns"])
    @d.get_summary(base="base_model_ci")
    @d.get_full_description(base="base_model_ci")
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
            Array of shape `(n_samples,)` used for confidence interval calculation. If `None`, use :paramref:`x_test`.
        kwargs
            Keyword arguments for underlying :paramref:`model`'s confidence method
            or for :meth:`default_confidence_interval`.

        Returns
        -------
        :class:`numpy.ndarray`
            Updates the following fields:

                - :paramref:`conf_int` - %(base_model_conf_int.summary)s
        """

    @d.dedent
    def default_confidence_interval(
        self,
        x_test: Optional[np.ndarray] = None,
        **kwargs,
    ) -> np.ndarray:
        """
        Calculate the confidence interval, if the underlying :paramref:`model` has no method for it.

        This formula is taken from [DeSalvo70]_, eq. 5.

        Parameters
        ----------
        %(base_model_ci.parameters)s

        Returns
        -------
        %(base_model_ci.returns)s
                - :paramref:`x_hat` - %(base_model_x_hat.summary)s
                - :paramref:`y_hat` - %(base_model_y_hat.summary)s

        References
        ----------
        .. [DeSalvo70] DeSalvo, J. S. (1970),
            *Standard Error of Forecast in Multiple Regression: Proof of a Useful Result.*,
            `RAND Corporation <https://www.rand.org/pubs/papers/P4365.html>`__.
        """

        use_ixs = self.w > 0
        x_hat = self.x[use_ixs]
        if x_test is None:
            x_test = self.x_test

        self._y_hat = self.predict(x_hat, key_added="_x_hat", **kwargs)
        self._y_test = self.predict(x_test, key_added="_x_test", **kwargs)

        n = np.sum(use_ixs)
        sigma_hat = np.sqrt(
            ((self.y_hat - self.y[use_ixs].squeeze()) ** 2).sum() / (n - 2)
        )
        mean = np.mean(self.x)

        stds = sigma_hat * np.sqrt(
            1 + 1 / n + ((self.x_test - mean) ** 2) / ((self.x - mean) ** 2).sum()
        )
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
        abs_prob_cmap: mcolors.ListedColormap = cm.viridis,
        cell_color: str = "black",
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
        dpi: int = None,
        fig: mpl.figure.Figure = None,
        ax: mpl.axes.Axes = None,
        return_fig: bool = False,
        save: Optional[str] = None,
        **kwargs,
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
        lineage_color
            Color for the lineage.
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
        cbar
            Whether to show colorbar.
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
            If :paramref:`self` is :class:`cellrank.ul.models.GAMR`, it can also specify the confidence level,
            the default is `0.95`. Only used when ``show_lineage_probability=True``.
        lineage_probability_color
            Color to use when plotting the smoothed ``lineage_probability``.
            If `None`, it's the same as ``lineage_color``. Only used when ``show_lineage_probability=True``.
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
        kwargs
            Keyword arguments for :meth:`matplotlib.axes.Axes.legend`, e.g. to disable the legend, specify ``loc=None``.
            Only available when ``show_lineage_probability=True``.

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
        hide_cells = (
            hide_cells or self.x_all is None or self.w_all is None or self.y_all is None
        )

        lineage_probability_color = (
            lineage_color
            if lineage_probability_color is None
            else lineage_probability_color
        )

        scaler = kwargs.pop(
            "scaler",
            self._create_scaler(
                lineage_probability,
                show_conf_int=conf_int,
            ),
        )

        if lineage_probability:
            if ylabel in ("expression", self._gene):
                ylabel = f"scaled {ylabel}"

        vmin, vmax = None, None
        if not hide_cells:
            vmin, vmax = _minmax(self.w_all, perc)
            _ = ax.scatter(
                self.x_all.squeeze(),
                scaler(self.y_all.squeeze()),
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
            title = (
                f"{self._gene} @ {self._lineage}"
                if self._lineage is not None
                else f"{self._gene}"
            )

        ax.plot(
            self.x_test, scaler(self.y_test), color=lineage_color, lw=lw, label=title
        )

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

        if (
            lineage_probability
            and not isinstance(self, FittedModel)
            and not np.allclose(self.w, 1.0)
        ):
            from cellrank.pl._utils import _is_any_gam_mgcv

            model = deepcopy(self)
            model._y = self._reshape_and_retype(self.w).copy()
            model = model.fit()

            if not lineage_probability_conf_int:
                y = model.predict()
            elif _is_any_gam_mgcv(model):
                y = model.predict(
                    level=lineage_probability_conf_int
                    if isinstance(lineage_probability_conf_int, float)
                    else 0.95
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

            if kwargs.get("loc", "best") not in (None, "none"):
                ax.legend(handles=handle, **kwargs)

        if (
            cbar
            and not hide_cells
            and not same_plot
            and not np.allclose(self.w_all, 1.0)
        ):
            norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="2%", pad=0.1)
            _ = mpl.colorbar.ColorbarBase(
                cax,
                norm=norm,
                cmap=abs_prob_cmap,
                ticks=np.linspace(norm.vmin, norm.vmax, 5),
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
            Array of shape `(n, 1)` with dtype set to :paramref:`_dtype`.
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
        self, attr_name: Optional[str], arr: Optional[np.ndarray], ndim: int = 2
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
            Expected number of dimensions of the ``arr``.

        Returns
        -------
        :class:`numpy.ndarray`
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
            setattr(dst, attr, _copy(getattr(self, attr)))

    def _shallowcopy_attributes(self, dst: "BaseModel") -> None:
        for attr in ["_gene", "_lineage", "_is_bulk"]:
            setattr(dst, attr, _copy(getattr(self, attr)))
        # user is not exposed to this
        dst._obs_names = self._obs_names
        # always deepcopy the model (we're manipulating it in multiple threads/processed)
        dst._model = deepcopy(self.model)

    @abstractmethod
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
            None
            if self.model is None
            else _dup_spaces.sub(" ", str(self.model).replace("\n", " ")).strip(),
        )

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
    """
    Model representing a failure of the original :paramref:`model`.

    Parameters
    ----------
    model
        The original model which has failed.
    exc
        The exception that caused the :paramref:`model` to fail or a :class:`str` containing the message.
        In the latter case, :meth:`cellrank.ul.models.FailedModel.reraise` a :class:`RuntimeError` with that message.
        If `None`, :class`UnknownModelError` will eventually be raised.
    """

    # in a functional programming language like Haskell essentially BaseModel would be a Maybe monad and this Nothing
    def __init__(
        self, model: BaseModel, exc: Optional[Union[BaseException, str]] = None
    ):
        if not isinstance(model, BaseModel):
            raise TypeError(
                f"Expected `model` to be of type `BaseModel`, found `{type(model).__name__!r}`."
            )
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
        **kwargs,
    ) -> "FailedModel":
        """Do nothing and return self."""
        return self

    def predict(
        self,
        x_test: Optional[np.ndarray] = None,
        key_added: Optional[str] = "_x_test",
        **kwargs,
    ) -> np.ndarray:
        """Do nothing."""

    def confidence_interval(
        self, x_test: Optional[np.ndarray] = None, **kwargs
    ) -> np.ndarray:
        """Do nothing."""

    def default_confidence_interval(
        self,
        x_test: Optional[np.ndarray] = None,
        **kwargs,
    ) -> np.ndarray:
        """Do nothing."""

    def reraise(self) -> None:
        """Raise the original exception with additional model information."""
        # retain the exception type and also the original exception
        raise type(self._exc)(f"Fatal model failure `{self}`.") from self._exc

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
    """
    Class representing an already fitted model. Useful when the smoothed gene expression was computed externally.

    Parameters
    ----------
    x_test
        %(base_model_x_test.summary)s
    y_test
        %(base_model_y_test.summary)s
    conf_int
        %(base_model_conf_int.summary)s
        If `None`, always set ``conf_int=False`` in :meth:`cellrank.ul.models.BaseModel.plot`.
    x_all
        %(base_model_x_all.summary)s
        If `None`, always sets ``hide_cells=True`` in :meth:`cellrank.ul.models.BaseModel.plot`.
    y_all
        %(base_model_y_all.summary)s
        If `None`, always sets `hide_cells=True`` in :meth:`cellrank.ul.models.BaseModel.plot`.
    w_all
        %(base_model_w_all.summary)s
        If `None` and :paramref:`x_all` and :paramref:`y_all` are present, it will be set an array of `1`.
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
            raise ValueError(
                f"Expected `x_test` to be of shape `(..., 1)`, found `{self.x_test.shape}`."
            )

        if self.y_test.shape != (self.x_test.shape[0],):
            raise ValueError(
                f"Expected `y_test` to be of shape `({self.x_test.shape[0]},)`, "
                f"found `{self.y_test.shape}`."
            )

        if conf_int is not None:
            self._conf_int = _densify_squeeze(conf_int, self._dtype)
            if self.conf_int.shape != (self.x_test.shape[0], 2):
                raise ValueError(
                    f"Expected `conf_int` to be of shape `({self.x_test.shape[0]}, 2)`."
                )
        else:
            logg.debug(
                "No `conf_int` have been supplied, will be ignored during plotting"
            )

        if x_all is not None and y_all is not None:
            self._x_all = _densify_squeeze(x_all, self._dtype)[:, np.newaxis]
            if self.x_all.ndim != 2 or self.x_all.shape[1] != 1:
                raise ValueError(
                    f"Expected `x_all` to be of shape `(..., 1)`, found `{self.x_all.shape}`."
                )

            self._y_all = _densify_squeeze(y_all, self._dtype)[:, np.newaxis]
            if self.y_all.shape != self.x_all.shape:
                raise ValueError(
                    f"Expected `y_all` to be of shape `({self.x_all.shape[0]}, 1)`, found `{self.y_all.shape}`."
                )

            if w_all is not None:
                self._w_all = _densify_squeeze(w_all, self._dtype)
                if self.w_all.ndim != 1 or self.w_all.shape[0] != self.x_all.shape[0]:
                    raise ValueError(
                        f"Expected `w_all` to be of shape `({self.x_all.shape[0]},)`, "
                        f"found `{self.w_all.shape}`."
                    )
            else:
                logg.debug("Setting `w_all` to an array of `1`")
                self._w_all = np.ones_like(self.x_all).squeeze(1)
        else:
            logg.debug(
                "None or partially incomplete `x_all` and `y_all` have been supplied, "
                "will be ignored during plotting"
            )

        self._prepared = True

    def prepare(self, *_args, **_kwargs) -> "FittedModel":
        """Do nothing and return self."""
        return self

    def fit(
        self,
        x: Optional[np.ndarray] = None,
        y: Optional[np.ndarray] = None,
        w: Optional[np.ndarray] = None,
        **kwargs,
    ) -> "FittedModel":
        """Do nothing and return self."""
        return self

    @d.dedent
    def predict(
        self,
        x_test: Optional[np.ndarray] = None,
        key_added: Optional[str] = "_x_test",
        **kwargs,
    ) -> np.ndarray:
        """%(base_model_y_test.summary)s"""  # noqa
        return self._y_test

    @d.dedent
    def confidence_interval(
        self, x_test: Optional[np.ndarray] = None, **kwargs
    ) -> np.ndarray:
        """%(base_model_conf_int.summary)s Raise a :class:`RuntimeError` if not present."""
        if self.conf_int is None:
            raise RuntimeError(
                "No confidence interval has been supplied. "
                "Use `conf_int=...` when instantiating this class."
            )
        return self.conf_int

    @d.dedent
    def default_confidence_interval(
        self,
        x_test: Optional[np.ndarray] = None,
        **kwargs,
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
        """Create a :class:`cellrank.ul.models.FittedModel` instance from :class:`cellrank.ul.models.BaseModel`."""
        if not isinstance(model, BaseModel):
            raise TypeError(
                f"Expected `model` to be of type `BaseModel`, found `{type(model).__name__!r}`."
            )

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
