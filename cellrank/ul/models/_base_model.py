# -*- coding: utf-8 -*-
"""Base class for all models."""
import re
from abc import ABC, abstractmethod
from copy import copy as _copy
from typing import Any, Tuple, Union, TypeVar, Optional

import numpy as np
from scipy.ndimage import convolve

import matplotlib as mpl
from matplotlib import cm as cm
from matplotlib import colors as mcolors
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from cellrank.tl import Lineage
from cellrank.ul._docs import d
from cellrank.tl._utils import save_fig, _densify_squeeze
from cellrank.ul._utils import _minmax
from cellrank.tl._constants import AbsProbKey

AnnData = TypeVar("AnnData")
_dup_spaces = re.compile(r" +")  # used on repr for underlying model's repr


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
        self._prepared = False

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
    def prepared(self):
        """Whether the model is prepared for fitting."""
        return self._prepared

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
        lineage: Optional[str],
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
            Gene in :paramref:`adata` ``.var_names`` or in :paramref:`adata` ``.raw.var_names``.
        lineage
            Name of a lineage in :paramref:`adata` ``.obsm[lineage_key]``.
            If `None`, all weights will be set to `1`.
        %(backward)s
        %(time_range)s
        data_key
            Key in :paramref:`adata` ``.layers`` or `'X'` for :paramref:`adata` ``.X``.
        time_key
            Key in :paramref:`adata` ``.obs`` where the pseudotime is stored.
        use_raw
            Whether to access :paramref:`adata` ``.raw`` or not.
        threshold
            Consider only cells with weights > ``threshold`` when estimating the test endpoint.
            If `None`, use the median of the weights.
        weight_threshold
            Set all weights below ``weight_threshold`` to either `0` if a :class:`float`,
            or if a :class:`tuple`, to the second value.
        filter_dropouts
            Filter out all cells with expression lower than this.
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
        """
        if use_raw and self.adata.raw is None:
            raise AttributeError("AnnData object has no attribute `.raw`.")

        if data_key not in ["X", "obs"] + list(self.adata.layers.keys()):
            raise KeyError(
                f"Data key must be a key of `adata.layers`: `{list(self.adata.layers.keys())}`, "
                f"`adata.X` or `adata.obs`."
            )
        if time_key not in self.adata.obs:
            raise KeyError(f"Time key `{time_key!r}` not found in `adata.obs`.")

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
            # use `>=` because weights can all be 1
            w_test = w[w >= threshold]
            n_window = n_test_points // 20
            tmp = convolve(w_test, np.ones(n_window) / n_window, mode="nearest")
            val_end = x[w >= threshold][-1 if lineage is None else np.nanargmax(tmp)]

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
        self._prepared = True

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
        lineage_color: str = "black",
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
            Percentile by which to clip the absorption probabilities./
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

        _ = ax.plot(self.x_test, self.y_test, color=lineage_color, lw=lw, label=title)

        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.set_xlabel(xlabel)

        ax.margins(margins)

        if show_conf_int and self.conf_int is not None:
            ax.fill_between(
                self.x_test.squeeze(),
                self.conf_int[:, 0],
                self.conf_int[:, 1],
                alpha=lineage_alpha,
                color=lineage_color,
                linestyle="--",
            )

        if (
            show_cbar
            and not hide_cells
            and not same_plot
            and not np.allclose(self.w_all, 1)
        ):
            norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="2.5%", pad=0.1)
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
            "_prepared",
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
