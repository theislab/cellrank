import collections
import copy
import itertools
import pathlib
from typing import (
    Any,
    Callable,
    Dict,
    List,
    Mapping,
    NamedTuple,
    Optional,
    Sequence,
    Tuple,
    TypeVar,
    Union,
)

import numpy as np
import pandas as pd
from pandas.api.types import infer_dtype, is_categorical_dtype, is_numeric_dtype

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

from anndata import AnnData

from cellrank import logging as logg
from cellrank._utils._colors import _create_categorical_colors
from cellrank._utils._docs import d
from cellrank._utils._enum import DEFAULT_BACKEND
from cellrank._utils._parallelize import parallelize
from cellrank._utils._utils import _unique_order_preserving, save_fig
from cellrank.models import GAMR, BaseModel, FailedModel, SKLearnModel
from cellrank.models._base_model import ColorType

__all__ = ["composition"]

Queue = TypeVar("Queue")
Graph = TypeVar("Graph")


_ERROR_INCOMPLETE_SPEC = "No options were specified for {}. " "Consider specifying a fallback model using '*'."
_time_range_type = Optional[Union[float, Tuple[Optional[float], Optional[float]]]]
_return_model_type = Mapping[str, Mapping[str, BaseModel]]
_input_model_type = Union[BaseModel, _return_model_type]
_callback_type = Optional[Union[Callable, Mapping[str, Mapping[str, Callable]]]]


class BulkRes(NamedTuple):
    x_test: np.ndarray
    y_test: np.ndarray


def _is_any_gam_mgcv(models: Union[BaseModel, Dict[str, Dict[str, BaseModel]]]) -> bool:
    """Return whether any models to be fit are from R's `mgcv` package.

    Parameters
    ----------
    models
        Model used for fitting.

    Returns
    -------
    :obj:`True` if any of the models is from R's mgcv package, else :obj:`False`.
    """
    return isinstance(models, GAMR) or (
        isinstance(models, dict) and any(isinstance(m, GAMR) for ms in models.values() for m in ms.values())
    )


def _create_models(
    model: _input_model_type, obs: Sequence[str], lineages: Sequence[Optional[str]]
) -> _return_model_type:
    """Create models for each gene and lineage.

    Parameters
    ----------
    obs
        Sequence of observations, such as genes.
    lineages
        Sequence of genes.

    Returns
    -------
    The created models.
    """

    def process_lineages(obs_name: str, lin_names: Union[BaseModel, Dict[Optional[str], Any]]):
        if isinstance(lin_names, BaseModel):
            # sharing the same models for all lineages
            for lin_name in lineages:
                models[obs_name][lin_name] = copy.copy(lin_names)
            return
        if not isinstance(lin_names, dict):
            raise TypeError(
                f"Expected the model to be either a lineage specific `dict` or a `BaseModel`, "
                f"found `{type(lin_names).__name__!r}`."
            )

        lin_rest_model = lin_names.get("*", None)  # do not pop
        if lin_rest_model is not None and not isinstance(lin_rest_model, BaseModel):
            raise TypeError(
                f"Expected the lineage fallback model for gene `{obs_name!r}` to be of type `BaseModel`, "
                f"found `{type(lin_rest_model).__name__!r}`."
            )

        for lin_name, mod in lin_names.items():
            if lin_name == "*":
                continue
            if not isinstance(mod, BaseModel):
                raise TypeError(
                    f"Expected the model for gene `{obs_name!r}` and lineage `{lin_name!r}` "
                    f"to be of type `BaseModel`, found `{type(mod).__name__!r}`."
                )
            models[obs_name][lin_name] = copy.copy(mod)

        if set(models[obs_name].keys()) & lineages == lineages:
            return

        if lin_rest_model is not None:
            for lin_name in lineages - set(models[obs_name].keys()):
                models[obs_name][lin_name] = copy.copy(lin_rest_model)
        else:
            raise ValueError(_ERROR_INCOMPLETE_SPEC.format(f"all lineages for gene `{obs_name!r}`"))

    if not len(lineages):
        raise ValueError("No lineages have been selected.")

    if not len(obs):
        raise ValueError("No genes have been selected.")

    if isinstance(model, BaseModel):
        return {
            o: {lin: copy.copy(model) for lin in _unique_order_preserving(lineages)}
            for o in _unique_order_preserving(obs)
        }

    lineages, obs = (
        set(_unique_order_preserving(lineages)),
        set(_unique_order_preserving(obs)),
    )
    models = collections.defaultdict(dict)

    if isinstance(model, dict):
        obs_rest_model = model.pop("*", None)
        if obs_rest_model is not None and not isinstance(obs_rest_model, BaseModel):
            raise TypeError(
                f"Expected the gene fallback model to be of type `BaseModel`, "
                f"found `{type(obs_rest_model).__name__!r}`."
            )

        for obs_name, lin_names in model.items():
            process_lineages(obs_name, lin_names)

        if obs_rest_model is not None:
            for obs_name in obs - set(model.keys()):
                process_lineages(obs_name, model.get(obs_name, obs_rest_model))
        elif set(model.keys()) != obs:
            raise ValueError(_ERROR_INCOMPLETE_SPEC.format(f"genes `{list(obs - set(model.keys()))}`."))
    else:
        raise TypeError(
            f"Class `{type(model).__name__!r}` must be of type `BaseModel` or "
            f"a gene and lineage specific `dict` of `BaseModel`.."
        )

    if set(models.keys()) & obs != obs:
        raise ValueError(f"Missing gene models for the following genes: `{list(obs - set(models.keys()))}`.")

    for gene, vs in models.items():
        if set(vs.keys()) & lineages != lineages:
            raise ValueError(f"Missing lineage models for the gene `{gene!r}`: `{list(lineages - set(vs.keys()))}`.")

    return models


def _fit_bulk_helper(
    genes: Sequence[str],
    models: _input_model_type,
    callbacks: _callback_type,
    lineages: Sequence[Optional[str]],
    time_range: Sequence[Union[float, Tuple[float, float]]],
    return_models: bool = False,
    queue: Optional[Queue] = None,
    **kwargs,
) -> Dict[str, Dict[str, BaseModel]]:
    """Fit model for given genes and lineages.

    Parameters
    ----------
    genes
        Genes for which to fit the models.
    models
        Gene and lineage specific models.
    callbacks
        Gene and lineage specific prepare callbacks.
    lineages
        Lineages for which to fit the models.
    time_range
        Minimum and maximum pseudotimes.
    return_models
        Whether to return the full models or just tuple ``(x_test,  y_test)``.
    queue
        Signalling queue in the parent process/thread used to update the progress bar.
    kwargs
        Keyword arguments for :func:`~cellrank.models.BaseModel.prepare`.

    Returns
    -------
    The fitted models, optionally containing the confidence interval in the form of
    ``{'gene1': {'lineage1': <model11>, ...}, ...}``.
    If any step has failed, the model will be of type :class:`~cellrank.models.FailedModel`.
    """
    if len(lineages) != len(time_range):
        raise ValueError(
            f"Expected `lineage` and `time_range` to be of same length, "
            f"found `{len(lineages)}` != `{len(time_range)}`."
        )

    conf_int = return_models and kwargs.pop("conf_int", False)
    res = {}

    for gene in genes:
        res[gene] = {}
        for ln, tr in zip(lineages, time_range):
            cb = callbacks[gene][ln]
            model = models[gene][ln]
            model._is_bulk = True

            model = cb(model, gene=gene, lineage=ln, time_range=tr, **kwargs)
            model = model.fit()
            # GAMR is a bit faster if we don't need the conf int
            # if it's needed, `.predict` will calculate it and `confidence_interval` will do nothing
            if not conf_int:
                model.predict()
            elif _is_any_gam_mgcv(model):
                model.predict(level=conf_int if isinstance(conf_int, float) else 0.95)
            else:
                model.predict()
                model.confidence_interval()

            res[gene][ln] = model if return_models else BulkRes(model.x_test, model.y_test)

        if queue is not None:
            queue.put(1)

    if queue is not None:
        queue.put(None)

    return res


def _fit_bulk(
    models: Mapping[str, Mapping[str, Callable]],
    callbacks: Mapping[str, Mapping[str, Callable]],
    genes: Union[str, Sequence[str]],
    lineages: Union[str, Sequence[str]],
    time_range: _time_range_type,
    parallel_kwargs: dict,
    return_models: bool = False,
    filter_all_failed: bool = True,
    **kwargs,
) -> Tuple[_return_model_type, _return_model_type, Sequence[str], Sequence[str]]:
    """Fit models for given genes and lineages.

    Parameters
    ----------
    models
        Gene and lineage specific estimators.
    callbacks
        Functions which are called to prepare the ``models``.
    genes
        Genes for which to fit the ``models``.
    lineages
        Lineages for which to fit the ``models``.
    time_range
        Possibly ``lineages`` specific start- and endtimes.
    parallel_kwargs
        Keyword arguments for :func:`~cellrank._utils._parallelize.parallelize`.
    return_models
        Whether to return the full models or just a dictionary of dictionaries of :class:`~typing.NamedTuple`,
        ``(x_test, y_test)``. This is highly discouraged because no meaningful error messages will be produced.
    filter_all_failed
        Whether to filter out all models which have failed.

    Returns
    -------
    All the models, including the failed ones. It is a nested dictionary where keys are the ``genes`` and the values
    is again a :class:`dict`, where keys are ``lineages`` and values are the failed or fitted models or
    the :class:`collections.namedtuple`, based on ``return_models = True``.

    Same as above, but can contain failed models if ``filter_all_failed=False``. In that case, it is guaranteed
    that this dictionary will contain only genes which have been successfully fitted for at least 1 lineage.
    If ``return_models = True``, the models are just a :class:`collections.namedtuple` of `(x_test, y_test)`.

    All the genes of the filtered models.

    All the lineage of the filtered models.
    """
    if isinstance(genes, str):
        genes = [genes]
    if isinstance(lineages, str):
        lineages = [lineages]

    if isinstance(time_range, (tuple, float, int, type(None))):
        time_range = [time_range] * len(lineages)
    elif len(time_range) != len(lineages):
        raise ValueError(f"Expected time ranges to be of length `{len(lineages)}`, found `{len(time_range)}`.")

    n_jobs = parallel_kwargs.pop("n_jobs", 1)

    start = logg.info(f"Computing trends using `{n_jobs}` core(s)")
    models = parallelize(
        _fit_bulk_helper,
        genes,
        unit="gene" if kwargs.get("data_key", "gene") != "obs" else "obs",
        n_jobs=n_jobs,
        extractor=lambda modelss: {k: v for m in modelss for k, v in m.items()},
    )(
        models=models,
        callbacks=callbacks,
        lineages=lineages,
        time_range=time_range,
        return_models=return_models,
        **kwargs,
    )
    logg.info("    Finish", time=start)

    return _filter_models(models, return_models=return_models, filter_all_failed=filter_all_failed)


def _filter_models(
    models, return_models: bool = False, filter_all_failed: bool = True
) -> Tuple[_return_model_type, _return_model_type, Sequence[str], Sequence[str]]:
    def is_valid(x: Union[BaseModel, BulkRes]) -> bool:
        if return_models:
            assert isinstance(x, BaseModel), f"Expected `BaseModel`, found `{type(x).__name__!r}`."
            return bool(x)

        return x.x_test is not None and x.y_test is not None and np.all(np.isfinite(x.y_test))

    modelmat = pd.DataFrame(models).T
    modelmask = modelmat.applymap(is_valid)
    to_keep = modelmask[modelmask.any(axis=1)]
    to_keep = to_keep.loc[:, to_keep.any(axis=0)].T

    filtered_models = {
        gene: {
            ln: models[gene][ln]
            for ln in (ln for ln in v if (is_valid(models[gene][ln]) if filter_all_failed else True))
        }
        for gene, v in to_keep.to_dict().items()
    }

    if not len(filtered_models):
        if not return_models:
            raise RuntimeError(
                "Fitting has failed for all gene/lineage combinations. "
                "Specify `return_models=True` for more information."
            )
        for ms in models.values():
            for model in ms.values():
                assert isinstance(model, FailedModel), f"Expected `FailedModel`, found `{type(model).__name__!r}`."
                model.reraise()

    if not np.all(modelmask.values):
        failed_models = modelmat.values[~modelmask.values]
        logg.warning(
            f"Unable to fit `{len(failed_models)}` models." + ""
            if return_models
            else "Consider specify `return_models=True` for further inspection."
        )
        logg.debug("The failed models were:\n`{}`".format("\n".join(f"    {m}" for m in failed_models)))

    # lineages is the max number of lineages
    return models, filtered_models, tuple(filtered_models.keys()), tuple(to_keep.index)


@d.dedent
def _trends_helper(
    models: Dict[str, Dict[str, Any]],
    gene: str,
    transpose: bool = False,
    lineage_names: Optional[Sequence[str]] = None,
    same_plot: bool = False,
    sharey: Union[str, bool] = False,
    show_ylabel: bool = True,
    show_lineage: Union[bool, np.ndarray] = True,
    show_xticks_and_label: Union[bool, np.ndarray] = True,
    lineage_cmap: Optional[Union[mpl.colors.ListedColormap, Sequence]] = None,
    lineage_probability_color: Optional[str] = None,
    fate_prob_cmap=cm.viridis,
    gene_as_title: bool = False,
    cell_color: Optional[str] = None,
    legend_loc: Optional[str] = "best",
    fig: mpl.figure.Figure = None,
    axes: Union[mpl.axes.Axes, Sequence[mpl.axes.Axes]] = None,
    **kwargs: Any,
) -> None:
    """Plot an expression gene for some lineages.

    Parameters
    ----------
    %(adata)s
    %(model)s
    gene
        Name of the gene in :attr:`~anndata.AnnData.var_names`.
    ln_key
        Key in :attr:`~anndata.AnnData.bosm` where to find the lineages.
    lineage_names
        Names of lineages to plot.
    same_plot
        Whether to plot all lineages in the same plot or separately.
    sharey
        Whether the y-axis is being shared.
    show_ylabel
        Whether to show y-label on the y-axis. Usually, only the first column will contain the y-label.
    show_lineage
        Whether to show the lineage as the title. Usually, only first row will contain the lineage names.
    show_xticks_and_label
        Whether to show x-ticks and x-label. Usually, only the last row will show this.
    lineage_cmap
        Colormap to use when coloring the lineage. When ``transpose = True``, this corresponds to the color of genes.
    lineage_probability_color
        Actual color of one ``lineage``. Only used when ``same_plot = True`` and ``transpose = True`` and
        ``lineage_probability = True``.
    fate_prob_cmap:
        Colormap to use when coloring in the fate probabilities, if they are being plotted.
    gene_as_title
        Whether to use the gene names as titles (with lineage names as well) or on the y-axis.
    legend_loc
        Location of the legend. If :obj:`None`, don't show any legend.
    fig
        Figure to use.
    ax
        Ax to use.
    kwargs
        Keyword arguments for :meth:`~cellrank.models.BaseModel.plot`.

    Returns
    -------
    %(just_plots)s
    """
    n_lineages = len(lineage_names)
    if same_plot:
        axes = [axes] * len(lineage_names)

    axes = np.ravel(axes)

    percs = kwargs.pop("perc", None)
    if percs is None or not isinstance(percs[0], (tuple, list)):
        percs = [percs]

    same_perc = False  # we need to show colorbar always if percs differ
    if len(percs) != n_lineages or n_lineages == 1:
        if len(percs) != 1:
            raise ValueError(f"Percentile must be a collection of size `1` or `{n_lineages}`, got `{len(percs)}`.")
        same_perc = True
        percs = percs * n_lineages

    hide_cells = kwargs.pop("hide_cells", False)
    show_cbar = kwargs.pop("cbar", True)
    show_prob = kwargs.pop("lineage_probability", False)

    if same_plot:
        if not transpose:
            lineage_colors = (
                lineage_cmap.colors if lineage_cmap is not None and hasattr(lineage_cmap, "colors") else lineage_cmap
            )
        else:
            # this should be fine w.r.t. to the missing genes, since they are in the same order AND
            # we're also passing the failed models (this is important)
            # these are actually gene colors
            if lineage_cmap is not None:
                lineage_colors = (
                    lineage_cmap.colors
                    if hasattr(lineage_cmap, "colors")
                    else [c for _, c in zip(lineage_names, lineage_cmap)]
                )
            else:
                lineage_colors = _create_categorical_colors(n_lineages)
    else:
        lineage_colors = (("black" if not colors.is_color_like(lineage_cmap) else lineage_cmap),) * n_lineages

    if n_lineages > len(lineage_colors):
        raise ValueError(f"Expected at least `{n_lineages}` colors, found `{len(lineage_colors)}`.")

    lineage_color_mapper = {ln: lineage_colors[i] for i, ln in enumerate(lineage_names)}
    successful_models = {ln: models[gene][ln] for ln in lineage_names if models[gene][ln]}

    if show_prob and same_plot:
        minns, maxxs = zip(
            *(
                models[gene][n]._return_min_max(
                    show_conf_int=kwargs.get("conf_int", False),
                )
                for n in lineage_names
            )
        )
        minn, maxx = min(minns), max(maxxs)
        kwargs["loc"] = legend_loc
        kwargs["scaler"] = lambda x: (x - minn) / (maxx - minn)
    else:
        kwargs["loc"] = None

    if isinstance(show_xticks_and_label, bool):
        show_xticks_and_label = [show_xticks_and_label] * len(lineage_names)
    elif len(show_xticks_and_label) != len(lineage_names):
        raise ValueError(
            f"Expected `show_xticks_label` to be the same length as `lineage_names`, "
            f"found `{len(show_xticks_and_label)}` != `{len(lineage_names)}`."
        )

    if isinstance(show_lineage, bool):
        show_lineage = [show_lineage] * len(lineage_names)
    elif len(show_lineage) != len(lineage_names):
        raise ValueError(
            f"Expected `show_lineage` to be the same length as `lineage_names`, "
            f"found `{len(show_lineage)}` != `{len(lineage_names)}`."
        )

    last_ax = None
    ylabel_shown = False
    cells_shown = False
    obs_legend_loc = kwargs.pop("obs_legend_loc", "best")

    for i, (name, ax, perc) in enumerate(zip(lineage_names, axes, percs)):
        model = models[gene][name]
        if isinstance(model, FailedModel):
            if not same_plot:
                ax.remove()
            continue

        if same_plot:
            if gene_as_title:
                title = gene
                ylabel = "expression" if show_ylabel else None
            else:
                title = ""
                ylabel = gene
        else:
            if gene_as_title:
                title = None
                ylabel = "expression" if not ylabel_shown else None
            else:
                title = (name if name is not None else "no lineage") if show_lineage[i] else ""
                ylabel = gene if not ylabel_shown else None

        model.plot(
            ax=ax,
            fig=fig,
            perc=perc,
            cell_color=cell_color,
            cbar=False,
            obs_legend_loc=None,
            title=title,
            hide_cells=True if hide_cells else cells_shown if same_plot else False,
            same_plot=same_plot,
            lineage_color=lineage_color_mapper[name],
            lineage_probability_color=lineage_probability_color,
            fate_prob_cmap=fate_prob_cmap,
            lineage_probability=show_prob,
            ylabel=ylabel,
            **kwargs,
        )
        if sharey in ("row", "all", True) and not ylabel_shown:
            plt.setp(ax.get_yticklabels(), visible=True)

        if show_xticks_and_label[i]:
            plt.setp(ax.get_xticklabels(), visible=True)
        else:
            ax.set_xlabel(None)

        last_ax = ax
        ylabel_shown = True
        cells_shown = True

    key, color, typp, mapper = model._get_colors(cell_color, same_plot=same_plot)
    if typp == ColorType.CAT:
        if not hide_cells:
            model._maybe_add_legend(fig, ax, mapper=mapper, title=key, loc=obs_legend_loc, is_line=False)
    elif typp == ColorType.CONT and same_perc and show_cbar and not hide_cells:
        if isinstance(color, np.ndarray):
            # plotting cont. observation other than lin. probs as a color
            vmin = np.min(color)
            vmax = np.max(color)
        else:
            vmin = np.min([model.w_all for model in successful_models.values()])
            vmax = np.max([model.w_all for model in successful_models.values()])
        norm = colors.Normalize(vmin=vmin, vmax=vmax)

        for ax in axes:
            children = [c for c in ax.get_children() if isinstance(c, mpl.collections.PathCollection)]
            if len(children):
                children[0].set_norm(norm)

        divider = make_axes_locatable(last_ax)
        cax = divider.append_axes("right", size="2%", pad=0.1)
        _ = mpl.colorbar.ColorbarBase(
            cax,
            norm=norm,
            cmap=fate_prob_cmap,
            label=key,
            ticks=np.linspace(norm.vmin, norm.vmax, 5),
        )

    if same_plot and lineage_names != [None]:
        model._maybe_add_legend(
            fig,
            ax,
            mapper={ln: lineage_color_mapper[ln] for ln in successful_models},
            loc=legend_loc,
        )


def _position_legend(ax: mpl.axes.Axes, legend_loc: str, **kwargs) -> mpl.legend.Legend:
    """Position legend inside or outside the figure.

    Parameters
    ----------
    ax
        Ax where to position the legend.
    legend_loc
        Position of legend.
    kwargs
        Keyword arguments for :meth:`matplotlib.axes.Axes.legend`.

    Returns
    -------
    The created legend.
    """
    if legend_loc == "center center out":
        raise ValueError("Invalid option: `'center center out'`.")
    if legend_loc == "best":
        return ax.legend(loc="best", **kwargs)

    tmp, loc = legend_loc.split(" "), ""

    if len(tmp) == 1:
        height, rest = tmp[0], []
        width = "right" if height in ("upper", "top", "center") else "left"
    else:
        height, width, *rest = legend_loc.split(" ")
        if rest:
            if len(rest) != 1:
                raise ValueError(f"Expected only 1 additional modifier ('in' or 'out'), found `{list(rest)}`.")
            if rest[0] not in ("in", "out"):
                raise ValueError(f"Invalid modifier `{rest[0]!r}`. Valid options are: `'in', 'out'`.")
            if rest[0] == "in":  # ignore in, it's default
                rest = []

    if height in ("upper", "top"):
        y = 1.55 if width == "center" else 1.025
        loc += "upper"
    elif height == "center":
        y = 0.5
        loc += "center"
    elif height in ("lower", "bottom"):
        y = -0.55 if width == "center" else -0.025
        loc += "lower"
    else:
        raise ValueError(
            f"Invalid legend position on y-axis: `{height!r}`. "
            f"Valid options are: `'upper', 'top', 'center', 'lower', 'bottom'`."
        )

    if width == "left":
        x = -0.05
        loc += " right" if rest else " left"
    elif width == "center":
        x = 0.5
        if height != "center":  # causes to be like top center
            loc += " center"
    elif width == "right":
        x = 1.05
        loc += " left" if rest else " right"
    else:
        raise ValueError(
            f"Invalid legend position on x-axis: `{width!r}`. " f"Valid options are: `'left', 'center', 'right'`."
        )

    if rest:
        kwargs["bbox_to_anchor"] = (x, y)

    return ax.legend(loc=loc, **kwargs)


def _get_backend(model, backend: str) -> str:
    return DEFAULT_BACKEND if _is_any_gam_mgcv(model) else backend


@d.dedent
def _create_callbacks(
    adata: AnnData,
    callback: Optional[Callable],
    obs: Sequence[str],
    lineages: Sequence[Optional[str]],
    perform_sanity_check: Optional[bool] = None,
    **kwargs,
) -> Dict[str, Dict[str, Callable]]:
    """Create models for each gene and lineage.

    Parameters
    ----------
    %(adata)s
    callback
        Gene and lineage specific prepare callbacks.
    obs
        Sequence of observations, such as genes.
    lineages
        Sequence of genes.
    perform_sanity_check
        Whether to check if all callbacks have the correct signature. This is done by instantiating
        dummy model and running the function. We're assuming that the callback isn't really a pricey operation.
        If :obj:`None`, it is only performed for non-default callbacks.
    kwargs
        Keyword arguments for ``callback`` when performing the sanity check.

    Returns
    -------
    The created callbacks.
    """

    def process_lineages(obs_name: str, lin_names: Optional[Union[Callable, Dict[Optional[str], Any]]]) -> None:
        if lin_names is None:
            lin_names = _default_model_callback

        if callable(lin_names):
            # sharing the same models for all lineages
            for lin_name in lineages:
                callbacks[obs_name][lin_name] = lin_names
            return
        if not isinstance(lin_names, dict):
            raise TypeError(
                f"Expected the lineage callback to be either `callable` or a dictionary of callables, "
                f"found `{type(lin_names).__name__!r}`."
            )

        lin_rest_callback = lin_names.get("*", _default_model_callback) or _default_model_callback  # do not pop
        if not callable(lin_rest_callback):
            raise TypeError(
                f"Expected the lineage fallback callback for gene `{obs_name!r}` to be `callable`, "
                f"found `{type(lin_rest_callback).__name__!r}`."
            )

        for lin_name, cb in lin_names.items():
            if lin_name == "*":
                continue
            if not callable(cb):
                raise TypeError(
                    f"Expected the callback for gene `{obs_name!r}` and lineage `{lin_name!r}` "
                    f"to be `callable`, found `{type(cb).__name__!r}`."
                )
            callbacks[obs_name][lin_name] = cb

        for lin_name in lineages - set(callbacks[obs_name].keys()):
            callbacks[obs_name][lin_name] = lin_rest_callback

    def maybe_sanity_check(callbacks: Dict[str, Dict[str, Callable]]) -> None:
        if not perform_sanity_check:
            return

        from sklearn.svm import SVR

        logg.debug("Performing callback sanity checks")
        for gene in callbacks:
            for lineage, cb in callbacks[gene].items():
                # create the model here because the callback can search the attribute
                dummy_model = SKLearnModel(adata, model=SVR())
                try:
                    model = cb(dummy_model, gene=gene, lineage=lineage, **kwargs)
                    assert model is dummy_model, (
                        "Creation of new models is not allowed. " "Ensure that callback returns the same model."
                    )
                    assert model.prepared, "Model is not prepared. Ensure that callback calls `.prepare()`."
                    assert model._gene == gene, f"Callback modified the gene from `{gene!r}` to `{model._gene!r}`."
                    assert (
                        model._lineage == lineage
                    ), f"Callback modified the lineage from `{lineage!r}` to `{model._lineage!r}`."
                    if isinstance(model, FailedModel):
                        model.reraise()
                except Exception as e:  # noqa: BLE001
                    raise RuntimeError(
                        f"Callback validation failed for gene `{gene!r}` and lineage `{lineage!r}`."
                    ) from e

    def all_callbacks_are_default(cbs: dict) -> bool:
        # this correctly implicitly handles '*': None
        for vs in cbs.values():
            if isinstance(vs, dict):
                for cb in vs.values():
                    if callable(cb) and cb is not _default_model_callback:
                        return False
            elif callable(vs) and vs is not _default_model_callback:
                return False

        return True

    if not len(lineages):
        raise ValueError("No lineages have been selected.")

    if not len(obs):
        raise ValueError("No genes have been selected.")

    if callback is None:
        callback = _default_model_callback

    if perform_sanity_check is None:
        perform_sanity_check = (
            not all_callbacks_are_default(callback)
            if isinstance(callback, dict)
            else callback is not _default_model_callback
        )

    if callable(callback):
        callbacks = {o: {lin: callback for lin in lineages} for o in obs}
        maybe_sanity_check(callbacks)
        return callbacks

    lineages, obs = (
        set(_unique_order_preserving(lineages)),
        set(_unique_order_preserving(obs)),
    )
    callbacks = collections.defaultdict(dict)

    if isinstance(callback, dict):
        # can be specified as None
        obs_rest_callback = callback.pop("*", _default_model_callback) or _default_model_callback

        for obs_name, lin_names in callback.items():
            process_lineages(obs_name, lin_names)

        if callable(obs_rest_callback):
            for obs_name in obs - set(callback.keys()):
                process_lineages(obs_name, callback.get(obs_name, obs_rest_callback))
        else:
            raise TypeError(
                f"Expected the gene fallback callback to be `callable`, "
                f"found `{type(obs_rest_callback).__name__!r}`."
            )
    else:
        raise TypeError(
            f"Class `{type(callback).__name__!r}` must be `callable` or "
            f"a gene and lineage specific `dict` of `callables`."
        )

    if set(callbacks.keys()) & obs != obs:
        raise ValueError(f"Missing gene callbacks for the following genes: `{list(obs - set(callbacks.keys()))}`.")

    for gene, vs in callbacks.items():
        if set(vs.keys()) & lineages != lineages:
            raise ValueError(f"Missing lineage callbacks for gene `{gene!r}`: `{list(lineages - set(vs.keys()))}`.")

    maybe_sanity_check(callbacks)

    return callbacks


def _default_model_callback(model: BaseModel, **kwargs) -> BaseModel:
    # we could filter kwargs, but it's better not to - this will detect if we pass useless stuff
    return model.prepare(**kwargs)


@d.dedent
def composition(
    adata: AnnData,
    key: str,
    fontsize: Optional[str] = None,
    figsize: Optional[Tuple[float, float]] = None,
    dpi: Optional[float] = None,
    save: Optional[Union[str, pathlib.Path]] = None,
    **kwargs: Any,
) -> None:
    """Plot a pie chart for categorical annotation.

    Parameters
    ----------
    %(adata)s
    key
        Key in :attr:`~anndata.AnnData.obs` containing categorical observation.
    fontsize
        Font size for the pie chart labels.
    %(plotting)s
    kwargs
        Keyword arguments for :meth:`~matplotlib.axes.Axes.pie`.

    Returns
    -------
    %(just_plots)s
    """
    if key not in adata.obs:
        raise KeyError(f"Data not found in `adata.obs[{key!r}]`.")
    if not is_categorical_dtype(adata.obs[key]):
        raise TypeError(f"Expected `adata.obs[{key!r}]` is not `categorical`, found `{infer_dtype(adata.obs[key])}`.")

    colors = adata.uns.get(f"{key}_colors", None)
    x = adata.obs[key].value_counts()

    # plot these fractions in a pie plot
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    ax.pie(
        x=x,
        labels=x.index,
        colors=colors,
        textprops={"fontsize": fontsize},
        **kwargs,
    )
    ax.set_title(f"composition by {key}")

    if save is not None:
        save_fig(fig, save)


# modified from: https://github.com/CarlEkerot/held-karp
def _held_karp(dists: np.ndarray) -> Tuple[float, np.ndarray]:
    """Held-Karp algorithm solves the Traveling Salesman Problem.

    This algorithm uses dynamic programming with memoization.

    Parameters
    ----------
    dists
        Distance matrix.

    Returns
    -------
    The cost and the path.
    """
    n = len(dists)

    # Maps each subset of the nodes to the cost to reach that subset, as well
    # as what node it passed before reaching this subset.
    # Node subsets are represented as set bits.
    C = {}

    # Set transition cost from initial state
    for k in range(1, n):
        C[1 << k, k] = (dists[0][k], 0)

    # Iterate subsets of increasing length and store intermediate results
    # in classic dynamic programming manner
    for subset_size in range(2, n):
        for subset in itertools.combinations(range(1, n), subset_size):
            # Set bits for all nodes in this subset
            bits = 0
            for bit in subset:
                bits |= 1 << bit

            # Find the lowest cost to get to this subset
            for k in subset:
                prev = bits & ~(1 << k)

                res = []
                for m in subset:
                    if m == 0 or m == k:
                        continue
                    res.append((C[prev, m][0] + dists[m][k], m))
                C[bits, k] = min(res)

    # We're interested in all bits but the least significant (the start state)
    bits = (2**n - 1) - 1

    # Calculate optimal cost
    res = []
    for k in range(1, n):
        res.append((C[bits, k][0] + dists[k][0], k))
    opt, parent = min(res)

    # Backtrack to find full path
    path = []
    for _ in range(n - 1):
        path.append(parent)
        new_bits = bits & ~(1 << parent)
        _, parent = C[bits, parent]
        bits = new_bits

    # Add implicit start state
    path.append(0)

    return opt, np.array(path)[::-1]


def _get_categorical_colors(adata: AnnData, cluster_key: str) -> Tuple[np.ndarray, Mapping[str, str]]:
    if cluster_key not in adata.obs:
        raise KeyError(f"Unable to find data in `adata.obs[{cluster_key!r}].`")
    if not is_categorical_dtype(adata.obs[cluster_key]):
        raise TypeError(
            f"Expected `adata.obs[{cluster_key!r}]` to be categorical, "
            f"found `{infer_dtype(adata.obs[cluster_key])}`."
        )

    color_key = f"{cluster_key}_colors"
    try:
        colors = adata.uns[color_key]
    except KeyError:
        adata.uns[color_key] = colors = _create_categorical_colors(len(adata.obs[cluster_key].cat.categories))
    mapper = dict(zip(adata.obs[cluster_key].cat.categories, colors))
    mapper[np.nan] = "grey"

    return colors, mapper


def _get_sorted_colors(
    adata: AnnData,
    cluster_key: Union[str, Sequence[str]],
    time_key: Optional[str] = None,
    tmin: float = -np.inf,
    tmax: float = np.inf,
) -> List[np.ndarray]:
    if time_key is not None:
        if time_key not in adata.obs:
            raise KeyError(f"Unable to find time in `adata.obs[{time_key!r}]`.")
        adata = adata[(adata.obs[time_key] >= tmin) & (adata.obs[time_key] <= tmax)]
        if not adata.n_obs:
            raise ValueError(f"Specified time range `{[tmin, tmax]}` does not contain any data.")
        order = np.argsort(adata.obs[time_key].values)
    else:
        order = np.arange(adata.n_obs)

    if isinstance(cluster_key, str):
        cluster_key = (cluster_key,)
    cluster_key = _unique_order_preserving(cluster_key)

    res = []
    for ck in cluster_key:
        try:
            cols, mapper = _get_categorical_colors(adata, ck)
            res.append(np.array([colors.to_hex(mapper[v]) for v in adata.obs[ck].values[order]]))
        except TypeError:
            if not is_numeric_dtype(adata.obs[ck]):
                raise TypeError(
                    f"Expected `adata.obs[{cluster_key!r}]` to be numeric, "
                    f"found `{infer_dtype(adata.obs[cluster_key])}`."
                ) from None
            res.append(np.asarray(adata.obs[ck])[order])

    return res
