# -*- coding: utf-8 -*-
"""Base properties used within the estimators."""
import sys
from abc import ABC, ABCMeta, abstractmethod
from typing import Any, Dict, List, Tuple, Union, Iterable, Optional
from inspect import isabstract

import scvelo as scv
from anndata import AnnData

import numpy as np
import pandas as pd
from scipy.sparse import issparse, spmatrix
from pandas.api.types import is_categorical_dtype

import matplotlib as mpl
from matplotlib import cm

import cellrank.logging as logg
from cellrank.ul._docs import d, _initial, _terminal
from cellrank.tl._utils import RandomKeys, _make_cat, _partition, _complex_warning
from cellrank.ul._utils import Pickleable
from cellrank.tl._colors import _create_categorical_colors
from cellrank.tl.kernels import PrecomputedKernel
from cellrank.tl._lineage import Lineage
from cellrank.tl._constants import Direction, DirPrefix, DirectionPlot
from cellrank.tl.kernels._utils import _filter_kwargs
from cellrank.tl.estimators._utils import (
    Metadata,
    Sequence,
    _create_property,
    singledispatchmethod,
    _delegate_method_dispatch,
)
from cellrank.tl.kernels._base_kernel import KernelExpression
from cellrank.tl.estimators._constants import META_KEY, A, F, P


def is_abstract(classname: str) -> bool:  # TODO: determine the necessity of this
    """
    Check whether class with a given name inside this module is abstract.

    Parameters
    ----------
    classname
        Name of the class.

    Returns
    -------
    bool
        `True` if the class is abstract, otherwise `False`.

    """

    cls = getattr(sys.modules[__name__], classname, None)
    return cls is not None and isabstract(cls)


class PropertyMeta(ABCMeta, type):
    """Metaclass for all the properties."""

    @staticmethod
    def update_attributes(md: Metadata, attributedict: Dict[str, Any]) -> Optional[str]:
        """
        Update ``attributedict`` with new attribute and property.

        Parameters
        ----------
        md
            Metadata object, containing attribute name, property names, etc.
        attributedict
            Dictionary with attributes.

        Returns
        -------
        str or None
            Old property name, if found, otherwise a newly constructed one, or `None` if no property is desired.
        """

        # TODO: determine whether supporting strings is such a good idea
        if not isinstance(md, Metadata):
            raise TypeError(
                f"Expected a `Metadata` object, found `{type(md).__name__}`."
            )

        if not isinstance(md.attr, (str, A)):
            raise TypeError(
                f"Attribute `{md.attr}` must be of type `A` or `str`, found `{type(md.attr).__name__!r}`."
            )
        if not str(md.attr).startswith("_"):
            raise ValueError(f"Attribute `{md.attr!r}` must start with `'_'`.")

        attributedict[str(md.attr)] = md.default

        if md.prop == P.NO_PROPERTY:
            return

        if not isinstance(md.prop, (str, P)):
            raise TypeError(
                f"Property must be of type `P` or `str`, found `{type(md.attr).__name__!r}`."
            )

        prop_name = str(md.attr).lstrip("_") if md.prop == P.EMPTY else str(md.prop)

        if not len(prop_name):
            raise ValueError(f"Property name for attribute `{md.attr}` is empty.")
        if prop_name.startswith("_"):
            raise ValueError(
                f"Property musn't start with an underscore: `{prop_name!r}`."
            )

        attributedict[prop_name] = _create_property(
            str(md.attr), prop_name, doc=md.doc, return_type=md.dtype
        )

        return prop_name

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

        compute_md, metadata = attributedict.pop(META_KEY, None), []

        if compute_md is None:
            return super().__new__(cls, clsname, superclasses, attributedict)

        if isinstance(compute_md, str):
            compute_md = Metadata(attr=compute_md)
        elif not isinstance(compute_md, (tuple, list)):
            raise TypeError(
                f"Expected property metadata to be `list` or `tuple`,"
                f"found `{type(compute_md).__name__!r}`."
            )
        elif len(compute_md) == 0:
            raise ValueError("No metadata found.")
        else:
            compute_md, *metadata = [
                Metadata(attr=md) if isinstance(md, str) else md for md in compute_md
            ]

        prop_name = PropertyMeta.update_attributes(compute_md, attributedict)
        plot_name = str(compute_md.plot_fmt).format(prop_name)

        if compute_md.compute_fmt != F.NO_FUNC:
            if "_compute" in attributedict:
                attributedict[
                    str(compute_md.compute_fmt).format(prop_name)
                ] = attributedict["_compute"]

        if (
            compute_md.plot_fmt != F.NO_FUNC
            and VectorPlottable in superclasses
            and plot_name not in attributedict
            and not is_abstract(clsname)
        ):
            raise TypeError(
                f"Method `{plot_name}` is not implemented for class `{clsname}`."
            )

        for md in metadata:
            PropertyMeta.update_attributes(md, attributedict)

        res = super().__new__(cls, clsname, superclasses, attributedict)

        if compute_md.plot_fmt != F.NO_FUNC and Plottable in res.mro():
            # _this is intended singledispatchmethod
            # unfortunately, `_plot` is not always in attributedict, so we can't just check for it
            # and res._plot is just a regular function
            # if this gets buggy in the future, consider switching from singlemethoddispatch
            setattr(
                res,
                plot_name,
                _delegate_method_dispatch(res._plot, "_plot", prop_name, skip=2),
            )

        return res


class Property(ABC, metaclass=PropertyMeta):
    """Base class for all the properties."""

    pass


class KernelHolder(ABC):
    """Base class which holds a :class:`cellrank.tool.kernels._kernel.KernelExpression`."""

    def __init__(
        self,
        obj: Union[AnnData, np.ndarray, spmatrix, KernelExpression],
        key: Optional[str] = None,
        obsp_key: Optional[str] = None,
        write_to_adata: bool = True,
    ):
        if isinstance(obj, KernelExpression):
            self._kernel = obj
        elif isinstance(obj, (np.ndarray, spmatrix)):
            self._kernel = PrecomputedKernel(obj)
        elif isinstance(obj, AnnData):
            if obsp_key is None:
                raise ValueError(
                    "Specify `obsp_key=...` when supplying an `AnnData` object."
                )
            elif obsp_key not in obj.obsp.keys():
                raise KeyError(f"Key `{obsp_key!r}` not found in `adata.obsp`.")
            self._kernel = PrecomputedKernel(obj.obsp[obsp_key], adata=obj)
        else:
            raise TypeError(
                f"Expected an object of type `KernelExpression`, `numpy.ndarray`, `scipy.sparse.spmatrix` "
                f"or `anndata.AnnData`, got `{type(obj).__name__!r}`."
            )

        if self.kernel._transition_matrix is None:
            # access the private attribute to avoid accidentally computing the transition matrix
            # in principle, it doesn't make a difference, apart from not seeing the message
            logg.warning("Computing transition matrix using the default parameters")
            self.kernel.compute_transition_matrix()

        if write_to_adata:
            self.kernel.write_to_adata(key=key)

    @property
    def _direction(self):
        return Direction.BACKWARD if self.kernel.backward else Direction.FORWARD

    @property
    def transition_matrix(self) -> Union[np.ndarray, spmatrix]:
        """Transition matrix."""
        return self.kernel.transition_matrix

    @property
    def issparse(self) -> bool:
        """Whether the transition matrix is sparse or not."""
        return issparse(self.transition_matrix)

    @property
    def kernel(self) -> KernelExpression:
        """Underlying kernel."""
        return self._kernel

    @property
    @d.dedent
    def adata(self) -> AnnData:
        """
        Annotated data object.

        Returns
        -------
        %(adata_ret)s
        """  # noqa
        return self.kernel.adata

    def __len__(self):
        return self.kernel.transition_matrix.shape[0]

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}[n={len(self)}, kernel={repr(self.kernel)}]"

    def __str__(self) -> str:
        return f"{self.__class__.__name__}[n={len(self)}, kernel={str(self.kernel)}]"


class VectorPlottable(KernelHolder, Property):
    """
    Injector class which plots vectors.

    To be used in conjunction with:

        - :class:`cellrank.tool.estimators._decomposition.Eig`
        - :class:`cellrank.tool.estimators._decomposition.Schur`
    """

    @d.get_sectionsf("plot_vectors", sections=["Parameters", "Returns"])
    @d.dedent
    def _plot_vectors(
        self,
        vectors: Optional[np.ndarray],
        prop: str,
        use: Optional[Union[int, Tuple[int], List[int]]] = None,
        abs_value: bool = False,
        cluster_key: Optional[str] = None,
        **kwargs,
    ) -> None:
        """
        Plot vectors in an embedding.

        Parameters
        ----------
        use
            Which or how many vectors are to be plotted.
        abs_value
            Whether to take the absolute value before plotting.
        cluster_key
            Key in :paramref:`adata` ``.obs`` for plotting categorical observations.
        %(basis)s
        **kwargs
            Keyword arguments for :func:`scvelo.pl.scatter`.

        Returns
        -------
        %(just_plots)s
        """

        if prop not in (P.EIG.v, P.SCHUR.v):
            raise ValueError(
                f"Invalid kind `{prop!r}`. Valid options are `{P.EIG!r}`, `{P.SCHUR!r}``."
            )
        if vectors is None:
            raise RuntimeError(f"Compute `.{prop}` first as `{F.COMPUTE.fmt(prop)}()`.")

        if prop == P.SCHUR.s:
            is_schur = True
            name = "Schur "
        else:
            is_schur = False
            name = "eigen"

        # check whether dimensions are consistent
        if self.adata.n_obs != vectors.shape[0]:
            raise ValueError(
                f"Number of cells ({self.adata.n_obs}) is inconsistent with the first"
                f"dimension of vectors ({vectors.shape[0]})."
            )

        if use is None:
            m = (
                getattr(self, P.EIG.s).get("eigengap", vectors.shape[1]) + 1
                if hasattr(self, P.EIG.s) and not is_schur
                else vectors.shape[1]
            )
            use = list(range(is_schur, m))
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
            raise ValueError("Nothing to plot.")

        muse = max(use)
        if muse >= vectors.shape[1]:
            raise ValueError(
                f"Maximum specified {name}vector ({muse}) is larger "
                f"than the number of computed {name}vectors ({vectors.shape[1]})."
            )
        V_ = vectors[:, use]

        if is_schur:
            title = [f"{name}vector {i}" for i in use]
        else:
            D = kwargs.pop("D")
            V_ = _complex_warning(V_, use, use_imag=kwargs.pop("use_imag", False))
            title = [fr"$\lambda_{i}$={d:.02f}" for i, d in zip(use, D[use])]

        if abs_value:
            V_ = np.abs(V_)

        color = list(V_.T)
        if cluster_key is not None:
            color = [cluster_key] + color
        cmap = kwargs.pop("cmap", "viridis")

        logg.debug(f"Plotting `{use}` {name}vectors")

        scv.pl.scatter(self.adata, color=color, title=title, cmap=cmap, **kwargs)


class Plottable(KernelHolder, Property):
    """
    Injector which plots macrostates or terminal states or absorption probabilities.

    To be used in conjunction with:

        - :class:`cellrank.tool.estimators._property.MacroStates`.
        - :class:`cellrank.tool.estimators._property.TermStates`.
        - :class:`cellrank.tool.estimators._property.AbsProbs`.
    """

    @d.get_sectionsf("plot_discrete")
    @d.dedent
    def _plot_discrete(
        self,
        data: pd.Series,
        prop: str,
        lineages: Optional[Union[str, Sequence[str]]] = None,
        cluster_key: Optional[str] = None,
        same_plot: bool = True,
        title: Optional[Union[str, List[str]]] = None,
        **kwargs,
    ) -> None:
        """
        Plot the states for each uncovered lineage.

        Parameters
        ----------
        lineages
            Plot only these lineages. If `None`, plot all lineages.
        cluster_key
            Key from :paramref:`adata` ``.obs`` for plotting categorical observations.
        same_plot
            Whether to plot the lineages on the same plot or separately.
        title
            The title of the plot.
        %(basis)s
        **kwargs
            Keyword arguments for :func:`scvelo.pl.scatter`.

        Returns
        -------
        %(just_plots)s
        """

        if data is None:
            raise RuntimeError(
                f"Compute `.{prop}` first as `.{F.COMPUTE.fmt(prop)}()`."
            )
        if not is_categorical_dtype(data):
            raise TypeError(
                f"Expected property `.{prop}` to be categorical, found `{type(data).__name__!r}`."
            )
        if prop in (P.ABS_PROBS.s, P.TERM.s):
            colors = getattr(self, A.TERM_COLORS.v, None)
        elif prop == P.MACRO.v:
            colors = getattr(self, A.MACRO_COLORS.v, None)
        else:
            logg.debug("No colors found. Creating new ones")
            colors = _create_categorical_colors(len(data.cat.categories))
        colors = dict(zip(data.cat.categories, colors))

        if (
            lineages is not None
        ):  # these are states per-se, but I want to keep the arg names for dispatch the same
            if isinstance(lineages, str):
                lineages = [lineages]
            for state in lineages:
                if state not in data.cat.categories:
                    raise ValueError(
                        f"Invalid state `{state!r}`. Valid options are `{list(data.cat.categories)}`."
                    )
            data = data.copy()
            to_remove = list(set(data.cat.categories) - set(lineages))

            if len(to_remove) == len(data.cat.categories):
                raise RuntimeError(
                    "Nothing to plot because empty subset has been selected."
                )

            for state in to_remove:
                data[data == state] = np.nan
            data.cat.remove_categories(to_remove, inplace=True)

        if cluster_key is None:
            cluster_key = []
        elif isinstance(cluster_key, str):
            cluster_key = [cluster_key]
        if not isinstance(cluster_key, list):
            cluster_key = list(cluster_key)

        same_plot = same_plot or len(data.cat.categories) == 1
        kwargs["legend_loc"] = kwargs.get("legend_loc", "on data")

        with RandomKeys(
            self.adata, None if same_plot else len(data.cat.categories), where="obs"
        ) as keys:
            if same_plot:
                key = keys[0]
                self.adata.obs[key] = data
                self.adata.uns[f"{key}_colors"] = [
                    colors[c] for c in data.cat.categories
                ]

                if title is None:
                    title = (
                        f"{prop.replace('_', ' ')} "
                        f"({Direction.BACKWARD if self.kernel.backward else Direction.FORWARD})"
                    )
                if isinstance(title, str):
                    title = [title]

                scv.pl.scatter(
                    self.adata,
                    title=cluster_key + title,
                    color=cluster_key + keys,
                    **_filter_kwargs(scv.pl.scatter, **kwargs),
                )
            else:
                for key, cat in zip(keys, data.cat.categories):
                    d = data.copy()
                    d[data != cat] = None
                    d.cat.set_categories([cat], inplace=True)

                    self.adata.obs[key] = d
                    self.adata.uns[f"{key}_colors"] = [colors[cat]]

                scv.pl.scatter(
                    self.adata,
                    color=cluster_key + keys,
                    title=(
                        cluster_key
                        + [
                            f"{_initial if self.kernel.backward else _terminal} state {c}"
                            for c in data.cat.categories
                        ]
                    )
                    if title is None
                    else title,
                    **_filter_kwargs(scv.pl.scatter, **kwargs),
                )

    @d.get_sectionsf("plot_continuous")
    @d.dedent
    def _plot_continuous(
        self,
        probs: Optional[Lineage],
        prop: str,
        diff_potential: Optional[pd.Series] = None,
        lineages: Optional[Union[str, Iterable[str]]] = None,
        cluster_key: Optional[str] = None,
        mode: str = "embedding",
        time_key: str = "latent_time",
        show_dp: bool = True,
        title: Optional[str] = None,
        same_plot: bool = True,
        cmap: Union[str, mpl.colors.ListedColormap] = cm.viridis,
        **kwargs,
    ) -> None:
        """
        Plot continuous observations such as macrostates memberships or lineages in an embedding.

        Parameters
        ----------
        lineages
            Plot only these lineages. If `None`, plot all lineages.
        cluster_key
            Key from :paramref:`adata` ``.obs`` for plotting categorical observations.
        %(time_mode)s
        time_key
            Key from :paramref:`adata` ``.obs`` to use as a pseudotime ordering of the cells.
        title
            Either `None`, in which case titles are ``'{to, from} {terminal, initial} {state}'``,
            or an array of titles, one per lineage.
        same_plot
            Whether to plot the lineages on the same plot using color gradients when ``mode='embedding'``.
        cmap
            Colormap to use.
        %(basis)s
        **kwargs
            Keyword arguments for :func:`scvelo.pl.scatter`.

        Returns
        -------
        %(just_plots)s
        """

        if probs is None:
            raise RuntimeError(
                f"Compute `.{prop}` first as `.{F.COMPUTE.fmt(prop)}()`."
            )

        if isinstance(lineages, str):
            lineages = [lineages]

        if lineages is None:
            lineages = probs.names
            A = probs
        else:
            A = probs[lineages]

        if not len(lineages):
            raise RuntimeError(
                "Nothing to plot because empty subset has been selected."
            )

        prefix = DirPrefix.BACKWARD if self.kernel.backward else DirPrefix.FORWARD
        same_plot = same_plot and mode == "embedding"  # set this silently

        diff_potential = (
            [diff_potential.values]
            if show_dp
            and not same_plot
            and diff_potential is not None
            and probs.shape[1] > 1
            else []
        )

        A = A.copy()  # the below code modifies stuff inplace
        X = A.X  # list(A.T) behaves differently, because it's Lineage

        if X.shape[1] == 1:
            same_plot = (
                False  # because color_gradients for 1 state is buggy (looks empty)
            )
            # this is the case for only 1 recurrent class - all cells have prob. 1 of going there
            # however, matplotlib's plotting really picks up the slightest differences in the colormap, here we set
            # everything to one, if applicable
            if np.allclose(X, 1.0):
                X = np.ones_like(X)

        for col in X.T:
            mask = ~np.isclose(col, 1.0)
            # change the maximum value - the 1 is artificial and obscures the color scaling
            if np.sum(mask):
                max_not_one = np.max(col[mask])
                col[~mask] = max_not_one

        if mode == "time":
            if time_key not in self.adata.obs.keys():
                raise KeyError(f"Time key `{time_key!r}` not found in `adata.obs`.")
            time = self.adata.obs[time_key]
            if cluster_key is not None:
                logg.warning(
                    f"Cluster key `{cluster_key!r}` is ignored when `mode='time'`"
                )
            cluster_key = None

        color = list(X.T) + diff_potential
        if title is None:
            if same_plot:
                title = [
                    f"{prop.replace('_', ' ')} "
                    f"({DirectionPlot.BACKWARD if self.kernel.backward else Direction.FORWARD})"
                ]
            else:
                title = [f"{prefix} {lin}" for lin in lineages] + (
                    ["differentiation potential"] if diff_potential else []
                )
        elif isinstance(title, str):
            title = [title]

        if isinstance(cluster_key, str):
            cluster_key = [cluster_key]
        elif cluster_key is None:
            cluster_key = []
        if not isinstance(cluster_key, list):
            cluster_key = list(cluster_key)

        if not same_plot:
            color = cluster_key + color
            title = cluster_key + title

        if mode == "embedding":
            if same_plot:
                kwargs["color_gradients"] = A
                if len(cluster_key):
                    logg.warning(
                        "Ignoring `cluster_key` when plotting continuous observations in the same plot"
                    )
                # kwargs["color"] = cluster_key  this results in a bug, cluster_key data is overwritten, will make a PR
            else:
                kwargs["color"] = color

            if probs.shape[1] == 1 and prop in (P.MACRO_MEMBER.s, P.TERM_PROBS.s):
                if "perc" not in kwargs:
                    logg.warning(
                        "Did not detect percentile for stationary distribution. Setting `perc=[0, 95]`"
                    )
                    kwargs["perc"] = [0, 95]
                kwargs["color"] = X
                kwargs.pop("color_gradients", None)

            scv.pl.scatter(
                self.adata,
                title=title,
                color_map=cmap,
                **_filter_kwargs(scv.pl.scatter, **kwargs),
            )
        elif mode == "time":
            scv.pl.scatter(
                self.adata,
                x=time,
                color_map=cmap,
                y=color,
                title=title,
                xlabel=[time_key] * len(title),
                ylabel=["probability"] * len(title),
                **_filter_kwargs(scv.pl.scatter, **kwargs),
            )
        else:
            raise ValueError(
                f"Invalid mode `{mode!r}`. Valid options are: `'embedding'` or `'time'`."
            )

    @singledispatchmethod
    @d.dedent
    def _plot(
        self,
        data,
        prop: str,
        discrete: bool = False,
        lineages: Optional[Union[str, Sequence[str]]] = None,
        cluster_key: Optional[str] = None,
        mode: str = "embedding",
        time_key: str = "latent_time",
        show_dp: bool = True,
        title: Optional[str] = None,
        same_plot: bool = False,
        cmap: Union[str, mpl.colors.ListedColormap] = "viridis",
        **kwargs,
    ) -> None:
        """
        Plot discrete states or probabilities in an embedding.

        Parameters
        ----------
        discrete
            Whether to plot in discrete or continuous mode.
        %(plot_continuous.parameters)s

        Returns
        -------
        %(just_plots)s
        """
        # we keep all the arguments because from this function, we get the docs + signature
        raise RuntimeError(f"Compute `.{prop}` first as `.{F.COMPUTE.fmt(prop)}()`.")

    @_plot.register(pd.Series)
    def _(self, data: pd.Series, prop: str, discrete: bool = False, **kwargs) -> None:
        if discrete and kwargs.get("mode", "embedding") == "time":
            logg.warning(
                "`mode='time'` is implemented in continuous case, plotting in continuous mode"
            )
            discrete = False

        if discrete:
            self._plot_discrete(data, prop, **kwargs)
        elif prop == P.MACRO.v:  # GPCCA
            prop = P.MACRO_MEMBER.v
            self._plot_continuous(getattr(self, prop, None), prop, None, **kwargs)
        elif prop == P.TERM.v:
            probs = getattr(self, A.TERM_ABS_PROBS.s, None)
            # we have this only in GPCCA
            if isinstance(probs, Lineage):
                self._plot_continuous(probs, P.TERM_PROBS.v, **kwargs)
            else:
                logg.warning(
                    f"Unable to plot continuous observations for `{prop!r}`, plotting in discrete mode"
                )
                self._plot_discrete(data, prop, **kwargs)
        else:
            raise NotImplementedError(
                f"Unable to plot property `.{prop}` in discrete mode."
            )

    @_plot.register(Lineage)
    def _(self, data: Lineage, prop: str, discrete: bool = False, **kwargs) -> None:
        if discrete and kwargs.get("mode", "embedding") == "time":
            logg.warning(
                "`mode='time'` is implemented in continuous case, plotting in continuous mode"
            )
            discrete = False

        if not discrete:
            diff_potential = getattr(self, P.DIFF_POT.v, None)
            self._plot_continuous(data, prop, diff_potential, **kwargs)
        elif prop == P.ABS_PROBS.v:
            # for discrete and abs. probs, plot the terminal states
            prop = P.TERM.v
            self._plot_discrete(getattr(self, prop, None), prop, **kwargs)
        else:
            raise NotImplementedError(
                f"Unable to plot property `.{prop}` in continuous mode."
            )


class Macrostates(Plottable):
    """Class dealing with macrostates."""

    __prop_metadata__ = [
        Metadata(attr=A.MACRO, prop=P.MACRO, dtype=pd.Series),
        Metadata(
            attr=A.MACRO_MEMBER,
            prop=P.MACRO_MEMBER,
            dtype=Lineage,
        ),
        Metadata(attr=A.MACRO_COLORS, prop=P.NO_PROPERTY, dtype=np.ndarray),
    ]

    @abstractmethod
    def compute_macrostates(self, *args, **kwargs) -> None:  # noqa
        pass


class TerminalStates(Plottable):
    """Class dealing with terminal states."""

    __prop_metadata__ = [
        Metadata(attr=A.TERM, prop=P.TERM, dtype=pd.Series, doc="Terminal states."),
        Metadata(
            attr=A.TERM_PROBS,
            prop=P.TERM_PROBS,
            dtype=pd.Series,
            doc="Terminal states probabilities.",
        ),
        Metadata(attr=A.TERM_COLORS, prop=P.NO_PROPERTY, dtype=np.ndarray),
    ]

    @abstractmethod
    def set_terminal_states(self, *args, **kwargs) -> None:  # noqa
        pass

    @abstractmethod
    def compute_terminal_states(self, *args, **kwargs) -> None:  # noqa
        pass

    # TODO: back compat

    @abstractmethod
    def _write_terminal_states(self, *args, **kwargs) -> None:
        pass


class AbsProbs(Plottable):
    """Class dealing with absorption probabilities."""

    __prop_metadata__ = [
        Metadata(
            attr=A.ABS_PROBS,
            prop=P.ABS_PROBS,
            dtype=Lineage,
            doc="Absorption probabilities.",
        ),
        Metadata(
            attr=A.DIFF_POT,
            prop=P.DIFF_POT,
            dtype=pd.Series,
            doc="Differentiation potential.",
        ),
        Metadata(attr=A.LIN_ABS_TIMES, prop=P.LIN_ABS_TIMES, dtype=pd.DataFrame),
    ]

    @abstractmethod
    def _write_absorption_probabilities(self, *args, **kwargs) -> None:
        pass


class LinDrivers(Plottable):  # noqa
    __prop_metadata__ = [
        Metadata(
            attr=A.LIN_DRIVERS,
            prop=P.LIN_DRIVERS,
            dtype=pd.DataFrame,
            doc="Lineage drivers.",
            plot_fmt=F.NO_FUNC,
            # in essence ignore Plottable (could be done by registering DataFrame, but it's ugly
        )
    ]


class Partitioner(KernelHolder, ABC):
    """Abstract base class for partitioning transition matrix into sets of recurrent and transient states."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._is_irreducible = None  # no need to refactor these
        self._rec_classes = None
        self._trans_classes = None

    def compute_partition(self) -> None:
        """
        Compute communication classes for the Markov chain.

        Returns
        -------
        None
            Nothing, but updates the following fields:

                - :paramref:`recurrent_classes`
                - :paramref:`transient_classes`
                - :paramref:`is_irreducible`
        """

        start = logg.info("Computing communication classes")
        n_states = len(self)

        rec_classes, trans_classes = _partition(self.transition_matrix)

        self._is_irreducible = len(rec_classes) == 1 and len(trans_classes) == 0

        if not self._is_irreducible:
            self._trans_classes = _make_cat(
                trans_classes, n_states, self.adata.obs_names
            )
            self._rec_classes = _make_cat(rec_classes, n_states, self.adata.obs_names)
            logg.info(
                f"Found `{(len(rec_classes))}` recurrent and `{len(trans_classes)}` transient classes\n"
                f"Adding `.recurrent_classes`\n"
                f"       `.transient_classes`\n"
                f"       `.is_irreducible`\n"
                f"    Finish",
                time=start,
            )
        else:
            logg.warning(
                "The transition matrix is irreducible, cannot further partition it\n    Finish",
                time=start,
            )

    @property
    def is_irreducible(self):
        """Whether the Markov chain is irreducible or not."""
        return self._is_irreducible

    @property
    def recurrent_classes(self):
        """Recurrent classes of the Markov chain."""  # noqa
        return self._rec_classes

    @property
    def transient_classes(self):
        """Transient classes of the Markov chain."""  # noqa
        return self._trans_classes


class LineageEstimatorMixin(TerminalStates, AbsProbs, LinDrivers, Pickleable, ABC):
    """Mixin containing terminal states and absorption probabilities."""

    pass
