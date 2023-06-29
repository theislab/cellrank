import abc
import contextlib
import copy as copy_
import inspect
from typing import (
    Any,
    Callable,
    Dict,
    List,
    Literal,
    Mapping,
    Optional,
    Sequence,
    Tuple,
    Union,
)

import numpy as np
import pandas as pd
import scipy.sparse as sp

from anndata import AnnData

from cellrank._utils._docs import d
from cellrank._utils._key import Key
from cellrank._utils._lineage import Lineage
from cellrank.estimators.mixins import KernelMixin
from cellrank.kernels import PrecomputedKernel
from cellrank.kernels._base_kernel import KernelExpression
from cellrank.kernels.mixins import AnnDataMixin, IOMixin

__all__ = ["BaseEstimator"]

Attr_t = (Literal["X", "raw", "layers", "obs", "var", "obsm", "varm", "obsp", "varp", "uns"],)


@d.get_sections(base="base_estimator", sections=["Parameters"])
class BaseEstimator(IOMixin, KernelMixin, AnnDataMixin, abc.ABC):
    """Base class for all estimators.

    Parameters
    ----------
    object
        Can be one of the following types:

        - :class:`~anndata.AnnData` - annotated data object.
        - :class:`~scipy.sparse.spmatrix`, :class:`~numpy.ndarray` - row-normalized transition matrix.
        - :class:`~cellrank.kernels.KernelExpression` - kernel expression.
        - :class:`str` - key in :attr:`~anndata.AnnData.obsp` where the transition matrix is stored and
          ``adata`` must be provided in this case.
        - :class:`bool` - directionality of the transition matrix that will be used to infer its storage location.
          If :obj:`None`, the directionality will be determined automatically and
          ``adata`` must be provided in this case.
    kwargs
        Keyword arguments for the :class:`~cellrank.kernels.PrecomputedKernel`.
    """

    def __init__(
        self,
        object: Union[str, bool, np.ndarray, sp.spmatrix, AnnData, KernelExpression],
        **kwargs: Any,
    ):
        if isinstance(object, KernelExpression):
            if object.transition_matrix is None:
                raise RuntimeError("Compute transition matrix first as `.compute_transition_matrix()`.")
        else:
            object = PrecomputedKernel(object, copy=False, **kwargs)
        super().__init__(kernel=object)

        self._params: Dict[str, Any] = {}
        self._shadow_adata = AnnData(
            X=sp.csr_matrix(self.adata.shape, dtype=self.adata.X.dtype),
            obs=self.adata.obs[[]].copy(),
            var=self.adata.var[[]].copy(),
            raw=None if self.adata.raw is None else self.adata.raw.to_adata(),
        )

    def __init_subclass__(cls, **kwargs: Any):
        super().__init_subclass__()

    @abc.abstractmethod
    def fit(self, *args: Any, **kwargs: Any) -> "BaseEstimator":
        """Fit the estimator.

        Parameters
        ----------
        args
            Positional arguments.
        kwargs
            Keyword arguments.

        Returns
        -------
        Self.
        """

    @abc.abstractmethod
    def predict(self, *args: Any, **kwargs: Any) -> "BaseEstimator":
        """Run the prediction.

        Parameters
        ----------
        args
            Positional arguments.
        kwargs
            Keyword arguments.

        Returns
        -------
        Self.
        """

    def _set(
        self,
        attr: Optional[str] = None,
        obj: Optional[Union[pd.DataFrame, Mapping[str, Any]]] = None,
        key: Optional[str] = None,
        value: Optional[Union[np.ndarray, pd.Series, pd.DataFrame, Lineage, AnnData, Dict[str, Any]]] = None,
        copy: bool = True,
        shadow_only: bool = False,
    ) -> None:
        """Set an attribute and optionally update ``obj['{key}']``.

        Parameters
        ----------
        attr
            Attribute to set. Only updated when we're not in the shadow. If :obj:`None`, don't update anything.
            See :attr:`_in_shadow` and ``obj`` for more information.
        obj
            Object which to update with ``value`` alongside the ``attr``.
            Usually, an attribute of :attr:`adata` is passed here.
        key
            Key in ``obj`` to update with ``value``. Only used when ``obj != None``.
        value
            Value to set. If :obj:`None` and ``key != None``, it removes the values under ``obj['{key}']``, if present.
        copy
            Whether to copy the ``value`` before setting it in ``obj``.
        shadow_only
            Whether to update the ``obj`` if we are not in the shadow.

        Returns
        -------
        Nothing, just optionally updates ``attr`` and/or ``obj[{key}]``.

        Raises
        ------
        AttributeError
            If ``attr`` doesn't exist.
        """
        if not self._in_shadow:
            if attr is not None:
                if not hasattr(self, attr):
                    raise AttributeError(attr)
                setattr(self, attr, value)
            if shadow_only:
                return

        if obj is None:
            return

        if key is not None:
            if value is None:
                with contextlib.suppress(KeyError):
                    del obj[key]

            else:
                obj[key] = copy_.copy(value) if copy else value

    def _get(
        self,
        *,
        obj: Union[pd.DataFrame, Mapping[str, Any]],
        key: str,
        shadow_attr: Optional[Literal["obs", "obsm", "var", "varm", "uns"]] = None,
        dtype: Optional[Union[type, Tuple[type, ...]]] = None,
        copy: bool = True,
        allow_missing: bool = False,
    ) -> Any:
        """Get data from an object and set an attribute.

        Parameters
        ----------
        obj
            Object from which to extract the data.
        key
            Key in ``obj`` where the data is stored.
        shadow_attr
            Attribute of :attr:`_shadow_adata` where to save the extracted data. If :obj:`None`, don't update it.
        dtype
            Valid type(s) of the extracted data.
        copy
            Copy the data before setting the ``self_attr``.
        allow_missing
            Whether to allow ``key`` to be missing in ``obj``.

        Returns
        -------
        The extracted values and optionally updates :attr:`_shadow_adata` ``.{shadow_attr}``.

        Raises
        ------
        AttributeError
            If ``attr`` doesn't exist.
        TypeError
            If ``dtype != None`` and the extracted values are not instances of ``dtype``.
        KeyError
            If ``allow_missing = False`` and ``key`` was not found in ``obj``.
        """
        if shadow_attr is not None and not hasattr(self._shadow_adata, shadow_attr):
            raise AttributeError(shadow_attr)

        try:
            data = obj[key]
            if dtype is not None and not isinstance(data, dtype):
                raise TypeError(f"Expected object to be of type `{dtype}`, found `{type(data).__name__}`.")
            if copy:
                data = copy_.copy(data)
            if shadow_attr is not None:
                getattr(self._shadow_adata, shadow_attr)[key] = data
            return data
        except KeyError:
            if not allow_missing:
                raise
            return None

    @property
    @contextlib.contextmanager
    def _shadow(self) -> None:
        """Temporarily set :attr:`adata` to :attr:`_shadow_adata`.

        Used to construct the serialization object in :meth:`to_adata`.
        """
        if self._in_shadow:
            yield
        else:
            adata = self.adata
            try:
                self.adata = self._shadow_adata
                yield
            finally:
                self.adata = adata

    @property
    def _in_shadow(self) -> bool:
        """Return `True` if :attr:`adata` is :attr:`_shadow_adata`."""
        return self.adata is self._shadow_adata

    def _create_params(
        self,
        locs: Optional[Mapping[str, Any]] = None,
        func: Optional[Callable] = None,
        remove: Sequence[str] = (),
    ) -> Dict[str, Any]:
        """Create parameters of interest from a function call.

        Parameters
        ----------
        locs
            Environment from which to get the parameters. If `None`, get the caller's environment.
        func
            Function of interest. If :obj:`None`, use the caller.
        remove
            Keys in ``locs`` which should not be included in the result.

        Returns
        -------
        The parameters as a :class:`dict`.

        Notes
        -----
        *args/**kwargs are always ignored and the values in ``locs`` are not copied.
        """
        frame = inspect.currentframe()
        try:
            if locs is None:
                locs = frame.f_back.f_locals
            if func is None:
                name = frame.f_back.f_code.co_name
                func = dict(inspect.getmembers(self)).get(name, None)
            if not callable(func):
                raise TypeError(f"Expected `func` to be `callable`, found `{type(func).__name__}`.")

            params = {}
            for name, param in inspect.signature(func).parameters.items():
                if param.kind in (inspect.Parameter.VAR_POSITIONAL, inspect.Parameter.VAR_KEYWORD):
                    continue
                if name in remove:
                    continue
                if name in locs:
                    params[name] = locs[name]

            return params
        except AttributeError:
            # frame can be None
            return {}
        except TypeError:
            return {}
        finally:
            del frame

    def _read_params(self, key: str) -> Dict[str, Any]:
        """Read ``key`` from estimator params in :attr:`adata`.

        Usually called in :meth:`_read_adata` during :meth:`from_adata`.
        """
        ekey = Key.uns.estimator(self.backward) + "_params"
        return dict(self.adata.uns.get(ekey, {}).get(key, {}))

    @d.dedent
    def to_adata(
        self,
        keep: Union[Literal["all"], Sequence[Attr_t]] = ("X", "raw"),
        *,
        copy: Union[bool, Sequence[Attr_t]] = True,
    ) -> AnnData:
        """%(to_adata.full_desc)s

        Parameters
        ----------
        keep
            Which attributes to keep from the underlying :attr:`adata`. Valid options are:

            - ``'all'`` - keep all attributes specified in the signature.
            - :class:`~typing.Sequence` - keep only subset of these attributes.
            - :class:`dict` - the keys correspond the attribute names and values to a subset of keys
              which to keep from this attribute. If the values are specified either as :obj:`True` or ``'all'``,
              everything from this attribute will be kept.
        copy
            Whether to copy the data. Can be specified on per-attribute basis.
            Useful for attributes that are array-like.

        Returns
        -------
        Annotated data object.
        """  # noqa: D400

        def handle_attribute(attr: Attr_t, keys: List[str], *, copy: bool) -> None:
            try:
                if attr == "X":
                    adata.X = copy_.deepcopy(self.adata.X) if copy else self.adata.X
                    return
                if attr == "raw":
                    adata.raw = self.adata.raw.to_adata() if self.adata.raw is not None else None
                    return

                old = getattr(self.adata, attr)
                new = getattr(adata, attr)
                if keys == ["all"]:
                    keys = list(old.keys())

                # fmt: off
                if isinstance(new, pd.DataFrame):
                    old = old[keys]
                    # avoid duplicates
                    old = old[old.columns.difference(new.columns)]
                    setattr(adata, attr, pd.merge(new, old, how="inner", left_index=True, right_index=True, copy=copy))
                elif isinstance(new, Mapping):
                    old = {k: old[k] for k in keys}
                    # old has preference, since it's user supplied
                    setattr(adata, attr, {**new, **(copy_.deepcopy(old) if copy else old)})
                else:
                    raise TypeError(f"Expected `adata.{attr}` to be either `Mapping` or `pandas. DataFrame`, "
                                    f"found `{type(new).__name__}`.")
                # fmt: on
            except KeyError:
                missing = sorted(k for k in keys if k not in old)
                raise KeyError(f"Unable to find key(s) `{missing}` in `adata.{attr}`.") from None

        adata = self._shadow_adata.copy()
        _adata = self.adata
        try:
            # kernel and estimator share the adata
            self.adata = adata
            self.kernel.write_to_adata()
        finally:
            self.adata = _adata
        key = Key.uns.estimator(self.backward) + "_params"
        adata.uns[key] = copy_.deepcopy(self.params)

        # fmt: off
        if isinstance(keep, str):
            if keep == "all":
                keep = ["X", "raw", "layers", "obs", "var", "obsm", "varm", "obsp", "varp", "uns"]
            else:
                keep = [keep]
        if not isinstance(keep, Mapping):
            keep = {attr: True for attr in keep}
        if isinstance(copy, bool):
            copy = {attr: copy for attr in keep}
        elif isinstance(copy, str):
            copy = [copy]
        if not isinstance(copy, Mapping):
            copy = {attr: True for attr in copy}
        # fmt: on

        for attr, keys in keep.items():
            if keys is True:
                keys = ["all"]
            elif isinstance(keys, str):
                keys = [keys]
            if keys is False or not len(keys):
                continue
            handle_attribute(attr, keys=list(keys), copy=copy.get(attr, False))

        return adata

    @classmethod
    @d.dedent
    def from_adata(cls, adata: AnnData, obsp_key: str) -> "BaseEstimator":
        """%(from_adata.full_desc)s

        Parameters
        ----------
        %(adata)s
        obsp_key
            Key in :attr:`~anndata.AnnData.obsp` where the transition matrix is stored.

        Returns
        -------
        %(from_adata.returns)s
        """  # noqa: D400
        return super().from_adata(adata, obsp_key=obsp_key)

    @d.dedent
    def copy(self, *, deep: bool = False) -> "BaseEstimator":
        """Return a copy of self.

        Parameters
        ----------
        deep
            Whether to return a deep copy or not. If :obj:`True`, this also copies the :attr:`adata`.

        Returns
        -------
        A copy of self.
        """
        k = copy_.deepcopy(self.kernel) if deep else copy_.copy(self.kernel)
        res = type(self)(k)
        for k, v in self.__dict__.items():
            if isinstance(v, Mapping):
                res.__dict__[k] = copy_.deepcopy(v)
            elif k != "_kernel":
                res.__dict__[k] = copy_.deepcopy(v) if deep else copy_.copy(v)

        return res

    def __copy__(self) -> "BaseEstimator":
        return self.copy(deep=False)

    def _format_params(self) -> str:
        return f"kernel={self.kernel!s}"

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}[{self._format_params()}]"

    def __str__(self) -> str:
        return repr(self)

    @property
    def params(self) -> Dict[str, Any]:
        """Estimator parameters."""
        return self._params
