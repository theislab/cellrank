"""Utility functions which deal with delegating methods."""
import inspect
from typing import *  # noqa
from functools import partial, singledispatch, update_wrapper
from collections import namedtuple

import wrapt

import cellrank.logging as logg
from cellrank.tl.estimators._constants import F

AnnData = TypeVar("AnnData")
try:
    Metadata = namedtuple(
        "Metadata",
        ["attr", "prop", "dtype", "default", "doc", "compute_fmt", "plot_fmt"],
        defaults=[None, None, None, None, None, F.COMPUTE, F.PLOT],
    )
except TypeError:
    # Python < 3.7
    Metadata = namedtuple(
        "Metadata",
        ["attr", "prop", "dtype", "default", "doc", "compute_fmt", "plot_fmt"],
    )
    Metadata.__new__.__defaults__ = defaults = (
        None,
        None,
        None,
        None,
        None,
        F.COMPUTE,
        F.PLOT,
    )


# copy of functools.singledispatchmethod (for Python < 3.8)
class singledispatchmethod:
    """
    Single-dispatch generic method descriptor.

    Supports wrapping existing descriptors and handles non-descriptor callables as instance methods.
    """

    def __init__(self, func):
        if not callable(func) and not hasattr(func, "__get__"):
            raise TypeError(f"{func!r} is not callable or a descriptor")

        self.dispatcher = singledispatch(func)
        self.func = func

    def register(self, cls, method=None):
        """
        Generic_method.register(cls, func) -> func

        Registers a new implementation for the given *cls* on a *generic_method*.
        """  # noqa
        return self.dispatcher.register(cls, func=method)

    def __get__(self, obj, cls=None):
        def _method(*args, **kwargs):
            method = self.dispatcher.dispatch(args[0].__class__)
            return method.__get__(obj, cls)(*args, **kwargs)

        _method.__isabstractmethod__ = self.__isabstractmethod__
        _method.register = self.register
        update_wrapper(_method, self.func)
        return _method

    @property
    def __isabstractmethod__(self):
        return getattr(self.func, "__isabstractmethod__", False)


def argspec_factory(
    wrapped: Callable,
    skip: int = 0,
    return_type: Optional[Type] = None,
    fn: Optional[Callable] = None,
) -> Callable:
    """
    Create a dummy adapter function with desired signature.

    Parameters
    ----------
    wrapped
        Function being wrapped whose signature we wish to replicate.
    skip
        Number of arguments from left to skip.
    return_type
        Return type of the function.
    fn
        If not `None`, take the signature from this function.

    Returns
    -------
    :class:`Callable`
        An adapter with the correct signature.
    """
    # for locals(), this is a whitelist of types we allow
    import matplotlib  # noqa  this one seems to be missing, so whitelist it

    NoneType = type(None)

    sig = inspect.signature(wrapped if fn is None else fn)
    if return_type is not None:
        sig._return_annotation = return_type

    keys = list(sig.parameters.keys())[skip:]
    sig = sig.replace(parameters=[sig.parameters[k] for k in keys])

    exec(f"def adapter{sig}: pass", globals(), locals())

    return locals()["adapter"]


def _delegate(
    *,
    prop_name: Optional[str] = None,
    return_type: Optional[Type] = None,
    skip: int = 2,
) -> Callable:
    @wrapt.decorator()
    def pass_through(wrapped, _instance, args, kwargs):
        return wrapped(*args, **kwargs)

    if prop_name is None:
        return pass_through

    adapter = wrapt.adapter_factory(
        partial(argspec_factory, skip=skip, return_type=return_type)
    )

    @wrapt.decorator(adapter=adapter)
    def wrapper(wrapped, instance, args, kwargs):
        return wrapped(getattr(instance, prop_name), prop_name, *args, **kwargs)

    return wrapper


def _delegate_method_dispatch(fn: Callable, attr: str, prop_name: str, skip: int = 2):

    adapter = wrapt.adapter_factory(partial(argspec_factory, fn=fn, skip=skip))

    @wrapt.decorator(adapter=adapter)
    def delegate(wrapped, instance, args, kwargs):
        # first, get the correct method based on the dispatcher (can't use wrapped.dispatch..., since it's a function)
        # then, call the actual function
        return getattr(instance, attr)(
            getattr(instance, prop_name), prop_name, *args, **kwargs
        )

    return delegate(fn)


def _create_property(
    attr_name: str,
    prop_name: str,
    doc: Optional[str] = None,
    return_type: Optional[Type] = None,
) -> property:
    def wrapper(self) -> return_type:
        return getattr(self, attr_name)

    if not doc:
        doc = f"{prop_name.replace('_', ' ').capitalize()}."
    wrapper.__doc__ = doc

    return property(wrapper)


def _print_insufficient_number_of_cells(groups: Iterable[Any], n_cells: int):
    if groups:
        logg.debug(
            f"The following groups have less than requested number of cells ({n_cells}): "
            f"`{', '.join(sorted(map(str, groups)))}`"
        )
