# -*- coding: utf-8 -*-
"""Module used for finding initial and terminal states."""
from typing import Union, TypeVar, Optional, Sequence

from cellrank.ul._docs import d, _initial, _terminal, inject_docs
from cellrank.tl._constants import FinalStatesKey, FinalStatesPlot
from cellrank.tl.estimators import GPCCA
from cellrank.tl.estimators._constants import P
from cellrank.tl.kernels._precomputed_kernel import DummyKernel

AnnData = TypeVar("AnnData")


_find_docs = """\
Plot {direction} states uncovered by :class:`cellrank.tl.{fn_name}`.

Parameters
----------
%(adata)s
estimator
    Estimator class that was used to compute the {direction} states.
discrete
    If `True`, plot probability distribution of {direction} states.
    Only available when {direction} were estimated by :class:`cellrank.tl.estimators.GPCCA`.
title
    Title of the figure. If `None`, it will be set to `{title!r}`.
**kwargs
    Keyword arguments for :meth:`cellrank.tl.estimators.BaseEstimator.plot_final_states`.

Returns
-------
%(just_plots)s
"""


def _initial_terminal(
    adata: AnnData,
    backward: bool = False,
    discrete: bool = False,
    title: Optional[Union[str, Sequence[str]]] = None,
    **kwargs,
) -> None:

    pk = DummyKernel(adata=adata, backward=backward)
    mc = GPCCA(pk, read_from_adata=True, write_to_adata=False)

    if title is None and kwargs.get("same_plot", True):
        title = FinalStatesPlot.BACKWARD.s if backward else FinalStatesPlot.FORWARD.s

    if mc._get(P.FIN) is None:
        raise RuntimeError(
            f"Compute {_initial if backward else _terminal} states first as "
            f"`cellrank.tl.compute_{FinalStatesKey.BACKWARD if backward else FinalStatesKey.FORWARD}()`."
        )

    mc.plot_final_states(discrete=discrete, title=title, **kwargs)


@d.dedent
@inject_docs(
    __doc__=_find_docs.format(
        direction=_initial,
        fn_name=FinalStatesKey.BACKWARD.s,
        title=FinalStatesPlot.BACKWARD.s,
    )
)
def initial_states(
    adata: AnnData,
    discrete: bool = True,
    title: Optional[Union[str, Sequence[str]]] = None,
    **kwargs,
) -> Optional[AnnData]:  # noqa

    return _initial_terminal(
        adata,
        backward=True,
        discrete=discrete,
        title=title,
        **kwargs,
    )


@d.dedent
@inject_docs(
    __doc__=_find_docs.format(
        direction=_terminal,
        fn_name=FinalStatesKey.FORWARD.s,
        title=FinalStatesPlot.FORWARD.s,
    )
)
def terminal_states(
    adata: AnnData,
    discrete: bool = True,
    title: Optional[Union[str, Sequence[str]]] = None,
    **kwargs,
) -> Optional[AnnData]:  # noqa

    return _initial_terminal(
        adata,
        backward=False,
        discrete=discrete,
        title=title,
        **kwargs,
    )
