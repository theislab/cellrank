# -*- coding: utf-8 -*-
"""Module used for finding initial and terminal states."""
from typing import TypeVar, Optional

from cellrank.ul._docs import d, _initial, _terminal, inject_docs
from cellrank.tl._utils import _check_estimator_type
from cellrank.tl._constants import FinalStatesKey
from cellrank.tl.estimators import GPCCA
from cellrank.tl.estimators._constants import P
from cellrank.tl.estimators._base_estimator import BaseEstimator
from cellrank.tl.kernels._precomputed_kernel import DummyKernel

AnnData = TypeVar("AnnData")


_find_docs = """\
Plot direction} states uncovered by :class:`cellrank.tl.{fn_name}`.

Parameters
----------
%(adata)s
estimator
    Estimator class that was used to compute the {direction} states.
discrete
    If `True`, plot probability distribution of {direction} states.
    Only available when ``estimator`` is :class:`cellrank.tl.estimators.GPCCA` type.
**kwargs
    Keyword arguments for :meth:`cellrank.tl.estimators.BaseEstimator.plot_final_states`.

Returns
-------
%(just_plots)s
"""


def _initial_terminal(
    adata: AnnData,
    estimator: type(BaseEstimator) = GPCCA,
    backward: bool = False,
    discrete: bool = False,
    **kwargs,
) -> None:

    _check_estimator_type(estimator)

    pk = DummyKernel(adata=adata, backward=backward)
    mc = estimator(pk, read_from_adata=True, write_to_adata=False)

    if mc._get(P.FIN) is None:
        raise RuntimeError(
            f"Compute {_initial if backward else _terminal} states first as "
            f"`cellrank.tl.compute_{FinalStatesKey.BACKWARD if backward else FinalStatesKey.FORWARD}`."
        )

    mc.plot_final_states(discrete=discrete, **kwargs)


@d.dedent
@inject_docs(
    __doc__=_find_docs.format(direction=_initial, fn_name=FinalStatesKey.BACKWARD.s)
)
def initial_states(
    adata: AnnData,
    estimator: type(BaseEstimator) = GPCCA,
    discrete: bool = True,
    **kwargs,
) -> Optional[AnnData]:  # noqa

    return _initial_terminal(
        adata,
        estimator=estimator,
        backward=True,
        discrete=discrete,
        **kwargs,
    )


@d.dedent
@inject_docs(
    __doc__=_find_docs.format(direction=_terminal, fn_name=FinalStatesKey.FORWARD.s)
)
def terminal_states(
    adata: AnnData,
    estimator: type(BaseEstimator) = GPCCA,
    discrete: bool = True,
    **kwargs,
) -> Optional[AnnData]:  # noqa

    return _initial_terminal(
        adata,
        estimator=estimator,
        backward=False,
        discrete=discrete,
        **kwargs,
    )
