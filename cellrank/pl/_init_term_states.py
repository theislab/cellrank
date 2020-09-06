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
cluster_key
    If given, plot cluster annotations left of the {direction} states.
%(time_mode)s
time_key
    Key from ``adata.obs`` to use as a pseudotime ordering of the cells.
discrete
    If `True`, plot probability distribution of {direction} states.
    Only available when {direction} were estimated by :class:`cellrank.tl.estimators.GPCCA`.
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
    states: Optional[Union[str, Sequence[str]]] = None,
    cluster_key: Optional[str] = None,
    mode: str = "embedding",
    time_key: str = "latent_time",
    **kwargs,
) -> None:

    pk = DummyKernel(adata=adata, backward=backward)
    mc = GPCCA(pk, read_from_adata=True, write_to_adata=False)

    if mc._get(P.FIN) is None:
        raise RuntimeError(
            f"Compute {_initial if backward else _terminal} states first as "
            f"`cellrank.tl.compute_{FinalStatesKey.BACKWARD if backward else FinalStatesKey.FORWARD}()`."
        )

    n_states = len(mc._get(P.FIN).cat.categories)
    if n_states == 1 or (
        states is not None and (isinstance(states, str) or len(states) == 1)
    ):
        kwargs["same_plot"] = True

    if kwargs.get("title", None) is None:
        if discrete:
            if kwargs.get("same_plot", True):
                kwargs["title"] = (
                    FinalStatesPlot.BACKWARD.s
                    if backward
                    else FinalStatesPlot.FORWARD.s
                )
        elif (
            mode == "embedding"
            and kwargs.get("title", None) is None
            and (
                kwargs.get("same_plot", True)
                and n_states > 1
                and (
                    states is None or (not isinstance(states, str) and len(states) > 1)
                )
            )
        ):
            kwargs["title"] = (
                FinalStatesPlot.BACKWARD.s if backward else FinalStatesPlot.FORWARD.s
            )

    _ = kwargs.pop("lineages", None)

    mc.plot_final_states(
        lineages=states,
        cluster_key=cluster_key,
        mode=mode,
        time_key=time_key,
        discrete=discrete,
        **kwargs,
    )


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
    discrete: bool = False,
    states: Optional[Union[str, Sequence[str]]] = None,
    cluster_key: Optional[str] = None,
    mode: str = "embedding",
    time_key: str = "latent_time",
    **kwargs,
) -> Optional[AnnData]:  # noqa

    return _initial_terminal(
        adata,
        backward=True,
        discrete=discrete,
        states=states,
        cluster_key=cluster_key,
        mode=mode,
        time_key=time_key,
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
    discrete: bool = False,
    states: Optional[Union[str, Sequence[str]]] = None,
    cluster_key: Optional[str] = None,
    mode: str = "embedding",
    time_key: str = "latent_time",
    **kwargs,
) -> Optional[AnnData]:  # noqa

    return _initial_terminal(
        adata,
        backward=False,
        discrete=discrete,
        states=states,
        cluster_key=cluster_key,
        mode=mode,
        time_key=time_key,
        **kwargs,
    )
