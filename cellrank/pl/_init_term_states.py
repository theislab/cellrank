"""Module used for finding initial and terminal states."""
from typing import Union, TypeVar, Optional, Sequence

from cellrank.ul._docs import d, _initial, _terminal, inject_docs
from cellrank.tl._constants import TermStatesKey, TerminalStatesPlot
from cellrank.tl.estimators import GPCCA
from cellrank.tl.estimators._constants import P
from cellrank.tl.kernels._precomputed_kernel import DummyKernel

AnnData = TypeVar("AnnData")


_find_docs = """\
Plot {direction} states uncovered by :class:`cellrank.tl.{fn_name}`.

Parameters
----------
%(adata)s
discrete
    If `True`, plot probability distribution of {direction} states.
    Only available when {direction} were estimated by :class:`cellrank.tl.estimators.GPCCA`.
states
    Subset of {direction} states to plot. If `None`, plot all {direction} states.
cluster_key
    If given, plot cluster annotations left of the {direction} states.
%(time_mode)s
time_key
    Key in ``adata.obs`` where the pseudotime is stored.
%(basis)s
kwargs
    Keyword arguments for :meth:`cellrank.tl.estimators.BaseEstimator.plot_terminal_states`.

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

    if mc._get(P.TERM) is None:
        raise RuntimeError(
            f"Compute {_initial if backward else _terminal} states first as "
            f"`cellrank.tl.compute_{TermStatesKey.BACKWARD if backward else TermStatesKey.FORWARD}()`."
        )

    n_states = len(mc._get(P.TERM).cat.categories)
    if n_states == 1 or (
        states is not None and (isinstance(states, str) or len(states) == 1)
    ):
        kwargs["same_plot"] = True

    if kwargs.get("title", None) is None:
        if discrete:
            if kwargs.get("same_plot", True):
                kwargs["title"] = (
                    TerminalStatesPlot.BACKWARD.s
                    if backward
                    else TerminalStatesPlot.FORWARD.s
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
                TerminalStatesPlot.BACKWARD.s
                if backward
                else TerminalStatesPlot.FORWARD.s
            )

    _ = kwargs.pop("lineages", None)

    mc.plot_terminal_states(
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
        fn_name=TermStatesKey.BACKWARD.s,
        title=TerminalStatesPlot.BACKWARD.s,
    )
)
def initial_states(  # noqa: D103
    adata: AnnData,
    discrete: bool = False,
    states: Optional[Union[str, Sequence[str]]] = None,
    cluster_key: Optional[str] = None,
    mode: str = "embedding",
    time_key: str = "latent_time",
    **kwargs,
) -> Optional[AnnData]:

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
        fn_name=TermStatesKey.FORWARD.s,
        title=TerminalStatesPlot.FORWARD.s,
    )
)
def terminal_states(  # noqa: D103
    adata: AnnData,
    discrete: bool = False,
    states: Optional[Union[str, Sequence[str]]] = None,
    cluster_key: Optional[str] = None,
    mode: str = "embedding",
    time_key: str = "latent_time",
    **kwargs,
) -> Optional[AnnData]:

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
