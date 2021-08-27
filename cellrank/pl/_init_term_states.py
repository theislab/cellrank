"""Module used for finding initial and terminal states."""
from typing import Union, Optional, Sequence
from typing_extensions import Literal

from anndata import AnnData
from cellrank.ul._docs import d, _initial, _terminal, inject_docs
from cellrank.tl.estimators import GPCCA
from cellrank.tl.kernels._precomputed_kernel import DummyKernel

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


# TODO: refactor me
def _initial_terminal(
    adata: AnnData,
    backward: bool = False,
    discrete: bool = False,
    states: Optional[Union[str, Sequence[str]]] = None,
    cluster_key: Optional[str] = None,
    mode: Literal["embedding", "time"] = "embedding",
    time_key: str = "latent_time",
    **kwargs,
) -> None:

    pk = DummyKernel(adata=adata, backward=backward)
    mc = GPCCA(pk, read_from_adata=True, write_to_adata=False)

    if mc.terminal_states is None:
        raise RuntimeError(
            f"Compute {_initial if backward else _terminal} states first as "
            f"`cellrank.tl.{'initial' if backward else 'terminal'}_states()`."
        )

    n_states = len(mc.terminal_states.cat.categories)
    if n_states == 1 or (
        states is not None and (isinstance(states, str) or len(states) == 1)
    ):
        kwargs["same_plot"] = True

    if kwargs.get("title", None) is None:
        if discrete:
            if kwargs.get("same_plot", True):
                kwargs["title"] = "initial states" if backward else "terminal states"
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
            kwargs["title"] = "initial states" if backward else "terminal states"

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
        fn_name="initial",
        title="initial states",
    )
)
def initial_states(  # noqa: D103
    adata: AnnData,
    discrete: bool = False,
    states: Optional[Union[str, Sequence[str]]] = None,
    cluster_key: Optional[str] = None,
    mode: Literal["embedding", "time"] = "embedding",
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
        fn_name="terminal",
        title="terminal states",
    )
)
def terminal_states(  # noqa: D103
    adata: AnnData,
    discrete: bool = False,
    states: Optional[Union[str, Sequence[str]]] = None,
    cluster_key: Optional[str] = None,
    mode: Literal["embedding", "time"] = "embedding",
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
