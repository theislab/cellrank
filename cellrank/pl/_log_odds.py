from typing import Tuple, Union, Optional, Sequence
from pathlib import Path

import scanpy as sc
from anndata import AnnData
from cellrank.ul._docs import d
from cellrank.tl._utils import save_fig

import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt


@d.dedent
def log_odds(
    adata: AnnData,
    fate: str,
    gene: Optional[Union[str, Sequence[str]]] = None,
    threshold: Optional[float] = None,  # TODO
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (10, 4),
    dpi: Optional[int] = None,
    save: Optional[Union[str, Path]] = None,
) -> None:
    """
    TODO Plot log-odds.

    Log-odds are plotted as a function of the experimental time-point. Optionally, gene
    expression can be overlaid to check for genes that show early fate bias.

    Parameters
    ----------
    %(adata)s
    fate
        TODO.
    gene
        TODO.
        Key from `adata.var_names`.
    title
        Title of the plot.
    %(plotting)s

    Returns
    -------
    %(just_plots)s
    """

    # define log-odds
    fate1 = adata.obsm["to_terminal_states"][fate].X
    fate2 = 1 - fate1
    log_odds = np.log(
        1e-9 + np.divide(fate1, fate2, out=np.zeros_like(fate1), where=fate2 != 0)
    )
    adata.obs["log_odds"] = log_odds

    # define keys to be pulled from adata
    keys = ["day", "log_odds"]

    # optionally, add gene expression
    if gene is not None:
        # threshold gene expression
        adata.obs[f"{gene}+"] = (adata[:, gene].X.A > threshold).flatten()
        keys.append(f"{gene}+")

    # get correponding dataframe
    plotting_df = sc.get.obs_df(adata, keys=keys)

    # filter out cells where either fate is zero
    zero_mask = np.logical_or(fate1 == 0, fate2 == 0)
    plotting_df = plotting_df[~zero_mask]

    # prepare for plotting
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    # use seaborn here for some jiiter
    if gene is not None:
        ax = sns.stripplot(
            "day",
            "log_odds",
            data=plotting_df[~plotting_df[f"{gene}+"]],
            jitter=0.4,
            color="k",
            ax=ax,
            s=0.5,
        )
        ax = sns.stripplot(
            "day",
            "log_odds",
            data=plotting_df[plotting_df[f"{gene}+"]],
            jitter=0.4,
            color="r",
            ax=ax,
            s=2,
        )
    else:
        ax = sns.stripplot(
            "day", "log_odds", data=plotting_df, jitter=0.4, color="k", ax=ax, s=2
        )

    ticks = np.arange(0, len(plotting_df["day"].cat.categories), 4)
    ax.set_xticks(ticks=ticks)

    # TODO
    ax.set_xlabel("Day")
    ax.set_ylabel("Log Odds")
    ax.set_title(title)

    if save is not None:
        save_fig(fig, save)
