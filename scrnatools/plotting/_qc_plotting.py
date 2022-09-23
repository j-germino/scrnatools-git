"""
Creates common preprocessing QC plots for an AnnData object
From scrnatools package

Created on Mon Jan 10 15:57:46 2022

@author: joe germino (joe.germino@ucsf.edu)
"""
# external package imports
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from anndata import AnnData
from typing import Optional, Tuple

# scrnatools package imports
from .._configs import configs
from .._utils import debug, check_path

logger = configs.create_logger(__name__.split('_', 1)[1])

# -------------------------------------------------------function----------------------------------------------------- #


@debug(logger, configs)
def qc_plotting(
        adata: AnnData,
        counts_thresholds: Tuple[int, int] = [1000, 30000],
        genes_thresholds: Tuple[int, int] = [100, 5000],
        mt_threshold: int = 10,
        show_thresholds: bool = True,
        batch_key: Optional[str] = None,
        show_legend: bool = True,
        figsize: Tuple[int, int] = (9, 3),
        dpi: int = 300,
        save_path: Optional[str] = None,
):
    """
    Creates common preprocessing QC plots for an AnnData object
    Parameters
    ----------
    adata
        The AnnData containing data to plot
    counts_thresholds
        The lower and upper thresholds to be used on total counts when filtering cells. Default (1000, 30000)
    genes_thresholds
        The lower and upper thresholds to be used on number of genes when filtering cells . Default (100, 5000)
    mt_threshold
        The threshold to be used on % mito reads when filtering cells. Default 10
    show_thresholds
        Whether to show the thresholds as dashed lines on each plot. Default True
    batch_key
        A column name in 'adata.obs' that annotates different batches of data for separate plotting. If 'None' treats
        all cells as coming from the same batch. Default None
    show_legend
        Whether to show the 'batch_key" labels as a legend. Default True
    figsize
        The size of the figure. Default (9, 3)
    dpi
        The resolution of the figure to save. Default 300
    save_path
        The path to save the figure to (/path/to/dir/filename). Default None

    Raises
    -------
    ValueError
        If batch_key is not a valid key in 'adata.obs.columns'
    """
    # Setup plots
    sns.set_style("ticks")
    sns.set_context("paper")
    fig = plt.figure(figsize=figsize)
    if batch_key is None:
        adata.obs["qc_plot_batch"] = "None"
        batch_key = "qc_plot_batch"
        show_legend = False
    else:
        if batch_key not in adata.obs.columns:
            raise ValueError(f"{batch_key} is not a valid column in 'adata.obs'")

    # Plot 1: Histogram of % mito counts
    plt.subplot(1, 3, 1)
    ax = sns.histplot(
        x=adata.obs.pct_counts_mt,
        hue=adata.obs[batch_key],
    )
    if not show_legend:
        ax.get_legend().remove()
    if show_thresholds:
        plt.axvline(
            x=mt_threshold,
            ymin=0,
            ymax=1,
            color="black",
            linestyle="--",
        )

    # Plot 2: total counts vs num genes
    plt.subplot(1, 3, 2)
    ax = sns.scatterplot(
        x=adata.obs.total_counts,
        y=adata.obs.n_genes_by_counts,
        s=2,
        hue=adata.obs[batch_key],
        linewidth=0,
    )
    if not show_legend:
        ax.get_legend().remove()
    if show_thresholds:
        for threshold in counts_thresholds:
            plt.axvline(
                x=threshold,
                ymin=0,
                ymax=1,
                color="black",
                linestyle="--",
            )
        for threshold in genes_thresholds:
            plt.axhline(
                y=threshold,
                xmin=0,
                xmax=1,
                color="black",
                linestyle="--",
            )

    # Plot 3: Rank ordered total counts
    plt.subplot(1, 3, 3)
    cell_data = pd.DataFrame()
    for category in adata.obs[batch_key].unique():
        category_data = adata[adata.obs[batch_key] == category].obs.copy()
        category_data["rank"] = category_data.total_counts.rank(method="first", ascending=False,)
        cell_data = pd.concat([cell_data, category_data])
    cell_data = cell_data.sort_values(by="rank")
    ax = sns.lineplot(
        x=cell_data["rank"],
        y=cell_data.total_counts,
        hue=cell_data[batch_key],
        hue_order=adata.obs[batch_key].unique(),
    )
    if not show_legend:
        ax.get_legend().remove()
    ax.set(yscale="log")
    ax.set(xscale="log")
    if show_thresholds:
        for threshold in counts_thresholds:
            plt.axhline(
                y=threshold,
                xmin=0,
                xmax=1,
                color="black",
                linestyle="--",
            )

    # Figure level adjustments
    fig.tight_layout()
    if save_path is not None:
        if "/" not in save_path:
            save_path = f"./{save_path}"
        check_path(save_path.rsplit("/", 1)[0])
        logger.info(f"Saving figure to {save_path}")
        plt.savefig(save_path, dpi=dpi, facecolor="white",)
    plt.show()
