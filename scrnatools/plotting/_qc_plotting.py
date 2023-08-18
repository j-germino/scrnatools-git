"""
Creates common preprocessing QC plots for an AnnData object.
From scrnatools package

Created on Mon Jan 10 15:57:46 2022

@author: joe germino (joe.germino@ucsf.edu)
"""
# external package imports
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from anndata import AnnData
from typing import Tuple

# scrnatools package imports
from .._configs import configs
from .._utils import check_path

logger = configs.create_logger(__name__.split('_', 1)[1])

# -------------------------------------------------------function----------------------------------------------------- #


def qc_plotting(
        adata: AnnData,
        counts_thresholds: Tuple[int, int] = (1000, 30000),
        genes_thresholds: Tuple[int, int] = (100, 5000),
        mt_threshold: int = 10,
        show_thresholds: bool = True,
        batch_key: str = None,
        show_legend: bool = True,
        figsize: Tuple[int, int] = (9, 3),
        dpi: int = 300,
        save_path: str = None,
):
    """Creates common preprocessing QC plots for an AnnData object.

    Args:
        adata (AnnData): The dataset containing data to plot
        counts_thresholds (Tuple[int, int], optional): The lower and upper thresholds to be used on total counts when filtering cells. Defaults to (1000, 30000).
        genes_thresholds (Tuple[int, int], optional): The lower and upper thresholds to be used on number of genes when filtering genes. Defaults to (100, 5000).
        mt_threshold (int, optional): The threshold to be used on % mito reads when filtering cells. Defaults to 10.
        show_thresholds (bool, optional): Whether to show the thresholds as dashed lines on each plot. Defaults to True.
        batch_key (str, optional):  A column name in 'adata.obs' that annotates different batches of data for separate plotting. If 'None' treats all cells as coming from the same batch. Defaults to None.
        show_legend (bool, optional):  Whether to show the 'batch_key' labels as a legend. Defaults to True.
        figsize (Tuple[int, int], optional): The size of the figure. Defaults to (9, 3).
        dpi (int, optional): The resolution of the figure to save. Defaults to 300.
        save_path (str, optional): The path to save the figure to (/path/to/dir/filename). Defaults to None.

    Raises:
        ValueError: If batch_key is not a valid key in 'adata.obs.columns'
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
            raise ValueError(
                f"{batch_key} is not a valid column in 'adata.obs'"
            )

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
        category_data["rank"] = category_data.total_counts.rank(
            method="first", ascending=False,)
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
