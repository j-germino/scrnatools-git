"""
Plots qc metrics for kallisto isoform data.
From scrnatools package

Created on Mon Jan 10 15:57:46 2022

@author: joe germino (joe.germino@ucsf.edu)
"""

# external imports
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import numpy as np
from anndata import AnnData
from matplotlib.patches import Patch

# scrnatools package imports
from scrnatools._configs import configs

logger = configs.create_logger(__name__.split('_', 1)[1])

# -------------------------------------------------------function----------------------------------------------------- #


def isoform_qc_plot(
        adata: AnnData,
        batch_key: str = None,
):
    """Plots qc metrics for kallisto isoform data.

    Args:
        adata (AnnData): The dataset containing isoform data to plot.
        batch_key (str, optional): An obs column identifying batches to plot separately. Defaults to None.

    Raises:
        ValueError: If the batch key provided is not a column in in 'adata.obs'
    """

    if batch_key not in adata.obs.columns:
        raise ValueError(f"{batch_key} is not a valid column in 'adata.obs'")

    sc.set_figure_params(figsize=(6, 6), dpi=80, dpi_save=300,
                         facecolor="white", frameon=False)
    sns.set_context("notebook")
    plt.rcParams["axes.grid"] = False

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22',
              '#17becf']
    if len(adata.obs[batch_key].unique()) > len(colors):
        logger.info(
            f"Too many batches in adata.obs.{batch_key}, colors will overlap")
    if batch_key is None:
        adata.obs["qc_plot_batch"] = "all_cells"
        batch_key = "qc_plot_batch"

    # Counts per ec plot
    fig = plt.figure(figsize=(12, 6))
    fig.add_subplot(1, 2, 1)
    plt.semilogy(np.sort(adata.var['total_counts'].values)[::-1])
    plt.ylabel("total_counts")
    plt.xlabel("ec")

    # Knee plot (all ec counts per barcode)
    fig.add_subplot(1, 2, 2)
    legend_elements = []
    for i, batch in enumerate(adata.obs[batch_key].unique()):
        q = np.sort(
            np.array(adata[adata.obs[batch_key] == batch].X.sum(axis=1)).reshape(-1)
        )[::-1]
        plt.loglog(q, linewidth=1)
        plt.ylabel("UMI Counts")
        plt.xlabel("Barcodes")
        legend_elements.append(Patch(facecolor=colors[i], label=batch))

    plt.legend(handles=legend_elements)
    plt.show()
