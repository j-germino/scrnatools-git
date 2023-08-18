"""
Plots the density of expression of genes on an embedding.
From scrnatools package

Created on Mon Jan 10 15:57:46 2022

@author: joe germino (joe.germino@ucsf.edu)
"""
# external package imports
from anndata import AnnData
from typing import List, Tuple
import seaborn as sns
import matplotlib.pyplot as plt
from math import ceil
from scipy.sparse import issparse
from statsmodels.nonparametric.kernel_density import KDEMultivariateConditional, EstimatorSettings
from sklearn.preprocessing import StandardScaler
import numpy as np

# scrnatools package imports
from .._configs import configs
from .._utils import check_path

logger = configs.create_logger(__name__.split('_', 1)[1])

# -------------------------------------------------------function----------------------------------------------------- #


def gene_density_plot(
        adata: AnnData,
        gene_list: List[str],
        data_loc: str = "X",
        thresh: int = 1,
        latent_rep: str = "X_umap",
        est_settings=None,
        cmap: str = "magma",
        s: int = None,
        ncols: int = 3,
        figsize: Tuple[int] = (3, 3),
        title: str = None,
        save_path: str = None,
        dpi: int = 300,
):
    """Plots the density of expression of genes on an embedding.

    Args:
        adata (AnnData): The dataset containing the gene expression and cell data.
        gene_list (List[str]): A list of genes to plot.
        data_loc (str, optional): The location of the expression data to use for density calculations, can be a layer in 'adata.layers' or 'X' to use the data stored in adata.X. Defaults to "X".
        thresh (int, optional): Kernel density threshold. Defaults to 1.
        latent_rep (str, optional): The 2D representation to plot gene expression for each cell on. Defaults to "X_umap".
        est_settings (optional): Custom settings for the kernel density estimator. Defaults to None.
        cmap (str, optional): The pyplot colormap to use for plotting. Defaults to "magma".
        s (int, optional): The size of data points to plot. Defaults to None.
        ncols (int, optional): The number of columns in the figure. Defaults to 3.
        figsize (Tuple[int], optional): The dimensions of each subplot for a single gene. Defaults to (3, 3).
        title (str, optional): The title of the figure. Defaults to None.
        save_path (str, optional): The path to save the figure. Defaults to None.
        dpi (int, optional): The resolution of the saved image. Defaults to 300.

    Raises:
        ValueError: If the 'data_loc' provided is not 'X' or a valid layer in 'adata.layers'
    """

    # set up figure
    if len(gene_list) < ncols:
        ncols = len(gene_list)
    nrows = ceil((len(gene_list)) / ncols)
    fig = plt.figure(figsize=(ncols * figsize[0], nrows * figsize[1]))
    sns.set_theme(context="paper", style="white", )

    # Iterate over genes to plot
    for index, gene in enumerate(gene_list):
        if gene not in adata.var_names:
            logger.info(f"{gene} not in 'adata.var_names', skipping")
        else:
            temp_data = adata[:, gene]
            if data_loc == "X":
                if issparse(temp_data.X):
                    x = temp_data.X.todense()
                else:
                    x = temp_data.X
            else:
                if data_loc in temp_data.layers:
                    if issparse(temp_data.X):
                        x = temp_data.layers[data_loc].todense()
                    else:
                        x = temp_data.layers[data_loc]
                else:
                    raise ValueError(
                        f"{data_loc} not 'X' or a valid layer in 'adata.layers'")
            scaler = StandardScaler()
            x_scaled = scaler.fit_transform(x).reshape(-1)
            density = KDEMultivariateConditional(
                endog=x_scaled,
                exog=adata.obsm[latent_rep],
                dep_type="c",
                indep_type="cc",
                bw="normal_reference",
                defaults=est_settings
            )
            z1 = density.cdf(
                thresh + np.zeros_like(x_scaled),
                adata.obsm[latent_rep]
            )

            cm = plt.get_cmap(cmap)
            s = int(80000 / adata.shape[0]) if s is None else s

            # Plot density of gene expression
            idx = np.argsort(x_scaled)
            plt.subplot(nrows, ncols, index + 1)
            ax = sns.scatterplot(
                x=adata.obsm[latent_rep][idx, 0],
                y=adata.obsm[latent_rep][idx, 1],
                hue=1 - z1[idx],
                palette=cm,
                s=s,
                linewidth=0
            )
            ax.set_aspect("equal")  # Make sure the aspect ratio is square
            ax.set(xticklabels=[], yticklabels=[], xlabel=None,
                   ylabel=None, title=gene)  # Get rid of x and y labels

            norm = plt.Normalize(1 - z1[idx].min(), 1 - z1[idx].max())
            sm = plt.cm.ScalarMappable(cmap=cm, norm=norm)
            sm.set_array([])

            # Remove the legend and add a colorbar
            ax.get_legend().remove()
            color_bar = ax.figure.colorbar(sm)
            color_bar.set_label("Probability of Expression")
    sns.despine(left=True, bottom=True)  # no borders
    fig.suptitle(title)  # figure title
    fig.tight_layout(rect=[0, 0.03, 1, 0.98])

    if save_path is not None:
        if "/" not in save_path:
            save_path = f"./{save_path}"
        check_path(save_path.rsplit("/", 1)[0])
        logger.info(f"Saving figure to {save_path}")
        plt.savefig(save_path, dpi=dpi, facecolor="white",)
    plt.show()
