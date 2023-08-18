"""
Creates a UMAP/TSNE plot of a gene's expression from a layer using scanpy.
From scrnatools package

Created on Mon Jan 10 15:57:46 2022

@author: joe germino (joe.germino@ucsf.edu)
"""
# external package imports
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import scipy
from typing import List, Tuple
from anndata import AnnData

# scrnatools package imports
from .._configs import configs

logger = configs.create_logger(__name__.split('_', 1)[1])

# -------------------------------------------------------function----------------------------------------------------- #


def gene_embedding(
        adata: AnnData,
        gene_list: List[str],
        layer: str,
        figsize: Tuple[int, int] = (4, 4),
        dpi: int = 80,
        dpi_save: int = 300,
        min_quantile: float = 0.01,
        max_quantile: float = 0.99,
        use_rep="X_umap",
        *args,
        **kwargs,
):
    """Creates a UMAP/TSNE plot of a gene's expression from a layer using scanpy.

    Args:
        adata (AnnData): The dataset containing gene expression data.
        gene_list (List[str]): The list of genes to plot.
        layer (str): The layer containing expression data, or 'X' to use the expression data stored in adata.X.
        figsize (Tuple[int, int], optional): The size of each panel. Defaults to (4, 4).
        dpi (int, optional): The dpi of plots shown inline. Defaults to 80.
        dpi_save (int, optional): The dpi of saved plots. Defaults to 300.
        min_quantile (float, optional):  The quantile of expression for a gene to set the minimum of the colorbar to. Defaults to 0.01.
        max_quantile (float, optional): The quantile of expression for a gene to set the maximum of the colorbar to. Defaults to 0.99.
        use_rep (str, optional): The embedding coordinates in 'adata.obsm' to use. Defaults to "X_umap".

    Raises:
        ValueError: If 'layer' is not a valid key in 'adata.layers'
    """

    # Setup figure
    sc.set_figure_params(figsize=figsize, dpi=dpi,
                         dpi_save=dpi_save, facecolor="white", frameon=False)
    sns.set_context("paper")
    plt.rcParams["axes.grid"] = False
    vmin = []
    vmax = []
    if layer != "X" and layer not in adata.layers.keys():
        raise ValueError(
            f"{layer} is not 'X' or a valid key in 'adata.layers'")
    # Get min/max expression values for the colorbars for each gene
    for gene in gene_list:
        if layer != "X":
            if scipy.sparse.issparse(adata.layers[layer]):
                min = np.quantile(
                    adata[:, gene].layers[layer].todense().tolist(), min_quantile
                )
                max = np.quantile(
                    adata[:, gene].layers[layer].todense().tolist(), max_quantile
                )
                vmin.append(min)
                if max < 1:
                    vmax.append(1)
                else:
                    vmax.append(max)
            else:
                min = np.quantile(adata[:, gene].layers[layer], min_quantile)
                max = np.quantile(adata[:, gene].layers[layer], max_quantile)
                vmin.append(min)
                vmax.append(max)
        else:
            min = np.quantile(adata[:, gene].X, min_quantile)
            max = np.quantile(adata[:, gene].X, max_quantile)
            vmin.append(min)
            vmax.append(max)
    # Plot the data
    adata_copy = adata.copy()
    adata_copy.obsm["X_umap"] = adata.obsm[use_rep]
    if layer == "X":
        sc.pl.umap(
            adata_copy,
            color=gene_list,
            vmin=vmin,
            vmax=vmax,
            *args,
            **kwargs,
        )
    else:
        sc.pl.umap(
            adata_copy,
            color=gene_list,
            layer=layer,
            vmin=vmin,
            vmax=vmax,
            *args,
            **kwargs,
        )
