"""
Plots a heatmap of expression of genes.
From scrnatools package

Created on Tue March 28 17:47:46 2023

@author: joe germino (joe.germino@ucsf.edu), nolan horner (nolan.horner@ucsf.edu)
"""
# external package imports
from anndata import AnnData
from typing import List
import seaborn as sns
import matplotlib.pyplot as plt
import scrnatools

# scrnatools package imports
from .._configs import configs
from .._utils import check_path

logger = configs.create_logger(__name__.split('_', 1)[1])

# -------------------------------------------------------function----------------------------------------------------- #


def gene_heatmap(
        adata: AnnData,
        gene_list: List[str],
        obs_key: str,
        layer: str = "X",
        obs_values: List[str] = ['All'],
        cbar_args: List[int] = None,
        swap_axes: bool = False,
        cell_size: int = 30,
        save_path: str = None,
        dpi: int = 300,
        *args, **kwargs
):
    """Plots a heatmap of expression of genes.

    Args:
        adata (AnnData): The dataset containing the gene expression and cell data.
        gene_list (List[str]):  A list of genes to plot.
        obs_key (str): The categorical grouping to display on the heatmap.
        layer (str, optional): The type of the expression data to use for the heatmap, can be a layer in 'adata.layers' or 'X' to use the data stored in adata.X. Defaults to "X".
        obs_values (List[str], optional): Values from obs key group to display on heatmap. Defaults to ['All'].
        cbar_args (List[int], optional): List of integers to position color bar on heatmap. [x position, y position, width, height]. 'None' automatically adjusts to right middle of heatmap.. Defaults to None.
        swap_axes (bool, optional): AWhether to swap x and y axes. Defaults to False.
        cell_size (int, optional): Integer specifying the size of each cell in the heatmap. Defaults to 30.
        save_path (str, optional): The path to save the figure. Defaults to None.
        dpi (int, optional): The resolution of the saved image. Defaults to 300.

    Raises:
        ValueError: If a gene in 'gene_list' provided is not in provided AnnData layer.
        ValueError: If 'gene_list' provided less than 2.
        ValueError: If the 'obs_values' provided less than 2.
        ValueError: If a value in 'obs_values' provided is not in provided AnnData layer obs key.
    """

    invalid_genes = [i for i in gene_list if i not in adata.var_names]
    if len(invalid_genes) > 0:
        raise ValueError(f"Genes not found in data: {invalid_genes}")
    if len(gene_list) < 2:
        raise ValueError(
            f"Please choose two or more genes. {len(gene_list)} genes were entered.")
    if layer == "raw":
        raw_adata = adata.raw.to_adata()
        expression_matrix = scrnatools.tl.get_expression_matrix(
            raw_adata[:, gene_list],
            gene_data="X"
        )
    else:
        expression_matrix = scrnatools.tl.get_expression_matrix(
            adata[:, gene_list],
            gene_data=layer
        )
    expression_matrix[obs_key] = adata.obs[obs_key]
    # matrix with mean expression levels of selected genes for each cell type
    expression_matrix = expression_matrix.groupby(obs_key).mean()
    if obs_values != ['All']:  # include only cell types from cell_type_list
        if len(obs_values) < 2:
            raise ValueError(
                f"Please choose more than one obs value. Or use '[All]' to see all {obs_key}s")
        for i in obs_values:
            if i not in expression_matrix.index:
                raise ValueError(
                    f"{i} is not in adata's obs_key {obs_key}.\nPossible {obs_key}s: {list(expression_matrix.index.categories)}"
                )
        for j in expression_matrix.index:
            if j not in obs_values:
                expression_matrix = expression_matrix.drop(j)

    if swap_axes:
        expression_matrix = expression_matrix.transpose()

    figdpi = plt.rcParams['figure.dpi']
    totalWidth = plt.rcParams['figure.subplot.right'] - \
        plt.rcParams['figure.subplot.left']
    totalHeight = plt.rcParams['figure.subplot.top'] - \
        plt.rcParams['figure.subplot.bottom']
    nrows, ncols = expression_matrix.shape
    figWidth = (ncols*cell_size/figdpi)/totalWidth
    figHeight = (nrows*cell_size/figdpi)/totalHeight

    sns.set_theme(context="paper", style="white", )
    cg = sns.clustermap(
        expression_matrix, figsize=(figWidth, figHeight),
        *args, **kwargs
    )

    # calculate the size of the heatmap axes
    axWidth = (ncols*cell_size)/(figWidth*figdpi)
    axHeight = (nrows*cell_size)/(figHeight*figdpi)

    # resize heatmap
    ax_heatmap_og_pos = cg.ax_heatmap.get_position()
    cg.ax_heatmap.set_position(
        [ax_heatmap_og_pos.x0, ax_heatmap_og_pos.y0, axWidth, axHeight]
    )

    # layout of colorbar
    if cbar_args is None:
        cbar_x_pos = 0.85 + ax_heatmap_og_pos.x0
        cbar_y_pos = ax_heatmap_og_pos.y0
        cbar_width = 25/(figWidth*figdpi)
        cbar_height = axHeight
        cbar_args = [cbar_x_pos, cbar_y_pos, cbar_width, cbar_height]

    # set the layout
    cg.ax_cbar.set_position(cbar_args)
    cg.ax_cbar.grid(False)
    cg.ax_row_dendrogram.set_visible(False)
    cg.ax_col_dendrogram.set_visible(False)
    cg.ax_heatmap.yaxis.tick_left()
    cg.ax_heatmap.set_ylabel("")
    cg.ax_heatmap.set_xlabel("")
    cg.ax_heatmap.grid(False)
    cg.ax_heatmap.tick_params(axis='y', labelrotation=0)
    cg.ax_heatmap.tick_params(axis='x', labelrotation=90)
    cg.ax_heatmap.tick_params(bottom=True, left=True)

    if save_path is not None:
        if "/" not in save_path:
            save_path = f"./{save_path}"
        check_path(save_path.rsplit("/", 1)[0])
        logger.info(f"Saving figure to {save_path}")
        plt.savefig(save_path, dpi=dpi, facecolor="white", bbox_inches="tight")
    plt.show()
