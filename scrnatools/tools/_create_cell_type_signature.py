"""
Creates a pseudobulked cell type gene signature from scRNAseq data.
From scrnatools package

Created on Mon Jan 10 15:57:46 2022

@author: joe germino (joe.germino@ucsf.edu)
"""

# external imports
from anndata import AnnData
from pandas import DataFrame
import numpy as np

# scrnatools package imports
from .._configs import configs
from .._utils import check_path

logger = configs.create_logger(__name__.split('_', 1)[1])


# -------------------------------------------------------function----------------------------------------------------- #


def create_cell_type_signature(
        adata: AnnData,
        save_path: str = None,
        cell_type_labels: str = "cell_type",
        data_loc: str = "raw",
) -> DataFrame:
    """Creates a pseudobulked cell type gene signature from scRNAseq data.

    Args:
        adata (AnnData): The AnnData object containing the cell type expression data, with cell labels as a key in 'adata.obs' and library size corrected and log-normalized counts in adata.X, adata.raw, or a layer in adata.layers.
        save_path (str, optional): The path to save a csv containing the cell type gene signatures to. Defaults to None.
        cell_type_labels (str, optional): The column name in 'adata.obs' containing the cell type labels for each cell. Defaults to "cell_type".
        data_loc (str, optional): The location of the library size corrected and log-normalized gene expression data in 'adata'. Default is 'raw' but can also be 'X' or a valid key from 'adata.layers'. Defaults to "raw".

    Raises:
        ValueError: If 'cell_type_labels' is not a valid column in 'adata.obs'
        ValueError: If 'data_loc' is not 'X', 'raw', or a valid layer in 'adata.layers'

    Returns:
        pd.DataFrame: A DataFrame containing the psuedobulked gene signatures for each cell type in adata.obs[cell_type_labels] with columns containing data each cell type and rows containing the average expression for each gene within that call
    """

    if cell_type_labels not in adata.obs:
        raise ValueError(
            f"{cell_type_labels} not a valid column in 'adata.obs'"
        )
    if data_loc == "raw":
        signatures = DataFrame(index=adata.raw.var_names)
        raw_adata = adata.raw.to_adata()
        raw_adata.obs[cell_type_labels] = adata.obs[cell_type_labels]
    else:
        signatures = DataFrame(index=adata.var_names)
    for cell_type in adata.obs[cell_type_labels].unique():
        if data_loc == "X":
            signatures[cell_type] = adata[adata.obs[cell_type_labels]
                                          == cell_type].X.mean(axis=0)
        elif data_loc == "raw":
            signatures[cell_type] = np.array(
                raw_adata[
                    raw_adata.obs[cell_type_labels] == cell_type
                ].X.todense()
            ).mean(axis=0)
        else:
            if data_loc in adata.layers:
                signatures[cell_type] = adata[
                    adata.obs[cell_type_labels] == cell_type
                ].layers[data_loc].mean(axis=0)
            else:
                raise ValueError(
                    f"{data_loc} not 'X', 'raw', or a valid layer in 'adata.layers'")
    if save_path is not None:
        check_path(save_path.rsplit("/")[0])
        logger.info(f"Saving signature DataFrame to {save_path}")
        signatures.to_csv(save_path)
    return signatures
