"""
Performs basic preprocessing of kallisto isoform data.
From scrnatools package

Created on Mon Jan 10 15:57:46 2022

@author: joe germino (joe.germino@ucsf.edu)
"""

# external imports
from anndata import AnnData
import scanpy as sc
from pandas import DataFrame
import numpy as np

# scrnatools package imports
from .._configs import configs

logger = configs.create_logger(__name__.split('_', 1)[1])


# -------------------------------------------------------function----------------------------------------------------- #


def isoform_preprocessing(
        adata: AnnData,
        ec: DataFrame,
):
    """Performs basic preprocessing of kallisto isoform data.

    Args:
        adata (AnnData): The AnnData containing kallisto isoform data.
        ec (DataFrame): The equivalence class matrix to rename vars with.
    """

    adata.obs_names_make_unique()
    # Filter out equivalence classes with no counts
    ngenes = adata.shape[1]
    sc.pp.filter_genes(adata, min_counts=1)
    logger.info(f"Filtered {ngenes - adata.shape[1]} ecs with no counts")
    # Filter cells with no ecs
    ncells = adata.shape[0]
    sc.pp.filter_cells(adata, min_genes=1)
    logger.info(f"Filtered {ncells - adata.shape[0]} cells with no ecs")
    # rename vars to be ecs
    ec = ec.iloc[np.array(adata.var_names, dtype=int)]
    adata.var_names = ec.loc[:, 1].values
