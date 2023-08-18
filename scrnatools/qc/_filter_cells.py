"""
Filters cells based on gene number, total counts, and % mitochondrial reads.
From scrnatools package

Created on Mon Jan 10 15:57:46 2022

@author: joe germino (joe.germino@ucsf.edu)
"""

# external imports
from anndata import AnnData
from typing import Tuple

# scrnatools package imports
from .._configs import configs

logger = configs.create_logger(__name__.split('_', 1)[1])

# -------------------------------------------------------function----------------------------------------------------- #


def filter_cells(
        adata: AnnData,
        genes_thresholds: Tuple[int, int],
        counts_thresholds: Tuple[int, int],
        mt_threshold: int = 10,
) -> AnnData:
    """Filters cells based on gene number, total counts, and % mitochondrial reads.

    Args:
        adata (AnnData): The AnnData with the data to filter.
        genes_thresholds (Tuple[int, int]): A Tuple of thresholds for the number of genes per cell with 'gene_thresholds[0]' being the lower bound and 'gene_thresholds[1]' being the upper bound (both exclusive).
        counts_thresholds (Tuple[int, int]): A Tuple of thresholds for the number of total counts per cell with 'count_thresholds[0]' being the lower bound and 'count_thresholds[1]' being the upper bound (both exclusive).
        mt_threshold (int, optional): The maximum percent mitochondrial reads per cell. Defaults to 10.

    Returns:
        AnnData:  The dataset filtered on cells that don't pass the thresholds.
    """

    logger.info(f"Number of cells before QC filtering: {len(adata.obs)}")
    filtered_adata = adata[adata.obs.pct_counts_mt < mt_threshold].copy()
    filtered_adata = filtered_adata[
        filtered_adata.obs.total_counts < counts_thresholds[1]
    ]
    filtered_adata = filtered_adata[
        filtered_adata.obs.total_counts > counts_thresholds[0]
    ]
    filtered_adata = filtered_adata[
        filtered_adata.obs.n_genes_by_counts < genes_thresholds[1]
    ]
    filtered_adata = filtered_adata[
        filtered_adata.obs.n_genes_by_counts > genes_thresholds[0]
    ].copy()
    logger.info(
        f"Number of cells after QC filtering: {len(filtered_adata.obs)}")
    return filtered_adata
