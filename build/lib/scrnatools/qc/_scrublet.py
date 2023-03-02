"""
Filters out doublets using scrublet
From scrnatools package

Created on Mon Jan 10 15:57:46 2022

@author: joe germino (joe.germino@ucsf.edu)
"""
# external package imports
from anndata import AnnData
from typing import Optional
import scrublet as scr
import pandas as pd

# scrnatools package imports
from .._configs import configs
from .._utils import debug

logger = configs.create_logger(__name__.split('_', 1)[1])

# -------------------------------------------------------function----------------------------------------------------- #


@debug(logger, configs)
def scrublet(
        adata: AnnData,
        raw_counts_layer: str,
        doublet_threshold: float = 0.2,
        batch_key: Optional[str] = None,
) -> AnnData:
    """
    Filters out doublets using scrublet

    Parameters
    ----------
    adata
        The AnnData to process
    doublet_threshold
        The doublet score threshold to call doublets/singlets on. Default 0.2
    batch_key
        A column in 'adata.obs' annotating different batches of data. If 'None' treats all cells as coming from the same
        batch. Default None
    raw_counts_layer
        A key in 'adata.layers' pointing to the raw counts data

    Returns
    -------
    The AnnData provided with doublets filtered out inplace based on doublet_threshold
    """
    # Setup
    scrublet_predictions = pd.DataFrame()
    if batch_key is None:
        adata.obs["scrublet_batch_key"] = "None"
        batch_key = "scrublet_batch_key"
    # Calculate doublet scores for each batch indepenently
    for key in adata.obs[batch_key].unique():
        subset_data = adata[adata.obs[batch_key] == key]
        counts_matrix = subset_data.layers[raw_counts_layer]
        scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(
            min_counts=2,
            min_cells=3,
            min_gene_variability_pctl=85,
            n_prin_comps=30
        )
        scrub.call_doublets(threshold=doublet_threshold)
        scrub.plot_histogram()
        scrublet_prediction = pd.DataFrame({"scrublet_score": doublet_scores}, index=subset_data.obs.index)
        scrublet_predictions = pd.concat([scrublet_predictions, scrublet_prediction])
    # Add scrublet scores to adata.obs and filter out cells with a score > doublet_threshold
    adata.obs["scrublet_score"] = scrublet_predictions.scrublet_score
    adata.obs["scrublet_called_doublet"] = adata.obs.scrublet_score > doublet_threshold
    num_doublets = adata.obs.scrublet_called_doublet.sum()
    pct_doublets = round(num_doublets / len(adata) * 100, 3)
    logger.info(f"{pct_doublets}% of cells classified as doublets ({num_doublets} cells)")
    return adata[~adata.obs.scrublet_called_doublet].copy()
