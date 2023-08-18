"""
Calculates the cosine similarity of cells to provided cell type signatures.
From scrnatools package

Created on Mon Jan 10 15:57:46 2022

@author: joe germino (joe.germino@ucsf.edu)
"""

# external imports
from anndata import AnnData
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.preprocessing import StandardScaler
import pandas as pd

# scrnatools package imports
from .._configs import configs

logger = configs.create_logger(__name__.split('_', 1)[1])


# -------------------------------------------------------function----------------------------------------------------- #


def cell_type_similarity(
        adata: AnnData,
        signatures: pd.DataFrame,
        normalize: str = "cell",
) -> AnnData:
    """Calculates the cosine similarity of cells to provided cell type signatures.

    Args:
        adata (AnnData): The AnnData object containing cells to score, with library size rescaled and log-normalized counts stored in 'adata.raw'.
        signatures (pd.DataFrame): A DataFrame containing transcriptome signatures for each cell type to be tested. Columns are cell types and rows are the log-normalized expression of each gene in that cell type. Can be created from cell type clusters of a scRNAseq dataset using 'create_cell_type_signature'
        normalize (str, optional): How to rescale the cosine similarity scores. Default is 'cell', which calculates a z-score normalized for each cell across all cell types in 'signatures'. Other possible values are 'cell type' which calculates a z-score normalized for each cell type in 'signatures' across all cells and 'none' which returns the raw cosine similarity values. Defaults to "cell".

    Returns:
        AnnData: The AnnData passed in with 'adata' with the cosine similarity scores for each immgen cell type added as a column to 'adata.obs'
    """

    # Raw attribute contains the log-normalized counts
    cell_data = pd.DataFrame(
        adata.raw.X.toarray(),
        columns=adata.raw.var_names,
        index=adata.obs.index
    ).T
    # Normalize per cell expression to sum to 1
    cell_data = cell_data / cell_data.sum(axis=0)
    # Join with gene x cell type signature matrix so that only shared genes are kept
    all_df = cell_data.join(
        signatures,
        how="inner"
    )
    # Normalize per cell expression to sum to 1
    all_df = all_df / all_df.sum(axis=0)

    # Calculate cosine similarity of single cells to each cell type
    sim = cosine_similarity(
        # Single-cell expression data
        all_df[all_df.columns[:cell_data.shape[1]]].values.T,
        # Cell type gene signature expression data
        all_df[all_df.columns[cell_data.shape[1]:]].values.T
    )

    # Create dataframe with per cell similarity scores for each cell type
    similarity = pd.DataFrame(
        sim,
        columns=all_df.columns[cell_data.shape[1]:],
        index=all_df.columns[:cell_data.shape[1]]
    )
    if normalize == "cell type":
        logger.info(
            f"Scaling cosine similarity scores to z-scores by cell type"
        )
        # Scale similarity scores
        scaler = StandardScaler()
        minmax_scale = scaler.fit(similarity)
        x_minmax = minmax_scale.transform(similarity)
        scale_by_row = pd.DataFrame(x_minmax)
    elif normalize == "cell":
        logger.info(f"Scaling cosine similarity scores to z-scores by cell")
        # Scale similarity scores
        scaler = StandardScaler()
        minmax_scale = scaler.fit(similarity.T)
        x_minmax = minmax_scale.transform(similarity.T).T
        scale_by_row = pd.DataFrame(x_minmax)
    else:
        logger.info(f"Leaving similarity scores unscaled")
        scale_by_row = pd.DataFrame(similarity)

    # Add cosine similarity scores to adata obs
    scale_by_row.index = similarity.index
    if normalize != "none":
        scale_by_row.columns = similarity.columns + \
            "_cosine_similarity_" + normalize + "_z_score"
    else:
        scale_by_row.columns = similarity.columns + "_cosine_similarity"
    adata_copy = adata.copy()
    adata_copy.obs = adata_copy.obs.join(scale_by_row)
    return adata_copy
