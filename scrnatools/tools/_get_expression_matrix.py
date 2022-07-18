"""
Extracts the cell x gene expression matrix from an AnnData object
From scrnatools package

Created on Mon Jan 10 15:57:46 2022

@author: joe germino (joe.germino@ucsf.edu)
"""
# external package imports
from typing import Optional
from scipy.sparse import issparse
from pandas import DataFrame
from anndata import AnnData

# scrnatools package imports
from .._configs import configs
from .._utils import debug

logger = configs.create_logger(__name__.split('_', 1)[1])

# -------------------------------------------------------function----------------------------------------------------- #


@debug(logger, configs)
def get_expression_matrix(
        anndata: AnnData,
        gene_data: Optional[str] = None,
) -> DataFrame:
    """
    Extracts the cell x gene expression matrix from an AnnData object

    Parameters
    ----------
    anndata
        The AnnData to pull the expression matrix from
    gene_data
        Where to get the expression data from the AnnData object. If a layer from the AnnData['layers'] is passed,
        that is used, otherwise 'X' or 'raw' can be used.

    Returns
    -------
    A Pandas DataFrame containing the expression matrix (cells x genes)

    Raises
    -------
    ValueError
        When 'layer' is not one of the layers in 'anndata'
    """
    if gene_data == "X":
        # Convert to a dense matrix if needed
        if issparse(anndata.X):
            matrix = DataFrame(anndata.X.todense(), index=anndata.obs.index, columns=anndata.var_names)
        else:
            matrix = DataFrame(anndata.X, index=anndata.obs.index, columns=anndata.var_names)
    # get a layer from the AnnData if specified
    elif gene_data == "raw":
        raw_data = anndata.raw.to_adata()
        # Convert to a dense matrix if needed
        if issparse(raw_data.X):
            matrix = DataFrame(raw_data.X.todense(), index=anndata.obs.index, columns=anndata.raw.var_names)
        else:
            matrix = DataFrame(raw_data.X, index=anndata.obs.index, columns=anndata.raw.var_names)
    elif gene_data in anndata.layers:
        if issparse(anndata.layers[gene_data]):
            matrix = DataFrame(anndata.layers[gene_data].todense(), index=anndata.obs.index, columns=anndata.var_names)
        else:
            matrix = DataFrame(anndata.layers[gene_data])
    else:
        raise ValueError(f"{gene_data} is not 'X', 'raw', or a valid layer name in '{anndata.layers}'")
    return matrix
