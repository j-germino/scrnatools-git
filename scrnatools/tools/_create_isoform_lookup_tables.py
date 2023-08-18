"""
Creates the lookup tables for isoform data transcripts and genes.
From scrnatools package

Created on Mon Jan 10 15:57:46 2022

@author: joe germino (joe.germino@ucsf.edu)
"""

# external imports
from anndata import AnnData
from pandas import DataFrame
import numpy as np
import pandas as pd
from typing import Tuple, Dict

# scrnatools package imports
from .._configs import configs

logger = configs.create_logger(__name__.split('_', 1)[1])


# -------------------------------------------------------function----------------------------------------------------- #


def create_isoform_lookup_tables(
        adata: AnnData,
        t2enst_path: DataFrame,
        t2g_path: DataFrame,
) -> Tuple[Dict[str, str], Dict[str, str], Dict[str, str]]:
    """Creates the lookup tables for isoform data transcripts and genes.

    Args:
        adata (AnnData): The AnnData containing kallisto isoform data.
        t2enst_path (DataFrame): Path to the transcript to ensembl transcript id mapping dataframe (from kallisto alignment).
        t2g_path (DataFrame): Path to the transcript to gene mapping dataframe (from the kallisto reference used for alignment).

    Returns:
        Tuple[Dict[str, str], Dict[str, str], Dict[str, str]]: The equivalence class to transcript dict, the equivalence class to gene dict, and the gene to equivalence class dict.
    """

    # Import the transcript list from kallisto
    t2enst = pd.read_csv(t2enst_path, header=None,)
    t2enst.columns = ['enst']

    # Import the transcript to gene mapping dataframe
    t2g = pd.read_csv(t2g_path, header=None, sep='\t')
    t2g.columns = ['enst', 'ensg', 'gname']
    t2g.index = t2g.enst

    # Get ecs from adata
    ecs = np.array(adata.var_names.values, dtype=str)
    # Convert ec string to a list of tx indices
    ecdict = {}
    for ec, txs in enumerate(ecs):
        ecdict[str(ec)] = list(np.array(txs.split(','), dtype=int))
    # Save ec list to adata
    adata.var['ecs'] = ecs
    # Set adata var_names to be new ec index
    # var names are set names
    adata.var_names = np.array(np.arange(adata.shape[1]), dtype=str)
    # Create dictionary mapping tx indices to tx ids for each ec
    ec2tx = {}
    for (ec, txs) in ecdict.items():
        ec2tx[ec] = list(t2enst.loc[txs].enst.values)
    # Create mapping of each ec to the genes whose transcripts are contained in it
    ec2g = {}
    key_errors = {}
    ec_error = []
    for (ec, txs) in ec2tx.items():
        try:
            # get unique gene names for each ec list of transcripts
            ec2g[ec] = list(np.unique(t2g.loc[txs].gname.values))
            ec_error.append(False)
        except KeyError:
            # if there's a tx that doesn't have a t2g mapping
            valid_txs = [i for i in txs if i in t2g.index]
            if len(valid_txs) > 0:
                # if there is at least one valid gene for the ec txs use them instead
                ec2g[ec] = list(np.unique(t2g.loc[valid_txs].gname.values))
            else:
                # Otherwise, fill in none for ec2g mapping
                ec2g[ec] = "None"
            ec_error.append(True)
            # If the ec list of transcripts contains a tx not in the tx to gene mapping dict, note the ec and txs that
            # are the problem
            errors = [i for i in txs if i not in t2g.index]
            key_errors[ec] = errors
    adata.var["ec_tx_error"] = ec_error
    # Create an inverse mapping dictionary getting all ecs and their txs for a given gene
    inv_map = {}
    for ec, genes in ec2g.items():
        for gene in genes:
            # Get the current dict of ec's for that gene
            inv_map[gene] = inv_map.get(gene, {})
            # Append the current ec
            inv_map[gene][ec] = ec2tx[ec]
    # save the dicts
    logger.info(f"Number of ec errors: {len(key_errors)}")
    adata.uns['key_errors'] = key_errors
    return ec2tx, ec2g, inv_map
