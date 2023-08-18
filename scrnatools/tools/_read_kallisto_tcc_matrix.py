"""
Reads in a kallisto bustools isoform matrix.
From scrnatools package

Created on Mon Jan 10 15:57:46 2022

@author: joe germino (joe.germino@ucsf.edu)
"""

# external imports
import scanpy as sc
import pandas as pd
from pandas import DataFrame
from typing import Tuple
from anndata import AnnData
import numpy as np

# scrnatools package imports
from .._configs import configs

logger = configs.create_logger(__name__.split('_', 1)[1])


# -------------------------------------------------------function----------------------------------------------------- #


def read_kallisto_tcc_matrix(
        sample_path: str,
) -> Tuple[AnnData, DataFrame]:
    """Reads in a kallisto bustools isoform matrix.

    Args:
        sample_path (str): The path to and sample name for the mtx, barcodes, and ec files from kallisto bustools output.

    Returns:
        Tuple[AnnData, DataFrame]: An AnnData with the tcc data and the ec matrix.
    """

    logger.info(f"Reading sample: {sample_path}")
    data = sc.read_mtx(f"{sample_path}.mtx")
    data.obs_names = pd.read_csv(
        f"{sample_path}.barcodes.txt",
        header=None
    ).values.T[0]
    ec = pd.read_csv(f"{sample_path}.ec.txt", sep='\t', header=None)
    data.var_names = np.array(np.arange(data.shape[1]), dtype=str)
    return data, ec
