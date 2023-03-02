"""
Calculates a differential density ratio of cells grouped by condition/control
From scrnatools package

Created on Mon Jan 10 15:57:46 2022

@author: joe germino (joe.germino@ucsf.edu)
"""

# external package imports
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from statsmodels.nonparametric.kernel_density import KDEMultivariateConditional,EstimatorSettings
from anndata import AnnData
from typing import Optional, Tuple

# scrnatools package imports
from .._configs import configs
from .._utils import debug

logger = configs.create_logger(__name__.split('_', 1)[1])

# -------------------------------------------------------function----------------------------------------------------- #


@debug(logger, configs)
def log_density_ratio(
        adata: AnnData,
        group_by: str,
        groups: Tuple[str, str],
        latent_rep: str = "X_umap",
        key_added: str = "log_density_ratio",
        est_settings: Optional = None,
):
    """
    Calculates a differential density ratio of cells grouped by condition/control
    Parameters
    ----------
    adata
        The AnnData containing labeled cell data
    group_by
        The column in 'adata.obs' containing condition and control labels
    groups
        The labels in 'adata.obs[group_by]' that distinguish condition and control cells, i.e. ["WT", "KO"]
    latent_rep
        The latent representation to calculate differential density across
    key_added
        The column in 'adata.obs' that will the density ratios will be added to
    est_settings

    Returns
    -------
        The original adata modified to contain the differential density ratios in a new column in 'adata.obs'

    Raises
    -------
    ValueError
        If 'group_by' is not a valid column name in 'adata.obs'
    ValueError
        If 'groups' are not valid labels in 'adata.obs[group_by]'
    """
    if group_by not in adata.obs:
        raise ValueError(f"{group_by} is not a valid column key in 'adata.obs'")
    for group in groups:
        if group not in adata.obs[group_by].unique():
            raise ValueError(
                f"{group} is not a valid category of 'adata.obs.{group_by}' ({adata.obs[group_by].unique()})"
            )
    density = KDEMultivariateConditional(
        endog=adata.obsm[latent_rep],
        exog=np.array(adata.obs[group_by].values == groups[1]),
        dep_type="cc",
        indep_type="u",
        bw="normal_reference",
        defaults=est_settings
    )
    logger.info(f"Calculating log density ratio of {groups[1]}/{groups[0]}")
    z1 = density.pdf(adata.obsm[latent_rep], 0 + np.zeros((adata.obsm[latent_rep].shape[0], 1)))
    z2 = density.pdf(adata.obsm[latent_rep], 1 + np.zeros((adata.obsm[latent_rep].shape[0], 1)))
    adata.obs[key_added] = np.log(z2/z1)  # group2 vs group1
    return adata
