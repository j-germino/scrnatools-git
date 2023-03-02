"""
Loads the lookup tables for an isoform AnnData object
From scrnatools package

Created on Mon Jan 10 15:57:46 2022

@author: joe germino (joe.germino@ucsf.edu)
"""

# external imports
import json
from typing import Tuple, Dict

# scrnatools package imports
from .._configs import configs
from .._utils import debug, check_path

logger = configs.create_logger(__name__.split('_', 1)[1])


# -------------------------------------------------------function----------------------------------------------------- #


@debug(logger, configs)
def load_isoform_lookup_tables(path: str) -> Tuple[Dict[str, str], Dict[str, str], Dict[str, str]]:
    """
    Loads the lookup tables for an isoform AnnData object

    Parameters
    ----------
    path
        The path and file name of the lookup dicts to load. (_ec2tx.json, _ec2g.json, and _inv_map.json all appended to
        this file name)

    Returns
    -------
    ec2tx
        The equivalence class to transcript dict
    ec2g
        The equivalence class to gene dict
    inv_map
        The gene to equivalence class dict
    """
    with open(f"{path}_ec2tx.json") as d:
        ec2tx = json.load(d)
    with open(f"{path}_ec2g.json") as d:
        ec2g = json.load(d)
    with open(f"{path}_inv_map.json") as d:
        inv_map = json.load(d)
    return ec2tx, ec2g, inv_map
