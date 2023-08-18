"""
Saves the lookup tables for an isoform AnnData object.
From scrnatools package

Created on Mon Jan 10 15:57:46 2022

@author: joe germino (joe.germino@ucsf.edu)
"""

# external imports
import json
from typing import Dict

# scrnatools package imports
from .._configs import configs

logger = configs.create_logger(__name__.split('_', 1)[1])


# -------------------------------------------------------function----------------------------------------------------- #


def save_isoform_lookup_tables(
        path: str,
        ec2tx: Dict[str, str],
        ec2g: Dict[str, str],
        inv_map: Dict[str, str],
):
    """Saves the lookup tables for an isoform AnnData object.

    Args:
        path (str): The path and file name to save the lookup dicts to. (_ec2tx.json, _ec2g.json, and _inv_map.json all appended to this file name).
        ec2tx (Dict[str, str]): The equivalence class to transcript dict to save.
        ec2g (Dict[str, str]): The equivalence class to gene dict to save.
        inv_map (Dict[str, str]): The gene to equivalence class dict to save.
    """

    with open(f"{path}_ec2tx.json", "w") as f:
        json_file = json.dumps(ec2tx)
        f.write(json_file)
    with open(f"{path}_ec2g.json", "w") as f:
        json_file = json.dumps(ec2g)
        f.write(json_file)
    with open(f"{path}_inv_map.json", "w") as f:
        json_file = json.dumps(inv_map)
        f.write(json_file)
