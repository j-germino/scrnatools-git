"""
Creates a csv with the immgen cell populations signature data.
From scrnatools package

Created on Mon Jan 10 15:57:46 2022

@author: joe germino (joe.germino@ucsf.edu)
"""

# external imports
import pandas as pd
import subprocess
import os
from shutil import which

# scrnatools package imports
from .._configs import configs
from .._utils import check_path

logger = configs.create_logger(__name__.split('_', 1)[1])


# -------------------------------------------------------function----------------------------------------------------- #


def get_immgen_similarity_signatures(
        save_path: str = "datasets"
) -> pd.DataFrame:
    """Creates a csv with the immgen cell populations signature data.

    Args:
        save_path (str, optional): The path to save the immgen data and signature file to. Defaults to "datasets".

    Raises:
        OSError: if wget is not installed on the system

    Returns:
        pd.DataFrame: A DataFrame containing the immgen cell type expression signatures, with cell types as columns and genes as rows
    """

    if save_path[-1] == "/":
        save_path = save_path[:-1]
    # Download data from github to datasets folder
    if which("wget") is None:
        raise OSError("wget is required to download files. Please install using 'sudo apt-get install wget' on UNIX or "
                      "download and install from 'http://gnuwin32.sourceforge.net/packages/wget.htm' on Windows")
    if not os.path.exists(f"{save_path}/immgen"):
        check_path(save_path)
        subprocess.run(["wget", "-P", f"{save_path}/",
                        "https://gist.github.com/vasilisNt/5e23eeefc188e1e772f428c74ef43277/raw/67f83d282b0b2180a8eeff74edf079d8826b12ba/immgen.tar.gz"],
                       stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL,)
        subprocess.run(
            ["tar", "-xzf", f"{save_path}/immgen.tar.gz", "-C", f"{save_path}/"]
        )
        subprocess.run(["rm", "-rf", f"{save_path}/immgen.tar.gz"])

    # Import immgen data
    immgen = pd.read_csv(
        f"{save_path}/immgen/Immgen_expression.txt",
        sep='\t',
        index_col=0
    )
    probes = pd.read_csv(
        f"{save_path}/immgen/Immgen_probes.txt",
        sep='\t',
        index_col=0,
        header=None
    )
    g2n = pd.read_csv(f"{save_path}/immgen/mart_export.txt", sep='\t', )

    # create DataFrame with signatures
    probes = probes.join(g2n.set_index('Gene stable ID')['Gene name'], on=1)
    probes.columns = ['gene_id', 'gene_name']
    immgen = immgen.join(probes['gene_name']).set_index('gene_name')
    immgen = immgen.groupby(immgen.index).aggregate('sum')

    # Save signatures csv
    immgen.to_csv(f"{save_path}/immgen/immgen_signatures.csv")

    return immgen
