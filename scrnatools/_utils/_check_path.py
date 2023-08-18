"""
Checks if a path exists and creates it if it doesn't.
From scrnatools package

Created on Mon Jan 10 15:57:46 2022

@author: joe germino (joe.germino@ucsf.edu)
"""
# external package imports
import os

# ----------------------function---------------------- #


def check_path(path: str):
    """Checks if a path exists and creates it if it doesn't.

    Args:
        path (str): The string of the path to check/create.
    """

    if not os.path.exists(path):
        os.makedirs(path)
