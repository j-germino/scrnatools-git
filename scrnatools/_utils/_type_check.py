"""
Checks the type of a variable.
From scrnatools package

Created on Mon Jan 10 15:57:46 2022

@author: joe germino (joe.germino@ucsf.edu)
"""
from typing import Any

# -------------------------------------------------------function----------------------------------------------------- #


def type_check(var: Any, varname: str, types: Any,):
    """Checks the type of a variable.

    Args:
        var (Any): the variable to check the type of.
        varname (str): the name of the variable.
        types (Any): the type the variable should be.

    Raises:
        TypeError: when 'var' is not one of 'types'.
    """

    if not isinstance(var, types):
        raise TypeError(f"{varname} must be of type {types}")
