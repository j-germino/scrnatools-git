"""
Checks the type of a variable
From scrnatools package

Created on Mon Jan 10 15:57:46 2022

@author: joe germino (joe.germino@ucsf.edu)
"""
from typing import Any

# -------------------------------------------------------function----------------------------------------------------- #


def type_check(var: Any, varname: str, types: Any,):
    """
    Checks the type of a variable

    Parameters
    ----------
    var
        the variable to check the type of
    varname
        the name of the variable
    types
        the type the variable should be

    Raises
    -------
    TypeError
        when 'var' is not one of 'types'
    """
    if not isinstance(var, types):
        raise TypeError(f"{varname} must be of type {types}")
