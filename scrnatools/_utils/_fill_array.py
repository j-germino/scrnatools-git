"""
Fills the default value for parameter that don't set values for all plots (keys).
From scrnatools package

Created on Mon Jan 10 15:57:46 2022

@author: joe germino (joe.germino@ucsf.edu)
"""
# external package imports
from typing import Any, List

# -------------------------------------------------------function----------------------------------------------------- #


def fill_array(variable: Any, default: Any, final_size: int) -> List[Any]:
    """Fills the default value for parameter that don't set values for all plots (keys).

    Args:
        variable (Any): The parameter to fill values for.
        default (Any): The default value for the argument on a single plot.
        final_size (int): The number of keys in the figure.

    Returns:
        List[Any]: A list of length 'final_size' with the default values for the parameter filling to the end of the list.
    """

    # ensure the parameter is a list
    if not isinstance(variable, list):
        variable = [variable]
    # get the current length of the parameter (i.e. the number of keys with user-specified values)
    num_vars = len(variable)
    if num_vars < final_size:
        for _ in range(final_size - num_vars):
            if num_vars == 1:
                # if only one value is provided for the parameter, use that by default for all keys
                variable.append(variable[0])
            else:
                # otherwise, add the default value for the parameter to fill the list
                variable.append(default)
    return variable
