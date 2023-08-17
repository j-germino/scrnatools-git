"""
Decorator function to handle debugging
From scrnatools package

Created on Mon Jan 10 15:57:46 2022

@author: joe germino (joe.germino@ucsf.edu)
"""
# external package imports
from logging import Logger
from time import perf_counter
from typing import Any

# -------------------------------------------------------function----------------------------------------------------- #


def debug(logger: Logger, configs,) -> Any:
    """A function decorator for creating various debugging logs"""

    def decorator(func):
        def wrapper(*args, **kwargs):
            # log the function call
            logger.debug(f"Function call of {func.__name__}")  # with args: {args}, and kwargs: {kwargs}")
            start = perf_counter()
            result = func(*args, **kwargs)
            end = perf_counter()
            # if tracking time to execute
            if configs.debug_timing:
                minutes = int((end - start) / 60)
                seconds = (end - start - minutes * 60)
                # log time to execute as XX min X.XXXX s
                logger.info(f"{func.__name__} took: {minutes:d} min {seconds:.4f} s")
            return result

        return wrapper

    return decorator
