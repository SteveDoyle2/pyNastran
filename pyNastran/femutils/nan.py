"""
Various numpy dependent functions for interfacing with nan data are defined in this file.
This includes:
 - is_array = isfinite(array)
 - is_array = isfinite_and_greater_than(array, value)
 - is_array = isfinite_and_nonzero(array)
 - is_array = isgreater_int(_array, value)

"""
import numpy as np

def isfinite(_array: np.ndarray):
    """are any of the values finite?"""
    return _array is not None and np.any(np.isfinite(_array))

def isfinite_and_greater_than(_array: np.ndarray, value: float):
    """are any of the values finite and greater than some value?"""
    return isfinite(_array) and abs(np.nanmax(_array) > value)

def isfinite_and_nonzero(_array: np.ndarray):
    """are any of the values finite and a value non-zero?"""
    return isfinite_and_greater_than(_array, 0.)

def isgreater_int(_array: np.ndarray, value: int):
    """is the max value greater than some integer value"""
    return _array is not None and _array.max() > value
