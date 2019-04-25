"""
Various numpy dependent functions for interfacing with nan data are defined in this file.
This includes:
 - is_array = isfinite(array)
 - is_array = isfinite_and_greater_than(array, value)
 - is_array = isfinite_and_nonzero(array)
"""
from __future__ import print_function, division
import numpy as np

def isfinite(_array):
    """are any of the values finite?"""
    return _array is not None and np.any(np.isfinite(_array))

def isfinite_and_greater_than(_array, value):
    """are any of the values finite and greater than some value?"""
    return isfinite(_array) and abs(np.nanmax(_array) > value)

def isfinite_and_nonzero(_array):
    """are any of the values finite and a value non-zero?"""
    return isfinite_and_greater_than(_array, 0.)
