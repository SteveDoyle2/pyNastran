"""
Default and unknown values, for bdf cards in particular.
If a value is blank or None coming from pyNastran, use the defaults.
If the h5Nastran developer doesn't know what the value should be, use the unknowns.
"""
import numpy as np


class DataHelper(object):
    default_double = np.nan
    default_int = -10
    default_str = ''

    unknown_double = -10101010.
    unknown_int = -20
    unknown_str = '?'
