from numpy import array, zeros, arange, concatenate, searchsorted, where, unique

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double_or_blank, integer_double_or_blank, blank)


class PCOMPG(object):
    type = 'PCOMPG'
    def __init__(self, model):
        pass

    def write_bdf(self, f, size=8, property_ids=None):
        pass