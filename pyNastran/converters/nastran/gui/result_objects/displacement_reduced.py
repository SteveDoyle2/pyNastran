import numpy as np
from pyNastran.op2.result_objects.table_object import (
    RealTableArray, ComplexTableArray)

class DisplacementReduced:
    def __init__(self, case: RealTableArray | ComplexTableArray,
                 nodal_disp: np.ndarray,
                 node_gridtype: np.ndarray):
        self.node_gridtype = node_gridtype
        self.data = nodal_disp
