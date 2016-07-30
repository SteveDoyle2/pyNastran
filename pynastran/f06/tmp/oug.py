#pylint: disable=C0301
from __future__ import print_function
from six.moves import zip

from pyNastran.op2.tables.oug.oug_displacements import RealDisplacement, ComplexDisplacement
from pyNastran.op2.tables.oug.oug_eigenvectors import (Eigenvector, ComplexEigenvector,
                                                       RealEigenvectorArray, ComplexEigenvectorArray)
from pyNastran.op2.tables.oug.oug_temperatures import RealTemperature


class OUG(object):
    def __init__(self):
        self.iSubcases = []
        self.i = 0
