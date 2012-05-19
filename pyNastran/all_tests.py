#bdf
from pyNastran.bdf.test.test_fieldWriter import TestFieldWriter
from pyNastran.bdf.test.bdf_unitTests import BDF_Test
from pyNastran.bdf.test.unit.test_coords import TestCoords

#op2
from pyNastran.op2.test.op2_unitTests import OP2_Test

#f06
#op4

import unittest
class AllTests(TestFieldWriter,BDF_Test,OP2_Test,TestCoords):
    pass

if __name__ == "__main__":
    unittest.main()
