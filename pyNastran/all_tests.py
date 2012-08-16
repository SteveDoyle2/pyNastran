import unittest

#bdf
from pyNastran.bdf.test.test_fieldWriter import TestFieldWriter
from pyNastran.bdf.test.bdf_unitTests import BDF_Test
from pyNastran.bdf.test.unit.test_coords import TestCoords

#op2
from pyNastran.op2.test.op2_unitTests import OP2_Test

#f06
from pyNastran.f06.test.f06_test import main as F06

#op4
from pyNastran.op4.test.op4_test import OP4_Test

#gui - just tests the imports
#import pyNastran.gui.gui


if __name__ == "__main__":
    unittest.main()
    #F06()
