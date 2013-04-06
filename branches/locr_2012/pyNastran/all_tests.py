import unittest

#bdf
from pyNastran.bdf.test.test_field_writer import TestFieldWriter
from pyNastran.bdf.test.bdf_unit_tests import TestBDF
from pyNastran.bdf.test.unit.test_coords import TestCoords

#op2
from pyNastran.op2.test.op2_unit_tests import TestOP2

#f06
from pyNastran.f06.test.f06_test import main as F06
from pyNastran.f06.test.f06_unit_tests import TestF06


#op4
from pyNastran.op4.test.op4_test import TestOP4

#gui - just tests the imports
#import pyNastran.gui.gui


if __name__ == "__main__":
    unittest.main()
    #F06()
