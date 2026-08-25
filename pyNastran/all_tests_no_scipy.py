# no pytables, no pandas, no vtk, no hdf5

import os
import sys

import pyNastran

pkg_path = pyNastran.__path__[0]

#bdf
#from pyNastran.bdf.test.all_tests import *

#op2
from pyNastran.op2.test.test_op2_no_scipy import *

#f06
#from pyNastran.f06.test.all_tests import *

#op4
from pyNastran.op4.test.op4_unit_tests import TestOP4 #, TestOP4Fast

#utils
from pyNastran.utils.test.all_tests import *
from pyNastran.femutils.test.all_tests import *


if __name__ == "__main__":  # pragma: no cover
    import unittest
    unittest.main()
