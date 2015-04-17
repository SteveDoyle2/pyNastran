#bdf
from pyNastran.bdf.test.all_tests import *

#op2
from pyNastran.op2.test.all_tests import *

#f06
try:
    from pyNastran.f06.test.all_tests import *
except ImportError:
    pass

#op4
from pyNastran.op4.test.op4_test import TestOP4

#utils
from pyNastran.utils.test.all_tests import *

#gui - just tests the imports
#import pyNastran.gui.gui


if __name__ == "__main__":
    import unittest
    unittest.main()
