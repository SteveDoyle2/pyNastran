#import warnings
#warnings.filterwarnings('ignore', 'missing __init__.py*')

import sys
import os
import pyNastran
pkg_path = pyNastran.__path__[0]

#bdf
from pyNastran.bdf.test.all_tests_no_matplotlib import *
from pyNastran.dev.bdf_vectorized2.bdfv_unit_tests import *

#op2
from pyNastran.op2.test.all_tests import *

#f06
#from pyNastran.f06.test.all_tests import *

#op4
from pyNastran.op4.test.op4_unit_tests import TestOP4

#utils
from pyNastran.utils.test.all_tests import *
from pyNastran.femutils.test.all_tests import *

# converters
from pyNastran.converters.test_formats import *

#gui - just tests the imports
from pyNastran.gui.test.all_tests import *


if __name__ == "__main__":  # pragma: no cover
    import unittest
    unittest.main()
