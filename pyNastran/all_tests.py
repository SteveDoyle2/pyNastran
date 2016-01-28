#import warnings
#warnings.filterwarnings('ignore', 'missing __init__.py*')

#bdf
from pyNastran.bdf.test.all_tests import *

#op2
from pyNastran.op2.test.all_tests import *

#f06
#try:
#    from pyNastran.f06.test.all_tests import *
#except ImportError:
#    pass

#op4
from pyNastran.op4.test.op4_test import TestOP4

#utils
from pyNastran.utils.test.all_tests import *

# converters
try:
    from pyNastran.converters.test_formats import *
except ImportError:
    pass

#gui - just tests the imports
#import pyNastran.gui.gui

import pyNastran
pkg_path = pyNastran.__path__[0]
manual_path = os.path.join(pkg_path, '..', 'docs_sphinx', 'manual') # , 'py_to_rst.py'
sys.path.append(manual_path)

from py_to_rst import create_rst_from_python_files
#create_rst_from_python_files()

if __name__ == "__main__":
    import unittest
    unittest.main()
