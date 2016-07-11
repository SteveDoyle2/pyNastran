#import warnings
#warnings.filterwarnings('ignore', 'missing __init__.py*')

import sys
import os
import pyNastran
pkg_path = pyNastran.__path__[0]
manual_path = os.path.abspath(os.path.join(pkg_path, '..', 'docs', 'html_docs', 'manual')) # , 'py_to_rst.py'
sys.path.append(manual_path)

#print(sys.path)
from notebook_to_markdown import create_rst_from_ipython_notebooks
#create_rst_from_ipython_notebooks()

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
from pyNastran.gui.test.all_tests import *
#on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
#if not on_rtd:
#    import pyNastran.gui.gui
#import pyNastran.gui.gui


if __name__ == "__main__":
    import unittest
    unittest.main()
