"""with pyQt5/pySide2 and vtk"""
#import os
import sys
from pyNastran.gui.menus.test.test_about import *
# if 'XVFB' in os.environ or sys.platform == 'win32':   # XVFB is for TravisCI and doesn't work

if sys.platform == 'win32' and len(sys.argv) == 1:
    from pyNastran.gui.menus.test.test_groups import *

if __name__ == "__main__":  # pragma: no cover
    import unittest
    unittest.main()
