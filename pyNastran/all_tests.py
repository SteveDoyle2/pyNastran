import numpy as np
#import warnings
#warnings.filterwarnings('ignore', 'missing __init__.py*')
from pyNastran.gui.qt_version import qt_version
if qt_version == 'pyqt4':
    import PyQt4
elif qt_version == 'pyqt5':
    import PyQt5
elif qt_version == 'pyside':
    import PySide
else:
    raise NotImplementedError(qt_version)

import pyNastran.gui.formats
import pyNastran.gui.gui_common
from pyNastran.all_tests_no_gui import *
from pyNastran.converters.test_gui_formats import *


if __name__ == "__main__":  # pragma: no cover
    import unittest
    with np.errstate(divide='raise', over='raise', under='raise', invalid='raise'):
        unittest.main()
