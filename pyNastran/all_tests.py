import numpy as np

try:
    import matplotlib  # type: ignore
    IS_MATPLOTLIB = True
except ModuleNotFoundError:
    IS_MATPLOTLIB = False

try:
    # missing-imports
    import h5py  # type: ignore
    IS_H5PY = True
except ModuleNotFoundError:
    IS_H5PY = False


if IS_MATPLOTLIB:
    from pyNastran.gui.matplotlib_backend import matplotlib_backend
    matplotlib.use(matplotlib_backend)

#import warnings
#warnings.filterwarnings('ignore', 'missing __init__.py*')
from pyNastran.gui.qt_version import qt_version
#if qt_version == 'pyqt5':
    #import PyQt5  # pylint: disable=unused-import
#elif qt_version == 'pyside2':
    #import PySide2  # pylint: disable=unused-import
#else:  # pragma: no cover
    #raise NotImplementedError(qt_version)

from pyNastran.gui import IS_DEV
import pyNastran.gui.formats
import pyNastran.gui.gui_common
#if not IS_DEV:
    #from pyNastran.gui.menus.test_groups import *
from pyNastran.all_tests_no_gui import *
from pyNastran.converters.test_gui_formats import *
from pyNastran.gui.test.all_tests_gui import *

#if IS_H5PY:
    #from pyNastran.dev.h5.test.test_h5 import *


if __name__ == "__main__":  # pragma: no cover
    import unittest
    with np.errstate(divide='raise', over='raise', under='raise', invalid='raise'):
        unittest.main()
