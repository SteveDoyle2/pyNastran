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

import pyNastran.gui.formats
import pyNastran.gui.gui_common
#if not IS_DEV:
    #from pyNastran.gui.menus.test_groups import *
from pyNastran.converters.test_gui_formats import *
from pyNastran.gui.test.all_tests_gui import *

#if IS_H5PY:
    #from pyNastran.dev.h5.test.test_h5 import *


def run_tests_cleanup():
    import os
    import unittest
    with np.errstate(divide='raise', over='raise', under='raise', invalid='raise'):
        #standard - has CI core dump issue
        #unittest.main()

        # 1. Run unittest without letting it invoke sys.exit()
        res = unittest.main(exit=False)

        # 2. Force coverage.py to save data before we kill the process
        try:
            import coverage
            cov = coverage.Coverage.current()
            if cov is not None:
                cov.stop()
                cov.save()
        except ImportError:
            pass  # Safe fallback if coverage isn't running locally

        # 3. Bypass VTK's crashing destructors using a hard exit
        if res.result.wasSuccessful():
            # Bypass VTK's crashing destructors and exit clean
            os._exit(0)
        else:
            # Properly fail the CI pipeline if an actual test failed
            os._exit(1)

if __name__ == "__main__":  # pragma: no cover
    run_tests_cleanup()
