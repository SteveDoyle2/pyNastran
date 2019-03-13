import os
import warnings
import unittest

import numpy as np

from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.converters.stl.stl_io import STL_IO
import pyNastran

warnings.simplefilter('always')
np.seterr(all='raise')

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, 'converters', 'stl')


class STL_GUI(FakeGUIMethods):
    def __init__(self):
        FakeGUIMethods.__init__(self)
        self.model = STL_IO(self)
        self.build_fmts(['stl'], stop_on_failure=True)


class STL_GUITest(unittest.TestCase):

    def test_stl_geometry(self):
        from cpylog import get_logger
        log = get_logger(level='warning')
        geometry_filename = os.path.join(MODEL_PATH, 'sphere.stl')

        test = STL_GUI()
        test.log = log
        #test.model.load_nastran_geometry(geometry_filename)
        #test.model.load_stl_geometry(geometry_filename)
        test.on_load_geometry(geometry_filename, geometry_format='stl', raise_error=True)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
