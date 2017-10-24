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


class STL_GUI(STL_IO, FakeGUIMethods):
    def __init__(self):
        FakeGUIMethods.__init__(self)
        STL_IO.__init__(self)


class STL_GUITest(unittest.TestCase):

    def test_stl_geometry(self):
        from pyNastran.utils.log import get_logger
        log = get_logger(level='warning')
        geometry_filename = os.path.join(MODEL_PATH, 'sphere.stl')

        test = STL_GUI()
        test.log = log
        #test.load_nastran_geometry(geometry_filename)
        test.load_stl_geometry(geometry_filename)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
