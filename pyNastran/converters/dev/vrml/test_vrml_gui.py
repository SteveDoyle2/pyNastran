import os
import unittest
from cpylog import get_logger

import pyNastran
from pyNastran.gui.testing_methods import FakeGUIMethods
#from pyNastran.bdf.bdf import BDF
from pyNastran.converters.dev.vrml.vrml_io import Vrml_io
#from pyNastran.converters.dev.vrml.vrml import read_vrml
PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, 'converters', 'dev', 'vrml')


PKG_PATH = pyNastran.__path__[0]


class VrmlGui(FakeGUIMethods):
    def __init__(self):
        FakeGUIMethods.__init__(self)
        self.model = Vrml_io(self)


class TestVrmlGui(unittest.TestCase):

    def test_vrml_io_1(self):
        """tests geometry"""
        vrml_filename = os.path.join(MODEL_PATH, 'pyramid_sphere.wrl')

        test = VrmlGui()
        test.log = get_logger(level='warning', encoding='utf-8')
        test.model.load_vrml_geometry(vrml_filename)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
