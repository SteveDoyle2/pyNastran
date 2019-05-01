"""
tests Usm3d
"""
import os
import unittest
from cpylog import get_logger

import pyNastran
#from pyNastran.bdf.bdf import BDF
from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.converters.usm3d.usm3d_io import Usm3dIO


PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, 'converters', 'usm3d', 'box')


class Usm3dGUI(FakeGUIMethods):
    def __init__(self):
        FakeGUIMethods.__init__(self)
        self.model = Usm3dIO(self)
        self.build_fmts(['usm3d'], stop_on_failure=True)


class TestUsm3dGUI(unittest.TestCase):

    def test_usm3d_geometry_01(self):
        """tests the box.cogwg/box.flo turbulent model"""
        log = get_logger(level='error', encoding='utf-8')
        geometry_filename = os.path.join(MODEL_PATH, 'box.cogsg')
        flo_filename = os.path.join(MODEL_PATH, 'box.flo')

        test = Usm3dGUI()
        test.log = log
        #test.on_load_geometry(geometry_filename, geometry_format='usm3d', raise_error=True)
        test.model.load_usm3d_geometry(geometry_filename)
        test.model.load_usm3d_results(flo_filename)
        test.model.on_reload_usm3d()

        test.model.load_usm3d_geometry(geometry_filename)
        test.model.load_usm3d_results(flo_filename)
        test.model.on_reload_usm3d()

if __name__ == '__main__':  # pragma: no cover
    unittest.main()

