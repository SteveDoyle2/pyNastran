import os
import unittest

import pyNastran
from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.converters.lawgs.wgs_reader import read_lawgs
from pyNastran.converters.lawgs.wgs_io import LaWGS_IO

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, 'converters', 'lawgs')


class LaWGS_GUI(LaWGS_IO, FakeGUIMethods):
    def __init__(self):
        FakeGUIMethods.__init__(self)
        self.model = LaWGS_IO(self)
        self.build_fmts(['lawgs'], stop_on_failure=True)


class TestLawgsGUI(unittest.TestCase):
    """tests LAWGS"""

    def test_tmx_geometry(self):
        """tests tmx1242.wgs"""
        geometry_filename = os.path.join(MODEL_PATH, 'tmx1242.wgs')
        test = LaWGS_GUI()
        test.model.load_lawgs_geometry(geometry_filename)
        test.on_load_geometry(geometry_filename, geometry_format='lawgs', raise_error=True)

    def test_tmd_geometry(self):
        """tests tnd6480.wgs"""
        geometry_filename = os.path.join(MODEL_PATH, 'tnd6480.wgs')
        test = LaWGS_GUI()
        test.model.load_lawgs_geometry(geometry_filename)

    def test_lawgs_to_plot3d(self):
        """tests write_as_plot3d"""
        lawgs_filename = os.path.join(MODEL_PATH, 'tnd6480.wgs')
        plot3d_filename = os.path.join(MODEL_PATH, 'tnd6480.p3d')
        model = read_lawgs(lawgs_filename, log=None, debug=False)
        model.write_as_plot3d(plot3d_filename)
        #os.remove(plot3d_filename)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
