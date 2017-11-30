import os
import unittest

from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.converters.lawgs.wgs_io import LaWGS_IO
import pyNastran

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, 'converters', 'lawgs')


class LaWGS_GUI(LaWGS_IO, FakeGUIMethods):
    def __init__(self):
        FakeGUIMethods.__init__(self)
        LaWGS_IO.__init__(self)


class TestLawgsGUI(unittest.TestCase):

    def test_tmx_geometry(self):
        geometry_filename = os.path.join(MODEL_PATH, 'tmx1242.wgs')
        test = LaWGS_GUI()
        test.load_lawgs_geometry(geometry_filename)

    def test_tmd_geometry(self):
        geometry_filename = os.path.join(MODEL_PATH, 'tnd6480.wgs')
        test = LaWGS_GUI()
        test.load_lawgs_geometry(geometry_filename)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
