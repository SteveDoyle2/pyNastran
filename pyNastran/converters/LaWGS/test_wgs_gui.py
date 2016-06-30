import os

from pyNastran.gui.testing_methods import GUIMethods
from pyNastran.converters.LaWGS.wgs_io import LaWGS_IO
import pyNastran

pkg_path = pyNastran.__path__[0]
model_path = os.path.join(pkg_path, 'converters', 'LaWGS')

import unittest

class LaWGS_GUI(LaWGS_IO, GUIMethods):
    def __init__(self):
        GUIMethods.__init__(self)
        LaWGS_IO.__init__(self)


class TestLawgsGUI(unittest.TestCase):

    def test_tmx_geometry(self):
        geometry_filename = os.path.join(model_path, 'tmx1242.wgs')
        dirname = None

        test = LaWGS_GUI()
        #test.load_nastran_geometry(geometry_filename, None)
        test.load_lawgs_geometry(geometry_filename, dirname)

    def test_tmd_geometry(self):
        geometry_filename = os.path.join(model_path, 'tnd6480.wgs')
        dirname = None

        test = LaWGS_GUI()
        test.load_lawgs_geometry(geometry_filename, dirname)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()

