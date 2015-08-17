import os

from pyNastran.gui.testing_methods import add_dummy_gui_functions
from pyNastran.converters.LaWGS.wgsIO import LaWGS_IO
import pyNastran

pkg_path = pyNastran.__path__[0]
model_path = os.path.join(pkg_path, 'converters', 'LaWGS')

import unittest

class LawgsGUI(unittest.TestCase):

    def test_tmx_geometry(self):
        geometry_filename = os.path.join(model_path, 'tmx1242.wgs')
        dirname = None

        test = LaWGS_IO()
        add_dummy_gui_functions(test)

        #test.load_nastran_geometry(geometry_filename, None)
        test.load_lawgs_geometry(geometry_filename, dirname)

    def test_tmd_geometry(self):
        geometry_filename = os.path.join(model_path, 'tnd6480.wgs')
        dirname = None

        test = LaWGS_IO()
        add_dummy_gui_functions(test)
        test.load_lawgs_geometry(geometry_filename, dirname)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()

