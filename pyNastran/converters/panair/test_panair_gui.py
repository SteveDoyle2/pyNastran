import os

from pyNastran.gui.testing_methods import GUIMethods
from pyNastran.converters.panair.panair_io import PanairIO
import pyNastran

pkg_path = pyNastran.__path__[0]
model_path = os.path.join(pkg_path, 'converters', 'panair', 'M100')

import unittest

class PanairGUI(PanairIO, GUIMethods):
    def __init__(self):
        GUIMethods.__init__(self)
        PanairIO.__init__(self)


class TestPanairGUI(unittest.TestCase):

    def test_m100_geometry(self):
        geometry_filename = os.path.join(model_path, 'M100.inp')
        #agps_filename = os.path.join(model_path, 'agps')
        #out_filename = os.path.join(model_path, 'panair.out')
        dirname = None

        test = PanairGUI()
        #test.load_nastran_geometry(geometry_filename, None)
        test.load_panair_geometry(geometry_filename, dirname)

    def test_m100_results(self):
        geometry_filename = os.path.join(model_path, 'M100.inp')
        agps_filename = os.path.join(model_path, 'agps')
        out_filename = os.path.join(model_path, 'panair.out')
        dirname = None

        test = PanairGUI()
        test.load_panair_geometry(geometry_filename, dirname)
        test.load_panair_results(agps_filename, dirname)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()

