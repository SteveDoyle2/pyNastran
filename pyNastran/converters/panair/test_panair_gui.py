import os

import pyNastran
from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.converters.panair.panair_io import PanairIO
from pyNastran.utils.log import get_logger

pkg_path = pyNastran.__path__[0]
model_path = os.path.join(pkg_path, 'converters', 'panair', 'M100')

import unittest

class PanairGUI(PanairIO, FakeGUIMethods):
    def __init__(self):
        FakeGUIMethods.__init__(self)
        PanairIO.__init__(self)


class TestPanairGUI(unittest.TestCase):

    def test_m100_geometry(self):
        log = get_logger(level='warning')
        geometry_filename = os.path.join(model_path, 'M100.inp')
        #agps_filename = os.path.join(model_path, 'agps')
        #out_filename = os.path.join(model_path, 'panair.out')
        dirname = None

        test = PanairGUI()
        test.log = log
        test.load_panair_geometry(geometry_filename, dirname)

    def test_m100_results(self):
        log = get_logger(level='warning')
        geometry_filename = os.path.join(model_path, 'M100.inp')
        agps_filename = os.path.join(model_path, 'agps')
        out_filename = os.path.join(model_path, 'panair.out')
        dirname = None

        test = PanairGUI()
        test.log = log
        test.load_panair_geometry(geometry_filename, dirname)
        test.load_panair_results(agps_filename, dirname)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()

