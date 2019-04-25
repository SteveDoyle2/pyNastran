"""tests the pyNastranGUI PANAIR (A502) interface"""
import os
import unittest
from cpylog import get_logger

import pyNastran
from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.converters.panair.panair_io import PanairIO
from pyNastran.converters.panair.panair_out import read_panair_out

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, 'converters', 'panair')


class PanairGUI(FakeGUIMethods):
    def __init__(self):
        FakeGUIMethods.__init__(self)
        self.model = PanairIO(self)
        self.build_fmts(['panair'], stop_on_failure=True)


class TestPanairGUI(unittest.TestCase):

    def test_m100_geometry(self):
        """simple panair geometry test"""
        log = get_logger(level='warning')
        geometry_filename = os.path.join(MODEL_PATH, 'M100', 'M100.inp')
        #agps_filename = os.path.join(MODEL_PATH, 'agps')
        #out_filename = os.path.join(MODEL_PATH, 'panair.out')

        test = PanairGUI()
        test.log = log
        #test.model.load_panair_geometry(geometry_filename)
        test.on_load_geometry(geometry_filename, geometry_format='panair', raise_error=True)

    def test_m100_results(self):
        """contains one AGPS results"""
        log = get_logger(level='warning')
        geometry_filename = os.path.join(MODEL_PATH, 'M100', 'M100.inp')
        agps_filename = os.path.join(MODEL_PATH, 'M100', 'agps')
        out_filename = os.path.join(MODEL_PATH, 'M100', 'panair.out')
        read_panair_out(out_filename)

        test = PanairGUI()
        test.log = log
        test.model.load_panair_geometry(geometry_filename)
        test.model.load_panair_results(agps_filename)

    def test_naca0012_results(self):
        """contains multiple AGPS results"""
        log = get_logger(level='warning')
        geometry_filename = os.path.join(MODEL_PATH, 'naca0012', 'multi', 'naca0012_ar9.inp')
        agps_filename = os.path.join(MODEL_PATH, 'naca0012', 'multi', 'agps_0012')
        out_filename = os.path.join(MODEL_PATH, 'naca0012', 'multi', 'panair.out')

        test = PanairGUI()
        test.log = log
        test.model.load_panair_geometry(geometry_filename)
        test.model.load_panair_results(agps_filename)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
