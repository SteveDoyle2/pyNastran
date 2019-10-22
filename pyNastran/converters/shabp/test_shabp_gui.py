import os
import unittest

from cpylog import get_logger
import pyNastran
from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.converters.shabp.shabp_io import ShabpIO
from pyNastran.converters.shabp.shabp import read_shabp

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, 'converters', 'shabp', 'models')


class ShabpGUI(ShabpIO, FakeGUIMethods):
    def __init__(self):
        FakeGUIMethods.__init__(self)
        self.model = ShabpIO(self)
        self.build_fmts(['shabp'], stop_on_failure=True)


class TestShabpGUI(unittest.TestCase):

    def test_shabp_results_1(self):
        log = get_logger(level='warning')
        test = ShabpGUI()
        test.log = log
        shabp_infilename = os.path.join(MODEL_PATH, 'flap', 'flap_inviscid.mk5')
        shabp_outfilename = os.path.join(MODEL_PATH, 'flap', 'SHABP.OUT')

        #test.model.load_shabp_geometry(shabp_infilename)
        test.on_load_geometry(shabp_infilename, geometry_format='shabp', raise_error=True)
        unused_model = read_shabp(shabp_infilename, log=None, debug=None)
        #model.get_area_by_patch()
        #model.get_area_by_component()
        #model.get_area_xlength_by_component()
        test.on_load_results(shabp_outfilename)

    def _test_shabp_geometry_1(self):
        test = ShabpGUI()
        shabp_infilename = os.path.join(MODEL_PATH, 'orbiter.mk5')
        test.on_load_geometry(shabp_infilename, geometry_format='shabp', raise_error=True)

    def _test_shabp_geometry_2(self):
        test = ShabpGUI()
        shabp_infilename = os.path.join(MODEL_PATH, 'shuttle.mk5')
        test.on_load_geometry(shabp_infilename, geometry_format='shabp', raise_error=True)

    def test_shabp_geometry_3(self):
        test = ShabpGUI()
        shabp_infilename = os.path.join(MODEL_PATH, 'nose', 'noseX_working.mk5')
        test.on_load_geometry(shabp_infilename, geometry_format='shabp', raise_error=True)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
