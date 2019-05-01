import os
import unittest

from cpylog import get_logger
import pyNastran
from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.converters.shabp.shabp_io import ShabpIO

PKG_PATH = pyNastran.__path__[0]
model_path = os.path.join(PKG_PATH, 'converters', 'shabp')


class ShabpGUI(ShabpIO, FakeGUIMethods):
    def __init__(self):
        FakeGUIMethods.__init__(self)
        self.model = ShabpIO(self)
        self.build_fmts(['shabp'], stop_on_failure=True)


class TestShabpGUI(unittest.TestCase):

    def test_shabp_geometry_01(self):
        log = get_logger(level='warning')
        test = ShabpGUI()
        test.log = log
        shabp_infilename = os.path.join(model_path, 'models', 'flap', 'flap_inviscid.mk5')
        shabp_outfilename = os.path.join(model_path, 'models', 'flap', 'SHABP.OUT')
        #test.model.load_shabp_geometry(shabp_infilename)
        test.on_load_geometry(shabp_infilename, geometry_format='shabp', raise_error=True)
        test.on_load_results(shabp_outfilename)

    def _test_shabp_geometry_02(self):
        test = ShabpGUI()
        shabp_infilename = os.path.join(model_path, 'models', 'orbiter.mk5')
        test.on_load_geometry(shabp_infilename, geometry_format='shabp', raise_error=True)

    def _test_shabp_geometry_03(self):
        test = ShabpGUI()
        shabp_infilename = os.path.join(model_path, 'models', 'shuttle.mk5')
        test.on_load_geometry(shabp_infilename, geometry_format='shabp', raise_error=True)

    def test_shabp_geometry_04(self):
        test = ShabpGUI()
        shabp_infilename = os.path.join(model_path, 'models', 'nose', 'noseX_working.mk5')
        test.on_load_geometry(shabp_infilename, geometry_format='shabp', raise_error=True)


    #def test_shabp_results(self):
        #pass
        #geometry_filename = os.path.join(model_path, 'M100.inp')
        #agps_filename = os.path.join(model_path, 'agps')
        #out_filename = os.path.join(model_path, 'panair.out')

        #test = ShabpGUI()
        #test.model.load_panair_geometry(geometry_filename)
        #test.model.load_panair_results(agps_filename)

        #model = ShabpOut()
        #model.read_shabp_out('SHABP.OUT')


if __name__ == '__main__':  # pragma: no cover
    unittest.main()

