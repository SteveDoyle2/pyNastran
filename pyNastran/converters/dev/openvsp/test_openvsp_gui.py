import os
import warnings
import unittest

import numpy as np

from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.converters.dev.openvsp.adb_io import ADB_IO
from pyNastran.converters.dev.openvsp.degen_geom_io import DegenGeomIO
from pyNastran.converters.dev.openvsp.degen_geom import read_degen_geom
from cpylog import get_logger
import pyNastran

warnings.simplefilter('always')
np.seterr(all='raise')

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, 'converters', 'dev', 'openvsp')


class OpenVSP_GUI(DegenGeomIO, ADB_IO, FakeGUIMethods):
    def __init__(self):
        FakeGUIMethods.__init__(self)
        DegenGeomIO.__init__(self, self)
        ADB_IO.__init__(self, self)
        self.build_fmts(['degen_geom'], stop_on_failure=True)


class TestOpenVSP_GUI(unittest.TestCase):

    def test_openvsp_geometry_1(self):
        log = get_logger(level='warning')
        geometry_filename = os.path.join(MODEL_PATH, 'wing_fuse.csv')

        test = OpenVSP_GUI()
        test.log = log
        #test.load_nastran_geometry(geometry_filename)
        #test.load_degen_geom_geometry(geometry_filename)
        test.on_load_geometry(geometry_filename, geometry_format='degen_geom', stop_on_failure=True)

    def test_openvsp_geometry_2(self):
        log = get_logger(level='warning')
        degen_geom_csv = os.path.join(MODEL_PATH, 'wing_fuse.csv')
        bdf_filename = os.path.join(MODEL_PATH, 'wing_fuse.bdf')
        panair_filename = os.path.join(MODEL_PATH, 'wing_fuse.inp')
        model = read_degen_geom(degen_geom_csv, log=log, debug=False)
        model.write_bdf(bdf_filename)

        #panair_filename = 'panair.inp'
        panair_case_filename = 'model_DegenGeom.vspaero'
        model.write_panair(panair_filename, panair_case_filename)
        os.remove(bdf_filename)
        os.remove(panair_filename)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
