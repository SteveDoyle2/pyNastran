import os
import warnings
import unittest

import numpy as np

from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.converters.dev.openvsp.adb_io import ADB_IO
from pyNastran.converters.dev.openvsp.degen_geom_io import DegenGeomIO
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


class TestOpenVSP_GUI(unittest.TestCase):

    def test_openvsp_geometry(self):
        from pyNastran.utils.log import get_logger
        log = get_logger(level='warning')
        geometry_filename = os.path.join(MODEL_PATH, 'wing_fuse.csv')

        test = OpenVSP_GUI()
        test.log = log
        #test.load_nastran_geometry(geometry_filename)
        test.load_degen_geom_geometry(geometry_filename)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
