import os
import warnings
import unittest

import numpy as np
from cpylog import get_logger
warnings.simplefilter('always')
np.seterr(all='raise')


import pyNastran
from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.converters.su2.su2_io import SU2_IO

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, 'converters', 'su2')


class SU2_GUI(FakeGUIMethods):
    def __init__(self):
        FakeGUIMethods.__init__(self)
        self.model = SU2_IO(self)
        self.build_fmts(['su2'], stop_on_failure=True)


class TestSU2GUI(unittest.TestCase):

    def test_su2_geometry(self):
        """tests mesh_naca0012_inv.su2"""
        log = get_logger(level='warning')
        geometry_filename = os.path.join(MODEL_PATH, 'mesh_naca0012_inv.su2')

        test = SU2_GUI()
        test.log = log
        #test.model.load_su2_geometry(geometry_filename)
        test.on_load_geometry(geometry_filename, geometry_format='su2', raise_error=True)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()

