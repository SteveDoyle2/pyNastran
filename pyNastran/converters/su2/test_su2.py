import os
import warnings
import unittest

import numpy as np
from cpylog import get_logger
warnings.simplefilter('always')
np.seterr(all='raise')


import pyNastran
from pyNastran.converters.su2.su2_reader import read_su2

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, 'converters', 'su2')

class TestSU2(unittest.TestCase):

    def test_su2_01(self):
        """tests mesh_naca0012_inv.su2"""
        log = get_logger(level='debug')
        geometry_filename = os.path.join(MODEL_PATH, 'mesh_naca0012_inv.su2')
        read_su2(geometry_filename, log=log)

    def test_su2_02(self):
        """tests sliding_interface_pipe.su2"""
        log = get_logger(level='debug')
        geometry_filename = os.path.join(MODEL_PATH, 'sliding_interface_pipe.su2')
        read_su2(geometry_filename, log=log)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()

