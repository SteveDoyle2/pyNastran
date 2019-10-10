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
        log = get_logger(level='info')
        geometry_filename = os.path.join(MODEL_PATH, 'mesh_naca0012_inv.su2')
        geometry_filename2 = os.path.join(MODEL_PATH, 'mesh_naca0012_inv_out.su2')
        model, unused_zones = read_su2(geometry_filename, log=log)
        model.write_su2(geometry_filename2)
        os.remove(geometry_filename2)

    def test_su2_02(self):
        """tests sliding_interface_pipe.su2"""
        log = get_logger(level='info')
        geometry_filename = os.path.join(MODEL_PATH, 'sliding_interface_pipe.su2')
        geometry_filename2 = os.path.join(MODEL_PATH, 'sliding_interface_pipe_out.su2')
        model, unused_zones = read_su2(geometry_filename, log=log)
        model.write_su2(geometry_filename2)
        os.remove(geometry_filename2)

    def test_su2_03(self):
        """tests fea_fsi_wall_channel_2d_mesh_fea.su2"""
        log = get_logger(level='info')
        geometry_filename = os.path.join(MODEL_PATH, 'fea_fsi_wall_channel_2d_mesh_fea.su2')
        geometry_filename2 = os.path.join(MODEL_PATH, 'fea_fsi_wall_channel_2d_mesh_fea_out.su2')
        model, unused_zones = read_su2(geometry_filename, log=log)
        model.write_su2(geometry_filename2)
        os.remove(geometry_filename2)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()

