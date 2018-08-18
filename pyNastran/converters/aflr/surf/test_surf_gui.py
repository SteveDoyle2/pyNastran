"""
Defines:
 - SURF tests
"""
from __future__ import print_function
import os
import unittest

import pyNastran
from pyNastran.converters.nastran.nastran_to_surf import nastran_to_surf, read_bdf
from pyNastran.converters.aflr.surf.surf_io import SurfIO
from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.utils.log import get_logger

PKG_PATH = pyNastran.__path__[0]
model_path = os.path.join(PKG_PATH, 'converters', 'tecplot', 'models')
nastran_path = os.path.join(PKG_PATH, '..', 'models')


class SurfGui(SurfIO, FakeGUIMethods):
    """defines the UGRID 2D/3D interface"""
    def __init__(self):
        FakeGUIMethods.__init__(self)
        SurfIO.__init__(self, self)
        self.build_fmts(['surf'], stop_on_failure=True)

class TestSurfGui(unittest.TestCase):
    """defines *.surf tests"""
    def test_surf_gui_01(self):
        """tests two_blade_wake_sym_extended.surf"""
        ugrid_filename = os.path.join(PKG_PATH, 'converters', 'aflr', 'ugrid', 'models',
                                      'two_blade_wake_sym_extended.surf')
        log = get_logger(level='warning')
        test = SurfGui()
        test.log = log
        test.on_load_geometry(ugrid_filename, geometry_format='surf', raise_error=True)
        #test.load_surf_geometry(ugrid_filename)

    def test_surf_01(self):
        """tests two_blade_wake_sym_extended.surf"""
        MODEL_PATH = os.path.join(PKG_PATH, '..', 'models')
        bdf_filename = os.path.join(MODEL_PATH, 'iSat', 'ISat_Launch_Sm_Rgd.dat')
        surf_filename = os.path.join(MODEL_PATH, 'iSat', 'ISat_Launch_Sm_Rgd.dat')
        pid_to_element_flags = {
        }
        bdf_model = read_bdf(bdf_filename)
        for pid, prop in bdf_model.properties.items():
            if prop.type in ['PSHELL', 'PCOMP']:
                # initial_normal_spacing, bl_thickness, grid_bc
                pid_to_element_flags[pid] = [0.01, 0.1, 1]

        #ugrid_filename = os.path.join(PKG_PATH, 'converters', 'aflr', 'ugrid', 'models',
                                      #'two_blade_wake_sym_extended.surf')
        #log = get_logger(level='warning')
        with self.assertRaises(RuntimeError):
            nastran_to_surf(bdf_model, pid_to_element_flags, surf_filename,
                            renumber_pids=None,
                            line_map=None, scale=1.0,
                            tol=1e-10, xref=True)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
