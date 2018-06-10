"""
Defines:
 - SURF tests
"""
from __future__ import print_function
import os
import unittest

import pyNastran
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


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
