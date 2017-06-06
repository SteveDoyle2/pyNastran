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

pkg_path = pyNastran.__path__[0]
model_path = os.path.join(pkg_path, 'converters', 'tecplot', 'models')
nastran_path = os.path.join(pkg_path, '..', 'models')


class SurfGui(SurfIO, FakeGUIMethods):
    """defines the UGRID 2D/3D interface"""
    def __init__(self):
        FakeGUIMethods.__init__(self)
        SurfIO.__init__(self)

class TestSurfGui(unittest.TestCase):
    """defines *.surf tests"""
    def test_surf_gui_01(self):
        """tests two_blade_wake_sym_extended.surf"""
        ugrid_filename = os.path.join(pkg_path, 'converters', 'aflr', 'ugrid', 'models',
                                      'two_blade_wake_sym_extended.surf')
        log = get_logger(level='warning')
        test = SurfGui()
        test.log = log
        test.load_surf_geometry(ugrid_filename, dirname=None)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
