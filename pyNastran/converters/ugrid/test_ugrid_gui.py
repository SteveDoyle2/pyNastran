"""
Defines:
 - UGRID 2D tests
 - UGRID 3D tests
 - SURF tests
"""
from __future__ import print_function
import os
import unittest

import pyNastran
from pyNastran.converters.ugrid.ugrid_io import UGRID_IO
from pyNastran.converters.ugrid.surf_io import SurfIO
from pyNastran.gui.testing_methods import GUIMethods
from pyNastran.converters.nastran.nastran_to_ugrid import nastran_to_ugrid

pkg_path = pyNastran.__path__[0]
model_path = os.path.join(pkg_path, 'converters', 'tecplot', 'models')
nastran_path = os.path.join(pkg_path, '..', 'models')


class UGRID_GUI(UGRID_IO, SurfIO, GUIMethods):
    """defines the UGRID 2D/3D interface"""
    def __init__(self):
        GUIMethods.__init__(self)
        UGRID_IO.__init__(self)
        SurfIO.__init__(self)

class TestUgridGUI(unittest.TestCase):
    """defines UGRID tests"""
    def test_ugrid_gui_01(self):
        """tests two_blade_wake_sym_extended.surf"""
        ugrid_filename = os.path.join(pkg_path, 'converters', 'ugrid', 'models',
                                      'two_blade_wake_sym_extended.surf')
        test = UGRID_GUI()
        test.load_surf_geometry(ugrid_filename, dirname=None)

    def test_ugrid_gui_02(self):
        """tests solid_bending.bdf"""
        nastran_filename1 = os.path.join(nastran_path, 'solid_bending', 'solid_bending.bdf')
        from pyNastran.bdf.mesh_utils.extract_free_faces import write_skin_solid_faces

        skin_filename = os.path.join(nastran_path, 'solid_bending', 'solid_bending_skin.bdf')
        write_skin_solid_faces(nastran_filename1, skin_filename,
                               write_solids=True, write_shells=True,
                               size=8, is_double=False, encoding=None,
                               punch=False)

        ugrid_filename = os.path.join(nastran_path, 'solid_bending', 'solid_bending.b8.ugrid')
        nastran_to_ugrid(skin_filename, ugrid_filename_out=ugrid_filename,
                         properties=None, check_shells=True, check_solids=True)
        assert os.path.exists(ugrid_filename), ugrid_filename
        test = UGRID_GUI()
        dirname = None
        test.load_ugrid_geometry(ugrid_filename, dirname, name='main',
                                 plot=True)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
