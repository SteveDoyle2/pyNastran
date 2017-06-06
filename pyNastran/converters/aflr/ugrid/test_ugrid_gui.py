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
from pyNastran.converters.aflr.ugrid.ugrid_io import UGRID_IO
from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.converters.nastran.nastran_to_ugrid import nastran_to_ugrid
from pyNastran.bdf.mesh_utils.extract_free_faces import write_skin_solid_faces
from pyNastran.utils.log import get_logger

pkg_path = pyNastran.__path__[0]
model_path = os.path.join(pkg_path, 'converters', 'tecplot', 'models')
nastran_path = os.path.join(pkg_path, '..', 'models')


class UGRID_GUI(UGRID_IO, FakeGUIMethods):
    """defines the UGRID 2D/3D interface"""
    def __init__(self):
        FakeGUIMethods.__init__(self)
        UGRID_IO.__init__(self)

class TestUgridGui(unittest.TestCase):
    """defines UGRID tests"""
    def test_ugrid_gui_01(self):
        """tests solid_bending.bdf"""
        nastran_filename1 = os.path.join(nastran_path, 'solid_bending', 'solid_bending.bdf')

        skin_filename = os.path.join(nastran_path, 'solid_bending', 'solid_bending_skin.bdf')
        log = get_logger(level='warning')
        write_skin_solid_faces(nastran_filename1, skin_filename,
                               write_solids=True, write_shells=True,
                               size=8, is_double=False, encoding=None,
                               punch=False, log=log)

        ugrid_filename = os.path.join(nastran_path, 'solid_bending', 'solid_bending.b8.ugrid')
        nastran_to_ugrid(skin_filename, ugrid_filename_out=ugrid_filename,
                         properties=None, check_shells=True, check_solids=True, log=log)
        assert os.path.exists(ugrid_filename), ugrid_filename
        test = UGRID_GUI()
        test.log = log
        dirname = None
        test.load_ugrid_geometry(ugrid_filename, dirname, name='main',
                                 plot=True)

    def test_ugrid_gui_02(self):
        """tests plate_with_circular_hole"""
        bdf_filename = os.path.join(nastran_path, 'plate_with_circular_hole', 'a101x.dat')
        ugrid_filename = os.path.join(nastran_path, 'plate_with_circular_hole', 'a101x.b8.ugrid')
        log = get_logger(level='warning')

        with self.assertRaises(RuntimeError):
            nastran_to_ugrid(bdf_filename, ugrid_filename_out=ugrid_filename,
                             properties=None, check_shells=True, check_solids=False, log=log)
        #assert os.path.exists(ugrid_filename), ugrid_filename
        #test = UGRID_GUI()
        #dirname = None
        #test.load_ugrid_geometry(ugrid_filename, dirname, name='main',
                                 #plot=True)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
