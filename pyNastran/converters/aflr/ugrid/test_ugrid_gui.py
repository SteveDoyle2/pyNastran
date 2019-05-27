"""
Defines:
 - UGRID 2D tests
 - UGRID 3D tests
 - SURF tests

"""
import os
import unittest
from cpylog import get_logger

import pyNastran
from pyNastran.converters.aflr.ugrid.ugrid_io import UGRID_IO
from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.converters.nastran.nastran_to_ugrid import nastran_to_ugrid
from pyNastran.bdf.mesh_utils.free_faces import write_skin_solid_faces

PKG_PATH = pyNastran.__path__[0]
UGRID_PATH = os.path.join(PKG_PATH, 'converters', 'aflr', 'ugrid', 'models')
#TECPLOT_PATH = os.path.join(PKG_PATH, 'converters', 'tecplot', 'models')
NASTRAN_PATH = os.path.join(PKG_PATH, '..', 'models')


class UGRID_GUI(UGRID_IO, FakeGUIMethods):
    """defines the UGRID 2D/3D interface"""
    def __init__(self):
        FakeGUIMethods.__init__(self)
        UGRID_IO.__init__(self, self)
        self.build_fmts(['ugrid', 'ugrid3d'], stop_on_failure=True)


class TestUgridGui(unittest.TestCase):
    """defines UGRID tests"""
    def test_ugrid_gui_01(self):
        """tests solid_bending.bdf"""
        nastran_filename1 = os.path.join(NASTRAN_PATH, 'solid_bending', 'solid_bending.bdf')

        skin_filename = os.path.join(NASTRAN_PATH, 'solid_bending', 'solid_bending_skin.bdf')
        log = get_logger(level='warning')
        write_skin_solid_faces(nastran_filename1, skin_filename,
                               write_solids=True, write_shells=True,
                               size=8, is_double=False, encoding=None,
                               punch=False, log=log)

        ugrid_filename = os.path.join(NASTRAN_PATH, 'solid_bending', 'solid_bending.b8.ugrid')
        nastran_to_ugrid(skin_filename, ugrid_filename_out=ugrid_filename,
                         properties=None, check_shells=True, check_solids=True, log=log)
        assert os.path.exists(ugrid_filename), ugrid_filename
        test = UGRID_GUI()
        test.log = log
        test.on_load_geometry(ugrid_filename, geometry_format='ugrid', raise_error=True)
        test.on_load_geometry(ugrid_filename, geometry_format='ugrid3d', raise_error=True)
        #test.load_ugrid_geometry(ugrid_filename, name='main', plot=True)
        #test.load_ugrid3d_geometry(ugrid_filename, name='main', plot=True)

    def test_ugrid_gui_02(self):
        """tests plate_with_circular_hole"""
        bdf_filename = os.path.join(NASTRAN_PATH, 'plate_with_circular_hole', 'a101x.dat')
        ugrid_filename = os.path.join(NASTRAN_PATH, 'plate_with_circular_hole', 'a101x.b8.ugrid')
        log = get_logger(level='warning')

        with self.assertRaises(RuntimeError):
            nastran_to_ugrid(bdf_filename, ugrid_filename_out=ugrid_filename,
                             properties=None, check_shells=True, check_solids=False, log=log)
        #assert os.path.exists(ugrid_filename), ugrid_filename
        #test = UGRID_GUI()
        #test.load_ugrid_geometry(ugrid_filename, name='main', plot=True)

    def test_ugrid2d_gui(self):
        """simple UGRID2D model"""
        ugrid_filename = 'quad_tri.ugrid'
        msg = (
            #(nnodes, ntrias, nquads), ntets, npyram5, npenta6, nhexas8s
            '5 1 1   0 0 0 0\n'
            '0. 0. 0.\n'
            '1. 0. 0.\n'
            '1. 1. 0.\n'
            '0. 1. 0.\n'
            '0. 2. 0.\n'
            '3 4 5\n'
            '1 2 3 4\n'
        )
        with open(ugrid_filename, 'w') as ugrid_file:
            ugrid_file.write(msg)

        log = get_logger(level='warning')
        test = UGRID_GUI()
        test.log = log
        test.on_load_geometry(ugrid_filename, geometry_format='ugrid', raise_error=True)
        #test.load_ugrid_geometry(ugrid_filename, name='main', plot=True)
        os.remove(ugrid_filename)

    def test_ugrid3d_gui_box(self):
        """simple UGRID3D box model"""
        ugrid_filename = os.path.join(UGRID_PATH, 'box.b8.ugrid')

        log = get_logger(level='warning')
        test = UGRID_GUI()
        test.log = log
        test.on_load_geometry(ugrid_filename, geometry_format='ugrid', raise_error=True)
        test.on_load_geometry(ugrid_filename, geometry_format='ugrid3d', raise_error=True)
        #test.load_ugrid_geometry(ugrid_filename, name='main', plot=True)
        #test.load_ugrid3d_geometry(ugrid_filename, name='main', plot=True)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
