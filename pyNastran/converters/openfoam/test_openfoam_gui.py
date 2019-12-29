import os
import unittest
from cpylog import get_logger

import pyNastran
from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.converters.openfoam.block_mesh import read_block_mesh, mirror_block_mesh
from pyNastran.converters.openfoam.face_file import FaceFile
from pyNastran.converters.openfoam.openfoam_io import OpenFoamIO
from pyNastran.utils import check_path

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, 'converters', 'openfoam', 'models')


class OpenFoamGUI(OpenFoamIO, FakeGUIMethods):
    def __init__(self):
        FakeGUIMethods.__init__(self)
        self.model = OpenFoamIO(self)
        self.build_fmts(['openfoam_hex', 'openfoam_shell', 'openfoam_faces'], stop_on_failure=True)

class TestOpenFoamGUI(unittest.TestCase):

    def test_openfoam_geometry_01(self):
        """tests the ascii three plugs model"""
        log = get_logger(level='warning', encoding='utf-8')
        geometry_filename = os.path.join(MODEL_PATH, 'SnakeRiverCanyon', 'system', 'blockMeshDict')
        bdf_filename = os.path.join(MODEL_PATH, 'SnakeRiverCanyon', 'system', 'blockMeshDict.bdf')
        face_filename = os.path.join(MODEL_PATH, 'SnakeRiverCanyon', 'system', 'faces')
        check_path(geometry_filename, 'geometry_filename')
        test = OpenFoamGUI()
        test.log = log
        test.on_load_geometry(geometry_filename, geometry_format='openfoam_shell', raise_error=True)
        test.on_load_geometry(geometry_filename, geometry_format='openfoam_hex', raise_error=True)
        os.remove('points.bdf')
        #test.load_openfoam_geometry_faces(geometry_filename)

        model = read_block_mesh(geometry_filename, log=log)
        block_mesh_name_out = 'blockMeshDict.out'
        model.write_block_mesh(
            block_mesh_name_out=block_mesh_name_out, make_symmetry=False)
        model.write_block_mesh(
            block_mesh_name_out=block_mesh_name_out, make_symmetry=True)
        model.write_bdf(bdf_filename, model.nodes, model.hexas)

        mirror_block_mesh(geometry_filename, block_mesh_name_out)
        os.remove(block_mesh_name_out)
        #nodes, hexas, quads, inames, bcs

    def test_openfoam_2(self):
        point_filename = 'points'
        with open(point_filename, 'w') as point_file:
            point_file.write('0. 0. 0.\n')

        face_filename = 'faces'
        with open(face_filename, 'w') as face_file:
            face_file.write('2\n')
            face_file.write('\n')
            face_file.write('3 1 2 3\n')
            face_file.write('3 1 3 4\n')

        log = get_logger(level='warning', encoding='utf-8')
        #test = OpenFoamGUI()
        #test.log = log
        #test.load_openfoam_faces_geometry(face_filename)
        faces = FaceFile(log=log, debug=False)
        faces.read_face_file(face_filename)

        faces.read_face_file(face_filename, ifaces_to_read=[1])
        faces.read_face_file(face_filename, ifaces_to_read=[0, 1])
        os.remove(point_filename)
        os.remove(face_filename)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()

