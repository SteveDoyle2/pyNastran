import os
import unittest
import numpy as np
from cpylog import get_logger

from pyNastran.converters.openfoam.block_mesh import BlockMesh, read_block_mesh, mirror_block_mesh
from pyNastran.converters.openfoam.points_file import read_points_file
from pyNastran.converters.openfoam.face_file import FaceFile
from pyNastran.converters.openfoam.boundary_file import read_boundary, read_boundary_file

from pyNastran.utils import check_path
DIRNAME = os.path.dirname(__file__)

class TestOpenFOAM(unittest.TestCase):
    def test_boundary_1(self):
        """tests the PointsFile, FaceFile, Boundary class using the Boundary"""
        #points_filename = os.path.join(DIRNAME, 'points.foam')
        boundary_filename = os.path.join(DIRNAME, 'boundary.foam')
        #with open(points_filename, 'w') as points_file:
            #points_file.write(
                #'4\n\n'
                #'(0. 0. 0.)\n'
                #'(1. 0. 0.)\n'
                #'(2. 0. 0.)\n'
                #'(3. 0. 0.)\n'
            #)
        #nfaces = 2
        boundary_msg = (
            '6\n'
            '(\n'
            '    inlet\n'
            '    {\n'
            '        type            patch;\n'
            '        nFaces          50;\n'
            '        startFace       10325;\n'
            '    }\n'
            '    outlet\n'
            '    {\n'
            '        type            patch;\n'
            '        nFaces          40;\n'
            '        startFace       10375;\n'
            '    }\n'
            '    bottom\n'
            '    {\n'
            '        type            symmetryPlane;\n'
            '        inGroups        1(symmetryPlane);\n'
            '        nFaces          25;\n'
            '        startFace       10415;\n'
            '    }\n'
            '    top\n'
            '    {\n'
            '        type            symmetryPlane;\n'
            '        inGroups        1(symmetryPlane);\n'
            '        nFaces          125;\n'
            '        startFace       10440;\n'
            '    }\n'
            '    obstacle\n'
            '    {\n'
            '        type            patch;\n'
            '        nFaces          110;\n'
            '        startFace       10565;\n'
            '    }\n'
            '    defaultFaces\n'
            '    {\n'
            '        type            empty;\n'
            '        inGroups        1(empty);\n'
            '        nFaces          10500;\n'
            '        startFace       10675;\n'
            '    }\n'
            ')\n'
            '\n'
            '// *************************************** //\n'
        )
        with open(boundary_filename, 'w') as boundary_file_obj:
            boundary_file_obj.write(boundary_msg)

        log = get_logger(level='warning', encoding='utf-8')
        boundary_file = read_boundary_file(
            boundary_filename, log=log, debug=False)
        #boundary = read_boundary(
            #point_filename, face_filename, boundary_filename,
            #log=None, debug=False)
        os.remove(boundary_filename)


    def test_points_1(self):
        """tests the PointsFile class"""
        points_filename = 'points.foam'
        with open(points_filename, 'w') as points_File:
            points_File.write(
                '4\n\n'
                '(0. 0. 0.)\n'
                '(1. 0. 0.)\n'
                '(2. 0. 0.)\n'
                '(3. 0. 0.)\n'
            )
        points = read_points_file(points_filename, ipoints_to_read=None, log=None,
                                  debug=None)
        #print(points)
        assert points.shape == (4, 3), points.shape

        points = read_points_file(points_filename, ipoints_to_read=[0, 3], log=None,
                                  debug=None)
        assert points.shape == (2, 3), points.shape
        #print(points)

        points = read_points_file(points_filename, ipoints_to_read=[3], log=None,
                                  debug=None)
        #print(points)
        assert points.shape == (1, 3), points.shape
        os.remove(points_filename)

    def test_blockmesh_1(self):
        """tests the BlockMesh class"""
        block = BlockMesh(log=None, debug=True)

        block.grading = [
            [4, 2, 6, None],
        ]
        block.nodes = np.array([
            [0., 0., 0.],
            [0., 1., 0.],
            [1., 1., 0.],
            [1., 0., 0.],

            [0., 0., 1.],
            [0., 1., 1.],
            [1., 1., 1.],
            [1., 0., 1.],
        ])
        block.hexas = [
            [1, 2, 3, 4, 5, 6, 7, 8],
        ]
        hex_id = 1
        block.make_hex_bars(hex_id)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
