import os
import unittest
import numpy as np
from cpylog import get_logger

from pyNastran.converters.openfoam.block_mesh import BlockMesh, read_block_mesh, mirror_block_mesh
from pyNastran.converters.openfoam.points_file import read_points_file
from pyNastran.converters.openfoam.face_file import FaceFile
from pyNastran.converters.openfoam.openfoam_io import OpenFoamIO
from pyNastran.utils import check_path


class TestOpenFOAM(unittest.TestCase):
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