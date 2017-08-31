from __future__ import print_function
import os
import unittest

import warnings
import numpy as np
warnings.simplefilter('always')
np.seterr(all='raise')

import pyNastran
from pyNastran.converters.aflr.aflr2.aflr2 import read_bedge, export_to_bedge
from pyNastran.converters.aflr.ugrid.ugrid2d_reader import UGRID2D_Reader

pkg_path = pyNastran.__path__[0]
TEST_PATH = os.path.join(pkg_path, 'converters', 'aflr', 'aflr2')

class TestBEdge(unittest.TestCase):
    """tests the bedge file format"""
    def test_bedge_1(self):
        """tests the m3 model"""
        bedge_filename = os.path.join(TEST_PATH, 'm3.bedge')
        bdf_filename = os.path.join(TEST_PATH, 'm3.bdf')
        fixed_points_filename = os.path.join(TEST_PATH, 'm3.xyz')

        model = read_bedge(bedge_filename, beta_reverse=179.7, log=None, debug=False)
        model2 = read_bedge(bedge_filename, beta_reverse=179.7, log=None, debug=False)
        assert len(model.nodes) == 858, 'nodes=%s' % len(model.nodes)
        assert len(model.bars) == 858, 'nbars=%s' % len(model.bars)

        model.write_nastran(bdf_filename)
        model.write_fixed_points(fixed_points_filename)

        bedge_filename_out = os.path.join(TEST_PATH, 'm3.bedge_out')
        model.merge_bedge(model2, bedge_filename_out)

        export_to_bedge(bedge_filename_out,
                        model.nodes, model.grid_bc, model.curves, model.subcurves,
                        axis=1, log=model.log)
        os.remove(bdf_filename)
        os.remove(fixed_points_filename)
        os.remove(bedge_filename_out)

    def test_ugrid2d(self):
        """simple UGRID2D model"""
        ugrid_filename = os.path.join(TEST_PATH, 'quad_tri.ugrid')
        #if not os.path.exists(ugrid_filename):
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
            #print('msg = %r' % msg)
            #print('----------------')
            ugrid_file.write(msg)

        model = UGRID2D_Reader(log=None, debug=True)
        model.read_ugrid(ugrid_filename)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
