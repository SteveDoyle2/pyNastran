import os
import warnings
from pathlib import Path
import unittest

import numpy as np
from cpylog import get_logger

import pyNastran
from pyNastran.converters.aflr.aflr2.aflr2 import (
    read_bedge, export_to_bedge, get_triangle_input)
from pyNastran.converters.aflr.ugrid.ugrid2d_reader import UGRID2D_Reader

PKG_PATH = Path(pyNastran.__path__[0])
TEST_PATH = PKG_PATH / 'converters' / 'aflr' / 'aflr2'
warnings.simplefilter('always')
np.seterr(all='raise')


class TestBEdge(unittest.TestCase):
    """tests the bedge file format"""

    def test_bedge_1(self):
        """tests the m3 model"""
        log = get_logger(level='warning')
        #log = get_logger(level='debug')
        bedge_filename = TEST_PATH / 'm3.bedge'
        bdf_filename = TEST_PATH / 'm3.bdf'
        esp_filename = TEST_PATH / 'm3.csm'
        tri_filename = TEST_PATH / 'm3.tri'
        fixed_points_filename = TEST_PATH / 'm3.xyz'

        holes = [
            [-0.06, -.1],  # fwd
            [0.4, 0],  # mid
            [0.9, 0],  # aft
        ]
        circles = [
            # R, x, y, Npoints
            #(9., 0., 0., 50),
            (1.5, 0.4, 0., 50),
            (2.0, 0.4, 0., 50),
        ]
        regions = [
            # x, y, pid, ???
            (0.4, 1.40, 1, 0),  # just inside inner circle
            (0.4, 1.75, 2, 0),  # just inside outer circle
        ]

        model = read_bedge(
            bedge_filename, beta_reverse=179.7,
            log=log, debug=False)

        triangle_input, options = get_triangle_input(
            model, curves_to_skip=[0],
            holes=holes, regions=regions, circles=circles,
            min_angle='20',
            max_area=0.03, tri_order=1)

        is_triangle = False
        if is_triangle:  # pragma: no cover
            temp_tags_map = {
                0: 100.0,
                1: 100.0,
                2: 100.0,
            }  # inner wall
            # 3: common wall
            bc_tags = [4]  # outer wall
            model.write_tri(
                tri_filename, curves_to_skip=[0],
                holes=holes, regions=regions, circles=circles,
                min_angle='20',
                max_area=0.03,
                bc_tags=bc_tags, temp_tags_map=temp_tags_map,
                plot_clear_regions=True,
                show=True)

        model.write_esp(esp_filename)
        model2 = read_bedge(bedge_filename, beta_reverse=179.7, log=log, debug=False)
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
        log = get_logger(level='warning')
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
            ugrid_file.write(msg)

        model = UGRID2D_Reader(log=log, debug=True)
        model.read_ugrid(ugrid_filename)
        os.remove(ugrid_filename)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
