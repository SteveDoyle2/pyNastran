import io
import os
import unittest
import warnings
from pathlib import Path

import numpy as np
#from cpylog import get_logger

import pyNastran
from pyNastran.converters.fluent.fluent import Fluent, read_vrt, read_cell, read_daten, read_fluent
from pyNastran.converters.fluent.nastran_to_fluent import nastran_to_fluent
from pyNastran.converters.fluent.fluent_to_tecplot import fluent_to_tecplot
from pyNastran.converters.tecplot.tecplot import read_tecplot
from pyNastran.converters.fluent.ugrid_to_fluent import ugrid_to_fluent_filename

warnings.simplefilter('always')
np.seterr(all='raise')

PKG_PATH = pyNastran.__path__[0]
BWB_PATH = Path(os.path.join(PKG_PATH, '..', 'models', 'bwb'))
TEST_PATH = Path(os.path.join(PKG_PATH, 'converters', 'fluent'))
UGRID_PATH = os.path.join(PKG_PATH, 'converters', 'aflr', 'ugrid', 'models')


class TestFluent(unittest.TestCase):
    def test_fluent_vrt_01(self):
        vrt_filename = TEST_PATH / 'test.cel'
        with open(vrt_filename, 'w') as vrt_file:
            vrt_file.write("""PROSTAR_VERTEX
     4000         0         0         0         0         0         0         0
             1     10.86546244027553      0.15709716073078      1.16760397346681
             2      8.91578650948743      1.88831278193532      1.27715047467296
             3      8.91578650948743      1.88831278193532      1.27715047467296
             4      8.91578650948743      1.88831278193532      1.27715047467296""")
        node_id, xyz = read_vrt(vrt_filename)
        assert len(xyz) == 4, xyz
        assert np.array_equal(node_id, [1, 2, 3, 4]), node_id

    def test_fluent_cell_01(self):
        cell_filename = TEST_PATH / 'test.cel'
        with open(cell_filename, 'w') as cell_file:
            cell_file.write("""PROSTAR_CELL
             4000         0         0         0         0         0         0         0
                 1          3          3          3          4
                 1          1          2          3""")
        (quads, tris), (element_ids, regions) = read_cell(cell_filename)
        assert len(quads) == 0, quads
        assert len(tris) == 1, tris
        assert np.array_equal(element_ids, [1]), element_ids
        assert np.array_equal(regions, [3]), regions

    def test_fluent_cell_02(self):
        cell_filename = TEST_PATH / 'test.cel'
        with open(cell_filename, 'w') as cell_file:
            cell_file.write("""PROSTAR_CELL
             4000         0         0         0         0         0         0         0
                 1          3          4         100         4
                 1          1          2          3          4""")
        (quads, tris), (element_ids, regions) = read_cell(cell_filename)
        assert len(quads) == 1, quads
        assert len(tris) == 0, tris
        assert np.array_equal(element_ids, [1]), element_ids
        assert np.array_equal(regions, [100]), regions

    def test_fluent_daten_01(self):
        daten_filename = TEST_PATH / 'test.daten'
        with open(daten_filename, 'w') as daten_file:
            daten_file.write("""# Shell Id, Cf Components X, Cf Components Y, Cf Components Z, Pressure Coefficient
         1       0.00158262927085      -0.00013470219276       0.00033890815696      -0.09045402854837
   4175456       0.00005536470364      -0.00050913009274       0.00226408640311      -0.44872838534457""")
        element_id, titles, results = read_daten(daten_filename, scale=2.0)
        assert len(element_id) == 2, element_id
        assert np.array_equal(element_id, [1, 4175456]), element_id

    def test_nastran_to_fluent(self):
        nastran_filename = BWB_PATH / 'bwb_saero.bdf'
        vrt_filename2 = BWB_PATH / 'bwb_saero2.vrt'
        vrt_filename = BWB_PATH / 'bwb_saero.vrt'
        cel_filename = BWB_PATH / 'bwb_saero.cel'
        daten_filename = BWB_PATH / 'bwb_saero.daten'
        h5_filename = BWB_PATH / 'bwb_saero.h5'
        tecplot_filename = BWB_PATH / 'bwb_saero.plt'
        if h5_filename.exists():
            os.remove(h5_filename)

        #model = Fluent()
        #is_loaded = model.read_h5(h5_filename)
        #assert is_loaded is True, h5_filename

        nastran_to_fluent(nastran_filename, vrt_filename)

        node_id, xyz = read_vrt(vrt_filename)
        assert len(xyz) == 10135, xyz.shape
        (quads, tris), (element_ids, regions) = read_cell(cel_filename)
        element_id, titles, results = read_daten(daten_filename, scale=2.0)
        model = read_fluent(vrt_filename)
        model.write_fluent(vrt_filename2)
        model2 = read_fluent(vrt_filename2)
        model2.get_filtered_data(
            regions_to_include=[1101, 1501, 1601])
        #print(np.unique(model2.region))

        tecplot = fluent_to_tecplot(vrt_filename, tecplot_filename)
        #read_tecplot(tecplot_filename)  # TODO: fix tecplot parsing; the file is correct...

    def test_ugrid3d_gui_box(self):
        """simple UGRID3D box model"""
        ugrid_filename = os.path.join(UGRID_PATH, 'box.b8.ugrid')
        fluent_filename = os.path.join(UGRID_PATH, 'box.vrt')
        fluent_model = ugrid_to_fluent_filename(ugrid_filename, fluent_filename)

    def test_fluent_area(self):
        model = Fluent()
        model.node_id = np.array([1, 2, 3, 4])
        model.xyz = np.array([
            [0., 0., 0.],
            [1., 0., 0.],
            [1., 1., 0.],
            [0., 1., 0.],
        ])
        model.tris = np.array([
            [1, 10, 1, 2, 3],
        ])
        model.quads = np.array([
            [2, 12, 1, 2, 3, 4],
        ])
        (tri_area, quad_area,
         tri_centroid, quad_centroid,
         tri_normal, quad_normal) = model.get_area_centroid_normal(
            model.tris, model.quads)
        assert np.allclose(tri_area[0], 0.5)
        assert np.allclose(quad_area[0], 1.0)
        assert np.allclose(tri_centroid, [2/3, 1/3, 0.]), tri_centroid
        assert np.allclose(quad_centroid, [0.5, 0.5, 0.]), quad_centroid
        assert np.allclose(tri_normal, [0., 0., 1.]), tri_normal
        assert np.allclose(quad_normal, [0., 0., 1.]), quad_normal

        (tri_area, quad_area,
         tri_centroid, quad_centroid) = model.get_area_centroid(
            model.tris, model.quads)
        assert np.allclose(tri_area[0], 0.5)
        assert np.allclose(quad_area[0], 1.0)
        assert np.allclose(tri_centroid, [2/3, 1/3, 0.]), tri_centroid
        assert np.allclose(quad_centroid, [0.5, 0.5, 0.]), quad_centroid


def main():  # pragma: no cover
    import time
    time0 = time.time()
    unittest.main()
    print("dt = %s" % (time.time() - time0))

if __name__ == '__main__':  # pragma: no cover
    main()
