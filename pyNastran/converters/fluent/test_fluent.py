import os
import unittest
import warnings
from pathlib import Path

import numpy as np
from cpylog import SimpleLogger

import pyNastran
from pyNastran.converters.fluent.fluent import Fluent, read_vrt, read_cell, read_daten, read_fluent
from pyNastran.converters.fluent.nastran_to_fluent import nastran_to_fluent
from pyNastran.converters.fluent.fluent_to_tecplot import fluent_to_tecplot
from pyNastran.converters.fluent.ugrid_to_fluent import ugrid_to_fluent_filename

try:
    import tables  # noqa: F401

    IS_TABLES = True
except ImportError:  # pragma: no cover
    IS_TABLES = False

warnings.simplefilter("always")
np.seterr(all="raise")

PKG_PATH = pyNastran.__path__[0]
BWB_PATH = Path(os.path.join(PKG_PATH, "..", "models", "bwb"))
TEST_PATH = Path(os.path.join(PKG_PATH, "converters", "fluent"))
UGRID_PATH = os.path.join(PKG_PATH, "converters", "aflr", "ugrid", "models")


class TestFluent(unittest.TestCase):
    def test_fluent_vrt_01(self):
        vrt_filename = TEST_PATH / "test.cel"
        with open(vrt_filename, "w") as vrt_file:
            vrt_file.write("""PROSTAR_VERTEX
     4000         0         0         0         0         0         0         0
             1     10.86546244027553      0.15709716073078      1.16760397346681
             2      8.91578650948743      1.88831278193532      1.27715047467296
             3      8.91578650948743      1.88831278193532      1.27715047467296
             4      8.91578650948743      1.88831278193532      1.27715047467296""")
        log = SimpleLogger(level="warning")
        node_id, xyz = read_vrt(vrt_filename, log)
        assert len(xyz) == 4, xyz
        assert np.array_equal(node_id, [1, 2, 3, 4]), node_id

    def test_fluent_cell_01(self):
        cell_filename = TEST_PATH / "test.cel"
        with open(cell_filename, "w") as cell_file:
            cell_file.write("""PROSTAR_CELL
             4000         0         0         0         0         0         0         0
                 1          3          3          3          4
                 1          1          2          3""")
        log = SimpleLogger(level="warning")
        (quads, tris), element_ids = read_cell(cell_filename, log)
        regions = tris[:, 1]
        assert len(quads) == 0, quads
        assert len(tris) == 1, tris
        assert np.array_equal(element_ids, [1]), element_ids
        assert np.array_equal(regions, [3]), regions

    def test_fluent_cell_02(self):
        cell_filename = TEST_PATH / "test.cel"
        with open(cell_filename, "w") as cell_file:
            cell_file.write("""PROSTAR_CELL
             4000         0         0         0         0         0         0         0
                 1          3          4         100         4
                 1          1          2          3          4""")
        log = SimpleLogger(level="warning")
        (quads, tris), element_ids = read_cell(cell_filename, log)
        assert len(quads) == 1, quads
        assert len(tris) == 0, tris
        regions = quads[:, 1]
        assert np.array_equal(element_ids, [1]), element_ids
        assert np.array_equal(regions, [100]), regions

    def test_fluent_daten_01(self):
        daten_filename = TEST_PATH / "test.daten"
        with open(daten_filename, "w") as daten_file:
            daten_file.write("""# Shell Id, Cf Components X, Cf Components Y, Cf Components Z, Pressure Coefficient
         1       0.00158262927085      -0.00013470219276       0.00033890815696      -0.09045402854837
   4175456       0.00005536470364      -0.00050913009274       0.00226408640311      -0.44872838534457""")
        element_id, titles, results = read_daten(daten_filename, scale=2.0)
        assert len(element_id) == 2, element_id
        assert np.array_equal(element_id, [1, 4175456]), element_id

    def test_nastran_to_fluent(self):
        nastran_filename = BWB_PATH / "bwb_saero.bdf"
        vrt_filename2 = BWB_PATH / "bwb_saero2.vrt"
        vrt_filename = BWB_PATH / "bwb_saero.vrt"
        cel_filename = BWB_PATH / "bwb_saero.cel"
        daten_filename = BWB_PATH / "bwb_saero.daten"
        h5_filename = BWB_PATH / "bwb_saero.h5"
        h5_filename2 = BWB_PATH / "bwb_saero2.h5"
        tecplot_filename = BWB_PATH / "bwb_saero.plt"
        if h5_filename.exists():
            os.remove(h5_filename)
        if h5_filename2.exists():
            os.remove(h5_filename2)

        # model = Fluent()
        # is_loaded = model.read_h5(h5_filename)
        # assert is_loaded is True, h5_filename

        log = SimpleLogger("warning")
        nastran_to_fluent(nastran_filename, vrt_filename, log=log)

        node_id, xyz = read_vrt(vrt_filename, log)
        assert len(xyz) == 10135, xyz.shape
        (quads, tris), element_ids = read_cell(cel_filename, log)
        element_id, titles, results = read_daten(daten_filename, scale=2.0)
        model = read_fluent(vrt_filename, log=log)
        model.write_fluent(vrt_filename2)
        model2 = read_fluent(vrt_filename2)
        model2.get_filtered_data(regions_to_include=[1101, 1501, 1601])
        # print(np.unique(model2.region))

        tecplot = fluent_to_tecplot(vrt_filename, tecplot_filename)
        # read_tecplot(tecplot_filename)  # TODO: fix tecplot parsing; the file is correct...

    def test_ugrid3d_gui_box(self):
        """simple UGRID3D box model"""
        ugrid_filename = os.path.join(UGRID_PATH, "box.b8.ugrid")
        fluent_filename = os.path.join(UGRID_PATH, "box.vrt")
        fluent_model = ugrid_to_fluent_filename(ugrid_filename, fluent_filename)

    def test_fluent_area(self):
        model = Fluent()
        model.node_id = np.array([1, 2, 3, 4])
        model.xyz = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 0.0],
            ]
        )
        model.tris = np.array(
            [
                [1, 10, 1, 2, 3],
            ]
        )
        model.quads = np.array(
            [
                [2, 12, 1, 2, 3, 4],
            ]
        )
        (tri_area, quad_area, tri_centroid, quad_centroid, tri_normal, quad_normal) = (
            model.get_area_centroid_normal(model.tris, model.quads)
        )
        assert np.allclose(tri_area[0], 0.5)
        assert np.allclose(quad_area[0], 1.0)
        assert np.allclose(tri_centroid, [2 / 3, 1 / 3, 0.0]), tri_centroid
        assert np.allclose(quad_centroid, [0.5, 0.5, 0.0]), quad_centroid
        assert np.allclose(tri_normal, [0.0, 0.0, 1.0]), tri_normal
        assert np.allclose(quad_normal, [0.0, 0.0, 1.0]), quad_normal

        (tri_area, quad_area, tri_centroid, quad_centroid) = model.get_area_centroid(
            model.tris, model.quads
        )
        assert np.allclose(tri_area[0], 0.5)
        assert np.allclose(quad_area[0], 1.0)
        assert np.allclose(tri_centroid, [2 / 3, 1 / 3, 0.0]), tri_centroid
        assert np.allclose(quad_centroid, [0.5, 0.5, 0.0]), quad_centroid

    def test_slice_by_element_id(self):
        """keep tri 5 and quad 3; everything else is removed in place"""
        model = _tri_quad_model()
        out = model.slice_by_element_id(tri_ids=[5], quad_ids=np.array([3]))
        assert out is None

        assert np.array_equal(model.tris, [[5, 9, 2, 5, 6]]), model.tris
        assert np.array_equal(model.quads, [[3, 7, 2, 5, 6, 3]]), model.quads
        assert np.array_equal(model.result_element_id, [3, 5]), model.result_element_id
        assert np.array_equal(model.results, [[30.0, -3.0], [50.0, -5.0]]), model.results
        assert np.array_equal(model.element_ids, [3, 5]), model.element_ids
        # nodes 1 and 4 were only used by the removed elements
        assert np.array_equal(model.node_id, [2, 3, 5, 6]), model.node_id
        assert np.array_equal(
            model.xyz, [[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 0.0, 0.0], [2.0, 1.0, 0.0]]
        ), model.xyz
        assert np.array_equal(model.region, [7, 9]), model.region

        # the sliced model still works
        (
            tri_area,
            quad_area,
            unused_tri_centroid,
            unused_quad_centroid,
            unused_tri_normal,
            unused_quad_normal,
        ) = model.get_area_centroid_normal(model.tris, model.quads)
        assert np.allclose(tri_area, [0.5]), tri_area
        assert np.allclose(quad_area, [1.0]), quad_area

    def test_slice_by_element_id_tris_or_quads_only(self):
        model = _tri_quad_model()
        model.slice_by_element_id(tri_ids=[2, 4])
        assert model.quads.shape == (0, 6), model.quads.shape
        assert np.array_equal(model.tris[:, 0], [2, 4]), model.tris
        assert np.array_equal(model.result_element_id, [2, 4])
        assert np.array_equal(model.node_id, [1, 2, 3, 4]), model.node_id

        model = _tri_quad_model()
        model.slice_by_element_id(quad_ids=[1])
        assert model.tris.shape == (0, 5), model.tris.shape
        assert np.array_equal(model.quads[:, 0], [1]), model.quads
        assert np.array_equal(model.results, [[10.0, -1.0]]), model.results
        assert np.array_equal(model.node_id, [1, 2, 3, 4]), model.node_id

    def test_slice_by_element_id_unsorted_results(self):
        """result_element_id doesn't need to be sorted"""
        model = _tri_quad_model()
        iorder = np.array([3, 0, 4, 2, 1])
        model.result_element_id = model.result_element_id[iorder]
        model.results = model.results[iorder, :]
        model.slice_by_element_id(tri_ids=[2, 4], quad_ids=[1])
        # results stay paired with their ids
        assert np.array_equal(model.results[:, 0], 10.0 * model.result_element_id)
        assert np.array_equal(np.sort(model.result_element_id), [1, 2, 4])

    def test_slice_by_element_id_errors(self):
        model = _tri_quad_model()
        with self.assertRaises(KeyError):
            model.slice_by_element_id(tri_ids=[1])  # 1 is a quad
        with self.assertRaises(KeyError):
            model.slice_by_element_id(quad_ids=[2])  # 2 is a tri
        with self.assertRaises(KeyError):
            model.slice_by_element_id(tri_ids=[99])
        with self.assertRaises(ValueError):
            model.slice_by_element_id()
        # failed calls don't change the model
        assert len(model.tris) == 3 and len(model.quads) == 2
        assert len(model.result_element_id) == 5

    def test_element_bbox(self):
        """element_bbox returns correct per-element min/max from node coords"""
        model = _tri_quad_model()
        quad_eids, tri_eids, quad_bbox, tri_bbox = model.element_bbox()
        assert np.array_equal(quad_eids, [1, 3])
        assert np.array_equal(tri_eids, [2, 4, 5])
        # quad 1: nodes 1,2,3,4 -> x=[0,1], y=[0,1]
        assert np.allclose(quad_bbox[0, 0], [0.0, 0.0, 0.0])
        assert np.allclose(quad_bbox[0, 1], [1.0, 1.0, 0.0])
        # quad 3: nodes 2,5,6,3 -> x=[1,2], y=[0,1]
        assert np.allclose(quad_bbox[1, 0], [1.0, 0.0, 0.0])
        assert np.allclose(quad_bbox[1, 1], [2.0, 1.0, 0.0])
        # tri 5: nodes 2,5,6 -> x=[1,2], y=[0,1]
        assert np.allclose(tri_bbox[2, 0], [1.0, 0.0, 0.0])
        assert np.allclose(tri_bbox[2, 1], [2.0, 1.0, 0.0])

    def test_filter_elements_by_xyz(self):
        """filter_elements_by_xyz with xmin/xmax/ymin/ymax/zmin/zmax.

        _tri_quad_model nodes: x=0 (nodes 1,4), x=1 (nodes 2,3), x=2 (nodes 5,6).
        All y in [0, 1], all z=0.
        Elements:
          quad 1 (nodes 1,2,3,4) x=[0,1] y=[0,1]
          tri  2 (nodes 1,2,3)   x=[0,1] y=[0,1]
          quad 3 (nodes 2,5,6,3) x=[1,2] y=[0,1]
          tri  4 (nodes 1,3,4)   x=[0,1] y=[0,1]
          tri  5 (nodes 2,5,6)   x=[1,2] y=[0,1]
        """
        # xmax=1.5 keeps elements with all nodes x <= 1.5 -> eids 1,2,4
        model = _tri_quad_model()
        model.filter_elements_by_xyz(xmax=1.5)
        kept = np.sort(model.result_element_id)
        assert np.array_equal(kept, [1, 2, 4]), kept

        # xmin=0.5 keeps elements with all nodes x >= 0.5 -> eids 3,5
        model = _tri_quad_model()
        model.filter_elements_by_xyz(xmin=0.5)
        kept = np.sort(model.result_element_id)
        assert np.array_equal(kept, [3, 5]), kept

        # xmin=0.5, xmax=1.5 -> nothing passes (no element fully in [0.5, 1.5])
        model = _tri_quad_model()
        with self.assertRaises(ValueError):
            model.filter_elements_by_xyz(xmin=0.5, xmax=1.5)

        # None means no restriction -> everything passes
        model = _tri_quad_model()
        model.filter_elements_by_xyz()
        kept = np.sort(model.result_element_id)
        assert np.array_equal(kept, [1, 2, 3, 4, 5]), kept

        # ymax=0.5 -> no element has all nodes y <= 0.5 (nodes 3,4,6 are at y=1)
        # actually quad 1 has node 3 at y=1, so it fails; but let's check
        # every element has at least one node at y=1, so nothing passes
        model = _tri_quad_model()
        with self.assertRaises(ValueError):
            model.filter_elements_by_xyz(ymax=0.5)

        # zmin=-1, zmax=1 -> all z=0 so everything passes
        model = _tri_quad_model()
        model.filter_elements_by_xyz(zmin=-1.0, zmax=1.0)
        kept = np.sort(model.result_element_id)
        assert np.array_equal(kept, [1, 2, 3, 4, 5]), kept

    def test_rotate_regions(self):
        """rotate 2 sets (flap + tail) about basic axes; region 1 doesn't move"""
        model = _three_region_model()
        xyz0 = model.xyz.copy()
        eye = np.eye(3)
        model.rotate_regions(
            [
                ([2], [2.0, 0.0, 0.0], eye, "z", 90.0),  # flap: x -> y
                ([3], np.zeros(3), eye, "x", 180.0),  # tail: (y, z) -> (-y, -z)
            ]
        )
        assert np.array_equal(model.xyz[:4], xyz0[:4])  # region 1
        # flap node (3, 1, 0): d=(1, 1, 0) -> (-1, 1, 0) -> (1, 1, 0)
        assert np.allclose(
            model.xyz[4:8], [[2.0, 0.0, 0.0], [2.0, 1.0, 0.0], [1.0, 1.0, 0.0], [1.0, 0.0, 0.0]]
        ), model.xyz[4:8]
        assert np.allclose(model.xyz[8:], xyz0[8:] * [1.0, -1.0, -1.0]), model.xyz[8:]
        # rigid: areas don't change
        _, quad_area, _, _ = model.get_area_centroid(model.tris, model.quads)
        assert np.allclose(quad_area, 1.0), quad_area

    def test_rotate_regions_local_axis(self):
        """local axes: x=basic y, y=-basic x, z=basic z"""
        hinge_axes = np.array([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
        origin = [2.0, 0.0, 0.0]

        # +90 about basic y turns z toward x, so x -> -z
        model = _three_region_model()
        model.rotate_regions([([2], origin, np.eye(3), "y", 90.0)])
        assert np.allclose(
            model.xyz[4:8], [[2.0, 0.0, 0.0], [2.0, 0.0, -1.0], [2.0, 1.0, -1.0], [2.0, 1.0, 0.0]]
        ), model.xyz[4:8]

        # local x is basic y, so it's the same rotation
        model2 = _three_region_model()
        model2.rotate_regions([([2], origin, hinge_axes, "x", 90.0)])
        assert np.allclose(model2.xyz, model.xyz), model2.xyz

        # local y is -basic x, so -90 about it == +90 about basic x
        model3 = _three_region_model()
        model3.rotate_regions([([2], origin, hinge_axes, "y", -90.0)])
        model4 = _three_region_model()
        model4.rotate_regions([([2], origin, np.eye(3), "x", 90.0)])
        assert np.allclose(model3.xyz, model4.xyz), model3.xyz

        # local z is basic z
        model5 = _three_region_model()
        model5.rotate_regions([([2], origin, hinge_axes, "z", 30.0)])
        model6 = _three_region_model()
        model6.rotate_regions([([2], origin, np.eye(3), "z", 30.0)])
        assert np.allclose(model5.xyz, model6.xyz), model5.xyz

    def test_rotate_regions_errors(self):
        model = _three_region_model()
        xyz0 = model.xyz.copy()
        eye = np.eye(3)
        zero = np.zeros(3)
        with self.assertRaises(KeyError):
            model.rotate_regions([([99], zero, eye, "x", 10.0)])
        with self.assertRaises(ValueError):  # not orthonormal
            model.rotate_regions([([2], zero, 2 * eye, "x", 10.0)])
        with self.assertRaises(ValueError):  # left-handed
            model.rotate_regions([([2], zero, np.diag([1.0, 1.0, -1.0]), "x", 10.0)])
        with self.assertRaises(ValueError):  # region 2 is in 2 sets
            model.rotate_regions([([2], zero, eye, "x", 10.0), ([2, 3], zero, eye, "y", 5.0)])
        with self.assertRaises(ValueError):  # bad axis
            model.rotate_regions([([2], zero, eye, "w", 10.0)])
        with self.assertRaises(ValueError):  # bad angle
            model.rotate_regions([([2], zero, eye, "x", np.nan)])
        with self.assertRaises(ValueError):  # missing axis/angle
            model.rotate_regions([([2], zero, eye)])
        # failed calls don't change the model
        assert np.array_equal(model.xyz, xyz0)

    def test_rotate_regions_split_shared_nodes(self):
        """tri 5 (region 9) shares nodes 2/5/6 with quad 3 (region 7) -> split"""
        model = _tri_quad_model()
        xyz0 = model.xyz.copy()
        model.rotate_regions([([9], np.zeros(3), np.eye(3), "x", 180.0)])

        # new nodes 7/8/9 are copies of 2/5/6; tri 5 uses them
        assert np.array_equal(model.node_id, [1, 2, 3, 4, 5, 6, 7, 8, 9]), model.node_id
        assert np.array_equal(model.tris[2], [5, 9, 7, 8, 9]), model.tris
        assert np.allclose(model.xyz[6:], xyz0[[1, 4, 5]] * [1.0, -1.0, -1.0]), model.xyz
        # the non-rotated elements are untouched
        assert np.array_equal(model.xyz[:6], xyz0)
        assert np.array_equal(model.quads, _tri_quad_model().quads)
        assert np.array_equal(model.tris[:2], _tri_quad_model().tris[:2])
        tri_area, quad_area, _, _ = model.get_area_centroid(model.tris, model.quads)
        assert np.allclose(quad_area, 1.0) and np.allclose(tri_area, 0.5)

    def test_rotate_regions_split_between_sets(self):
        """regions 8 and 9 share node 2 -> each set ends up with its own copy"""
        model = _tri_quad_model()
        xyz0 = model.xyz.copy()
        eye = np.eye(3)
        model.rotate_regions(
            [([8], np.zeros(3), eye, "x", 180.0), ([9], np.zeros(3), eye, "z", 0.0)]
        )
        # every node is used by exactly 1 set or the fixed region 7
        tri_nodes = model.tris[:, 2:]
        quad_nodes = model.quads[:, 2:]
        set8 = np.unique(tri_nodes[model.tris[:, 1] == 8])
        set9 = np.unique(tri_nodes[model.tris[:, 1] == 9])
        assert len(np.intersect1d(set8, set9)) == 0
        assert len(np.intersect1d(set8, quad_nodes)) == 0
        assert len(np.intersect1d(set9, quad_nodes)) == 0
        assert np.array_equal(model.xyz[:6], xyz0)  # originals stay with region 7

    def test_rotate_regions_cf_warning(self):
        messages = []
        log = SimpleLogger(
            level="debug", log_func=lambda typ, filename, lineno, msg: messages.append((typ, msg))
        )
        model = _three_region_model()
        model.log = log
        rotation = ([2], [2.0, 0.0, 0.0], np.eye(3), "z", 5.0)  # origin on region 2
        model.titles = np.array(["ElementID", "Pressure Coefficient"])
        model.rotate_regions([rotation])
        assert not any(typ == "WARNING" for typ, msg in messages), messages

        model.titles = np.array(["ElementID", "Cf Components X", "Pressure Coefficient"])
        model.rotate_regions([rotation])
        warnings_ = [msg for typ, msg in messages if typ == "WARNING"]
        assert len(warnings_) == 1 and "Cf Components X" in warnings_[0], messages

    def test_rotate_regions_far_origin(self):
        """region 2 is 1.4 across; an origin ~1000 away warns, the hinge edge doesn't"""
        messages = []
        log = SimpleLogger(
            level="debug", log_func=lambda typ, filename, lineno, msg: messages.append((typ, msg))
        )
        model = _three_region_model()
        model.log = log
        model.rotate_regions([([2], [2.0, 0.5, 0.0], np.eye(3), "x", 5.0)])  # on the edge
        assert not any(typ == "WARNING" for typ, msg in messages), messages

        model.rotate_regions([([2], [916.0, 401.0, 0.0], np.eye(3), "z", 5.0)])
        warnings_ = [msg for typ, msg in messages if typ == "WARNING"]
        assert len(warnings_) == 1 and "far away" in warnings_[0], messages

    def test_write_fluent_interleaved_ids(self):
        """quads 1/3 and tris 2/4/5 are interleaved; the daten must be written
        in element id order or get_filtered_data reads the wrong rows
        """
        model = _tri_quad_model()
        model.rotate_regions([([9], [2.0, 0.0, 0.0], np.eye(3), "x", 10.0)])
        fluent_filename = TEST_PATH / "interleaved.cel"
        model.write_fluent(fluent_filename)
        try:
            model2 = read_fluent(
                fluent_filename, auto_read_write_h5=False, log=SimpleLogger(level="warning")
            )
            assert np.array_equal(model2.result_element_id, [1, 2, 3, 4, 5])
            element_id, tris, quads, region, results = model2.get_filtered_data()
            # results are [quads, tris]; results[:, 0] = 10 * eid
            assert np.array_equal(results[:, 0], [10.0, 30.0, 20.0, 40.0, 50.0]), results
            assert np.allclose(model2.xyz, model.xyz)
        finally:
            for ext in (".cel", ".vrt", ".daten"):
                os.remove(fluent_filename.with_suffix(ext))

    def test_get_filtered_data_twice(self):
        """get_filtered_data must not reorder self.results; results[:, 0] = 10 * eid"""
        model = _tri_quad_model()
        results0 = model.results.copy()
        out1 = model.get_filtered_data()
        out2 = model.get_filtered_data()
        assert np.array_equal(model.results, results0), model.results
        for array1, array2 in zip(out1, out2):
            assert np.array_equal(array1, array2)
        # [quads, tris] order
        assert np.array_equal(out2[4][:, 0], [10.0, 30.0, 20.0, 40.0, 50.0]), out2[4]

    def test_get_filtered_data_region_model(self):
        """the filtered model's results line up with its result_element_id"""
        model = _tri_quad_model()
        element_id, tris, quads, region, results = model.get_filtered_data(
            regions_to_include=[7, 9]
        )
        assert np.array_equal(element_id, [1, 3, 5]), element_id
        assert np.array_equal(region, [7, 7, 9]), region
        assert np.array_equal(results[:, 0], [10.0, 30.0, 50.0]), results  # [quads, tris]

        model2 = model.get_filtered_data(regions_to_include=[7, 9], return_model=True)
        assert np.array_equal(model2.results[:, 0], 10.0 * model2.result_element_id), model2.results
        # filtering the filtered model gives the same answer
        out = model2.get_filtered_data()
        assert np.array_equal(out[4][:, 0], [10.0, 30.0, 50.0]), out[4]

    @unittest.skipIf(not IS_TABLES, "pytables is not installed")
    def test_read_fluent_stale_h5(self):
        """an h5 older than the vrt/cel/daten is ignored and rewritten"""
        model = _tri_quad_model()
        fluent_filename = TEST_PATH / "stale_h5.cel"
        h5_filename = fluent_filename.with_suffix(".h5")
        messages = []
        log = SimpleLogger(
            level="warning", log_func=lambda typ, filename, lineno, msg: messages.append((typ, msg))
        )
        try:
            model.write_fluent(fluent_filename)
            read_fluent(fluent_filename, auto_read_write_h5=True, log=log)  # writes the h5
            assert h5_filename.exists()

            model.results[:, 0] *= 2.0
            model.write_fluent(fluent_filename)
            # make the h5 older than the sources regardless of the clock resolution
            old_mtime = os.path.getmtime(fluent_filename) - 100.0
            os.utime(h5_filename, (old_mtime, old_mtime))

            model2 = read_fluent(fluent_filename, auto_read_write_h5=True, log=log)
            assert np.array_equal(model2.results[:, 0], 20.0 * model2.result_element_id)
            assert os.path.getmtime(h5_filename) > old_mtime  # rewritten
            warnings_ = [msg for typ, msg in messages if typ == "WARNING"]
            assert len(warnings_) == 1 and "older than" in warnings_[0], messages

            # the rewritten h5 is current, so it's used without a warning
            model3 = read_fluent(fluent_filename, auto_read_write_h5=True, log=log)
            assert np.array_equal(model3.results, model2.results)
            assert len(messages) == 1, messages
        finally:
            for ext in (".cel", ".vrt", ".daten", ".h5"):
                filename = fluent_filename.with_suffix(ext)
                if filename.exists():
                    os.remove(filename)


def _three_region_model() -> Fluent:
    """3 unit quads in regions 1/2/3 that don't share nodes"""
    model = Fluent(auto_read_write_h5=False, log=SimpleLogger(level="warning"))
    model.node_id = np.arange(1, 13)
    square = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
    model.xyz = np.vstack([square, square + [2.0, 0.0, 0.0], square + [4.0, 0.0, 1.0]])
    model.quads = np.array(
        [
            [1, 1, 1, 2, 3, 4],
            [2, 2, 5, 6, 7, 8],
            [3, 3, 9, 10, 11, 12],
        ]
    )
    model.tris = np.zeros((0, 5), dtype="int32")
    return model


def _tri_quad_model() -> Fluent:
    """2 quads (1, 3) and 3 tris (2, 4, 5); results = [10*eid, -eid]"""
    model = Fluent(auto_read_write_h5=False, log=SimpleLogger(level="warning"))
    model.node_id = np.array([1, 2, 3, 4, 5, 6])
    model.xyz = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
        ]
    )
    model.quads = np.array(
        [
            [1, 7, 1, 2, 3, 4],
            [3, 7, 2, 5, 6, 3],
        ]
    )
    model.tris = np.array(
        [
            [2, 8, 1, 2, 3],
            [4, 8, 1, 3, 4],
            [5, 9, 2, 5, 6],
        ]
    )
    model.titles = np.array(["ElementID", "Pressure Coefficient", "Other"])
    model.result_element_id = np.array([1, 2, 3, 4, 5])
    model.results = np.column_stack(
        [10.0 * model.result_element_id, -1.0 * model.result_element_id]
    )
    model.element_ids = np.array([1, 2, 3, 4, 5])
    return model


def main():  # pragma: no cover
    import time

    time0 = time.time()
    unittest.main()
    print("dt = %s" % (time.time() - time0))


if __name__ == "__main__":  # pragma: no cover
    main()
