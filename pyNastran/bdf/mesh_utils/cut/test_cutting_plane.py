"""defines cutting plane tests"""
import os
from pathlib import Path
import unittest
import numpy as np
from pyNastran.gui.matplotlib_backend import matplotlib_backend

try:
    import matplotlib  # pylint: disable=unused-import
    IS_MATPLOTLIB = True
except ModuleNotFoundError:  # pragma: no cover
    IS_MATPLOTLIB = False

if IS_MATPLOTLIB:
    matplotlib.use(matplotlib_backend)

import pyNastran
from pyNastran.bdf.bdf import read_bdf, BDF, CORD2R
from cpylog import SimpleLogger

from pyNastran.bdf.mesh_utils.cut.cutting_plane_plotter import cut_and_plot_model
from pyNastran.bdf.mesh_utils.cut.cut_model_by_plane import (
    cut_edge_model_by_coord, cut_face_model_by_coord,
    connect_face_rows, split_to_trias)
#from pyNastran.bdf.mesh_utils.bdf_merge import bdf_merge
from pyNastran.op2.op2_geom import read_op2_geom

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = Path(os.path.join(PKG_PATH, '..', 'models'))


class TestCuttingPlane(unittest.TestCase):
    """various cutting plane tests"""
    def test_cut_plate(self):
        """mode 10 is a sine wave"""
        log = SimpleLogger(level='warning', encoding='utf-8')
        dirname = MODEL_PATH / 'plate_py'
        bdf_filename = dirname / 'plate_py.dat'
        op2_filename = dirname / 'plate_py.op2'
        model = read_bdf(bdf_filename, log=log)
        op2_model = read_op2_geom(op2_filename, log=log)

        title = 'Mode 10 Eigenvector'
        p1 = None
        p2 = None
        zaxis = None
        coord = CORD2R(1, rid=0, origin=[0., 0., 0.],
                       zaxis=[0., 0., 1], xzplane=[1., 0., 0.],
                       comment='')
        model.coords[1] = coord
        ytol = 2.

        # no result
        nodal_result = None
        cut_and_plot_model(title, p1, p2, zaxis,
                           model, coord, nodal_result, model.log, ytol,
                           plane_atol=1e-5, csv_filename=None, invert_yaxis=False,
                           cut_type='edge', plot=False, show=False)

        # real
        nodal_result = op2_model.eigenvectors[1].data[9, :, 2]
        cut_and_plot_model(title, p1, p2, zaxis,
                           model, coord, nodal_result, model.log, ytol,
                           plane_atol=1e-5, csv_filename='real_result.csv', invert_yaxis=False,
                           cut_type='edge', plot=IS_MATPLOTLIB, show=False)

        # complex
        nodal_result2 = np.asarray(nodal_result, dtype='complex64')
        nodal_result2.imag = -nodal_result.real
        cut_and_plot_model(
            title, p1, p2, zaxis,
            model, coord, nodal_result2, model.log, ytol,
            plane_atol=1e-5,
            csv_filename='complex_result.csv', invert_yaxis=True,
            cut_type='edge', plot=IS_MATPLOTLIB, show=False)
        dirname = Path('.')
        os.remove(dirname / 'real_result.csv')
        os.remove(dirname / 'complex_result.csv')

    def test_cut_plate_eids(self):
        """recover element ids"""
        log = SimpleLogger(level='warning', encoding='utf-8')
        bdf_filename = MODEL_PATH / 'plate_py' / 'plate_py.dat'
        model = read_bdf(bdf_filename, log=log)
        nnodes = len(model.nodes)
        nodal_result = np.ones(nnodes)

        coord = CORD2R(1, rid=0, origin=[0., 0., 0.], zaxis=[0., 0., 1], xzplane=[1., 0., 0.],
                       comment='')
        model.coords[1] = coord

        found_cut, unique_geometry_array, unique_results_array, unused_rods = cut_face_model_by_coord(
            bdf_filename, coord,
            nodal_result, plane_atol=1e-5, skip_cleanup=True,
            csv_filename='',
            plane_bdf_filename1='',
            plane_bdf_filename2='',
            face_data=None,
            debug_vectorize=False,
        )
        found_cut, unique_geometry_array, unique_results_array, unused_rods = cut_face_model_by_coord(
            bdf_filename, coord,
            nodal_result, plane_atol=1e-5, skip_cleanup=True,
            csv_filename='cut_face.csv',
            plane_bdf_filename1='plane_face.bdf',
            face_data=None,
        )
        #print(unique_geometry_array)
        #print(unique_results_array)
        unique_geometry_array = np.array(unique_geometry_array)
        unique_results_array = np.array(unique_results_array)
        assert unique_geometry_array.shape == (1, 40, 4), unique_geometry_array.shape
        assert unique_results_array.shape == (1, 40, 7), unique_results_array.shape
        unique_geometry_array = unique_geometry_array[0, :, :]
        unique_results_array = unique_results_array[0, :, :]

        assert unique_geometry_array.shape == (40, 4), unique_geometry_array.shape
        assert unique_results_array.shape == (40, 7), unique_results_array.shape
        #print(unique_geometry_array)
        #print(unique_results_array)
        os.remove('cut_face.csv')
        os.remove('plane_face.bdf')

    def test_cut_shell_model_edge_1(self):
        """
        tests:
         - cut_edge_model_by_coord
         - cut_face_model_by_coord
        """
        model, nodal_result = _cut_shell_model_quads()
        coord = CORD2R(1, rid=0, origin=[0.5, 0., 0.], zaxis=[0.5, 0., 1], xzplane=[1.5, 0., 0.],
                       comment='')
        model.coords[1] = coord
        tol = 2.
        #-------------------------------------------------------------------------
        title = 'result'
        p1 = None
        p2 = None
        zaxis = None
        cut_and_plot_model(
            title, p1, p2, zaxis,
            model, coord, nodal_result, model.log, tol,
            plane_atol=1e-5,
            csv_filename=None,
            invert_yaxis=False,
            cut_type='edge', plot=IS_MATPLOTLIB, show=False)
        #=========================================================================
        out = cut_edge_model_by_coord(
            model, coord, tol, nodal_result,
            plane_atol=1e-5)
        unused_local_points_array, unused_global_points_array, result_array = out
        assert len(result_array) == 16, len(result_array)

        found_cut, unused_geometry_array, result_array, unused_rods = cut_face_model_by_coord(
            model, coord, nodal_result,
            plane_atol=1e-5,
            face_data=None,
        )
        result_array = np.array(result_array)
        assert result_array.shape == (1, 8, 7), result_array.shape
        #os.remove('plane_edge.bdf')
        os.remove('plane_face1.bdf')

    def test_cut_shell_model_edge_2(self):
        """
        tests:
         - cut_edge_model_by_coord
         - cut_face_model_by_coord
        """
        tol = 2.
        coord = CORD2R(1, rid=0, origin=[0.5, 0., 0.],
                       zaxis=[0.5, 0., 1], xzplane=[1.5, 0., 0.],
                       comment='')
        model, nodal_result = _cut_shell_model_quads()
        #-------------------------------------------------------------------------
        # triangles
        split_to_trias(model)
        model.coords[1] = coord
        model.write_bdf('tris.bdf')

        #print('----------------------------')
        title = 'result'
        p1 = None
        p2 = None
        zaxis = None
        cut_and_plot_model(
            title, p1, p2, zaxis,
            model, coord, nodal_result, model.log, tol,
            plane_atol=1e-5,
            csv_filename=None,
            invert_yaxis=False,
            cut_type='edge', plot=IS_MATPLOTLIB, show=False)

        out = cut_edge_model_by_coord(
            model, coord, tol, nodal_result,
            plane_atol=1e-5, csv_filename='cut_edge_2.csv')
        unused_local_points_array, unused_global_points_array, result_array = out
        assert len(result_array) == 20, len(result_array)

        found_cut, unused_geometry_arrays, result_arrays, unused_rods = cut_face_model_by_coord(
            model, coord, nodal_result,
            plane_atol=1e-5,
            csv_filename='cut_face_2.csv',
            face_data=None)
        assert len(result_arrays[0]) == 8, len(result_arrays)
        os.remove('tris.bdf')
        os.remove('cut_edge_2.csv')
        os.remove('cut_face_2.csv')
        #os.remove('plane_edge.bdf')
        os.remove('plane_face1.bdf')

    def test_cut_shell_model_face_1(self):
        """
        tests:
         - cut_edge_model_by_coord
         - cut_face_model_by_coord
        """
        tol = 2.
        coord = CORD2R(1, rid=0, origin=[0.5, 0., 0.], zaxis=[0.5, 0., 1], xzplane=[1.5, 0., 0.],
                       comment='')
        model, nodal_result = _cut_shell_model_quads()
        #-------------------------------------------------------------------------
        # triangles
        split_to_trias(model)
        model.coords[1] = coord
        model.write_bdf('tris.bdf')

        #print('----------------------------')
        title = 'result'
        p1 = None
        p2 = None
        zaxis = None
        #print(nodal_result)
        with self.assertRaises(TypeError):
            cut_and_plot_model(title, p1, p2, zaxis,
                               model, coord, nodal_result, model.log, tol,
                               plane_atol=1e-5,
                               csv_filename=None,
                               invert_yaxis=False,
                               cut_type='face', plot=IS_MATPLOTLIB, show=False)
        os.remove('tris.bdf')

    def test_connect_face_rows(self):
        """in order"""
        geometry_array = np.array([
            [1, 1, 2],
            [2, 2, 3],
            [3, 3, 4],
            [4, 4, 5],
            [5, 5, 6],
        ])
        nedges = geometry_array.shape[0]
        results_array = np.arange(0, nedges)
        #print(results_array)
        iedges, unused_geometry_arrays2, unused_results_arrays2 = connect_face_rows(
            geometry_array, results_array, skip_cleanup=False)
        assert np.array_equal(iedges, [[0, 1, 2, 3, 4]]), 'iedges=%s' % iedges
        #-----------------------------------------------------------------------

        # out of order
        geometry_array = np.array([
            [1, 1, 2], # 0
            [2, 4, 5], # 3
            [3, 5, 6], # 4
            [4, 3, 4], # 2
            [5, 2, 3], # 1
        ])
        nedges = geometry_array.shape[0]
        results_array = np.arange(0, nedges)
        iedges, unused_geometry_arrays2, unused_results_arrays2 = connect_face_rows(
            geometry_array, results_array, skip_cleanup=False)
        assert np.array_equal(iedges, [[0, 4, 3, 1, 2]]), 'iedges=%s' % iedges
        #print(geometry_array2)

        #-----------------------------------------------------------------------

        # in order, two blocks
        #print('*****************')
        geometry_array = np.array([
            # block 1
            [1, 1, 2],
            [2, 2, 3],
            [3, 3, 4],

            # block 2
            [10, 10, 20],
            [20, 20, 30],
            [30, 30, 40],
        ])
        nedges = geometry_array.shape[0]
        results_array = np.arange(0, nedges)
        #print(results_array)
        iedges, unused_geometry_array2, unused_results_array2 = connect_face_rows(
            geometry_array, results_array, skip_cleanup=False)
        assert np.array_equal(iedges, [[0, 1, 2], [3, 4, 5]]), 'iedges=%s' % iedges

    def test_connect_face_rows_ring_1(self):
        """in order, one ring"""
        geometry_array = np.array([
            [1, 1, 2],
            [2, 2, 3],
            [3, 3, 4],
            [4, 1, 4],
        ])
        nedges = geometry_array.shape[0]
        results_array = np.arange(0, nedges)
        #print(results_array)
        iedges, unused_geometry_array2, unused_results_array2 = connect_face_rows(
            geometry_array, results_array, skip_cleanup=False)
        assert np.array_equal(iedges, [[0, 1, 2, 3, 0]]), 'iedges=%s' % iedges

    def test_connect_face_rows_ring_2(self):
        """in order, two rings"""
        geometry_array = np.array([
            [1, 1, 2],
            [2, 2, 3],
            [3, 3, 4],
            [4, 1, 4],

            [10, 10, 20],
            [20, 20, 30],
            [30, 30, 40],
            [40, 10, 40],
        ])
        nedges = geometry_array.shape[0]
        results_array = np.arange(0, nedges)
        #print(results_array)
        iedges, unused_geometry_array2, unused_results_array2 = connect_face_rows(
            geometry_array, results_array, skip_cleanup=False)
        assert np.array_equal(iedges, [[0, 1, 2, 3, 0], [4, 5, 6, 7, 4]]), 'iedges=%s' % iedges


def _cut_shell_model_quads():
    """helper method"""
    log = SimpleLogger(level='error')
    pid = 10
    mid1 = 100
    model = BDF(log=log)

    # intersects (min)
    model.add_grid(1, [0., 0., 0.])
    model.add_grid(2, [1., 0., 0.])
    model.add_grid(3, [1., 1., 0.])
    model.add_grid(4, [0., 1., 0.])
    model.add_cquad4(1, pid, [1, 2, 3, 4])

    # intersects (max)
    model.add_grid(5, [0., 0., 1.])
    model.add_grid(6, [1., 0., 1.])
    model.add_grid(7, [1., 1., 1.])
    model.add_grid(8, [0., 1., 1.])
    model.add_cquad4(2, pid, [5, 6, 7, 8])

    # intersects (mid)
    model.add_grid(9, [0., 0., 0.5])
    model.add_grid(10, [1., 0., 0.5])
    model.add_grid(11, [1., 1., 0.5])
    model.add_grid(12, [0., 1., 0.5])
    model.add_cquad4(3, pid, [9, 10, 11, 12])

    # doesn't intersect
    model.add_grid(13, [10., 0., 0.])
    model.add_grid(14, [11., 0., 0.])
    model.add_grid(15, [11., 1., 0.])
    model.add_grid(16, [10., 1., 0.])
    model.add_cquad4(4, pid, [13, 14, 15, 16])

    model.add_pshell(pid, mid1=mid1, t=2.)
    E = 1.0
    G = None
    nu = 0.3
    model.add_mat1(mid1, E, G, nu, rho=1.0)
    model.validate()

    model.cross_reference()

    #xyz_points = [
        #[0.4, 0.6, 0.], [-1., -1, 0.],]
    nodal_result = np.linspace(0., 1., num=16)
    return model, nodal_result


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
