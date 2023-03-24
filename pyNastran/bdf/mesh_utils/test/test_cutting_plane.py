"""defines cutting plane tests"""
import os
from itertools import count
from typing import Any
import unittest
import numpy as np
#import PySide
from pyNastran.gui.matplotlib_backend import matplotlib_backend

try:
    import matplotlib  # pylint: disable=unused-import
    IS_MATPLOTLIB = True
except ImportError:  # pragma: no cover
    IS_MATPLOTLIB = False

if IS_MATPLOTLIB:
    matplotlib.use(matplotlib_backend)
    import matplotlib.pyplot as plt  # pylint: disable=unused-import

import pyNastran
from pyNastran.bdf.bdf import read_bdf, BDF, CORD2R
from cpylog import SimpleLogger

from pyNastran.bdf.mesh_utils.cut_model_by_plane import (
    cut_edge_model_by_coord, cut_face_model_by_coord, connect_face_rows,
    split_to_trias, calculate_area_moi, _get_shell_inertia)
from pyNastran.bdf.mesh_utils.cutting_plane_plotter import cut_and_plot_model
#from pyNastran.bdf.mesh_utils.bdf_merge import bdf_merge
from pyNastran.op2.op2_geom import read_op2_geom

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, '..', 'models')


class TestCuttingPlane(unittest.TestCase):
    """various cutting plane tests"""
    def test_shell_inertia(self):
        model = BDF(debug=False, log=None, mode='msc')

        # x-axis is at 0 degrees
        cid1 = 1
        origin = [0., 0., 0.]
        zaxis = [0., 0., 1.]
        xzplane = [0., 1., 0.]
        model.add_cord2r(cid1, origin, zaxis, xzplane)

        # x-axis is at 90 degrees
        cid2 = 2
        origin = [0., 0., 0.]
        zaxis = [0., 0., 1.]
        xzplane = [1., 0., 0.]
        model.add_cord2r(cid2, origin, zaxis, xzplane)

        # x-axis is at 45 degrees
        cid3 = 3
        origin = [0., 0., 0.]
        zaxis = [0., 0., 1.]
        xzplane = [1., 1., 0.]
        model.add_cord2r(cid3, origin, zaxis, xzplane)


        tply = 0.007
        thicknesses = [tply, tply]
        mids_ud = [1, 1]
        mids_45 = [1, 1]
        model.add_pcomp(1, [mids_ud[0]], [thicknesses[0]], thetas=[0.])
        model.add_pcomp(2, [mids_ud[0]], [thicknesses[0]], thetas=[90.])
        model.add_pcomp(3, [mids_ud[0]], [thicknesses[0]], thetas=[45.])

        model.add_pcomp(4, mids_ud, thicknesses, thetas=[45., -45.])
        model.add_pcomp(11, mids_45, thicknesses, thetas=[0., 90.])
        model.add_pcomp(12, mids_45, thicknesses, thetas=[45., -45.])

        # fibers are in the y direction
        nids = [1, 2, 3, 4]
        element1 = model.add_cquad4(1, 1, nids, theta_mcid=cid1)
        element2 = model.add_cquad4(2, 2, nids, theta_mcid=cid2)
        element3 = model.add_cquad4(3, 3, nids, theta_mcid=cid3)
        element4 = model.add_cquad4(4, 4, nids, theta_mcid=cid3)

        element11 = model.add_cquad4(11, 11, nids, theta_mcid=cid1)
        element12 = model.add_cquad4(12, 12, nids, theta_mcid=cid3)

        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])

        # fabric
        mid = 1
        e11 = 7600555.0
        e22 = 7029000.0
        nu12 = 0.042
        g12 = 360471.0
        model.add_mat8(mid, e11, e22, nu12, g12=g12, g1z=1e8, g2z=1e8,
                       rho=0., a1=0., a2=0., tref=0., Xt=0., Xc=None, Yt=0., Yc=None,
                       S=0., ge=0., F12=0., strn=0., comment='')
        lengthi = 1.
        normal_plane = np.array([0., 1., 0.])
        normal_plane_vector = normal_plane.copy().reshape((3, 1))

        model.cross_reference()

        # not rotated: 0 deg
        thicknessi, areai, imat_rotation_angle_deg, Ex, Ey, Gxy, nu_xy = _get_shell_inertia(
            element1, normal_plane, normal_plane_vector, lengthi)
        assert np.allclose(thicknessi, tply)
        assert np.allclose(areai, thicknessi*lengthi)
        assert np.allclose(imat_rotation_angle_deg, 0.)
        assert np.allclose(e11, Ex)
        assert np.allclose(e22, Ey)
        assert np.allclose(g12, Gxy)
        assert np.allclose(nu_xy, 0.)

        # rotate by 90 degrees: 90 deg
        thicknessi, areai, imat_rotation_angle_deg, Ex, Ey, Gxy, nu_xy = _get_shell_inertia(
            element2, normal_plane, normal_plane_vector, lengthi)
        assert np.allclose(thicknessi, tply)
        assert np.allclose(areai, thicknessi*lengthi)
        assert np.allclose(imat_rotation_angle_deg, 90.)
        #assert np.allclose(e11, Ey)
        #assert np.allclose(e22, Ex)
        #assert np.allclose(g12, Gxy)
        #assert np.allclose(nu_xy, 0.)

        # rotate by 45 degrees: +45 deg
        thicknessi, areai, imat_rotation_angle_deg, Ex, Ey, Gxy, nu_xy = _get_shell_inertia(
            element3, normal_plane, normal_plane_vector, lengthi)
        assert np.allclose(thicknessi, tply)
        assert np.allclose(areai, thicknessi*lengthi)
        assert np.allclose(imat_rotation_angle_deg, 45.)
        #assert np.allclose(Ex, 1317118.060260035)
        #assert np.allclose(Ey, 1317118.0602600349)
        #assert np.allclose(Gxy, 3510140.123933055)
        #assert np.allclose(nu_xy, 4.766776571060659e-13)

        # rotate by 45 degrees: +/-45 deg
        thicknessi, areai, imat_rotation_angle_deg, Ex, Ey, Gxy, nu_xy = _get_shell_inertia(
            element4, normal_plane, normal_plane_vector, lengthi)
        assert np.allclose(thicknessi, sum(thicknesses))
        assert np.allclose(areai, thicknessi*lengthi)
        assert np.allclose(imat_rotation_angle_deg, 45.)
        #assert np.allclose(Ex, 1317292.3250658484)
        #assert np.allclose(Ey, 1317292.3250658484)
        #assert np.allclose(Gxy, 3515514.7806900046)
        #assert np.allclose(nu_xy, 4.766908439759332e-13)

        # fabric - rotate by 45 degrees: +/- 0/90 deg
        thicknessi, areai, imat_rotation_angle_deg, Ex, Ey, Gxy, nu_xy = _get_shell_inertia(
            element11, normal_plane, normal_plane_vector, lengthi)
        assert np.allclose(thicknessi, sum(thicknesses))
        assert np.allclose(areai, thicknessi*lengthi)
        assert np.allclose(imat_rotation_angle_deg, 0.)
        #assert np.allclose(Ex, e11)
        #assert np.allclose(Ey, e22)
        #assert np.allclose(Gxy, g12)
        #assert np.allclose(nu_xy, nu12)

        # fabric - rotate by 45 degrees: +/- 45 deg
        #thicknessi, areai, imat_rotation_angle_deg, Ex, Ey, Gxy, nu_xy = _get_shell_inertia(
            #element12, normal_plane, normal_plane_vector, lengthi)
        #assert np.allclose(thicknessi, sum(thicknesses))
        #assert np.allclose(areai, thicknessi*lengthi)
        #assert np.allclose(imat_rotation_angle_deg, 45.)
        #assert np.allclose(Ex, 1317292.3250658484)
        #assert np.allclose(Ey, 1317292.3250658484)
        #assert np.allclose(Gxy, 3515514.7806900046)
        #assert np.allclose(nu_xy, 4.766908439759332e-13)

        x = 1

    def test_cut_plate(self):
        """mode 10 is a sine wave"""
        log = SimpleLogger(level='warning', encoding='utf-8', log_func=None)
        bdf_filename = os.path.join(MODEL_PATH, 'plate_py', 'plate_py.dat')
        op2_filename = os.path.join(MODEL_PATH, 'plate_py', 'plate_py.op2')
        model = read_bdf(bdf_filename, log=log)
        op2_model = read_op2_geom(op2_filename, log=log)

        title = 'Mode 10 Eigenvector'
        p1 = None
        p2 = None
        zaxis = None
        coord = CORD2R(1, rid=0, origin=[0., 0., 0.], zaxis=[0., 0., 1], xzplane=[1., 0., 0.],
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
        cut_and_plot_model(title, p1, p2, zaxis,
                           model, coord, nodal_result2, model.log, ytol,
                           plane_atol=1e-5, csv_filename='complex_result.csv', invert_yaxis=True,
                           cut_type='edge', plot=IS_MATPLOTLIB, show=False)
        os.remove('real_result.csv')
        os.remove('complex_result.csv')

    def _test_cut_box(self):  # pragma: no cover
        """recover element ids"""
        log = SimpleLogger(level='warning', encoding='utf-8', log_func=None)
        #bdf_filename = r'SEction_1_box.bdf'  # x-axis
        #normal_plane = np.array([1., 0., 0.])
        bdf_filename = 'Section_1_box_4.bdf'  # y-axis
        normal_plane = np.array([0., 1., 0.])
        dys, coords = get_coords_box()
        y, A, I, J, EI, GJ, avg_centroid, plane_bdf_filenames, plane_bdf_filenames2 = cut_and_plot_moi(
            bdf_filename, normal_plane, log,
            dys, coords,
            ytol=0.0,
            plot=True, show=True)

        show = True
        if IS_MATPLOTLIB:
            plot_inertia(y, A, I, J, EI, GJ, avg_centroid, show=show)

    def test_cut_bwb(self):
        """recover element ids"""
        log = SimpleLogger(level='warning', encoding='utf-8', log_func=None)
        is_bwb = True
        if is_bwb:
            bdf_filename = os.path.join(MODEL_PATH, 'bwb', 'bwb_saero.bdf')  # ymax~=1262.0
            dys, coords = get_coords_bwb()
        else:  #  pragma: no cover
            bdf_filename = r'C:\NASA\asm\all_modes_mach_0.85\flutter.bdf'  # ymax=1160.601
            dys, coords = get_coords_crm()
        normal_plane = np.array([0., 1., 0.])

        y, A, I, J, EI, J, avg_centroid, plane_bdf_filenames, plane_bdf_filenames2 = cut_and_plot_moi(
            bdf_filename, normal_plane, log,
            dys, coords,
            ytol=2.0,
            plot=True, show=True)

        show = True
        #show = False
        if IS_MATPLOTLIB:
            GJ = J
            plot_inertia(y, A, I, J, EI, GJ, avg_centroid, show=show)
            os.remove('normalized_inertia_vs_span.png')
            os.remove('area_vs_span.png')
            os.remove('amoi_vs_span.png')
            os.remove('e_amoi_vs_span.png')
            os.remove('cg_vs_span.png')

        #bdf_merge(plane_bdf_filenames, bdf_filename_out='merge.bdf', renumber=True,
                  #encoding=None, size=8, is_double=False, cards_to_skip=None,
                  #log=None, skip_case_control_deck=False)
        for plane_bdf_filename in plane_bdf_filenames:
            os.remove(plane_bdf_filename)
        os.remove('thetas.csv')
        os.remove('equivalent_beam_model.bdf')
        os.remove('cut_data_vs_span.csv')
        #os.remove('cut_face.csv')
        #if IS_MATPLOTLIB:
            #os.remove('area_vs_span.png')
            #os.remove('amoi_vs_span.png')
            #os.remove('normalized_inertia_vs_span.png')
            #os.remove('cg_vs_span.png')
            #os.remove('e_amoi_vs_span.png')

    def test_cut_plate_eids(self):
        """recover element ids"""
        log = SimpleLogger(level='warning', encoding='utf-8', log_func=None)
        bdf_filename = os.path.join(MODEL_PATH, 'plate_py', 'plate_py.dat')
        model = read_bdf(bdf_filename, log=log)
        nnodes = len(model.nodes)
        nodal_result = np.ones(nnodes)

        coord = CORD2R(1, rid=0, origin=[0., 0., 0.], zaxis=[0., 0., 1], xzplane=[1., 0., 0.],
                       comment='')
        model.coords[1] = coord
        ytol = 2.

        unique_geometry_array, unique_results_array, unused_rods = cut_face_model_by_coord(
            bdf_filename, coord, ytol,
            nodal_result, plane_atol=1e-5, skip_cleanup=True,
            csv_filename='cut_face.csv',
            plane_bdf_filename='plane_face.bdf',
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
        cut_and_plot_model(title, p1, p2, zaxis,
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

        unused_geometry_array, result_array, unused_rods = cut_face_model_by_coord(
            model, coord, tol, nodal_result,
            plane_atol=1e-5)
        result_array = np.array(result_array)
        assert result_array.shape == (1, 8, 7), result_array.shape
        #os.remove('plane_edge.bdf')
        os.remove('plane_face.bdf')

    def test_cut_shell_model_edge_2(self):
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
        cut_and_plot_model(title, p1, p2, zaxis,
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

        unused_geometry_arrays, result_arrays, unused_rods = cut_face_model_by_coord(
            model, coord, tol, nodal_result,
            plane_atol=1e-5, csv_filename='cut_face_2.csv')
        assert len(result_arrays[0]) == 8, len(result_arrays)
        os.remove('tris.bdf')
        os.remove('cut_edge_2.csv')
        os.remove('cut_face_2.csv')
        #os.remove('plane_edge.bdf')
        os.remove('plane_face.bdf')

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


def get_coords_bwb(ncuts: int=2000) -> tuple[list[float], list[CORD2R]]:  # pragma: no cover
    dys = []
    coords = []
    for i in range(ncuts):
        dy = 100. * i + 1.  #  bwb
        coord = CORD2R(1, rid=0, origin=[0., dy, 0.], zaxis=[0., dy, 1], xzplane=[1., dy, 0.])
        dys.append(dy)
        coords.append(coord)
    return dys, coords

def get_coords_crm(ncuts: int=2000) -> tuple[list[float], list[CORD2R]]:  # pragma: no cover
    dys = []
    coords = []
    for i in range(ncuts):
        dy = 4. * i + 1.  #  CRM
        coord = CORD2R(1, rid=0, origin=[0., dy, 0.], zaxis=[0., dy, 1], xzplane=[1., dy, 0.])
        dys.append(dy)
        coords.append(coord)
    return dys, coords

def get_coords_box(ncuts: int) -> tuple[list[float], list[CORD2R]]:  # pragma: no cover
    dys = []
    coords = []
    for i in range(ncuts):
        dy = -0.1 * i - 0.1  #  box
        #coord = CORD2R(1, rid=0, origin=[0., dy, 0.], zaxis=[0., dy, 1], xzplane=[1., dy, 0.])
        coord = CORD2R(1, rid=0, origin=[0., dy, 0.], zaxis=[0., dy, 1], xzplane=[1., dy, 0.])

        #if dy < -5:
            #print('break', dy)
            #break

        #origin = np.array([0., dy, 0.])
        #xzplane = origin + dx
        #xzplane = np.array([1., dy, 0.])
        #coord = CORD2R.add_axes(cid, rid=0, origin=p1, xaxis=p2-p1, yaxis=None, zaxis=None,
                                 #xyplane=None, yzplane=None, xzplane=None, comment='')
        #print(coord)
        dys.append(dy)
        coords.append(coord)
    return dys, coords

def cut_and_plot_moi(bdf_filename: str, normal_plane: np.ndarray, log: SimpleLogger,
                     dys: list[float],
                     coords: list[CORD2R],
                     ytol: float=2.0,
                     dirname: str='',
                     plot: bool=True, show: bool=False) -> tuple[Any, Any, Any, Any, Any]: # y, A, I, EI, avg_centroid
    model = read_bdf(bdf_filename, log=log)
    model2 = read_bdf(bdf_filename, log=log)

    out = _get_station_data(
        model, model2,
        dys, coords, normal_plane,
        ytol, dirname)
    thetas, y, dx, dz, A, I, J, EI, GJ, avg_centroid, plane_bdf_filenames, plane_bdf_filenames2 = out

    assert len(y) > 0, y
    thetas_csv_filename = os.path.join(dirname, 'thetas.csv')

    with open(thetas_csv_filename, 'w') as csv_filename:
        csv_filename.write('# eid(%i),theta,Ex,Ey,Gxy\n')
        for eid, (theta, Ex, Ey, Gxy) in sorted(thetas.items()):
            csv_filename.write('%d,%f,%f,%f,%f\n' % (eid, theta, Ex, Ey, Gxy))

    inid = 1
    beam_model = BDF(debug=False)
    avg_centroid[:, 1] = y

    # wrong
    mid = 1
    E = 3.0e7
    G = None
    nu = 0.3
    model.add_mat1(mid, E, G, nu, rho=0.1)

    #   0    1    2    3    4    5
    # [Ixx, Iyy, Izz, Ixy, Iyz, Ixz]
    Ix = I[:, 0]
    Iy = I[:, 1]
    Iz = I[:, 2]
    Ixz = I[:, 5]

    ExIx = EI[:, 0]
    ExIy = EI[:, 1]
    ExIz = EI[:, 2]
    ExIxz = EI[:, 5]

    J = Ix + Iz
    #i1, i2, i12 = Ix, Iy, Ixy
    for inid, xyz in enumerate(avg_centroid):
        beam_model.add_grid(inid+1, xyz)
    for eid in range(1, len(A)):
        pid = eid
        nids = [eid, eid + 1]
        x = [1., 0., 0.]
        g0 = None
        beam_model.add_cbeam(eid, pid, nids, x, g0, offt='GGG', bit=None,
                             pa=0, pb=0, wa=None, wb=None, sa=0, sb=0, comment='')

        # j = i1 + i2
        so = ['YES', 'YES']
        xxb = [0., 1.]
        area = [A[eid-1], A[eid]]
        i1 = [Ix[eid-1], Ix[eid]]
        i2 = [Iz[eid-1], Iz[eid]]
        i12 = [Ixz[eid-1], Ixz[eid]]
        j = [J[eid-1], J[eid]]
        beam_model.add_pbeam(pid, mid, xxb, so, area, i1, i2, i12, j, nsm=None,
                             c1=None, c2=None, d1=None, d2=None, e1=None, e2=None, f1=None, f2=None,
                             k1=1., k2=1., s1=0., s2=0., nsia=0., nsib=None, cwa=0., cwb=None,
                             m1a=0., m2a=0., m1b=None, m2b=None,
                             n1a=0., n2a=0., n1b=None, n2b=None,
                             comment='')

    beam_model_bdf_filename = os.path.join(dirname, 'equivalent_beam_model.bdf')
    beam_model.write_bdf(beam_model_bdf_filename)

    X = np.vstack([y, dx, dz, A, Ix, Iz, Ixz, ExIx, ExIz, ExIxz]).T
    Y = np.hstack([X, avg_centroid])
    header = 'y, dx, dz, A, Ix, Iz, Ixz, Ex*Ix, Ex*Iz, Ex*Ixz, xcentroid, ycentroid, zcentroid'
    cut_data_span_filename = os.path.join(dirname, 'cut_data_vs_span.csv')
    np.savetxt(cut_data_span_filename, Y, header=header, delimiter=',')

    plot_inertia(y, A, I, J, EI, GJ, avg_centroid, show=show, dirname=dirname)
    return y, A, I, J, EI, GJ, avg_centroid, plane_bdf_filenames, plane_bdf_filenames2

def _get_station_data(model: BDF, model2: BDF,
                     dys, coords, normal_plane: np.ndarray,
                     ytol: float, dirname: str, ) -> tuple[
                         dict[int, tuple[float, float, float, float]],  # thetas
                         #y, dx, dz,
                         #A, I, J,
                         #EI, GJ, avg_centroid
                         Any, Any, Any,
                         Any, Any, Any,
                         Any, Any, Any,
                         #plane_bdf_filenames, plane_bdf_filenames2,
                         list[str], list[str]]:
    # initialize theta
    thetas = {}
    for eid in model.elements:
        #  theta, Ex, Ey, Gxy
        thetas[eid] = (0., 0., 0., 0.)

    #p1 = np.array([466.78845, 735.9053, 0.0])
    #p2 = np.array([624.91345, 639.68896, -0.99763656])
    #dx = p2 - p1
    nodal_result = None
    plane_bdf_filenames = []
    plane_bdf_filenames2 = []
    y = []
    dx = []
    dz = []
    A = []
    I = []
    J = []
    EI = []
    GJ = []
    avg_centroid = []

    assert len(dys) > 0, dys
    for i, dy, coord in zip(count(), dys, coords):
        model.coords[1] = coord
        plane_bdf_filename = os.path.join(dirname, f'plane_face_{i:d}.bdf')
        plane_bdf_filename2 = os.path.join(dirname, f'plane_face2_{i:d}.bdf')
        cut_face_filename = os.path.join(dirname, f'cut_face_{i:d}.csv')
        if os.path.exists(cut_face_filename):
            os.remove(cut_face_filename)
        try:
            out = cut_face_model_by_coord(
                model2, coord, ytol,
                nodal_result, plane_atol=1e-5, skip_cleanup=True,
                #csv_filename=cut_face_filename,
                csv_filename=None,
                #plane_bdf_filename=None)
                plane_bdf_filename=plane_bdf_filename,
                plane_bdf_filename2=plane_bdf_filename2,
                plane_bdf_offset=dy)
        except PermissionError:
            print(f'failed to delete {plane_bdf_filename}')
            continue
        except RuntimeError:
            # incorrect ivalues=[0, 1, 2]; dy=771. for CRM
            continue
        unused_unique_geometry_array, unused_unique_results_array, rods = out

        if not os.path.exists(plane_bdf_filename):
            break
        plane_bdf_filenames.append(plane_bdf_filename)
        plane_bdf_filenames2.append(plane_bdf_filename2)
        # eid, nid, inid1, inid2
        #print(unique_geometry_array)
        #moi_filename = 'amoi_%i.bdf' % i
        moi_filename = None
        dxi, dzi, Ai, Ii, EIi, avg_centroidi = calculate_area_moi(
            model, rods, normal_plane, thetas, moi_filename=moi_filename)

        #print(out)
        Ji = GJi = 1.0
        y.append(dy)
        dx.append(dxi)  # length
        dz.append(dzi)  # height
        A.append(Ai)
        I.append(Ii)
        J.append(Ji)
        EI.append(EIi)
        GJ.append(GJi)
        avg_centroid.append(avg_centroidi)
        #break
    assert len(y) > 0, y

    y = np.array(y, dtype='float64')
    A = np.array(A, dtype='float64')
    dx = np.array(dx, dtype='float64')
    dz = np.array(dz, dtype='float64')
    I = np.array(I, dtype='float64')
    J = np.array(J, dtype='float64')
    EI = np.array(EI, dtype='float64')
    GJ = np.array(GJ, dtype='float64')
    avg_centroid = np.array(avg_centroid, dtype='float64')
    return thetas, y, dx, dz, A, I, J, EI, GJ, avg_centroid, plane_bdf_filenames, plane_bdf_filenames2


def plot_inertia(y, A, I, J, EI, GJ, avg_centroid, ifig: int=1, show: bool=True, dirname: str=''):
    """helper method for test"""
    #plt.plot(y, I[:, 0] / I[:, 0].max(), 'ro-', label='Qxx')
    #plt.plot(y, I[:, 1] / I[:, 1].max(), 'bo-', label='Qyy')
    #plt.plot(y, I[:, 2] / I[:, 2].max(), 'go-', label='Qxy')
    aI = np.abs(I)
    aEI = np.abs(EI)
    aGJ = np.abs(GJ)

    fig = plt.figure(ifig)
    ax = fig.gca()
    ax.plot(y, I[:, 0] / aI[:, 0].max(), 'ro-', label='Ixx')
    ax.plot(y, I[:, 1] / aI[:, 1].max(), 'bo-', label='Izz')
    ax.plot(y, I[:, 2] / aI[:, 2].max(), 'go-', label='Ixz')

    ax.plot(y, EI[:, 0] / aEI[:, 0].max(), 'ro', label='EIxx', linestyle='--')
    ax.plot(y, EI[:, 1] / aEI[:, 1].max(), 'bo', label='EIzz', linestyle='--')
    ax.plot(y, EI[:, 2] / aEI[:, 2].max(), 'go', label='EIxz', linestyle='--')
    #ax.plot(y, GJ / aGJ.max(), 'go-', label='GJ', linestyle='--')

    ax.grid(True)
    ax.set_xlabel('Span, y')
    ax.set_ylabel('Normalized Area MOI, I')
    ax.legend()
    fig.savefig('normalized_inertia_vs_span.png')
    #-------------------------------------------------------

    fig = plt.figure(ifig + 1)
    ax = fig.gca()
    ax.plot(y, A, 'ro', label='Area', linestyle='-')

    ax.grid(True)
    ax.set_xlabel('Span, y')
    ax.set_ylabel('Area, A')
    ax.legend()
    fig.savefig('area_vs_span.png')
    #-------------------------------------------------------

    fig = plt.figure(ifig + 2)
    ax = fig.gca()
    ax.plot(y, I[:, 0], 'ro-', label='Ixx')
    ax.plot(y, I[:, 1], 'bo-', label='Izz')
    ax.plot(y, I[:, 2], 'go-', label='Ixz')
    ax.grid(True)
    ax.set_xlabel('Span, y')
    ax.set_ylabel('Area MOI, I')
    ax.legend()
    fig.savefig('amoi_vs_span.png')
    #-------------------------------------------------------


    fig = plt.figure(ifig + 3)
    ax = fig.gca()
    ax.plot(y, EI[:, 0], 'ro-', label='EIxx')
    #ax.plot(y, I[:, 0], 'bo-', label='Ixx')
    ax.grid(True)
    ax.set_xlabel('Span, y')
    ax.set_ylabel('Exx*Area MOI, Exx*I')
    ax.legend()
    fig.savefig('e_amoi_vs_span.png')
    #-------------------------------------------------------

    fig = plt.figure(ifig + 4)
    ax = fig.gca()
    ax.plot(y, avg_centroid[:, 0], 'ro-', label='xcg')
    ax.plot(y, avg_centroid[:, 2], 'bo-', label='zcg')
    ax.grid(True)
    ax.set_xlabel('Span, y')
    ax.set_ylabel('CG')
    ax.legend()
    fig.savefig('cg_vs_span.png')
    #-------------------------------------------------------

    if show:
        plt.show()
    ifig += 4
    return ifig

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

    #tol = 2.
    nodal_result = np.linspace(0., 1., num=16)
    return model, nodal_result


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
