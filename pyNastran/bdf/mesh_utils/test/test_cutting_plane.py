"""defines cutting plane tests"""
import os
import copy
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

from pyNastran.bdf.mesh_utils.moi_plotter import (
    cut_and_plot_moi, plot_inertia)
from pyNastran.bdf.mesh_utils.cut_model_by_plane import (
    cut_edge_model_by_coord, cut_face_model_by_coord,
    connect_face_rows, split_to_trias,
    _get_shell_inertia, _setup_faces,
)
from pyNastran.bdf.mesh_utils.cutting_plane_plotter import cut_and_plot_model
#from pyNastran.bdf.mesh_utils.bdf_merge import bdf_merge
from pyNastran.op2.op2_geom import read_op2_geom

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = Path(os.path.join(PKG_PATH, '..', 'models'))


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
        # assert np.allclose(Ex, e11)
        # assert np.allclose(Ey, e22)
        # assert np.allclose(Gxy, g12)
        # assert np.allclose(nu_xy, nu12)

        # fabric - rotate by 45 degrees: +/- 45 deg
        # thicknessi, areai, imat_rotation_angle_deg, Ex, Ey, Gxy, nu_xy = _get_shell_inertia(
        #     element12, normal_plane, normal_plane_vector, lengthi)
        # assert np.allclose(thicknessi, sum(thicknesses))
        # assert np.allclose(areai, thicknessi*lengthi)
        # assert np.allclose(imat_rotation_angle_deg, 45.)
        # assert np.allclose(Ex, 1317292.3250658484)
        # assert np.allclose(Ey, 1317292.3250658484)
        # assert np.allclose(Gxy, 3515514.7806900046)
        # assert np.allclose(nu_xy, 4.766908439759332e-13)
        x = 1

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

    def test_cut_bwb(self):
        """cut_and_plot_moi"""
        log = SimpleLogger(level='warning', encoding='utf-8')
        dirname = MODEL_PATH / 'bwb'
        bdf_filename = dirname / 'bwb_saero.bdf'  # ymax~=1262.0
        model = read_bdf(bdf_filename, log=log)
        dys, coords = get_coords_bwb()
        normal_plane = np.array([0., 1., 0.])

        face_data = _setup_faces(model)
        # y0, A0, I0, J0, EI0, J0, avg_centroid0, plane_bdf_filenames10, plane_bdf_filenames20 = cut_and_plot_moi(
        #     model, normal_plane, log,
        #     dys, coords,
        #     ytol=2.0,
        #     dirname=dirname,
        #     plot=False, show=False, face_data=face_data)

        y, A, I, J, EI, J, avg_centroid, plane_bdf_filenames1, plane_bdf_filenames2 = cut_and_plot_moi(
            bdf_filename, normal_plane, log,
            dys, coords,
            ytol=2.0,
            dirname=dirname,
            plot=True, show=False, face_data=face_data)
        # assert np.allclose(avg_centroid, avg_centroid0)

        # show = True
        show = False
        if IS_MATPLOTLIB:
            GJ = J
            plot_inertia(y, A, I, J, EI, GJ, avg_centroid, show=show)
            os.remove(dirname / 'normalized_inertia_vs_span.png')
            os.remove(dirname / 'area_vs_span.png')
            os.remove(dirname / 'amoi_vs_span.png')
            os.remove(dirname / 'e_amoi_vs_span.png')
            os.remove(dirname / 'cg_vs_span.png')

        # bdf_merge(plane_bdf_filenames, bdf_filename_out='merge.bdf', renumber=True,
        #           encoding=None, size=8, is_double=False, cards_to_skip=None,
        #           log=None, skip_case_control_deck=False)
        for plane_bdf_filename in plane_bdf_filenames1:
            os.remove(plane_bdf_filename)
        os.remove(dirname / 'thetas.csv')
        # os.remove(dirname / 'equivalent_beam_model.bdf')
        os.remove(dirname / 'cut_data_vs_span.csv')
        # os.remove('cut_face.csv')
        # if IS_MATPLOTLIB:
        #     os.remove('area_vs_span.png')
        #     os.remove('amoi_vs_span.png')
        #     os.remove('normalized_inertia_vs_span.png')
        #     os.remove('cg_vs_span.png')
        #     os.remove('e_amoi_vs_span.png')

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
        ytol = 2.

        unique_geometry_array, unique_results_array, unused_rods = cut_face_model_by_coord(
            bdf_filename, coord, ytol,
            nodal_result, plane_atol=1e-5, skip_cleanup=True,
            csv_filename='',
            plane_bdf_filename1='',
            plane_bdf_filename2='',
            face_data=None,
            debug_vectorize=False,
        )
        unique_geometry_array, unique_results_array, unused_rods = cut_face_model_by_coord(
            bdf_filename, coord, ytol,
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


def get_coords_bwb(ncuts: int=2000) -> tuple[list[float], list[CORD2R]]:  # pragma: no cover
    """gets coords from y=0 to y=100*ncuts"""
    dys = []
    coords = []
    for i in range(ncuts):
        dy = 100. * i + 1.  #  bwb
        coord = CORD2R(1, rid=0, origin=[0., dy, 0.], zaxis=[0., dy, 1], xzplane=[1., dy, 0.])
        dys.append(dy)
        coords.append(coord)
    return dys, coords

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
