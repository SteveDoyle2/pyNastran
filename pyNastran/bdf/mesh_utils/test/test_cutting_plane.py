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
        # show = True
        show = False

        # log = SimpleLogger(level='warning', encoding='utf-8')
        log = SimpleLogger(level='debug', encoding='utf-8')
        dirname = MODEL_PATH / 'bwb'
        bdf_filename = dirname / 'bwb_saero.bdf'  # ymax~=1262.0
        model = read_bdf(bdf_filename, log=log)

        ymax = 1401.
        ncut = 10
        dy = (ymax-1.) / ncut
        i = np.arange(ncut, dtype='int32')
        # xstations = 50 * i + 1.
        ystations = dy * i + 1.

        origin = [0., 0., 0.]
        zaxis = [0., 0., 1.]
        xzplane = [1., 0., 0.]
        coords = get_coords_bwb(
            ystations, cid=1,
            origin=origin,
            zaxis=zaxis,
            xzplane=xzplane)

        coord0 = coords[0]
        normal_plane = coord0.j
        # normal_plane = np.array([0., 1., 0.])
        # assert np.allclose(normal_plane, normal_plane2)

        face_data = _setup_faces(model)
        # nids, xyz_cid0, elements = face_data
        # y0, A0, I0, J0, EI0, J0, avg_centroid0, plane_bdf_filenames10, plane_bdf_filenames20 = cut_and_plot_moi(
        #     model, normal_plane, log,
        #     dys, coords,
        #     ytol=2.0,
        #     dirname=dirname,
        #     plot=False, show=False, face_data=face_data)

        moi_data = cut_and_plot_moi(
            bdf_filename, normal_plane, log,
            ystations, coords,
            ytol=2.0,
            dirname=dirname,
            plot=True, show=False, face_data=face_data,
            debug_vectorize=True,
        )
        # y, A, I, J, EI, GJ, avg_centroid, plane_bdf_filenames, plane_bdf_filenames2
        (y, A, I, J, EI, GJ, avg_centroid,
         plane_bdf_filenames1, plane_bdf_filenames2) = moi_data
        # assert np.allclose(avg_centroid, avg_centroid0)
        # print(f'y = {y.tolist()}')
        # print(f'A = {A.tolist()}')
        # print(f'I = {I.tolist()}')
        # print(f'J = {J.tolist()}')
        # print(f'EI = {EI.tolist()}')
        # print(f'GJ = {GJ.tolist()}')

        y_expected = [1.0, 141.0, 281.0, 421.0, 561.0, 701.0, 841.0, 981.0, 1121.0, 1261.0]
        A_expected = [27325.7035460869, 16791.75582131023, 3439.6593514378583, 1885.2005923508825, 1112.9466242866476, 902.0982620165437, 692.8361798071458, 528.0741288550782, 397.36522577804004, 286.59151276901025]
        I_expected = [
            [3403933746.398361, 8.911280281243883e-27, 87431348.83639666, -1.1524729838475616e-09, -8.115831147516347e-11, -113043322.10390155],
            [1333461947.3875391, 2.2380375734653757e-27, 41162522.72557803, -3.816419553672272e-11, 2.2490616251215358e-11, -16958977.338277496],
            [213026144.23155284, 2.421936632273299e-28, 4799720.5299459305, -2.789842444913242e-12, -3.8705946989251635e-14, -5341199.828407566],
            [41040943.73186442, 1.963256730304734e-28, 545179.3969024371, -1.2228253392113826e-11, -4.381822725941327e-13, 449897.1318028397],
            [10103921.858063448, 2.8359530263005473e-28, 145368.04628856102, 5.970341815480735e-12, 8.713118441540724e-13, 523741.2250703718],
            [5374113.109890341, 1.135309608747665e-28, 87417.81394571027, 1.069894873862433e-12, -1.6083252621051996e-13, 380076.317217249],
            [2845327.0112954406, 6.312458138991356e-29, 53886.33108173666, 2.2219427152776314e-12, 3.100104645122421e-13, 238803.7329846358],
            [1483124.33122053, 8.203848533496161e-29, 30610.86321437863, 1.167927493310214e-12, -8.845064798186607e-14, 142410.19408779056],
            [684110.0163668194, 8.231340803540231e-29, 16725.831058479114, 5.647390666346783e-12, 7.090225668624896e-13, 78866.11262783181],
            [250269.1763145865, 3.010055271280282e-29, 8678.51145346907, 1.2739523923193558e-12, 1.1959210151601072e-13, 38470.74439627658]]
        J_expected = [3491365095.234758, 1374624470.1131172, 217825864.76149878, 41586123.12876686, 10249289.90435201, 5461530.9238360515, 2899213.3423771774, 1513735.1944349087, 700835.8474252985, 258947.6877680556]
        EI_expected = [
            [3.295171322240443e+16, 1.0055942705118451e-19, 834725104132483.5, -0.016966902163427177, -0.0008386804640027322, -1305402543493497.0],
            [1.1424869043875884e+16, 2.436246336255071e-20, 346881503958272.25, -0.0014203565392383113, 0.00012928505137646288, -158076619489430.62],
            [1376220583763065.0, 1.4987514310082692e-21, 30847567187819.86, -6.739996447080615e-05, 2.686140388554683e-07, -37232864334729.03],
            [232849135638355.66, 1.163610376266595e-21, 3370933367866.2837, -6.256159445564982e-05, -2.130227964903126e-06, 2725803693395.542],
            [53434265948650.6, 1.5070800170642208e-21, 756376205291.4624, 2.8350427118543883e-05, 3.4861744237313086e-06, 2421036127784.6895],
            [28448320392572.47, 5.8563062886418695e-22, 452376143171.62006, 1.233289212310551e-05, -4.952364981663227e-07, 1838842284503.897],
            [11922707419038.45, 3.128149669510187e-22, 255514290711.5685, 1.4411389995587915e-05, 1.8297740043860598e-06, 1043116272629.821],
            [5081816000830.39, 3.04143608378661e-22, 110098452419.56664, 4.1518689145376174e-06, -2.915248889290899e-07, 483517837043.9404],
            [2127609200637.4373, 2.5769541338500123e-22, 52962312533.78477, 1.7402752844301667e-05, 2.1839472821258388e-06, 244954003117.0613],
            [778925018736.649, 9.542065160328333e-23, 27450908134.014584, 3.8466670979195605e-06, 3.395408960754233e-07, 119686895556.34503]]
        GJ_expected = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        assert np.allclose(y, y_expected)
        assert np.allclose(A, A_expected)
        assert np.allclose(I, I_expected)
        assert np.allclose(J, J_expected)
        assert np.allclose(EI, EI_expected)
        assert np.allclose(GJ, GJ_expected)

        for plane_bdf_filename in plane_bdf_filenames1:
            os.remove(plane_bdf_filename)
        for plane_bdf_filename in plane_bdf_filenames2:
            os.remove(plane_bdf_filename)

        if IS_MATPLOTLIB:
            plot_inertia(y, A, I, J, EI, GJ, avg_centroid, show=show)
            asdf
            os.remove(dirname / 'normalized_inertia_vs_span.png')
            os.remove(dirname / 'area_vs_span.png')
            os.remove(dirname / 'amoi_vs_span.png')
            os.remove(dirname / 'e_amoi_vs_span.png')
            os.remove(dirname / 'cg_vs_span.png')

        # bdf_merge(plane_bdf_filenames, bdf_filename_out='merge.bdf', renumber=True,
        #           encoding=None, size=8, is_double=False, cards_to_skip=None,
        #           log=None, skip_case_control_deck=False)
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


def get_coords_bwb(ystations: np.ndarray,
                   cid: int=1,
                   base_coord: CORD2R=None,
                   origin: np.ndarray=None,
                   zaxis: np.ndarray=None,
                   xzplane: np.ndarray=None) -> list[CORD2R]:  # pragma: no cover
    """gets coords from y=0 to y=100*ncuts"""
    if base_coord:
        raise NotImplementedError('base_coord is not yet implemented')
    else:
        origin = np.asarray(origin, dtype='float64')
        zaxis = np.asarray(zaxis, dtype='float64')
        xzplane = np.asarray(xzplane, dtype='float64')

    coords = []
    nstation = len(ystations)
    dxyz = np.zeros((nstation, 3), dtype='float64')
    dxyz[:, 1] = ystations

    for dxyzi in dxyz:
        coord = CORD2R(1, rid=0,
                       origin=origin+dxyzi,
                       zaxis=zaxis+dxyzi,
                       xzplane=xzplane+dxyzi)
        coords.append(coord)
    return coords

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
