"""defines cutting plane tests"""
import os
# import copy
import time
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

from pyNastran.bdf.mesh_utils.cut.moi_plotter import (
    cut_and_plot_moi, plot_inertia, _get_shell_inertia, load_moi_data)
from pyNastran.bdf.mesh_utils.cut.cut_model_by_plane import (
    _setup_faces)

PKG_PATH = pyNastran.__path__[0]
TEST_PATH = Path(__file__).parent
MODEL_PATH = Path(os.path.join(PKG_PATH, '..', 'models'))


class TestStiffnessPlot(unittest.TestCase):
    def test_shell_inertia(self):
        log = SimpleLogger(level='warning', encoding='utf-8')
        model = BDF(debug=False, log=log, mode='msc')

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

    def test_cut_quad_mat8(self):
        """cut_and_plot_moi"""
        dirname = TEST_PATH
        log = SimpleLogger(level='warning', encoding='utf-8')
        # log = SimpleLogger(level='debug', encoding='utf-8')
        dy = 0.5
        t = 0.1
        e11 = 100.
        e22 = 300.
        nu12 = 0.3
        g12 = 400.
        ystations = [0.]
        normal_plane = np.array([0., 1., 0.])

        cut_length = 3.0
        # skipped this...1/12 * cut_length * t**3
        i_expected = 0.09375  # A*d^2 of 2 triangles
        y_expected = [0.0]
        A_expected = [t * cut_length]
        I_expected = [[i_expected, 0.0, 0.0, 0.0, 0.0, 0.0]]
        J_expected = [i_expected]
        ExI_expected = [[e22 * i_expected, 0.0, 0.0, 0.0, 0.0, 0.0]]  # 2812500.0
        EyI_expected = [[e11 * i_expected, 0.0, 0.0, 0.0, 0.0, 0.0]]  # 2812500.0
        GJ_expected = [37.5]
        centroid_expected = [[1.5, 0.0, 0.0]]
        E_expected = [e11, e22, g12]
        Ex_expected = [e22]
        Ey_expected = [e11]
        G_expected = [g12]
        pid = 11
        mid = 12
        for type in {'PSHELL', 'PCOMP', 'PCOMPG'}:
            model, coord = _build_quad(log, dy, zoffset=0.0)
            model.add_mat8(mid, e11, e22, nu12, g12=g12, g1z=1e8, g2z=1e8)
            coords = [coord]
            if type == 'PSHELL':
                prop = model.add_pshell(pid, mid, t=t)
            elif type == 'PCOMP':
                prop = model.add_pcomp(pid, mid, [t])
            elif type == 'PCOMPG':
                prop = model.add_pcompg(pid, [8], mid, [t])
            else:  # pragma: no cover
                raise RuntimeError(type)
            model.cross_reference()
            A, B, D = prop.get_individual_ABD_matrices(theta_offset=0.)
            imat_rotation_angle = 0.0
            Ex, Ey, Gxy, nu_xy = prop.get_Ainv_equivalent_pshell(imat_rotation_angle, 0.1)
            assert np.allclose(e11, Ex)
            assert np.allclose(e22, Ey)
            assert np.allclose(g12, Gxy)
            moi_data = cut_and_plot_moi(
                model, normal_plane, log,
                ystations, coords,
                dirname=dirname,
                plot=False, show=False, face_data=None,
                stop_on_failure=True,
                cut_data_span_filename='',
                beam_model_bdf_filename='',
                thetas_csv_filename='y_thetas.csv',
                # normalized_inertia_png_filename='y_normalized_inertia_vs_span.png',
                # area_span_png_filename='y_area_vs_span.png',
                # amoi_span_png_filename='y_amoi_vs_span.png',
                # e_amoi_span_png_filename='y_e_amoi_vs_span.png',
                # cg_span_png_filename='y_cg_vs_span.png',
                debug_vectorize=True,
            )
            (y, A, I, J,
             ExI, EyI, GJ, avg_centroid,
             plane_bdf_filenames1, plane_bdf_filenames2, ifig) = moi_data
            Ex = ExI[:, 0] / I[:, 0]
            Ey = EyI[:, 0] / I[:, 0]
            G = GJ / J
            # print(f'y = {y.tolist()}')
            # print(f'A = {A.tolist()}')
            # print(f'I = {I.tolist()}')
            # print(f'J = {J.tolist()}')
            # print(f'EI = {EI.tolist()}')
            # print(f'GJ = {GJ.tolist()}')
            # print(f'avg_centroid = {avg_centroid.tolist()}')
            assert np.allclose(G, G_expected)
            assert np.allclose(Ex, Ex_expected), (type, Ex, Ex_expected)
            assert np.allclose(Ey, Ey_expected)

            assert np.allclose(y, y_expected)
            assert np.allclose(A, A_expected)
            assert np.allclose(I, I_expected), (I, I_expected)
            assert np.allclose(J, J_expected)
            assert np.allclose(ExI, ExI_expected), ExI.tolist()
            assert np.allclose(EyI, EyI_expected)
            assert np.allclose(GJ, GJ_expected), GJ.tolist()
            assert np.allclose(avg_centroid, centroid_expected), avg_centroid.tolist()
            del model.properties[pid]

    def test_cut_quad_mat1(self):
        """cut_and_plot_moi"""
        dirname = TEST_PATH
        log = SimpleLogger(level='warning', encoding='utf-8')
        # log = SimpleLogger(level='debug', encoding='utf-8')
        dy = 0.5
        t = 0.1
        E = 3.0e7

        pid = 11
        mid = 12
        ystations = [0.]
        normal_plane = np.array([0., 1., 0.])

        cut_length = 3.0
        # skipped this...1/12 * cut_length * t**3
        i_expected = 0.09375  # A*d^2 of 2 triangles
        y_expected = [0.0]
        A_expected = [t*cut_length]
        I_expected = [[i_expected, 0.0, 0.0, 0.0, 0.0, 0.0]]
        J_expected = [i_expected]
        ExI_expected = [[E*i_expected, 0.0, 0.0, 0.0, 0.0, 0.0]]  # 2812500.0
        EyI_expected = [[E*i_expected, 0.0, 0.0, 0.0, 0.0, 0.0]]  # 2812500.0
        GJ_expected = [1081730.7692307692]
        centroid_expected = [[1.5, 0.0, 0.0]]

        cut_data_span_filename = dirname / 'test_cut_quad_shell_mat1.csv'

        for type in {'PSHELL', 'PCOMP', 'PCOMPG'}:
            model, coord = _build_quad(log, dy, zoffset=0.0)
            model.add_mat1(mid, E=E, G=None, nu=0.3)
            coords = [coord]
            if type == 'PSHELL':
                model.add_pshell(pid, mid, t=t)
            elif type == 'PCOMP':
                model.add_pcomp(pid, mid, [t])
            elif type == 'PCOMPG':
                model.add_pcompg(pid, [8], mid, [t])
            else:  # pragma: no cover
                raise RuntimeError(type)
        model.cross_reference()
        moi_data = cut_and_plot_moi(
            model, normal_plane, log,
            ystations, coords,
            dirname=dirname,
            plot=False, show=False, face_data=None,
            stop_on_failure=True,
            cut_data_span_filename=cut_data_span_filename,
            beam_model_bdf_filename='',
            thetas_csv_filename='y_thetas.csv',
            # normalized_inertia_png_filename='y_normalized_inertia_vs_span.png',
            # area_span_png_filename='y_area_vs_span.png',
            # amoi_span_png_filename='y_amoi_vs_span.png',
            # e_amoi_span_png_filename='y_e_amoi_vs_span.png',
            # cg_span_png_filename='y_cg_vs_span.png',
            debug_vectorize=True,
        )
        (y, A, I, J,
         ExI, EyI, GJ, avg_centroid,
         plane_bdf_filenames1, plane_bdf_filenames2, ifig) = moi_data

        y1, A1, I1, J1, ExI1, EyI1, GJ1, avg_centroid1 = load_moi_data(cut_data_span_filename)
        assert np.allclose(y, y1)
        assert np.allclose(A, A1)
        assert np.allclose(I, I1)
        assert np.allclose(J, J1)
        assert np.allclose(ExI, ExI1)
        assert np.allclose(EyI, EyI1)
        assert np.allclose(GJ, GJ1)
        assert np.allclose(avg_centroid, avg_centroid1)

        # print(f'y = {y.tolist()}')
        # print(f'A = {A.tolist()}')
        # print(f'I = {I.tolist()}')
        # print(f'J = {J.tolist()}')
        # print(f'EI = {EI.tolist()}')
        # print(f'GJ = {GJ.tolist()}')
        # print(f'avg_centroid = {avg_centroid.tolist()}')

        assert np.allclose(y, y_expected)
        assert np.allclose(A, A_expected)
        assert np.allclose(I, I_expected), (I, I_expected)
        assert np.allclose(J, J_expected)
        assert np.allclose(ExI, ExI_expected)
        assert np.allclose(EyI, EyI_expected)
        assert np.allclose(GJ, GJ_expected), GJ.tolist()
        assert np.allclose(avg_centroid, centroid_expected), avg_centroid.tolist()
        del model.properties[pid]

    def test_cut_quad_shell_mat1_zoffset(self):
        """cut_and_plot_moi"""
        dirname = TEST_PATH
        log = SimpleLogger(level='warning', encoding='utf-8')
        # log = SimpleLogger(level='debug', encoding='utf-8')
        dy = 0.5
        model, coord = _build_quad(log, dy, zoffset=10.0)
        t = 0.1
        E = 3.0e7
        model.add_pshell(11, 12, t=t)
        model.add_mat1(12, E=E, G=None, nu=0.3)
        model.cross_reference()
        coords = [coord]
        ystations = [0.]
        normal_plane = np.array([0., 1., 0.])

        moi_data = cut_and_plot_moi(
            model, normal_plane, log,
            ystations, coords,
            dirname=dirname,
            plot=False, show=False, face_data=None,
            stop_on_failure=True,
            cut_data_span_filename='',
            beam_model_bdf_filename='',
            thetas_csv_filename='y_thetas.csv',
            # normalized_inertia_png_filename='y_normalized_inertia_vs_span.png',
            # area_span_png_filename='y_area_vs_span.png',
            # amoi_span_png_filename='y_amoi_vs_span.png',
            # e_amoi_span_png_filename='y_e_amoi_vs_span.png',
            # cg_span_png_filename='y_cg_vs_span.png',
            debug_vectorize=True,
        )
        (y, A, I, J,
         ExI, EyI, GJ, avg_centroid,
         plane_bdf_filenames1, plane_bdf_filenames2, ifig) = moi_data
        # print(f'y = {y.tolist()}')
        # print(f'A = {A.tolist()}')
        # print(f'I = {I.tolist()}')
        # print(f'J = {J.tolist()}')
        # print(f'EI = {EI.tolist()}')
        # print(f'GJ = {GJ.tolist()}')
        # print(f'avg_centroid = {avg_centroid.tolist()}')
        cut_length = 3.0
        # skipped this...1/12 * cut_length * t**3
        i_expected = 0.09375  # A*d^2 of 2 triangles
        y_expected = [0.0]
        A_expected = [t*cut_length]
        I_expected = [[i_expected, 0.0, 0.0, 0.0, 0.0, 0.0]]
        J_expected = [i_expected]
        ExI_expected = [[E*i_expected, 0.0, 0.0, 0.0, 0.0, 0.0]]  # 2812500.0
        EyI_expected = [[E*i_expected, 0.0, 0.0, 0.0, 0.0, 0.0]]  # 2812500.0
        GJ_expected = [1081730.7692307692]
        centroid_expected = [[1.5, 0.0, 10.0]]

        assert np.allclose(y, y_expected)
        assert np.allclose(A, A_expected)
        assert np.allclose(I, I_expected), (I, I_expected)
        assert np.allclose(J, J_expected)
        assert np.allclose(ExI, ExI_expected)
        assert np.allclose(EyI, EyI_expected)
        assert np.allclose(GJ, GJ_expected), GJ.tolist()
        with self.assertRaises(AssertionError):
            assert np.allclose(avg_centroid, centroid_expected), avg_centroid.tolist()

    def test_cut_quad_shell_zoffset(self):
        """cut_and_plot_moi"""
        dirname = TEST_PATH
        log = SimpleLogger(level='warning', encoding='utf-8')
        # log = SimpleLogger(level='debug', encoding='utf-8')
        dy = 0.5
        model, coord = _build_quad(log, dy, zoffset=10.0)
        t = 0.1
        E = 3.0e7
        model.add_pshell(11, 12, t=t)
        model.add_mat1(12, E=E, G=None, nu=0.3)
        model.cross_reference()
        coords = [coord]
        ystations = [0.]
        normal_plane = np.array([0., 1., 0.])

        moi_data = cut_and_plot_moi(
            model, normal_plane, log,
            ystations, coords,
            dirname=dirname,
            plot=False, show=False, face_data=None,
            stop_on_failure=True,
            cut_data_span_filename='',
            beam_model_bdf_filename='',
            thetas_csv_filename='y_thetas.csv',
            # normalized_inertia_png_filename='y_normalized_inertia_vs_span.png',
            # area_span_png_filename='y_area_vs_span.png',
            # amoi_span_png_filename='y_amoi_vs_span.png',
            # e_amoi_span_png_filename='y_e_amoi_vs_span.png',
            # cg_span_png_filename='y_cg_vs_span.png',
            debug_vectorize=True,
        )
        (y, A, I, J,
         ExI, EyI, GJ, avg_centroid,
         plane_bdf_filenames1, plane_bdf_filenames2, ifig) = moi_data
        # print(f'y = {y.tolist()}')
        # print(f'A = {A.tolist()}')
        # print(f'I = {I.tolist()}')
        # print(f'J = {J.tolist()}')
        # print(f'EI = {EI.tolist()}')
        # print(f'GJ = {GJ.tolist()}')
        # print(f'avg_centroid = {avg_centroid.tolist()}')
        cut_length = 3.0
        # skipped this...1/12 * cut_length * t**3
        i_expected = 0.09375  # A*d^2 of 2 triangles
        y_expected = [0.0]
        A_expected = [t*cut_length]
        I_expected = [[i_expected, 0.0, 0.0, 0.0, 0.0, 0.0]]
        J_expected = [i_expected]
        ExI_expected = [[E*i_expected, 0.0, 0.0, 0.0, 0.0, 0.0]]  # 2812500.0
        EyI_expected = [[E*i_expected, 0.0, 0.0, 0.0, 0.0, 0.0]]  # 2812500.0
        GJ_expected = [1081730.7692307692]
        centroid_expected = [[1.5, 0.0, 0.0]]

        assert np.allclose(y, y_expected)
        assert np.allclose(A, A_expected)
        assert np.allclose(I, I_expected), (I, I_expected)
        assert np.allclose(J, J_expected)
        assert np.allclose(ExI, ExI_expected)
        assert np.allclose(EyI, EyI_expected)
        assert np.allclose(GJ, GJ_expected), GJ.tolist()
        assert np.allclose(avg_centroid, centroid_expected), avg_centroid.tolist()

    def test_cut_tet(self):
        """cut_and_plot_moi"""
        dirname = TEST_PATH
        log = SimpleLogger(level='warning', encoding='utf-8')
        # log = SimpleLogger(level='debug', encoding='utf-8')
        dy = 0.5
        model, coord = _build_tet(log, dy)
        t = 0.1
        E = 3.0e7
        model.add_psolid(11, 12)
        model.add_mat1(12, E=E, G=None, nu=0.3)
        model.cross_reference()
        coords = [coord]
        ystations = [0.]
        normal_plane = np.array([0., 1., 0.])

        with self.assertRaises(NotImplementedError):
            moi_data = cut_and_plot_moi(
                model, normal_plane, log,
                ystations, coords,
                dirname=dirname,
                plot=False, show=False, face_data=None,
                include_solids=True,
                stop_on_failure=True,
                cut_data_span_filename='',
                beam_model_bdf_filename='',
                thetas_csv_filename='y_thetas.csv',
                # normalized_inertia_png_filename='y_normalized_inertia_vs_span.png',
                # area_span_png_filename='y_area_vs_span.png',
                # amoi_span_png_filename='y_amoi_vs_span.png',
                # e_amoi_span_png_filename='y_e_amoi_vs_span.png',
                # cg_span_png_filename='y_cg_vs_span.png',
                debug_vectorize=True,
            )
        # (y, A, I, J,
        #  ExI, EyI, GJ, avg_centroid,
        #  plane_bdf_filenames1, plane_bdf_filenames2, ifig) = moi_data
        # print(f'y = {y.tolist()}')
        # print(f'A = {A.tolist()}')
        # print(f'I = {I.tolist()}')
        # print(f'J = {J.tolist()}')
        # print(f'ExI = {ExI.tolist()}')
        # print(f'EyI = {EyI.tolist()}')
        # print(f'GJ = {GJ.tolist()}')
        # # print(f'avg_centroid = {avg_centroid.tolist()}')
        # cut_area = 3.0
        # # skipped this...1/12 * cut_length * t**3
        # i_expected = 0.09375  # A*d^2 of 2 triangles
        # y_expected = [0.0]
        # A_expected = [cut_area]
        # I_expected = [[i_expected, 0.0, 0.0, 0.0, 0.0, 0.0]]
        # J_expected = [i_expected]
        # ExI_expected = [[E*i_expected, 0.0, 0.0, 0.0, 0.0, 0.0]]  # 2812500.0
        # EyI_expected = [[E*i_expected, 0.0, 0.0, 0.0, 0.0, 0.0]]  # 2812500.0
        # GJ_expected = [1081730.7692307692]
        # centroid_expected = [[1.5, 0.0, 0.0]]
        #
        # assert np.allclose(y, y_expected)
        # assert np.allclose(A, A_expected)
        # assert np.allclose(I, I_expected), (I, I_expected)
        # assert np.allclose(J, J_expected)
        # assert np.allclose(ExI, ExI_expected)
        # assert np.allclose(EyI, EyI_expected)
        # assert np.allclose(GJ, GJ_expected), GJ.tolist()
        # assert np.allclose(avg_centroid, centroid_expected), avg_centroid.tolist()

    def test_cut_bwb(self):
        """cut_and_plot_moi"""
        # show = True
        t0 = time.time()
        show = False
        run_y_cuts = True
        run_x_cuts = True

        log = SimpleLogger(level='warning', encoding='utf-8')
        # log = SimpleLogger(level='debug', encoding='utf-8')
        dirname = MODEL_PATH / 'bwb'
        bdf_filename = dirname / 'bwb_saero.bdf'  # ymax~=1262.0
        model = read_bdf(bdf_filename, log=log)
        # model.log.level = 'debug'

        ymax = 1401.
        ncut_y = 10
        dy = (ymax-1.) / ncut_y
        i = np.arange(ncut_y, dtype='int32')
        # xstations = 50 * i + 1.
        ystations = dy * i + 1.

        # y is outboard
        origin = np.array([0., 0., 0.])
        zaxis = np.array([0., 0., 1.])   # z
        xzplane = np.array([1., 0., 0.]) # x

        cid0, coords = get_coords_bwb(
            model, ystations, cid0=-1,
            origin=origin,
            zaxis=zaxis,
            xzplane=xzplane)
        add_coords = True
        bdf_filename_out = dirname / 'y_bwb_saero.bdf'
        fadd_coords(model, coords, bdf_filename_out,
                    add_coords=add_coords)

        normal_plane = coords[0].j
        # normal_plane = np.array([0., 1., 0.])
        # assert np.allclose(normal_plane, normal_plane2)

        log, *face_data = _setup_faces(model)
        # nids, xyz_cid0, elements = face_data
        # y0, A0, I0, J0, EI0, J0, avg_centroid0, plane_bdf_filenames10, plane_bdf_filenames20 = cut_and_plot_moi(
        #     model, normal_plane, log,
        #     dys, coords,
        #     dirname=dirname,
        #     plot=False, show=False, face_data=face_data)

        if run_y_cuts:
            log.info('working on y-cuts')
            moi_data = cut_and_plot_moi(
                bdf_filename, normal_plane, log,
                ystations, coords,
                dirname=dirname,
                plot=True, show=False, face_data=face_data,
                cut_data_span_filename='y_cut_data_vs_span.csv',
                beam_model_bdf_filename='y_equivalent_beam_model.bdf',
                thetas_csv_filename='y_thetas.csv',
                normalized_inertia_png_filename='y_normalized_inertia_vs_span.png',
                area_span_png_filename='y_area_vs_span.png',
                amoi_span_png_filename='y_amoi_vs_span.png',
                e_amoi_span_png_filename='y_e_amoi_vs_span.png',
                cg_span_png_filename='y_cg_vs_span.png',
                debug_vectorize=True,
            )
            (y, A, I, J,
             ExI, EyI, GJ, avg_centroid,
             plane_bdf_filenames1, plane_bdf_filenames2, ifig) = moi_data
            # assert np.allclose(avg_centroid, avg_centroid0)
            # print(f'y = {y.tolist()}')
            # print(f'A = {A.tolist()}')
            # print(f'I = {I.tolist()}')
            # print(f'J = {J.tolist()}')
            # print(f'EI = {EI.tolist()}')
            # print(f'GJ = {GJ.tolist()}')

            y_expected = [1.0, 141.0, 281.0, 421.0, 561.0, 701.0, 841.0, 981.0, 1121.0, 1261.0]
            A_expected = [27325.7035460869, 16791.75582131023, 3439.6593514378583, 1885.2005923508825, 1112.9466242866476,
                          902.0982620165437, 692.8361798071458, 528.0741288550782, 397.36522577804004, 286.59151276901025]
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
            J_expected = [3491365095.234758, 1374624470.1131172, 217825864.76149878, 41586123.12876686, 10249289.90435201, 5461530.9238360515,
                          2899213.3423771774, 1513735.1944349087, 700835.8474252985, 258947.6877680556]
            ExI_expected = [
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
            EyI_expected = [
                [3.30291794969561e+16, 9.914430634459927e-20, 829639166631982.6, -0.016724795295494543, -0.0008199629361946648, -1311475594974912.8],
                [1.1367775459177206e+16, 2.4110535055375546e-20, 345163741056031.75, -0.0013782381954762866, 0.00012794305713784447, -157741311102661.12],
                [1363575430159224.5, 1.487036422484434e-21, 30598238977826.594, -7.104785875813303e-05, -5.969107229406736e-08, -38430882185175.41],
                [234021754167591.72, 1.1694514144944916e-21, 3426073080525.612, -6.280806449176725e-05, -2.2836120177074253e-06, 2854885803354.4224],
                [53700287868470.336, 1.4738734024992309e-21, 740959869292.9684, 2.8334698365116492e-05, 3.58910621827633e-06, 2411945681754.337],
                [28577869538088.812, 5.814147705733206e-22, 443705495395.6413, 1.3012154621287177e-05, -4.816673998844504e-07, 1834538624017.6018],
                [11881165360517.701, 3.0568127168517245e-22, 251438302297.99298, 1.4080660253727768e-05, 1.7827618937643696e-06, 1039150935120.2333],
                [5008238252044.59, 3.0875799633053617e-22, 110181565188.11455, 4.540787324493114e-06, -2.5512325871018163e-07, 476060331894.8157],
                [2128603897250.3176, 2.5772461523789405e-22, 52995377672.16867, 1.7416397212813894e-05, 2.186036081758313e-06, 245172127763.47583],
                [779303809782.374, 9.540748307180537e-23, 27469421334.232525, 3.843636626470701e-06, 3.38601732516192e-07, 119778660481.12492]]
            GJ_expected = [1.703986879714322e+16, 5994675543453522.0, 715043894381704.6, 118078327598852.66, 28896697595393.2, 15363440774671.04,
                           5753383451231.989, 2920241713409.1665, 1214164704037.6692, 447803875533.6692]
            assert np.allclose(y, y_expected)
            assert np.allclose(A, A_expected)
            assert np.allclose(I, I_expected)
            assert np.allclose(J, J_expected)
            assert np.allclose(ExI, ExI_expected)
            assert np.allclose(EyI, EyI_expected), EyI.tolist()
            assert np.allclose(GJ, GJ_expected), GJ.tolist()
            for plane_bdf_filename in plane_bdf_filenames1:
                os.remove(plane_bdf_filename)
            for plane_bdf_filename in plane_bdf_filenames2:
                os.remove(plane_bdf_filename)

            if IS_MATPLOTLIB:
                plot_inertia(log, y, A, I, J, ExI, EyI, GJ, avg_centroid, show=show)

        if run_x_cuts:
            # xmax = 1001.
            ncut_x = 21
            # dx = (xmax-1.) / ncut
            i = np.arange(ncut_x, dtype='int32')
            # xstations = 10 * i + 1
            xstations = np.linspace(0.14312, 1614.17, num=ncut_x)[1:-1]

            # y is outboard
            origin = np.array([0., 0., 0.])
            zaxis = np.array([0., 0., 1.])   # z
            xzplane = np.array([0., 1., 0.]) # y

            cid0, xcoords = get_coords_bwb(
                model, xstations, axis=0, cid0=-1,
                origin=origin,
                zaxis=zaxis,
                xzplane=xzplane)

            log.info('working on x-cuts')
            # xmin=0.14312  xmax=1614.17 dx=1614.027
            # ymin=-0.01102 ymax=1262.0  dy=1262.011
            # zmin=-105.05  zmax=282.209 dz=387.25903
            normal_plane = coords[0].i
            log.debug(f'normal_plane = {normal_plane}')
            moi_data = cut_and_plot_moi(
                bdf_filename, normal_plane, log,
                xstations, xcoords,
                dirname=dirname, ifig=10,
                plot=True, show=False, face_data=face_data,
                cut_data_span_filename='x_cut_data_vs_span.csv',
                beam_model_bdf_filename='x_equivalent_beam_model.bdf',
                thetas_csv_filename='x_thetas.csv',
                normalized_inertia_png_filename='x_normalized_inertia_vs_span.png',
                area_span_png_filename='x_area_vs_span.png',
                amoi_span_png_filename='x_amoi_vs_span.png',
                e_amoi_span_png_filename='x_e_amoi_vs_span.png',
                cg_span_png_filename='x_cg_vs_span.png',
                debug_vectorize=True,
            )
            (x, A, I, J,
             ExI, EyI, GJ, avg_centroid,
             plane_bdf_filenames1, plane_bdf_filenames2, ifig) = moi_data
            # log.warning(f'x = {x.tolist()}')
            # log.warning(f'A = {A.tolist()}')
            x_expected = [80.844464, 161.54580800000002, 242.24715200000003, 322.94849600000003, 403.64984000000004, 484.35118400000005,
                          565.0525279999999, 645.753872, 726.4552160000001, 807.15656, 887.857904, 968.559248, 1049.260592, 1129.961936,
                          1210.66328, 1291.364624, 1372.065968, 1452.7673120000002, 1533.468656]
            ax_expected = [392.32365100458486, 1093.07584073101, 1648.1132186236794, 2048.09990888604, 2418.1562457012187, 2774.3435792781092,
                           4338.529928442088, 4973.310993528411, 5695.986120156011, 6028.017865844101, 6406.784071105117, 5143.303096908391,
                           5604.80749223895, 6056.25268430806, 4670.299171635897, 4201.107609278721, 3235.886100101357, 2083.159010726759, 1441.6045438028273]
            I_expected = [
                [221897.53083113694, 3.7983639915375733e-29, 490585.4080952155, -1.3634463940838429e-12, -2.2362747388142087e-13, -19311.70957664522],
                [1135660.8112785984, 1.373539641117107e-28, 1952549.4967115684, -5.438462476795749e-12, -4.811809236665373e-12, 564774.7917087536],
                [2250752.9886096143, 4.26935073673555e-28, 3139518.92069423, -4.609031580526869e-13, -1.5938656038357523e-11, 793751.846454918],
                [3769691.8313080794, 5.829553049402392e-28, 4404891.178114804, -1.499956998435433e-11, -3.5346722713561794e-12, 1066016.2123976718],
                [5823377.837602431, 4.167896798995401e-28, 5730123.807352846, 3.3068081084315233e-12, 6.312571154282228e-12, 1363768.3947602392],
                [8571967.00998755, 3.6990385173288806e-28, 7017210.056155724, -1.0354438534431915e-12, 1.723521682412084e-11, 1715442.5637448356],
                [17683660.38590234, 4.715842100428592e-28, 14848167.185298493, -1.1210711198113354e-11, -1.2411601462789255e-11, 7108422.4829052845],
                [24932189.990431517, 7.39345153195075e-28, 17015656.141243584, -2.4796820516447402e-11, -2.4356887460267413e-11, 8976769.662816338],
                [38698154.80266441, 1.0840939814800784e-27, 19796331.194030203, 8.304407185107903e-12, 4.849429924380216e-12, 13308830.968513666],
                [51304288.515147164, 8.318971749911368e-28, 21493956.676934727, 1.623710787739812e-11, 1.2600753099333035e-12, 16706547.975627016],
                [75349425.20702761, 8.187071433949908e-28, 23310243.883227274, -1.877845800040756e-11, -1.2856240689643125e-12, 22806540.732509054],
                [94527894.60177897, 1.1747030447106125e-27, 13431125.329021016, 8.331331649168944e-11, 3.3419442992674405e-11, 17074852.872649502],
                [175704618.55822277, 1.71744539566443e-27, 14176002.586274857, 2.7877534629651125e-10, 4.335586873204759e-11, 27765440.09937934],
                [305988717.9499091, 1.6728844724723693e-27, 14652825.334953113, -1.6677078996340593e-10, -4.563688685161472e-11, 41635691.16114099],
                [362118508.1307739, 6.039432586887828e-28, 11645206.668794185, -8.365132003151304e-11, -5.707375026785586e-12, 36289871.977354966],
                [472358953.592614, 4.6507529070138925e-28, 8795869.026929501, 1.9249370289199633e-12, -7.925298419066906e-12, 39000766.40437172],
                [529230564.55807215, 4.282371898532258e-28, 6413078.724857627, 4.20476652404626e-11, -3.491202371845323e-13, 34152362.65613742],
                [32243503.639903653, 2.7056944874312897e-28, 7683345.754385823, -6.436079031265966e-12, -1.6989270200404068e-12, 10858173.343513252],
                [23430036.049582027, 1.4352984163003176e-28, 12452844.40630308, -2.1036973800523316e-12, -1.8083012004535386e-14, 14853002.6361171]]
            J_expected = [712482.9389263524, 3088210.307990167, 5390271.909303844, 8174583.009422883, 11553501.644955277, 15589177.066143274, 32531827.571200833, 41947846.1316751, 58494485.99669461, 72798245.1920819, 98659669.09025489, 107959019.93079999, 189880621.14449763, 320641543.2848622, 373763714.79956806, 481154822.6195435, 535643643.2829298, 39926849.39428948, 35882880.455885105]
            ExI_expected = [
                [1637184286334.9233, 2.8021158468839396e-22, 3618876781150.281, -1.0052182112457289e-05, -1.6494873953252972e-06, -142809814651.8627],
                [10095614240289.502, 1.2089151762018814e-21, 16496285015070.463, -4.8028779838279166e-05, -4.87899394117178e-05, 5698885665116.398],
                [21038398805027.05, 4.2605944309416014e-21, 25320302308082.0, -1.6636888485889658e-06, -0.00014606196456040527, 7332025815600.03],
                [36411031132090.94, 6.665201807409123e-21, 34418535949989.42, -0.00020786705924220097, -2.7303889612632887e-05, 9407594048506.8],
                [57550645552937.71, 4.2241651394029154e-21, 43990464726705.516, 5.998891792727669e-05, 4.7342033772719154e-05, 11659186213235.516],
                [85424002799352.7, 3.49951384710969e-21, 53251192618309.97, -3.537154912603066e-06, 0.00013590754235504182, 14354555880835.223],
                [176622053933508.28, 5.307940297781259e-21, 135751786495349.58, -4.738874066985743e-05, -0.00010242628300806387, 69346253740960.26],
                [249466387285724.7, 8.202013634455497e-21, 155757465096236.12, -0.0002686478510943583, -0.0002168431928989825, 87073833903341.84],
                [367383746166928.56, 1.2503811897092048e-20, 180540861457562.78, 8.621118050627118e-05, 2.728447857832072e-05, 123836287114791.4],
                [461655146472030.1, 9.18615560520371e-21, 195549780487692.34, 0.00012846367513981688, -1.904471194880576e-05, 150508668022991.8],
                [648341516638356.9, 9.011864363061857e-21, 211251206531269.34, -0.00020018763986618358, -1.7537203521728115e-05, 200039129450810.88],
                [782483071726457.9, 1.1717120651666672e-20, 101993829187559.33, 0.0007305614145337446, 0.00027863842544280165, 137298401830296.78],
                [1355225209430539.0, 1.7334507908580536e-20, 107818763599914.86, 0.0026674695413966697, 0.00040174651955560234, 217383757877501.8],
                [2242515811484036.8, 1.7260998496012703e-20, 111755788030440.06, -0.0014687192381861167, -0.0004387903428892533, 319396170865319.44],
                [1944932705274887.0, 3.7519449378182506e-21, 70445422696078.83, -0.0004753615127818212, -3.51835918567862e-05, 198332953987211.25],
                [2333753367080524.0, 2.4407048275683007e-21, 48905243375276.664, 3.9620758039479674e-05, -4.687931897931743e-05, 197555430503770.53],
                [2217886393895651.2, 2.4240054793082022e-21, 32223616895968.605, 0.00017244028794740013, -5.857456697942972e-06, 152931577653961.5],
                [180576346079041.34, 1.5916574070506534e-21, 30966361158264.15, -4.149281975540677e-05, -3.7689020243736565e-06, 48738001502955.47],
                [110920646839548.8, 8.300475585830254e-22, 48766620679325.81, -1.3918000555511405e-05, -3.0903813983376974e-06, 61924176288314.54]]
            EyI_expected = [
                [1637184297585.3806, 2.802115852694404e-22, 3618876827638.434, -1.0052182167397996e-05, -1.6494875823069066e-06, -142809816431.40533],
                [9991302026715.37, 1.1983283412620032e-21, 16403363880205.479, -4.756739789982026e-05, -4.81992375815353e-05, 5608924542301.52],
                [20890605271618.906, 4.255811046271339e-21, 25284120946604.812, -1.9920787201381187e-06, -0.0001459527857664239, 7287269766076.512],
                [36100357532413.55, 6.5927821972141556e-21, 34372927448596.66, -0.00020475070151550618, -2.6677403218795737e-05, 9330999569973.676],
                [56903902821919.0, 4.175129272569716e-21, 43921865995788.82, 5.8346537359325534e-05, 4.717956084182895e-05, 11548387469985.654],
                [84804917203632.2, 3.498201229832784e-21, 53238237617047.72, -2.791573804861473e-06, 0.00013583064410259906, 14306962105111.543],
                [178688595056977.03, 5.380644804163906e-21, 136885663434360.75, -4.5172413495305427e-05, -0.00010230740073449735, 70130806791742.47],
                [252283953176704.94, 8.312623952590479e-21, 157068848483454.7, -0.0002728251621293441, -0.00021831029579883485, 88009093211224.89],
                [372389984622246.75, 1.2722463270796538e-20, 182156799559390.28, 8.501437439866044e-05, 2.6144620705329416e-05, 125298737103408.47],
                [468579684669571.25, 9.349364164357739e-21, 197479491034590.1, 0.00013080614028360634, -2.055079067295667e-05, 152488067111204.53],
                [653559695539177.8, 9.166935802571101e-21, 213195118563393.2, -0.00020462862630587442, -1.832444451599298e-05, 201634971676723.44],
                [785813254542123.6, 1.187349412526408e-20, 102157203311272.17, 0.0007316192146983297, 0.0002796743179911789, 137691465576039.52],
                [1353767746854281.5, 1.7568021144828676e-20, 107896536459465.67, 0.002695866376948589, 0.00040507963646525396, 217208614483178.03],
                [2272047618206551.5, 1.7555759478234918e-20, 112595348240389.12, -0.00149801851443868, -0.00044664835442235346, 324228566309275.25],
                [1971742270713136.0, 3.7745462803104216e-21, 70527356841834.64, -0.0004841126157461029, -3.7046020751910726e-05, 200662496883935.7],
                [2315493416347223.5, 2.4221646518754087e-21, 48177774317139.5, 4.370369752080637e-05, -4.591235823461671e-05, 194462725305966.84],
                [2211671016460806.8, 2.4204112315448614e-21, 31903448463639.715, 0.00017144782851788044, -6.110556166269301e-06, 151788127523731.66],
                [181571584942580.62, 1.597938734326795e-21, 31088654693132.008, -4.1771709693110584e-05, -3.678096539091203e-06, 48978902253935.82],
                [111308913025634.05, 8.327247466975484e-22, 48874758600414.11, -1.4205997157181862e-05, -3.273198197739828e-06, 62118990365954.03]]

            GJ_expected = [2743207688607.3154, 13519468742297.812, 23529818190215.14, 35837638612444.36, 50967502623265.56,
                  69697812451200.37, 155195580780968.66, 200997256485291.75, 270141129771216.94, 326173799506531.4,
                  425198906188080.0, 437298336390693.6, 735240713516634.2, 1127198943532038.0, 1066322381238653.2,
                  1236997531280349.5, 1172560428793695.0, 112820404543552.34, 85435975493538.38]

            # print(f'x = {x.tolist()}')
            # print(f'A = {A.tolist()}')
            # print(f'I = {I.tolist()}')
            # print(f'J = {J.tolist()}')
            # print(f'EI = {EI.tolist()}')
            # print(f'GJ = {GJ.tolist()}')
            assert np.allclose(x, x_expected)
            assert np.allclose(A, ax_expected)
            assert np.allclose(I, I_expected)
            assert np.allclose(J, J_expected)
            assert np.allclose(ExI, ExI_expected)
            assert np.allclose(EyI, EyI_expected), EyI.tolist()
            assert np.allclose(GJ, GJ_expected), GJ.tolist()

            # assert np.allclose(x, x_expected)
            # assert np.allclose(A, ax_expected)
            log.info('cleanup')
            for plane_bdf_filename in plane_bdf_filenames1:
                os.remove(plane_bdf_filename)
            for plane_bdf_filename in plane_bdf_filenames2:
                os.remove(plane_bdf_filename)
        if run_x_cuts:
            _cleanup_moi_files(dirname, 'x_')
        if run_y_cuts:
            _cleanup_moi_files(dirname, 'y_')
        print(f'dt = {time.time() - t0}')


def _cleanup_moi_files(dirname: Path, tag: str) -> None:
    if IS_MATPLOTLIB:
        os.remove(dirname / f'{tag}normalized_inertia_vs_span.png')
        os.remove(dirname / f'{tag}area_vs_span.png')
        os.remove(dirname / f'{tag}amoi_vs_span.png')
        os.remove(dirname / f'{tag}e_amoi_vs_span.png')
        os.remove(dirname / f'{tag}cg_vs_span.png')

    # bdf_merge(plane_bdf_filenames, bdf_filename_out='merge.bdf', renumber=True,
    #           encoding=None, size=8, is_double=False, cards_to_skip=None,
    #           log=None, skip_case_control_deck=False)
    os.remove(dirname / f'{tag}thetas.csv')
    # os.remove(dirname / 'equivalent_beam_model.bdf')
    os.remove(dirname / f'{tag}cut_data_vs_span.csv')
    # os.remove('cut_face.csv')
    # if IS_MATPLOTLIB:
    #     os.remove('area_vs_span.png')
    #     os.remove('amoi_vs_span.png')
    #     os.remove('normalized_inertia_vs_span.png')
    #     os.remove('cg_vs_span.png')
    #     os.remove('e_amoi_vs_span.png')


def get_coords_bwb(model: BDF,
                   ystations: np.ndarray,
                   axis: int=1,
                   cid0: int=-1,
                   base_coord: CORD2R=None,
                   origin: np.ndarray=None,
                   zaxis: np.ndarray=None,
                   xzplane: np.ndarray=None) -> tuple[int, list[CORD2R]]:
    """gets coords from y=0 to y=100*ncuts"""
    if cid0 == -1:
        cid0 = max(model.coords) + 1

    if base_coord:
        raise NotImplementedError('base_coord is not yet implemented')
    else:
        origin = np.asarray(origin, dtype='float64')
        zaxis = np.asarray(zaxis, dtype='float64')
        xzplane = np.asarray(xzplane, dtype='float64')

    coords = []
    nstation = len(ystations)
    dxyz = np.zeros((nstation, 3), dtype='float64')
    dxyz[:, axis] = ystations

    for icid, dxyzi in enumerate(dxyz):
        cid = cid0 + icid
        coord = CORD2R(cid, rid=0,
                       origin=origin+dxyzi,
                       zaxis=zaxis+dxyzi,
                       xzplane=xzplane+dxyzi)
        coords.append(coord)
    return cid0, coords

def _build_quad(log: SimpleLogger,
                dy: float,
                zoffset: float=0.0) -> tuple[BDF, CORD2R]:
    """
    ^ y
    4-----3
    |     |
    |     |
    1-----2--->x
    """
    model = BDF(log=log)
    model.add_grid(101, [0., 0., 0.])
    model.add_grid(102, [3., 0., 0.])
    model.add_grid(103, [3., 3., 0.])
    model.add_grid(104, [0., 3., 0.])
    model.add_cquad4(10, 11, [101, 102, 103, 104], zoffset=zoffset)
    coord = model.add_cord2r(
        cid=1,
        origin=[0., dy, 0.],
        zaxis=[0., dy, 1.],
        xzplane=[1., dy, 0.])
    return model, coord

def _build_tet(log: SimpleLogger, dy: float):
    model = BDF(log=log)
    model.add_grid(101, [0., 0., 0.])
    model.add_grid(102, [3., 0., 0.])
    model.add_grid(103, [3., 3., 0.])
    model.add_grid(104, [0., 0., 3.])
    model.add_ctetra(10, 11, [101, 102, 103, 104])
    coord = model.add_cord2r(
        cid=1,
        origin=[0., dy, 0.],
        zaxis=[0., dy, 1.],
        xzplane=[1., dy, 0.])
    return model, coord

def fadd_coords(model: BDF, coords: list,
                bdf_filename_out: Path,
                add_coords: bool=True) -> None:
    if not add_coords:
        return
    for coord in coords:
        cid = coord.cid
        assert cid not in model.coords
        model.coords[cid] = coord
    model.write_bdf(bdf_filename_out)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
