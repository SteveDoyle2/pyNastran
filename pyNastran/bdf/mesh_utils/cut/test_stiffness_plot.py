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

try:
    import pandas as pd
    IS_PANDAS = True
except ModuleNotFoundError:
    IS_PANDAS = False

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

#: the nine entries the tests below want out of ``cut_and_plot_moi``'s dict
CORE_KEYS = ('stations', 'L', 'A', 'I', 'J', 'ExI', 'EyI', 'GJ', 'avg_centroid')


def unpack_moi(out_dict: dict) -> tuple:
    """
    ``(stations, L, A, I, J, ExI, EyI, GJ, avg_centroid)``, by name.

    The dict has grown entries over time (neutral axis, shear center), so
    unpacking ``out_dict.values()`` positionally breaks every test at once
    the next time something is added.
    """
    missing = [key for key in CORE_KEYS if key not in out_dict]
    assert not missing, f'cut_and_plot_moi stopped returning {missing}'
    return tuple(out_dict[key] for key in CORE_KEYS)


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
        L_expected = [cut_length]
        A_expected = [t * cut_length]
        I_expected = [[i_expected, 0.0, 0.0, 0.0, 0.0, 0.0]]
        J_expected = [i_expected]
        ExI_expected = [[e22 * i_expected, 0.0, 0.0, 0.0, 0.0, 0.0]]  # 2812500.0
        EyI_expected = [[e11 * i_expected, 0.0, 0.0, 0.0, 0.0, 0.0]]  # 2812500.0
        # A single cut quad is one straight wall: an OPEN section, so the
        # torsion constant is the thin-strip value j = s*t^3/3, not the polar
        # moment.  The old baseline was G*(Ix+Iz) = 400*0.09375 = 37.5, which
        # is 93.75x too stiff -- Ix here is the in-plane second moment of a
        # 3-wide strip and describes bending, not twist.
        j_open = cut_length * t ** 3 / 3.  # 0.001
        GJ_expected = [g12 * j_open]  # 0.4
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
            # expected:
            # percent=0.16666666666666666 (a,b)=(1,2) -> (102,103)
            #   avg_local = [ 3.00000000e+00 -5.55111512e-17  0.00000000e+00]
            #   p1_local  = [ 3.  -0.5  0. ]
            #   p2_local  = [3.  2.5 0. ]
            # percent=0.16666666666666666 (a,b)=(0,2) -> (101,103)
            #   avg_local = [ 5.00000000e-01 -5.55111512e-17  0.00000000e+00]
            #   p1_local  = [ 0.  -0.5  0. ]
            #   p2_local  = [3.  2.5 0. ]
            # percent=0.16666666666666666 (a,b)=(0,2) -> (101,103)
            #   avg_local = [ 5.00000000e-01 -5.55111512e-17  0.00000000e+00]
            #   p1_local  = [ 0.  -0.5  0. ]
            #   p2_local  = [3.  2.5 0. ]
            # percent=0.16666666666666666 (a,b)=(0,3) -> (101,104)
            #   avg_local = [ 0.00000000e+00 -5.55111512e-17  0.00000000e+00]
            #   p1_local  = [ 0.  -0.5  0. ]
            #   p2_local  = [0.  2.5 0. ]
            # ------
            # rods = (rod_eid_nodes, rod_nids, rod_interp_local_xyz)
            # rod_eid_nodes:
            # [[ 10   1   2]
            #  [-10   5   6]]
            # rod_nids:
            # [1 2 5 6]
            # rod_xyzs:
            # [[ 3.00000000e+00 -5.55111512e-17  0.00000000e+00]
            #  [ 5.00000000e-01 -5.55111512e-17  0.00000000e+00]
            #  [ 5.00000000e-01 -5.55111512e-17  0.00000000e+00]
            #  [ 0.00000000e+00 -5.55111512e-17  0.00000000e+00]]
            #
            #--------------------------------------------------------
            # some stuff seems wrong...
            # [source_eid, new_nid, source_nid1, source_nid2]
            # geometry - seems like it's missing a -10 line?
            # array([[ 10,   2, 101, 103],
            #        [-10,   6, 101, 104],
            #        [ 10,   1, 102, 103]], dtype=int32)]
            # rod_eids=[
            #  [ 10   1   2]
            #  [-10   5   6]]
            # rod_xyz_global:
            # [1 2 5 6]
            # rod_xyz_local:
            # [[ 3.00000000e+00 -5.55111512e-17  0.00000000e+00]
            #  [ 5.00000000e-01 -5.55111512e-17  0.00000000e+00]
            #  [ 5.00000000e-01 -5.55111512e-17  0.00000000e+00]
            #  [ 0.00000000e+00 -5.55111512e-17  0.00000000e+00]]

            x_vector = [0., 0., 1.]
            moi_data = cut_and_plot_moi(
                model, normal_plane, log,
                ystations, coords, x_vector,
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
                debug_v3=False,
                # debug_v3=False,
            )
            out_dict, plane_bdf_filenames1, plane_bdf_filenames2, ifig = moi_data
            (y, L, A, I, J,
             ExI, EyI, GJ, avg_centroid) = unpack_moi(out_dict)

            Ex = ExI[:, 0] / I[:, 0]
            Ey = EyI[:, 0] / I[:, 0]
            # GJ is now a real torsion constant, so it is NOT G*J with J the
            # polar moment; dividing by J no longer recovers G.  Divide by the
            # open-section constant that GJ is actually built from instead.
            G = GJ / j_open
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
            assert np.allclose(L, L_expected), f'L={L} expected={L_expected}'
            assert np.allclose(A, A_expected), f'A={A} expected={A_expected}'
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
        # open section: one straight wall, so j = s*t^3/3 rather than the polar
        # moment.  Old baseline was G*(Ix+Iz) = 1081730.77, 93.75x too stiff.
        G = E / (2. * (1. + 0.3))  # 11538461.54
        GJ_expected = [G * cut_length * t ** 3 / 3.]  # 11538.46
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

        x_vector = [0., 1., 0.]
        moi_data = cut_and_plot_moi(
            model, normal_plane, log,
            ystations, coords, x_vector,
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
        (out_dict, plane_bdf_filenames1, plane_bdf_filenames2, ifig) = moi_data
        (y, L, A, I, J,
         ExI, EyI, GJ, avg_centroid) = unpack_moi(out_dict)

        if IS_PANDAS:
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
        assert np.allclose(GJ, GJ_expected), (GJ_expected, GJ.tolist())
        assert np.allclose(avg_centroid, centroid_expected), avg_centroid.tolist()
        del model.properties[pid]

    def test_cut_ellipse_constant_area(self):
        """
        Prismatic 2:1 elliptical tube extruded along +y; ``cut_and_plot_moi``
        plus the equivalent beam model it writes.

        The section does not change along the span, so every station must
        return the *same* area / inertia / centroid / GJ.  That is a strong
        self-check that needs no closed-form value at all.

        The 2:1 aspect ratio then makes the 1-2 axis mapping unambiguous.  The
        major axis is along global x and the CBEAM v-vector is [1, 0, 0], so
        ``y_elem`` lies along +x and ``I1 = int(y_elem^2 dA)`` must be the
        LARGER of the two.  A circle or a square could not tell these apart.

        Finally, GJ must be the Bredt-Batho closed-cell value; the polar
        moment ``G*(Ix+Iz)`` is ~1.5x too stiff for a 2:1 ellipse.
        """
        dirname = TEST_PATH
        tag = 'ellipse_'
        log = SimpleLogger(level='warning', encoding='utf-8')

        a, b, t = 20., 10., 0.1
        E, nu = 1.0e7, 0.3
        G = E / (2. * (1. + nu))
        span, nspan, ntheta = 100., 40, 120

        model, pts = _build_ellipse_tube(
            log, a, b, t, span, nspan, ntheta, E, nu)
        model.write_bdf(dirname / 'ellipse.bdf')
        exact = _thin_wall_section(pts, t)

        # deliberately off the node planes (span/nspan = 2.5), so the cut has
        # to interpolate rather than land on coincident grids
        ystations = np.array([21.3, 33.7, 46.1, 58.5, 70.9])
        coords = [CORD2R(1000 + i, rid=0, origin=[0., ys, 0.],
                         zaxis=[0., ys, 1.], xzplane=[1., ys, 0.])
                  for i, ys in enumerate(ystations)]
        normal_plane = coords[0].j
        assert np.allclose(normal_plane, [0., 1., 0.]), normal_plane

        beam_bdf_filename = tag + 'equivalent_beam_model.bdf'
        x_vector = [0., 0., 1.]
        moi_data = cut_and_plot_moi(
            model, normal_plane, log, ystations, coords, x_vector,
            dirname=dirname, plot=False, show=False, stop_on_failure=True,
            xyz_round=3,
            area_round=3,
            inertia_round=2,
            cut_data_span_filename='',
            beam_model_bdf_filename=beam_bdf_filename,
            thetas_csv_filename=tag + 'thetas.csv')
        (out_dict, plane_bdf_filenames1, plane_bdf_filenames2, unused_ifig) = moi_data
        (y, L, A, I, J, ExI, EyI, GJ, avg_centroid) = unpack_moi(out_dict)

        assert np.isfinite(A).all(), f'missing cuts: A={A}'
        assert np.allclose(y, ystations), y

        # ------------------------------------------------------------------
        # 1) prismatic: nothing varies along the span
        # ------------------------------------------------------------------
        # the cut interpolates along the element edges, so a station that does
        # not land on a node plane picks up a little roundoff; ~1e-5 is the
        # observed spread, so 1e-4 flags a real span dependence
        for name, value in [('L', L), ('A', A), ('J', J), ('GJ', GJ)]:
            assert np.allclose(value, value[0], rtol=1e-4), \
                f'{name} varies vs span: {value}'
        for name, value in [('I', I), ('ExI', ExI), ('EyI', EyI)]:
            atol = 1e-4 * np.abs(value).max()
            assert np.allclose(value, value[0, :], rtol=1e-4, atol=atol), \
                f'{name} varies vs span:\n{value}'
        # only the y column of the centroid varies; it *is* the station
        assert np.allclose(avg_centroid[:, 1], ystations), avg_centroid

        # ------------------------------------------------------------------
        # 2) section integrals vs the closed-form thin-wall values
        # ------------------------------------------------------------------
        # area and perimeter are integrated exactly, whatever the subdivision
        assert np.allclose(A, exact['A'], rtol=1e-6), (A[0], exact['A'])
        assert np.allclose(L, exact['perimeter'], rtol=1e-6), (L[0], exact['perimeter'])

        # lump <= cut <= strip; see _thin_wall_section
        eps = 1e-9
        assert (exact['int_x2_lump'] * (1. - eps) <= I[:, 0]).all() and \
               (I[:, 0] <= exact['int_x2_strip'] * (1. + eps)).all(), \
            (I[:, 0], exact['int_x2_lump'], exact['int_x2_strip'])
        assert (exact['int_z2_lump'] * (1. - eps) <= I[:, 2]).all() and \
               (I[:, 2] <= exact['int_z2_strip'] * (1. + eps)).all(), \
            (I[:, 2], exact['int_z2_lump'], exact['int_z2_strip'])

        # doubly symmetric, so the product of inertia is numerical noise
        assert np.abs(I[:, 5]).max() < 1e-6 * exact['int_x2_lump'], I[:, 5]
        assert np.allclose(avg_centroid[:, 0], 0., atol=1e-6 * a), avg_centroid
        assert np.allclose(avg_centroid[:, 2], 0., atol=1e-6 * a), avg_centroid
        assert np.allclose(ExI[:, 0] / I[:, 0], E, rtol=1e-6), ExI[0, 0] / I[0, 0]

        # the whole point of a 2:1 section: int(x^2) must dominate int(z^2).
        # for a thin elliptical shell the ratio is ~2.9, not (a/b)**2 = 4,
        # because the wall is not uniformly distributed in x
        ratio = I[0, 0] / I[0, 2]
        assert 2.5 < ratio < 3.5, ratio
        assert np.all(I[:, 0] > I[:, 2]), (I[:, 0], I[:, 2])

        # ------------------------------------------------------------------
        # 3) torsion is Bredt-Batho, not the polar moment
        # ------------------------------------------------------------------
        assert np.allclose(GJ, G * exact['J'], rtol=1e-3), (GJ[0], G * exact['J'])
        gj_polar = G * (exact['int_x2_lump'] + exact['int_z2_lump'])
        assert gj_polar > 1.3 * GJ[0], (gj_polar, GJ[0])

        # ------------------------------------------------------------------
        # 4) the equivalent beam deck reproduces every stiffness
        # ------------------------------------------------------------------
        beam_model = read_bdf(dirname / beam_bdf_filename, punch=True, debug=None)
        assert len(beam_model.nodes) == len(ystations), beam_model.nodes
        assert len(beam_model.elements) == len(ystations) - 1, beam_model.elements

        mat = beam_model.materials[1]
        # the MAT1 carries real moduli so that rho*A is a meaningful mass and
        # K*G*A is a meaningful shear stiffness
        assert np.allclose(mat.e, 2. * mat.g * (1. + mat.nu)), (mat.e, mat.g, mat.nu)
        assert np.allclose(mat.e, E, rtol=1e-6), mat.e
        assert np.allclose(mat.g, G, rtol=1e-6), mat.g

        for pid, prop in sorted(beam_model.properties.items()):
            assert np.allclose(prop.A[0], exact['A'], rtol=1e-4), (pid, prop.A)
            assert np.allclose(mat.e * prop.A[0], E * exact['A'], rtol=1e-4)
            # I1 = int(y_e^2 dA), I2 = int(z_e^2 dA) in the *element* frame.
            # The beam runs along +y and x_vector = [0, 0, 1], so
            # y_e = +z_global and z_e = +x_global: I1 picks up int(z^2) and
            # I2 picks up int(x^2), i.e. the opposite of the cut-plane order.
            assert np.allclose(mat.e * prop.i1[0], ExI[0, 2], rtol=1e-4), \
                (pid, mat.e * prop.i1[0], ExI[0, 2])
            assert np.allclose(mat.e * prop.i2[0], ExI[0, 0], rtol=1e-4), \
                (pid, mat.e * prop.i2[0], ExI[0, 0])
            assert np.allclose(mat.g * prop.j[0], G * exact['J'], rtol=1e-3)
            # the ellipse is 2:1 with the long axis along x = z_e, so I2 is
            # the strong one.  Swapping x_vector to [1,0,0] swaps these.
            assert prop.i2[0] > prop.i1[0], (pid, prop.i1, prop.i2)
            k1 = 1.0 if prop.k1 is None else prop.k1
            assert np.allclose(k1 * mat.g * prop.A[0], G * exact['A'], rtol=1e-4)

        # plot=False and cut_data_span_filename='', so only these two exist;
        # _cleanup_moi_files() would trip over the missing plots
        for fname in plane_bdf_filenames1 + plane_bdf_filenames2:
            os.remove(fname)
        # os.remove(dirname / beam_bdf_filename)
        os.remove(dirname / (tag + 'thetas.csv'))

    def test_cut_ellipse_fuselage_frame(self):
        """
        The equivalent beam model must be written in the BASIC frame, even
        when the cut coord is rotated relative to it.

        ``avg_centroid`` comes out of the cutter in the cut coord's LOCAL
        frame (the plots and the csv want in-plane coordinates), and column 1
        is then overwritten with the station.  For a wing cut that happens to
        be the basic frame, because the coord is built so its axes coincide
        with the global ones -- so the bug is invisible there.  A fuselage cut
        marches along +x with a rotated coord, and the GRIDs used to come out
        permuted as ``[ycg, station, zcg]``.

        The section is deliberately centered off-axis so that a permutation
        cannot hide behind a zero.
        """
        dirname = TEST_PATH
        tag = 'fuse_ellipse_'
        log = SimpleLogger(level='warning', encoding='utf-8')

        a, b, t = 20., 10., 0.1
        E, nu = 1.0e7, 0.3
        ycg, zcg = 3., 7.
        span, nspan, ntheta = 100., 40, 80

        model, unused_pts = _build_ellipse_tube(
            log, a, b, t, span, nspan, ntheta, E, nu,
            axis=0, center=(ycg, zcg))

        # the fuselage recipe: march along +x, cut in the global yz plane.
        # the local axes are i=+y, j=-x, k=+z, so the coord IS rotated.
        xstations = np.array([21.3, 46.1, 70.9])
        origin = np.array([0., 0., 0.])
        zaxis = np.array([0., 0., 1.])
        xzplane = np.array([0., 1., 0.])
        coords = []
        for i, xs in enumerate(xstations):
            dxyz = np.array([xs, 0., 0.])
            coords.append(CORD2R(2000 + i, rid=0, origin=origin + dxyz,
                                 zaxis=zaxis + dxyz, xzplane=xzplane + dxyz))
        normal_plane = coords[0].j
        assert np.allclose(normal_plane, [-1., 0., 0.]), normal_plane

        beam_bdf_filename = tag + 'equivalent_beam_model.bdf'
        x_vector = [0., 0., 1.]
        moi_data = cut_and_plot_moi(
            model, normal_plane, log, xstations, coords, x_vector,
            dirname=dirname, plot=False, show=False, stop_on_failure=True,
            cut_data_span_filename='',
            beam_model_bdf_filename=beam_bdf_filename,
            thetas_csv_filename=tag + 'thetas.csv')
        (out_dict, plane_bdf_filenames1, plane_bdf_filenames2, unused_ifig) = moi_data
        (unused_x, unused_L, A, unused_I, unused_J, unused_ExI, unused_EyI,
         unused_GJ, avg_centroid) = unpack_moi(out_dict)
        assert np.isfinite(A).all(), f'missing cuts: A={A}'

        # the beam GRIDs are in the basic frame
        beam_model = read_bdf(dirname / beam_bdf_filename, punch=True, debug=None)
        xyz = np.array([beam_model.nodes[nid].xyz
                        for nid in sorted(beam_model.nodes)])
        xyz_expected = np.column_stack([
            xstations,
            np.full(len(xstations), ycg),
            np.full(len(xstations), zcg)])
        assert np.allclose(xyz, xyz_expected, atol=1e-6), \
            f'beam GRIDs are not in the basic frame:\n{xyz}\nexpected\n{xyz_expected}'

        # ...while the reported avg_centroid stays in the local frame, which
        # is what plot_inertia and the csv header assume
        local_expected = np.column_stack([
            np.full(len(xstations), ycg),
            xstations,
            np.full(len(xstations), zcg)])
        assert np.allclose(avg_centroid, local_expected, atol=1e-6), avg_centroid

        # the CBEAM axis must run down the fuselage, not across it
        for eid in sorted(beam_model.elements):
            elem = beam_model.elements[eid]
            n1, n2 = elem.node_ids
            dxyz = beam_model.nodes[n2].xyz - beam_model.nodes[n1].xyz
            assert abs(dxyz[0]) > 1e-6, (eid, dxyz)
            assert np.allclose(dxyz[1:], 0., atol=1e-6), (eid, dxyz)

        for fname in plane_bdf_filenames1 + plane_bdf_filenames2:
            os.remove(fname)
        os.remove(dirname / beam_bdf_filename)
        os.remove(dirname / (tag + 'thetas.csv'))

    def test_cut_ellipse_beam_grid_offset(self):
        """
        cut_and_plot_moi with prescribed beam GRID locations.

        The beam nodes are put on an arbitrary straight reference axis (stand-in
        for a set of load control points) rather than on the section centroids,
        and the difference is carried on the CBEAM WA/WB offsets.  A CBEAM's
        element axis runs from GA+WA to GB+WB, so the emitted model has to be
        geometrically identical to the un-offset one -- same element axis, same
        section properties -- with only the GRIDs moved.

        """
        dirname = TEST_PATH
        log = SimpleLogger(level='warning', encoding='utf-8')
        tag = 'bgo_'
        a, b, t = 20., 10., 0.1
        span, nspan, ntheta = 100., 20, 40
        E, nu = 1.0e7, 0.3
        # section deliberately off the y-axis so the offsets are nonzero
        xcg, zcg = 4., -6.
        model, unused_pts = _build_ellipse_tube(
            log, a, b, t, span, nspan, ntheta, E, nu,
            axis=1, center=(xcg, zcg))

        ystations = np.array([20., 50., 80.])
        coords = [CORD2R(4000 + i, rid=0, origin=[0., ys, 0.],
                         zaxis=[0., ys, 1.], xzplane=[1., ys, 0.])
                  for i, ys in enumerate(ystations)]
        normal_plane = coords[0].j
        x_vector = [0., 0., 1.]
        nstation = len(ystations)

        # a straight reference axis that is nowhere near the centroid
        ref_xyz = np.column_stack([
            np.full(nstation, 15.),
            ystations,
            np.full(nstation, 30.)])
        ref_ids = np.arange(100001, 100001 + nstation)

        kwargs = dict(
            dirname=dirname, plot=False, show=False, face_data=None,
            stop_on_failure=True, cut_data_span_filename='',
            thetas_csv_filename=tag + 'thetas.csv')

        base_bdf_filename = tag + 'base.bdf'
        off_bdf_filename = tag + 'offset.bdf'
        out_base = cut_and_plot_moi(
            model, normal_plane, log, ystations, coords, x_vector,
            beam_model_bdf_filename=base_bdf_filename, **kwargs)
        out_off = cut_and_plot_moi(
            model, normal_plane, log, ystations, coords, x_vector,
            beam_grid_xyz=ref_xyz, beam_grid_ids=ref_ids, beam_id0=500,
            beam_model_bdf_filename=off_bdf_filename, **kwargs)

        base = read_bdf(dirname / base_bdf_filename, punch=True, debug=None)
        off = read_bdf(dirname / off_bdf_filename, punch=True, debug=None)

        # the GRIDs went exactly where they were told, under the given ids
        assert sorted(off.nodes) == ref_ids.tolist(), sorted(off.nodes)
        xyz = np.array([off.nodes[nid].xyz for nid in ref_ids])
        assert np.allclose(xyz, ref_xyz, atol=0., rtol=0.), xyz

        # ...and the default is still a GRID on each centroid, numbered 1..n
        assert sorted(base.nodes) == [1, 2, 3], sorted(base.nodes)
        xyz_base = np.array([base.nodes[nid].xyz for nid in sorted(base.nodes)])
        centroid_expected = np.column_stack([
            np.full(nstation, xcg), ystations, np.full(nstation, zcg)])
        assert np.allclose(xyz_base, centroid_expected, atol=1e-6), xyz_base

        # beam_id0 moved the elements/properties/material out of the way
        assert sorted(off.elements) == [500, 501], sorted(off.elements)
        assert sorted(off.properties) == [500, 501], sorted(off.properties)
        assert sorted(off.materials) == [500], sorted(off.materials)
        assert sorted(base.elements) == [1, 2], sorted(base.elements)

        for eid_base, eid_off in zip(sorted(base.elements), sorted(off.elements)):
            elem_base = base.elements[eid_base]
            elem_off = off.elements[eid_off]

            # GA+WA / GB+WB reproduce the un-offset element axis exactly
            for iend in (0, 1):
                w = elem_off.wa if iend == 0 else elem_off.wb
                xyz_off = (np.asarray(off.nodes[elem_off.node_ids[iend]].xyz) +
                           np.asarray(w))
                xyz_cen = np.asarray(base.nodes[elem_base.node_ids[iend]].xyz)
                assert np.allclose(xyz_off, xyz_cen, atol=1e-9), \
                    f'eid={eid_off:d} end={iend:d}: {xyz_off} != {xyz_cen}'

            # cutting at the station means the offset is purely in-plane;
            # an axial component would stretch the element
            assert abs(elem_off.wa[1]) < 1e-10, elem_off.wa
            assert abs(elem_off.wb[1]) < 1e-10, elem_off.wb

            # moving the nodes must not touch the section properties
            prop_base = base.properties[elem_base.pid]
            prop_off = off.properties[elem_off.pid]
            for field in ('A', 'i1', 'i2', 'i12', 'j', 'k1', 'k2'):
                assert np.allclose(getattr(prop_base, field),
                                   getattr(prop_off, field)), field

        # a prescribed point that is the wrong shape is a caller error
        with self.assertRaises(ValueError):
            cut_and_plot_moi(
                model, normal_plane, log, ystations, coords, x_vector,
                beam_grid_xyz=ref_xyz[:-1],
                beam_model_bdf_filename=off_bdf_filename, **kwargs)
        with self.assertRaises(ValueError):
            cut_and_plot_moi(
                model, normal_plane, log, ystations, coords, x_vector,
                beam_grid_xyz=ref_xyz, beam_grid_ids=np.array([7, 7, 8]),
                beam_model_bdf_filename=off_bdf_filename, **kwargs)

        # both runs write the same plane_face_* names, so dedupe before unlink
        for fname in set(out_base[1] + out_base[2] + out_off[1] + out_off[2]):
            os.remove(fname)
        os.remove(dirname / base_bdf_filename)
        os.remove(dirname / off_bdf_filename)
        os.remove(dirname / (tag + 'thetas.csv'))

    def test_cut_ellipse_element_frame(self):
        """
        I1/I2/I12 are written in the CBEAM element frame, not the cut frame.

        MSC defines I1 = int(y_e^2 dA), I2 = int(z_e^2 dA) and
        I12 = int(y_e*z_e dA) about the *element* axes, and those axes come
        from the orientation vector::

            x_e = GA+WA -> GB+WB        y_e = v normal to x_e        z_e = x_e cross y_e

        For a beam along +y that means v = [1,0,0] gives y_e = +x_global and
        z_e = -z_global, while v = [0,0,1] gives y_e = +z_global and
        z_e = +x_global.  Same section, same cut, same everything else -- so
        the two runs must come out with I1 and I2 exchanged and I12 negated.
        This used to be written straight out of the cut-plane integrals with
        no regard for v at all, which silently transposed the bending axes
        and flipped the product of inertia.

        A tilted ellipse is used so that all three moments are distinct and
        I12 is comfortably non-zero.
        """
        dirname = TEST_PATH
        log = SimpleLogger(level='warning', encoding='utf-8')
        tag = 'efr_'
        a, b, t = 20., 10., 0.1
        span, nspan, ntheta = 100., 20, 40
        E, nu = 1.0e7, 0.3
        model, unused_pts = _build_ellipse_tube(
            log, a, b, t, span, nspan, ntheta, E, nu, axis=1, tilt=0.4)

        ystations = np.array([30., 50., 70.])
        coords = [CORD2R(4200 + i, rid=0, origin=[0., ys, 0.],
                         zaxis=[0., ys, 1.], xzplane=[1., ys, 0.])
                  for i, ys in enumerate(ystations)]
        normal_plane = coords[0].j
        kwargs = dict(
            dirname=dirname, plot=False, show=False, face_data=None,
            stop_on_failure=True, cut_data_span_filename='',
            thetas_csv_filename=tag + 'thetas.csv')

        outs, props = {}, {}
        for name, x_vector in (('x', [1., 0., 0.]), ('z', [0., 0., 1.])):
            bdf_filename = f'{tag}{name}.bdf'
            outs[name] = cut_and_plot_moi(
                model, normal_plane, log, ystations, coords, x_vector,
                beam_model_bdf_filename=bdf_filename, **kwargs)
            beam_model = read_bdf(dirname / bdf_filename, punch=True, debug=None)
            props[name] = beam_model.properties[min(beam_model.properties)]

        px, pz = props['x'], props['z']

        # the section is tilted, so nothing here is degenerate
        assert not np.allclose(px.i1[0], px.i2[0]), (px.i1, px.i2)
        assert abs(px.i12[0]) > 0.05 * abs(px.i1[0]), px.i12

        # v = [1,0,0] vs v = [0,0,1]: the 1 and 2 axes trade places.
        # rtol is 1e-5 rather than exact because these come back through an
        # 8-character small-field BDF, and the leading minus sign on I12
        # costs it a significant digit relative to its positive twin.
        assert np.allclose(pz.i1[0], px.i2[0], rtol=1e-5), (pz.i1, px.i2)
        assert np.allclose(pz.i2[0], px.i1[0], rtol=1e-5), (pz.i2, px.i1)
        # ...and the product of inertia changes sign with the handedness
        assert np.allclose(pz.i12[0], -px.i12[0], rtol=1e-5), (pz.i12, px.i12)

        # everything that does not depend on the orientation vector is unmoved
        for field in ('A', 'j', 'k1', 'k2'):
            assert np.allclose(getattr(px, field), getattr(pz, field)), \
                (field, getattr(px, field), getattr(pz, field))

        # I1 + I2 is the polar moment, which no rotation can change
        assert np.allclose(px.i1[0] + px.i2[0], pz.i1[0] + pz.i2[0],
                           rtol=1e-10)

        # v parallel to the beam axis leaves plane 1 undefined
        with self.assertRaises(ValueError):
            cut_and_plot_moi(
                model, normal_plane, log, ystations, coords, [0., 1., 0.],
                beam_model_bdf_filename=f'{tag}bad.bdf', **kwargs)

        for fname in set(sum((outs[k][1] + outs[k][2] for k in outs), [])):
            os.remove(fname)
        for name in ('x', 'z'):
            os.remove(dirname / f'{tag}{name}.bdf')
        os.remove(dirname / (tag + 'thetas.csv'))

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

        x_vector = [0., 0., 1.]
        moi_data = cut_and_plot_moi(
            model, normal_plane, log,
            ystations, coords, x_vector,
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
        out_dict, plane_bdf_filenames1, plane_bdf_filenames2, ifig = moi_data
        (y, L, A, I, J,
         ExI, EyI, GJ, avg_centroid) = unpack_moi(out_dict)
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
        # open section: one straight wall, so j = s*t^3/3 rather than the polar
        # moment.  Old baseline was G*(Ix+Iz) = 1081730.77, 93.75x too stiff.
        G = E / (2. * (1. + 0.3))  # 11538461.54
        GJ_expected = [G * cut_length * t ** 3 / 3.]  # 11538.46
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

        x_vector = [0., 0., 1.]
        moi_data = cut_and_plot_moi(
            model, normal_plane, log,
            ystations, coords, x_vector,
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
        out_dict, plane_bdf_filenames1, plane_bdf_filenames2, ifig = moi_data
        (y, L, A, I, J,
         ExI, EyI, GJ, avg_centroid) = unpack_moi(out_dict)
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
        # open section: one straight wall, so j = s*t^3/3 rather than the polar
        # moment.  Old baseline was G*(Ix+Iz) = 1081730.77, 93.75x too stiff.
        G = E / (2. * (1. + 0.3))  # 11538461.54
        GJ_expected = [G * cut_length * t ** 3 / 3.]  # 11538.46
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
        x_vector = [0., 0., 1.]

        with self.assertRaises(NotImplementedError):
            moi_data = cut_and_plot_moi(
                model, normal_plane, log,
                ystations, coords, x_vector,
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

        x_vector = [0., 0., 1.]
        if run_y_cuts:
            log.info('working on y-cuts')
            moi_data = cut_and_plot_moi(
                bdf_filename, normal_plane, log,
                ystations, coords, x_vector,
                dirname=dirname,
                plot=IS_MATPLOTLIB, show=False, face_data=face_data,
                cut_data_span_filename='y_cut_data_vs_span.csv',
                beam_model_bdf_filename='y_equivalent_beam_model.bdf',
                thetas_csv_filename='y_thetas.csv',
                normalized_inertia_png_filename='y_normalized_inertia_vs_span.png',
                area_span_png_filename='y_area_vs_span.png',
                amoi_span_png_filename='y_amoi_vs_span.png',
                e_amoi_span_png_filename='y_e_amoi_vs_span.png',
                centroid_span_png_filename='y_centroid_vs_span.png',
                debug_vectorize=False,
            )
            out_dict, plane_bdf_filenames1, plane_bdf_filenames2, ifig = moi_data
            (y, L, A, I, J,
             ExI, EyI, GJ, avg_centroid) = unpack_moi(out_dict)
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
            # Bredt-Batho multi-cell torsion, replacing the old G*(Ix+Iz) polar
            # baseline.  The wing box is a closed multi-cell section, so the
            # polar moment badly overstated it -- 42x at the root, 27x at the
            # tip.  These are regression values, not closed-form ones; the
            # closed-form check lives in test_cut_ellipse_constant_area.
            GJ_expected = [406405373109638.2, 86481383450659.14, 55375491645159.75,
                           6269981882883.581, 1178856065913.1702, 601070280743.0763,
                           284866251513.561, 122035533137.60596, 46849456050.85022,
                           16802569331.901482]
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
                xstations, xcoords, x_vector,
                dirname=dirname, ifig=10,
                plot=IS_MATPLOTLIB, show=False, face_data=face_data,
                cut_data_span_filename='x_cut_data_vs_span.csv',
                beam_model_bdf_filename='x_equivalent_beam_model.bdf',
                thetas_csv_filename='x_thetas.csv',
                normalized_inertia_png_filename='x_normalized_inertia_vs_span.png',
                area_span_png_filename='x_area_vs_span.png',
                amoi_span_png_filename='x_amoi_vs_span.png',
                e_amoi_span_png_filename='x_e_amoi_vs_span.png',
                centroid_span_png_filename='x_centroid_vs_span.png',
                debug_vectorize=True,
            )
            out_dict, plane_bdf_filenames1, plane_bdf_filenames2, ifig = moi_data
            (x, L, A, I, J,
             ExI, EyI, GJ, avg_centroid) = unpack_moi(out_dict)
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

            # Bredt-Batho, replacing the old G*(Ix+Iz) polar baseline.  This
            # assertion was previously unreachable -- the y-cut GJ check above
            # failed first, so this list was never exercised and silently went
            # stale too.  Regression values, not closed-form.
            GJ_expected = [1909997028.546258, 143618441147.58948, 253466248971.01608,
                  326852460412.9585, 385726808364.47565, 452019448660.6323,
                  915364684102.5784, 2018834400232.94, 4172528234519.123,
                  7309103592385.41, 10844932249461.475, 13291126771296.396,
                  14782088705379.807, 14372117610030.074, 32006475322096.867,
                  22739867489634.14, 13707972167928.145, 11723023385219.502,
                  2400531926446.8584]

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
        # os.remove(dirname / f'{tag}area_vs_span.png')
        os.remove(dirname / f'{tag}amoi_vs_span.png')
        os.remove(dirname / f'{tag}e_amoi_vs_span.png')
        os.remove(dirname / f'{tag}centroid_vs_span.png')
        os.remove(dirname / f'{tag}equivalent_beam_model.bdf')
        # os.remove(dirname / f'{tag}thetas.csv')
        # 'y_amoi_vs_span.png', 'y_bwb_saero.bdf',
        #
        # 'y_centroid_vs_span.png',
        # 'y_e_amoi_vs_span.png',
        # 'y_normalized_inertia_vs_span.png',

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
    # [101, 102, 103]
    #
    model.add_cquad4(10, 11, [101, 102, 103, 104], zoffset=zoffset)
    coord = model.add_cord2r(
        cid=1,
        origin=[0., dy, 0.],
        zaxis=[0., dy, 1.],
        xzplane=[1., dy, 0.])
    return model, coord

def _ellipse_pts(a: float, b: float, ntheta: int) -> np.ndarray:
    """midline of an ellipse; ``a`` is along x, ``b`` is along z"""
    theta = np.linspace(0., 2*np.pi, ntheta, endpoint=False)
    return np.column_stack([a*np.cos(theta), b*np.sin(theta)])


def _thin_wall_section(pts: np.ndarray, t: float) -> dict[str, float]:
    """
    Closed-form thin-wall section properties for a faceted closed midline.

    Two versions of the second moments are returned, and the cutter must land
    between them:

    - ``*_lump`` lumps each facet's area at its own centroid
    - ``*_strip`` integrates x^2 continuously along each facet

    The cutter triangulates every CQUAD4 before cutting, so each facet comes
    back as two sub-segments lumped at their own centroids.  Refining a
    midpoint rule on the convex integrand x^2 always moves the answer up
    toward the exact strip value without overshooting it, so
    ``lump <= cut <= strip`` is guaranteed and is a much sharper statement
    than any hand-tuned tolerance.
    """
    p1 = pts
    p2 = np.roll(pts, -1, axis=0)
    x1, z1 = p1[:, 0], p1[:, 1]
    x2, z2 = p2[:, 0], p2[:, 1]
    ell = np.hypot(x2 - x1, z2 - z1)
    dA = t * ell

    area = dA.sum()
    xc = (dA * (x1 + x2) / 2).sum() / area
    zc = (dA * (z1 + z2) / 2).sum() / area
    xm = (x1 + x2) / 2 - xc
    zm = (z1 + z2) / 2 - zc

    # Bredt-Batho for the single closed cell, plus the open-section term
    area_enclosed = 0.5 * abs((x1 * z2 - x2 * z1).sum())
    j_bredt = 4. * area_enclosed ** 2 / (ell / t).sum()
    j_open = (ell * t ** 3).sum() / 3.
    return {
        'A': area, 'xc': xc, 'zc': zc,
        'int_x2_lump': (dA * xm ** 2).sum(),
        'int_z2_lump': (dA * zm ** 2).sum(),
        'int_x2_strip': (dA * ((x1 - xc) ** 2 + (x1 - xc) * (x2 - xc) +
                               (x2 - xc) ** 2) / 3).sum(),
        'int_z2_strip': (dA * ((z1 - zc) ** 2 + (z1 - zc) * (z2 - zc) +
                               (z2 - zc) ** 2) / 3).sum(),
        'int_xz': (dA * xm * zm).sum(),
        'perimeter': ell.sum(),
        'J': j_bredt + j_open,
    }


def _build_ellipse_tube(log: SimpleLogger,
                        a: float, b: float, t: float,
                        span: float, nspan: int, ntheta: int,
                        E: float, nu: float,
                        pid: int=11, mid: int=12,
                        axis: int=1,
                        center: tuple[float, float]=(0., 0.),
                        tilt: float=0.,
                        ) -> tuple[BDF, np.ndarray]:
    """
    Prismatic elliptical tube extruded along ``+axis``.

    ^ z
    |    _____
    |  /       \\
    | (    +    )  --> x     a by b, constant along the span
    |  \\ _____ /

    ``axis=1`` extrudes along +y and puts ``a`` along x and ``b`` along z
    (a wing cut).  ``axis=0`` extrudes along +x and puts ``a`` along y and
    ``b`` along z (a fuselage cut).  ``center`` offsets the section in those
    same two in-plane directions, which is what makes a frame error visible.
    ``tilt`` (radians) rotates the section about the extrusion axis, which
    gives it a non-zero product of inertia -- needed to see an I12 error at
    all, since an untilted ellipse has I12 = 0 whatever the convention.
    """
    model = BDF(log=log)
    model.add_mat1(mid, E=E, G=None, nu=nu)
    model.add_pshell(pid, mid1=mid, t=t, mid2=mid, mid3=mid)

    iaxes = [i for i in range(3) if i != axis]
    pts = _ellipse_pts(a, b, ntheta)
    if tilt:
        c, s = np.cos(tilt), np.sin(tilt)
        pts = pts @ np.array([[c, s], [-s, c]])
    pts = pts + np.asarray(center, dtype='float64')
    stations = np.linspace(0., span, nspan + 1)

    def nid(itheta: int, istation: int) -> int:
        return istation * ntheta + (itheta % ntheta) + 1

    for istation, station in enumerate(stations):
        for itheta, (p, q) in enumerate(pts):
            xyz = np.zeros(3, dtype='float64')
            xyz[axis] = station
            xyz[iaxes[0]] = p
            xyz[iaxes[1]] = q
            model.add_grid(nid(itheta, istation), xyz)

    eid = 1
    for istation in range(nspan):
        for itheta in range(ntheta):
            model.add_cquad4(eid, pid, [
                nid(itheta, istation), nid(itheta+1, istation),
                nid(itheta+1, istation+1), nid(itheta, istation+1)])
            eid += 1
    model.cross_reference()
    return model, pts


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
