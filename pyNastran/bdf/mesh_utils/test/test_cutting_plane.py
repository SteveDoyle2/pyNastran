"""defines cutting plane tests"""
import os
from itertools import count
from typing import Tuple, Any
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
    split_to_trias, calculate_area_moi)
from pyNastran.bdf.mesh_utils.cutting_plane_plotter import cut_and_plot_model
#from pyNastran.bdf.mesh_utils.bdf_merge import bdf_merge
from pyNastran.op2.op2_geom import read_op2_geom

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, '..', 'models')


class TestCuttingPlane(unittest.TestCase):
    """various cutting plane tests"""
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

    def test_cut_bwb(self):
        """recover element ids"""
        log = SimpleLogger(level='warning', encoding='utf-8', log_func=None)
        is_bwb = True
        if is_bwb:
            bdf_filename = os.path.join(MODEL_PATH, 'bwb', 'bwb_saero.bdf')  # ymax~=1262.0
        else:  #  pragma: no cover
            bdf_filename = r'C:\NASA\asm\all_modes_mach_0.85\flutter.bdf'  # ymax=1160.601
        normal_plane = np.array([0., 1., 0.])
        y, A, I, EI, avg_centroid, plane_bdf_filenames = cut_and_plot_moi(
            bdf_filename, normal_plane, log,
            plot=True, show=True)

        show = True
        #show = False
        if IS_MATPLOTLIB:
            plot_inertia(y, A, I, EI, avg_centroid, show=show)
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


def cut_and_plot_moi(bdf_filename: str, normal_plane: np.ndarray, log: SimpleLogger,
                     ytol: float=2.0, ncuts: int=2000, dirname: str='',
                     plot: bool=True, show: bool=False) -> Tuple[Any, Any, Any, Any, Any]: # y, A, I, EI, avg_centroid
    model = read_bdf(bdf_filename, log=log)
    model2 = read_bdf(bdf_filename, log=log)

    # initialize theta
    thetas = {}
    for eid in model.elements:
        #  theta, Ex, Ey, Gxy
        thetas[eid] = (0., 0., 0., 0.)

    #p1 = np.array([466.78845, 735.9053, 0.0])
    #p2 = np.array([624.91345, 639.68896, -0.99763656])
    #dx = p2 - p1
    ytol = 2.
    nodal_result = None
    plane_bdf_filenames = []
    y = []
    A = []
    I = []
    EI = []
    avg_centroid = []

    is_bwb = True
    for i in range(ncuts):
        if is_bwb:
            dy = 100. * i + 1.  #  bwb
            coord = CORD2R(1, rid=0, origin=[0., dy, 0.], zaxis=[0., dy, 1], xzplane=[1., dy, 0.])
        else:  #  pragma: no cover
            dy = 4. * i + 1.  #  CRM
            coord = CORD2R(1, rid=0, origin=[0., dy, 0.], zaxis=[0., dy, 1], xzplane=[1., dy, 0.])
            #origin = np.array([0., dy, 0.])
            #xzplane = origin + dx
            #xzplane = np.array([1., dy, 0.])
            #coord = CORD2R.add_axes(cid, rid=0, origin=p1, xaxis=p2-p1, yaxis=None, zaxis=None,
                                     #xyplane=None, yzplane=None, xzplane=None, comment='')
            print(coord)
        model.coords[1] = coord
        plane_bdf_filename = os.path.join(dirname, f'plane_face_{i:d}.bdf')
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
                plane_bdf_filename=plane_bdf_filename, plane_bdf_offset=dy)
        except RuntimeError:
            # incorrect ivalues=[0, 1, 2]; dy=771. for CRM
            continue
        unused_unique_geometry_array, unused_unique_results_array, rods = out

        if not os.path.exists(plane_bdf_filename):
            break
        plane_bdf_filenames.append(plane_bdf_filename)
        # eid, nid, inid1, inid2
        #print(unique_geometry_array)
        #moi_filename = 'amoi_%i.bdf' % i
        moi_filename = None
        out = calculate_area_moi(model, rods, normal_plane, thetas, moi_filename=moi_filename)
        #print(out)
        Ai, Ii, EIi, avg_centroidi = out
        y.append(dy)
        A.append(Ai)
        I.append(Ii)
        EI.append(EIi)
        avg_centroid.append(avg_centroidi)
        #break

    thetas_csv_filename = os.path.join(dirname, 'thetas.csv')

    with open(thetas_csv_filename, 'w') as csv_filename:
        csv_filename.write('# eid(%i),theta,Ex,Ey,Gxy\n')
        for eid, (theta, Ex, Ey, Gxy) in sorted(thetas.items()):
            csv_filename.write('%i,%f,%f,%f,%f\n' % (eid, theta, Ex, Ey, Gxy))

    y = np.array(y, dtype='float64')
    A = np.array(A, dtype='float64')
    I = np.array(I, dtype='float64')
    EI = np.array(EI, dtype='float64')
    avg_centroid = np.array(avg_centroid, dtype='float64')

    inid = 1
    beam_model = BDF(debug=False)
    avg_centroid[:, 1] = y

    # wrong
    mid = 1
    E = 3.0e7
    G = None
    nu = 0.3
    model.add_mat1(mid, E, G, nu, rho=0.1)

    Ix = I[:, 0]
    Iy = I[:, 1]
    Ixy = I[:, 2]
    J = Ix + Iy
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
        i2 = [Iy[eid-1], Iy[eid]]
        i12 = [Ixy[eid-1], Ixy[eid]]
        j = [J[eid-1], J[eid]]
        beam_model.add_pbeam(pid, mid, xxb, so, area, i1, i2, i12, j, nsm=None,
                             c1=None, c2=None, d1=None, d2=None, e1=None, e2=None, f1=None, f2=None,
                             k1=1., k2=1., s1=0., s2=0., nsia=0., nsib=None, cwa=0., cwb=None,
                             m1a=0., m2a=0., m1b=None, m2b=None,
                             n1a=0., n2a=0., n1b=None, n2b=None,
                             comment='')

    beam_model_bdf_filename = os.path.join(dirname, 'equivalent_beam_model.bdf')
    beam_model.write_bdf(beam_model_bdf_filename)

    X = np.vstack([y, A]).T
    Y = np.hstack([X, I, EI, avg_centroid])
    header = 'y, A, Ix, Iz, Ixz, Ex*Ix, Ex*Iz, Ex*Ixz, xcentroid, ycentroid, zcentroid'
    cut_data_span_filename = os.path.join(dirname, 'cut_data_vs_span.csv')
    np.savetxt(cut_data_span_filename, Y, header=header, delimiter=',')

    if IS_MATPLOTLIB and (plot or show):
        plot_inertia(y, A, I, EI, avg_centroid, show=show, dirname=dirname)
    else:
        plane_bdf_filenames = []
    return y, A, I, EI, avg_centroid, plane_bdf_filenames

def plot_inertia(y, A, I, EI, avg_centroid, ifig: int=1, show: bool=True, dirname: str=''):
    """hepler method for test"""
    #plt.plot(y, I[:, 0] / I[:, 0].max(), 'ro-', label='Qxx')
    #plt.plot(y, I[:, 1] / I[:, 1].max(), 'bo-', label='Qyy')
    #plt.plot(y, I[:, 2] / I[:, 2].max(), 'go-', label='Qxy')
    aI = np.abs(I)
    aEI = np.abs(EI)

    fig = plt.figure(ifig)
    ax = fig.gca()
    ax.plot(y, I[:, 0] / aI[:, 0].max(), 'ro-', label='Ixx')
    ax.plot(y, I[:, 1] / aI[:, 1].max(), 'bo-', label='Izz')
    ax.plot(y, I[:, 2] / aI[:, 2].max(), 'go-', label='Ixz')

    ax.plot(y, EI[:, 0] / aEI[:, 0].max(), 'ro-', label='EIxx', linestyle='--')
    ax.plot(y, EI[:, 1] / aEI[:, 1].max(), 'bo-', label='EIzz', linestyle='--')
    ax.plot(y, EI[:, 2] / aEI[:, 2].max(), 'go-', label='EIxz', linestyle='--')

    ax.grid(True)
    ax.set_xlabel('Span, y')
    ax.set_ylabel('Normalized Area MOI, I')
    ax.legend()
    fig.savefig('normalized_inertia_vs_span.png')
    #-------------------------------------------------------

    fig = plt.figure(ifig + 1)
    ax = fig.gca()
    ax.plot(y, A, 'ro-', label='Area', linestyle='-')

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
