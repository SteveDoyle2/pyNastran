"""
Calculate EI(y) and GJ(y)
"""
import os
import copy
from pathlib import Path
from typing import Any
from itertools import count

import numpy as np
try:
    import matplotlib.pyplot as plt  # pylint: disable=unused-import
    IS_MATPLOTLIB = True
except ModuleNotFoundError:  # pragma: no cover
    IS_MATPLOTLIB = False

from cpylog import SimpleLogger
from pyNastran.utils import PathLike
from pyNastran.bdf.bdf import BDF, read_bdf, CORD2R
from pyNastran.bdf.mesh_utils.cut_model_by_plane import (
    cut_face_model_by_coord,
    calculate_area_moi,
)


def cut_and_plot_moi(bdf_filename: PathLike | BDF,
                     normal_plane: np.ndarray,
                     log: SimpleLogger,
                     dys: list[float] | np.ndarray,
                     coords: list[CORD2R],
                     ytol: float=2.0,
                     face_data=None,
                     dirname: PathLike='',
                     debug_vectorize: bool=True,
                     plot: bool=True,
                     cut_data_span_filename: PathLike='cut_data_vs_span.csv',
                     beam_model_bdf_filename: PathLike='equivalent_beam_model.bdf',
                     show: bool=False) -> tuple[Any, Any, Any, Any, Any]: # y, A, I, EI, avg_centroid
    """
    For a shell structure, cut and plot
     - moments of inertia
     - stiffness

    Parameters
    ----------
    bdf_filename : PathLike
        the path to the bdf
    normal_plane : (3,) float ndarray
        the plane normal that defines ???
    log : SimpleLogger
        the logger
    dys : list[float] or ystations
        the y-stations to march down
    coords : list[CORD2R]
        coords to take cuts at; cutting plane normal is the y-axis?
        x:   defines axial direction (E1*A)
        y/z: defines transverse directions (E1*Iy)
    ytol : float; default=2.0
        ???
    face_data : ???
        nids : np.ndarray
            node ids
        xyz_cid0 : (nnode, 3) np.ndarray
            xyz values
        elements = dict[key, value]
            key = line, shell?
            value = (eids, ??)

    dirname : PathLike; default=''
        directory for output plots/csv/bdfs
    plot : bool; default=True
        not used
    show : bool; default=False
        show the plots at the end
    """
    if isinstance(dirname, str):
        dirname = Path(dirname)
    if isinstance(bdf_filename, PathLike):
        model = read_bdf(bdf_filename, log=log)
        model_static = copy.deepcopy(model)
    else:
        model = bdf_filename
        model_static = bdf_filename

    out = _get_station_data(
        model, model_static,
        dys, coords, normal_plane,
        ytol, dirname, face_data=face_data,
        debug_vectorize=debug_vectorize,
    )
    thetas, y, dx, dz, A, I, J, EI, GJ, avg_centroid, plane_bdf_filenames, plane_bdf_filenames2 = out

    assert len(y) > 0, y
    thetas_csv_filename = dirname / 'thetas.csv'

    with open(thetas_csv_filename, 'w') as csv_filename:
        csv_filename.write('# eid(%d),theta,Ex,Ey,Gxy\n')
        for eid, (theta, ex, ey, gxy) in sorted(thetas.items()):
            csv_filename.write(f'{eid:d},{theta},{ex},{ey},{gxy}\n')

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

    if beam_model_bdf_filename:
        beam_model_bdf_filename = dirname / beam_model_bdf_filename
        _write_beam_model(
            avg_centroid,
            A, Ix, Iz, Ixz,
            beam_model_bdf_filename,
        )

    if cut_data_span_filename:
        cut_data_span_filename = dirname / cut_data_span_filename
        X = np.vstack([y, dx, dz, A, Ix, Iz, Ixz, ExIx, ExIz, ExIxz]).T
        Y = np.hstack([X, avg_centroid])
        header = 'y, dx, dz, A, Ix, Iz, Ixz, Ex*Ix, Ex*Iz, Ex*Ixz, xcentroid, ycentroid, zcentroid'
        np.savetxt(cut_data_span_filename, Y, header=header, delimiter=',')

    plot_inertia(y, A, I, J, EI, GJ, avg_centroid, show=show, dirname=dirname)
    return y, A, I, J, EI, GJ, avg_centroid, plane_bdf_filenames, plane_bdf_filenames2


def _write_beam_model(avg_centroid: np.ndarray,
                      A, I, J, EI, GJ,
                      bdf_filename: PathLike=''):
    if isinstance(bdf_filename, str) and len(bdf_filename) == 0:
        return

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

    mid = 1
    beam_model = BDF(debug=False)
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
    beam_model.write_bdf(bdf_filename)


def _get_station_data(model: BDF,
                      model_static: BDF,
                      dys: list[float],
                      coords: list[CORD2R],
                      normal_plane: np.ndarray,
                      ytol: float, dirname: Path,
                      debug_vectorize: bool=True,
                      face_data=None) -> tuple[
                         dict[int, tuple[float, float, float, float]],  # thetas
                         #y, dx, dz,
                         #A, I, J,
                         #EI, GJ, avg_centroid
                         Any, Any, Any,
                         Any, Any, Any,
                         Any, Any, Any,
                         #plane_bdf_filenames, plane_bdf_filenames2,
                         list[str], list[str]]:
    """
    Helper for ``cut_and_plot_moi``

    Parameters
    ----------
    model : BDF()
        ???
    model_static : BDF()
        ???
    dys : list[float]
        the y values to make cuts at
    coords : list[CORD2R]
    normal_plane :
    ytol : float; default=2.0
    dirname :
    face_data : ???
        ???
    """
    log = model.log

    # initialize theta
    thetas = {}
    for eid in model.elements:
        #  theta, Ex, Ey, Gxy
        thetas[eid] = (0., 0., 0., 0.)

    #p1 = np.array([466.78845, 735.9053, 0.0])
    #p2 = np.array([624.91345, 639.68896, -0.99763656])
    #dx = p2 - p1
    nodal_result = None
    plane_bdf_filenames1 = []
    plane_bdf_filenames2 = []

    ny = len(dys)
    assert ny > 0, dys
    y = np.full(ny, np.nan, dtype='float64')
    dx = np.full(ny, np.nan, dtype='float64')
    dz = np.full(ny, np.nan, dtype='float64')
    A = np.full(ny, np.nan, dtype='float64')
    I = np.full((ny, 6), np.nan, dtype='float64')
    J = np.full(ny, np.nan, dtype='float64')
    EI = np.full((ny, 6), np.nan, dtype='float64')
    GJ = np.full(ny, np.nan, dtype='float64')
    avg_centroid = np.full((ny, 3), np.nan, dtype='float64')

    for icut, dy, coord in zip(count(), dys, coords):
        model.coords[1] = coord
        plane_bdf_filename1 = dirname / f'plane_face1_{icut:d}.bdf'
        plane_bdf_filename2 = dirname / f'plane_face2_{icut:d}.bdf'
        cut_face_filename = dirname / f'cut_face_{icut:d}.csv'
        if os.path.exists(cut_face_filename):
            os.remove(cut_face_filename)
        try:
            out = cut_face_model_by_coord(
                model_static, coord, ytol,
                nodal_result, plane_atol=1e-5,
                skip_cleanup=True,
                #csv_filename=cut_face_filename,
                csv_filename=None,
                #plane_bdf_filename=None)
                plane_bdf_filename1=plane_bdf_filename1,
                plane_bdf_filename2=plane_bdf_filename2,
                plane_bdf_offset=dy, face_data=face_data,
                debug_vectorize=debug_vectorize,
            )
        except PermissionError:
            print(f'failed to delete {plane_bdf_filename1}')
            raise
            # continue
        except RuntimeError:
            # incorrect ivalues=[0, 1, 2]; dy=771. for CRM
            raise
            # continue
        unused_unique_geometry_array, unused_unique_results_array, rods = out

        if not os.path.exists(plane_bdf_filename1):
            break
        plane_bdf_filenames1.append(plane_bdf_filename1)
        plane_bdf_filenames2.append(plane_bdf_filename2)
        # eid, nid, inid1, inid2
        #print(unique_geometry_array)
        #moi_filename = 'amoi_%i.bdf' % i
        moi_filename = None
        log.info(f'calculate_area_moi {icut:d}')
        dxi, dzi, Ai, Ii, EIi, avg_centroidi = calculate_area_moi(
            model, rods, normal_plane, thetas,
            moi_filename=moi_filename)

        #print(out)
        Ji = GJi = 1.0
        y[icut] = dy
        dx[icut] = dxi  # length
        dz[icut] = dzi  # height
        A[icut] = Ai
        I[icut, :] = Ii
        # print(Ji, EIi, GJi)
        # print(len(Ji), len(EIi), len(GJi))
        J[icut] = Ji
        EI[icut, :] = EIi
        GJ[icut] = GJi
        avg_centroid[icut, :] = avg_centroidi
        #break

    out = (
        thetas, y, dx, dz,
        A, I, J, EI, GJ,
        avg_centroid,
        plane_bdf_filenames1, plane_bdf_filenames2
    )
    return out


def plot_inertia(y, A, I, J, EI, GJ, avg_centroid,
                 ifig: int=1, show: bool=True,
                 dirname: PathLike='',
                 normalized_inertia_png_filename: PathLike='normalized_inertia_vs_span.png',
                 area_span_png_filename: PathLike='area_vs_span.png',
                 amoi_span_png_filename: PathLike = 'amoi_vs_span.png',
                 e_amoi_span_png_filename: PathLike = 'e_amoi_vs_span.png',
                 cg_span_png_filename: PathLike = 'cg_vs_span.png') -> None:
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
    png_filename = os.path.join(dirname, normalized_inertia_png_filename)
    fig.savefig(png_filename)
    #-------------------------------------------------------

    fig = plt.figure(ifig + 1)
    ax = fig.gca()
    ax.plot(y, A, 'ro', label='Area', linestyle='-')

    ax.grid(True)
    ax.set_xlabel('Span, y')
    ax.set_ylabel('Area, A')
    ax.legend()
    png_filename = os.path.join(dirname, area_span_png_filename)
    fig.savefig(png_filename)
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
    png_filename = os.path.join(dirname, amoi_span_png_filename)
    fig.savefig(png_filename)
    #-------------------------------------------------------


    fig = plt.figure(ifig + 3)
    ax = fig.gca()
    ax.plot(y, EI[:, 0], 'ro-', label='EIxx')
    #ax.plot(y, I[:, 0], 'bo-', label='Ixx')
    ax.grid(True)
    ax.set_xlabel('Span, y')
    ax.set_ylabel('Exx*Area MOI, Exx*I')
    ax.legend()
    png_filename = os.path.join(dirname, e_amoi_span_png_filename)
    fig.savefig(png_filename)
    #-------------------------------------------------------

    fig = plt.figure(ifig + 4)
    ax = fig.gca()
    ax.plot(y, avg_centroid[:, 0], 'ro-', label='xcg')
    ax.plot(y, avg_centroid[:, 2], 'bo-', label='zcg')
    ax.grid(True)
    ax.set_xlabel('Span, y')
    ax.set_ylabel('CG')
    ax.legend()
    png_filename = os.path.join(dirname, cg_span_png_filename)
    fig.savefig(png_filename)
    #-------------------------------------------------------

    if show:
        plt.show()
    ifig += 4
    return ifig
