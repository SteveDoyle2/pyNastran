"""
Calculate EI(y) and GJ(y)
"""
from __future__ import annotations
import os
import copy
from pathlib import Path
from typing import Any, TYPE_CHECKING
from itertools import count

import numpy as np
try:
    import matplotlib.pyplot as plt  # pylint: disable=unused-import
    IS_MATPLOTLIB = True
except ModuleNotFoundError:  # pragma: no cover
    IS_MATPLOTLIB = False

from cpylog import SimpleLogger
from pyNastran.utils import PathLike
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.cards.coordinate_systems import (
    CORD2R, Coord,
    xyz_to_rtz_array, rtz_to_xyz_array)
from pyNastran.bdf.bdf import BDF, read_bdf
from pyNastran.bdf.mesh_utils.cut.cut_model_by_plane import (
    cut_face_model_by_coord,
    fis_tri_cut, _setup_faces,
)
if TYPE_CHECKING:
    from pyNastran.bdf.cards.elements.shell import CTRIA3, CQUAD4
Rods = tuple[np.ndarray, np.ndarray, np.ndarray]


def cut_and_plot_moi(bdf_filename: PathLike | BDF,
                     normal_plane: np.ndarray,
                     log: SimpleLogger,
                     dys: list[float] | np.ndarray,
                     coords: list[CORD2R],
                     face_data=None,
                     dirname: PathLike='',
                     ifig: int=1,
                     debug_vectorize: bool=True,
                     stop_on_failure: bool=False,
                     cut_data_span_filename: PathLike='cut_data_vs_span.csv',
                     beam_model_bdf_filename: PathLike='equivalent_beam_model.bdf',
                     thetas_csv_filename: PathLike='thetas.csv',
                     normalized_inertia_png_filename: PathLike='normalized_inertia_vs_span.png',
                     area_span_png_filename: PathLike='area_vs_span.png',
                     amoi_span_png_filename: PathLike='amoi_vs_span.png',
                     e_amoi_span_png_filename: PathLike='e_amoi_vs_span.png',
                     cg_span_png_filename: PathLike='cg_vs_span.png',
                     plot: bool=True,
                     show: bool=False) -> tuple[Any, Any, Any, Any,     # y, A, I, J
                                                Any, Any, Any,          # EI, GJ, avg_centroid,
                                                list[str], list[str]]:  # plane_bdf_filenames1, plane_bdf_filenames2
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
        dirname, face_data=face_data,
        debug_vectorize=debug_vectorize,
        stop_on_failure=stop_on_failure,
    )
    (thetas, y, dx, dz, A, I, J, EI, GJ, avg_centroid,
     plane_bdf_filenames, plane_bdf_filenames2) = out

    assert len(y) > 0, y
    thetas_csv_filename = dirname / thetas_csv_filename

    with open(thetas_csv_filename, 'w') as csv_filename:
        csv_filename.write('# eid(%d),theta,Ex,Ey,Gxy\n')
        for eid, (theta, ex, ey, gxy) in sorted(thetas.items()):
            csv_filename.write(f'{eid:d},{theta},{ex},{ey},{gxy}\n')

    avg_centroid[:, 1] = y

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

        # wrong
        # model.add_mat1(mid=1, E=3.0e7, G=None, nu=0.3, rho=0.1)

        _write_beam_model(
            avg_centroid,
            A, I, J,
            EI, GJ,
            beam_model_bdf_filename,
        )

    if cut_data_span_filename:
        cut_data_span_filename = dirname / cut_data_span_filename
        X = np.vstack([y, dx, dz, A, Ix, Iz, Ixz, ExIx, ExIz, ExIxz]).T
        Y = np.hstack([X, avg_centroid])
        header = 'y, dx, dz, A, Ix, Iz, Ixz, Ex*Ix, Ex*Iz, Ex*Ixz, xcentroid, ycentroid, zcentroid'
        np.savetxt(cut_data_span_filename, Y, header=header, delimiter=',')

    # if plot:
    ifig = plot_inertia(
        y, A, I, J, EI, GJ, avg_centroid, show=show,
        dirname=dirname, ifig=ifig,
        normalized_inertia_png_filename=normalized_inertia_png_filename,
        area_span_png_filename=area_span_png_filename,
        amoi_span_png_filename=amoi_span_png_filename,
        e_amoi_span_png_filename=e_amoi_span_png_filename,
        cg_span_png_filename=cg_span_png_filename,
    )
    return y, A, I, J, EI, GJ, avg_centroid, plane_bdf_filenames, plane_bdf_filenames2


def _write_beam_model(avg_centroid: np.ndarray,
                      A: np.ndarray,
                      I: np.ndarray,
                      J: np.ndarray,
                      EI, GJ,
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
                      dirname: Path,
                      plane_atol: float=1e-5,
                      debug_vectorize: bool=True,
                      stop_on_failure: bool=False,
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
    dirname : Path | str
        base directory for output files/pictures
    face_data : ???
        ???
    """
    log = model.log

    # initialize theta
    thetas = {}
    for eid in model.elements:
        #  theta, Ex, Ey, Gxy
        thetas[eid] = (0., 0., 0., 0.)

    if face_data is None:
        # TODO: could filter out unused nodes
        _log, *face_data = _setup_faces(model)
    nodes, xyz_cid0, elements = face_data
    tri_eids, tri_nodes = elements['tri3']
    nnode = len(nodes)
    ntri = len(tri_eids)

    #p1 = np.array([466.78845, 735.9053, 0.0])
    #p2 = np.array([624.91345, 639.68896, -0.99763656])
    #dx = p2 - p1
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

    log.debug(f'dys={dys}; n={len(dys):d}')
    assert len(dys) == len(coords), (len(dys), len(coords))
    ncuts_found = 0
    for icut, dy, coord in zip(count(), dys, coords):
        itri_nodes = np.searchsorted(nodes, tri_nodes)
        xyz_cid = coord.transform_node_to_local_array(xyz_cid0)
        y_cid = xyz_cid[:, 1]
        # is_tri_cut = fis_tri_cut(y_cid, itri_nodes, ntri)

        model.coords[1] = coord
        plane_bdf_filename1 = dirname / f'plane_face1_{icut:d}.bdf'
        plane_bdf_filename2 = dirname / f'plane_face2_{icut:d}.bdf'
        cut_face_filename = dirname / f'cut_face_{icut:d}.csv'
        if os.path.exists(cut_face_filename):
            os.remove(cut_face_filename)

        found_cut, rods = _get_station_datai(
            model, model_static,
            dy, coord, plane_atol=plane_atol,
            debug_vectorize=debug_vectorize,
            stop_on_failure=stop_on_failure,
            plane_bdf_filename1=plane_bdf_filename1,
            plane_bdf_filename2=plane_bdf_filename2,
            face_data=face_data)

        # if not os.path.exists(plane_bdf_filename1) or len(rods) == 0:
        if not found_cut:
            log.debug(coord)
            log.debug(f'skipping calculate_area_moi {icut:d} (station={dy:g})')
            continue
            # break
        plane_bdf_filenames1.append(plane_bdf_filename1)
        plane_bdf_filenames2.append(plane_bdf_filename2)
        # eid, nid, inid1, inid2
        #print(unique_geometry_array)
        #moi_filename = 'amoi_%i.bdf' % i
        moi_filename = None
        log.info(f'calculate_area_moi {icut:d} (station={dy})')
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
        ncuts_found += 1
        #break
    if ncuts_found == 0:
        raise RuntimeError('no cuts found...')

    out = (
        thetas, y, dx, dz,
        A, I, J, EI, GJ,
        avg_centroid,
        plane_bdf_filenames1, plane_bdf_filenames2
    )
    return out


def _get_station_datai(model: BDF,
                       model_static: BDF,
                       dy: float,
                       coord: CORD2R,
                       plane_atol: float=1e-5,
                       debug_vectorize: bool=True,
                       stop_on_failure: bool=False,
                       plane_bdf_filename1: PathLike='',
                       plane_bdf_filename2: PathLike='',
                       face_data=None):
    nodal_result = None
    try:
        out = cut_face_model_by_coord(
            model_static, coord,
            nodal_result, plane_atol=plane_atol,
            skip_cleanup=True,
            # csv_filename=cut_face_filename,
            csv_filename=None,
            # plane_bdf_filename=None)
            plane_bdf_filename1=plane_bdf_filename1,
            plane_bdf_filename2=plane_bdf_filename2,
            plane_bdf_offset=dy, face_data=face_data,
            debug_vectorize=debug_vectorize,
            stop_on_failure=stop_on_failure,
        )
    except PermissionError:
        print(f'failed to delete {plane_bdf_filename1}')
        raise
        # continue
    except RuntimeError:
        # incorrect ivalues=[0, 1, 2]; dy=771. for CRM
        raise
        # continue
    found_cut, unused_unique_geometry_array, unused_unique_results_array, rods = out
    return found_cut, rods


def plot_inertia(y, A, I, J, EI, GJ, avg_centroid,
                 ifig: int=1, show: bool=True,
                 dirname: PathLike='',
                 normalized_inertia_png_filename: PathLike='normalized_inertia_vs_span.png',
                 area_span_png_filename: PathLike='area_vs_span.png',
                 amoi_span_png_filename: PathLike='amoi_vs_span.png',
                 e_amoi_span_png_filename: PathLike='e_amoi_vs_span.png',
                 cg_span_png_filename: PathLike='cg_vs_span.png') -> int:
    """helper method for test"""
    #plt.plot(y, I[:, 0] / I[:, 0].max(), 'ro-', label='Qxx')
    #plt.plot(y, I[:, 1] / I[:, 1].max(), 'bo-', label='Qyy')
    #plt.plot(y, I[:, 2] / I[:, 2].max(), 'go-', label='Qxy')
    aI = np.abs(I)
    aEI = np.abs(EI)
    aGJ = np.abs(GJ)

    fig = plt.figure(ifig)
    ax = fig.gca()
    ai_max = aI[:, :3].max(axis=0)
    aei_max = aEI[:, :3].max(axis=0)
    ai_max[ai_max == 0] = 1.
    aei_max[aei_max == 0] = 1.
    assert len(ai_max) == 3, (ai_max.shape, aI)
    ax.plot(y, I[:, 0] / ai_max[0], 'ro-', label='Ixx')
    ax.plot(y, I[:, 1] / ai_max[1], 'bo-', label='Izz')
    ax.plot(y, I[:, 2] / ai_max[2], 'go-', label='Ixz')

    ax.plot(y, EI[:, 0] / aei_max[0], 'ro', label='EIxx', linestyle='--')
    ax.plot(y, EI[:, 1] / aei_max[1], 'bo', label='EIzz', linestyle='--')
    ax.plot(y, EI[:, 2] / aei_max[2], 'go', label='EIxz', linestyle='--')
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


def calculate_area_moi(model: BDF,
                       rods: Rods,
                       normal_plane: np.ndarray,
                       thetas: dict[int, tuple[float, float, float, float]],
                       moi_filename: PathLike='',
                       eid_filename: PathLike='eid_file.csv',
                       ) -> tuple[Any, Any, Any, Any]:
    """
    The inertia of a square plate about the midplane is:
     Ixx = 1/12*b*h^3
     Iyy = 1/12*h*b^3
     Izz = 0.
     Ixy = Ixz = Iyz = 0.
    These terms are small for a real structure
    and the math gets harder for odd shapes,
    so we calculate just the A*d^2 terms.

    TODO: nevermind...this is just a 2d inertial formula
          of a flat plat that's been rotated

    Parameters
    ----------
    model : BDF
        the model object
    rods : (eids, nids, xyzs)
        eids : (nelements,) int ndarray
            the element id that was split
        nids : (nelements, 2) int ndarray
            the n1, n2 in xyzs that define the cut shell element
        xyzs : (nnodes, 3) float ndarray
            the xyz of the nodes
    normal_plane : (3,) float ndarray
        the direction of the cut plane
    thetas : dict[eid] = (thetad, Ex, Ey, Gxy)???
        thetas[eid] = (thetad, Ex, Ey, Gxy)
    moi_filename : str; default=None
        writes a csv file

    Returns
    -------
    total_area
    Isum
    Jsum
    EIsum
    GJsum
    avg_centroid
    """
    assert isinstance(rods, tuple), type(rods)
    assert isinstance(thetas, dict), type(thetas)
    rod_elements, rod_nids, rod_xyzs = rods
    assert isinstance(rod_elements, np.ndarray), type(rod_elements)
    assert isinstance(rod_nids, np.ndarray), type(rod_nids)
    assert isinstance(rod_xyzs, np.ndarray), type(rod_xyzs)

    eids = np.abs(rod_elements[:, 0])
    neids = len(eids)
    all_nids = rod_nids
    n1 = rod_elements[:, 1]
    n2 = rod_elements[:, 2]
    inid1 = np.searchsorted(all_nids, n1)
    inid2 = np.searchsorted(all_nids, n2)
    xyz1 = rod_xyzs[inid1, :]
    xyz2 = rod_xyzs[inid2, :]
    centroid = (xyz1 + xyz2) / 2.
    length = np.linalg.norm(xyz2 - xyz1, axis=1)
    assert len(length) == neids

    centroid, area, thickness, E = get_element_inertias(
        model, normal_plane, thetas,
        eids, length, centroid)

    # [Ixx, Iyy, Izz, Ixy, Iyz, Ixz]
    I: np.ndarray = np.zeros((len(area), 6), dtype='float64')

    # (Ex, Ey, Gxy)
    Ex = E[:, 0]

    total_area = area.sum()
    avg_centroid = (centroid * area[:, np.newaxis]) .sum(axis=0) / total_area
    assert len(avg_centroid) == 3, len(avg_centroid)
    # y corresponds to the station in the plane of the coordinate system
    # and is 0. because we're in the local plane
    x = centroid[:, 0] - avg_centroid[0]
    y = centroid[:, 1] - avg_centroid[1]
    z = centroid[:, 2] - avg_centroid[2]

    xmin = x.min()
    xmax = x.max()
    ixmin = np.where(x == xmin)[0][0]
    ixmax = np.where(x == xmax)[0][0]
    xyz_min = centroid[ixmin, :]
    xyz_max = centroid[ixmax, :]
    d = xyz_max - xyz_min
    dx = d[0]
    dz = d[2]
    theta = np.arctan2(dx, dz)
    thetad = np.degrees(theta)

    nnodes = len(x)
    delta = np.zeros((nnodes, 3))
    delta[:, 1] = thetad

    xyz = np.zeros((nnodes, 3))
    # we're swapping what axes we have to make the transform easier
    xyz[:, 0] = x
    xyz[:, 1] = z
    xyz[:, 2] = 0.
    rtz = xyz_to_rtz_array(xyz)
    rtz2 = rtz + delta
    xyz2 = rtz_to_xyz_array(rtz2)
    x2 = xyz2[:, 0]
    y2 = xyz2[:, 1]
    #z2 = xyz2[:, 2]

    #origin = d
    #zaxis = np.array([0., 1., 0.])
    #xzplane = d
    dxi = x2.max() - x2.min()
    dyi = y2.max() - y2.min()
    #dzi = z2.max() - z2.min()  # zero by definition

    I[:, 0] = area * (x * x)  # Ixx
    I[:, 1] = area * (y * y)  # Iyy
    I[:, 2] = area * (z * z)  # Izz
    I[:, 3] = area * (x * y)  # Ixy
    I[:, 4] = area * (y * z)  # Iyz
    I[:, 5] = area * (x * z)  # Ixz

    Isum = I.sum(axis=0)
    ExIsum = (Ex[:, np.newaxis] * I).sum(axis=0)
    assert len(Isum) == 6, len(Isum)

    if moi_filename is not None:
        dirname = os.path.dirname(moi_filename)
        eid_filename = os.path.join(dirname, eid_filename)
        _write_moi_file(
            moi_filename, eid_filename,
            eids, n1, n2, xyz1, xyz2, length, thickness, area,
            centroid, avg_centroid, I, E
        )
    return dxi, dyi, total_area, Isum, ExIsum, avg_centroid


def _write_moi_file(moi_filename: PathLike,
                    eid_filename: PathLike,
                    eids, n1, n2, xyz1, xyz2,
                    length, thickness, area,
                    centroid, avg_centroid, I, E) -> None:
    eidi = 1
    mid = 1
    nid0 = max(n1.max(), n2.max()) + 1
    with open(moi_filename, 'w') as bdf_file, open(eid_filename, 'w') as eid_file:
        bdf_file.write('$ pyNastran: punch=True\n')
        bdf_file.write('MAT1,1,3.0e7,,0.3\n')
        grid = ['GRID', nid0, 0, avg_centroid[0], avg_centroid[2], 0.]
        bdf_file.write(print_card_8(grid))
        bdf_file.write(f'CONM2   {1:8d}{nid0:8d}\n')

        fmt = ('%s,' * 7)[:-1] + '\n'
        eid_file.write('# eid(%i),pid(%i),area,thickness,Ixx,Izz,Ixz\n')
        for eid, n1i, n2i, xyz1i, xyz2i, lengthi, thicknessi, areai, centroidi, Ii, Ei in zip(
                eids, n1, n2, xyz1, xyz2, length, thickness, area, centroid, I, E):
            actual_eid = abs(eid)

            assert nid0 not in [n1i, n2i], (n1i, n2i)
            pidi = actual_eid
            #pid = eidi
            grid1 = ['GRID', n1i, None] + xyz1i.tolist()
            grid2 = ['GRID', n2i, None] + xyz2i.tolist()
            #crod = ['CROD', eidi, pid, n1i, n2i]
            A, J, nsm = Ii
            #prod = ['PROD', pid, mid, A, J, 0., nsm]
            assert eidi > 0, eidi
            conrod = ['CONROD', eidi, n1i, n2i, mid, A, J, 0., nsm]
            bdf_file.write(print_card_8(grid1))
            bdf_file.write(print_card_8(grid2))
            #bdf_file.write(print_card_8(crod))
            #bdf_file.write(print_card_8(prod))
            bdf_file.write(print_card_8(conrod))
            eidi += 1
            #PID | MID |  A  |  J  |  C  | NSM
            eid_file.write(fmt % (eidi, pidi, areai, thicknessi, Ii[0], Ii[1], Ii[2]))


def get_element_inertias(model: BDF,
                         normal_plane: np.ndarray,
                         thetas: dict[int, tuple[float, float, float, float]],
                         eids: list[int],
                         length: list[float],
                         centroid: list[np.ndarray],
                         ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    normal_plane_vector = normal_plane.copy().reshape((3, 1))
    cg_list: list[np.ndarray] = []
    area_list: list[float] = []
    thickness_list: list[float] = []
    E_list: list[tuple[float, float, float]] = []

    log = model.log
    for eid, lengthi, centroidi in zip(eids, length, centroid):
        #print(eid, lengthi)
        element = model.elements[eid]
        if element.type in ['CTRIA3', 'CQUAD4']:
            thicknessi, areai, thetad, Ex, Ey, Gxy, nu_xy = _get_shell_inertia(
                element, normal_plane, normal_plane_vector, lengthi)
            thetas[eid] = (thetad, Ex, Ey, Gxy)
            thickness_list.append(thicknessi)
            area_list.append(areai)
            cg_list.append(centroidi)
            E_list.append((Ex, Ey, Gxy))
        else:
            log.warning(element)

    centroid = np.array(cg_list, dtype='float64')
    area = np.array(area_list, dtype='float64')
    thickness = np.array(thickness_list, dtype='float64')
    E = np.array(E_list, dtype='float64')
    return centroid, area, thickness, E

def _get_shell_inertia(element: CTRIA3 | CQUAD4,
                       normal_plane: np.ndarray,
                       normal_plane_vector: np.ndarray,
                       lengthi: float,) -> tuple[float, float, float,
                                                 float, float, float, float]:
    """
    Parameters
    ----------
    element : CTRIA3 / CQUAD4
        the object to cut
    normal_plane : (3,) float ndarray
        the normal vector of the cutting plane (should be roughly normal to the element face)
    normal_plane_vector : (3,1) float ndarray
        the normal vector of the cutting plane (should be roughly normal to the element face)
    lengthi : float
        the length the cutting plane makes with the element

    Returns
    -------
    thicknessi : float
        the total thickness of the element
    areai : float
        the cut area of the element
    imat_rotation_angle_deg : float
        the angle between the cutting plane and the normal_plane / normal_plane_vector
        this is NOT the angle of the fiber
    Ex : float
        the moduli normal to the cut plane
    Ey : float
        the moduli parallel to the cut plane (normal to Ex)
    Gxy : float
        the inplane shear moduli
    nu_xy : float
        the correlary to in-plane nu12

    """
    pid_ref = element.pid_ref
    thicknessi = element.Thickness()
    dxyz, centroid, imat, unused_jmat, element_normal = element.material_coordinate_system()
    #print('imat = ', imat)
    #print('normal = ', normal)
    n1, n2, n3 = element_normal
    n12 = n1 * n2
    n13 = n1 * n3
    n23 = n2 * n3
    # http://scipp.ucsc.edu/~haber/ph216/rotation_12.pdf
    # expanding eq 20 into
    # R(n,theta) = R0 + R1*sin(theta) + R2*cos(theta)
    # R0 = np.array([
    #     [n1 ** 2, n12, n13],
    #     [n12, n2 ** 2, n23],
    #     [n13, n23, n3 ** 2],
    # ], dtype='float64')
    R1 = np.array([
        [0., -n3, n2],
        [n3, 0., -n1],
        [-n2, n1, 0.],
    ], dtype='float64')
    R2 = np.array([
        [1 - n1 ** 2, -n12, -n13],
        [-n12, 1 - n2 ** 2, -n23],
        [-n13, -n23, 1 - n3 ** 2],
    ])
    imat = imat.reshape(3, 1)
    #print(normal_plane.shape, R1.shape, imat.shape)
    #a = np.linalg.multi_dot([normal_plane.T, R0, imat])
    b = np.linalg.multi_dot([normal_plane_vector.T, R1, imat])
    c = np.linalg.multi_dot([normal_plane_vector.T, R2, imat])

    #  maximize m' dot p = p.T dot m
    # m' = R dot m
    #    = a + b*sin(theta) + c*cos(theta)
    #  d/d(theta) = b*cos(theta)*sin(theta) = 0
    #
    #  d/d(theta) = b*cos(theta) - c*sin(theta) = 0
    #  b*cos(theta) = c*sin(theta)
    #  tan(theta) = b/c
    #
    # the theta to rotate by in order to orient imat with the normal
    #print(b, c)
    imat_rotation_angle = np.arctan2(b, c).item()
    imat_rotation_angle_deg = np.degrees(imat_rotation_angle)
    if imat_rotation_angle_deg <= -90.:
        imat_rotation_angle_deg += 180.
    elif imat_rotation_angle_deg > 90.:
        imat_rotation_angle_deg -= 180.

    #element_normal = element.Normal()
    # cos(theta) = a o b / (|a| * |b|)
    # |a| = length of normal vector = 1.0
    # |b| = length of normal_plane vector = 1.0
    #
    # cos(theta) = a o b
    # then we take the absolute value because we don't care if the element is +/- theta off

    abs_cos_theta = abs(normal_plane @ element_normal)
    assert isinstance(imat_rotation_angle, float), imat_rotation_angle
    if abs_cos_theta > 0.9:  # <25.8 degrees
        # filter out elements that are in-plane
        thicknessi = 0.
        areai = 0.
        Ex = 0.
        Ey = 0.
        Gxy = 0.
        nu_xy = 0.
    else:
        Ex, Ey, Gxy, nu_xy = pid_ref.get_Ainv_equivalent_pshell(
            imat_rotation_angle_deg, thicknessi, # degrees=True,
        )

        #thicknessi = prop.Thickness()
        areai = thicknessi * lengthi

    # pid = pid_ref.pid
    # if pid == 10:
    #     import copy
    #     pid_ref45 = copy.deepcopy(pid_ref)
    #     pid_ref45.mids_ref = [copy.deepcopy(pid_ref.mids_ref[0])]
    #     pid_ref45.thetas = [copy.deepcopy(pid_ref.thetas[0])]
    #     pid_ref45.thicknesses = [copy.deepcopy(pid_ref.thicknesses[0])]
    #     pid_ref45.mids = [copy.deepcopy(pid_ref.mids[0])]
    #     pid_ref45.get_thetas()
    #     Ex45, Ey45, Gxy45, nu_xy45 = pid_ref45.get_Ainv_equivalent_pshell(
    #         imat_rotation_angle_deg, thicknessi)
    #
    #     pid_ref0 = copy.deepcopy(pid_ref)
    #     pid_ref0.mids_ref = [copy.deepcopy(pid_ref.mids_ref[1])]
    #     pid_ref0.thetas = [copy.deepcopy(pid_ref.thetas[1])]
    #     pid_ref0.thicknesses = [copy.deepcopy(pid_ref.thicknesses[1])]
    #     pid_ref0.mids = [copy.deepcopy(pid_ref.mids[1])]
    #     pid_ref0.get_thetas()
    #     Ex0, Ey0, Gxy0, nu_xy0 = pid_ref0.get_Ainv_equivalent_pshell(
    #         imat_rotation_angle_deg, thicknessi)
    return thicknessi, areai, imat_rotation_angle_deg, Ex, Ey, Gxy, nu_xy
