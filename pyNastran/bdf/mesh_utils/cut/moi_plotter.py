"""
Calculate EI(y) and GJ(y)
"""
from __future__ import annotations
import os
import copy
from pathlib import Path
from typing import Any, Optional, TYPE_CHECKING
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
    cut_face_model_by_coord, _setup_faces,
    # is_element_cut,
)
if TYPE_CHECKING:
    from pyNastran.bdf.cards.elements.shell import CTRIA3, CQUAD4
Rods = tuple[np.ndarray, np.ndarray, np.ndarray]


def cut_and_plot_moi(bdf_filename: PathLike | BDF,
                     normal_plane: np.ndarray,
                     log: SimpleLogger,
                     stations: list[float] | np.ndarray,
                     coords: list[CORD2R],
                     include_lines: bool=False,
                     include_solids: bool=False,
                     face_data: Optional[Any]=None,
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
                     show: bool=False) -> tuple[Any, Any, Any, Any,     # y, A, I, J,
                                                Any, Any, Any, Any,     # ExI, EyI, GJ, avg_centroid,
                                                list[str], list[str], int]:  # plane_bdf_filenames1, plane_bdf_filenames2, ifig
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
    stations : list[float] or ystations
        the y-stations to march down
    coords : list[CORD2R]
        coords to take cuts at; cutting plane normal is the y-axis?
        x:   defines axial direction (E1*A)
        y/z: defines transverse directions (E1*Iy)
    face_data : tuple
        nids : np.ndarray
            node ids
        xyz_cid0 : (nnode, 3) np.ndarray
            xyz values
        elements = dict[key, value]
            key: str
                'line', 'tri3', 'tet4'
            value: tuple[eids, node_ids]
                The element ids and associated node ids
                Tri3s have negative element ids, which correspond to the split CQUAD4s?
    include_lines: bool=False
        unused
    include_solids: bool=False
        unused
    dirname : PathLike; default=''
        directory for output plots/csv/bdfs
    plot : bool; default=True
        not used
    show : bool; default=False
        show the plots at the end
    ifig : int; default=1
        lets you change the figure ID, useful when you do multiple cuts
    stop_on_failure : bool; default=False
        useful for debugging or things you know should be cut
    debug_vectorize : bool; default=True
        test faster cutting method
     cut_data_span_filename : PathLike; default='cut_data_vs_span.csv'
        changes the filename
     beam_model_bdf_filename : PathLike; default='equivalent_beam_model.bdf'
        changes the filename
     thetas_csv_filename : PathLike; default='thetas.csv'
        changes the filename
     normalized_inertia_png_filename : PathLike; default='normalized_inertia_vs_span.png'
        changes the filename
     area_span_png_filename : PathLike; default='area_vs_span.png'
        changes the filename
     amoi_span_png_filename : PathLike; default='amoi_vs_span.png'
        changes the filename
     e_amoi_span_png_filename : PathLike; default='e_amoi_vs_span.png'
        changes the filename
     cg_span_png_filename : PathLike; default='cg_vs_span.png'
        changes the filename

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
        stations, coords, normal_plane,
        dirname, face_data=face_data,
        include_lines=include_lines, include_solids=include_solids,
        debug_vectorize=debug_vectorize,
        stop_on_failure=stop_on_failure,
    )
    (thetas, stations, dx, dz, A, I, J, ExI, EyI, GJ, avg_centroid,
     plane_bdf_filenames, plane_bdf_filenames2) = out

    assert len(stations) > 0, stations
    thetas_csv_filename = dirname / thetas_csv_filename
    with open(thetas_csv_filename, 'w') as csv_filename:
        csv_filename.write('# eid(%d),theta,Ex,Ey,Gxy\n')
        for eid, (theta, ex, ey, gxy) in sorted(thetas.items()):
            csv_filename.write(f'{eid:d},{theta},{ex},{ey},{gxy}\n')

    avg_centroid[:, 1] = stations

    #   0    1    2    3    4    5
    # [Ixx, Iyy, Izz, Ixy, Iyz, Ixz]
    Ix = I[:, 0]
    # Iy = I[:, 1]
    Iz = I[:, 2]
    Ixz = I[:, 5]

    ExIx = ExI[:, 0]
    # ExIy = ExI[:, 1]
    ExIz = ExI[:, 2]
    ExIxz = ExI[:, 5]
    # Ex = ExIx / Ix
    # Ey = ExIz / Iz
    J = Ix + Iz
    G = GJ / J
    #i1, i2, i12 = Ix, Iy, Ixy

    if beam_model_bdf_filename:
        beam_model_bdf_filename = dirname / beam_model_bdf_filename
        # wrong
        # model.add_mat1(mid=1, E=3.0e7, G=None, nu=0.3, rho=0.1)
        _write_beam_model(
            avg_centroid,
            A, I, J,
            ExI, GJ,
            beam_model_bdf_filename,
        )

    if cut_data_span_filename:
        cut_data_span_filename = dirname / cut_data_span_filename
        inotnan = np.isfinite(stations)
        X = np.column_stack([stations, dx, dz, A,
                             I, J, #Ex, Ey, G,
                             ExI,
                             EyI,
                             GJ,
                             avg_centroid])[inotnan, :]

        header = (
            'station,dx,dz,A,'
            'Ix,Iy,Iz,Ixy,Ixz,Iyz,J,'
            'Ex*Ix,Ex*Iy,Ex*Iz,Ex*Ixy,Ex*Ixz,Ex*Iyz,'
            'Ey*Ix,Ey*Iy,Ey*Iz,Ey*Ixy,Ey*Ixz,Ey*Iyz,'
            'GJ,'
            'xcentroid,ycentroid,zcentroid')
        np.savetxt(cut_data_span_filename, X, header=header, delimiter=',')

    if plot:
        ifig = plot_inertia(
            log, stations, A, I, J,
            ExI, EyI, GJ, avg_centroid, show=show,
            dirname=dirname, ifig=ifig,
            normalized_inertia_png_filename=normalized_inertia_png_filename,
            area_span_png_filename=area_span_png_filename,
            amoi_span_png_filename=amoi_span_png_filename,
            e_amoi_span_png_filename=e_amoi_span_png_filename,
            cg_span_png_filename=cg_span_png_filename,
        )

    return stations, A, I, J, ExI, EyI, GJ, avg_centroid, plane_bdf_filenames, plane_bdf_filenames2, ifig


def load_moi_data(csv_filename: PathLike) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray,
                                                   np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    import pandas as pd
    # 'station, dx, dz, A, '
    # 'Ix, Iy, Iz, Ixy, Ixz, Iyz, J, '
    # 'Ex*Ix, Ex*Iy, Ex*Iz, Ex*Ixy, Ex*Ixz, Ix*Iyz, '
    # 'Ey*Ix, Ey*Iy, Ey*Iz, Ey*Ixy, Ey*Ixz, Iy*Iyz, '
    # 'GJ,'
    # 'xcentroid, ycentroid, zcentroid')
    df = pd.read_csv(csv_filename)
    df.columns = df.columns.str.strip(' #')
    # print(df.columns)
    y = df['station'].to_numpy()
    A = df['A'].to_numpy()
    I = df[['Ix', 'Iy', 'Iz', 'Ixy', 'Iyz', 'Ixz']].to_numpy()
    ExI = df[['Ex*Ix', 'Ex*Iy', 'Ex*Iz', 'Ex*Ixy', 'Ex*Iyz', 'Ex*Ixz']].to_numpy()
    EyI = df[['Ey*Ix', 'Ey*Iy', 'Ey*Iz', 'Ey*Ixy', 'Ey*Iyz', 'Ey*Ixz']].to_numpy()
    GJ = df['GJ'].to_numpy()
    J = df['J'].to_numpy()
    avg_centroid = df[['xcentroid', 'ycentroid', 'zcentroid']].to_numpy()
    return y, A, I, J, ExI, EyI, GJ, avg_centroid


def _write_beam_model(avg_centroid: np.ndarray,
                      A: np.ndarray,
                      I: np.ndarray,
                      J: np.ndarray,
                      EI: np.ndarray, GJ: np.ndarray,
                      bdf_filename: PathLike=''):
    if isinstance(bdf_filename, str) and len(bdf_filename) == 0:
        return

    #   0    1    2    3    4    5
    # [Ixx, Iyy, Izz, Ixy, Iyz, Ixz]
    Ix = I[:, 0]
    # Iy = I[:, 1]
    Iz = I[:, 2]
    Ixz = I[:, 5]

    # ExIx = EI[:, 0]
    # ExIy = EI[:, 1]
    # ExIz = EI[:, 2]
    # ExIxz = EI[:, 5]

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
                      include_lines: bool=False,
                      include_solids: bool=False,
                      debug_vectorize: bool=True,
                      stop_on_failure: bool=False,
                      face_data=None) -> tuple[
                         dict[int, tuple[float, float, float, float]],  # thetas
                         #y, dx, dz,
                         Any, Any, Any,
                         #A, I, J,
                         Any, Any, Any,
                         #ExI, EyI, GJ,
                         Any, Any, Any,
                         #avg_centroid, plane_bdf_filenames, plane_bdf_filenames2,
                         Any, list[str], list[str]]:
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
    normal_plane : np.ndarray
    dirname : Path | str
        base directory for output files/pictures
    face_data : ???
        ???
    stop_on_failure : bool; default=False
        useful for debugging or things you know should be cut
    """
    log = model.log

    # initialize theta
    thetas = {}
    for eid in model.elements:
        #  theta, Ex, Ey, Gxy
        thetas[eid] = (0., 0., 0., 0.)

    if face_data is None:
        # TODO: could filter out unused nodes
        _log, *face_data = _setup_faces(
            model,
            include_lines=include_lines, include_solids=include_solids)
    nodes, xyz_cid0, elements = face_data
    tri_eids, tri_nodes, _zoffset = elements['tri3']
    # nnode = len(nodes)
    # ntri = len(tri_eids)

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
    area = np.full(ny, np.nan, dtype='float64')
    inertia = np.full((ny, 6), np.nan, dtype='float64')
    J = np.full(ny, np.nan, dtype='float64')
    ExI = np.full((ny, 6), np.nan, dtype='float64')
    EyI = np.full((ny, 6), np.nan, dtype='float64')
    GJ = np.full(ny, np.nan, dtype='float64')
    avg_centroid = np.full((ny, 3), np.nan, dtype='float64')

    log.debug(f'dys={dys}; n={len(dys):d}')
    assert len(dys) == len(coords), (len(dys), len(coords))
    ncuts_found = 0
    for icut, dy, coord in zip(count(), dys, coords):
        # itri_nodes = np.searchsorted(nodes, tri_nodes)
        # xyz_cid = coord.transform_node_to_local_array(xyz_cid0)
        # y_cid = xyz_cid[:, 1]
        # is_tri_cut = fis_tri_cut(y_cid, itri_nodes, ntri)

        model.coords[coord.cid] = coord
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
        (dxi, dzi, areai,
         inertiai, Ji,
         ExIi, EyIi, GJi, avg_centroidi) = calculate_area_moi(
            model, rods, normal_plane, thetas,
            moi_filename=moi_filename)

        #print(out)
        y[icut] = dy
        dx[icut] = dxi  # length
        dz[icut] = dzi  # height
        area[icut] = areai
        inertia[icut, :] = inertiai
        # print(Ji, EIi, GJi)
        # print(len(Ji), len(EIi), len(GJi))
        J[icut] = Ji
        ExI[icut, :] = ExIi
        EyI[icut, :] = EyIi
        GJ[icut] = GJi
        avg_centroid[icut, :] = avg_centroidi
        ncuts_found += 1
        #break
    if ncuts_found == 0:
        raise RuntimeError('no cuts found...')

    out = (
        thetas, y, dx, dz,
        area, inertia, J,
        ExI, EyI, GJ,
        avg_centroid, plane_bdf_filenames1, plane_bdf_filenames2
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


def plot_inertia(log: SimpleLogger,
                 y: np.ndarray, A: np.ndarray,
                 I: np.ndarray, J: np.ndarray,
                 ExI: np.ndarray, EyI: np.ndarray, GJ: np.ndarray,
                 avg_centroid: np.ndarray,
                 ifig: int=1, show: bool=True,
                 dirname: PathLike='',
                 normalized_inertia_png_filename: PathLike='normalized_inertia_vs_span.png',
                 area_span_png_filename: PathLike='area_vs_span.png',
                 amoi_span_png_filename: PathLike='amoi_vs_span.png',
                 e_amoi_span_png_filename: PathLike='e_amoi_vs_span.png',
                 cg_span_png_filename: PathLike='cg_vs_span.png') -> int:
    """helper method for test

    inertia: [Ixx, Iyy, Izz, Ixy, Iyz, Ixz]
    """
    absI = np.abs(I)
    absExI = np.abs(ExI)
    # absGJ = np.abs(GJ)

    assert isinstance(ifig, int), ifig
    fig = plt.figure(ifig)
    ax = fig.gca()
    ai_max = absI[:, :3].max(axis=0)
    aei_max = absExI[:, :3].max(axis=0)
    ai_max[ai_max == 0] = 1.
    aei_max[aei_max == 0] = 1.
    assert len(ai_max) == 3, (ai_max.shape, absI)
    ax.plot(y, I[:, 0] / ai_max[0], 'ro-', label='Ixx')
    ax.plot(y, I[:, 1] / ai_max[1], 'bo-', label='Izz')
    ax.plot(y, I[:, 2] / ai_max[2], 'go-', label='Ixz')

    ax.plot(y, ExI[:, 0] / aei_max[0], 'ro', label='ExIxx', linestyle='--')
    ax.plot(y, ExI[:, 1] / aei_max[1], 'bo', label='ExIzz', linestyle='--')
    ax.plot(y, ExI[:, 2] / aei_max[2], 'go', label='ExIxz', linestyle='--')
    #ax.plot(y, GJ / aGJ.max(), 'go-', label='GJ', linestyle='--')

    ax.grid(True)
    ax.set_xlabel('Span, y')
    ax.set_ylabel('Normalized Area MOI, I')
    ax.legend()
    png_filename = os.path.join(dirname, normalized_inertia_png_filename)
    log.info(f'saving {png_filename}')
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
    log.debug(f'saving {png_filename}')
    fig.savefig(png_filename)
    #-------------------------------------------------------

    fig = plt.figure(ifig + 2)
    ax = fig.gca()
    ax.plot(y, I[:, 0], 'ro-', label='Ixx')
    ax.plot(y, I[:, 1], 'bo-', label='Iyy')
    ax.plot(y, I[:, 2], 'go-', label='Izz')
    ax.grid(True)
    ax.set_xlabel('Span, y')
    ax.set_ylabel('Area MOI, I')
    ax.legend()
    png_filename = os.path.join(dirname, amoi_span_png_filename)
    log.debug(f'saving {png_filename}')
    fig.savefig(png_filename)
    #-------------------------------------------------------


    fig = plt.figure(ifig + 3)
    ax = fig.gca()
    ax.plot(y, ExI[:, 0], 'ro-', label='EIxx')
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
    log.debug(f'saving {png_filename}')
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
                       ) -> tuple[np.ndarray, np.ndarray, np.ndarray,               # dxi, dyi, total_area,
                                  np.ndarray, np.ndarray,                           # Isum, Jsum,
                                  np.ndarray, np.ndarray, np.ndarray, np.ndarray]:  # ExIsum, EyIsum, GJsum, avg_centroid
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
    inertia: np.ndarray = np.zeros((len(area), 6), dtype='float64')

    # (Ex, Ey, Gxy)
    ex = E[:, 0]
    ey = E[:, 1]
    gxy = E[:, 2]

    total_area = area.sum()
    if total_area == 0.0:
        avg_centroid = centroid.mean(axis=0)
    else:
        avg_centroid = (centroid * area[:, np.newaxis]) .sum(axis=0) / total_area
    assert len(avg_centroid) == 3, len(avg_centroid)
    # y corresponds to the station in the plane of the coordinate system
    # and is 0. because we're in the local plane
    x = centroid[:, 0] - avg_centroid[0]
    y = centroid[:, 1] - avg_centroid[1]
    z = centroid[:, 2] - avg_centroid[2]

    xmin = x.min()
    xmax = x.max()
    ixmins = np.where(x == xmin)[0]
    ixmaxs = np.where(x == xmax)[0]
    try:
        ixmin = ixmins[0]
        ixmax = ixmaxs[0]
    except IndexError:
        print(f'x = {x}')
        print(f'ixmins = {ixmins}')
        print(f'ixmaxs = {ixmaxs}')
        raise RuntimeError('bad cut?')

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

    inertia[:, 0] = area * (x * x)  # Ixx
    inertia[:, 1] = area * (y * y)  # Iyy
    inertia[:, 2] = area * (z * z)  # Izz
    inertia[:, 3] = area * (x * y)  # Ixy
    inertia[:, 4] = area * (y * z)  # Iyz
    inertia[:, 5] = area * (x * z)  # Ixz

    # cut is in xz plane
    ix = inertia[:, 0]
    iz = inertia[:, 2]
    J = ix + iz

    Isum = inertia.sum(axis=0)
    Jsum = J.sum()
    ExIsum = (ex[:, np.newaxis] * inertia).sum(axis=0)
    EyIsum = (ey[:, np.newaxis] * inertia).sum(axis=0)
    GJsum = (gxy * J).sum()
    assert len(Isum) == 6, len(Isum)

    if moi_filename is not None:
        dirname = os.path.dirname(moi_filename)
        eid_filename = os.path.join(dirname, eid_filename)
        _write_moi_file(
            moi_filename, eid_filename,
            eids, n1, n2, xyz1, xyz2, length, thickness, area,
            centroid, avg_centroid, inertia, E,
        )
    out = (
        dxi, dyi, total_area,
        Isum, Jsum,
        ExIsum, EyIsum, GJsum, avg_centroid,
    )
    return out


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


def plot_compare_inertia(log: SimpleLogger,
                         csv_filenames: list[tuple[Path, str, str]],
                         # ifig: int=1,
                         x: str='x',
                         y: str='y',
                         z: str='z',
                         yrange=None,
                         span_label = 'Span, y (in)',
                         dirname='', save: bool=True,
                         show: bool=True) -> int:
    """helper method for test

    inertia: [Ixx, Iyy, Izz, Ixy, Iyz, Ixz]
    """
    xx = f'{x}{x}'
    # xy = f'{x}{y}'
    xz = f'{x}{z}'
    # yz = f'{y}{z}'
    zz = f'{z}{z}'
    yy = f'{y}{y}'

    marker = ''
    # absI = np.abs(I)
    # absExI = np.abs(ExI)
    # absGJ = np.abs(GJ)

    save0 = save
    save = False
    ilast_file = len(csv_filenames) - 1
    data = []
    for ifile, (cut_data_span_filename, tag, linestyle) in enumerate(csv_filenames):
        is_last_file = (ifile == ilast_file)
        if is_last_file:
            # last file
            save = save0

        log.info(f'loading {str(cut_data_span_filename)}')
        datai = load_moi_data(cut_data_span_filename)
        data.append(datai)
        station, A, I, J, ExI, EyI, GJ, avg_centroid = datai
        # assert station.max() < 2000., station.max()
        Ixx = I[:, 0]
        Iyy = I[:, 1]
        Izz = I[:, 2]
        Ixy = I[:, 3]
        Iyz = I[:, 4]
        Ixz = I[:, 5]
        Ex = ExI[:, 0]/I[:,0]
        Ey = EyI[:, 0]/I[:,0]
        G = GJ/J
        ExA =  Ex * A

        # assert isinstance(ifig, int), ifig
        #-------------------------------------------------------
        ifig = 1
        fig = plt.figure(ifig)
        ax = fig.gca()
        ax.plot(station, A, 'r', marker=marker, label=f'{tag}Area', linestyle=linestyle)
        ax.grid(True)
        ax.set_xlabel(span_label)
        ax.set_ylabel('Area, A ($in^2$)')
        ax.legend()

        if save:
            png_filename = dirname / 'area_vs_span.png'
            log.debug(f'saving {png_filename}')
            fig.savefig(png_filename)
        ifig += 1
        #-------------------------------------------------------
        fig = plt.figure(ifig)
        ax = fig.gca()
        ax.semilogy(station, Ixx, 'C0', linestyle=linestyle, marker=marker, label=f'{tag}I{xx}')
        # ax.semilogy(station, Iyy, 'C1', linestyle=linestyle, marker=marker, label='Iyy')
        ax.semilogy(station, Izz, 'C2', linestyle=linestyle, marker=marker, label=f'{tag}I{zz}')
        # ax.semilogy(station, Ixy, 'C3', linestyle=linestyle, marker=marker, label='Ixy')
        # ax.semilogy(station, Iyz, 'C4', linestyle=linestyle, marker=marker, label='Iyz')
        ax.semilogy(station, np.abs(Ixz), 'C5', linestyle=linestyle, marker=marker, label=f'{tag}I{xz}')
        ax.semilogy(station, J, 'k', linestyle=linestyle, marker=marker, label=f'{tag}J')
        ax.grid(True)
        ax.set_xlabel(span_label)
        ax.set_ylabel('Area MOI, I ($in^4$)')
        ax.legend()
        png_filename = dirname / 'Ixx_Izz_Ixz_J_log.png'
        log.debug(f'saving {png_filename}')
        fig.savefig(png_filename)
        ifig += 1

        fig = plt.figure(ifig)
        ax = fig.gca()
        ax.plot(station, Ixx, 'C0', linestyle=linestyle, marker=marker, label=f'{tag}I{xx}')
        # ax.plot(station, Iyy, 'C1', linestyle=linestyle, marker=marker, label='Iyy')
        ax.plot(station, Izz, 'C2', linestyle=linestyle, marker=marker, label=f'{tag}I{zz}')
        # ax.plot(station, Ixy, 'C3', linestyle=linestyle, marker=marker, label='Ixy')
        # ax.plot(station, Iyz, 'C4', linestyle=linestyle, marker=marker, label='Iyz')
        ax.plot(station, Ixz, 'C5', linestyle=linestyle, marker=marker, label=f'{tag}I{xz}')
        ax.plot(station, J, 'k', linestyle=linestyle, marker=marker, label=f'{tag}J')
        if save:
            ax.grid(True)
            ax.set_xlabel(span_label)
            ax.set_ylabel('Area MOI, I ($in^4$)')
            ax.legend()
            png_filename = dirname / 'Ixx_Izz_Ixz_J.png'
            log.debug(f'saving {png_filename}')
            fig.savefig(png_filename)
        ifig += 1

        #-------------------------------------------------------
        ifig_EyIzz = ifig
        fig = plt.figure(ifig)
        ax = fig.gca()
        ax.plot(station, ExI[:, 2], color='r', linestyle=linestyle, marker=marker, label=f'{tag}E{y}*I{zz}')  # Ey*Izz
        #ax.plot(station, I[:, 0], 'bo-', label='Ixx')
        if save:
            ax.grid(True)
            ax.set_xlabel(span_label)
            ax.set_ylabel(f'Stiffness: E{y}*Area')
            ax.legend()
            png_filename = dirname / f'stiffness_E{y}I{zz}.png'
            fig.savefig(png_filename)
        ifig += 1

        #---------------------------------------------------
        ifig_GJ = ifig
        fig = plt.figure(ifig)
        ax = fig.gca()
        ax.plot(station, GJ, color='k', linestyle=linestyle, marker=marker, label=f'{tag}G{xz}*J')
        #ax.plot(station, I[:, 0], 'b-', marker=marker, label='Ixx')
        if save:
            ax.grid(True)
            ax.set_xlabel(span_label)
            ax.set_ylabel(f'Stiffness: G{xz}*J')
            ax.legend()
            png_filename = dirname / 'stiffness_GJ.png'
            fig.savefig(png_filename)
        ifig += 1

        #---------------------------------------------------
        fig = plt.figure(ifig)
        ax = fig.gca()
        ax.plot(station, GJ, 'k', linestyle=linestyle, marker=marker, label=f'{tag}G{xz}*J')
        ax.plot(station, ExI[:, 0], 'r', linestyle=linestyle, marker=marker, label=f'{tag}E{y}*I{xx}')
        ax.plot(station, ExA, 'b', linestyle=linestyle, marker=marker, label=f'{tag}E{y}*A')
        if save:
            ax.grid(True)
            ax.set_xlabel(span_label)
            ax.set_ylabel(f'Stiffness: G{xz}*J')
            ax.legend()
            png_filename = dirname / 'stiffness_ExIxx_GJ.png'
            fig.savefig(png_filename)
        ifig += 1
        #---------------------------------------------------
        fig = plt.figure(ifig)
        ax = fig.gca()
        ax.plot(station, G, color='k', linestyle=linestyle, marker=marker, label=f'G{xz}')
        ax.plot(station, Ex, color='r', linestyle=linestyle, marker=marker, label=f'E{y}')
        ax.plot(station, Ey, color='b', linestyle=linestyle, marker=marker, label=f'E{x}')
        if save:
            ax.grid(True)
            ax.set_xlabel(span_label)
            ax.set_ylabel(f'Effective Modulus: E{x}, E{z}, G{xz}')
            ax.legend()
            png_filename = dirname / 'stiffness_Ex_Ey_G.png'
            fig.savefig(png_filename)
        ifig += 1

        #-------------------------------------------------------
        fig = plt.figure(ifig)
        ax = fig.gca()
        ax.plot(station, avg_centroid[:, 0], color='r', linestyle=linestyle, marker=marker, label=f'{tag}xcg')
        ax.plot(station, avg_centroid[:, 2], color='b', linestyle=linestyle, marker=marker, label=f'{tag}zcg')
        if save:
            ax.grid(True)
            ax.set_xlabel(span_label)
            ax.set_ylabel('CG (in)')
            ax.legend()
            png_filename = dirname / 'cg_vs_span.png'
            log.debug(f'saving {png_filename}')
            fig.savefig(png_filename)
        #-------------------------------------------------------
        ifig += 4

    #----------------------------------------------------------
    ifig = 1
    if len(data) == 2:
        tag1 = csv_filenames[0][1]
        tag2 = csv_filenames[1][1]
        station1, A1, I1, J1, ExI1, EyI1, GJ1, avg_centroid1 = data[0]
        station2, A2, I2, J2, ExI2, EyI2, GJ2, avg_centroid2 = data[1]
        ustation = np.unique(np.hstack([station1, station2]))
        common_station = np.intersect1d(station1, station2)
        istation1 = np.searchsorted(station1, common_station)
        istation2 = np.searchsorted(station2, common_station)

        # assert np.allclose(station1, station2)
        # Ixx1 = I1[:, 0]
        # Izz1 = I1[:, 2]
        # Ixz1 = I1[:, 5]
        Ex1 = ExI1[:, 0]/I1[:,0]
        # Ey1 = EyI1[:, 0]/I1[:,0]
        # G1 = GJ1/J1
        ExA1 =  Ex1 * A1

        # Ixx2 = I2[:, 0]
        # Izz2 = I2[:, 2]
        # Ixz2 = I2[:, 5]
        Ex2 = ExI2[:, 0]/I2[:,0]
        # Ey2 = EyI2[:, 0]/I2[:,0]
        # G = GJ2/J2
        ExA2 =  Ex2 * A2

        # fig = plt.figure(ifig)
        # ax = fig.gca()

        ax_yscale = 'linear'
        ifig = 10
        fig = plt.figure(ifig)
        ax = fig.gca()
        ax.set_yscale(ax_yscale)
        ax2 = ax.twinx()
        # 1 / 2 will give us less
        ExA_ratio = ExA1[istation1]/ExA2[istation2] -1
        ax.plot(station1, ExA1, color='r', linestyle='-', marker=marker, label=f'{tag1}E{y}*A')  # Ey*Izz
        ax.plot(station2, ExA2, color='b', linestyle='--', marker=marker, label=f'{tag2}E{y}*A')  # Ey*Izz
        ax2.plot(common_station, ExA_ratio*100, color='k', linestyle='-', marker=marker, label=f'ratio')  # Ey*Izz
        ax.grid()
        ax.legend()
        ax.set_xlabel(span_label)
        ax.set_ylabel(f'Stiffness ($lb_f*in^2$): E{y}*I{zz}')
        ax2.set_ylabel('%Difference')

        ifig += 1
        fig = plt.figure(ifig)
        ax = fig.gca()
        ax.set_yscale(ax_yscale)
        ax2 = ax.twinx()
        # 1 / 2 will give us less
        ax.grid()
        ExI_ratio = ExI1[istation1, 2]/ExI2[istation2, 2] - 1
        ax.plot(station1, ExI1[:, 2], color='r', linestyle='-', marker=marker, label=f'{tag1}E{y}I{zz}')  # Ey*Izz
        ax.plot(station2, ExI2[:, 2], color='b', linestyle='--', marker=marker, label=f'{tag2}E{y}I{zz}')  # Ey*Izz
        ax2.plot(common_station, ExI_ratio*100, color='k', linestyle='-', marker=marker, label=f'ratio')  # Ey*Izz
        ax.set_ylabel(f'Stiffness ($lb_f*in^2$): E{y}*I{zz}')
        ax2.set_ylabelf'%Difference')
        ax.legend()
        ifig += 1

        fig = plt.figure(ifig)
        ax = fig.gca()
        ax.set_yscale(ax_yscale)
        ax2 = ax.twinx()
        GJ_ratio = GJ1[istation1]/GJ2[istation2]-1
        ax.plot(station1, GJ1, color='r', linestyle='-', marker=marker, label=f'{tag1}GJ')
        ax.plot(station2, GJ2, color='b', linestyle='--', marker=marker, label=f'{tag2}GJ')
        ax2.plot(common_station, GJ_ratio*100, color='k', linestyle='-', marker=marker, label=f'ratio')
        ax.set_ylabel('Stiffness ($lb_f*in^2$): G*J')
        ax2.set_ylabel('%Difference')
        ax.grid()
        ax.legend()
        plt.show()
    ifig += 1

    if show:
        plt.show()
    return ifig
