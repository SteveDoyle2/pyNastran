"""
defines:
 - nids, nid_cd, icd_transform, xyz_cid0 = get_nid_cd_xyz_cid0(model)
 - nids, nid_cd, xyz_cid0, icd_transform, eids, element_centroids_cid0 = smt_setup(model)
 - xyz1, xyz2, xyz3, i, k, origin, xzplane, dim_max, stations = setup_coord_from_plane(
        model, xyz_cid0,
        p1, p2, p3, zaxis,
        method='Z-Axis Projection',
        cid_p1=0, cid_p2=0, cid_p3=0, cid_zaxis=0,
        nplanes=11)
  - plot_shear_moment_torque(model, gpforce, coord,
        idir=0, itime=0,
        nplanes=11, show=True)
  - plot_smt(x, force_sum, moment_sum, show=True)

"""
from __future__ import annotations
import numpy as np
from typing import Union, cast, TYPE_CHECKING

try:
    import matplotlib.pyplot as plt
    IS_MATPLOTLIB = True
except:
    IS_MATPLOTLIB = False

from pyNastran.bdf.mesh_utils.cut_model_by_plane import (
    get_nid_cd_xyz_cid0, get_element_centroids, get_stations)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF, CORD2R
    from pyNastran.op2.op2_geom import OP2Geom
    from pyNastran.op2.tables.ogf_gridPointForces.ogf_objects import RealGridPointForcesArray

from pyNastran.nptyping_interface import (
    NDArrayNint, NDArrayN2int, NDArray3float, NDArrayN3float)

def create_shear_moment_torque(model: BDF,
                               gpforce: RealGridPointForcesArray,
                               p1: np.ndarray,
                               p2: np.ndarray,
                               p3: np.ndarray,
                               zaxis: np.ndarray,
                               method: str='Vector',
                               station_location: str='End-Origin',
                               cid_p1: int=0, cid_p2: int=0, cid_p3: int=0, cid_zaxis: int=0,
                               nplanes: int=20,
                               root_filename=None,
                               csv_filename=None,
                               length_scale: float=1.0, length_unit: str='',
                               force_scale: float=1.0, force_unit: str='',
                               moment_scale: float=1.0, moment_unit: str='',
                               show: bool=True,
                               ) -> tuple[np.ndarray, np.ndarray]:
    """
    Creates a shear moment torque plot for the active plot result
    Plane Actor is drawn in the i-k plane

    Parameters
    ----------
    model : BDF
        the model object
    gpforce : RealGridPointForcesArray
        the results object
    p1: (3,) float ndarray
        defines the starting point for the shear, moment, torque plot
    p3: (3,) float ndarray
        defines the end point for the shear, moment, torque plot
    p2: (3,) float ndarray
        defines the XZ plane for the shears/moments
    zaxis: (3,) float ndarray
        the direction of the z-axis
    station_location : str; default='End-Origin'
        'End-Origin': magnitude along p3-p1
        'X' : Global X location
        'Y' : Global Y location
        'Z' : Global Z location
    cid_p1 / cid_p2 / cid_p3; default=0
        the coordinate systems for p1, p2, and p3
    method : str
        'coord id':
           zaxis:  N/A
           p2:     use the coord id from p2 as the output frame
        'CORD2R':
           zaxis:  point on the z-axis
           p2:     point on the xz-plane
        'Vector':
           zaxis:  k vector
           p2:     xz-plane vector
         'Z-Axis Projection':
           zaxis:  point on the z-axis
           p2:     p2 is a point on the xz-plane
    show: bool; default=True
        shows the plots

    Returns
    -------
    force_sum / moment_sum : (nstations, 3) float ndarray
        the forces/moments at the station

    """
    log = model.log
    nids, nid_cd, icd_transform, xyz_cid0 = get_nid_cd_xyz_cid0(model)
    eids, element_centroids_cid0 = get_element_centroids(model)

    xyz1, xyz2, xyz3, i, k, coord, iaxis_march, dim_max, stations = setup_coord_from_plane(
        model, xyz_cid0,
        p1, p2, p3, zaxis,
        method=method,
        cid_p1=cid_p1, cid_p2=cid_p2, cid_p3=cid_p3, cid_zaxis=cid_zaxis,
        nplanes=nplanes,
    )
    #is_failed, stations, coord, iaxis_march = self.plot_plane(
        #model, xyz_cid0,
        #p1, p2, p3,
        #zaxis, method=method,
        #cid_p1=cid_p1, cid_p2=cid_p2, cid_p3=cid_p3, cid_zaxis=cid_zaxis,
        #nplanes=nplanes,
        #stop_on_failure=stop_on_failure)
    #if is_failed:
        #force_sum = np.zeros((0, 3))
        #moment_sum = np.zeros((0, 3))
        #return force_sum, moment_sum

    coord = cast(CORD2R, coord)
    force_sum, moment_sum, new_coords, nelems, nnodes = gpforce.shear_moment_diagram(
        nids, xyz_cid0, nid_cd, icd_transform,
        eids, element_centroids_cid0,
        stations, model.coords, coord,
        iaxis_march=iaxis_march,
        itime=0, idir=0,
        nodes_tol=None, debug=False, log=log)

    force_sum *= force_scale
    moment_sum *= moment_scale

    cids, origins, xyz_stations, xlabel = get_xyz_stations(
        stations*length_scale, new_coords,
        station_location=station_location)
    origins *= length_scale

    if csv_filename:
        write_smt_to_csv(
            csv_filename,
            stations, nelems, nnodes, cids, origins,
            force_sum, moment_sum,
            length_unit=length_unit,
            force_unit=force_unit,
            moment_unit=moment_unit)
    plot_smt(
        xyz_stations,
        force_sum, moment_sum,
        nelems, nnodes, show=show,
        xtitle=xlabel, xlabel=xlabel,
        length_unit=length_unit,
        force_unit=force_unit,
        moment_unit=moment_unit,
        root_filename=root_filename,
        plot_force_components=False,
        plot_moment_components=False,
    )
    return force_sum, moment_sum

def smt_setup(model: BDF) -> tuple[NDArrayNint, NDArrayN2int, NDArrayN3float,
                                   dict[int, NDArrayNint], NDArrayNint, NDArrayN3float]:
    nids, nid_cd, icd_transform, xyz_cid0 = get_nid_cd_xyz_cid0(model)
    eids, element_centroids_cid0 = get_element_centroids(model, fdtype='float64')
    return nids, nid_cd, xyz_cid0, icd_transform, eids, element_centroids_cid0

def setup_coord_from_plane(model: Union[BDF, OP2Geom], xyz_cid0: NDArrayN3float,
                           p1: NDArray3float, p2: NDArray3float, p3: NDArray3float,
                           zaxis: NDArray3float,
                           method: str='Z-Axis Projection',
                           cid_p1: int=0, cid_p2: int=0, cid_p3: int=0, cid_zaxis: int=0,
                           nplanes: int=11, ) -> tuple[
                               np.ndarray, np.ndarray, np.ndarray,
                               np.ndarray, np.ndarray,
                               CORD2R, CORD2R, float, np.ndarray]:
    """
    Parameters
    ----------
    model : BDF
        the geometry model
    xyz_cid0: (3,nnodes) float ndarray
        the nodes locations in cid=0
    p1: (3,) float ndarray
        defines the starting point for the shear, moment, torque plot
    p3: (3,) float ndarray
        defines the end point for the shear, moment, torque plot
    p2: (3,) float ndarray
        defines the XZ plane for the shears/moments
    zaxis: (3,) float ndarray
        the direction of the z-axis
    cid_p1 / cid_p2 / cid_p3 : int
        the coordinate systems for p1, p2, and p3
    method : str
       'CORD2R':
          zaxis: point on the z-axis
          p2:     point on the xz-plane
       'Vector':
          zaxis:  k vector
          p2:     xz-plane vector
        'Z-Axis Projection':
          zaxis:  point on the z-axis
          p2:     p2 is a point on the xz-plane
    nplanes : int; default=11
        the number of planes

    Returns
    -------
    xyz1 / xyz2 / xyz3 : (3,) float ndarray
        the 1=starting 2=ending, 3=normal coordinates of the
        coordinate frames to create in the cid=0 frame
    i / k : (3,) float ndarray
        the i and k vectors of the coordinate system
    coord_out : Coord
        the generated coordinate system where the x-axis defines
        the direction to be marched
    stations : (n,) float ndarray
        the coordinates in the x-axis that will be marched down

    """
    #xyz_min, xyz_max = model.xyz_limits
    xyz_min = xyz_cid0.min(axis=0)
    xyz_max = xyz_cid0.max(axis=0)
    assert len(xyz_min) == 3, xyz_min

    dxyz = np.abs(xyz_max - xyz_min)
    dim_max = dxyz.max()
    izero = np.where(dxyz == 0)
    dxyz[izero] = dim_max
    xyz1, xyz2, xyz3, i, k, coord_out, coord_march, stations = get_stations(
        model, p1, p2, p3, zaxis,
        method=method, cid_p1=cid_p1, cid_p2=cid_p2, cid_p3=cid_p3,
        cid_zaxis=cid_zaxis, nplanes=nplanes)
    return xyz1, xyz2, xyz3, i, k, coord_out, coord_march, dim_max, stations

def plot_shear_moment_torque(model: OP2Geom,
                             gpforce: RealGridPointForcesArray,
                             coord: CORD2R,
                             itime: int=0,
                             nplanes: int=11, show: bool=True,
                             #xtitle: str='x', xlabel: str='xlabel',
                             #force_unit: str='', moment_unit: str=''
                             ) -> None:
    nids, nid_cd, xyz_cid0, icd_transform, eids, element_centroids_cid0 = smt_setup(model)
    element_centroids_coord = coord.transform_node_to_local_array(element_centroids_cid0)
    idir = 0
    x = element_centroids_coord[idir]
    xmin = x.min()
    xmax = x.max()
    dx = xmax - xmin
    assert abs(dx) > 0., f'dx={dx} xmin={xmin} xmax={xmax}'
    stations = np.linspace(0., dx, num=nplanes, endpoint=True)
    force_sum, moment_sum, new_coords, nelems, nnodes = gpforce.shear_moment_diagram(
        nids, xyz_cid0, nid_cd, icd_transform,
        eids, element_centroids_cid0, stations,
        model.coords, coord,
        iaxis_march=None,
        itime=itime, debug=False, log=model.log)
    plot_smt(stations,
             force_sum,
             moment_sum,
             nelems, nnodes, show=show)
    return

def plot_smt(x: np.ndarray,
             force_sum: np.ndarray, moment_sum: np.ndarray,
             nelems: np.ndarray, nnodes: np.ndarray,
             plot_force_components: bool=True,
             plot_moment_components: bool=True,
             root_filename: str='',
             show: bool=True,
             xtitle: str='x', xlabel: str='xlabel',
             length_unit: str='',
             force_unit: str='',
             moment_unit: str='') -> None:
    """plots the shear, moment, torque plots"""
    plt.close()

    #xtitle = 'Y'
    #xlabel = 'Spanwise Location, Y (in)'
    #moment_unit = 'in-kip'
    #force_unit = 'kip'

    length_unit2 = f' ({length_unit})' if length_unit else ''
    force_unit2 = f' ({force_unit})' if force_unit else ''
    moment_unit2 = f' ({moment_unit})' if moment_unit else ''

    #f, ax = plt.subplots()
    # ax = fig.subplots()
    if plot_force_components:
        fig = plt.figure(1)
        ax = fig.gca()
        ax.plot(x, force_sum[:, 0], '-*')
        ax.set_title(f'{xtitle} vs. Axial')
        ax.set_xlabel(f'{xlabel}{length_unit2}')
        ax.set_ylabel(f'Axial{force_unit2}')
        ax.grid(True)
        fig.tight_layout()

        fig = plt.figure(2)
        ax = fig.gca()
        ax.plot(x, force_sum[:, 1], '-*')
        ax.set_title(f'{xtitle} vs. Shear Y')
        ax.set_xlabel(f'{xlabel}{length_unit2}')
        ax.set_ylabel(f'Shear Y{force_unit2}')
        ax.grid(True)
        fig.tight_layout()

        fig = plt.figure(3)
        ax = fig.gca()
        ax.plot(x, force_sum[:, 2], '-*')
        ax.set_title(f'{xtitle} vs. Shear Z')
        ax.set_xlabel(f'{xlabel}{length_unit2}')
        ax.set_ylabel(f'Shear Z{force_unit2}')
        ax.grid(True)
        fig.tight_layout()

    if plot_moment_components:
        fig = plt.figure(4)
        ax = fig.gca()
        ax.plot(x, moment_sum[:, 0], '-*')
        ax.set_title(f'{xtitle} vs. Torque')
        ax.set_xlabel(f'{xlabel}{length_unit2}')
        ax.set_ylabel(f'Torque{moment_unit2}')
        ax.grid(True)
        fig.tight_layout()

        fig = plt.figure(5)
        ax = fig.gca()
        ax.plot(x, moment_sum[:, 1], '-*')
        ax.set_title(f'{xtitle} vs. Moment Y')
        ax.set_xlabel(f'{xlabel}{length_unit2}')
        ax.set_ylabel(f'Moment Y{moment_unit2}')
        ax.grid(True)
        fig.tight_layout()

        fig = plt.figure(6)
        ax = fig.gca()
        ax.plot(x, moment_sum[:, 2], '-*')
        ax.set_title(f'{xtitle} vs. Moment Z')
        ax.set_xlabel(f'{xlabel}{length_unit2}')
        ax.set_ylabel(f'Moment Z{moment_unit2}')
        ax.grid(True)
        fig.tight_layout()

    #-----------------------------------------------
    fig = plt.figure(7)
    ax = fig.gca()
    ax.plot(x, force_sum[:, 0], '-*', label='Force X')
    ax.plot(x, force_sum[:, 1], '-*', label='Force Y')
    ax.plot(x, force_sum[:, 2], '-*', label='Force Z')
    #ax.set_title(f'{xtitle} vs. Force')
    ax.set_xlabel(f'{xlabel}{length_unit2}')
    ax.set_ylabel(f'Force{force_unit2}')
    ax.legend()
    ax.grid(True)
    fig.tight_layout()
    if root_filename:
        fig.savefig(f'{root_filename}_forces.png')

    fig = plt.figure(8)
    ax = fig.gca()
    ax.plot(x, moment_sum[:, 0], '-*', label='Torque')
    ax.plot(x, moment_sum[:, 1], '-*', label='Moment Y')
    ax.plot(x, moment_sum[:, 2], '-*', label='Moment Z')
    #ax.set_title(f'{xtitle} vs. Moment')
    ax.set_xlabel(f'{xlabel}{length_unit2}')
    ax.set_ylabel(f'Moment{moment_unit2}')
    ax.legend()
    ax.grid(True)
    fig.tight_layout()
    if root_filename:
        fig.savefig(f'{root_filename}_moments.png')
    #-----------------------------------------------
    fig = plt.figure(9)
    ax = fig.gca()
    ax.plot(x, nnodes / nnodes.max(), '-*', label=f'n_nodes (N={nnodes.max()})')
    ax.plot(x, nelems / nelems.max(), '-*', label=f'n_elems (N={nelems.max()})')
    ax.set_title('Monotonic Nodes/Elements')
    ax.set_xlabel(f'{xlabel}{length_unit2}')
    ax.set_ylabel('Fraction of Nodes, Elements')
    ax.legend()
    ax.grid(True)
    fig.tight_layout()
    if root_filename:
        fig.savefig(f'{root_filename}_nodes_elements.png')

    if show:
        plt.show()

def get_xyz_stations(stations: np.ndarray,
                     new_coords: dict[int, CORD2R],
                     station_location: str='End-Origin',
                     ) -> tuple[list[int], np.ndarray, np.ndarray, str]:
    """helper to simplify xlabel on SMT output"""
    cids_list = []
    origins_list = []
    for cid, coord in sorted(new_coords.items()):
        cids_list.append(cid)
        origins_list.append(coord.origin)
    origins = np.vstack(origins_list)

    station_location_lower = station_location.lower()
    if station_location_lower == 'end-origin':
        xyz_stations = stations
        xlabel = 'i Station'
    elif station_location_lower == 'x':
        xyz_stations = origins[:, 0]
        xlabel = 'X'
    elif station_location_lower == 'y':
        xyz_stations = origins[:, 1]
        xlabel = 'Y'
    elif station_location_lower == 'z':
        xyz_stations = origins[:, 2]
        xlabel = 'Z'
    else:  # pragma: no cover
        raise RuntimeError(station_location_lower)
    return cids_list, origins, xyz_stations, xlabel

def write_smt_to_csv(csv_filename: str,
                     stations: np.ndarray,
                     nelems: np.ndarray, nnodes: np.ndarray,
                     cids: list[int],
                     origins: np.ndarray,
                     force_sum: np.ndarray,
                     moment_sum: np.ndarray,
                     length_unit: str='',
                     force_unit: str='',
                     moment_unit: str='') -> None:
    """writes the shear, moment, torque data"""
    length_label = f'({length_unit})' if length_unit else ''
    force_label = f'({force_unit})' if force_unit else ''
    moment_label = f'({moment_unit})' if moment_unit else ''
    with open(csv_filename, 'w') as csv_file:
        header = (
            f'Station{length_label},nelements,nnodes,'
            f'coord_id,origin_x{length_label},origin_y{length_label},origin_z{length_label},'
            f'Fx{force_label},Fy{force_label},Fz{force_label},'
            f'Mx{moment_label},My{moment_label},Mz{moment_label}\n')
        csv_file.write(header)
        for station, nelem, nnode, coord_id, origin, force_sumi, moment_sumi in zip(
            stations, nelems, nnodes, cids, origins, force_sum, moment_sum):
            csv_file.write(
                f'{station},{nelem:d},{nnode:d},{coord_id:d},'
                f'{origin[0]},{origin[1]},{origin[2]},'
                f'{force_sumi[0]},{force_sumi[1]},{force_sumi[2]},'
                f'{moment_sumi[0]},{moment_sumi[1]},{moment_sumi[2]}\n')
