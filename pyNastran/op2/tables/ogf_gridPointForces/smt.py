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
from typing import TYPE_CHECKING

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

from pyNastran.nptyping_interface import NDArrayNint, NDArrayN2int, NDArray3float, NDArrayN3float

def smt_setup(model: BDF) -> tuple[NDArrayNint, NDArrayN2int, NDArrayN3float,
                                   dict[int, NDArrayNint], NDArrayNint, NDArrayN3float]:
    nids, nid_cd, icd_transform, xyz_cid0 = get_nid_cd_xyz_cid0(model)
    eids, element_centroids_cid0 = get_element_centroids(model, fdtype='float64')
    return nids, nid_cd, xyz_cid0, icd_transform, eids, element_centroids_cid0

def setup_coord_from_plane(model: tuple[BDF, OP2Geom], xyz_cid0: NDArrayN3float,
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
        'Z-Axis Projection'
           p1-p2 defines the x-axis
           k is defined by the z-axis
       'CORD2R' : typical
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
             force_unit: str='', moment_unit: str='') -> None:
    """plots the shear, moment, torque plots"""
    plt.close()

    #xtitle = 'Y'
    #xlabel = 'Spanwise Location, Y (in)'
    #moment_unit = 'in-kip'
    #force_unit = 'kip'

    force_unit2 = f' ({force_unit})' if force_unit else ''
    moment_unit2 = f' ({moment_unit})' if moment_unit else ''

    #f, ax = plt.subplots()
    # ax = fig.subplots()
    if plot_force_components:
        fig = plt.figure(1)
        ax = fig.gca()
        ax.plot(x, force_sum[:, 0], '-*')
        ax.set_title(f'{xtitle} vs. Axial')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(f'Axial{force_unit2}')
        ax.grid(True)
        fig.tight_layout()

        fig = plt.figure(2)
        ax = fig.gca()
        ax.plot(x, force_sum[:, 1], '-*')
        ax.set_title(f'{xtitle} vs. Shear Y')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(f'Shear Y{force_unit2}')
        ax.grid(True)
        fig.tight_layout()

        fig = plt.figure(3)
        ax = fig.gca()
        ax.plot(x, force_sum[:, 2], '-*')
        ax.set_title(f'{xtitle} vs. Shear Z')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(f'Shear Z{force_unit2}')
        ax.grid(True)
        fig.tight_layout()

    if plot_moment_components:
        fig = plt.figure(4)
        ax = fig.gca()
        ax.plot(x, moment_sum[:, 0], '-*')
        ax.set_title(f'{xtitle} vs. Torque')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(f'Torque{moment_unit2}')
        ax.grid(True)
        fig.tight_layout()

        fig = plt.figure(5)
        ax = fig.gca()
        ax.plot(x, moment_sum[:, 1], '-*')
        ax.set_title(f'{xtitle} vs. Moment Y')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(f'Moment Y{moment_unit2}')
        ax.grid(True)
        fig.tight_layout()

        fig = plt.figure(6)
        ax = fig.gca()
        ax.plot(x, moment_sum[:, 2], '-*')
        ax.set_title(f'{xtitle} vs. Moment Z')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(f'Moment Z{moment_unit2}')
        ax.grid(True)
        fig.tight_layout()

    #-----------------------------------------------
    fig = plt.figure(7)
    ax = fig.gca()
    ax.plot(x, force_sum[:, 0], '-*', label=f'Force X{force_unit2}')
    ax.plot(x, force_sum[:, 1], '-*', label=f'Force Y{force_unit2}')
    ax.plot(x, force_sum[:, 2], '-*', label=f'Force Z{force_unit2}')
    #ax.set_title(f'{xtitle} vs. Force')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(f'Force{force_unit2}')
    ax.legend()
    ax.grid(True)
    fig.tight_layout()
    if root_filename:
        fig.savefig(f'{root_filename}_forces.png')

    fig = plt.figure(8)
    ax = fig.gca()
    ax.plot(x, moment_sum[:, 0], '-*', label=f'Torque{moment_unit2}')
    ax.plot(x, moment_sum[:, 1], '-*', label=f'Moment Y{moment_unit2}')
    ax.plot(x, moment_sum[:, 2], '-*', label=f'Moment Z{moment_unit2}')
    #ax.set_title(f'{xtitle} vs. Moment')
    ax.set_xlabel(xlabel)
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
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Fraction of Nodes, Elements')
    ax.legend()
    ax.grid(True)
    fig.tight_layout()
    if root_filename:
        fig.savefig(f'{root_filename}_nodes_elements.png')

    if show:
        plt.show()
