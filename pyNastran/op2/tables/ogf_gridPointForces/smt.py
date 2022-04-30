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

from pyNastran.bdf.mesh_utils.cut_model_by_plane import (
    get_nid_cd_xyz_cid0, get_element_centroids, get_stations)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF, CORD2R
    from pyNastran.op2.op2_geom import OP2Geom
    from pyNastran.op2.tables.ogf_gridPointForces.ogf_objects import RealGridPointForcesArray

from pyNastran.nptyping_interface import NDArrayNint, NDArrayN2int, NDArray3float, NDArrayN3float

def smt_setup(model: BDF):
    nids, nid_cd, icd_transform, xyz_cid0 = get_nid_cd_xyz_cid0(model)
    eids, element_centroids_cid0 = get_element_centroids(model)
    return nids, nid_cd, xyz_cid0, icd_transform, eids, element_centroids_cid0

def setup_coord_from_plane(model: Tuple[BDF, OP2Geom], xyz_cid0: NDArrayN3float,
                           p1: NDArray3float, p2: NDArray3float, p3: NDArray3float,
                           zaxis: NDArray3float,
                           method: str='Z-Axis Projection',
                           cid_p1: int=0, cid_p2: int=0, cid_p3: int=0, cid_zaxis: int=0,
                           nplanes: int=11, ):
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
    xyz1, xyz2, xyz3, i, k, coord_out, stations = get_stations(
        model, p1, p2, p3, zaxis,
        method=method, cid_p1=cid_p1, cid_p2=cid_p2, cid_p3=cid_p3,
        cid_zaxis=cid_zaxis, idir=0, nplanes=nplanes)
    return xyz1, xyz2, xyz3, i, k, coord_out, dim_max, stations

def plot_shear_moment_torque(model: OP2Geom,
                             gpforce: RealGridPointForcesArray,
                             coord: CORD2R,
                             idir: int=0, itime: int=0,
                             nplanes: int=11, show: bool=True):
    nids, nid_cd, xyz_cid0, icd_transform, eids, element_centroids_cid0 = smt_setup(model)
    element_centroids_cid = coord.transform_node_to_local_array(element_centroids_cid0)
    x = element_centroids_cid[idir]
    xmin = x.min()
    xmax = x.max()
    dx = xmax - xmin
    assert abs(dx) > 0., f'dx={dx} xmin={xmin} xmax={xmax}'
    stations = np.linspace(0., dx, num=nplanes, endpoint=True)
    force_sum, moment_sum = gpforce.shear_moment_diagram(
        xyz_cid0, eids, nids, icd_transform,
        element_centroids_cid0,
        model.coords, nid_cd, stations, coord,
        idir=idir, itime=itime, debug=False, log=model.log)
    plot_smt(stations, force_sum, moment_sum, show=show)

def plot_smt(x, force_sum, moment_sum, show=True):
    """plots the shear, moment, torque plots"""
    import matplotlib.pyplot as plt
    plt.close()
    #f, ax = plt.subplots()
    # ax = fig.subplots()
    fig = plt.figure(1)
    ax = fig.gca()
    ax.plot(x, force_sum[:, 0], '-*')
    ax.set_xlabel('X')
    ax.set_ylabel('Axial')
    ax.grid(True)

    fig = plt.figure(2)
    ax = fig.gca()
    ax.plot(x, force_sum[:, 1], '-*')
    ax.set_xlabel('X')
    ax.set_ylabel('Shear Y')
    ax.grid(True)

    fig = plt.figure(3)
    ax = fig.gca()
    ax.plot(x, force_sum[:, 2], '-*')
    ax.set_xlabel('X')
    ax.set_ylabel('Shear Z')
    ax.grid(True)

    fig = plt.figure(4)
    ax = fig.gca()
    ax.plot(x, moment_sum[:, 0], '-*')
    ax.set_xlabel('X')
    ax.set_ylabel('Torque')
    ax.grid(True)

    fig = plt.figure(5)
    ax = fig.gca()
    ax.plot(x, moment_sum[:, 1], '-*')
    ax.set_xlabel('X')
    ax.set_ylabel('Moment Y')
    ax.grid(True)

    fig = plt.figure(6)
    ax = fig.gca()
    ax.plot(x, moment_sum[:, 2], '-*')
    ax.set_xlabel('X')
    ax.set_ylabel('Moment Z')
    ax.grid(True)

    if show:
        plt.show()
