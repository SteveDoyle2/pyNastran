from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.bdf.cards.coordinate_systems import (
    CORD2R, # Coord,
    # xyz_to_rtz_array, rtz_to_xyz_array,
)
if TYPE_CHECKING:
    from pyNastran.bdf.bdf import BDF
    from pyNastran.nptyping_interface import NDArray3float


def get_stations(model: BDF,
                 p1: NDArray3float, p2: NDArray3float, p3: NDArray3float,
                 zaxis: NDArray3float,
                 method: str='Vector',
                 cid_p1: int=0, cid_p2: int=0, cid_p3: int=0, cid_zaxis: int=0,
                 nplanes: int=20) -> tuple[NDArray3float, NDArray3float, NDArray3float,
                                           NDArray3float, NDArray3float,
                                           CORD2R,
                                           NDArray3float, NDArrayNfloat]:
    """
    Gets the axial stations

    Parameters
    ----------
    p1: (3,) float ndarray
        defines the starting point for the shear, moment, torque plot
    p3: (3,) float ndarray
        defines the end point for the shear, moment, torque plot
    p2: (3,) float ndarray
        defines the XZ plane for the shears/moments (depends on method)
        'Vectors':
           i = p2
        'Z-Axis Projection':
           i = p2 - p1
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

    Returns
    -------
    xyz1 / xyz2 / xyz3 : (3,) float ndarray
        the 1=starting 2=ending, 3=normal coordinates of the
        coordinate frames to create in the cid=0 frame
    i / k : (3,) float ndarray
        the i and k vectors of the coordinate system
    coord_out : Coord
        the output coordinate system
    iaxis_march : (3,) float ndarray
        the normalized x-axis that defines the direction to march
    stations : (n,) float ndarray
        the coordinates in the x-axis that will be marched down

    Example
    -------
    For the BWB example, we to calculate an SMT down the global x-axis

    1--------> y
    |
    |
    |
    2, 3
    |
    v x

    # axial
    p1 = np.array([0., 0., 0.]) # origin
    p2 = np.array([1600., 0., 0.]) # xaxis
    p3 = np.array([1600., 0., 0.]) # end
    zaxis = np.array([0., 0., 1.])
    method = 'Z-Axis Projection'

    xyz1, xyz2, xyz3, i, k, coord_out, stations = get_stations(
        model, p1, p2, p3, zaxis,
        method=method, cid_p1=0, cid_p2=0, cid_p3=0,
        cid_zaxis=0, nplanes=100)
    print(stations)

    """
    p1 = np.asarray(p1).astype('float64') # start
    p2 = np.asarray(p2).astype('float64') # xz-plane
    p3 = np.asarray(p3).astype('float64') # end point
    zaxis = np.asarray(zaxis).astype('float64')

    # define a local coordinate system
    xyz1, xyz2, unused_z_global, i, k, origin, zaxis2, xzplane = p1_p2_zaxis_to_cord2r(
        model, p1, p2, zaxis,
        cid_p1=cid_p1, cid_p2=cid_p2, cid_zaxis=cid_zaxis,
        method=method)
    xyz3 = model.coords[cid_p3].transform_node_to_global(p3)

    try:
        coord_out = CORD2R(-1, origin=origin, zaxis=zaxis2, xzplane=xzplane)
    except Exception:
        msg = f'Cannot create ouput coordinate system.  origin={origin} zaxis={zaxis} xzplane={xzplane}\n'
        #msg += coord_out.get_stats()
        raise ValueError(msg)

    #coord_march = coord_out
    xaxis_march = xyz3 - xyz1
    xaxis_march_norm = np.linalg.norm(xaxis_march)
    if xaxis_march_norm == 0.:
        msg = f'Coincident starting and end points.  dx={xaxis_march_norm} xyz1={xyz1} xyz3={xyz3}\n'
        #msg += coord_out.get_stats()
        raise ValueError(msg)

    iaxis_march = xaxis_march / xaxis_march_norm
    # k has been rotated into the output coordinate frame, so we'll maintain that
    # k is length=1
    assert np.allclose(np.linalg.norm(k), 1.0)
    jaxis_march = np.cross(k, iaxis_march)
    jaxis_march_norm = np.linalg.norm(jaxis_march)
    if jaxis_march_norm == 0.:
        msg = f'Equal k axis and iaxis.  k={str(k)} iaxis_march={str(iaxis_march)}\n'
        #msg += coord_out.get_stats()
        raise ValueError(msg)

    kaxis_march = np.cross(iaxis_march, jaxis_march)
    kaxis_march_norm = np.linalg.norm(kaxis_march)
    if kaxis_march_norm == 0.:
        msg = f'Equal iaxis and jaxis.  k={str(k)} iaxis_march={str(iaxis_march)} jaxis_march={str(jaxis_march)}\n'
        #msg += coord_out.get_stats()
        raise ValueError(msg)

    coord_march = CORD2R(-1, origin=origin, zaxis=origin+kaxis_march, xzplane=origin+iaxis_march)
    #coord_march = CORD2R(None, origin=origin, zaxis=axis_march, xzplane=xzplane)
    #print(coord_out.get_stats())

    xyz1p = coord_march.transform_node_to_local(xyz1) # start
    xyz3p = coord_march.transform_node_to_local(xyz3) # end
    xaxis = xyz3p - xyz1p

    # we want to give the dx the right sign in the coord_out frame
    #i_abs = np.abs(coord_march.i)
    #i_abs_max = i_abs.max()
    #idir = np.where(i_abs == i_abs_max)[0][0]
    #isign = np.sign(coord_march.i[idir])
    dx = xaxis[0] # * isign

    if abs(dx) == 0.:
        msg = f'Coincident starting and end points.  dx={dx} xyz1={xyz1} xyz3={xyz3}\n'
        msg += coord_out.get_stats()
        raise ValueError(msg)
    x_stations_march = np.linspace(0., dx, num=nplanes, endpoint=True)
    assert x_stations_march.shape == (nplanes, ), x_stations_march.shape
    #stations.sort()
    return xyz1, xyz2, xyz3, i, k, coord_out, iaxis_march, x_stations_march


def p1_p2_zaxis_to_cord2r(model: BDF,
                          p1: NDArray3float, p2: NDArray3float, zaxis: NDArray3float,
                          method: str='Z-Axis Projection',
                          cid_p1: int=0, cid_p2: int=0,
                          cid_zaxis: int=0) -> tuple[NDArray3float, NDArray3float, NDArray3float,
                                                      NDArray3float, NDArray3float,
                                                      NDArray3float, NDArray3float, NDArray3float]:
    """
    Creates the coordinate system that will define the cutting plane

    Parameters
    ----------
    model : BDF
        the bdf object
    method: str; default='Z-Axis Projection'
        method = 'CORD2R'
           p1:    origin
           zaxis: origin + zaxis (zaxis)
           p3:    origin + xzplane
        method = 'Z-Axis Projection'
           p1:    origin
           zaxis: zaxis
           p2:    origin + xzplane
    cid_p1/p2/zaxis: int; default=0
        the coordinate system for the points

    Returns
    -------
    p1 / p2 / p3 : (3,) float ndarray
        A, B, C in the CORD2R in the rid=0 frame
    i, k : the i/k vectors for the coord_out frame
    origin, zaxis, xzplane : (3,) float ndarray
        the CORD2R-ready coordinate system
    """
    p1 = np.asarray(p1).astype('float64') # origin
    p2 = np.asarray(p2).astype('float64') # xz-plane
    zaxis = np.asarray(zaxis).astype('float64')
    #print("coord:")
    #print('  p1 =', p1)
    #print('  p2 =', p2)
    #print('  zaxis =', zaxis)

    coords = model.coords
    xyz1 = coords[cid_p1].transform_node_to_global(p1)
    method_lower = method.lower().strip()
    if method_lower == 'vector':
        origin = xyz1
        xyz2 = origin + p2
        z_global = origin + zaxis
        method_lower = 'cord2r'
    elif method_lower == 'coord id':
        coord = coords[cid_p2]
        origin = xyz1
        xyz2 = origin + coord.i
        z_global = origin + coord.k
        method_lower = 'cord2r'
    else:
        xyz2 = coords[cid_p2].transform_node_to_global(p2)
        z_global = coords[cid_zaxis].transform_node_to_global(zaxis)

    if method_lower == 'cord2r':
        origin = xyz1
        xzplane_ = xyz2 - xyz1
        zaxis_ = z_global - xyz1
        i, k = _determine_cord2r(origin, zaxis_, xzplane_)
        origin_zaxis = z_global # origin + i
        origin_xzplane = xyz2 # origin + k
    #elif method == 'Vectors':
        #xyz2, i, k, origin, origin_zaxis, origin_xzplane = _project_vectors(xyz1, xyz2, z_global)
    elif method_lower == 'z-axis projection':
        i, k, origin, origin_zaxis, origin_xzplane = project_z_axis(xyz1, xyz2, z_global)
    else:
        raise NotImplementedError(f"method={method!r}; valid_methods=['CORD2R', 'Vector', 'Z-Axis Projection']")
    #print(f'origin={origin}')
    #print(f'origin_zaxis={origin_zaxis}')
    #print(f'origin_xzplane={origin_xzplane}')
    return xyz1, xyz2, z_global, i, k, origin, origin_zaxis, origin_xzplane


def project_z_axis(p1: NDArray3float,
                   p2: NDArray3float,
                   z_global: NDArray3float) -> tuple[NDArray3float, NDArray3float, NDArray3float,
                                                     NDArray3float, NDArray3float]:
    """
    p1-p2 defines the x-axis
    k is defined by the z-axis
    k = z / |z|
    xz = p2 - p1
    xz /= |xz|
    j = z × xz
    j /= |j|
    i = k × j
    """
    p1 = np.asarray(p1).astype('float64') # origin
    p2 = np.asarray(p2).astype('float64') # xz-plane
    x = p2 - p1
    norm_x = np.linalg.norm(x)
    if norm_x == 0.:
        raise RuntimeError(f'p1={p1} and p2={p2} are coincident; distance={norm_x}')
    iprime = x / norm_x
    k = z_global / np.linalg.norm(z_global)
    j = np.cross(k, iprime)
    jhat = j / np.linalg.norm(j)
    i = np.cross(jhat, k)

    origin = p1
    zaxis = p1 + k
    xzplane = p1 + i
    return i, k, origin, zaxis, xzplane


def _determine_cord2r(origin: np.ndarray,
                      zaxis: np.ndarray,
                      xzplane: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    k = zaxis / np.linalg.norm(zaxis)
    iprime = xzplane / np.linalg.norm(xzplane)
    j = np.cross(k, iprime)
    j /= np.linalg.norm(j)
    i = np.cross(j, k)
    return i, k
