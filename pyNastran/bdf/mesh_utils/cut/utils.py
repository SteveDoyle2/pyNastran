from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

if TYPE_CHECKING:
    from pyNastran.bdf.bdf import BDF
    from pyNastran.nptyping_interface import NDArray3float


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
