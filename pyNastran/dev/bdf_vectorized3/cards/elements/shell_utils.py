from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

from .solid_volume import volume_chexa, volume_cpenta

if TYPE_CHECKING:
    #from pyNastran.dev.bdf_vectorized3.bdf import BDF
    #from pyNastran.bdf.cards.materials import MAT1, MAT8
    from pyNastran.dev.bdf_vectorized3.cards.grid import GRID


def tri_area(grid: GRID, nodes: np.ndarray) -> np.ndarray:
    xyz = grid.xyz_cid0()
    nid = grid.node_id
    nelements, nnodes = nodes.shape
    assert nnodes == 3, nodes.shape
    inode = np.searchsorted(nid, nodes)
    actual_nodes = nid[inode]
    assert np.array_equal(actual_nodes, nodes)
    in1 = inode[:, 0]
    in2 = inode[:, 1]
    in3 = inode[:, 2]
    xyz1 = xyz[in1, :]
    xyz2 = xyz[in2, :]
    xyz3 = xyz[in3, :]
    a = xyz1 - xyz2
    b = xyz1 - xyz3
    area = 0.5 * np.linalg.norm(np.cross(a, b), axis=1)
    assert len(area) == nelements
    return area

def tri_centroid(grid: GRID, nodes: np.ndarray) -> np.ndarray:
    xyz = grid.xyz_cid0()
    nid = grid.node_id
    inode = np.searchsorted(nid, nodes)
    assert np.array_equal(nid[inode], nodes)
    in1 = inode[:, 0]
    in2 = inode[:, 1]
    in3 = inode[:, 2]
    xyz1 = xyz[in1, :]
    xyz2 = xyz[in2, :]
    xyz3 = xyz[in3, :]
    centroid = (xyz1 + xyz2 + xyz3) / 3.
    return centroid

def tri_area_centroid_normal(grid: GRID, nodes: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    xyz = grid.xyz_cid0()
    nid = grid.node_id
    nelements, nnodes = nodes.shape
    assert nnodes == 3, nodes.shape
    inode = np.searchsorted(nid, nodes)
    actual_nodes = nid[inode]
    assert np.array_equal(actual_nodes, nodes)
    in1 = inode[:, 0]
    in2 = inode[:, 1]
    in3 = inode[:, 2]
    xyz1 = xyz[in1, :]
    xyz2 = xyz[in2, :]
    xyz3 = xyz[in3, :]
    a = xyz1 - xyz2
    b = xyz1 - xyz3

    normal = np.cross(a, b)
    assert normal.shape[0] == nelements

    norm = np.linalg.norm(normal, axis=1)
    area = 0.5 * norm
    assert len(area) == nelements

    centroid = (xyz1 + xyz2 + xyz3) / 3.
    assert centroid.shape[0] == nelements

    assert normal.shape == (nelements, 3)
    ipos = norm > 0.
    ibad = ~ipos
    if ibad.sum() == 0:
        unit_normal = normal / norm[:, np.newaxis]
    else:
        unit_normal = np.full(normal.shape, np.nan, normal.dtype)
        unit_normal[ipos, :] = normal[ipos, :] / norm[ipos, np.newaxis]
        #print('norm =', norm)
        #print('unit_normal =', unit_normal)
        print('a[bad,:] =', a[ibad, :])
        print('  n1 =', nodes[ibad, 0])
        print('  n2 =', nodes[ibad, 1])
        print('  n3 =', nodes[ibad, 2])
        print('  xyz1 =', xyz1[ibad, :])
        print('  xyz2 =', xyz2[ibad, :])
        print('  xyz3 =', xyz3[ibad, :])
        print('b[bad,:] =', b[ibad, :])

    assert unit_normal.shape == (nelements, 3)
    return area, centroid, unit_normal

def quad_area(grid: GRID, nodes: np.ndarray) -> np.ndarray:
    nelements, nnodes = nodes.shape
    assert nnodes == 4, nnodes
    nid = grid.node_id
    xyz = grid.xyz_cid0()
    inode = np.searchsorted(nid, nodes)
    assert np.array_equal(nid[inode], nodes)
    in1 = inode[:, 0]
    in2 = inode[:, 1]
    in3 = inode[:, 2]
    in4 = inode[:, 3]
    xyz1 = xyz[in1, :]
    xyz2 = xyz[in2, :]
    xyz3 = xyz[in3, :]
    xyz4 = xyz[in4, :]
    #area = 0.5 * norm(cross(n3-n1, n4-n2))
    a = xyz3 - xyz1
    b = xyz4 - xyz2
    area = 0.5 * np.linalg.norm(np.cross(a, b), axis=1)
    assert len(area) == nelements
    return area

def quad_centroid(grid: GRID, nodes: np.ndarray) -> np.ndarray:
    nelements, nnodes = nodes.shape
    assert nnodes == 4, nnodes
    xyz = grid.xyz_cid0()
    nid = grid.node_id
    inode = np.searchsorted(nid, nodes)
    assert np.array_equal(nid[inode], nodes)
    in1 = inode[:, 0]
    in2 = inode[:, 1]
    in3 = inode[:, 2]
    in4 = inode[:, 3]
    xyz1 = xyz[in1, :]
    xyz2 = xyz[in2, :]
    xyz3 = xyz[in3, :]
    xyz4 = xyz[in4, :]
    centroid = (xyz1 + xyz2 + xyz3 + xyz4) / 4.
    assert centroid.shape[0] == nelements
    return centroid

def quad_area_centroid_normal(grid: GRID, nodes: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    nelements, nnodes = nodes.shape
    assert nnodes == 4, nnodes
    nid = grid.node_id
    xyz = grid.xyz_cid0()
    inode = np.searchsorted(nid, nodes)
    assert np.array_equal(nid[inode], nodes)
    in1 = inode[:, 0]
    in2 = inode[:, 1]
    in3 = inode[:, 2]
    in4 = inode[:, 3]
    xyz1 = xyz[in1, :]
    xyz2 = xyz[in2, :]
    xyz3 = xyz[in3, :]
    xyz4 = xyz[in4, :]
    #area = 0.5 * norm(cross(n3-n1, n4-n2))
    a = xyz3 - xyz1
    b = xyz4 - xyz2
    normal = np.cross(a, b)
    assert normal.shape[0] == nelements

    norm = np.linalg.norm(normal, axis=1)
    area = 0.5 * norm
    assert len(area) == nelements

    centroid = (xyz1 + xyz2 + xyz3 + xyz4) / 4.
    assert centroid.shape[0] == nelements

    assert normal.shape == (nelements, 3)
    unit_normal = normal / norm[:, np.newaxis]

    assert unit_normal.shape == (nelements, 3)
    return area, centroid, unit_normal


def tri_volume(grid: GRID,
               nodes: np.ndarray,
               dthickness: np.ndarray) -> np.ndarray:
    """
    Exact quad volume considering differential thicknesses

    Parameters
    ----------
    dthickness : (nelement, 3) float ndarray
        differential thickness on the quad
    """

    nelements, nnodes = nodes.shape
    assert nnodes == 3, nnodes
    nid = grid.node_id
    xyz = grid.xyz_cid0()
    inode = np.searchsorted(nid, nodes)
    assert np.array_equal(nid[inode], nodes)
    in1 = inode[:, 0]
    in2 = inode[:, 1]
    in3 = inode[:, 2]
    xyz1 = xyz[in1, :]
    xyz2 = xyz[in2, :]
    xyz3 = xyz[in3, :]

    a = xyz2 - xyz1
    b = xyz3 - xyz1
    normal = np.cross(a, b)
    assert normal.shape[0] == nelements

    norm = np.linalg.norm(normal, axis=1)
    assert normal.shape == (nelements, 3)
    unit_normal = normal / norm[:, np.newaxis]

    iequal = (dthickness.max(axis=0) == dthickness.min(axis=0))
    inot_equal = ~iequal
    nequal = iequal.sum()
    nnot_equal = inot_equal.sum()
    volume = np.full(nelements, np.nan, dtype='float64')
    if nequal:
        thicknessi = dthickness[iequal, 0]
        areai = 0.5 * norm[iequal]
        volumei = areai * thicknessi
        volume[iequal] = volumei

    if nnot_equal:
        unit_normali = unit_normal[inot_equal, :]
        xyz4 = xyz1[inot_equal, :] + unit_normali * dthickness[inot_equal, 0]
        xyz5 = xyz2[inot_equal, :] + unit_normali * dthickness[inot_equal, 1]
        xyz6 = xyz3[inot_equal, :] + unit_normali * dthickness[inot_equal, 2]
        volumei = volume_cpenta(
            xyz1[inot_equal, :], xyz2[inot_equal, :], xyz3[inot_equal, :],
            xyz4, xyz5, xyz6)
        volume[iequal] = volumei
    return volume

def quad_volume(grid: GRID,
                nodes: np.ndarray,
                dthickness: np.ndarray) -> np.ndarray:
    """
    Exact quad volume considering differential thicknesses

    Parameters
    ----------
    dthickness : (nelement, 4) float ndarray
        differential thickness on the quad
    """

    nelements, nnodes = nodes.shape
    assert nnodes == 4, nnodes
    nid = grid.node_id
    xyz = grid.xyz_cid0()
    inode = np.searchsorted(nid, nodes)
    assert np.array_equal(nid[inode], nodes)
    in1 = inode[:, 0]
    in2 = inode[:, 1]
    in3 = inode[:, 2]
    in4 = inode[:, 3]
    xyz1 = xyz[in1, :]
    xyz2 = xyz[in2, :]
    xyz3 = xyz[in3, :]
    xyz4 = xyz[in4, :]

    a = xyz3 - xyz1
    b = xyz4 - xyz2
    normal = np.cross(a, b)
    assert normal.shape[0] == nelements

    norm = np.linalg.norm(normal, axis=1)
    assert normal.shape == (nelements, 3)
    unit_normal = normal / norm[:, np.newaxis]

    iequal = (dthickness.max(axis=0) == dthickness.min(axis=0))
    inot_equal = ~iequal
    nequal = iequal.sum()
    nnot_equal = inot_equal.sum()
    volume = np.full(nelements, np.nan, dtype='float64')
    if nequal:
        thicknessi = dthickness[iequal, 0]
        areai = 0.5 * norm[iequal]
        volumei = areai * thicknessi
        volume[iequal] = volumei

    if nnot_equal:
        unit_normali = unit_normal[inot_equal, :]
        xyz5 = xyz1[inot_equal, :] + unit_normali * dthickness[inot_equal, 0]
        xyz6 = xyz2[inot_equal, :] + unit_normali * dthickness[inot_equal, 1]
        xyz7 = xyz3[inot_equal, :] + unit_normali * dthickness[inot_equal, 2]
        xyz8 = xyz4[inot_equal, :] + unit_normali * dthickness[inot_equal, 3]
        volumei = volume_chexa(xyz1[inot_equal, :], xyz2[inot_equal, :],
                               xyz3[inot_equal, :], xyz4[inot_equal, :],
                               xyz5, xyz6, xyz7, xyz8)
        volume[iequal] = volumei
    return volume
