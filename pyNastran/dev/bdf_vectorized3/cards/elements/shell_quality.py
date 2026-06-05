from __future__ import annotations
import warnings
from typing import TYPE_CHECKING
import numpy as np

if TYPE_CHECKING:  # pragma: no cover
    #from pyNastran.dev.bdf_vectorized3.bdf import BDF
    #from pyNastran.bdf.cards.materials import MAT1, MAT8
    from pyNastran.dev.bdf_vectorized3.cards.grid import GRID


Quality = tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray,
                np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray,
                np.ndarray, np.ndarray, np.ndarray]


def _cross3(a: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Cross product returning separate x,y,z components (avoids column_stack)."""
    return (a[:, 1]*b[:, 2] - a[:, 2]*b[:, 1],
            a[:, 2]*b[:, 0] - a[:, 0]*b[:, 2],
            a[:, 0]*b[:, 1] - a[:, 1]*b[:, 0])


def _norm3(cx: np.ndarray, cy: np.ndarray, cz: np.ndarray) -> np.ndarray:
    """Norm of vector given as separate x,y,z components."""
    return np.sqrt(cx*cx + cy*cy + cz*cz)


def _cross_norm(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """||a x b|| without allocating the full cross product array."""
    cx, cy, cz = _cross3(a, b)
    return _norm3(cx, cy, cz)


def _vec_norm(v: np.ndarray) -> np.ndarray:
    """Row-wise norm of (n,3) array — 2x faster than np.linalg.norm(axis=1)."""
    return np.sqrt(np.einsum('ij,ij->i', v, v))


def _dot(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Row-wise dot product of (n,3) arrays."""
    return np.einsum('ij,ij->i', a, b)


def _fast_searchsorted(nid: np.ndarray, nodes: np.ndarray) -> np.ndarray:
    """Map node IDs to indices; fast path when nid = 1..N (contiguous)."""
    n = len(nid)
    if n > 0 and nid[0] == 1 and nid[-1] == n:
        # Contiguous 1-based IDs — direct index
        inode = nodes - 1
    else:
        inode = np.searchsorted(nid, nodes)
        assert np.array_equal(nid[inode], nodes)
    return inode


def tri_quality_nodes(grid: GRID, nodes: np.ndarray) -> Quality:
    """
    gets the quality metrics for a tri

    area, max_skew, aspect_ratio, min_theta, max_theta, dideal_theta, min_edge_length
    """
    xyz = grid.xyz_cid0()
    nid = grid.node_id
    inode = _fast_searchsorted(nid, nodes)
    return tri_quality_xyz(xyz, inode)

def tri_quality_xyz(xyz: np.ndarray, inode: np.ndarray) -> Quality:
    """
    gets the quality metrics for a tri

    area, max_skew, aspect_ratio, min_theta, max_theta, dideal_theta, min_edge_length
    """
    nelements, nnodes = inode.shape
    assert nnodes == 3, inode.shape
    in1 = inode[:, 0]
    in2 = inode[:, 1]
    in3 = inode[:, 2]
    p1 = xyz[in1, :]
    p2 = xyz[in2, :]
    p3 = xyz[in3, :]
    return tri_quality_xyz0(p1, p2, p3)
    # ---------------------------------------------------------

def tri_quality_xyz0(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> Quality:
    """
    gets the quality metrics for a tri

    area, max_skew, aspect_ratio, min_theta, max_theta, dideal_theta, min_edge_length
    """
    nelements = len(p1)

    e1 = (p1 + p2) / 2.
    e2 = (p2 + p3) / 2.
    e3 = (p3 + p1) / 2.

    e21 = e2 - e1
    e31 = e3 - e1
    e32 = e3 - e2

    e3_p2 = e3 - p2
    e2_p1 = e2 - p1
    e1_p3 = e1 - p3

    v21 = p2 - p1
    v32 = p3 - p2
    v13 = p1 - p3
    length21 = _vec_norm(v21)
    length32 = _vec_norm(v32)
    length13 = _vec_norm(v13)
    lengths = np.array([length21, length32, length13])  # (3, n)
    min_edge_length = lengths.min(axis=0)
    area = 0.5 * _cross_norm(v21, v13)

    ne31 = _vec_norm(e31)
    ne21 = _vec_norm(e21)
    ne32 = _vec_norm(e32)
    ne2_p1 = _vec_norm(e2_p1)
    ne3_p2 = _vec_norm(e3_p2)
    ne1_p3 = _vec_norm(e1_p3)

    ne2_p1__ne31 = ne2_p1 * ne31
    ne3_p2__ne21 = ne3_p2 * ne21
    ne1_p3__ne32 = ne1_p3 * ne32

    cos_skew1 = _dot(e2_p1,  e31) / ne2_p1__ne31
    cos_skew2 = -cos_skew1
    cos_skew3 = _dot(e3_p2,  e21) / ne3_p2__ne21
    cos_skew4 = -cos_skew3
    cos_skew5 = _dot(e1_p3,  e32) / ne1_p3__ne32
    cos_skew6 = -cos_skew5

    skews = np.abs(np.arccos(np.clip([
        cos_skew1, cos_skew2, cos_skew3,
        cos_skew4, cos_skew5, cos_skew6], -1., 1.)))
    max_skew = np.pi / 2. - skews.min(axis=0)

    length_max = lengths.max(axis=0)
    length_min = lengths.min(axis=0)

    izero = (length_min == 0.0)
    ipos = ~izero

    aspect_ratio = np.full(nelements, np.nan, dtype='float64')
    min_theta = np.full(nelements, np.nan, dtype='float64')
    max_theta = np.full(nelements, np.nan, dtype='float64')
    dideal_theta = np.full(nelements, np.nan, dtype='float64')
    PIOVER3 = np.pi / 3.

    length_ratio = length_max / length_min
    try:
        aspect_ratio[ipos] = length_ratio[ipos]
    except ValueError:
        warnings.warn(f'len(aspect_ratio)={len(aspect_ratio)}; len(ipos)={len(ipos)}; ipos.sum()={ipos.sum()} len(length_max/min)={len(length_ratio)}')
        raise

    cos_theta1 = _dot(v21, -v13) / (length21 * length13)
    cos_theta2 = _dot(v32, -v21) / (length32 * length21)
    cos_theta3 = _dot(v13, -v32) / (length13 * length32)
    thetas = np.arccos(np.clip([cos_theta1, cos_theta2, cos_theta3], -1., 1.))

    min_theta[ipos] = thetas.min(axis=0)[ipos]
    max_theta[ipos] = thetas.max(axis=0)[ipos]
    dideal_theta[ipos] = np.maximum(max_theta[ipos] - PIOVER3, PIOVER3 - min_theta[ipos])

    taper_ratio = np.full(nelements, np.nan, dtype=area.dtype)
    area_ratio = np.full(nelements, np.nan, dtype=area.dtype)
    max_warp = np.full(nelements, np.nan, dtype=area.dtype)
    nastran_skew = np.degrees(min_theta)
    nastran_taper = np.full(nelements, np.nan, dtype=area.dtype)
    nastran_warp = np.full(nelements, np.nan, dtype=area.dtype)
    out = (area, taper_ratio, area_ratio, np.degrees(max_skew), aspect_ratio,
           np.degrees(min_theta), np.degrees(max_theta), np.degrees(dideal_theta),
           min_edge_length, np.degrees(max_warp),
           nastran_skew, nastran_taper, nastran_warp)
    return out


def quad_quality_nodes(grid: GRID, nodes: np.ndarray) -> Quality:
    """
    gets the quality metrics for a quad

    area, max_skew, aspect_ratio, min_theta, max_theta, dideal_theta, min_edge_length
    """
    xyz = grid.xyz_cid0()
    nid = grid.node_id
    nelements, nnodes = nodes.shape
    assert nnodes == 4, nodes.shape
    inode = _fast_searchsorted(nid, nodes)
    in1 = inode[:, 0]
    in2 = inode[:, 1]
    in3 = inode[:, 2]
    in4 = inode[:, 3]
    p1 = xyz[in1, :]
    p2 = xyz[in2, :]
    p3 = xyz[in3, :]
    p4 = xyz[in4, :]
    out = quad_quality_xyz(p1, p2, p3, p4)
    return out

def quad_quality_xyz(p1: np.ndarray, p2: np.ndarray,
                     p3: np.ndarray, p4: np.ndarray) -> Quality:
    nelement = p1.shape[0]
    PIOVER2 = np.pi / 2.
    v21 = p2 - p1
    v32 = p3 - p2
    v43 = p4 - p3
    v14 = p1 - p4
    length21 = _vec_norm(v21)
    length32 = _vec_norm(v32)
    length43 = _vec_norm(v43)
    length14 = _vec_norm(v14)
    lengths = np.array([length21, length32, length43, length14])  # (4, n)
    min_edge_length = lengths.min(axis=0)

    p12 = (p1 + p2) / 2.
    p23 = (p2 + p3) / 2.
    p34 = (p3 + p4) / 2.
    p14 = (p4 + p1) / 2.
    v31 = p3 - p1
    v42 = p4 - p2
    nx, ny, nz = _cross3(v31, v42)
    area = 0.5 * _norm3(nx, ny, nz)
    # Keep normal as stacked array for later dot products
    normal = np.column_stack([nx, ny, nz])

    # Sub-triangle areas (used for area_ratio AND taper_ratio)
    # Compute cross-product norms once, reuse for both
    cn1 = _cross_norm(-v14, v21)  # |v41 x v21|
    cn2 = _cross_norm(-v21, v32)  # |v12 x v32|
    cn3 = _cross_norm(v43, -v32)  # |v43 x v23|
    cn4 = _cross_norm(v14, v43)   # |v14 x v43|

    # area_ratio
    areas = np.array([cn1, cn2, cn3, cn4])  # (4, n)
    max_area = areas.max(axis=0)
    min_area = areas.min(axis=0)

    is_nan_area_ratio = False
    min_area_min = min_area.min()
    if min_area_min > 0.0:
        area_ratioi1 = area / min_area
    else:
        area_ratioi1 = np.full(nelement, np.nan, dtype=area.dtype)
        iarea = np.where(min_area != 0.)[0]
        area_ratioi1[iarea] = area[iarea] / min_area[iarea]
        is_nan_area_ratio = True
        warnings.warn(f'invalid area_ratio for a quad/hexa/penta because min_area=0')

    area_min = area.min()
    if area_min > 0.0:
        area_ratioi2 = max_area / area
    else:
        area_ratioi2 = np.full(nelement, np.nan, dtype=area.dtype)
        iarea = np.where(area != 0.)[0]
        area_ratioi2[iarea] = max_area[iarea] / area[iarea]
        is_nan_area_ratio = True
        warnings.warn(f'invalid area_ratio for a quad/hexa/penta because area=0')

    area_ratio = np.maximum(area_ratioi1, area_ratioi2)
    if is_nan_area_ratio:
        inan = np.isnan(area_ratio)
        inot_nan = ~inan
        not_nan_area = area_ratio[inot_nan]
        if inot_nan.sum() > 0:
            area_ratio[inan] = not_nan_area.max() * 2

    # taper_ratio — reuse cn1..cn4 (half-areas = cn/2)
    area1 = 0.5 * cn1
    area2 = 0.5 * cn2
    area3 = 0.5 * cn3
    area4 = 0.5 * cn4
    aavg = (area1 + area2 + area3 + area4) / 4.
    taper_ratio = np.full(nelement, np.nan, dtype=area.dtype)
    num = (np.abs(area1 - aavg) + np.abs(area2 - aavg) +
           np.abs(area3 - aavg) + np.abs(area4 - aavg))
    itaper = (aavg != 0.0)
    taper_ratio[itaper] = num[itaper] / aavg[itaper]

    # skew
    e13 = p34 - p12
    e42 = p23 - p14
    ne42 = _vec_norm(e42)
    ne13 = _vec_norm(e13)
    ne13_ne42 = ne13 * ne42
    ne13_ne42_min = ne13_ne42.min()
    if ne13_ne42_min > 0.0:
        cos_skew1 = _dot(e13, e42) / ne13_ne42
        cos_skew2 = -cos_skew1
    else:
        ine = np.where(ne13_ne42 != 0.)[0]
        cos_skew1 = np.full(nelement, np.nan, dtype=area.dtype)
        cos_skew2 = np.full(nelement, np.nan, dtype=area.dtype)
        cos_skew1[ine] = _dot(e13[ine, :], e42[ine, :]) / ne13_ne42[ine]
        cos_skew2[ine] = -cos_skew1[ine]
        warnings.warn(f'invalid skew for a quad/hexa/penta because there are collapsed edges')

    skews = np.abs(np.arccos(np.clip([cos_skew1, cos_skew2], -1., 1.)))
    max_skew = np.pi / 2. - skews.min(axis=0)

    # aspect ratio (reuse lengths already computed)
    length_max = lengths.max(axis=0)
    length_min = lengths.min(axis=0)
    aspect_ratio = np.full(nelement, np.nan, dtype=area.dtype)
    iaspect = (length_min > 0)
    aspect_ratio[iaspect] = length_max[iaspect] / length_min[iaspect]

    # corner angles
    cos_theta1 = np.full(nelement, np.nan, dtype=area.dtype)
    cos_theta2 = np.full(nelement, np.nan, dtype=area.dtype)
    cos_theta3 = np.full(nelement, np.nan, dtype=area.dtype)
    cos_theta4 = np.full(nelement, np.nan, dtype=area.dtype)

    icos = (length21 * length14) > 0
    cos_theta1[icos] = _dot(v21[icos, :], -v14[icos, :]) / (length21[icos] * length14[icos])
    icos = (length32 * length21) > 0
    cos_theta2[icos] = _dot(v32[icos, :], -v21[icos, :]) / (length32[icos] * length21[icos])
    icos = (length43 * length32) > 0
    cos_theta3[icos] = _dot(v43[icos, :], -v32[icos, :]) / (length43[icos] * length32[icos])
    icos = (length14 * length43) > 0
    cos_theta4[icos] = _dot(v14[icos, :], -v43[icos, :]) / (length14[icos] * length43[icos])

    # Determine winding sign from cross products at each corner
    # Use _dot with normal instead of full np.cross
    v21_v32_x, v21_v32_y, v21_v32_z = _cross3(v21, v32)
    v32_v43_x, v32_v43_y, v32_v43_z = _cross3(v32, v43)
    v43_v14_x, v43_v14_y, v43_v14_z = _cross3(v43, v14)
    v14_v21_x, v14_v21_y, v14_v21_z = _cross3(v14, v21)

    normal2 = np.sign(v21_v32_x*nx + v21_v32_y*ny + v21_v32_z*nz)
    normal3 = np.sign(v32_v43_x*nx + v32_v43_y*ny + v32_v43_z*nz)
    normal4 = np.sign(v43_v14_x*nx + v43_v14_y*ny + v43_v14_z*nz)
    normal1 = np.sign(v14_v21_x*nx + v14_v21_y*ny + v14_v21_z*nz)
    n = np.array([normal1, normal2, normal3, normal4])  # (4, n)
    theta_additional = np.where(n < 0, 2*np.pi, 0.)

    cos_thetas = np.array([cos_theta1, cos_theta2, cos_theta3, cos_theta4])
    theta = n * np.arccos(np.clip(cos_thetas, -1., 1.)) + theta_additional
    min_theta = theta.min(axis=0)
    max_theta = theta.max(axis=0)
    dideal_theta = np.maximum(max_theta - PIOVER2, PIOVER2 - min_theta)

    # warp angle
    cos_warp1 = np.full(nelement, np.nan, dtype=area.dtype)
    cos_warp2 = np.full(nelement, np.nan, dtype=area.dtype)
    v41 = -v14
    n123_x, n123_y, n123_z = _cross3(v21, v31)
    n134_x, n134_y, n134_z = _cross3(v31, v41)
    length123 = _norm3(n123_x, n123_y, n123_z)
    length134 = _norm3(n134_x, n134_y, n134_z)
    length_warp1 = length123 * length134
    iwarp = (length_warp1 != 0.0)
    cos_warp1[iwarp] = ((n123_x[iwarp]*n134_x[iwarp] +
                         n123_y[iwarp]*n134_y[iwarp] +
                         n123_z[iwarp]*n134_z[iwarp]) / length_warp1[iwarp])

    n124_x, n124_y, n124_z = _cross3(v21, v41)
    n234_x, n234_y, n234_z = _cross3(v32, v42)
    length_124 = _norm3(n124_x, n124_y, n124_z)
    length_234 = _norm3(n234_x, n234_y, n234_z)
    length_warp2 = length_124 * length_234
    iwarp2 = (length_warp2 != 0.0)
    cos_warp2[iwarp2] = ((n124_x[iwarp2]*n234_x[iwarp2] +
                          n124_y[iwarp2]*n234_y[iwarp2] +
                          n124_z[iwarp2]*n234_z[iwarp2]) / length_warp2[iwarp2])

    warps = np.abs(np.arccos(np.clip([cos_warp1, cos_warp2], -1., 1.)))
    max_warp = warps.max(axis=0)

    # --- NX Nastran GEOMCHECK metrics ---
    # Skew (Q4_SKEW): angle between midpoint vectors (same geometry as Altair skew)
    # Already computed above as max_skew — convert to NX convention (the angle itself)
    nastran_skew = np.pi / 2.0 - max_skew  # Altair: pi/2 - min_angle; NX reports the angle

    # Taper (Q4_TAPER): (A_max - Q) / Q
    # Corner triangles: ABD, BCA, CDB, DAC (skip one vertex each)
    area_abd = 0.5 * _cross_norm(p2 - p1, p4 - p1)
    area_bca = 0.5 * _cross_norm(p3 - p2, p1 - p2)
    area_cdb = 0.5 * _cross_norm(p4 - p3, p2 - p3)
    area_dac = 0.5 * _cross_norm(p1 - p4, p3 - p4)
    corner_areas = np.array([area_abd, area_bca, area_cdb, area_dac])  # (4, n)
    a_max = corner_areas.max(axis=0)
    q_half = 0.5 * area  # Q = 0.5 * total quad area
    nastran_taper = np.full(nelement, np.nan, dtype=area.dtype)
    iq = (q_half > 0.0)
    nastran_taper[iq] = (a_max[iq] - q_half[iq]) / q_half[iq]

    # Warp (Q4_WARP): W = |AB . PK| / (D_AC + D_BD)
    # PK = (AC x BD) / |AC x BD|  (out-of-plane unit vector from diagonals)
    # AC = v31, BD = v42 (already computed)
    pk_x, pk_y, pk_z = _cross3(v31, v42)
    pk_norm = _norm3(pk_x, pk_y, pk_z)
    d_ac = _vec_norm(v31)
    d_bd = _vec_norm(v42)
    denom_warp = d_ac + d_bd

    nastran_warp = np.zeros(nelement, dtype=area.dtype)
    iwarp_nx = (pk_norm > 0.0) & (denom_warp > 0.0)
    if iwarp_nx.any():
        pk_unit_x = pk_x[iwarp_nx] / pk_norm[iwarp_nx]
        pk_unit_y = pk_y[iwarp_nx] / pk_norm[iwarp_nx]
        pk_unit_z = pk_z[iwarp_nx] / pk_norm[iwarp_nx]
        ab = v21[iwarp_nx]
        hh = ab[:, 0]*pk_unit_x + ab[:, 1]*pk_unit_y + ab[:, 2]*pk_unit_z
        nastran_warp[iwarp_nx] = np.abs(hh) / denom_warp[iwarp_nx]

    out = (area, taper_ratio, area_ratio, np.degrees(max_skew), aspect_ratio,
           np.degrees(min_theta), np.degrees(max_theta), np.degrees(dideal_theta),
           min_edge_length, np.degrees(max_warp),
           np.degrees(nastran_skew), nastran_taper, nastran_warp)
    return out
