from __future__ import annotations
import warnings
from typing import TYPE_CHECKING
import numpy as np

if TYPE_CHECKING:  # pragma: no cover
    #from pyNastran.dev.bdf_vectorized3.bdf import BDF
    #from pyNastran.bdf.cards.materials import MAT1, MAT8
    from pyNastran.dev.bdf_vectorized3.cards.grid import GRID


Quality = tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray,
                np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]


def tri_quality_nodes(grid: GRID, nodes: np.ndarray) -> Quality:
    """
    gets the quality metrics for a tri

    area, max_skew, aspect_ratio, min_theta, max_theta, dideal_theta, min_edge_length
    """
    xyz = grid.xyz_cid0()
    nid = grid.node_id
    inode = np.searchsorted(nid, nodes)
    actual_nodes = nid[inode]
    assert np.array_equal(actual_nodes, nodes)
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
    #a = xyz1 - xyz2
    #b = xyz1 - xyz3

    #normal = np.cross(a, b)
    #assert normal.shape[0] == nelements

    #norm = np.linalg.norm(normal, axis=1)
    #area = 0.5 * norm
    #assert len(area) == nelements

    e1 = (p1 + p2) / 2.
    e2 = (p2 + p3) / 2.
    e3 = (p3 + p1) / 2.

    #    3
    #    / \
    # e3/   \ e2
    #  /    /\
    # /    /  \
    # 1---/----2
    #    e1
    e21 = e2 - e1
    e31 = e3 - e1
    e32 = e3 - e2

    e3_p2 = e3 - p2
    e2_p1 = e2 - p1
    e1_p3 = e1 - p3

    v21 = p2 - p1
    v32 = p3 - p2
    v13 = p1 - p3
    length21 = np.linalg.norm(v21, axis=1)
    length32 = np.linalg.norm(v32, axis=1)
    length13 = np.linalg.norm(v13, axis=1)
    lengths = np.vstack([length21,
                         length32,
                         length13]).T
    assert lengths.shape == (nelements, 3), lengths.shape
    min_edge_length = lengths.min(axis=1)
    normal = np.cross(v21, v13)
    area = 0.5 * np.linalg.norm(normal, axis=1)

    ne31 = np.linalg.norm(e31, axis=1)
    ne21 = np.linalg.norm(e21, axis=1)
    ne32 = np.linalg.norm(e32, axis=1)
    ne2_p1 = np.linalg.norm(e2_p1, axis=1)
    ne3_p2 = np.linalg.norm(e3_p2, axis=1)
    ne1_p3 = np.linalg.norm(e1_p3, axis=1)

    ne2_p1__ne31 = ne2_p1 * ne31
    ne3_p2__ne21 = ne3_p2 * ne21
    ne1_p3__ne32 = ne1_p3 * ne32

    # (32, ) = (32, 3) o (32 o 3)
    #te2_p1 = np.einsum('ij,ij->i', e2_p1, e31)
    #assert len(te2_p1) == nelements

    cos_skew1 = np.einsum('ij,ij->i', e2_p1,  e31) / ne2_p1__ne31
    cos_skew2 = np.einsum('ij,ij->i', e2_p1, -e31) / ne2_p1__ne31
    cos_skew3 = np.einsum('ij,ij->i', e3_p2,  e21) / ne3_p2__ne21
    cos_skew4 = np.einsum('ij,ij->i', e3_p2, -e21) / ne3_p2__ne21
    cos_skew5 = np.einsum('ij,ij->i', e1_p3,  e32) / ne1_p3__ne32
    cos_skew6 = np.einsum('ij,ij->i', e1_p3, -e32) / ne1_p3__ne32
    assert len(cos_skew1) == nelements

    skews = np.abs(np.arccos(np.clip([
        cos_skew1, cos_skew2, cos_skew3,
        cos_skew4, cos_skew5, cos_skew6], -1., 1.)))
    assert skews.shape == (6, nelements), skews.shape
    #max_skew = np.pi / 2. - np.abs(np.arccos(np.clip([
        #cos_skew1, cos_skew2, cos_skew3,
        #cos_skew4, cos_skew5, cos_skew6], -1., 1.))).min(axis=0)
    max_skew = np.pi / 2. - skews.min(axis=0)

    l21 = np.linalg.norm(v21, axis=1).reshape(nelements, 1)
    l32 = np.linalg.norm(v32, axis=1).reshape(nelements, 1)
    l13 = np.linalg.norm(v13, axis=1).reshape(nelements, 1)
    lengths = np.hstack([l21, l32, l13])
    #lengths = np.linalg.norm(lengthsi, axis=1)
    assert lengths.shape == (nelements, 3), lengths.shape
    #assert lengths.shape == (nelements, ), lengths.shape
    #del lengthsi

    #assert len(lengths) == 3, lengths
    length_max = lengths.max(axis=1)
    length_min = lengths.min(axis=1)
    assert length_min.shape == (nelements, ), length_min.shape

    izero = (length_min == 0.0)
    ipos = ~izero

    aspect_ratio = np.full(nelements, np.nan, dtype='float64')
    # assume length_min = length21 = nan, so:
    #   cos_theta1 = nan
    #   thetas = [nan, b, c]
    #   min_theta = max_theta = dideal_theta = nan
    min_theta = np.full(nelements, np.nan, dtype='float64')
    max_theta = np.full(nelements, np.nan, dtype='float64')
    dideal_theta = np.full(nelements, np.nan, dtype='float64')
    #if izero:
        #aspect_ratio[izero] = np.nan
        ## assume length_min = length21 = nan, so:
        ##   cos_theta1 = nan
        ##   thetas = [nan, b, c]
        ##   min_theta = max_theta = dideal_theta = nan
        #min_theta[izero] = np.nan
        #max_theta[izero] = np.nan
        #dideal_theta[izero] = np.nan
    #else:
    PIOVER3 = np.pi / 3.

    ## TODO: make this not fail... this to not fail
    length_ratio = length_max / length_min
    try:
        aspect_ratio[ipos] = length_ratio[ipos]
    except ValueError:
        warnings.warn(f'len(aspect_ratio)={len(aspect_ratio)}; len(ipos)={len(ipos)}; ipos.sum()={ipos.sum()} len(length_max/min)={len(length_ratio)}')
        raise

    length_cos1 = length21 * length13
    length_cos2 = length32 * length21
    length_cos3 = length13 * length32
    cos_theta1 = np.einsum('ij,ij->i', v21, -v13) / length_cos1
    cos_theta2 = np.einsum('ij,ij->i', v32, -v21) / length_cos2
    cos_theta3 = np.einsum('ij,ij->i', v13, -v32) / length_cos3
    thetas = np.arccos(np.clip([cos_theta1, cos_theta2, cos_theta3], -1., 1.))

    # https://scicomp.stackexchange.com/questions/25012/calculate-jacobian-of-triangular-element-given-coordinates-of-vertices-and-displ
    #j11 = xb - xa
    #j12 = xc - xa

    #j21 = yb - ya
    #j22 = yc - ya
    #det_j = j11 * y22 - j12 * j21
    #jacobian = det_j

    # old
    #min_theta[ipos] = thetas.min(axis=0)
    #max_theta[ipos] = thetas.max(axis=0)
    #dideal_theta[ipos] = np.maximum(max_theta[ipos] - PIOVER3, PIOVER3 - min_theta[ipos])

    # new
    ## TODO: why do we need to use ipos twice?
    min_theta[ipos] = thetas.min(axis=0)[ipos]
    max_theta[ipos] = thetas.max(axis=0)[ipos]
    #dideal_theta[ipos] = np.maximum(max_theta - PIOVER3, PIOVER3 - min_theta)[ipos]       # more work
    dideal_theta[ipos] = np.maximum(max_theta[ipos] - PIOVER3, PIOVER3 - min_theta[ipos])  # less work

    #theta_deg = np.degrees(np.arccos(max_cos_theta))
    #if theta_deg < 60.:
        #print('p1=%s' % xyz_cid0[p1, :])
        #print('p2=%s' % xyz_cid0[p2, :])
        #print('p3=%s' % xyz_cid0[p3, :])
        #print('theta1=%s' % np.degrees(np.arccos(cos_theta1)))
        #print('theta2=%s' % np.degrees(np.arccos(cos_theta2)))
        #print('theta3=%s' % np.degrees(np.arccos(cos_theta3)))
        #print('max_theta=%s' % theta_deg)
        #asdf

    test = False
    if test:
        from pyNastran.bdf.mesh_utils.delete_bad_elements import tri_quality as tri_quality_old
        for p1i, p2i, p3i, areai, max_skewi, aspect_ratioi, min_thetai, max_thetai, dideal_thetai, min_edge_lengthi in zip(
            p1, p2, p3, area, max_skew, aspect_ratio, min_theta, max_theta, dideal_theta, min_edge_length):
            (areai_old, max_skewi_old, aspect_ratioi_old, min_thetai_old, max_thetai_old,
             dideal_thetai_old, min_edge_lengthi_old) = tri_quality_old(p1i, p2i, p3i)
            assert np.allclose(min_edge_lengthi, min_edge_lengthi_old)
            assert np.allclose(areai, areai_old)

            if np.isnan(aspect_ratioi_old):
                assert np.isnan(aspect_ratioi)
            else:
                assert np.allclose(aspect_ratioi, aspect_ratioi_old)

            assert np.allclose(max_skewi, max_skewi_old), f'skews={skews}'
            assert np.allclose(min_thetai, min_thetai_old)
            assert np.allclose(max_thetai, max_thetai_old)
            assert np.allclose(dideal_thetai, dideal_thetai_old)

    taper_ratio = np.full(nelements, np.nan, dtype=area.dtype)
    area_ratio = np.full(nelements, np.nan, dtype=area.dtype)
    max_warp = np.full(nelements, np.nan, dtype=area.dtype)
    out = (area, taper_ratio, area_ratio, np.degrees(max_skew), aspect_ratio,
           np.degrees(min_theta), np.degrees(max_theta), np.degrees(dideal_theta),
           min_edge_length, np.degrees(max_warp))
    return out
    #return area, max_skew, aspect_ratio, min_theta, max_theta, dideal_theta, min_edge_length


def quad_quality_nodes(grid: GRID, nodes: np.ndarray) -> Quality:
    """
    gets the quality metrics for a quad

    area, max_skew, aspect_ratio, min_theta, max_theta, dideal_theta, min_edge_length
    """
    xyz = grid.xyz_cid0()
    nid = grid.node_id
    nelements, nnodes = nodes.shape
    assert nnodes == 4, nodes.shape
    inode = np.searchsorted(nid, nodes)
    actual_nodes = nid[inode]
    assert np.array_equal(actual_nodes, nodes)
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
    length21 = np.linalg.norm(v21, axis=1)
    length32 = np.linalg.norm(v32, axis=1)
    length43 = np.linalg.norm(v43, axis=1)
    length14 = np.linalg.norm(v14, axis=1)
    assert length21.shape == (nelement, )
    lengths = np.array([length21, length32, length43, length14])
    min_edge_length = lengths.min(axis=0)
    assert min_edge_length.shape == (nelement, )

    p12 = (p1 + p2) / 2.
    p23 = (p2 + p3) / 2.
    p34 = (p3 + p4) / 2.
    p14 = (p4 + p1) / 2.
    v31 = p3 - p1
    v42 = p4 - p2
    normal = np.cross(v31, v42)
    area = 0.5 * np.linalg.norm(normal, axis=1)

    # still kind of in development
    #
    # the ratio of the ideal area to the actual area
    # this is an hourglass check
    areas = np.vstack([
        np.linalg.norm(np.cross(-v14, v21), axis=1), # v41 x v21
        np.linalg.norm(np.cross(v32, -v21), axis=1), # v32 x v12
        np.linalg.norm(np.cross(v43, -v32), axis=1), # v43 x v23
        np.linalg.norm(np.cross(v14,  v43), axis=1),  # v14 x v43
    ]).T
    assert areas.shape == (nelement, 4)
    #
    # for:
    #   area=1; area1=0.5 -> area_ratioi1=2.0; area_ratio=2.0
    #   area=1; area1=2.0 -> area_ratioi2=2.0; area_ratio=2.0
    max_area = areas.max(axis=1)
    min_area = areas.min(axis=1)
    assert len(min_area) == nelement

    is_nan_area_ratio = False
    min_area_min = min_area.min()
    if min_area_min > 0.0:
        area_ratioi1 = area / min_area
    else:
        assert min_area_min == 0., min_area
        area_ratioi1 = np.full(nelement, np.nan, dtype=area.dtype)
        iarea = np.where(min_area != 0.)[0]
        area_ratioi1[area] = area[area] / min_area[iarea]
        is_nan_area_ratio = True
        warnings.warn(f'invalid area_ratio for a quad/hexa/penta because min_area=0')

    area_min = area.min()
    if area_min > 0.0:
        area_ratioi2 = max_area / area
    else:
        assert area_min == 0.0, area
        area_ratioi2 = np.full(nelement, np.nan, dtype=area.dtype)
        iarea = np.where(area != 0.)[0]
        area_ratioi2[iarea] = max_area[iarea] / area[iarea]
        is_nan_area_ratio = True
        warnings.warn(f'invalid area_ratio for a quad/hexa/penta because area=0')

    area_ratios = np.array([area_ratioi1, area_ratioi2])
    area_ratio = area_ratios.max(axis=0)
    if is_nan_area_ratio:
        inan = np.isnan(area_ratio)
        area_ratio[inan] = area_ratio[~inan].max() * 2
    assert area_ratio.shape == (nelement, )
    del area_ratios

    area1 = 0.5 * np.linalg.norm(np.cross(-v14, v21), axis=1) # v41 x v21
    area2 = 0.5 * np.linalg.norm(np.cross(-v21, v32), axis=1) # v12 x v32
    area3 = 0.5 * np.linalg.norm(np.cross(v43, v32), axis=1)  # v43 x v32
    area4 = 0.5 * np.linalg.norm(np.cross(v14, -v43), axis=1) # v14 x v34
    aavg = (area1 + area2 + area3 + area4) / 4.
    taper_ratio = (abs(area1 - aavg) + abs(area2 - aavg) +
                   abs(area3 - aavg) + abs(area4 - aavg)) / aavg
    assert len(taper_ratio) == nelement

    #    e3
    # 4-------3
    # |       |
    # |e4     |  e2
    # 1-------2
    #     e1
    e13 = p34 - p12
    e42 = p23 - p14
    ne42 = np.linalg.norm(e42, axis=1)
    ne13 = np.linalg.norm(e13, axis=1)
    ne13_ne42 = ne13 * ne42
    ne13_ne42_min = ne13_ne42.min()
    if ne13_ne42_min > 0.0:
        cos_skew1 = np.einsum('ij,ij->i', e13,  e42) / ne13_ne42
        cos_skew2 = np.einsum('ij,ij->i', e13, -e42) / ne13_ne42
    #except FloatingPointError:
    else:
        assert ne13_ne42_min == 0.0, ne13_ne42_min
        ine = np.where(ne13_ne42 != 0.)[0]
        cos_skew1 = np.full(nelement, np.nan, dtype=area.dtype)
        cos_skew2 = np.full(nelement, np.nan, dtype=area.dtype)
        cos_skew1[ine] = np.einsum('ij,ij->i', e13[ine, :],  e42[ine, :]) / ne13_ne42[ine]
        cos_skew2[ine] = np.einsum('ij,ij->i', e13[ine, :], -e42[ine, :]) / ne13_ne42[ine]
        warnings.warn(f'invalid skew for a quad/hexa/penta because there are collapsed edges')

    skews = np.abs(np.arccos(np.clip([
        cos_skew1, cos_skew2], -1., 1.)))
    assert skews.shape == (2, nelement), skews.shape
    max_skew = np.pi / 2. - skews.min(axis=0)
    assert len(max_skew) == nelement

    #aspect_ratio = max(p12, p23, p34, p14) / max(p12, p23, p34, p14)
    #lengths = np.linalg.norm([v21, v32, v43, v14], axis=1)
    l21 = np.linalg.norm(v21, axis=1).reshape(nelement, 1)
    l32 = np.linalg.norm(v32, axis=1).reshape(nelement, 1)
    l43 = np.linalg.norm(v43, axis=1).reshape(nelement, 1)
    l14 = np.linalg.norm(v14, axis=1).reshape(nelement, 1)
    lengths = np.hstack([l21, l32, l43, l14])
    #lengths = np.linalg.norm(lengthsi, axis=1)
    assert lengths.shape == (nelement, 4), lengths.shape

    #assert len(lengths) == 3, lengths
    length_max = lengths.max(axis=1)
    length_min = lengths.min(axis=1)
    aspect_ratio = length_max / length_min
    assert len(aspect_ratio) == nelement

    cos_theta1 = np.einsum('ij,ij->i', v21, -v14) / (length21 * length14)
    cos_theta2 = np.einsum('ij,ij->i', v32, -v21) / (length32 * length21)
    cos_theta3 = np.einsum('ij,ij->i', v43, -v32) / (length43 * length32)
    cos_theta4 = np.einsum('ij,ij->i', v14, -v43) / (length14 * length43)
    #max_thetai = np.arccos([cos_theta1, cos_theta2, cos_theta3, cos_theta4]).max()

    # dot the local normal with the normal vector
    # then take the norm of that to determine the angle relative to the normal
    # then take the sign of that to see if we're pointing roughly towards the normal
    #
    # np.sign(np.linalg.norm(np.dot(
    # a x b = ab sin(theta)
    # a x b / ab = sin(theta)
    # sin(theta) < 0. -> normal is flipped
    v21_v32 = np.cross(v21, v32)
    v32_v43 = np.cross(v32, v43)
    v43_v14 = np.cross(v43, v14)
    v14_v21 = np.cross(v14, v21)
    normal2 = np.sign(np.einsum('ij,ij->i', v21_v32, normal))
    normal3 = np.sign(np.einsum('ij,ij->i', v32_v43, normal))
    normal4 = np.sign(np.einsum('ij,ij->i', v43_v14, normal))
    normal1 = np.sign(np.einsum('ij,ij->i', v14_v21, normal))
    n = np.array([normal1, normal2, normal3, normal4]) # (4,n)
    theta_additional = np.where(n < 0, 2*np.pi, 0.)    # (4,n)

    cos_thetas = np.array([cos_theta1, cos_theta2, cos_theta3, cos_theta4])
    theta = n * np.arccos(np.clip(
        cos_thetas, -1., 1.)) + theta_additional
    min_theta = theta.min(axis=0)
    max_theta = theta.max(axis=0)

    dideal_theta = np.maximum(max_theta - PIOVER2, PIOVER2 - min_theta)
    assert len(min_theta) == nelement
    assert len(max_theta) == nelement
    assert len(dideal_theta) == nelement
    #print('theta_max = ', theta_max)


    # warp angle
    # split the quad and find the normals of each triangl
    # find the angle between the two triangles
    #
    # 4---3
    # | / |
    # |/  |
    # 1---2
    #
    v41 = -v14
    n123 = np.cross(v21, v31)
    n134 = np.cross(v31, v41)
    length123 = np.linalg.norm(n123, axis=1)
    length134 = np.linalg.norm(n134, axis=1)
    #v1 o v2 = v1 * v2 cos(theta)
    length_warp1 = length123 * length134
    cos_warp1 = np.einsum('ij,ij->i', n123, n134) / length_warp1

    # split the quad in the order direction and take the maximum of the two splits
    # 4---3
    # | \ |
    # |  \|
    # 1---2
    n124 = np.cross(v21, v41)
    n234 = np.cross(v32, v42)
    legth_124 = np.linalg.norm(n124, axis=1)
    legth_234 = np.linalg.norm(n234, axis=1)
    length_warp2 = legth_124 * legth_234
    cos_warp2 = np.einsum('ij,ij->i', n124, n234) / length_warp2

    warps = np.abs(np.arccos(
        np.clip([cos_warp1, cos_warp2], -1., 1.)))
    max_warp = warps.max(axis=0)
    assert len(max_warp) == nelement

    #https://math.stackexchange.com/questions/2430691/jacobian-determinant-for-bi-linear-quadrilaterals
    #https://scicomp.stackexchange.com/questions/35881/fem-shape-functions-on-triangular-elements-transition-from-2d-to-3d
    #dx_d

    if np.isnan(max_warp.max()):
        print(f'NAN max_warp = {max_warp}')
    #for warp in warps:
        #print(warp.tolist())

    test = False
    if test:  # pragma: no cover
        from pyNastran.bdf.mesh_utils.delete_bad_elements import quad_quality as quad_quality_old
        for p1i, p2i, p3i, p4i, areai, taper_ratioi, area_ratioi, max_skewi, aspect_ratioi, \
            min_thetai, max_thetai, dideal_thetai, min_edge_lengthi, max_warpi in zip(
            p1, p2, p3, p4, area, taper_ratio, area_ratio, max_skew, aspect_ratio,
            min_theta, max_theta, dideal_theta, min_edge_length, max_warp):
            out_old = quad_quality_old('dummy', p1i, p2i, p3i, p4i)

            (areai_old, taper_ratioi_old, area_ratioi_old, max_skewi_old, aspect_ratioi_old,
             min_thetai_old, max_thetai_old, dideal_thetai_old, min_edge_lengthi_old, max_warpi_old) = out_old
            assert np.allclose(min_edge_lengthi, min_edge_lengthi_old)
            assert np.allclose(areai, areai_old)

            if np.isnan(aspect_ratioi_old):
                assert np.isnan(aspect_ratioi)
            else:
                assert np.allclose(aspect_ratioi, aspect_ratioi_old)
            assert np.allclose(max_warpi, max_warpi_old, atol=0.0001)
            assert np.allclose(taper_ratioi, taper_ratioi_old)
            assert np.allclose(area_ratioi, area_ratioi_old)

            assert np.allclose(max_skewi, max_skewi_old), f'skews={skews}'
            assert np.allclose(min_thetai, min_thetai_old)
            assert np.allclose(max_thetai, max_thetai_old)
            assert np.allclose(dideal_thetai, dideal_thetai_old)
        #return area, max_skew, aspect_ratio, min_theta, max_theta, dideal_theta, min_edge_length

    out = (area, taper_ratio, area_ratio, np.degrees(max_skew), aspect_ratio,
           np.degrees(min_theta), np.degrees(max_theta), np.degrees(dideal_theta),
           min_edge_length, np.degrees(max_warp))
    return out
