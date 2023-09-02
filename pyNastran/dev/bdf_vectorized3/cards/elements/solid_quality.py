from __future__ import annotations
from typing import TYPE_CHECKING

import numpy as np
from .shell_quality import tri_quality_xyz, quad_quality_xyz, Quality # tri_quality_xyz0,
if TYPE_CHECKING:
    from .solid import CTETRA, CHEXA, CPENTA, CPYRAM


def penta_quality(self: CPENTA) -> Quality:
    """
    4-5-6

    1-2-3"""
    grid = self.model.grid
    nodes = self.base_nodes

    xyz = grid.xyz_cid0()
    nid = grid.node_id
    inode = np.searchsorted(nid, nodes)
    actual_nodes = nid[inode]
    assert np.array_equal(actual_nodes, nodes)

    tri_faces = [
        (3, 4, 5),
        (0, 1, 2),
    ]
    quad_faces = [
        (0, 1, 4, 3),
        (1, 2, 5, 4),
        (2, 0, 3, 5),
    ]
    out = _tri_quad_quality(tri_faces, quad_faces, inode, xyz)
    return out


def pyram_quality(self: CPYRAM) -> Quality:
    """
       5

    1-2-3-4
    """
    grid = self.model.grid
    nodes = self.base_nodes

    xyz = grid.xyz_cid0()
    nid = grid.node_id
    inode = np.searchsorted(nid, nodes)
    actual_nodes = nid[inode]
    assert np.array_equal(actual_nodes, nodes)

    tri_faces = [
        (0, 1, 4),
        (1, 2, 4),
        (2, 3, 4),
        (3, 0, 4),
    ]
    quad_faces = [(0, 1, 2, 3)]
    out = _tri_quad_quality(tri_faces, quad_faces, inode, xyz)
    return out

def tetra_quality(self: CTETRA) -> Quality:
    grid = self.model.grid
    nodes = self.base_nodes

    xyz = grid.xyz_cid0()
    nid = grid.node_id
    inode = np.searchsorted(nid, nodes)
    actual_nodes = nid[inode]
    assert np.array_equal(actual_nodes, nodes)
    faces = [
        (0, 1, 2),
        (0, 1, 3),
        (1, 2, 3),
        (2, 0, 3),
    ]
    #out = _tri_quality(faces, inode, xyz)

    nelements = inode.shape[0]
    area = np.full(nelements, np.nan, dtype='float64')
    taper_ratio = np.full((nelements, 4), np.nan, dtype='float64')
    area_ratio = np.full((nelements, 4), np.nan, dtype='float64')
    max_skew = np.full((nelements, 4), np.nan, dtype='float64')
    aspect_ratio = np.full((nelements, 4), np.nan, dtype='float64')
    min_theta = np.full((nelements, 4), np.nan, dtype='float64')
    max_theta = np.full((nelements, 4), np.nan, dtype='float64')
    dideal_theta = np.full((nelements, 4), np.nan, dtype='float64')
    min_edge_length = np.full((nelements, 4), np.nan, dtype='float64')
    max_warp = np.full((nelements, 4), np.nan, dtype='float64')
    for i, face in enumerate(faces):
        qualityi = tri_quality_xyz(xyz, inode[:, face])
        (areai, taper_ratioi, area_ratioi, max_skewi, aspect_ratioi,
         min_thetai, max_thetai, dideal_thetai, min_edge_lengthi, max_warpi) = qualityi
        taper_ratio[:, i] = taper_ratioi
        area_ratio[:, i] = area_ratioi
        max_skew[:, i] = max_skewi
        aspect_ratio[:, i] = aspect_ratioi
        min_theta[:, i] = min_thetai
        max_theta[:, i] = max_thetai
        dideal_theta[:, i] = dideal_thetai
        min_edge_length[:, i] = min_edge_lengthi
        max_warp[:, i] = max_warpi

    taper_ratio = np.max(taper_ratio, axis=1)
    area_ratio = np.max(area_ratio, axis=1)
    max_skew = np.max(max_skew, axis=1)
    aspect_ratio = np.max(aspect_ratio, axis=1)
    min_theta = np.min(min_theta, axis=1)
    max_theta = np.max(max_theta, axis=1)
    dideal_theta = np.max(dideal_theta, axis=1)
    min_edge_length = np.max(min_edge_length, axis=1)
    max_warp = np.max(max_warp, axis=1)
    assert len(taper_ratio) == nelements

    out = (area, taper_ratio, area_ratio, max_skew, aspect_ratio,
           min_theta, max_theta, dideal_theta, min_edge_length, max_warp)
    return out


def chexa_quality(self: CHEXA) -> Quality:
    grid = self.model.grid
    nodes = self.base_nodes

    xyz = grid.xyz_cid0()
    nid = grid.node_id
    inode = np.searchsorted(nid, nodes)
    actual_nodes = nid[inode]
    assert np.array_equal(actual_nodes, nodes)

    quad_faces = [
        (0, 1, 2, 3),  # 1, 2, 3, 4 # inward
        (0, 1, 5, 4),  # 1, 2, 6, 5
        (1, 2, 6, 5),  # 2, 3, 7, 6
        (2, 3, 7, 6),  # 3, 4, 8, 7
        (3, 0, 4, 7),  # 4, 1, 5, 8
        (4, 5, 6, 7),  # 5, 6, 7, 8
    ]
    out = _quad_quality(quad_faces, inode, xyz)
    #nelements = inode.shape[0]
    #area = np.full(nelements, np.nan, dtype='float64')
    #taper_ratio = np.full((nelements, 6), np.nan, dtype='float64')
    #area_ratio = np.full((nelements, 6), np.nan, dtype='float64')
    #max_skew = np.full((nelements, 6), np.nan, dtype='float64')
    #aspect_ratio = np.full((nelements, 6), np.nan, dtype='float64')
    #min_theta = np.full((nelements, 6), np.nan, dtype='float64')
    #max_theta = np.full((nelements, 6), np.nan, dtype='float64')
    #dideal_theta = np.full((nelements, 6), np.nan, dtype='float64')
    #min_edge_length = np.full((nelements, 6), np.nan, dtype='float64')
    #max_warp = np.full((nelements, 6), np.nan, dtype='float64')
    #for i, face in enumerate(faces):
        ##if len(face) == 3:
        #in1 = inode[:, face[0]]
        #in2 = inode[:, face[1]]
        #in3 = inode[:, face[2]]
        #in4 = inode[:, face[3]]
        #n1 = xyz[in1, :]
        #n2 = xyz[in2, :]
        #n3 = xyz[in3, :]
        #n4 = xyz[in4, :]
        #qualityi = quad_quality_xyz(n1, n2, n3, n4)
        #(areai, taper_ratioi, area_ratioi, max_skewi, aspect_ratioi,
         #min_thetai, max_thetai, dideal_thetai, min_edge_lengthi, max_warpi) = qualityi
        #taper_ratio[:, i] = taper_ratioi
        #area_ratio[:, i] = area_ratioi
        #max_skew[:, i] = max_skewi
        #aspect_ratio[:, i] = aspect_ratioi
        #min_theta[:, i] = min_thetai
        #max_theta[:, i] = max_thetai
        #dideal_theta[:, i] = dideal_thetai
        #min_edge_length[:, i] = min_edge_lengthi
        #max_warp[:, i] = max_warpi

    #taper_ratio = np.max(taper_ratio, axis=1)
    #area_ratio = np.max(area_ratio, axis=1)
    #max_skew = np.max(max_skew, axis=1)
    #aspect_ratio = np.max(aspect_ratio, axis=1)
    #min_theta = np.min(min_theta, axis=1)
    #max_theta = np.max(max_theta, axis=1)
    #dideal_theta = np.max(dideal_theta, axis=1)
    #min_edge_length = np.max(min_edge_length, axis=1)
    #max_warp = np.max(max_warp, axis=1)
    #assert len(taper_ratio) == nelements

    #out = (area, taper_ratio, area_ratio, max_skew, aspect_ratio,
           #min_theta, max_theta, dideal_theta, min_edge_length, max_warp)
    return out


def _tri_quad_quality(tri_faces: list[tuple[int, int, int]],
                      quad_faces: list[tuple[int, int, int, int]],
                      inode: np.ndarray, xyz: np.ndarray) -> Quality:
    quality3 = _tri_quality(tri_faces, inode, xyz)
    quality4 = _quad_quality(quad_faces, inode, xyz)
    (area3, taper_ratio3, area_ratio3, max_skew3, aspect_ratio3,
     min_theta3, max_theta3, dideal_theta3, min_edge_length3, unused_max_warp3) = quality3
    (area4, taper_ratio4, area_ratio4, max_skew4, aspect_ratio4,
     min_theta4, max_theta4, dideal_theta4, min_edge_length4, max_warp4) = quality4
    area = (area3 + area4) / 2
    taper_ratio     = np.column_stack([taper_ratio3, taper_ratio4]).max(axis=1)
    area_ratio      = np.column_stack([area_ratio3, area_ratio4]).max(axis=1)
    max_skew        = np.column_stack([max_skew3, max_skew4]).max(axis=1)
    aspect_ratio    = np.column_stack([aspect_ratio3, aspect_ratio4]).max(axis=1)
    min_theta       = np.column_stack([min_theta3, min_theta4]).min(axis=1)
    max_theta       = np.column_stack([max_theta3, max_theta4]).max(axis=1)
    dideal_theta    = np.column_stack([dideal_theta3, dideal_theta4]).max(axis=1)
    min_edge_length = np.column_stack([min_edge_length3, min_edge_length4]).max(axis=1)
    max_warp = max_warp4
    assert len(min_edge_length) == len(min_edge_length3)

    out = (area, taper_ratio, area_ratio, max_skew, aspect_ratio,
           min_theta, max_theta, dideal_theta, min_edge_length, max_warp)
    return out


def _tri_quality(tri_faces: list[tuple[int, int, int]],
                 inode: np.ndarray,
                 xyz: np.ndarray) -> Quality:
    nelements = len(inode)
    nfaces = len(tri_faces)

    area = np.full(nelements, np.nan, dtype='float64')
    taper_ratio = np.full((nelements, nfaces), np.nan, dtype='float64')
    area_ratio = np.full((nelements, nfaces), np.nan, dtype='float64')
    max_skew = np.full((nelements, nfaces), np.nan, dtype='float64')
    aspect_ratio = np.full((nelements, nfaces), np.nan, dtype='float64')
    min_theta = np.full((nelements, nfaces), np.nan, dtype='float64')
    max_theta = np.full((nelements, nfaces), np.nan, dtype='float64')
    dideal_theta = np.full((nelements, nfaces), np.nan, dtype='float64')
    min_edge_length = np.full((nelements, nfaces), np.nan, dtype='float64')
    max_warp = np.full((nelements, nfaces), np.nan, dtype='float64')

    for i, face in enumerate(tri_faces):
        #in1 = inode[:, face[0]]
        #in2 = inode[:, face[1]]
        #in3 = inode[:, face[2]]
        #n1 = xyz[in1, :]
        #n2 = xyz[in2, :]
        #n3 = xyz[in3, :]
        #qualityi = tri_quality_xyz(n1, n2, n3)
        qualityi = tri_quality_xyz(xyz, inode[:, face])
        (areai, taper_ratioi, area_ratioi, max_skewi, aspect_ratioi,
         min_thetai, max_thetai, dideal_thetai, min_edge_lengthi, max_warpi) = qualityi
        taper_ratio[:, i] = taper_ratioi
        area_ratio[:, i] = area_ratioi
        max_skew[:, i] = max_skewi
        aspect_ratio[:, i] = aspect_ratioi
        min_theta[:, i] = min_thetai
        max_theta[:, i] = max_thetai
        dideal_theta[:, i] = dideal_thetai
        min_edge_length[:, i] = min_edge_lengthi
        max_warp[:, i] = max_warpi

    taper_ratio = np.max(taper_ratio, axis=1)
    area_ratio = np.max(area_ratio, axis=1)
    max_skew = np.max(max_skew, axis=1)
    aspect_ratio = np.max(aspect_ratio, axis=1)
    min_theta = np.min(min_theta, axis=1)
    max_theta = np.max(max_theta, axis=1)
    dideal_theta = np.max(dideal_theta, axis=1)
    min_edge_length = np.max(min_edge_length, axis=1)
    max_warp = np.max(max_warp, axis=1)
    assert len(taper_ratio) == nelements

    out = (area, taper_ratio, area_ratio, max_skew, aspect_ratio,
           min_theta, max_theta, dideal_theta, min_edge_length, max_warp)
    return out

def _quad_quality(quad_faces: list[tuple[int, int, int, int]],
                  inode: np.ndarray,
                  xyz: np.ndarray) -> Quality:
    nelements = len(inode)
    nfaces = len(quad_faces)

    area = np.full(nelements, np.nan, dtype='float64')
    taper_ratio = np.full((nelements, nfaces), np.nan, dtype='float64')
    area_ratio = np.full((nelements, nfaces), np.nan, dtype='float64')
    max_skew = np.full((nelements, nfaces), np.nan, dtype='float64')
    aspect_ratio = np.full((nelements, nfaces), np.nan, dtype='float64')
    min_theta = np.full((nelements, nfaces), np.nan, dtype='float64')
    max_theta = np.full((nelements, nfaces), np.nan, dtype='float64')
    dideal_theta = np.full((nelements, nfaces), np.nan, dtype='float64')
    min_edge_length = np.full((nelements, nfaces), np.nan, dtype='float64')
    max_warp = np.full((nelements, nfaces), np.nan, dtype='float64')

    for i, face in enumerate(quad_faces):
        in1 = inode[:, face[0]]
        in2 = inode[:, face[1]]
        in3 = inode[:, face[2]]
        in4 = inode[:, face[3]]
        n1 = xyz[in1, :]
        n2 = xyz[in2, :]
        n3 = xyz[in3, :]
        n4 = xyz[in4, :]
        qualityi = quad_quality_xyz(n1, n2, n3, n4)
        (areai, taper_ratioi, area_ratioi, max_skewi, aspect_ratioi,
             min_thetai, max_thetai, dideal_thetai, min_edge_lengthi, max_warpi) = qualityi
        taper_ratio[:, i] = taper_ratioi
        area_ratio[:, i] = area_ratioi
        max_skew[:, i] = max_skewi
        aspect_ratio[:, i] = aspect_ratioi
        min_theta[:, i] = min_thetai
        max_theta[:, i] = max_thetai
        dideal_theta[:, i] = dideal_thetai
        min_edge_length[:, i] = min_edge_lengthi
        max_warp[:, i] = max_warpi

    taper_ratio = np.max(taper_ratio, axis=1)
    area_ratio = np.max(area_ratio, axis=1)
    max_skew = np.max(max_skew, axis=1)
    aspect_ratio = np.max(aspect_ratio, axis=1)
    min_theta = np.min(min_theta, axis=1)
    max_theta = np.max(max_theta, axis=1)
    dideal_theta = np.max(dideal_theta, axis=1)
    min_edge_length = np.max(min_edge_length, axis=1)
    max_warp = np.max(max_warp, axis=1)
    assert len(taper_ratio) == nelements

    out = (area, taper_ratio, area_ratio, max_skew, aspect_ratio,
           min_theta, max_theta, dideal_theta, min_edge_length, max_warp)
    return out
