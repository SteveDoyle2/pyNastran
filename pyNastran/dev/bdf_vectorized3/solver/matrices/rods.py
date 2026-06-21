from __future__ import annotations
import copy
from itertools import count
from typing import Any, TYPE_CHECKING

import numpy as np
from scipy.sparse import dok_matrix

from cpylog import SimpleLogger

from ..utils import lambda1d, DOF_MAP


if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.nptyping_interface import NDArrayNNfloat
    from pyNastran.dev.bdf_vectorized3.bdf import (
        BDF,
    )
    from pyNastran.dev.bdf_vectorized3.bdf_interface.bdf_attributes import (
        CROD, PROD, CONROD, CTUBE, PTUBE,)


def _build_kbb_rod(
        model: BDF,
        element_name: str,
        Kbb: dok_matrix,
        dof_map: DOF_MAP,
        xyz_cid0: np.ndarray) -> int:
    """fill the CROD Kbb matrix"""
    if element_name == 'CROD':
        elem = model.crod
        if len(elem) == 0:
            return 0
        prop = model.prod
        pids = elem.property_id
        prop2 = prop.slice_card_by_id(pids, assume_sorted=True)
        A = prop2.area()
        J = prop2.J
    elif element_name == 'CTUBE':
        elem = model.ctube
        if len(elem) == 0:
            return 0
        prop = model.ptube
        pids = elem.property_id
        prop2: PTUBE = prop.slice_card_by_id(pids, assume_sorted=True)
        A = prop2.area()
        J = prop2.J()
    elif element_name == 'CONROD':
        elem = model.conrod
        if len(elem) == 0:
            return 0
        prop2 = elem
        A = elem.area()
        J = elem.J
    else:  # pragma: no cover
        raise NotImplementedError(element_name)
    assert isinstance(A, np.ndarray), (element_name, A)
    assert isinstance(J, np.ndarray), (element_name, J)

    mat1 = model.mat1

    nodes = elem.nodes
    inid = model.grid.index(nodes)
    inid1 = inid[:, 0]
    inid2 = inid[:, 1]
    xyz1 = xyz_cid0[inid1, :]
    xyz2 = xyz_cid0[inid2, :]
    dxyz = xyz2 - xyz1
    L = np.linalg.norm(dxyz, axis=1)

    #L = elem.length()
    assert len(L) == elem.n
    if L.min() <= 0.0:
        ifailed = (L <= 0.0)
        eids_failed = eids[ifailed]
        raise FatalError(f'{elem.type} length must be greater than 0 for eids={eids}')

    material_id = prop2.material_id
    imat = mat1.index(material_id)
    assert elem.n == prop2.n
    assert elem.n == len(imat)
    G = mat1.G[imat]
    E = mat1.E[imat]
    
    neids = len(elem)
    i = np.arange(neids)

    # 1) primary direction
    vx = dxyz / L[:, np.newaxis]

    # 2) find the max(abs(vx)) in each row
    vx_abs = np.abs(vx)
    imax0 = np.argmax(vx_abs, axis=1)
    imax = (i, imax0)
    assert len(imax0) == neids

    # 3) define a 0 matrix (delta) and add the max
    #    value to it. Offset the column entry by 1
    jmax = (i, imax0-1)
    delta = np.zeros((neids, 3), dtype=imax0.dtype)
    delta[jmax] = vx_abs[imax]
    vy_og = vx + delta

    # 4) cross to get vz
    vz = np.cross(vx, vy_og, axis=1)
    assert vx.shape == vz.shape, (vx.shape, vz.shape)
    vz_norm = np.linalg.norm(vz, axis=1)
    assert len(vz_norm) == neids
    vz /= vz_norm[:, np.newaxis]

    # 5) cross to get vy
    vy = np.cross(vz, vx, axis=1)
    vy_norm = np.linalg.norm(vy, axis=1)
    vy /= vy_norm[:, np.newaxis]

    k_axial = A * E / L
    k_torsion = G * J / L
    T3 = np.dstack([vx, vy, vz])
    assert T3.shape == (neids, 3, 3), T3.shape
    for eid, nodes, T3i, k_axiali, k_torsioni in zip(elem.element_id, elem.nodes, T3, k_axial, k_torsion):
        assert T3i.shape == (3, 3), T3i.shape
        _build_kbbi_conrod_crod(eid, Kbb, dof_map, nodes, T3i, k_axiali, k_torsioni)
    return len(elem)


def _build_kbbi_conrod_crod(
    eid: int,
    Kbb: dok_matrix,
    dof_map: DOF_MAP,
    nodes,
    T3: np.ndarray,
    k_axial: float,
    k_torsion: float,
    fdtype: str = "float64",) -> None:
    """
    fill the ith rod Kbb matrix
    # k_axial = A * E / L
    # k_torsion = G * J / L
    """
    nid1, nid2 = nodes
    log = SimpleLogger(level='debug')
    
    # 1D rod; element coordinate system
    k = np.array([[1.0, -1.0], [-1.0, 1.0]])
    k3 = np.array([
        [1.0, -1.0, 0.0],
        [-1.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
    ])  # 1D rod; element coordinate system
    k6 = np.zeros((6, 6))
    k6[0, 0] = k6[3, 3] = 1.0
    k6[0, 3] = k6[3, 0] = -1.0

    # Tr = np.row_stack([vx, vy, vz])
    # Tv = np.vstack([vx, vy, vz])
    z = np.zeros((3, 3))
    T6 = np.block([[T3, z], [z, T3]])
    # K = Lambda.T @ k @ Lambda
    K3 = T3.T @ k3 @ T3
    K6 = T6 @ k6 @ T6.T
    # print(K)

    # K = np.array([
    # [1., -1., 0.],
    # [-1., 1., 0.],
    # [0., 0., 0.]])

    # axial + torsion; assume 3D
    # u1fx, u1fy, u1fz, u2fx, u2fy, u2fz
    Ka = K6 * k_axial

    # u1mx, u1my, u1mz, u2mx, u2my, u2mz
    Kt = K6 * k_torsion

    # print(K)
    np.set_printoptions(linewidth=180)
    # print(Ka)

    # index in Ka
    idofs = np.array([
        0, 1, 2,
        3, 4, 5,
    ], dtype="int32",)

    # mapping of dofs
    ni1 = dof_map[(nid1, 1)]
    nj2 = dof_map[(nid2, 1)]
    n_ijv = np.array([
        ni1,
        ni1 + 1,
        ni1 + 2,  # node 1
        nj2,
        nj2 + 1,
        nj2 + 2,  # node 2
    ], dtype="int32")
 
    # axial
    for dof1, i1 in zip(idofs, n_ijv):
        for dof2, i2 in zip(idofs, n_ijv):
            ki = Ka[dof1, dof2]
            if abs(ki) > 0.0:
                log.debug(f"axial {ni1} {nj2} dof1/2=({dof1}, {dof2}); k={ki}")
                #print("axial", ni1, nj2, f"({dof1}, {dof2});", (dof1, dof2), ki)
                Kbb[i1, i2] += ki
    # print(Kbb.todense())

    # torsion
    n_ijv += 3
    for dof1, i1 in zip(idofs, n_ijv):
        for dof2, i2 in zip(idofs, n_ijv):
            ki = Kt[dof1, dof2]
            if abs(ki) > 0.0:
                log.debug(f"torsion {ni1} {nj2} i1/2=({i1}, {i2}); dof1/2=({dof1}, {dof2}); k={ki}")
                Kbb[i1, i2] += ki
        # print(K2)
    #print(Kbb.todense())
    return

