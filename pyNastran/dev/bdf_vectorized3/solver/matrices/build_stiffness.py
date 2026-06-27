from __future__ import annotations
import copy
from itertools import count
from typing import Any, TYPE_CHECKING

import numpy as np
from scipy.sparse import coo_matrix

from cpylog import SimpleLogger
from pyNastran.utils.scipy_utils import dok_matrix, csc_matrix

from pyNastran.dev.bdf_vectorized3.solver.elements.solids import (
    build_kbb_chexa, build_kbb_ctetra, build_kbb_cpenta)
from pyNastran.dev.bdf_vectorized3.cards.base_card import searchsorted_filter

from pyNastran.dev.bdf_vectorized3.solver.elements.shells import build_kbb_cquad4, build_kbb_ctria3
from pyNastran.dev.bdf_vectorized3.solver.elements.beam import (
    timoshenko_stiffness,
    beam_transforms,
    thermal_load_beam,
    geometric_stiffness,
)
from .springs import (
    _build_kbb_celas1, _build_kbb_celas2,
    _build_kbb_celas3, _build_kbb_celas4,)

from .bush import assemble_cbush_matrices
from .rods import _build_kbb_rod
from ..utils import lambda1d, DOF_MAP


if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.nptyping_interface import NDArrayNNfloat
    from pyNastran.dev.bdf_vectorized3.bdf import (
        BDF,
        # CELAS1, CELAS2, CELAS3, CELAS4,
        # CBAR, PBAR, PBARL, PBEAM, PBEAML, # , CBEAM
        # MAT1,
    )
    from pyNastran.dev.bdf_vectorized3.bdf_interface.bdf_attributes import (
        CELAS1, CELAS2,
        CBAR, CTUBE, PTUBE,
    )


class _COOAccumulator:
    """Collects COO triplets for fast sparse matrix assembly."""

    __slots__ = ('_rows', '_cols', '_vals', '_ndof')

    def __init__(self, ndof: int) -> None:
        self._rows: list[np.ndarray] = []
        self._cols: list[np.ndarray] = []
        self._vals: list[np.ndarray] = []
        self._ndof = ndof

    def add_matrix(self, n_ijv: list | np.ndarray, K: np.ndarray) -> None:
        """Add a dense element matrix at the given DOF indices."""
        n_ijv = np.asarray(n_ijv, dtype='int32')
        ndof_e = len(n_ijv)
        rows = np.repeat(n_ijv, ndof_e)
        cols = np.tile(n_ijv, ndof_e)
        vals = K.ravel()
        mask = vals != 0.0
        self._rows.append(rows[mask])
        self._cols.append(cols[mask])
        self._vals.append(vals[mask])

    def add_scalar(self, i: int, j: int, value: float) -> None:
        """Add a single scalar entry."""
        self._rows.append(np.array([i], dtype='int32'))
        self._cols.append(np.array([j], dtype='int32'))
        self._vals.append(np.array([value], dtype='float64'))

    def to_csc(self):
        """Convert accumulated entries to CSC sparse matrix."""
        if not self._rows:
            return csc_matrix((self._ndof, self._ndof))
        rows = np.concatenate(self._rows)
        cols = np.concatenate(self._cols)
        vals = np.concatenate(self._vals)
        return coo_matrix((vals, (rows, cols)), shape=(self._ndof, self._ndof)).tocsc()


def build_Kgg(
    model: BDF,
    dof_map: DOF_MAP,
    ndof: int,
    ngrid: int,
    ndof_per_grid: int,
    idtype: str = "int32",
    fdtype: str = "float32",) -> tuple[NDArrayNNfloat, Any]:
    """[K] = d{P}/dx"""
    model.log.debug("starting build_Kgg")
    Kbb = dok_matrix((ndof, ndof), dtype=fdtype)

    out = model.get_xyz_in_coord_array(cid=0, fdtype=fdtype, idtype=idtype)
    nid_cp_cd, xyz_cid0, xyz_cp, unused_icd_transform, unused_icp_transform = out
    all_nids = nid_cp_cd[:, 0]
    del xyz_cp, nid_cp_cd
    all_nids

    nelements = 0
    xyz_cid0 = model.grid.xyz_cid0()
    #for name in ['CELAS1', 'CELAS2', 'CELAS3', 'CELAS4']:
    #    nelements += _build_kbb_celas(model, name, Kbb, dof_map)
    nelements += _build_kbb_celas1(model, Kbb, dof_map)
    nelements += _build_kbb_celas2(model, Kbb, dof_map)
    nelements += _build_kbb_celas3(model, Kbb, dof_map)
    nelements += _build_kbb_celas4(model, Kbb, dof_map)

    for name in ['CROD', 'CTUBE', 'CONROD']:
        nelements += _build_kbb_rod(model, name, Kbb, dof_map, xyz_cid0)

    # Beam elements use COO batch assembly for performance
    coo = _COOAccumulator(ndof)
    nelements += _build_kbb_cbar_coo(model, coo, dof_map)
    nelements += _build_kbb_cbeam_coo(model, coo, dof_map)

    nelements += build_kbb_cquad4(
        model, Kbb, dof_map, all_nids, xyz_cid0, idtype="int32", fdtype="float64")
    nelements += build_kbb_ctria3(
        model, Kbb, dof_map, all_nids, xyz_cid0, idtype="int32", fdtype="float64")
    #nelements += _build_kbb_cbush(model, Kbb, dof_map)
    ndof = Kbb.shape[0]
    nelements += assemble_cbush_matrices(
        model, dof_map,
        Kbb, Mgg=None,
        include_damping=False, include_mass=False,
        extra_mass=None)
    nelements += _build_kbb_cshear(model, Kbb, dof_map, xyz_cid0)

    # Solid elements (CHEXA, CTETRA, CPENTA)
    nelements += build_kbb_chexa(model, coo, dof_map, xyz_cid0)
    nelements += build_kbb_ctetra(model, coo, dof_map, xyz_cid0)
    nelements += build_kbb_cpenta(model, coo, dof_map, xyz_cid0)

    assert nelements > 0, [elem for elem in model.element_cards if elem.n]

    # Merge dok (springs/rods/shells) + COO (beams/solids) into single CSC
    Kbb_csc = Kbb.tocsc()
    Kbb_beams = coo.to_csc()
    Kbb2 = Kbb_csc + Kbb_beams

    Kgg = Kbb_to_Kgg(model, Kbb2, ngrid, ndof_per_grid)

    model.log.debug("end of build_Kgg")
    return Kgg


def _build_kbb_cbar(model: BDF, Kbb: dok_matrix, dof_map: DOF_MAP,
                    fdtype: str = "float64") -> int:
    """Fill the CBAR Kbb matrix using a Timoshenko beam."""
    elem: CBAR = model.cbar
    nelements = elem.n
    if nelements == 0:
        return nelements

    LAIJEG = elem.stiffness_info()
    L = LAIJEG[:, 0]
    A = LAIJEG[:, 1]
    I = LAIJEG[:, [2, 3, 4]]
    J = LAIJEG[:, 5]

    k1 = LAIJEG[:, 6]
    k2 = LAIJEG[:, 7]

    E = LAIJEG[:, 8]
    G = LAIJEG[:, 9]
    assert LAIJEG.shape[1] == 10, LAIJEG.shape
    Nu = 0.3

    xyz1, xyz2 = elem.get_xyz()
    v, ihat, yhat, zhat, wa, wb = elem.get_axes(xyz1, xyz2)

    Teb = beam_transforms(ihat, yhat, zhat)
    for (eid, (nid1, nid2), areai, i1, i2, i12, j, pa, pb,
         Tebi, lengthi, k1i, k2i, e, g, nu) in zip(
         elem.element_id, elem.nodes,
         L, A, I, J, k1, k2, E, G, Nu, Teb, elem.pa, elem.pb):

        Ke = timoshenko_stiffness(
            areai, e, g, lengthi, i1, i2, j,
            k1i, k2i, pa=pa, pb=pb)
        Kebb = Teb.T @ Ke @ Teb

        gi1 = dof_map[(nid1, 1)]
        gi2 = dof_map[(nid2, 1)]
        n_ijv = [
            gi1, gi1 + 1, gi1 + 2, gi1 + 3, gi1 + 4, gi1 + 5,
            gi2, gi2 + 1, gi2 + 2, gi2 + 3, gi2 + 4, gi2 + 5,
        ]
        for i, ii in enumerate(n_ijv):
            for j, jj in enumerate(n_ijv):
                kij = Kebb[i, j]
                if kij != 0.0:
                    Kbb[ii, jj] += kij
    return nelements


def _build_kbb_cbar_coo(model: BDF, coo: _COOAccumulator,
                        dof_map: DOF_MAP, fdtype: str = "float64") -> int:
    """Fill CBAR stiffness using COO batch assembly."""
    elem = model.cbar
    nelements = elem.n
    if nelements == 0:
        return 0

    LAIJEG = elem.stiffness_info()
    L = LAIJEG[:, 0]
    A = LAIJEG[:, 1]
    I = LAIJEG[:, [2, 3, 4]]
    J = LAIJEG[:, 5]

    k1 = LAIJEG[:, 6]
    k2 = LAIJEG[:, 7]

    E = LAIJEG[:, 8]
    G = LAIJEG[:, 9]
    #A = elem.area()
    #inertia = elem.inertia()
    xyz1, xyz2 = elem.get_xyz()
    length = np.linalg.norm(xyz2 - xyz1, axis=1)
    v, ihat, yhat, zhat, wa, wb = elem.get_axes(xyz1, xyz2)
    k = elem.k()
    e_g_nus = elem.e_g_nu()

    Teb = beam_transforms(ihat, yhat, zhat)
    for ((nid1, nid2), Li, Ai, inertiai, ji, pa, pb,
          Tebi, k1i, k2i, e_g_nu) in zip(
          elem.nodes, L, A, I, J, elem.pa, elem.pb,
          Teb, k1, k2, e_g_nus,):
        i1, i2, i12 = inertiai
        e, g, nu = e_g_nu

        Ke = timoshenko_stiffness(Ai, e, g, Li,
                                  i1, i2, ji, k1i, k2i, pa=pa, pb=pb)
        Kebb = Tebi.T @ Ke @ Tebi

        gi1 = dof_map[(nid1, 1)]
        gi2 = dof_map[(nid2, 1)]
        n_ijv = [
            gi1, gi1 + 1, gi1 + 2, gi1 + 3, gi1 + 4, gi1 + 5,
            gi2, gi2 + 1, gi2 + 2, gi2 + 3, gi2 + 4, gi2 + 5,
        ]
        coo.add_matrix(n_ijv, Kebb)
    return nelements


def _build_kbb_cbeam_coo(model: BDF,
                         coo: _COOAccumulator,
                         dof_map: DOF_MAP,
                         fdtype: str = "float64") -> int:
    """Fill CBEAM stiffness using COO batch assembly."""
    elem = model.cbeam
    nelements = len(elem)
    if nelements == 0:
        return 0

    A = elem.area()
    I = elem.inertia()
    xyz1, xyz2 = elem.get_xyz()
    L = np.linalg.norm(xyz2 - xyz1, axis=1)
    v, ihat, yhat, zhat, wa, wb = elem.get_axes(xyz1, xyz2)
    k = elem.k()
    e_g_nus = elem.e_g_nu()

    Teb = beam_transforms(ihat, yhat, zhat)
    for ((nid1, nid2), Li, Ai, inertiai, pa, pb,
            Tebi, (k1i, k2i), e_g_nu) in zip(
            elem.nodes, L, A, I, elem.pa, elem.pb,
            Teb, k, e_g_nus,):

        i1, i2, i12, j = inertiai
        e, g, nu = e_g_nu

        Ke = timoshenko_stiffness(
            Ai, e, g, Li, i1, i2, j, k1i, k2i,
            pa=pa, pb=pb)
        Kebb = Tebi.T @ Ke @ Tebi

        gi1 = dof_map[(nid1, 1)]
        gi2 = dof_map[(nid2, 1)]
        n_ijv = [
            gi1, gi1 + 1, gi1 + 2, gi1 + 3, gi1 + 4, gi1 + 5,
            gi2, gi2 + 1, gi2 + 2, gi2 + 3, gi2 + 4, gi2 + 5,
        ]
        coo.add_matrix(n_ijv, Kebb)
    return nelements


def ke_cbar(
    model: BDF,
    xyz1: np.ndarray,
    xyz2: np.ndarray,
    length: float,
    Teb: np.ndarray,
    area: float,
    i1: float,
    i2: float,
    i12: float,
    j: float,
    pa: int,
    pb: int,
    k1: float,
    k2: float,
    E: float,
    G: float,
    fdtype: str = "float64",) -> tuple[bool, np.ndarray]:
    """Get the elemental stiffness matrix in the basic frame."""
    assert length > 0, length
    Ke = timoshenko_stiffness(area, E, G, length, i1, i2, j, k1, k2, pa, pb)
    K = Teb.T @ Ke @ Teb
    return True, K


def _build_kbb_cbeam(
    model: BDF,
    Kbb: dok_matrix,
    dof_map: DOF_MAP,
    all_nids: np.ndarray,
    xyz_cid0: np.ndarray,
    idtype: str = "int32",
    fdtype: str = "float64",) -> int:
    """Fill the CBEAM Kbb matrix using a Timoshenko beam."""
    elem = model.cbeam
    nelements = len(elem)
    if nelements == 0:
        return nelements

    LAIJEG = elem.stiffness_info()
    L = LAIJEG[:, 0]
    A = LAIJEG[:, 1]
    I = LAIJEG[:, [2, 3, 4]]
    J = LAIJEG[:, 5]

    k1 = LAIJEG[:, 6]
    k2 = LAIJEG[:, 7]

    xyz1, xyz2 = elem.get_xyz()
    assert len(length) == nelements
    v, ihat, yhat, zhat, wa, wb = elem.get_axes(xyz1, xyz2)

    e_g_nus = elem.e_g_nu()

    Teb = beam_transforms(ihat, yhat, zhat)
    for (eid, (nid1, nid2), Li, Ai, Ii, pa, pb, Tebi,
         k1i, k2i, (e, g, nu)) in zip(
            elem.element_id,
            elem.nodes,
            L, A, inertia,
            elem.pa, elem.pb,
            Teb, k1, k2, e_g_nus):
        i1i, i2i, i12i, ji = Ii

        Ke = timoshenko_stiffness(
            Ai, e, g, Li, i1i, i2i, ji,
            k1=k1i, k2=k2i, pa=pa, pb=pb)
        K = Teb.T @ Ke @ Teb

        gi1 = dof_map[(nid1, 1)]
        gi2 = dof_map[(nid2, 1)]
        n_ijv = [
            gi1, gi1 + 1, gi1 + 2, gi1 + 3, gi1 + 4, gi1 + 5,
            gi2, gi2 + 1, gi2 + 2, gi2 + 3, gi2 + 4, gi2 + 5,
        ]
        for i, ii in enumerate(n_ijv):
            for j, jj in enumerate(n_ijv):
                kij = K[i, j]
                if kij != 0.0:
                    Kbb[ii, jj] += kij
    return nelements


def _beami_stiffness(
    A: float,
    E: float,
    G: float,
    L: float,
    Iy: float,
    Iz: float,
    Iyz: float,
    J: float,
    pa: int = 0,
    pb: int = 0,
    k1: float = 1e8,
    k2: float = 1e8,):
    """Wrapper for backward compatibility; delegates to timoshenko_stiffness."""
    return timoshenko_stiffness(A, E, G, L, Iy, Iz, J, k1, k2, pa, pb)


def _ke_cshear(p1: np.ndarray,
               p2: np.ndarray,
               p3: np.ndarray,
               p4: np.ndarray,
               G: float, t: float,
               area: float,
               normal: np.ndarray,
               lxv: np.ndarray,
               lyv: np.ndarray,
               detJv: float,
               dN_dxv: float,
               dN_dyv: float) -> tuple[np.ndarray, float]:
    """Build the 12x12 CSHEAR element stiffness in global coordinates.

    Uses the isoparametric bilinear formulation with a single shear strain
    DOF (constant strain), consistent with the Garvey formulation in the
    NX Nastran Theoretical Manual Section 5.3.

    Parameters
    ----------
    p1, p2, p3, p4 : np.ndarray
        Node coordinates in global frame (3,).
    G : float
        Shear modulus.
    t : float
        Panel thickness.

    Returns
    -------
    Ke : np.ndarray
        12x12 element stiffness in global translational DOFs.
    area : float
        Element area.
    """
    if 0:  # pragma: no cover
        # Local coordinate system: x along edge 1-2, z = normal
        e12 = p2 - p1
        lx = e12 - np.dot(e12, normal) * normal
        lx_norm = np.linalg.norm(lx)
        lx /= lx_norm
        ly = np.cross(normal, lx)

        # Project nodes to local 2D
        pts = np.array([p1, p2, p3, p4])
        origin = p1
        x_local = np.dot(pts - origin, lx)
        y_local = np.dot(pts - origin, ly)

        # Bilinear shape function derivatives at element center (xi=0, eta=0)
        # Reference nodes: (-1,-1), (1,-1), (1,1), (-1,1)
        xi_n = np.array([-1., 1., 1., -1.])
        eta_n = np.array([-1., -1., 1., 1.])
        dN_dxi = xi_n / 4.0
        dN_deta = eta_n / 4.0

        # Jacobian at center
        dxdxi = np.dot(dN_dxi, x_local)
        dydxi = np.dot(dN_dxi, y_local)
        dxdeta = np.dot(dN_deta, x_local)
        dydeta = np.dot(dN_deta, y_local)
        detJ = dxdxi * dydeta - dydxi * dxdeta
        assert np.allclose(detJ, detJv)
        inv_detJ = 1.0 / detJ
        print('  dxdxi', dxdxi)
        print('  dxdeta', dxdeta)
        print('  detJ', detJ)
        print('  inv_detJ', inv_detJ)

        # Shape function derivatives in physical coordinates
        dN_dx = inv_detJ * (dydeta * dN_dxi - dydxi * dN_deta)  # correct
        #dN_dx = dydeta * dN_dxi - dydxi * dN_deta  # debugging
        dN_dy = inv_detJ * (-dxdeta * dN_dxi + dxdxi * dN_deta)
        #dN_dxi *= inv_detJ
        #dN_dyi *= inv_detJ
        print('  dN_dx', dN_dx)
        assert np.allclose(dN_dx, dN_dxv), (dN_dx, dN_dxv)
        assert np.allclose(dN_dy, dN_dyv), (dN_dy, dN_dyv)
    else:
        lx = lxv
        ly = lyv
        detJ = detJv
        dN_dx = dN_dxv #* np.random.random(dN_dxv.shape)
        dN_dy = dN_dyv #* np.random.random(dN_dxv.shape)

    # B_shear (1x8): gamma_xy = du/dy + dv/dx
    B = np.zeros(8)
    i = np.arange(4)
    B[2 * i] = dN_dy[i]
    B[2 * i + 1] = dN_dx[i]

    # 8x8 local stiffness: K_local = G*t * B^T*B * detJ * 4
    K_local_8 = (G * t * detJ * 4.0) * np.outer(B, B)

    # Transformation from local 2D (u_local, v_local) to global 3D (ux, uy, uz)
    # For each node: [u_local, v_local] = [[lx], [ly]]^T . [ux, uy, uz]
    # T_node = [[lx[0], lx[1], lx[2]],
    #           [ly[0], ly[1], ly[2]]]  0(2x3)
    T_node = np.array([lx, ly])  # (2, 3)

    # Full transformation: T (8x12) block-diagonal with T_node
    T = np.zeros((8, 12))
    for i in range(4):
        T[2 * i:2 * i + 2, 3 * i:3 * i + 3] = T_node

    # K_global (12x12) = T^T @ K_local_8 @ T
    Ke = T.T @ K_local_8 @ T
    return Ke


def _build_kbb_cshear(model: BDF,
                      Kbb: dok_matrix,
                      dof_map: DOF_MAP,
                      xyz_cid0: np.ndarray) -> int:
    """Fill CSHEAR stiffness using isoparametric bilinear shear formulation."""
    elem = model.cshear
    if elem.n == 0:
        return 0
    nelem = elem.n
    log = model.log

    grid = model.grid
    pshear = model.pshear
    mat1 = model.mat1
    
    eids = elem.element_id
    pids = elem.property_id
    ipids = pshear.index(pids)

    # [N1, N2, N3, N4]
    nidss = elem.nodes
    inids = grid.index(elem.nodes)
    ts = pshear.t[ipids]
    mids = pshear.material_id[ipids]
    imids = mat1.index(mids)
    Gs = mat1.G[imids]
    
    p1 = xyz_cid0[inids[:, 0], :]
    p2 = xyz_cid0[inids[:, 1], :]
    p3 = xyz_cid0[inids[:, 2], :]
    p4 = xyz_cid0[inids[:, 3], :]

    d13 = p3 - p1
    d24 = p4 - p2
    normals = np.cross(d13, d24, axis=1)
    assert normals.shape == p1.shape, f'normals.shape={normals.shape}; nelem={nelem}'
    normi = np.linalg.norm(normals, axis=1)
    areas = 0.5 * normi
    assert areas.shape == (nelem,), f'areas.shape={areas.shape} nelem={nelem}'
    if areas.min() <= 0.0:
        ifailed = (areas <= 0.0)
        eids_failed = eids[ifailed]
        raise FatalError(f'CSHEAR area must be greater than 0 for eids={eids}')
    normals /= normi[:, np.newaxis]

    # Local coordinate system: x along edge 1-2, z = normal
    e12 = p2 - p1
    e12_normal = np.einsum('ij,ij->i', e12, normals) # np.dot(e12, normal)
    lx = e12 - e12_normal[:, np.newaxis] * normals
    lx_norm = np.linalg.norm(lx, axis=1)
    assert lx_norm.shape == (nelem,)
    if lx_norm.min() <= 1e-14:
        ifailed = (lx_norm <= 1e-14)
        eids_failed = eids[ifailed]
        raise FatalError(f'Bad CSHEAR for eids={eids}')

    lx /= lx_norm[:, np.newaxis]
    ly = np.cross(normals, lx, axis=1)
    assert ly.shape == p1.shape

    #pts = np.array([p1, p2, p3, p4])
    #origin = p1
    #x_local = np.dot(pts - origin, lx)
    #y_local = np.dot(pts - origin, ly)
    dp1 = p1 - p1
    dp2 = p2 - p1
    dp3 = p3 - p1
    dp4 = p4 - p1
    p1_xlocal = np.einsum('ij,ij->i', dp1, lx)
    p2_xlocal = np.einsum('ij,ij->i', dp2, lx)
    p3_xlocal = np.einsum('ij,ij->i', dp3, lx)
    p4_xlocal = np.einsum('ij,ij->i', dp4, lx)
    xlocals = np.column_stack([
        p1_xlocal, p2_xlocal, p3_xlocal, p4_xlocal])

    p1_ylocal = np.einsum('ij,ij->i', dp1, ly)
    p2_ylocal = np.einsum('ij,ij->i', dp2, ly)
    p3_ylocal = np.einsum('ij,ij->i', dp3, ly)
    p4_ylocal = np.einsum('ij,ij->i', dp4, ly)
    ylocals = np.column_stack([
        p1_ylocal, p2_ylocal, p3_ylocal, p4_ylocal])

    # Bilinear shape function derivatives at element center (xi=0, eta=0)
    # Reference nodes: (-1,-1), (1,-1), (1,1), (-1,1)
    xi_n = np.array([-1., 1., 1., -1.])
    eta_n = np.array([-1., -1., 1., 1.])
    dN_dxi = xi_n / 4.0
    dN_deta = eta_n / 4.0

    # Jacobian at center
    #dxdxi = np.dot(dN_dxi, x_local)
    #dydxi = np.dot(dN_dxi, y_local)
    dxdxi = np.einsum('f,nf->n', dN_dxi, xlocals)
    dydxi = np.einsum('f,nf->n', dN_dxi, ylocals)
    #print('dxdxi', dxdxi)
    assert dxdxi.shape == (nelem,), f'dxdxi.shape={dxdxi.shape}; nelem={nelem}'

    #dxdeta = np.dot(dN_deta, x_local)
    #dydeta = np.dot(dN_deta, y_local)
    dxdeta = np.einsum('f,nf->n', dN_deta, xlocals)
    dydeta = np.einsum('f,nf->n', dN_deta, ylocals)
    assert dydeta.shape == (nelem,), f'dydeta.shape={dydeta.shape}; nelem={nelem}'

    detJ = dxdxi * dydeta - dydxi * dxdeta
    assert detJ.shape == (nelem,), f'detJ.shape={detJ.shape}; nelem={nelem}'
    if detJ.min() <= 1e-30:
        ifailed = (detJ <= 1e-30)
        eids_failed = eids[ifailed]
        raise FatalError(f'Bad CSHEAR Jacobian for eids={eids}')
    inv_detJ = 1.0 / detJ
    #print(f'inv_detJ = {inv_detJ}')

    # Shape function derivatives in physical coordinates
    #dN_dx = inv_detJ * (dydeta * dN_dxi - dydxi * dN_deta)
    #dN_dy = inv_detJ * (-dxdeta * dN_dxi + dxdxi * dN_deta)

    # Shape function derivatives in physical coordinates
    #dN_dx = inv_detJ[:, np.newaxis] * (
    #    dydeta[:, np.newaxis] * dN_dxi[np.newaxis, :]
    #    - dydxi[:, np.newaxis] * dN_deta[np.newaxis, :])

    dN_dx = inv_detJ[:, np.newaxis] * (
        dydeta[:, np.newaxis] * dN_dxi[np.newaxis, :]
        - dydxi[:, np.newaxis] * dN_deta[np.newaxis, :])
    dN_dy = inv_detJ[:, np.newaxis] * (
        - dxdeta[:, np.newaxis] * dN_dxi[np.newaxis, :]
        + dxdxi[:, np.newaxis] * dN_deta[np.newaxis, :])
    #print(f'dN_dx = \n{dN_dx}')
    assert dN_dx.shape == (nelem,4), f'dN_dx.shape={dN_dx.shape}; nelem={nelem}'

    for eid, pid, ipid, nids, inid, t, mid, G, area, normal, lxi, lyi, detJi, dN_dxi, dN_dyi in zip(eids, pids, ipids, nidss, inids, ts, mids, Gs, areas, normals, lx, ly, detJ, dN_dx, dN_dy):
        if ipid >= pshear.n or pshear.property_id[ipid] != pid:
            log.error(f'invalid pid; skipping eid={eid}')
            continue

        p1 = xyz_cid0[inid[0], :]
        p2 = xyz_cid0[inid[1], :]
        p3 = xyz_cid0[inid[2], :]
        p4 = xyz_cid0[inid[3], :]
        Ke = _ke_cshear(p1, p2, p3, p4, G, t, area, normal,
            lxi, lyi, detJi, dN_dxi, dN_dyi)

        # Assemble 12x12 into global Kbb (T1, T2, T3 at each node)
        dof_indices = []
        for nid in nids:
            dof_indices.append(dof_map[(nid, 1)])

        for ia in range(4):
            for ib in range(4):
                for da in range(3):
                    for db in range(3):
                        val = Ke[3 * ia + da, 3 * ib + db]
                        if abs(val) > 0.0:
                            row = dof_indices[ia] + da
                            col = dof_indices[ib] + db
                            Kbb[row, col] += val
    return elem.n


def Kbb_to_Kgg(model: BDF,
               Kbb: np.ndarray | csc_matrix,
               ngrid: int, ndof_per_grid: int,
               inplace: bool=True) -> NDArrayNNfloat:
    """Transform Kbb (basic frame) to Kgg (displacement coordinate frame).

    For grids with CD != 0, rotates stiffness rows/columns from
    basic to the grid's displacement coordinate system.

    Vectorized: builds all T6 matrices at once, then applies per unique CD.
    """
    assert isinstance(Kbb, (np.ndarray, csc_matrix)), type(Kbb)
    #Kbb = tolil(Kbb)
    if not isinstance(Kbb, np.ndarray):
        Kbb = Kbb.tolil()

    ndof = Kbb.shape[0]
    assert ndof > 0, f"ngrid={ngrid} card_count={model.card_count}"

    Kgg = Kbb
    if not inplace:
        Kgg = copy.deepcopy(Kgg)

    grid = model.grid
    cd_array = grid.cd
    if np.all(cd_array == 0):
        return Kgg

    coord = model.coord
    unique_cds = np.unique(cd_array)
    unique_cds = unique_cds[unique_cds != 0]
    if len(unique_cds) == 0:
        return Kgg

    # Pre-build T6 for each unique CD
    T6_map: dict[int, np.ndarray] = {}
    for cd in unique_cds:
        cd_slice = coord.slice_card_by_id(np.array([cd]))
        R = np.column_stack([cd_slice.i[0], cd_slice.j[0], cd_slice.k[0]])
        T6 = np.zeros((6, 6))
        T6[0:3, 0:3] = R
        T6[3:6, 3:6] = R
        T6_map[cd] = T6

    # Find all grid indices with CD != 0
    icd_nonzero = np.where(cd_array != 0)[0]

    is_dense = isinstance(Kgg, np.ndarray)
    for i in icd_nonzero:
        T6 = T6_map[cd_array[i]]
        T6T = T6.T
        i1 = i * ndof_per_grid
        i2 = i1 + ndof_per_grid

        if is_dense:
            row_block = Kgg[i1:i2, :].copy()
            Kgg[i1:i2, :] = T6T @ row_block
            col_block = Kgg[:, i1:i2].copy()
            Kgg[:, i1:i2] = col_block @ T6
        else:
            row_block = Kgg[i1:i2, :].toarray()
            Kgg[i1:i2, :] = T6T @ row_block
            col_block = Kgg[:, i1:i2].toarray()
            Kgg[:, i1:i2] = col_block @ T6

    return Kgg
