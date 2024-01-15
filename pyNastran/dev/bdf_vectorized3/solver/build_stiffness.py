from __future__ import annotations
import copy
from itertools import count
from typing import Union, TYPE_CHECKING

import numpy as np
from scipy.sparse._dok import dok_matrix
from scipy.sparse.csc import csc_matrix

#from pyNastran.dev.solver.stiffness.shells import build_kbb_cquad4, build_kbb_cquad8
from pyNastran.dev.solver.utils import lambda1d, DOF_MAP
#from pyNastran.bdf.cards.elements.bars import get_bar_vector, get_bar_yz_transform
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.nptyping_interface import NDArrayNNfloat
    from pyNastran.dev.bdf_vectorized3.bdf import (
        BDF,
        #CELAS1, CELAS2, CELAS3, CELAS4,
        #CBAR, PBAR, PBARL, PBEAM, PBEAML, # , CBEAM
        #MAT1,
    )
    from pyNastran.dev.bdf_vectorized3.bdf_interface.bdf_attributes import (
        CELAS1, CELAS2,
        CBAR,
        CTUBE, PTUBE,
    )


def build_Kgg(model: BDF,
              dof_map: DOF_MAP,
              ndof: int,
              ngrid: int,
              ndof_per_grid: int,
              idtype: str='int32', fdtype: str='float32') -> tuple[NDArrayNNfloat, Any]:
    """[K] = d{P}/dx"""
    model.log.debug(f'starting build_Kgg')
    Kbb = dok_matrix((ndof, ndof), dtype=fdtype)
    #print(dof_map)

    #_get_loadid_ndof(model, subcase_id)
    out = model.get_xyz_in_coord_array(cid=0, fdtype=fdtype, idtype=idtype)
    nid_cp_cd, xyz_cid0, xyz_cp, unused_icd_transform, unused_icp_transform = out
    all_nids = nid_cp_cd[:, 0]
    del xyz_cp, nid_cp_cd
    all_nids

    nelements = 0
    # TODO: I think these need to be in the global frame
    nelements += _build_kbb_celas1(model, Kbb, dof_map)
    nelements += _build_kbb_celas2(model, Kbb, dof_map)
    nelements += _build_kbb_celas3(model, Kbb, dof_map)
    nelements += _build_kbb_celas4(model, Kbb, dof_map)

    nelements += _build_kbb_conrod(model, Kbb, dof_map)
    nelements += _build_kbb_crod(model, Kbb, dof_map)
    nelements += _build_kbb_ctube(model, Kbb, dof_map)
    nelements += _build_kbb_cbar(model, Kbb, dof_map)
    nelements += _build_kbb_cbeam(model, Kbb, dof_map,
                                  all_nids, xyz_cid0, idtype='int32', fdtype='float64')
    #nelements += build_kbb_cquad4(model, Kbb, dof_map,
                                  #all_nids, xyz_cid0, idtype='int32', fdtype='float64')
    #nelements += build_kbb_cquad8(model, Kbb, dof_map,
                                  #all_nids, xyz_cid0, idtype='int32', fdtype='float64')
    assert nelements > 0, [elem for elem in model.element_cards if elem.n]
    Kbb2 = Kbb.tocsc()

    #Kgg = Kbb_to_Kgg(model, Kbb, ngrid, ndof_per_grid, inplace=False)
    Kgg = Kbb_to_Kgg(model, Kbb2, ngrid, ndof_per_grid)

    model.log.debug(f'end of build_Kgg')
    return Kgg


def _build_kbb_celas1(model: BDF,
                      Kbb: dok_matrix,
                      dof_map: DOF_MAP) -> None:
    """fill the CELAS1 Kbb matrix"""
    celas = model.celas1
    pelas = model.pelas
    if celas.n == 0:
        return celas.n
    pelas = pelas.slice_card_by_id(celas.property_id, assume_sorted=True)
    nids1 = celas.nodes[:, 0]
    nids2 = celas.nodes[:, 1]
    c1s = celas.components[:, 0]
    c2s = celas.components[:, 1]
    #pelas = model.celas2
    ks = pelas.k
    for nid1, nid2, c1, c2, ki in zip(nids1, nids2, c1s, c2s, ks):
        i = dof_map[(nid1, c1)]
        j = dof_map[(nid2, c2)]
        k = ki * np.array([[1, -1,],
                           [-1, 1]])
        ibe = [
            (i, 0),
            (j, 1),
        ]
        for ib1, ie1 in ibe:
            for ib2, ie2 in ibe:
                Kbb[ib1, ib2] += k[ie1, ie2]
    return celas.n

def _build_kbb_celas2(model: BDF,
                      Kbb: dok_matrix,
                      dof_map: DOF_MAP) -> int:
    """fill the CELAS2 Kbb matrix"""
    celas = model.celas2
    if celas.n == 0:
        return celas.n
    pelas = celas
    nids1 = celas.nodes[:, 0]
    nids2 = celas.nodes[:, 1]
    c1s = celas.components[:, 0]
    c2s = celas.components[:, 1]
    pelas = model.celas2
    ks = pelas.k
    k_unscaled = np.array([[1, -1,],
                           [-1, 1]])
    for nid1, nid2, c1, c2, ki in zip(nids1, nids2, c1s, c2s, ks):
        i = dof_map[(nid1, c1)]
        j = dof_map[(nid2, c2)]
        k = ki * k_unscaled
        ibe = [
            (i, 0),
            (j, 1),
        ]
        for ib1, ie1 in ibe:
            for ib2, ie2 in ibe:
                Kbb[ib1, ib2] += k[ie1, ie2]
    return celas.n

def _build_kbb_celas3(model: BDF,
                      Kbb: dok_matrix,
                      dof_map: DOF_MAP) -> None:
    """fill the CELAS3 Kbb matrix"""
    element = model.celas3
    nelement = len(element)
    if nelement == 0:
        return nelement
    eids = element.element_id
    pelas = model.pelas.slice_card_by_id(element.property_id, assume_sorted=True)
    nids1 = element.spoints[:, 0]
    nids2 = element.spoints[:, 1]

    ks = pelas.k
    Ke = np.full((nelement, 2, 2), np.nan, dtype='float64')
    for ielement, nid1, nid2, ki in zip(count(), nids1, nids2, ks):
        #i = dof_map[(nid1, c1)]
        #j = dof_map[(nid2, c2)]
    #for eid in eids:
        #elem = model.elements[eid]
        #ki = elem.K()
        #print(elem, ki)
        #print(elem.get_stats())
        ke = _build_kbbi_celas34(Kbb, dof_map, nid1, nid2, ki)
        Ke[ielement, :, :] = ke
    return nelement

def _build_kbb_celas4(model: BDF,
                      Kbb: dok_matrix,
                      dof_map: DOF_MAP) -> int:
    """fill the CELAS4 Kbb matrix"""
    element = model.celas4
    nelement = len(element)
    if nelement == 0:
        return nelement

    eids = element.element_id
    ks = element.k
    nids1 = element.spoints[:, 0]
    nids2 = element.spoints[:, 1]

    Ke = np.full((nelement, 2, 2), np.nan, dtype='float64')
    for ielement, eid, nid1, nid2, ki in zip(count(), eids, nids1, nids2, ks):
        ke = _build_kbbi_celas34(Kbb, dof_map, nid1, nid2, ki)
        Ke[ielement, :, :] = ke
    return nelement

def _build_kbbi_celas12(Kbb: dok_matrix,
                        dof_map: DOF_MAP,
                        elem: Union[CELAS1, CELAS2],
                        ki: float) -> np.ndarray:
    """fill the CELASx Kbb matrix"""
    nid1, nid2 = elem.nodes
    c1, c2 = elem.c1, elem.c2
    i = dof_map[(nid1, c1)]
    j = dof_map[(nid2, c2)]
    ke = ki * np.array([[1, -1,],
                        [-1, 1]])
    ibe = [
        (i, 0),
        (j, 1),
    ]
    for ib1, ie1 in ibe:
        for ib2, ie2 in ibe:
            Kbb[ib1, ib2] += ke[ie1, ie2]
    #Kbb[j, i] += ki
    #Kbb[i, j] += ki
    #del i, j, ki, nid1, nid2, c1, c2
    return ke

def _build_kbbi_celas34(Kbb: dok_matrix,
                        dof_map: DOF_MAP,
                        nid1: int, nid2: int,
                        ki: float) -> np.ndarray:
    """fill the CELASx Kbb matrix"""
    #print(dof_map)
    i = dof_map[(nid1, 0)]
    j = dof_map[(nid2, 0)]
    ke = ki * np.array([[1, -1,],
                        [-1, 1]])
    ibe = [
        (i, 0),
        (j, 1),
    ]
    for ib1, ie1 in ibe:
        for ib2, ie2 in ibe:
            Kbb[ib1, ib2] += ke[ie1, ie2]
    #Kbb[j, i] += ki
    #Kbb[i, j] += ki
    #del i, j, ki, nid1, nid2, c1, c2
    return ke

def _build_kbb_cbar(model: BDF,
                    Kbb: dok_matrix,
                    dof_map: DOF_MAP,
                    fdtype: str='float64') -> int:
    """fill the CBAR Kbb matrix using an Euler-Bernoulli beam"""
    elem: CBAR = model.cbar
    nelements = elem.n
    if nelements == 0:
        return nelements

    #area = np.zeros(nelements, dtype='float64')
    #i = np.zeros((nelements, 3), dtype='float64')
    #j = np.zeros(nelements, dtype='float64')
    #nsm = np.zeros(nelements, dtype='float64')

    area = elem.area()
    inertia = elem.inertia()
    xyz1, xyz2 = elem.get_xyz()
    length = np.linalg.norm(xyz2 - xyz1, axis=1)
    assert len(length) == nelements
    v, ihat, yhat, zhat, wa, wb = elem.get_axes(xyz1, xyz2)
    k = elem.k()
    e_g_nus = elem.e_g_nu()

    Ke = np.full((nelements, 12, 12), np.nan, dtype=fdtype)
    for ielement, eid, (nid1, nid2), areai, inertiai, \
        pa, pb, ihati, jhati, khati, lengthi, \
        ki, e_g_nu in zip(count(),
                          elem.element_id, elem.nodes,
                          area, inertia,
                          elem.pa, elem.pb,
                          ihat, yhat, zhat, length,
                          k, e_g_nus):
        i1, i2, i12, j = inertiai
        e, g, nu = e_g_nu
        #nid1, nid2 = elem.nodes
        #is_passed, K = ke_cbar(model, elem, fdtype=fdtype)
        k1, k2 = ki
        is_passed, K = ke_cbar(model, xyz1, xyz2, lengthi,
                               ihati, jhati, khati,
                               areai, i1, i2, i12, j,
                               pa, pb, k1, k2, e, g,)
        assert is_passed

        i1 = dof_map[(nid1, 1)]
        i2 = dof_map[(nid2, 1)]
        n_ijv = [
            i1,
            i1 + 1, i1 + 2,
            i1 + 3, i1 + 4, i1 + 5,

            i2,
            i2 + 1, i2 + 2,
            i2 + 3, i2 + 4, i2 + 5,
        ]
        for i, i1 in enumerate(n_ijv):
            for j, i2 in enumerate(n_ijv):
                ki = K[i, j]
                if abs(ki) > 0.:
                    Kbb[i1, i2] += ki
        Ke[ielement, :, :] = K
    return nelements

def ke_cbar(model: BDF,
            #v, ihat, jhat, khat, wa, wb
            xyz1: np.ndarray, xyz2: np.ndarray, length: float,
            ihat: np.ndarray, jhat: np.ndarray, khat: np.ndarray,
            area: float,
            i1: float, i2: float, i12: float, j: float,
            pa: int, pb: int,
            k1: float, k2: float,
            E: float, G: float,
            fdtype: str='float64') -> tuple[bool, np.ndarray]:
    """get the elemental stiffness matrix in the basic frame"""
    #pid_ref = elem.pid_ref
    #mat = pid_ref.mid_ref

    #is_passed, (v, ihat, jhat, khat, wa, wb) = elem.get_axes(model)
    #T = np.vstack([ihat, jhat, khat])
    #z = np.zeros((3, 3), dtype='float64')
    #prop = elem.pid_ref
    #mat = prop.mid_ref
    #I1 = prop.I11()
    #I2 = prop.I22()
    #unused_I12 = prop.I12()
    #pa = elem.pa
    #pb = elem.pb
    #J = prop.J()
    #E = mat.E()
    #G = mat.G()
    z = np.zeros((3, 3), dtype='float64')
    T = z
    #unused_Teb = np.block([
        #[T, z],
        #[z, T],
    #])
    #is_failed, (v, ihat, jhat, khat, wa, wb) = elem.get_axes(model)
    #assert is_failed is False
    #print(wa, wb)
    #xyz1 = elem.nodes_ref[0].get_position() + wa
    #xyz2 = elem.nodes_ref[1].get_position() + wb
    dxyz = xyz2 - xyz1
    #L = np.linalg.norm(dxyz)
    assert length > 0, length
    #pid_ref = elem.pid_ref
    T = np.vstack([ihat, jhat, khat])
    z = np.zeros((3, 3), dtype=fdtype)
    Teb = np.block([
        [T, z, z, z],
        [z, T, z, z],
        [z, z, T, z],
        [z, z, z, T],
    ])
    Ke = _beami_stiffness(
        area, E, G, length,
        i1, i2, i12, j,
        k1=k1, k2=k2,
        pa=pa, pb=pb)
    K = Teb.T @ Ke @ Teb
    is_passed = True
    return is_passed, K

def _build_kbb_crod(model: BDF,
                    Kbb: dok_matrix,
                    dof_map: DOF_MAP) -> int:
    """fill the CROD Kbb matrix"""
    elem = model.crod
    if elem.n == 0:
        return elem.n
    prop = model.prod
    mat = model.mat1

    xyz_cid0 = model.grid.xyz_cid0()
    nodes = elem.nodes
    pids = elem.property_id
    inid = model.grid.index(nodes)
    inid1 = inid[:, 0]
    inid2 = inid[:, 1]
    xyz1 = xyz_cid0[inid1, :]
    xyz2 = xyz_cid0[inid2, :]
    dxyz = xyz2 - xyz1
    length = np.linalg.norm(dxyz, axis=1)
    assert len(length) == elem.n

    prop2 = prop.slice_card_by_id(pids, assume_sorted=True)
    area = prop2.area()
    J = prop2.J
    material_id = prop2.material_id
    mat1 = mat.slice_card_by_material_id(material_id)
    assert elem.n == prop.n
    assert elem.n == mat1.n

    G = mat1.G
    E = mat1.E
    L = length
    A = area
    k_axial = A * E / L
    k_torsion = G * J / L
    for nodes, dxyzi, k_axiali, k_torsioni in zip(elem.nodes, dxyz, k_axial, k_torsion):
        _build_kbbi_conrod_crod(Kbb, dof_map, nodes, dxyzi, k_axiali, k_torsioni)
    return len(elem)

def _build_kbb_ctube(model: BDF,
                     Kbb: dok_matrix,
                     dof_map: DOF_MAP) -> None:
    """fill the CTUBE Kbb matrix"""
    elem = model.ctube
    if elem.n == 0:
        return elem.n
    prop: PTUBE = model.ptube
    mat = model.mat1

    xyz_cid0 = model.grid.xyz_cid0()
    nodes = elem.nodes
    pids = elem.property_id
    inid = model.grid.index(nodes)
    inid1 = inid[:, 0]
    inid2 = inid[:, 1]
    xyz1 = xyz_cid0[inid1, :]
    xyz2 = xyz_cid0[inid2, :]
    dxyz = xyz2 - xyz1
    length = np.linalg.norm(dxyz, axis=1)
    assert len(length) == elem.n

    prop2: PTUBE = prop.slice_card_by_id(pids, assume_sorted=True)
    area = prop2.area()
    J = prop2.J()
    material_id = prop2.material_id
    mat1 = mat.slice_card_by_material_id(material_id)
    assert elem.n == prop.n
    assert elem.n == mat1.n

    G = mat1.G
    E = mat1.E
    L = length
    A = area
    k_axial = A * E / L
    k_torsion = G * J / L
    for nodes, dxyzi, k_axiali, k_torsioni in zip(elem.nodes, dxyz, k_axial, k_torsion):
        _build_kbbi_conrod_crod(Kbb, dof_map, nodes, dxyzi, k_axiali, k_torsioni)
    return len(elem)

def _build_kbb_conrod(model: BDF,
                      Kbb: dok_matrix,
                      dof_map: DOF_MAP) -> int:
    """fill the CONROD Kbb matrix"""
    elem = model.conrod
    if elem.n == 0:
        return elem.n
    prop = model.conrod
    mat = model.mat1

    xyz_cid0 = model.grid.xyz_cid0()
    nodes = elem.nodes
    inid = model.grid.index(nodes)
    inid1 = inid[:, 0]
    inid2 = inid[:, 1]
    xyz1 = xyz_cid0[inid1, :]
    xyz2 = xyz_cid0[inid2, :]
    dxyz = xyz2 - xyz1
    length = np.linalg.norm(dxyz, axis=1)
    assert len(length) == elem.n

    area = prop.area()
    material_id = prop.material_id
    mat1 = mat.slice_card_by_material_id(material_id)
    assert elem.n == prop.n
    assert elem.n == mat1.n
    J = prop.J
    G = mat1.G
    E = mat1.E
    L = length
    A = area
    k_axial = A * E / L
    k_torsion = G * J / L
    for nodes, dxyzi, k_axiali, k_torsioni in zip(elem.nodes, dxyz, k_axial, k_torsion):
        _build_kbbi_conrod_crod(Kbb, dof_map, nodes, dxyzi, k_axiali, k_torsioni)
    return len(elem)

def _build_kbbi_conrod_crod(Kbb: dok_matrix,
                            dof_map: DOF_MAP,
                            nodes,
                            dxyz12: np.ndarray,
                            k_axial: float,
                            k_torsion: float,
                            fdtype: str='float64') -> None:
    """fill the ith rod Kbb matrix"""
    nid1, nid2 = nodes
    #print(f'A = {A}')
    #k_axial = A * E / L
    #k_torsion = G * J / L

    assert isinstance(k_axial, float), k_axial
    assert isinstance(k_torsion, float), k_torsion
    k = np.array([[1., -1.],
                  [-1., 1.]])  # 1D rod
    Lambda = lambda1d(dxyz12, debug=False)
    K = Lambda.T @ k @ Lambda
    #i11 = dof_map[(n1, 1)]
    #i12 = dof_map[(n1, 2)]

    #i21 = dof_map[(n2, 1)]
    #i22 = dof_map[(n2, 2)]

    nki, nkj = K.shape
    K2 = np.zeros((nki*2, nkj*2), dtype=fdtype)

    # axial + torsion; assume 3D
    # u1fx, u1fy, u1fz, u2fx, u2fy, u2fz
    K2[:nki, :nki] = K * k_axial

    # u1mx, u1my, u1mz, u2mx, u2my, u2mz
    K2[nki:, nki:] = K * k_torsion

    i1 = 0
    i2 = 3 # dof_map[(n1, 2)]
    dofs = np.array([
        i1, i1+1, i1+2,
        i2, i2+1, i2+2,

        i1+3, i1+4, i1+5,
        i2+3, i2+4, i2+5,
    ], dtype='int32')

    ni1 = dof_map[(nid1, 1)]
    nj2 = dof_map[(nid2, 1)]
    n_ijv = [
        # axial
        ni1, ni1 + 1, ni1 + 2,  # node 1
        nj2, nj2 + 1, nj2 + 2,  # node 2

        # torsion
        ni1 + 3, ni1 + 4, ni1 + 5,  # node 1
        nj2 + 3, nj2 + 4, nj2 + 5,  # node 2
    ]
    for dof1, i1 in zip(dofs, n_ijv):
        for dof2, i2 in zip(dofs, n_ijv):
            ki = K2[dof1, dof2]
            if abs(ki) > 0.:
                #print(nij1, nij2, f'({i1}, {i2});', (dof1, dof2), ki)
                Kbb[i1, i2] += ki
        #print(K2)
    #print(Kbb)
    return

def _build_kbb_cbeam(model: BDF,
                     Kbb: dok_matrix,
                     dof_map: DOF_MAP,
                     all_nids: np.ndarray,
                     xyz_cid0: np.ndarray,
                     idtype: str='int32',
                     fdtype: str='float64') -> int:
    """TODO: Timoshenko beam, warping, I12"""
    str(all_nids)
    str(xyz_cid0)
    elem = model.cbeam
    nelements = len(elem)
    if nelements == 0:
        return nelements

    area = elem.area()
    inertia = elem.inertia()
    xyz1, xyz2 = elem.get_xyz()
    length = np.linalg.norm(xyz2 - xyz1, axis=1)
    assert len(length) == nelements
    v, ihat, yhat, zhat, wa, wb = elem.get_axes(xyz1, xyz2)
    k = elem.k()
    e_g_nus = elem.e_g_nu()

    Ke = np.full((nelements, 12, 12), np.nan, dtype=fdtype)
    for ielement, eid, (nid1, nid2), areai, inertiai, \
        pa, pb, ihati, jhati, khati, lengthi, \
        ki, e_g_nu in zip(count(), elem.element_id, elem.nodes,
                          area, inertia,
                          elem.pa, elem.pb,
                          ihat, yhat, zhat, length,
                          k, e_g_nus):
        i1, i2, i12, j = inertiai
        e, g, nu = e_g_nu
        #nid1, nid2 = elem.nodes
        #is_passed, K = ke_cbar(model, elem, fdtype=fdtype)
        k1, k2 = ki

        dxyz = xyz2 - xyz1
        L = np.linalg.norm(dxyz)
        assert L > 0, L
        is_passed, K = ke_cbar(model, xyz1, xyz2, lengthi,
                               ihati, jhati, khati,
                               areai, i1, i2, i12, j,
                               pa, pb, k1, k2, e, g,)
        #is_failed, (v, ihat, jhat, khat, wa, wb) = elem.get_axes(model)
        #print(wa, wb, ihat, jhat, khat)
        #assert is_failed is False

        T = np.vstack([ihati, jhati, khati])
        z = np.zeros((3, 3), dtype=fdtype)
        Teb = np.block([
            [T, z, z, z],
            [z, T, z, z],
            [z, z, T, z],
            [z, z, z, T],
        ])
        Kei = _beami_stiffness(
            areai, e, g, lengthi,
            i1, i2, i12, j,
            pa, pb, k1, k2)
        K = Teb.T @ Kei @ Teb

        i1 = dof_map[(nid1, 1)]
        j1 = dof_map[(nid2, 1)]
        n_ijv = [
            i1, i1 + 1, i1 + 2, i1 + 3, i1 + 4, i1 + 5, # node 1
            j1, j1 + 1, j1 + 2, j1 + 3, j1 + 4, j1 + 5, # node 2
        ]
        for i, i1 in enumerate(n_ijv):
            for j, i2 in enumerate(n_ijv):
                ki = K[i, j]
                if abs(ki) > 0.:
                    Kbb[i1, i2] += ki
        Ke[ielement, :, :] = Kei
    return nelements

def _beami_stiffness(A: float, E: float, G: float, L: float,
                     Iy: float, Iz: float, Iyz: float, J: float,
                     pa: int, pb: int,
                     k1: float, k2: float):
    """gets the ith Euler-Bernoulli beam stiffness"""
    kaxial = E * A / L
    ktorsion = G * J / L

    L2 = L * L
    if k1 is not None:
        phiy = 12. * E * Iy / (k1 * G * A * L2)
    if k2 is not None:
        phiz = 12. * E * Iy / (k2 * G * A * L2)

    phiy = 1.0
    phiz = 1.0
    ky = E * Iy / (L * phiy)
    kz = E * Iz / (L * phiz)

    K = np.zeros((12, 12), dtype='float64')
    # axial
    K[0, 0] = K[6, 6] = kaxial
    K[6, 0] = K[0, 6] = -kaxial

    # torsion
    K[3, 3] = K[9, 9] = ktorsion
    K[9, 3] = K[3, 9] = -ktorsion

    #Fx - 0, 6
    #Fy - 1, 7**
    #Fz - 2, 8
    #Mx - 3, 9
    #My - 4, 10
    #Mz - 5, 11**
    # bending z (Fy/1/7, Mz/5/11)
    #      1     5       7   11
    # 1  [12  & 6L   & -12 & 6L
    # 5  [6L  & 4L^2 & -6L & 2L^2
    # 7  [-12 &-6L   &  12 & -6L
    # 11 [6L  & 2L^2 & -6L & 4L^2
    K[1, 1] = K[7, 7] = 12. * kz
    K[1, 7] = K[1, 7] = -12. * kz
    K[1, 5] = K[5, 1] = K[11, 1] = K[1, 11] = 6. * L * kz

    K[5, 7] = K[7, 5] = K[7, 11] = K[11, 7] = -6. * L * kz
    K[5, 11] = K[11, 5] = 2. * L2 * kz * (2 - phiz)
    K[5, 5] = K[11, 11] = 4. * L2 * kz * (4 + phiz)

    #Fx - 0, 6
    #Fy - 1, 7
    #Fz - 2, 8**
    #Mx - 3, 9
    #My - 4, 10**
    #Mz - 5, 11
    # bending y (Fz/2/8, My/4/10)
    #      2     4       8   10
    # 2  [12  & 6L   & -12 & 6L
    # 4  [6L  & 4L^2 & -6L & 2L^2
    # 8  [-12 &-6L   &  12 & -6L
    # 10 [6L  & 2L^2 & -6L & 4L^2
    K[2, 2] = K[8, 8] = 12. * ky
    K[2, 8] = K[2, 8] = -12. * ky
    K[2, 4] = K[4, 2] = K[10, 2] = K[2, 10] = 6. * L * ky

    K[4, 8] = K[8, 4] = K[8, 10] = K[10, 8] = -6. * L * ky
    K[4, 10] = K[10, 4] = 2. * L * L * ky * (2. - phiy)
    K[4, 4] = K[10, 10] = 4. * L * L * ky * (4. + phiy)

    if pa != 0:
        assert pa > 0
        for pas in str(pa): # 123456
            pai = int(pas) - 1 # 012345 (0-5)
            K[pai, :] = 0
            K[:, pai] = 0
    if pb != 0:
        assert pb > 0
        for pbs in str(pb): # 123456
            pbi = int(pbs) + 5 #  (6 - 11)
            K[pbi, :] = 0
            K[:, pbi] = 0
    return K

def Kbb_to_Kgg(model: BDF,
               Kbb: Union[np.ndarray, csc_matrix],
               ngrid: int,
               ndof_per_grid: int,
               inplace=True) -> NDArrayNNfloat:
    """does an in-place transformation"""
    assert isinstance(Kbb, (np.ndarray, csc_matrix)), type(Kbb)
    #assert isinstance(Kbb, (np.ndarray, csc_matrix, dok_matrix)), type(Kbb)
    if not isinstance(Kbb, np.ndarray):
        Kbb = Kbb.tolil()

    ndof = Kbb.shape[0]
    assert ndof > 0, f'ngrid={ngrid} card_count={model.card_count}'
    nids = model._type_to_id_map['GRID']

    Kgg = Kbb
    if not inplace:
        Kgg = copy.deepcopy(Kgg)

    for i, nid in enumerate(nids):
        node = model.nodes[nid]
        if node.cd:
            model.log.debug(f'node {nid} has a CD={node.cd}')
            cd_ref = node.cd_ref
            T = cd_ref.beta_n(n=2)
            i1 = i * ndof_per_grid
            i2 = (i+1) * ndof_per_grid
            Ki = Kbb[i1:i2, i1:i2]
            Kgg[i1:i2, i1:i2] = T.T @ Ki @ T
    return Kgg
