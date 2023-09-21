from __future__ import annotations
import copy
from typing import Union, Optional, Any, TYPE_CHECKING

import numpy as np
import scipy.sparse as sci_sparse

#from pyNastran.dev.solver.stiffness.shells import build_kbb_cquad4, build_kbb_cquad8
from .utils import lambda1d, DOF_MAP
#from pyNastran.bdf.cards.elements.bars import get_bar_vector, get_bar_yz_transform
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.nptyping_interface import NDArrayNNfloat
    from pyNastran.dev.bdf_vectorized3.bdf import (
        BDF,
        #CELAS1, CELAS2, CELAS3, CELAS4,
        #CBAR, PBAR, PBARL, PBEAM, PBEAML, # , CBEAM
        #MAT1,
    )


def build_Kgg(model: BDF, dof_map: DOF_MAP,
              ndof: int,
              ngrid: int,
              ndof_per_grid: int,
              idtype: str='int32', fdtype: str='float32') -> tuple[NDArrayNNfloat, Any]:
    """[K] = d{P}/dx"""
    model.log.debug(f'starting build_Kgg')
    Kbb = sci_sparse.dok_matrix((ndof, ndof), dtype=fdtype)
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
    assert nelements > 0, [elem for elem in model.elements if elem.n]
    Kbb2 = Kbb.tocsc()

    #Kgg = Kbb_to_Kgg(model, Kbb, ngrid, ndof_per_grid, inplace=False)
    Kgg = Kbb_to_Kgg(model, Kbb2, ngrid, ndof_per_grid)

    model.log.debug(f'end of build_Kgg')
    return Kgg


def _build_kbb_celas1(model: BDF, Kbb, dof_map: DOF_MAP) -> None:
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

def _build_kbb_celas2(model: BDF, Kbb, dof_map: DOF_MAP) -> int:
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

def _build_kbb_celas3(model: BDF, Kbb, dof_map: DOF_MAP) -> None:
    """fill the CELAS3 Kbb matrix"""
    celas = model.celas3
    if celas.n == 0:
        return celas.n
    eids = celas.element_id
    pelas = model.pelas.slice_card_by_id(celas.property_id, assume_sorted=True)
    nids1 = celas.spoints[:, 0]
    nids2 = celas.spoints[:, 1]

    ks = pelas.k
    for nid1, nid2, ki in zip(nids1, nids2, ks):
        #i = dof_map[(nid1, c1)]
        #j = dof_map[(nid2, c2)]
    #for eid in eids:
        #elem = model.elements[eid]
        #ki = elem.K()
        #print(elem, ki)
        #print(elem.get_stats())
        _build_kbbi_celas34(Kbb, dof_map, nid1, nid2, ki)
    return len(eids)

def _build_kbb_celas4(model: BDF, Kbb, dof_map: DOF_MAP) -> None:
    """fill the CELAS4 Kbb matrix"""
    celas = model.celas4
    if celas.n == 0:
        return celas.n
    eids = model._type_to_id_map['CELAS4']
    for eid in eids:
        elem = model.elements[eid]
        ki = elem.K()
        #print(elem, ki)
        #print(elem.get_stats())
        _build_kbbi_celas34(Kbb, dof_map, elem, ki)
    return len(eids)

def _build_kbbi_celas12(Kbb, dof_map: DOF_MAP,
                        elem: Union[CELAS1, CELAS2], ki: float) -> None:
    """fill the CELASx Kbb matrix"""
    nid1, nid2 = elem.nodes
    c1, c2 = elem.c1, elem.c2
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
    #Kbb[j, i] += ki
    #Kbb[i, j] += ki
    #del i, j, ki, nid1, nid2, c1, c2

def _build_kbbi_celas34(Kbb, dof_map: DOF_MAP,
                        nid1: int, nid2: int, ki: float) -> None:
    """fill the CELASx Kbb matrix"""
    #print(dof_map)
    i = dof_map[(nid1, 0)]
    j = dof_map[(nid2, 0)]
    k = ki * np.array([[1, -1,],
                       [-1, 1]])
    ibe = [
        (i, 0),
        (j, 1),
    ]
    for ib1, ie1 in ibe:
        for ib2, ie2 in ibe:
            Kbb[ib1, ib2] += k[ie1, ie2]
    #Kbb[j, i] += ki
    #Kbb[i, j] += ki
    #del i, j, ki, nid1, nid2, c1, c2

def _build_kbb_cbar(model, Kbb, dof_map: DOF_MAP, fdtype: str='float64') -> int:
    """fill the CBAR Kbb matrix using an Euler-Bernoulli beam"""
    eids = model._type_to_id_map['CBAR']
    nelements = len(eids)
    if nelements == 0:
        return nelements

    for eid in eids:
        elem = model.elements[eid]  # type: CBAR
        nid1, nid2 = elem.nodes
        is_passed, K = ke_cbar(model, elem, fdtype=fdtype)
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
    return nelements

def ke_cbar(model: BDF, elem: CBAR, fdtype: str='float64'):
    """get the elemental stiffness matrix in the basic frame"""
    pid_ref = elem.pid_ref
    mat = pid_ref.mid_ref

    #is_passed, (wa, wb, ihat, jhat, khat) = elem.get_axes(model)
    #T = np.vstack([ihat, jhat, khat])
    #z = np.zeros((3, 3), dtype='float64')
    prop = elem.pid_ref
    mat = prop.mid_ref
    I1 = prop.I11()
    I2 = prop.I22()
    unused_I12 = prop.I12()
    pa = elem.pa
    pb = elem.pb
    #J = prop.J()
    #E = mat.E()
    #G = mat.G()
    z = np.zeros((3, 3), dtype='float64')
    T = z
    #unused_Teb = np.block([
        #[T, z],
        #[z, T],
    #])
    is_failed, (wa, wb, ihat, jhat, khat) = elem.get_axes(model)
    assert is_failed is False
    #print(wa, wb)
    xyz1 = elem.nodes_ref[0].get_position() + wa
    xyz2 = elem.nodes_ref[1].get_position() + wb
    dxyz = xyz2 - xyz1
    L = np.linalg.norm(dxyz)
    pid_ref = elem.pid_ref
    mat = pid_ref.mid_ref
    T = np.vstack([ihat, jhat, khat])
    z = np.zeros((3, 3), dtype=fdtype)
    Teb = np.block([
        [T, z, z, z],
        [z, T, z, z],
        [z, z, T, z],
        [z, z, z, T],
    ])
    k1 = pid_ref.k1
    k2 = pid_ref.k2
    Ke = _beami_stiffness(pid_ref, mat, L, I1, I2, k1=k1, k2=k2, pa=pa, pb=pb)
    K = Teb.T @ Ke @ Teb
    is_passed = not is_failed
    return is_passed, K

def _build_kbb_crod(model: BDF, Kbb, dof_map: DOF_MAP) -> None:
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
    Js = prop2.J
    material_id = prop2.material_id
    mat1 = mat.slice_card_by_material_id(material_id)
    assert elem.n == prop.n
    assert elem.n == mat1.n
    for nodes, dxyzi, L, A, E, G, J in zip(elem.nodes, dxyz, length, area, mat1.E, mat1.G, Js):
        _build_kbbi_conrod_crod(Kbb, dof_map, nodes, dxyzi, L, A, E, G, J)
    return len(elem)

def _build_kbb_ctube(model: BDF, Kbb, dof_map: DOF_MAP) -> None:
    """fill the CTUBE Kbb matrix"""
    elem = model.ctube
    if elem.n == 0:
        return elem.n
    prop = model.ptube
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
    Js = area * 0.
    material_id = prop2.material_id
    mat1 = mat.slice_card_by_material_id(material_id)
    assert elem.n == prop.n
    assert elem.n == mat1.n
    for nodes, dxyzi, L, A, E, G, J in zip(elem.nodes, dxyz, length, area, mat1.E, mat1.G, Js):
        _build_kbbi_conrod_crod(Kbb, dof_map, nodes, dxyzi, L, A, E, G, J)
    return len(elem)

def _build_kbb_conrod(model: BDF, Kbb, dof_map: DOF_MAP) -> int:
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
    Js = prop.J
    material_id = prop.material_id
    mat1 = mat.slice_card_by_material_id(material_id)
    assert elem.n == prop.n
    assert elem.n == mat1.n
    for nodes, dxyzi, L, A, E, G, J in zip(elem.nodes, dxyz, length, area, mat1.E, mat1.G, Js):
        _build_kbbi_conrod_crod(Kbb, dof_map, nodes, dxyzi, L, A, E, G, J)
    return len(elem)

def _build_kbbi_conrod_crod(Kbb, dof_map: DOF_MAP, nodes,
                            dxyz12, L, A, E, G, J, fdtype='float64') -> None:
    """fill the ith rod Kbb matrix"""
    nid1, nid2 = nodes
    #mat = elem.mid_ref
    #xyz1 = elem.nodes_ref[0].get_position()
    #xyz2 = elem.nodes_ref[1].get_position()
    #dxyz12 = xyz1 - xyz2
    #L = np.linalg.norm(dxyz12)
    #E = mat.E
    #G = mat.G()
    #J = elem.J()
    #A = elem.Area()
    #print(f'A = {A}')
    #L = elem.Length()
    k_axial = A * E / L
    k_torsion = G * J / L

    assert isinstance(k_axial, float), k_axial
    assert isinstance(k_torsion, float), k_torsion
    #Kbb[i, i] += ki[0, 0]
    #Kbb[i, j] += ki[0, 1]
    #Kbb[j, i] = ki[1, 0]
    #Kbb[j, j] = ki[1, 1]
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

    ni1 = dof_map[(nid1, 1)]
    nj1 = dof_map[(nid2, 1)]

    i1 = 0
    i2 = 3 # dof_map[(n1, 2)]
    if k_torsion == 0.0: # axial; 2D or 3D
        K2 = K * k_axial
        n_ijv = [
            # axial
            ni1, ni1 + 1, ni1 + 2,  # node 1
            nj1, nj1 + 1, nj1 + 2,  # node 2
        ]
        dofs = np.array([
            i1, i1+1, i1+2,
            i2, i2+1, i2+2,
        ], dtype='int32')
    elif k_axial == 0.0: # torsion; assume 3D
        K2 = K * k_torsion
        n_ijv = [
            # torsion
            ni1 + 3, ni1 + 4, ni1 + 5,  # node 1
            nj1 + 3, nj1 + 4, nj1 + 5,  # node 2
        ]
        dofs = np.array([
            i1, i1+1, i1+2,
            i2, i2+1, i2+2,
        ], dtype='int32')
    else:  # axial + torsion; assume 3D
        # u1fx, u1fy, u1fz, u2fx, u2fy, u2fz
        K2[:nki, :nki] = K * k_axial

        # u1mx, u1my, u1mz, u2mx, u2my, u2mz
        K2[nki:, nki:] = K * k_torsion

        dofs = np.array([
            i1, i1+1, i1+2,
            i2, i2+1, i2+2,

            i1+3, i1+4, i1+5,
            i2+3, i2+4, i2+5,
        ], dtype='int32')
        n_ijv = [
            # axial
            ni1, ni1 + 1, ni1 + 2,  # node 1
            nj1, nj1 + 1, nj1 + 2,  # node 2

            # torsion
            ni1 + 3, ni1 + 4, ni1 + 5,  # node 1
            nj1 + 3, nj1 + 4, nj1 + 5,  # node 2
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

def _build_kbb_cbeam(model: BDF, Kbb, dof_map: DOF_MAP,
                     all_nids, xyz_cid0, idtype='int32', fdtype='float64') -> int:
    """TODO: Timoshenko beam, warping, I12"""
    str(all_nids)
    str(xyz_cid0)
    eids = np.array(model._type_to_id_map['CBEAM'], dtype=idtype)
    nelements = len(eids)
    if nelements == 0:
        return nelements

    for eid in eids:
        elem = model.elements[eid]
        nid1, nid2 = elem.nodes
        xyz1 = elem.nodes_ref[0].get_position()
        xyz2 = elem.nodes_ref[1].get_position()
        dxyz = xyz2 - xyz1
        L = np.linalg.norm(dxyz)
        pid_ref = elem.pid_ref
        mat = pid_ref.mid_ref
        is_failed, (wa, wb, ihat, jhat, khat) = elem.get_axes(model)
        #print(wa, wb, ihat, jhat, khat)
        assert is_failed is False
        T = np.vstack([ihat, jhat, khat])
        z = np.zeros((3, 3), dtype=fdtype)
        Teb = np.block([
            [T, z, z, z],
            [z, T, z, z],
            [z, z, T, z],
            [z, z, z, T],
        ])
        Iy = pid_ref.I11()
        Iz = pid_ref.I22()
        k1 = pid_ref.k1
        k2 = pid_ref.k2
        pa = elem.pa
        pb = elem.pb
        Ke = _beami_stiffness(pid_ref, mat, L, Iy, Iz, pa, pb, k1=k1, k2=k2)
        K = Teb.T @ Ke @ Teb
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
    return nelements

def _beami_stiffness(prop: Union[PBAR, PBARL, PBEAM, PBEAML],
                     mat: MAT1,
                     L: float, Iy: float, Iz: float,
                     pa: int, pb: int,
                     k1: Optional[float]=None,
                     k2: Optional[float]=None):
    """gets the ith Euler-Bernoulli beam stiffness"""
    E = mat.E()
    G = mat.G()
    A = prop.Area()
    J = prop.J()

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

def Kbb_to_Kgg(model: BDF, Kbb: NDArrayNNfloat,
               ngrid: int, ndof_per_grid: int, inplace=True) -> NDArrayNNfloat:
    """does an in-place transformation"""
    assert isinstance(Kbb, (np.ndarray, sci_sparse.csc.csc_matrix)), type(Kbb)
    #assert isinstance(Kbb, (np.ndarray, sci_sparse.csc.csc_matrix, sci_sparse.dok.dok_matrix)), type(Kbb)
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
