from __future__ import annotations
from typing import Tuple, Any, TYPE_CHECKING
#from typing import List, Dict, Tuple, Union, Any

import numpy as np
import scipy.sparse as sci_sparse

from .utils import lambda1d
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


def build_Kbb(model: BDF, dof_map, ndof, dtype='float32') -> Tuple[np.array, Any]:
    """[K] = d{P}/dx"""
    Kbb = np.zeros((ndof, ndof), dtype=dtype)
    Kbbs = sci_sparse.dok_matrix((ndof, ndof), dtype=dtype)
    #print(dof_map)

    #_get_loadid_ndof(model, subcase_id)
    nelements = 0
    nelements += _build_kbb_celas1(model, Kbb, Kbbs, dof_map)
    nelements += _build_kbb_celas2(model, Kbb, Kbbs, dof_map)
    nelements += _build_kbb_conrod(model, Kbb, Kbbs, dof_map)
    nelements += _build_kbb_crod(model, Kbb, Kbbs, dof_map)
    nelements += _build_kbb_ctube(model, Kbb, Kbbs, dof_map)
    nelements += _build_kbb_cbar(model, Kbb, Kbbs, dof_map)
    assert nelements > 0, nelements
    Kbbs2 = Kbbs.tocsc()
    Kbb2 = Kbbs2.toarray()
    error = np.linalg.norm(Kbb - Kbb2)
    if error > 1e-12:
        model.log.warning(f'error = {error}')
    return Kbb, Kbbs2


def _build_kbb_celas1(model: BDF, Kbb, Kbbs, dof_map):
    """fill the CELAS1 Kbb matrix"""
    celas1s = model._type_to_id_map['CELAS1']
    for eid in celas1s:
        elem = model.elements[eid]
        ki = elem.K()
        #print(elem, ki)
        #print(elem.get_stats())
        _build_kbbi_celas12(Kbb, Kbbs, dof_map, elem, ki)
    return len(celas1s)

def _build_kbb_celas2(model: BDF, Kbb, Kbbs, dof_map):
    """fill the CELAS2 Kbb matrix"""
    celas2s = model._type_to_id_map['CELAS2']
    #celas3s = model._type_to_id_map['CELAS3']
    #celas4s = model._type_to_id_map['CELAS4']
    for eid in celas2s:
        elem = model.elements[eid]
        ki = elem.K()
        #print(elem, ki)
        #print(elem.get_stats())
        _build_kbbi_celas12(Kbb, Kbbs, dof_map, elem, ki)
    return len(celas2s)

def _build_kbbi_celas12(Kbb, Kbbs, dof_map, elem, ki):
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
            Kbbs[ib1, ib2] += k[ie1, ie2]
    #Kbb[j, i] += ki
    #Kbb[i, j] += ki
    #del i, j, ki, nid1, nid2, c1, c2

def _build_kbb_cbar(model, Kbb, Kbbs, dof_map):
    """fill the CBAR Kbb matrix"""
    cbars = model._type_to_id_map['CBAR']
    for eid in cbars:
        elem = model.elements[eid]
        pid_ref = elem.pid_ref
        mat = pid_ref.mid_ref
        _build_kbbi_conrod_crod(Kbb, Kbbs, dof_map, elem, mat)
    return len(cbars)

def _build_kbb_crod(model, Kbb, Kbbs, dof_map):
    """fill the CROD Kbb matrix"""
    crods = model._type_to_id_map['CROD']
    for eid in crods:
        elem = model.elements[eid]
        pid_ref = elem.pid_ref
        mat = pid_ref.mid_ref
        _build_kbbi_conrod_crod(Kbb, Kbbs, dof_map, elem, mat)
    return len(crods)

def _build_kbb_ctube(model: BDF, Kbb, Kbbs, dof_map):
    """fill the CTUBE Kbb matrix"""
    ctubes = model._type_to_id_map['CTUBE']
    for eid in ctubes:
        elem = model.elements[eid]
        pid_ref = elem.pid_ref
        mat = pid_ref.mid_ref
        _build_kbbi_conrod_crod(Kbb, Kbbs, dof_map, elem, mat)
    return len(ctubes)

def _build_kbb_conrod(model: BDF, Kbb, Kbbs, dof_map):
    """fill the CONROD Kbb matrix"""
    conrods = model._type_to_id_map['CONROD']
    for eid in conrods:
        elem = model.elements[eid]
        mat = elem.mid_ref
        _build_kbbi_conrod_crod(Kbb, Kbbs, dof_map, elem, mat)
    return len(conrods)

def _build_kbbi_conrod_crod(Kbb, Kbbs, dof_map, elem, mat, fdtype='float64'):
    """fill the ith rod Kbb matrix"""
    nid1, nid2 = elem.nodes
    #mat = elem.mid_ref
    xyz1 = elem.nodes_ref[0].get_position()
    xyz2 = elem.nodes_ref[1].get_position()
    dxyz12 = xyz1 - xyz2
    L = np.linalg.norm(dxyz12)
    E = mat.E
    G = mat.G()
    J = elem.J()
    A = elem.Area()
    E = elem.E()
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

    i1 = 0
    i2 = 3 # dof_map[(n1, 2)]
    if k_torsion == 0.0: # axial; 2D or 3D
        K2 = K * k_axial
        n_ijv = [
            # axial
            dof_map[(nid1, 1)], dof_map[(nid1, 2)], dof_map[(nid1, 3)],
            dof_map[(nid2, 1)], dof_map[(nid2, 2)], dof_map[(nid2, 3)],
        ]
        dofs = np.array([
            i1, i1+1, i1+2,
            i2, i2+1, i2+2,
        ], dtype='int32')
    elif k_axial == 0.0: # torsion; assume 3D
        K2 = K * k_torsion
        n_ijv = [
            # torsion
            dof_map[(nid1, 4)], dof_map[(nid1, 5)], dof_map[(nid1, 6)],
            dof_map[(nid2, 4)], dof_map[(nid2, 5)], dof_map[(nid2, 6)],
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
            dof_map[(nid1, 1)], dof_map[(nid1, 2)], dof_map[(nid1, 3)],
            dof_map[(nid2, 1)], dof_map[(nid2, 2)], dof_map[(nid2, 3)],

            # torsion
            dof_map[(nid1, 4)], dof_map[(nid1, 5)], dof_map[(nid1, 6)],
            dof_map[(nid2, 4)], dof_map[(nid2, 5)], dof_map[(nid2, 6)],
        ]
    for dof1, i1 in zip(dofs, n_ijv):
        for dof2, i2 in zip(dofs, n_ijv):
            ki = K2[dof1, dof2]
            if abs(ki) > 0.:
                #print(nij1, nij2, f'({i1}, {i2});', (dof1, dof2), ki)
                Kbb[i1, i2] = ki
                Kbbs[i1, i2] = ki
        #print(K2)
    #print(Kbb)
    return
