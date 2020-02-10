from __future__ import annotations
import numpy as np
from typing import TYPE_CHECKING

from pyNastran.dev.solver.build_stiffness import _lambda1d
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


def get_ieids_eids(model: BDF, etype, eids_str, ncols: int=1, idtype='int32', fdtype='float32'):
    """helper for the stress/strain/displacment recovery"""
    eids = model._type_to_id_map[etype]
    if len(eids) == 0:
        return 0, None, None, None

    if eids_str == 'ALL':
        neids = len(eids)
        ieids = np.arange(neids, dtype=idtype)
    else:
        ieids = np.searchsorted(eids_str, eids)
        neids = len(ieids)
    shape = (neids, ncols)
    empty_array = np.full(shape, np.nan, dtype=fdtype)
    return neids, ieids, eids, empty_array

def recover_strain_101(model, xg, dof_map, fdtype='float32'):
    """recovers the strains"""
    eids = 'ALL'
    nelements = 0
    nelements += _recover_strain_celas1(model, dof_map, xg, eids, fdtype=fdtype)
    nelements += _recover_strain_celas2(model, dof_map, xg, eids, fdtype=fdtype)
    nelements += _recover_strain_conrod(model, dof_map, xg, eids, fdtype=fdtype)
    #assert nelements > 0, nelements

def _recover_strain_celas1(model: BDF, dof_map, xg, eids, fdtype='float32'):
    """
    recovers static strain
    TODO: write the OP2/F06

    """
    neids, ielas, celas1s, strains = get_ieids_eids(model, 'CELAS1', eids, fdtype=fdtype)
    if not neids:
        return neids
    for ieid, eid in zip(ielas, celas1s):
        elem = model.elements[eid]
        strain = _recover_straini_celas12(xg, dof_map, elem)
        strains[ielas] = strain
    #spring_strain = RealSpringStrainArray.add_static_case(
        #table_name, node_gridtype, data, isubcase,
        #is_sort1=True, is_random=False, is_msc=True,
        #random_code=0, title=title, subtitle=subtitle, label=label)
    return neids # , spring_strain

def _recover_strain_celas2(model: BDF, dof_map, xg, eids, fdtype='float32'):
    """
    recovers static strain
    TODO: write the OP2/F06

    """
    neids, ielas, celas2s, strains = get_ieids_eids(model, 'CELAS2', eids, fdtype=fdtype)
    if not neids:
        return neids
    #celas3s = model._type_to_id_map['CELAS3']
    #celas4s = model._type_to_id_map['CELAS4']
    for ieid, eid in zip(ielas, celas2s):
        elem = model.elements[eid]
        strain = _recover_straini_celas12(xg, dof_map, elem)
        strains[ielas] = strain
    return neids

def _recover_strain_conrod(model: BDF, dof_map, xg, eids, fdtype='float32'):
    """
    recovers static strain
    TODO: write the OP2/F06

    """
    neids, irod, conrods, strains = get_ieids_eids(model, 'CONROD', eids, ncols=2, fdtype=fdtype)
    if not neids:
        return neids
    for ieid, eid in zip(irod, conrods):
        elem = model.elements[eid]
        axiali, torsioni = _recover_straini_rod(xg, dof_map, elem)
        strains[ieid, :] = axiali, torsioni
    return neids

def _recover_straini_celas12(xg, dof_map, elem):
    """get the static strain"""
    nid1, nid2 = elem.nodes
    c1, c2 = elem.c1, elem.c2
    i = dof_map[(nid1, c1)]
    j = dof_map[(nid2, c2)]
    strain = xg[i] - xg[j]  # TODO: check the sign
    return strain

def _recover_straini_rod(xb, dof_map, elem):
    """get the static strain"""
    nid1, nid2 = elem.nodes

    # axial
    i11 = dof_map[(nid1, 1)]
    i12 = dof_map[(nid1, 2)]
    i13 = dof_map[(nid1, 3)]

    i21 = dof_map[(nid2, 1)]
    i22 = dof_map[(nid2, 2)]
    i23 = dof_map[(nid2, 3)]

    # torsion
    i14 = dof_map[(nid1, 4)]
    i15 = dof_map[(nid1, 5)]
    i16 = dof_map[(nid1, 6)]

    i24 = dof_map[(nid2, 4)]
    i25 = dof_map[(nid2, 5)]
    i26 = dof_map[(nid2, 6)]

    q_axial = np.array([
        xb[i11], xb[i12], xb[i13],
        xb[i21], xb[i22], xb[i23]
    ])
    q_torsion = np.array([
        xb[i14], xb[i15], xb[i16],
        xb[i24], xb[i25], xb[i26]
    ])
    #print("type=%s n1=%s n2=%s" % (self.type, n1, n2))
    #print("n11=%s n12=%s n21=%s n22=%s" %(n11,n12,n21,n22))

    #print("q2[%s] = %s" % (self.eid, q2))
    #print("Lambda = \n"+str(Lambda))

    #print("Lsize = ", Lambda.shape)
    #print("qsize = ", q.shape)
    xyz1 = elem.nodes_ref[0].get_position()
    xyz2 = elem.nodes_ref[1].get_position()
    dxyz12 = xyz1 - xyz2
    Lambda = _lambda1d(dxyz12, debug=False)

    u_axial = Lambda @ q_axial
    u_torsion = Lambda @ q_torsion
    du_axial = u_axial[0] - u_axial[1]
    du_torsion = u_torsion[0] - u_torsion[1]

    return du_axial, du_torsion

