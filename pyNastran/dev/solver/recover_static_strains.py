from __future__ import annotations
import numpy as np
from typing import TYPE_CHECKING

from pyNastran.dev.solver.build_stiffness import lambda1d
from pyNastran.op2.op2_interface.hdf5_interface import (
    #RealSpringForceArray,
    RealSpringStrainArray, # RealSpringStressArray,
    #RealRodForceArray,
    RealRodStrainArray, # RealRodStressArray,
)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


def get_ieids_eids(model: BDF, etype, eids_str, ncols: int=1, idtype='int32', fdtype='float32'):
    """helper for the stress/strain/displacment recovery"""
    eids = np.array(model._type_to_id_map[etype], dtype=idtype)
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

def recover_strain_101(f06_file, op2,
                       model: BDF, dof_map, isubcase: int, xb, fdtype: str='float32',
                       title: str='', subtitle: str='', label: str='',
                       page_num: int=1, page_stamp: str='PAGE %s'):
    """recovers the strains"""
    eid_str = 'ALL'
    nelements = 0
    nelements += _recover_strain_celas(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CELAS1', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += _recover_strain_celas(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CELAS2', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    #celas3s = model._type_to_id_map['CELAS3']
    #celas4s = model._type_to_id_map['CELAS4']

    nelements += _recover_strain_rod(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CROD', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += _recover_strain_rod(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CTUBE', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += _recover_strain_rod(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CONROD', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    #assert nelements > 0, nelements

def _recover_strain_celas(f06_file, op2,
                          model: BDF, dof_map, isubcase, xg, eids_str,
                          element_name: str, fdtype='float32',
                          title: str='', subtitle: str='', label: str='',
                          page_num: int=1, page_stamp='PAGE %s') -> None:
    """
    recovers static strain
    TODO: write the OP2/F06

    """
    neids, ielas, eids, strains = get_ieids_eids(model, element_name, eids_str, fdtype=fdtype)
    if not neids:
        return neids
    for ieid, eid in zip(ielas, eids):
        elem = model.elements[eid]
        strain = _recover_straini_celas12(xg, dof_map, elem)
        strains[ielas] = strain

    data = strains.reshape(1, *strains.shape)
    table_name = 'OSTR1'
    spring_strain = RealSpringStrainArray.add_static_case(
        table_name, element_name, eids, data, isubcase, is_stress=False,
        is_sort1=True, is_random=False, is_msc=True,
        random_code=0, title=title, subtitle=subtitle, label=label)
    op2.celas1_strain[isubcase] = spring_strain
    spring_strain.write_f06(f06_file, header=None, page_stamp=page_stamp,
                            page_num=page_num, is_mag_phase=False, is_sort1=True)
    return neids

def _recover_strain_rod(f06_file, op2,
                        model: BDF, dof_map, isubcase, xb, eids_str,
                        element_name, fdtype='float32',
                        title: str='', subtitle: str='', label: str='',
                        page_num: int=1, page_stamp='PAGE %s') -> None:
    """
    recovers static strain
    TODO: write the OP2/F06

    """
    neids, irod, eids, strains = get_ieids_eids(model, element_name, eids_str, ncols=4, fdtype=fdtype)
    if not neids:
        return neids
    for ieid, eid in zip(irod, eids):
        elem = model.elements[eid]
        strains[ieid, :] = _recover_straini_rod(xb, dof_map, elem)

    data = strains.reshape(1, *strains.shape)
    table_name = 'OSTR1'
    spring_strain = RealRodStrainArray.add_static_case(
        table_name, element_name, eids, data, isubcase, is_stress=False,
        is_sort1=True, is_random=False, is_msc=True,
        random_code=0, title=title, subtitle=subtitle, label=label)
    op2.celas1_strain[isubcase] = spring_strain
    spring_strain.write_f06(f06_file, header=None, page_stamp=page_stamp,
                            page_num=page_num, is_mag_phase=False, is_sort1=True)

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
    Lambda = lambda1d(dxyz12, debug=False)

    u_axial = Lambda @ q_axial
    u_torsion = Lambda @ q_torsion
    du_axial = u_axial[0] - u_axial[1]
    du_torsion = u_torsion[0] - u_torsion[1]
    #headers = ['axial', 'SMa', 'torsion', 'SMt']

    return du_axial, 0., du_torsion, 0.

