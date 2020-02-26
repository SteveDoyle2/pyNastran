from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.dev.solver.utils import lambda1d, get_ieids_eids
from pyNastran.op2.op2_interface.hdf5_interface import (
    RealSpringForceArray, RealRodForceArray,
)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


def recover_force_101(f06_file, op2,
                       model: BDF, dof_map, isubcase: int, xb, fdtype: str='float32',
                       title: str='', subtitle: str='', label: str='',
                       page_num: int=1, page_stamp: str='PAGE %s'):
    """
    recovers the forces from:
     - FORCE = ALL

    """
    eid_str = 'ALL'
    nelements = 0
    nelements += _recover_force_celas(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CELAS1', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += _recover_force_celas(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CELAS2', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += _recover_force_celas(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CELAS3', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += _recover_force_celas(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CELAS4', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)

    nelements += _recover_force_rod(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CROD', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += _recover_force_rod(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CONROD', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += _recover_force_rod(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CTUBE', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    if nelements == 0:
        model.log.warning(f'no force output...{model.card_count}')

def _recover_force_celas(f06_file, op2,
                         model: BDF, dof_map, isubcase, xg, eids_str,
                         element_name: str, fdtype='float32',
                         title: str='', subtitle: str='', label: str='',
                         page_num: int=1, page_stamp='PAGE %s') -> None:
    """recovers static spring force"""
    neids, ielas, eids, forces = get_ieids_eids(model, element_name, eids_str, fdtype=fdtype)
    if not neids:
        return neids
    if element_name in {'CELAS1', 'CELAS2'}:
        for ieid, eid in zip(ielas, eids):
            elem = model.elements[eid]
            ki = elem.K()
            force = _recover_forcei_celas12(xg, dof_map, elem, ki)
            forces[ieid] = force
    elif element_name in {'CELAS3', 'CELAS4'}:
        for ieid, eid in zip(ielas, eids):
            elem = model.elements[eid]
            ki = elem.K()
            force = _recover_forcei_celas34(xg, dof_map, elem, ki)
            forces[ieid] = force
    else:  # pragma: no cover
        raise NotImplementedError(element_name)

    data = forces.reshape(1, *forces.shape)
    table_name = 'OEF1'
    spring_force = RealSpringForceArray.add_static_case(
        table_name, element_name, eids, data, isubcase,
        is_sort1=True, is_random=False, is_msc=True,
        random_code=0, title=title, subtitle=subtitle, label=label)
    if element_name == 'CELAS1':
        op2.celas1_force[isubcase] = spring_force
    elif element_name == 'CELAS2':
        op2.celas2_force[isubcase] = spring_force
    elif element_name == 'CELAS3':
        op2.celas3_force[isubcase] = spring_force
    elif element_name == 'CELAS4':
        op2.celas4_force[isubcase] = spring_force
    else:  # pragma: no cover
        raise NotImplementedError(element_name)
    spring_force.write_f06(f06_file, header=None, page_stamp=page_stamp,
                            page_num=page_num, is_mag_phase=False, is_sort1=True)
    return neids

def _recover_force_rod(f06_file, op2,
                       model: BDF, dof_map, isubcase, xb, eids_str,
                       element_name, fdtype='float32',
                       title: str='', subtitle: str='', label: str='',
                       page_num: int=1, page_stamp='PAGE %s') -> None:
    """recovers static rod force"""
    neids, irod, eids, forces = get_ieids_eids(model, element_name, eids_str, ncols=2, fdtype=fdtype)
    if not neids:
        return neids
    if element_name == 'CONROD':
        for ieid, eid in zip(irod, eids):
            elem = model.elements[eid]
            forces[ieid, :] = _recover_forcei_rod(xb, dof_map, elem, elem)
    elif element_name == 'CROD':
        for ieid, eid in zip(irod, eids):
            elem = model.elements[eid]
            forces[ieid, :] = _recover_forcei_rod(xb, dof_map, elem, elem.pid_ref)
    elif element_name == 'CTUBE':
        for ieid, eid in zip(irod, eids):
            elem = model.elements[eid]
            forces[ieid, :] = _recover_forcei_rod(xb, dof_map, elem, elem.pid_ref)
    else:  # pragma: no cover
        raise NotImplementedError(element_name)

    data = forces.reshape(1, *forces.shape)
    table_name = 'OEF1'
    force_obj = RealRodForceArray.add_static_case(
        table_name, element_name, eids, data, isubcase,
        is_sort1=True, is_random=False, is_msc=True,
        random_code=0, title=title, subtitle=subtitle, label=label)

    if element_name == 'CONROD':
        op2.conrod_force[isubcase] = force_obj
    elif element_name == 'CROD':
        op2.crod_force[isubcase] = force_obj
    elif element_name == 'CTUBE':
        op2.ctube_force[isubcase] = force_obj
    else:  # pragma: no cover
        raise NotImplementedError(element_name)

    force_obj.write_f06(f06_file, header=None, page_stamp=page_stamp,
                        page_num=page_num, is_mag_phase=False, is_sort1=True)
    return neids

def _recover_forcei_celas12(xg, dof_map, elem, ki: float):
    """get the static spring force"""
    # F = kx
    nid1, nid2 = elem.nodes
    c1, c2 = elem.c1, elem.c2
    i = dof_map[(nid1, c1)]
    j = dof_map[(nid2, c2)]
    strain = xg[j] - xg[i]  # TODO: check the sign
    force = ki * strain
    return force


def _recover_forcei_celas34(xg, dof_map, elem, ki: float):
    """get the static spring force"""
    # F = kx
    nid1, nid2 = elem.nodes
    i = dof_map[(nid1, 0)]
    j = dof_map[(nid2, 0)]
    strain = xg[j] - xg[i]  # TODO: check the sign
    force = ki * strain
    return force

def _recover_forcei_rod(xb, dof_map, elem, prop):
    """get the static rod force"""
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
    xyz1 = elem.nodes_ref[0].get_position()
    xyz2 = elem.nodes_ref[1].get_position()
    dxyz12 = xyz1 - xyz2
    Lambda = lambda1d(dxyz12, debug=False)

    u_axial = Lambda @ q_axial
    u_torsion = Lambda @ q_torsion
    du_axial = u_axial[0] - u_axial[1]
    du_torsion = u_torsion[0] - u_torsion[1]
    #headers = ['axial', 'SMa', 'torsion', 'SMt']

    #C = prop.c
    mat = prop.mid_ref

    L = np.linalg.norm(dxyz12)
    G = mat.G()
    J = elem.J()
    A = elem.Area()
    E = elem.E()

    axial_strain = du_axial / L
    #torsional_strain = du_torsion * C / L

    axial_stress = E * axial_strain
    #torsional_stress = G * torsional_strain
    axial_force = axial_stress * A
    torsional_moment = du_torsion * G * J / L

    return axial_force, torsional_moment
