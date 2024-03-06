r"""
Hooke law:
   {\sigma^0} = [D][B]{u_0}

"""
from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.dev.solver.utils import lambda1d, get_ieids_eids
from pyNastran.op2.op2_interface.op2_classes import (
    RealRodStressArray,
)
from .static_spring import _recover_stress_celas
from .utils import get_plot_request

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF, Subcase


def recover_stress_101(f06_file, op2,
                       model: BDF, dof_map, subcase: Subcase, xb, fdtype: str='float32',
                       title: str='', subtitle: str='', label: str='',
                       page_num: int=1, page_stamp: str='PAGE %s'):
    """
    recovers the stresses from:
     - STRESS = ALL

    """
    eid_str = 'ALL'
    unused_eids_write, write_f06, write_op2, quick_return = get_plot_request(
        subcase, 'STRESS')
    if quick_return:
        return page_num
    isubcase = subcase.id

    nelements = 0
    nelements += _recover_stress_celas(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CELAS1', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += _recover_stress_celas(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CELAS2', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += _recover_stress_celas(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CELAS3', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += _recover_stress_celas(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CELAS4', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)

    nelements += _recover_stress_rod(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CROD', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += _recover_stress_rod(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CONROD', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += _recover_stress_rod(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CTUBE', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    #assert nelements > 0, nelements
    if nelements == 0:
        model.log.warning(f'no stress output...{model.card_count}; {model.bdf_filename}')


def _recover_stress_rod(f06_file, op2,
                        model: BDF, dof_map, isubcase, xb, eids_str,
                        element_name, fdtype='float32',
                        title: str='', subtitle: str='', label: str='',
                        page_num: int=1, page_stamp='PAGE %s') -> None:
    """recovers static rod stress"""
    neids, irod, eids = get_ieids_eids(model, element_name, eids_str)
    if not neids:
        return neids
    stresses = np.full((neids, 4), np.nan, dtype=fdtype)
    if element_name == 'CONROD':
        for ieid, eid in zip(irod, eids):
            elem = model.elements[eid]
            stresses[ieid, :] = _recover_stressi_rod(xb, dof_map, elem, elem)
    elif element_name == 'CROD':
        for ieid, eid in zip(irod, eids):
            elem = model.elements[eid]
            stresses[ieid, :] = _recover_stressi_rod(xb, dof_map, elem, elem.pid_ref)
    elif element_name == 'CTUBE':
        for ieid, eid in zip(irod, eids):
            elem = model.elements[eid]
            stresses[ieid, :] = _recover_stressi_ctube(xb, dof_map, elem, elem.pid_ref)
    else:  # pragma: no cover
        raise NotImplementedError(element_name)

    data = stresses.reshape(1, *stresses.shape)
    table_name = 'OSTR1'
    stress_obj = RealRodStressArray.add_static_case(
        table_name, element_name, eids, data, isubcase,
        is_sort1=True, is_random=False, is_msc=True,
        random_code=0, title=title, subtitle=subtitle, label=label)

    stress = op2.op2_results.stress
    if element_name == 'CONROD':
        stress.conrod_stress[isubcase] = stress_obj
    elif element_name == 'CROD':
        stress.crod_stress[isubcase] = stress_obj
    elif element_name == 'CTUBE':
        stress.ctube_stress[isubcase] = stress_obj
    else:  # pragma: no cover
        raise NotImplementedError(element_name)

    stress_obj.write_f06(f06_file, header=None, page_stamp=page_stamp,
                         page_num=page_num, is_mag_phase=False, is_sort1=True)
    return neids

def _recover_stressi_rod(xb, dof_map, elem, prop):
    """get the static rod stress"""
    nid1, nid2 = elem.nodes

    # axial
    i11 = dof_map[(nid1, 1)]
    i12 = i11 + 1
    i13 = i11 + 2

    i21 = dof_map[(nid2, 1)]
    i22 = i21 + 1
    i23 = i21 + 2

    # torsion
    i14 = i11 + 3
    i15 = i11 + 4
    i16 = i11 + 5

    i24 = i21 + 3
    i25 = i21 + 4
    i26 = i21 + 5

    q_axial = np.array([
        xb[i11], xb[i12], xb[i13],
        xb[i21], xb[i22], xb[i23],
    ])
    q_torsion = np.array([
        xb[i14], xb[i15], xb[i16],
        xb[i24], xb[i25], xb[i26],
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

    C = prop.c
    mat = prop.mid_ref

    L = np.linalg.norm(dxyz12)
    G = mat.G()
    E = elem.E()

    axial_strain = du_axial / L
    torsional_strain = du_torsion * C / L

    axial_stress = E * axial_strain
    torsional_stress = G * torsional_strain

    return axial_stress, np.nan, torsional_stress, np.nan

def _recover_stressi_ctube(xb, dof_map, elem, prop):
    """get the static ctube stress"""
    nid1, nid2 = elem.nodes

    # axial
    i11 = dof_map[(nid1, 1)]
    i12 = i11 + 1
    i13 = i11 + 2

    i21 = dof_map[(nid2, 1)]
    i22 = i21 + 1
    i23 = i21 + 2

    # torsion
    i14 = i11 + 3
    i15 = i11 + 4
    i16 = i11 + 5

    i24 = i21 + 3
    i25 = i21 + 4
    i26 = i21 + 5

    q_axial = np.array([
        xb[i11], xb[i12], xb[i13],
        xb[i21], xb[i22], xb[i23],
    ])
    q_torsion = np.array([
        xb[i14], xb[i15], xb[i16],
        xb[i24], xb[i25], xb[i26],
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

    mat = prop.mid_ref

    L = np.linalg.norm(dxyz12)
    G = mat.G()
    E = elem.E()

    axial_strain = du_axial / L
    torsional_strain = du_torsion / L

    axial_stress = E * axial_strain
    torsional_stress = G * torsional_strain

    return axial_stress, np.nan, torsional_stress, np.nan
