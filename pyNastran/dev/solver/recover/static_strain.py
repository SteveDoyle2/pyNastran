from __future__ import annotations
from typing import Union, TYPE_CHECKING
import numpy as np

from pyNastran.nptyping_interface import NDArrayNfloat
from pyNastran.dev.solver.utils import lambda1d, get_ieids_eids
from pyNastran.dev.solver.build_stiffness import ke_cbar
from .static_spring import recover_celas
from pyNastran.op2.op2_interface.op2_classes import (
    #RealStrainEnergyArray,
    RealRodStrainArray,
    RealBarStrainArray,
)
from .utils import get_plot_request

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import (
        BDF, Subcase,
        CTUBE, PTUBE,
        CBAR, PBAR, PBARL,
    )


def recover_strain_101(f06_file, op2,
                       model: BDF, dof_map, subcase: Subcase, xb, fdtype: str='float32',
                       title: str='', subtitle: str='', label: str='',
                       page_num: int=1, page_stamp: str='PAGE %s'):
    """
    recovers the strains from:
     - STRAIN= ALL

    """
    eid_str = 'ALL'
    unused_eids_write, write_f06_strain, write_op2_strain, quick_return1 = get_plot_request(
        subcase, 'STRAIN')
    unused_eids_write, write_f06_stress, write_op2_stress, quick_return2 = get_plot_request(
        subcase, 'STRESS')
    unused_eids_write, write_f06_force, write_op2_force, quick_return3 = get_plot_request(
        subcase, 'FORCE')
    unused_eids_write, write_f06_ese, write_op2_ese, quick_return4 = get_plot_request(
        subcase, 'ESE')
    if all([quick_return1, quick_return2, quick_return3, quick_return4]):
        return page_num
    isubcase = subcase.id

    nelements = 0
    nelements += recover_celas(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CELAS1',
        write_f06_strain, write_op2_strain,
        write_f06_stress, write_op2_stress,
        write_f06_force, write_op2_force,
        write_f06_ese, write_op2_ese,
        fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += recover_celas(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CELAS2',
        write_f06_strain, write_op2_strain,
        write_f06_stress, write_op2_stress,
        write_f06_force, write_op2_force,
        write_f06_ese, write_op2_ese,
        fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += recover_celas(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CELAS3',
        write_f06_strain, write_op2_strain,
        write_f06_stress, write_op2_stress,
        write_f06_force, write_op2_force,
        write_f06_ese, write_op2_ese,
        fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += recover_celas(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CELAS4',
        write_f06_strain, write_op2_strain,
        write_f06_stress, write_op2_stress,
        write_f06_force, write_op2_force,
        write_f06_ese, write_op2_ese,
        fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)

    nelements += _recover_strain_rod(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CROD', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += _recover_strain_rod(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CONROD', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += _recover_strain_rod(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CTUBE', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += _recover_strain_bar(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CBAR', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)


    #assert nelements > 0, nelements
    if nelements == 0:
        model.log.warning(f'no strain output...{model.card_count}; {model.bdf_filename}')

def _recover_strain_rod(f06_file, op2,
                        model: BDF, dof_map, isubcase, xb, eids_str,
                        element_name, fdtype='float32',
                        title: str='', subtitle: str='', label: str='',
                        page_num: int=1, page_stamp='PAGE %s') -> None:
    """recovers static rod strain"""
    neids, irod, eids = get_ieids_eids(model, element_name, eids_str)
    if not neids:
        return neids
    strains = np.full((neids, 4), np.nan, dtype=fdtype)

    if element_name == 'CONROD':
        for ieid, eid in zip(irod, eids):
            elem = model.elements[eid]
            strains[ieid, :] = _recover_straini_rod(xb, dof_map, elem, elem)
    elif element_name == 'CROD':
        for ieid, eid in zip(irod, eids):
            elem = model.elements[eid]
            strains[ieid, :] = _recover_straini_rod(xb, dof_map, elem, elem.pid_ref)
    elif element_name == 'CTUBE':
        for ieid, eid in zip(irod, eids):
            elem = model.elements[eid]
            strains[ieid, :] = _recover_straini_ctube(xb, dof_map, elem, elem.pid_ref)
    else:  # pragma: no cover
        raise NotImplementedError(element_name)

    data = strains.reshape(1, *strains.shape)
    table_name = 'OSTR1'
    strain_obj = RealRodStrainArray.add_static_case(
        table_name, element_name, eids, data, isubcase,
        is_sort1=True, is_random=False, is_msc=True,
        random_code=0, title=title, subtitle=subtitle, label=label)

    if element_name == 'CONROD':
        op2.conrod_strain[isubcase] = strain_obj
    elif element_name == 'CROD':
        op2.crod_strain[isubcase] = strain_obj
    elif element_name == 'CTUBE':
        op2.ctube_strain[isubcase] = strain_obj
    else:  # pragma: no cover
        raise NotImplementedError(element_name)

    strain_obj.write_f06(f06_file, header=None, page_stamp=page_stamp,
                         page_num=page_num, is_mag_phase=False, is_sort1=True)
    return neids

def _recover_straini_rod(xb, dof_map, elem, prop):
    """get the static rod strain"""
    nid1, nid2 = elem.nodes

    # axial
    i1 = dof_map[(nid1, 1)]
    i2 = dof_map[(nid2, 1)]

    q_axial = np.array([
        xb[i1], xb[i1+1], xb[i1+2],
        xb[i2], xb[i2+2], xb[i2+2],
    ])
    q_torsion = np.array([
        xb[i1+3], xb[i1+4], xb[i1+5],
        xb[i2+3], xb[i2+4], xb[i2+5],
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

    xyz1 = elem.nodes_ref[0].get_position()
    xyz2 = elem.nodes_ref[1].get_position()
    dxyz12 = xyz1 - xyz2
    L = np.linalg.norm(dxyz12)

    axial_strain = du_axial / L
    torsional_strain = du_torsion * C / L

    return axial_strain, np.nan, torsional_strain, np.nan

def _recover_straini_ctube(xb, dof_map, elem: CTUBE, prop: PTUBE):
    """get the static ctube strain"""
    nid1, nid2 = elem.nodes

    # axial
    i1 = dof_map[(nid1, 1)]
    i2 = dof_map[(nid2, 1)]

    q_axial = np.array([
        xb[i1], xb[i1+1], xb[i1+2],
        xb[i2], xb[i2+2], xb[i2+2],
    ])
    q_torsion = np.array([
        xb[i1+3], xb[i1+4], xb[i1+5],
        xb[i2+3], xb[i2+4], xb[i2+5],
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

    xyz1 = elem.nodes_ref[0].get_position()
    xyz2 = elem.nodes_ref[1].get_position()
    dxyz12 = xyz1 - xyz2
    L = np.linalg.norm(dxyz12)

    axial_strain = du_axial / L
    torsional_strain = du_torsion / L

    return axial_strain, np.nan, torsional_strain, np.nan

def _recover_strain_bar(f06_file, op2,
                        model: BDF, dof_map, isubcase, xb, eids_str,
                        element_name: str, fdtype='float32',
                        title: str='', subtitle: str='', label: str='',
                        page_num: int=1, page_stamp='PAGE %s') -> None:
    """recovers static rod strain"""
    neids, irod, eids = get_ieids_eids(model, element_name, eids_str)
    if not neids:
        return neids
    strains = np.full((neids, 15), np.nan, dtype=fdtype)

    if element_name == 'CBAR':
        for ieid, eid in zip(irod, eids):
            elem = model.elements[eid]
            #[s1a, s2a, s3a, s4a, axial, smaxa, smina, MS_tension,
            # s1b, s2b, s3b, s4b,        sminb, sminb, MS_compression] - 15
            strains[ieid, :] = _recover_straini_cbar(model, xb, dof_map, elem, elem.pid_ref)
    else:  # pragma: no cover
        raise NotImplementedError(element_name)

    data = strains.reshape(1, *strains.shape)
    table_name = 'OSTR1'

    strain_obj = RealBarStrainArray.add_static_case(
        table_name, element_name, eids, data, isubcase,
        is_sort1=True, is_random=False, is_msc=True,
        random_code=0, title=title, subtitle=subtitle, label=label)

    if element_name == 'CBAR':
        op2.cbar_strain[isubcase] = strain_obj
    else:  # pragma: no cover
        raise NotImplementedError(element_name)

    strain_obj.write_f06(f06_file, header=None, page_stamp=page_stamp,
                         page_num=page_num, is_mag_phase=False, is_sort1=True)
    return neids

def _recover_straini_cbar(model: BDF, xb: NDArrayNfloat,
                          dof_map,
                          elem: CBAR, prop: Union[PBAR, PBARL], fdtype='float64'):
    """get the static bar strain"""
    nid1, nid2 = elem.nodes

    # axial
    i1 = dof_map[(nid1, 1)]
    j1 = dof_map[(nid2, 1)]

    q_all = np.hstack([
        xb[i1:i1+6],
        xb[j1:j1+6],
    ])
    #print(len(xb[i1:i1+3],))
    q_axial = np.hstack([
        xb[i1:i1+3],
        xb[j1:j1+3],
    ])
    #print(q_axial)
    #q_torsion = np.hstack([
        #xb[i1+3:i1+6],
        #xb[j1+3:j1+6],
    #])
    #{F} = [K]{u}

    nid1, nid2 = elem.nodes
    is_passed, Ke = ke_cbar(model, elem, fdtype=fdtype)
    assert is_passed

    pid_ref = elem.pid_ref
    mat = pid_ref.mid_ref

    #is_passed, (wa, wb, ihat, jhat, khat) = elem.get_axes(model)
    #T = np.vstack([ihat, jhat, khat])
    #z = np.zeros((3, 3), dtype='float64')
    prop = elem.pid_ref
    #mat = prop.mid_ref
    I1 = prop.I11()
    I2 = prop.I22()
    A = prop.Area()
    J = prop.J()
    unused_I12 = prop.I12()

    Fe = Ke @ q_all
    (Fx1, Fy1, Fz1, Mx1, My1, Mz1,
    Fx2, Fy2, Fz2, Mx2, My2, Mz2) = Fe
    s_axial = Fx1 / A
    s_torsion = Mx1 / J
    str([s_axial, s_torsion])

    xyz1 = elem.nodes_ref[0].get_position()
    xyz2 = elem.nodes_ref[1].get_position()
    #dxyz12 = xyz1 - xyz2
    #Lambda = lambda1d(dxyz12, debug=False)

    #u_axial = Lambda @ q_axial
    #u_torsion = Lambda @ q_torsion
    #du_axial = u_axial[0] - u_axial[1]
    #du_torsion = u_torsion[0] - u_torsion[1]

    xyz1 = elem.nodes_ref[0].get_position()
    xyz2 = elem.nodes_ref[1].get_position()
    dxyz12 = xyz1 - xyz2
    L = np.linalg.norm(dxyz12)
    Lambda = lambda1d(dxyz12, debug=False)

    u_axial = Lambda @ q_axial
    #u_torsion = Lambda @ q_torsion
    du_axial = u_axial[0] - u_axial[1]
    #du_torsion = u_torsion[0] - u_torsion[1]


    axial = du_axial / L

    # cdef = prop.get_cdef()
    if prop.type in ['PBARL', 'PBAR']:
        cdef = prop.get_cdef()
    else:
        raise NotImplementedError(prop.get_stats())
    model.log.info(f'cdef: {cdef}\n')

    Iy = I1
    Iz = I2
    G = mat.G()
    E = mat.E()
    strains = []
    for (T, My, Mz) in [(Mx1, My1, Mz1), (Mx2, My2, Mz2)]:
        for yz in cdef:
            y, z = yz
            radius = np.linalg.norm(yz)
            unused_torsional_stress = (T * radius) / J
            unused_phi = (T * L) / (G * J)

            stress = My * y / Iy + Mz * z / Iz
            strain = stress / E
            strains.append(strain)

    s1a = strains[0]
    s2a = strains[1]
    s3a = strains[2]
    s4a = strains[3]
    s1b = strains[4]
    s2b = strains[5]
    s3b = strains[6]
    s4b = strains[7]
    smaxa = max(s1a, s2a, s3a, s4a)
    smaxb = max(s1b, s2b, s3b, s4b)

    smina = min(s1a, s2a, s3a, s4a)
    sminb = min(s1b, s2b, s3b, s4b)
    MS_tension = MS_compression = np.nan
    out = (
        s1a, s2a, s3a, s4a, axial, smaxa, smina, MS_tension,
        s1b, s2b, s3b, s4b,        smaxb, sminb, MS_compression) # 15
    return out
