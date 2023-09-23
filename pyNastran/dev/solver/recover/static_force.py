from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.dev.solver.utils import lambda1d, get_ieids_eids
from pyNastran.op2.op2_interface.op2_classes import (
    RealSpringForceArray, RealRodForceArray, RealCBarForceArray,
)
from pyNastran.dev.solver.build_stiffness import ke_cbar
from .static_spring import _recover_force_celas
from .utils import get_plot_request

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF, Subcase, CBAR, PBAR, PBARL


def recover_force_101(f06_file, op2,
                       model: BDF, dof_map, subcase: Subcase, xb, fdtype: str='float32',
                       title: str='', subtitle: str='', label: str='',
                       page_num: int=1, page_stamp: str='PAGE %s'):
    """
    recovers the forces from:
     - FORCE = ALL

    """
    eid_str = 'ALL'
    unused_eids_write, write_f06, write_op2, quick_return = get_plot_request(
        subcase, 'FORCE')
    if quick_return:
        return page_num
    isubcase = subcase.id

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
    nelements += _recover_force_cbar(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CBAR', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    if nelements == 0:
        model.log.warning(f'no force output...{model.card_count}; {model.bdf_filename}')

def _recover_force_rod(f06_file, op2,
                       model: BDF, dof_map, isubcase, xb, eids_str,
                       element_name, fdtype='float32',
                       title: str='', subtitle: str='', label: str='',
                       page_num: int=1, page_stamp='PAGE %s') -> None:
    """recovers static rod force"""
    neids, irod, eids = get_ieids_eids(model, element_name, eids_str)
    if not neids:
        return neids
    forces = np.full((neids, 2), np.nan, dtype=fdtype)
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

def _recover_forcei_rod(xb, dof_map, elem, prop):
    """get the static rod force"""
    nid1, nid2 = elem.nodes

    i1 = dof_map[(nid1, 1)]
    i2 = dof_map[(nid2, 1)]

    q_axial = np.array([
        xb[i1], xb[i1+1], xb[i1+2],
        xb[i2], xb[i2+1], xb[i2+2]
    ])
    q_torsion = np.array([
        xb[i1+3], xb[i1+4], xb[i1+5],
        xb[i2+3], xb[i2+4], xb[i2+5]
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

def _recover_force_cbar(f06_file, op2,
                        model: BDF, dof_map, isubcase, xb, eids_str,
                        element_name, fdtype='float32',
                        title: str='', subtitle: str='', label: str='',
                        page_num: int=1, page_stamp='PAGE %s') -> None:
    """
    Recovers static CBAR force.

    .. todo:: doesn't support CBAR-100

    """
    neids, irod, eids = get_ieids_eids(model, element_name, eids_str)
    if not neids:
        return neids
    forces = np.full((neids, 8), np.nan, dtype=fdtype)

    for ieid, eid in zip(irod, eids):
        elem = model.elements[eid]
        forces[ieid, :] = _recover_forcei_cbar(model, xb, dof_map, elem, elem.pid_ref)

    data = forces.reshape(1, *forces.shape)
    table_name = 'OEF1'
    force_obj = RealCBarForceArray.add_static_case(
        table_name, 'CBAR', eids, data, isubcase,
        is_sort1=True, is_random=False, is_msc=True,
        random_code=0, title=title, subtitle=subtitle, label=label)

    op2.cbar_force[isubcase] = force_obj

    force_obj.write_f06(f06_file, header=None, page_stamp=page_stamp,
                        page_num=page_num, is_mag_phase=False, is_sort1=True)
    return neids

def _recover_forcei_cbar(model: BDF,
                         xb, dof_map, elem: CBAR,
                         prop: Union[PBAR, PBARL], fdtype: str='float64'):
    """get the static CBAR force"""
    #words = ['                                 F O R C E S   I N   B A R   E L E M E N T S         ( C B A R )\n',
             #'0    ELEMENT         BEND-MOMENT END-A            BEND-MOMENT END-B                - SHEAR -               AXIAL\n',
             #'       ID.         PLANE 1       PLANE 2        PLANE 1       PLANE 2        PLANE 1       PLANE 2         FORCE         TORQUE\n']
    nid1, nid2 = elem.nodes
    prop = elem.pid_ref
    mat = prop.mid_ref

    i1 = dof_map[(nid1, 1)]
    i2 = dof_map[(nid2, 1)]

    q_all = np.hstack([
        xb[i1:i1+6],
        xb[i2:i2+6],
    ])
    #q_axial = np.array([
        #xb[i1], xb[i1+1], xb[i1+2],
        #xb[i2], xb[i2+1], xb[i2+2]
    #])
    #q_torsion = np.array([
        #xb[i1+3], xb[i1+4], xb[i1+5],
        #xb[i2+3], xb[i2+4], xb[i2+5]
    #])

    #u_axial = Lambda @ q_axial
    #u_torsion = Lambda @ q_torsion

    nid1, nid2 = elem.nodes
    is_passed, Ke = ke_cbar(model, elem, fdtype=fdtype)
    assert is_passed

    #pid_ref = elem.pid_ref
    #mat = pid_ref.mid_ref


    # ------------------
    is_failed, (v, ihat, jhat, khat, wa, wb) = elem.get_axes(model)
    assert is_failed is False
    #print(wa, wb)
    #xyz1 = elem.nodes_ref[0].get_position() + wa
    #xyz2 = elem.nodes_ref[1].get_position() + wb
    #dxyz = xyz2 - xyz1
    #L = np.linalg.norm(dxyz)
    #pid_ref = elem.pid_ref
    #mat = pid_ref.mid_ref
    T = np.vstack([ihat, jhat, khat])
    z = np.zeros((3, 3), dtype=fdtype)
    Teb = np.block([
        [T, z, z, z],
        [z, T, z, z],
        [z, z, T, z],
        [z, z, z, T],
    ]) # 12x12
    q_element = Teb @ q_all
    u_e = q_element.reshape(12, 1)
    Fe = Ke @ q_element

    # ---------------------------------
    f_e = Fe
    #c, d, e, f = prop.get_cdef()
    #C1, C2 = c
    #D1, D2 = d
    #E1, E2 = e
    #F1, F2 = f
    C1, C2, D1, D2, E1, E2, F1, F2 = prop.get_cdef().ravel()
    A = prop.A
    E = mat.E()
    I1 = prop.I11()
    I2 = prop.I22()

    #f_e = obj.k_e * u_e
    force2stress = np.array([
        [1/A, 0, 0, 0, C2/I2, -C1/I1],
        [1/A, 0, 0, 0, D2/I2, -D1/I1],
        [1/A, 0, 0, 0, E2/I2, -E1/I1],
        [1/A, 0, 0, 0, F2/I2, -F1/I1],
    ])

    #print(force2stress.shape)
    #print(f_e.shape)
    stress = np.hstack([
        -force2stress @ f_e[:6],
        force2stress @ f_e[6:],
    ])
    #print(E, stress)

    # [End A Long. Stress or Strain at Point C;
    #  End A Long. Stress or Strain at Point D;
    #  End A Long. Stress or Strain at Point E;
    #  End A Long. Stress or Strain at Point F;
    #  End B Long. Stress or Strain at Point C;
    #  End B Long. Stress or Strain at Point D;
    #  End B Long. Stress or Strain at Point E;
    #  End B Long. Stress or Strain at Point F]
    strain = 1 / E * stress

    strain_energy = 0.5 * np.diag(u_e.T @ f_e).T
    #print('strain_energy =', strain_energy)
    # ---------------------------------
    #k1 = pid_ref.k1
    #k2 = pid_ref.k2
    #Ke = _beami_stiffness(pid_ref, mat, L, I1, I2, k1=k1, k2=k2, pa=pa, pb=pb)
    #K = Teb.T @ Ke @ Teb

    #is_passed, (v, ihat, jhat, khat, wa, wb) = elem.get_axes(model)
    #T = np.vstack([ihat, jhat, khat])
    #z = np.zeros((3, 3), dtype='float64')
    #I1 = prop.I11()
    #I2 = prop.I22()
    #A = prop.Area()
    #J = prop.J()
    #unused_I12 = prop.I12()

    (Fx1, Fy1, Fz1, Mx1, My1, Mz1,
    Fx2, Fy2, Fz2, Mx2, My2, Mz2) = Fe

    axial = Fx1
    torque = Mx1
    shear1 = Fy1
    shear2 = Fz1
    bending_moment_a1 = My1
    bending_moment_a2 = Mz1

    bending_moment_b1 = My2
    bending_moment_b2 = Mz2

    out = (
        bending_moment_a1, bending_moment_a2,
        bending_moment_b1, bending_moment_b2,
        shear1, shear2,
        axial, torque)
    return out
