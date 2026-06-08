from __future__ import annotations
from typing import TextIO, Optional, TYPE_CHECKING
import numpy as np

from pyNastran.dev.solver.utils import lambda1d
from pyNastran.dev.bdf_vectorized3.solver.utils import get_ieids_eids, get_element
from pyNastran.op2.op2_interface.op2_classes import (
    RealSpringForceArray, RealRodForceArray, RealCBarForceArray,
)
# from pyNastran.dev.solver.build_stiffness import ke_cbar
from .static_spring import _recover_force_celas
from .utils import get_plot_request

from pyNastran.dev.bdf_vectorized3.solver.elements.beam import (
    timoshenko_stiffness, beam_transform, recover_beam_force,
)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF, Subcase, CBAR, PBAR, PBARL
    DOF_MAP = dict[tuple[int, int], int]


def recover_force_101(f06_file: TextIO, op2: OP2,
                      model: BDF, dof_map: DOF_MAP,
                      subcase: Subcase, xb: np.ndarray, fdtype: str='float32',
                      title: str='', subtitle: str='', label: str='',
                      page_num: int=1, page_stamp: str='PAGE %s'):
    """
    recovers the forces from:
     - FORCE = ALL

    """
    if 'FORCE' not in subcase:
        return

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
    nelements += _recover_force_cbar(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CBEAM', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    if nelements == 0:
        model.log.warning(f'no force output...{model.card_count}; {model.bdf_filename}')
        # raise RuntimeError(f'nelements={nelements} cards={model.card_count}')


def _recover_force_rod(f06_file: TextIO, op2: OP2,
                       model: BDF, dof_map: DOF_MAP, isubcase: int,
                       xb: np.ndarray, eids_str,
                       element_name: str, fdtype='float32',
                       title: str='', subtitle: str='', label: str='',
                       page_num: int=1, page_stamp='PAGE %s') -> int:
    """recovers static rod force"""
    neids, ieids, eids = get_ieids_eids(model, element_name, eids_str)
    if not neids:
        return neids

    elem = get_element(model, element_name, ieids, eids)
    xyz1 = model.grid.get_position_by_node_id(elem.nodes[:, 0])
    xyz2 = model.grid.get_position_by_node_id(elem.nodes[:, 1])

    forces = np.full((neids, 2), np.nan, dtype=fdtype)

    if element_name == 'CONROD':
        prop = elem
    elif element_name == 'CROD':
        pid = elem.property_id
        prop = model.prod.slice_card_by_property_id(pid)
    elif element_name == 'CTUBE':
        pid = elem.property_id
        prop = model.ptube.slice_card_by_property_id(pid)
    else:  # pragma: no cover
        raise NotImplementedError(element_name)

    mat1 = model.mat1.slice_card_by_material_id(prop.material_id)
    E = mat1.E
    G = mat1.G
    if prop.type == 'PTUBE':
        A = prop.area()
        J = prop.J()
    else:
        A = prop.A
        J = prop.J

    assert isinstance(J, np.ndarray), (prop.type, J)
    assert isinstance(A, np.ndarray), (prop.type, A)
    assert isinstance(E, np.ndarray), (prop.type, E)
    for ieid, eid, nodes, xyz1i, xyz2i, Gi, Ji, Ai, Ei in zip(
            ieids, eids, elem.nodes, xyz1, xyz2, G, J, A, E):
        forces[ieid, :] = _recover_forcei_rod(
            xb, dof_map, nodes, xyz1i, xyz2i, Gi, Ji, Ai, Ei)

    data = forces.reshape(1, *forces.shape)
    table_name = 'OEF1'
    force_obj = RealRodForceArray.add_static_case(
        table_name, element_name, eids, data, isubcase,
        is_sort1=True, is_random=False, is_msc=True,
        random_code=0, title=title, subtitle=subtitle, label=label)

    force = op2.op2_results.force
    if element_name == 'CONROD':
        force.conrod_force[isubcase] = force_obj
    elif element_name == 'CROD':
        force.crod_force[isubcase] = force_obj
    elif element_name == 'CTUBE':
        force.ctube_force[isubcase] = force_obj
    else:  # pragma: no cover
        raise NotImplementedError(element_name)

    force_obj.write_f06(f06_file, header=None, page_stamp=page_stamp,
                        page_num=page_num, is_mag_phase=False, is_sort1=True)
    return neids

def _recover_forcei_rod(xb: np.ndarray,
                        dof_map: DOF_MAP,
                        nodes: np.ndarray,
                        xyz1: np.ndarray, xyz2: np.ndarray,
                        G, J, A, E):
    """get the static rod force"""
    nid1, nid2 = nodes
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
    dxyz12 = xyz1 - xyz2
    Lambda = lambda1d(dxyz12, debug=False)

    u_axial = Lambda @ q_axial
    u_torsion = Lambda @ q_torsion
    du_axial = u_axial[0] - u_axial[1]
    du_torsion = u_torsion[0] - u_torsion[1]
    #headers = ['axial', 'SMa', 'torsion', 'SMt']
    #C = prop.c
    L = np.linalg.norm(dxyz12)
    axial_strain = du_axial / L
    #torsional_strain = du_torsion * C / L

    axial_stress = E * axial_strain
    #torsional_stress = G * torsional_strain
    axial_force = axial_stress * A
    torsional_moment = du_torsion * G * J / L
    return axial_force, torsional_moment

def _recover_force_cbar(f06_file: TextIO, op2,
                        model: BDF, dof_map, isubcase, xb, eids_str,
                        element_name, fdtype='float32',
                        title: str='', subtitle: str='', label: str='',
                        page_num: int=1, page_stamp='PAGE %s') -> int:
    """
    Recovers static CBAR force.

    .. todo:: doesn't support CBAR-100

    """
    neids, ieids, eids = get_ieids_eids(model, element_name, eids_str)
    if not neids:
        return neids

    elem = get_element(model, element_name, ieids, eids)
    xyz1 = model.grid.get_position_by_node_id(elem.nodes[:, 0])
    xyz2 = model.grid.get_position_by_node_id(elem.nodes[:, 1])

    # mat1 = model.mat1.slice_card_by_material_id(prop.material_id)
    # E = mat1.E
    # G = mat1.G
    # if prop.type == 'PTUBE':
    #     A = prop.area()
    #     J = prop.J()
    # else:
    #     A = prop.A
    #     J = prop.J

    AIJEG = elem.stiffness_info()
    # columns: [length, area, I1, I2, I12, J, E, G]
    A = AIJEG[:, 1]
    I = AIJEG[:, [2, 3, 4]]
    J = AIJEG[:, 5]
    E = AIJEG[:, 6]
    G = AIJEG[:, 7]
    assert isinstance(J, np.ndarray), (elem.type, J)
    assert isinstance(A, np.ndarray), (elem.type, A)
    assert isinstance(E, np.ndarray), (elem.type, E)

    forces = np.full((neids, 8), np.nan, dtype=fdtype)

    v, ihat, yhat, zhat, wa, wb = elem.get_axes(xyz1, xyz2)
    for (ieid, eid, nodes, xyz1i, xyz2i,
         Ai, I1, Ji, Ei, Gi,
         vi, ihati, yhati, zhati, wai, wbi) in zip(
            ieids, eids, elem.nodes, xyz1, xyz2, A, I, J, E, G,
            v, ihat, yhat, zhat, wa, wb):
        forces[ieid, :] = _recover_forcei_cbar(model, xb, dof_map, nodes,
                                               xyz1i, xyz2i, Ai, I1, Ji, Ei, Gi,
                                               vi, ihati, yhati, zhati, wai, wbi)

    data = forces.reshape(1, *forces.shape)
    table_name = 'OEF1'
    force_obj = RealCBarForceArray.add_static_case(
        table_name, 'CBAR', eids, data, isubcase,
        is_sort1=True, is_random=False, is_msc=True,
        random_code=0, title=title, subtitle=subtitle, label=label)

    force = op2.op2_results.force
    force.cbar_force[isubcase] = force_obj

    force_obj.write_f06(f06_file, header=None, page_stamp=page_stamp,
                        page_num=page_num, is_mag_phase=False, is_sort1=True)
    return len(elem)

def ke_cbar(xyz1, xyz2,
            A, I, J, E, G,
            ihat, jhat, khat, wa, wb,
            fdtype: str='float64'):
    """Get the elemental stiffness matrix in the basic frame."""
    I1, I2, I12 = I
    dxyz = xyz2 - xyz1
    L = np.linalg.norm(dxyz)
    k1 = k2 = 1e8
    Ke = timoshenko_stiffness(A, E, G, L, I1, I2, J, k1, k2, pa=0, pb=0)
    Teb = beam_transform(ihat, jhat, khat)
    K = Teb.T @ Ke @ Teb
    return True, K

def _recover_forcei_cbar(model: BDF,
                         xb: np.ndarray,
                         dof_map: DOF_MAP,
                         nodes: np.ndarray,
                         xyz1: np.ndarray,
                         xyz2: np.ndarray,
                         A, I, J, E, G,
                         v, ihat, jhat, khat, wa, wb,
                         fdtype: str='float64'):
    """Get the static CBAR force."""
    nid1, nid2 = nodes
    I1, I2, I12 = I

    i1 = dof_map[(nid1, 1)]
    i2 = dof_map[(nid2, 1)]

    q_all = np.hstack([xb[i1:i1 + 6], xb[i2:i2 + 6]])

    dxyz = xyz2 - xyz1
    L = np.linalg.norm(dxyz)
    k1 = k2 = 1e8
    Ke = timoshenko_stiffness(A, E, G, L, I1, I2, J, k1, k2, pa=0, pb=0)
    Teb = beam_transform(ihat, jhat, khat)
    Fe = recover_beam_force(Ke, Teb, q_all)

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
