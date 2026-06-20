from __future__ import annotations
from typing import TextIO, Optional, TYPE_CHECKING
import numpy as np

from pyNastran.dev.bdf_vectorized3.solver.utils import (
    get_ieids_eids, get_element, lambda1d)
from pyNastran.op2.op2_interface.op2_classes import (
    RealRodForceArray,
)
from .utils import get_plot_request, save_strain_energy

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf_vectorized3.bdf import BDF, Subcase
DOF_MAP = dict[tuple[int, int], int]


def get_rod_pid_prop(model: BDF,
                     elem,
                     element_name: str) -> tuple:
    if element_name == 'CONROD':
        prop = elem
        pid = None
        A = prop.A
        J = prop.J
    elif element_name == 'CROD':
        pid = elem.property_id
        prop = model.prod.slice_card_by_property_id(pid)
        A = prop.A
        J = prop.J
    elif element_name == 'CTUBE':
        pid = elem.property_id
        prop = model.ptube.slice_card_by_property_id(pid)
        A = prop.area()
        J = prop.J()
    else:  # pragma: no cover
        raise NotImplementedError(element_name)

    assert isinstance(J, np.ndarray), (prop.type, J)
    assert isinstance(A, np.ndarray), (prop.type, A)

    mat1 = elem.model.mat1.slice_card_by_material_id(prop.material_id)
    E = mat1.E
    G = mat1.G
    assert isinstance(E, np.ndarray), (prop.type, E)
    return prop, pid, A, J, E, G

def get_rod_du(xb: np.ndarray,
               dof_map: DOF_MAP,
               elem, ieids, eids, dxyz, fdtype='float64'):
    neids = len(eids)
    du_axial = np.full(neids, np.nan, dtype=fdtype)
    du_torsion = np.full(neids, np.nan, dtype=fdtype)
    for ieid, eid, nodes, dxyzi in zip(ieids, eids, elem.nodes, dxyz):
        du_axiali, du_torsioni = _recover_dx_rod(
            xb, dof_map, nodes, dxyzi)
        du_axial[ieid] = du_axiali
        du_torsion[ieid] = du_torsioni
    return du_axial, du_torsion

def _recover_dx_rod(xb: np.ndarray,
                    dof_map: DOF_MAP,
                    nodes: np.ndarray,
                    dxyz: np.ndarray):
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
    Lambda = lambda1d(dxyz, debug=False)

    u_axial = Lambda @ q_axial
    u_torsion = Lambda @ q_torsion
    du_axial = u_axial[0] - u_axial[1]
    du_torsion = u_torsion[0] - u_torsion[1]
    return du_axial, du_torsion


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
    dxyz = xyz2 - xyz1
    L = np.linalg.norm(dxyz, axis=1)
    assert len(L) == neids

    _prop, _pid, A, J, E, G = get_rod_pid_prop(model, elem, element_name)
    du_axial, du_torsion = get_rod_du(
        xb, dof_map, elem, ieids, eids, dxyz, fdtype=fdtype)

    axial_force = A * E * du_axial / L
    torsional_moment = du_torsion * G * J / L

    forces = np.column_stack([axial_force, torsional_moment])

    data = forces.reshape(1, *forces.shape)
    table_name = 'OEF1'
    force_obj = RealRodForceArray.add_static_case(
        table_name, element_name, eids, data, isubcase,
        is_sort1=True, is_random=False, is_msc=True,
        random_code=0, title=title, subtitle=subtitle, label=label)

    force = op2.op2_results.force
    etype = element_name.lower()
    obj_dict = getattr(force, f'{etype}_force')
    obj_dict[isubcase] = force_obj

    force_obj.write_f06(f06_file, header=None, page_stamp=page_stamp,
                        page_num=page_num, is_mag_phase=False, is_sort1=True)
    return neids


def _recover_stress_rod(
    f06_file: TextIO,
    op2,
    model: BDF,
    dof_map: DOF_MAP,
    isubcase: int,
    xb: np.ndarray,
    eids_str: str,
    element_name: str,
    fdtype: str = "float32",
    title: str = "",
    subtitle: str = "",
    label: str = "",
    page_num: int = 1,
    page_stamp: str = "PAGE %s",) -> int:
    """Recovers static rod stress."""
    neids, ieids, eids = get_ieids_eids(model, element_name, eids_str)
    if not neids:
        return neids

    elem = get_element(model, element_name, ieids, eids)
    xyz1 = model.grid.get_position_by_node_id(elem.nodes[:, 0])
    xyz2 = model.grid.get_position_by_node_id(elem.nodes[:, 1])

    dxyz = xyz2 - xyz1
    L = np.linalg.norm(dxyz, axis=1)
    assert len(L) == neids

    prop, pid, A, J, E, G = get_rod_pid_prop(model, elem, element_name)
    du_axial, du_torsion = get_rod_du(
        xb, dof_map, elem, ieids, eids, dxyz, fdtype=fdtype)

    #axial_strain = np.zeros(neids, dtype=dtype)
    #torsional_moment = np.zeros(neids, dtype=dtype)
    #torsional_stress = np.zeros(neids, dtype=dtype)
    # jpos = (J > 0)

    axial_strain = du_axial / L
    torsional_strain = du_torsion / L
    
    # stress1 = E*e11

    # F = kx = kaxial*du = EA/L * du = EA * e11
    # e11 = du/L
    axial_force = E * A * axial_strain

    # stress = F/A = E*e11
    axial_stress = E * axial_strain
    
    # M = k*theta = ktor*du = GJ/L * du = GJ * e44
    torsional_moment = torsional_strain * G * J

    # stress = M*r/J = GJ * e44 * r / J
    #        = G * e44 * r
    radius = 1.0
    torsional_stress = G * torsional_strain * radius

    #headers = ['axial', 'SMa', 'torsion', 'SMt']
    # SM fields are NaN (safety margin not computed)
    nan = np.full(neids, np.nan, dtype=fdtype)
    stresses = np.column_stack([axial_stress, nan, torsional_stress, nan])

    data = stresses.reshape(1, *stresses.shape)
    table_name = "OES1"
    stress_obj = RealRodStressArray.add_static_case(
        table_name,
        element_name,
        eids,
        data,
        isubcase,
        is_sort1=True,
        is_random=False,
        is_msc=True,
        random_code=0,
        title=title,
        subtitle=subtitle,
        label=label,
    )

    stress = op2.op2_results.stress
    etype = element_name.lower()
    obj_dict = getattr(stress, f'{etype}_stress')
    obj_dict[isubcase] = stress_obj

    stress_obj.write_f06(
        f06_file,
        header=None,
        page_stamp=page_stamp,
        page_num=page_num,
        is_mag_phase=False,
        is_sort1=True,
    )
    return neids

def _recover_strain_rod(
    f06_file: TextIO,
    op2,
    model: BDF,
    dof_map: DOF_MAP,
    isubcase: int,
    xb: np.ndarray,
    eids_str: str,
    element_name: str,
    fdtype: str = "float32",
    title: str = "",
    subtitle: str = "",
    label: str = "",
    page_num: int = 1,
    page_stamp: str = "PAGE %s",
) -> int:
    """Recovers static rod strain."""
    neids, ieids, eids = get_ieids_eids(model, element_name, eids_str)
    if not neids:
        return neids

    elem = get_element(model, element_name, ieids, eids)
    xyz1 = model.grid.get_position_by_node_id(elem.nodes[:, 0])
    xyz2 = model.grid.get_position_by_node_id(elem.nodes[:, 1])

    dxyz = xyz2 - xyz1
    L = np.linalg.norm(dxyz, axis=1)
    assert len(L) == neids

    _prop, _pid, A, J, E, G = get_rod_pid_prop(model, elem, element_name)
    du_axial, du_torsion = get_rod_du(
        xb, dof_map, elem, ieids, eids, dxyz, fdtype=fdtype)

    axial_strain = du_axial / L
    torsional_strain = du_torsion / L

    # SM fields are NaN (safety margin not computed)
    #return axial_strain, np.nan, torsional_strain, np.nan

    #headers = ['axial', 'SMa', 'torsion', 'SMt']
    # SM fields are NaN (safety margin not computed)
    nan = np.full(neids, np.nan, dtype=fdtype)
    strains = np.column_stack([axial_strain, nan, torsional_strain, nan])

    data = strains.reshape(1, *strains.shape)
    table_name = "OSTR1"
    strain_obj = RealRodStrainArray.add_static_case(
        table_name,
        element_name,
        eids,
        data,
        isubcase,
        is_sort1=True,
        is_random=False,
        is_msc=True,
        random_code=0,
        title=title,
        subtitle=subtitle,
        label=label,
    )

    strain = op2.op2_results.strain
    etype = element_name.lower()
    obj_dict = getattr(strain, f'{etype}_strain')
    obj_dict[isubcase] = strain_obj

    strain_obj.write_f06(
        f06_file,
        header=None,
        page_stamp=page_stamp,
        page_num=page_num,
        is_mag_phase=False,
        is_sort1=True,
    )
    return neids

def _recover_strain_energy_rod(
    f06_file: TextIO, op2,
    model: BDF, dof_map: DOF_MAP,
    isubcase: int, xb: np.ndarray, eids_str: str,
    element_name: str, fdtype: str = 'float32',
    title: str = '', subtitle: str = '', label: str = '',
    page_num: int = 1, page_stamp: str = 'PAGE %s',) -> int:
    """Recover strain energy for CROD/CONROD/CTUBE: SE = 0.5 * u^T @ K @ u."""
    neids, ieids, eids = get_ieids_eids(model, element_name, eids_str)
    if not neids:
        return 0

    elem = get_element(model, element_name, ieids, eids)
    xyz1 = model.grid.get_position_by_node_id(elem.nodes[:, 0])
    xyz2 = model.grid.get_position_by_node_id(elem.nodes[:, 1])
    dxyz = xyz2 - xyz1
    L = np.linalg.norm(dxyz, axis=1)
    assert len(L) == neids

    _prop, _pid, A, J, E, G = get_rod_pid_prop(model, elem, element_name)
    du_axial, du_torsion = get_rod_du(
        xb, dof_map, elem, ieids, eids, dxyz, fdtype=fdtype)

    k_axial = E * A / L
    k_torsion = G * J / L
    strain_energy = 0.5 * (k_axial * du_axial**2 + k_torsion * du_torsion**2)

    save_strain_energy(
        op2, f06_file, page_num, page_stamp, element_name,
        strain_energy, eids, isubcase, title, subtitle, label,
    )
    return neids
