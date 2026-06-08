from __future__ import annotations
from typing import TextIO, TYPE_CHECKING
import numpy as np

from pyNastran.dev.solver.utils import lambda1d
from pyNastran.dev.bdf_vectorized3.solver.utils import get_ieids_eids, get_element
from pyNastran.dev.bdf_vectorized3.solver.elements.beam import (
    timoshenko_stiffness,
    beam_transform,
)
from pyNastran.op2.op2_interface.op2_classes import (
    RealStrainEnergyArray,
)
from .static_spring import _recover_strain_energy_celas
from .utils import get_plot_request

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF, Subcase

    DOF_MAP = dict[tuple[int, int], int]


def recover_strain_energy_101(f06_file: TextIO, op2,
                              model: BDF, dof_map,
                              subcase: Subcase, xb, fdtype: str='float32',
                              title: str='', subtitle: str='', label: str='',
                              page_num: int=1, page_stamp: str='PAGE %s'):
    """Recovers element strain energy from ESE = ALL."""
    if 'ESE' not in subcase:
        return
    eid_str = 'ALL'
    unused_eids_write, write_f06, write_op2, quick_return = get_plot_request(
        subcase, 'ESE')
    if quick_return:
        return page_num
    isubcase = subcase.id

    nelements = 0
    nelements += _recover_strain_energy_celas(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CELAS1', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += _recover_strain_energy_celas(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CELAS2', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += _recover_strain_energy_celas(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CELAS3', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += _recover_strain_energy_celas(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CELAS4', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)

    nelements += _recover_strain_energy_rod(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CROD', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += _recover_strain_energy_rod(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CONROD', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += _recover_strain_energy_rod(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CTUBE', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += _recover_strain_energy_beam(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CBAR', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += _recover_strain_energy_beam(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CBEAM', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    if nelements == 0:
        model.log.warning(f'no strain energy output...{model.card_count}; {model.bdf_filename}')


def _recover_strain_energy_rod(
    f06_file: TextIO, op2,
    model: BDF, dof_map: DOF_MAP,
    isubcase: int, xb: np.ndarray, eids_str: str,
    element_name: str, fdtype: str = 'float32',
    title: str = '', subtitle: str = '', label: str = '',
    page_num: int = 1, page_stamp: str = 'PAGE %s',
) -> int:
    """Recover strain energy for CROD/CONROD/CTUBE: SE = 0.5 * u^T @ K @ u."""
    neids, ieids, eids = get_ieids_eids(model, element_name, eids_str)
    if not neids:
        return 0

    elem = get_element(model, element_name, ieids, eids)
    xyz1 = model.grid.get_position_by_node_id(elem.nodes[:, 0])
    xyz2 = model.grid.get_position_by_node_id(elem.nodes[:, 1])

    if element_name == 'CONROD':
        prop = elem
    elif element_name == 'CROD':
        pid = elem.property_id
        prop = model.prod.slice_card_by_property_id(pid)
    elif element_name == 'CTUBE':
        pid = elem.property_id
        prop = model.ptube.slice_card_by_property_id(pid)
    else:
        raise NotImplementedError(element_name)

    mat1 = model.mat1.slice_card_by_material_id(prop.material_id)
    E = mat1.E
    G = mat1.G
    if hasattr(prop, 'area'):
        A = prop.area()
    else:
        A = prop.A
    if hasattr(prop, 'J'):
        if callable(prop.J):
            J_arr = prop.J()
        else:
            J_arr = prop.J
    else:
        J_arr = np.zeros(neids)

    strain_energies = np.full((neids, 1), np.nan, dtype=fdtype)

    for ieid, eid, nodes, xyz1i, xyz2i, Ei, Gi, Ai, Ji in zip(
        ieids, eids, elem.nodes, xyz1, xyz2, E, G, A, J_arr,
    ):
        nid1, nid2 = nodes
        dxyz = xyz2i - xyz1i
        L = np.linalg.norm(dxyz)

        i1 = dof_map[(nid1, 1)]
        i2 = dof_map[(nid2, 1)]

        # Rod element: axial (k=EA/L) and torsion (k=GJ/L)
        # SE = 0.5 * k * (du)^2
        Lambda = lambda1d(dxyz, debug=False)
        q_axial = np.array([
            xb[i1], xb[i1 + 1], xb[i1 + 2],
            xb[i2], xb[i2 + 1], xb[i2 + 2],
        ])
        q_torsion = np.array([
            xb[i1 + 3], xb[i1 + 4], xb[i1 + 5],
            xb[i2 + 3], xb[i2 + 4], xb[i2 + 5],
        ])
        u_axial = Lambda @ q_axial
        u_torsion = Lambda @ q_torsion
        du_axial = u_axial[0] - u_axial[1]
        du_torsion = u_torsion[0] - u_torsion[1]

        k_axial = Ei * Ai / L
        k_torsion = Gi * Ji / L if Ji > 0 else 0.0
        se = 0.5 * k_axial * du_axial**2 + 0.5 * k_torsion * du_torsion**2
        strain_energies[ieid] = se

    _save_strain_energy(
        op2, f06_file, page_num, page_stamp, element_name,
        strain_energies, eids, isubcase, title, subtitle, label,
    )
    return neids


def _recover_strain_energy_beam(
    f06_file: TextIO, op2,
    model: BDF, dof_map: DOF_MAP,
    isubcase: int, xb: np.ndarray, eids_str: str,
    element_name: str, fdtype: str = 'float32',
    title: str = '', subtitle: str = '', label: str = '',
    page_num: int = 1, page_stamp: str = 'PAGE %s',
) -> int:
    """Recover strain energy for CBAR/CBEAM: SE = 0.5 * u^T @ K @ u."""
    neids, ieids, eids = get_ieids_eids(model, element_name, eids_str)
    if not neids:
        return 0

    elem = get_element(model, element_name, ieids, eids)
    xyz1 = model.grid.get_position_by_node_id(elem.nodes[:, 0])
    xyz2 = model.grid.get_position_by_node_id(elem.nodes[:, 1])

    AIJEG = elem.stiffness_info()
    # columns: [length, area, I1, I2, I12, J, E, G]
    Avec = AIJEG[:, 1]
    Ivec = AIJEG[:, [2, 3, 4]]
    Jvec = AIJEG[:, 5]
    Evec = AIJEG[:, 6]
    Gvec = AIJEG[:, 7]

    v, ihat, yhat, zhat, wa, wb = elem.get_axes(xyz1, xyz2)
    k_arr = elem.k()

    strain_energies = np.full((neids, 1), np.nan, dtype=fdtype)

    for ieid, eid, nodes, xyz1i, xyz2i, Ai, Ii, Ji, Ei, Gi, ihati, jhati, khati, ki in zip(
        ieids, eids, elem.nodes, xyz1, xyz2, Avec, Ivec, Jvec, Evec, Gvec,
        ihat, yhat, zhat, k_arr,
    ):
        nid1, nid2 = nodes
        I1, I2, I12 = Ii
        k1, k2 = ki
        L = np.linalg.norm(xyz2i - xyz1i)

        Ke = timoshenko_stiffness(Ai, Ei, Gi, L, I1, I2, Ji, k1, k2)
        Teb = beam_transform(ihati, jhati, khati)

        gi1 = dof_map[(nid1, 1)]
        gi2 = dof_map[(nid2, 1)]
        q_basic = np.hstack([xb[gi1:gi1 + 6], xb[gi2:gi2 + 6]])
        q_elem = Teb @ q_basic
        se = 0.5 * q_elem @ Ke @ q_elem
        strain_energies[ieid] = se

    _save_strain_energy(
        op2, f06_file, page_num, page_stamp, element_name,
        strain_energies, eids, isubcase, title, subtitle, label,
    )
    return neids


def _save_strain_energy(
    op2, f06_file, page_num, page_stamp,
    element_name: str,
    strain_energies: np.ndarray,
    eids: np.ndarray,
    isubcase: int,
    title: str, subtitle: str, label: str,
) -> None:
    """Save strain energy results to OP2 and write F06."""
    if strain_energies.sum() == 0.0:
        return
    data = strain_energies.reshape(1, *strain_energies.shape)
    table_name = 'ONRGY1'
    try:
        se_obj = RealStrainEnergyArray.add_static_case(
            table_name, element_name, eids, data, isubcase,
            is_sort1=True, is_random=False, is_msc=True,
            random_code=0, title=title, subtitle=subtitle, label=label)
    except (KeyError, NotImplementedError, AssertionError):
        return

    ese = op2.op2_results.strain_energy
    attr = f'{element_name.lower()}_strain_energy'
    if hasattr(ese, attr):
        getattr(ese, attr)[isubcase] = se_obj

    se_obj.write_f06(
        f06_file, header=None, page_stamp=page_stamp,
        page_num=page_num, is_mag_phase=False, is_sort1=True,
    )
