from __future__ import annotations
from typing import TextIO, TYPE_CHECKING
import numpy as np

from pyNastran.dev.bdf_vectorized3.solver.utils import (
    get_ieids_eids, get_element, lambda1d)
from pyNastran.dev.bdf_vectorized3.solver.elements.beam import (
    timoshenko_stiffness,
    beam_transform,
    beam_transforms,
)
from pyNastran.op2.op2_interface.op2_classes import (
    RealStrainEnergyArray,)
from .static_spring import _recover_strain_energy_celas
from .rod import _recover_strain_energy_rod
#from .bar import _recover_strain_energy_beam
from .utils import get_plot_request, save_strain_energy, fix_xb_shape


if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.bdf import BDF, Subcase
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
    for name in ['CELAS1', 'CELAS2', 'CELAS3', 'CELAS4']:
        nelements += _recover_strain_energy_celas(
            f06_file, op2, model, dof_map, isubcase, xb, eid_str,
            name, fdtype=fdtype,
            title=title, subtitle=subtitle, label=label,
            page_num=page_num, page_stamp=page_stamp)

    for name in ['CROD', 'CTUBE', 'CONROD']:
        nelements += _recover_strain_energy_rod(
            f06_file, op2, model, dof_map, isubcase, xb, eid_str,
            name, fdtype=fdtype,
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


def _recover_strain_energy_beam(
        f06_file: TextIO,
        op2: OP2,
        model: BDF,
        dof_map: DOF_MAP,
        isubcase: int, xb: np.ndarray, eids_str: str,
        element_name: str, fdtype: str = 'float32',
        title: str = '', subtitle: str = '', label: str = '',
        page_num: int = 1,
        page_stamp: str = 'PAGE %s',) -> int:
    """Recover strain energy for CBAR/CBEAM: SE = 0.5 * u^T @ K @ u."""
    neids, ieids, eids = get_ieids_eids(model, element_name, eids_str)
    if not neids:
        return 0

    xb, nmode = fix_xb_shape(xb)

    elem = get_element(model, element_name, ieids, eids)
    xyz1 = model.grid.get_position_by_node_id(elem.nodes[:, 0])
    xyz2 = model.grid.get_position_by_node_id(elem.nodes[:, 1])

    LAIJEG = elem.stiffness_info()
    # columns: [length, area, I1, I2, I12, J, E, G]
    L = LAIJEG[:, 0]
    A = LAIJEG[:, 1]
    I = LAIJEG[:, [2, 3, 4]]
    J = LAIJEG[:, 5]
    E = LAIJEG[:, 6]
    G = LAIJEG[:, 7]

    v, ihat, yhat, zhat, wa, wb = elem.get_axes(xyz1, xyz2)
    ks = elem.k()

    strain_energy2 = np.full((neids, 1), np.nan, dtype=fdtype)

    Tebs = beam_transforms(ihat, yhat, zhat)
    Ke = np.zeros((neids, 12, 12), dtype=fdtype)
    q_basic = np.zeros((neids, 12), dtype=fdtype)
    q_elem2 = np.zeros((neids, 12), dtype=fdtype)
    for (ieid, nodes, Li, Ai, Ii, Ji, Ei, Gi, ki, Tebi,
        ) in zip(ieids, elem.nodes, L, A, I, J, E, G, ks, Tebs):
        nid1, nid2 = nodes
        I1, I2, I12 = Ii
        k1, k2 = ki
   
        Kei = timoshenko_stiffness(Ai, Ei, Gi, Li, I1, I2, Ji, k1, k2)
        Ke[ieid, :, :] = Kei

        gi1 = dof_map[(nid1, 1)]
        gi2 = dof_map[(nid2, 1)]
        q_basici = np.hstack([xb[gi1:gi1 + 6], xb[gi2:gi2 + 6]])
        q_basic[ieid, :] = q_basici
        #q_elemi = Teb2 @ q_basici
        #q_elem2[ieid, :] = q_elemi
        #strain_energy2[ieid] = 0.5 * q_elemi @ Kei @ q_elemi
    q_elem = np.einsum('ntu,nu->nt', Tebs, q_basic)
    strain_energy = 0.5 * np.einsum('nt,ntu,nu->n', q_elem, Ke, q_elem)
    #assert np.allclose(q_elem, q_elem2)
    #assert np.allclose(strain_energy, strain_energy2)

    strain_energy = strain_energy.reshape(nmode, neids, 1)
    save_strain_energy(
        op2, f06_file, page_num, page_stamp, element_name,
        strain_energy, eids, isubcase, title, subtitle, label)
    return neids
