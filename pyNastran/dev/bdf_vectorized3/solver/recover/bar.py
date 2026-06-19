from __future__ import annotations
from typing import TextIO, Optional, TYPE_CHECKING
import numpy as np

from pyNastran.dev.bdf_vectorized3.solver.utils import (
    get_ieids_eids, get_element, lambda1d)
from pyNastran.op2.op2_interface.op2_classes import (
    RealCBarForceArray,
)
# from pyNastran.dev.solver.build_stiffness import ke_cbar

from pyNastran.dev.bdf_vectorized3.solver.elements.beam import (
    timoshenko_stiffness, beam_transform, recover_beam_force,
)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf_vectorized3.bdf import BDF, Subcase
    from pyNastran.bdf_vectorized3.cards import CBAR, PBAR, PBARL
DOF_MAP = dict[tuple[int, int], int]


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
