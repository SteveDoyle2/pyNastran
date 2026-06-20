from __future__ import annotations
from typing import TextIO, Optional, TYPE_CHECKING
import numpy as np

from pyNastran.dev.bdf_vectorized3.cards.base_card import searchsorted_filter
from pyNastran.dev.bdf_vectorized3.solver.utils import (
    get_ieids_eids, get_element)
from pyNastran.op2.op2_interface.op2_classes import (
    RealCBarForceArray,
    RealBarStrainArray,
    RealBarStressArray,
)
from .utils import fix_xb_shape
from .beam import beam_stress_at_points, _set_min_max_stress

#from .static_stress import _get_bar_recovery_points
from pyNastran.dev.bdf_vectorized3.solver.elements.beam import (
    timoshenko_stiffness,
    beam_transforms,
)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf_vectorized3.bdf import BDF, Subcase
    from pyNastran.bdf_vectorized3.cards import CBAR, PBAR, PBARL
DOF_MAP = dict[tuple[int, int], int]


def _get_bar_recovery_points(model: BDF, elem) -> np.ndarray:
    """Get C,D,E,F stress recovery points per element from properties.

    Returns
    -------
    cdef : np.ndarray, shape (neids, 8)
        [C1, C2, D1, D2, E1, E2, F1, F2] per element.
    """
    neids = elem.n
    cdef = np.zeros((neids, 8), dtype="float64")
    pids = elem.property_id

    for prop in elem.allowed_properties:
        i_lookup, i_all = searchsorted_filter(prop.property_id, pids, msg='')
        if len(i_lookup) == 0:
            continue
        cdefi = np.column_stack([
            prop.c[i_all,:], prop.d[i_all,:],
            prop.e[i_all,:], prop.f[i_all,:]])
        cdef[i_lookup, :] = cdefi
    return cdef


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

    xb, nmode = fix_xb_shape(xb)
    elem = get_element(model, element_name, ieids, eids)

    LAIJEG = elem.stiffness_info()
    # columns: [length, area, I1, I2, I12, J, k1, k2, E, G]
    L = LAIJEG[:, 0]
    A = LAIJEG[:, 1]
    I = LAIJEG[:, [2, 3, 4]]
    J = LAIJEG[:, 5]
    #k1 = LAIJEG[:, 6] 
    #k2 = LAIJEG[:, 7]
    E = LAIJEG[:, 8]
    G = LAIJEG[:, 9]
    assert isinstance(J, np.ndarray), (elem.type, J)
    assert isinstance(A, np.ndarray), (elem.type, A)
    assert isinstance(E, np.ndarray), (elem.type, E)

    force2, Fe = get_bar_force_fe(
        model, xb, dof_map, elem, ieids, neids,
        LAIJEG, fdtype=fdtype)

    Fx1 = Fe[:, 0]
    Fy1 = Fe[:, 1]
    Fz1 = Fe[:, 2]
    Mx1 = Fe[:, 3]
    My1 = Fe[:, 4]
    Mz1 = Fe[:, 5]

    Fx2 = Fe[:, 6]
    Fy2 = Fe[:, 7]
    Fz2 = Fe[:, 8]
    #Mx2 = Fe[:, 9]
    My2 = Fe[:, 10]
    Mz2 = Fe[:, 11]

    axial = Fx1
    torque = Mx1
    shear1 = Fy1
    shear2 = Fz1
    bending_moment_a1 = My1
    bending_moment_a2 = Mz1

    bending_moment_b1 = My2
    bending_moment_b2 = Mz2

    force = np.column_stack([
        bending_moment_a1, bending_moment_a2,
        bending_moment_b1, bending_moment_b2,
        shear1, shear2, axial, torque])
    assert np.allclose(force, force2)

    data = force.reshape(nmode, neids, 8)
    table_name = 'OEF1'
    force_obj = RealCBarForceArray.add_static_case(
        table_name, 'CBAR', eids, data, isubcase,
        is_sort1=True, is_random=False, is_msc=True,
        random_code=0, title=title, subtitle=subtitle, label=label)

    force_dict = op2.op2_results.force
    force_dict.cbar_force[isubcase] = force_obj

    force_obj.write_f06(f06_file, header=None, page_stamp=page_stamp,
                        page_num=page_num, is_mag_phase=False, is_sort1=True)
    return len(elem)

def get_bar_force_fe(model: BDF,
                     xb: np.ndarray,
                     dof_map: DOF_MAP,
                     elem: CBAR,
                     ieids, neids: int,
                     LAIJEG: np.ndarray, fdtype='float32'):
    xyz1 = model.grid.get_position_by_node_id(elem.nodes[:, 0])
    xyz2 = model.grid.get_position_by_node_id(elem.nodes[:, 1])

    L = LAIJEG[:, 0]
    A = LAIJEG[:, 1]
    I = LAIJEG[:, [2, 3, 4]]
    J = LAIJEG[:, 5]
    k1 = LAIJEG[:, 6] 
    k2 = LAIJEG[:, 7]
    E = LAIJEG[:, 8]
    G = LAIJEG[:, 9]

    forces = np.full((neids, 8), np.nan, dtype=fdtype)
    Fe = np.full((neids, 12), np.nan, dtype=fdtype)
    Ke = np.full((neids, 12, 12), np.nan, dtype=fdtype)
    q = np.full((neids, 12), np.nan, dtype=fdtype)
    v, ihat, yhat, zhat, wa, wb = elem.get_axes(xyz1, xyz2)
    Teb = beam_transforms(ihat, yhat, zhat)
    for (ieid, nodes, Li, Ai, I1, Ji, Ei, Gi,
         vi, Tebi, wai, wbi, k1i, k2i) in zip(
            ieids, elem.nodes, L, A, I, J, E, G,
            v, Teb, wa, wb, k1, k2):
        qi, Kei, Fei = _recover_forcei_cbar(
            model, xb, dof_map, nodes,
            Li, Ai, I1, Ji, Ei, Gi,
            vi, Tebi, wai, wbi,
            k1i, k2i)

        (Fx1, Fy1, Fz1, Mx1, My1, Mz1,
         Fx2, Fy2, Fz2, Mx2, My2, Mz2) = Fei
        axial = Fx1
        torque = Mx1
        shear1 = Fy1
        shear2 = Fz1
        bending_moment_a1 = My1
        bending_moment_a2 = Mz1

        bending_moment_b1 = My2
        bending_moment_b2 = Mz2

        force = (
            bending_moment_a1, bending_moment_a2,
            bending_moment_b1, bending_moment_b2,
            shear1, shear2,
            axial, torque)
        q[ieid, :] = qi
        Ke[ieid, :, :] = Kei
        forces[ieid, :] = force
        Fe[ieid, :] = Fei

    #Fe2 = Ke @ Teb @ q
    Fe2 = np.einsum('ntu,nuv,nv->nt', Ke, Teb, q)
    assert np.allclose(Fe, Fe2)
    return forces, Fe


#def ke_cbar(L, A, I, J, E, G,
#            Teb, wa, wb, k1, k2,
#            fdtype: str='float64'):
#    """Get the elemental stiffness matrix in the basic frame."""
#    I1, I2, I12 = I
#    Ke = timoshenko_stiffness(A, E, G, L, I1, I2, J, k1, k2,
#                              pa=0, pb=0)
#    K = Teb.T @ Ke @ Teb
#    return True, K

def _recover_forcei_cbar(model: BDF,
                         xb: np.ndarray,
                         dof_map: DOF_MAP,
                         nodes: np.ndarray,
                         L, A, I, J, E, G,
                         v, Teb, wa, wb,
                         k1, k2,
                         fdtype: str='float64'):
    """Get the static CBAR force."""
    nid1, nid2 = nodes
    I1, I2, I12 = I

    i1 = dof_map[(nid1, 1)]
    i2 = dof_map[(nid2, 1)]

    q_all = np.hstack([xb[i1:i1 + 6], xb[i2:i2 + 6]])

    Ke = timoshenko_stiffness(A, E, G, L, I1, I2, J, k1, k2,
                              pa=0, pb=0)
    Fe = Ke @ Teb @ q_all  # recover beam force

    #(Fx1, Fy1, Fz1, Mx1, My1, Mz1,
    # Fx2, Fy2, Fz2, Mx2, My2, Mz2) = Fe
    return q_all, Ke, Fe


def _recover_strain_cbar(
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
    """Recovers static CBAR strain."""
    neids, ieids, eids = get_ieids_eids(model, element_name, eids_str)
    if not neids:
        return neids

    elem = get_element(model, element_name, ieids, eids)

    LAIJEG = elem.stiffness_info()
    forces, Fe = get_bar_force_fe(
        model, xb, dof_map, elem, ieids, neids,
        LAIJEG, fdtype=fdtype)

    # columns: [length, area, I1, I2, I12, J, k1, k2, E, G]
    L = LAIJEG[:, 0]
    A = LAIJEG[:, 1]
    I = LAIJEG[:, [2, 3, 4]]
    J = LAIJEG[:, 5]
    #k1 = LAIJEG[:, 6] 
    #k2 = LAIJEG[:, 7]
    E = LAIJEG[:, 8]
    G = LAIJEG[:, 9]

    cdef = _get_bar_recovery_points(model, elem)

    strain = np.full((neids, 15), np.nan, dtype=fdtype)
    for (ieid, eid, Ai, Ii, Ji, Ei, Gi, Fei, cdefi) in zip(
        ieids, eids, A, I, J, E, G, Fe, cdef,):
        strain[ieid, :] = _recover_straini_cbar(
            Ai, Ii, Ji, Ei, Gi, Fei, cdefi, fdtype=fdtype,)

    _set_min_max_stress(strain)
    data = strain.reshape(1, *strain.shape)
    table_name = "OSTR1"
    strain_obj = RealBarStrainArray.add_static_case(
        table_name, "CBAR", eids, data, isubcase,
        is_sort1=True, is_random=False, is_msc=True,
        random_code=0, title=title, subtitle=subtitle, label=label,
    )

    strain_dict = op2.op2_results.strain
    strain_dict.cbar_strain[isubcase] = strain_obj

    strain_obj.write_f06(
        f06_file, header=None, page_stamp=page_stamp,
        page_num=page_num, is_mag_phase=False, is_sort1=True,
    )
    return neids


def _recover_straini_cbar(A: float, I: np.ndarray, J: float,
                          E: float, G: float, Fe,
                          cdef: np.ndarray,
                          fdtype: str = "float64",):
    """Get static CBAR strain at points C, D, E, F for ends A and B."""
    I1, I2, I12 = I
    stress_a, stress_b = beam_stress_at_points(
        Fe, A, I1, I2, J, cdef)

    # strain = stress / E
    strain_a = stress_a / E
    strain_b = stress_b / E

    s1a, s2a, s3a, s4a = strain_a
    s1b, s2b, s3b, s4b = strain_b

    # axial strain from element forces
    axial = Fe[0] / (A * E)

    #smaxa = max(s1a, s2a, s3a, s4a)
    #smina = min(s1a, s2a, s3a, s4a)
    #smaxb = max(s1b, s2b, s3b, s4b)
    #sminb = min(s1b, s2b, s3b, s4b)
    #MS_tension = np.nan
    #MS_compression = np.nan

    #out = (
    #    s1a, s2a, s3a, s4a, axial, smaxa, smina, MS_tension,
    #    s1b, s2b, s3b, s4b, smaxb, sminb, MS_compression,
    #)
    out = (
        s1a, s2a, s3a, s4a, axial, np.nan, np.nan, np.nan,
        s1b, s2b, s3b, s4b, np.nan, np.nan, np.nan,
    )
    return out


def _recover_stress_cbar(
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
    """Recovers static CBAR stress."""
    neids, ieids, eids = get_ieids_eids(model, element_name, eids_str)
    if not neids:
        return neids

    elem = get_element(model, element_name, ieids, eids)

    LAIJEG = elem.stiffness_info()
    forces, Fe = get_bar_force_fe(
        model, xb, dof_map, elem, ieids, neids,
        LAIJEG, fdtype=fdtype)

    # columns: [length, area, I1, I2, I12, J, E, G]
    L = LAIJEG[:, 0]
    A = LAIJEG[:, 1]
    I = LAIJEG[:, [2, 3, 4]]
    J = LAIJEG[:, 5]
    E = LAIJEG[:, 6]
    G = LAIJEG[:, 7]
    #nu = LAIJEG[:, 7]
    nu = E / (2 * G) - 1

    cdef = _get_bar_recovery_points(model, elem)

    stress = np.full((neids, 15), np.nan, dtype=fdtype)

    for (ieid, Li, Ai, Ii, Ji, Ei, Gi, nui, Fei, cdefi,
         ) in zip(ieids, L, A, I, J, E, G, nu, Fe, cdef):
        stress[ieid, :] = _recover_stressi_cbar(
            Li, Ai, Ii, Ji, Ei, Gi,
            Fei, cdefi, fdtype=fdtype,)

    _set_min_max_stress(stress)
    data = stress.reshape(1, *stress.shape)
    table_name = "OES1"
    stress_obj = RealBarStressArray.add_static_case(
        table_name, "CBAR", eids, data, isubcase,
        is_sort1=True, is_random=False, is_msc=True,
        random_code=0, title=title, subtitle=subtitle, label=label,)

    stress_dict = op2.op2_results.stress
    stress_dict.cbar_stress[isubcase] = stress_obj

    stress_obj.write_f06(
        f06_file, header=None, page_stamp=page_stamp,
        page_num=page_num, is_mag_phase=False, is_sort1=True,
    )
    return neids


def _recover_stressi_cbar(
    L: float,
    A: float,
    I: np.ndarray,
    J: float,
    E: float,
    G: float,
    Fe,
    cdef: np.ndarray,
    fdtype: str = "float64",):
    """Get static CBAR stress at points C, D, E, F for ends A and B."""
    I1, I2, I12 = I

    stress_a, stress_b = beam_stress_at_points(
        Fe, A, I1, I2, J, cdef)

    s1a, s2a, s3a, s4a = stress_a
    s1b, s2b, s3b, s4b = stress_b
    axial = Fe[0] / A

    #smaxa = max(s1a, s2a, s3a, s4a)
    #smina = min(s1a, s2a, s3a, s4a)
    #smaxb = max(s1b, s2b, s3b, s4b)
    #sminb = min(s1b, s2b, s3b, s4b)
    #MS_tension = np.nan
    #MS_compression = np.nan

    #out = (
    #    s1a, s2a, s3a, s4a, axial, smaxa, smina, MS_tension,
    #    s1b, s2b, s3b, s4b, smaxb, sminb, MS_compression,
    #)
    out = (
        s1a, s2a, s3a, s4a, axial, np.nan, np.nan, np.nan,
        s1b, s2b, s3b, s4b, np.nan, np.nan, np.nan,
    )
    return out
