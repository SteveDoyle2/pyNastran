from __future__ import annotations
from typing import TextIO, Optional, TYPE_CHECKING
import numpy as np

from pyNastran.dev.bdf_vectorized3.cards.base_card import searchsorted_filter
from pyNastran.dev.bdf_vectorized3.solver.utils import (
    get_ieids_eids, get_element, lambda1d)
from pyNastran.op2.op2_interface.op2_classes import (
    RealCBarForceArray,
    RealBarStrainArray,
    RealBarStressArray,
)
# from pyNastran.dev.solver.build_stiffness import ke_cbar
from .utils import fix_xb_shape
from .beam import beam_stress_at_points
#from .static_stress import _get_bar_recovery_points
from pyNastran.dev.bdf_vectorized3.solver.elements.beam import (
    timoshenko_stiffness, beam_transform, recover_beam_force,
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
        cdef[i_lookup, 0] = prop.c[i_all, 0]
        cdef[i_lookup, 1] = prop.c[i_all, 1]
        cdef[i_lookup, 2] = prop.d[i_all, 0]
        cdef[i_lookup, 3] = prop.d[i_all, 1]
        cdef[i_lookup, 4] = prop.e[i_all, 0]
        cdef[i_lookup, 5] = prop.e[i_all, 1]
        cdef[i_lookup, 6] = prop.f[i_all, 0]
        cdef[i_lookup, 7] = prop.f[i_all, 1]
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
    # columns: [length, area, I1, I2, I12, J, E, G]
    L = LAIJEG[:, 0]
    A = LAIJEG[:, 1]
    I = LAIJEG[:, [2, 3, 4]]
    J = LAIJEG[:, 5]
    E = LAIJEG[:, 6]
    G = LAIJEG[:, 7]
    assert isinstance(J, np.ndarray), (elem.type, J)
    assert isinstance(A, np.ndarray), (elem.type, A)
    assert isinstance(E, np.ndarray), (elem.type, E)

    forces, Fe = get_bar_force_fe(
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

    forces2 = np.column_stack([
        bending_moment_a1, bending_moment_a2,
        bending_moment_b1, bending_moment_b2,
        shear1, shear2, axial, torque])
    assert np.allclose(forces, forces2)

    data = forces.reshape(nmode, neids, 8)
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
    E = LAIJEG[:, 6]
    G = LAIJEG[:, 7]

    forces = np.full((neids, 8), np.nan, dtype=fdtype)
    Fe = np.full((neids, 12), np.nan, dtype=fdtype)
    v, ihat, yhat, zhat, wa, wb = elem.get_axes(xyz1, xyz2)
    for (ieid, nodes, xyz1i, xyz2i,
         Ai, I1, Ji, Ei, Gi,
         vi, ihati, yhati, zhati, wai, wbi) in zip(
            ieids, elem.nodes, xyz1, xyz2, A, I, J, E, G,
            v, ihat, yhat, zhat, wa, wb):
        Fei = _recover_forcei_cbar(
            model, xb, dof_map, nodes,
            xyz1i, xyz2i, Ai, I1, Ji, Ei, Gi,
            vi, ihati, yhati, zhati, wai, wbi)

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
        forces[ieid, :] = force
        Fe[ieid, :] = Fei
    return forces, Fe


def ke_cbar(L, A, I, J, E, G,
            ihat, jhat, khat, wa, wb,
            fdtype: str='float64'):
    """Get the elemental stiffness matrix in the basic frame."""
    I1, I2, I12 = I
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

    #(Fx1, Fy1, Fz1, Mx1, My1, Mz1,
    # Fx2, Fy2, Fz2, Mx2, My2, Mz2) = Fe
    return Fe


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

    # columns: [length, area, I1, I2, I12, J, E, G]
    Lvec = LAIJEG[:, 0]
    Avec = LAIJEG[:, 1]
    Ivec = LAIJEG[:, [2, 3, 4]]
    Jvec = LAIJEG[:, 5]
    Evec = LAIJEG[:, 6]
    Gvec = LAIJEG[:, 7]

    cdef = _get_bar_recovery_points(model, elem)

    strains = np.full((neids, 15), np.nan, dtype=fdtype)
    for (ieid, eid, Ai, Ii, Ji, Ei, Gi, Fei, cdefi) in zip(
        ieids, eids, Avec, Ivec, Jvec, Evec, Gvec, Fe, cdef,):
        strains[ieid, :] = _recover_straini_cbar(
            Ai, Ii, Ji, Ei, Gi, Fei, cdefi, fdtype=fdtype,)

    data = strains.reshape(1, *strains.shape)
    table_name = "OSTR1"
    strain_obj = RealBarStrainArray.add_static_case(
        table_name, "CBAR", eids, data, isubcase,
        is_sort1=True, is_random=False, is_msc=True,
        random_code=0, title=title, subtitle=subtitle, label=label,
    )

    strain = op2.op2_results.strain
    strain.cbar_strain[isubcase] = strain_obj

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
    C1, C2, D1, D2, E1, E2, F1, F2 = cdef
    stress_a, stress_b = beam_stress_at_points(
        Fe, A, I1, I2, J,
        C1, C2, D1, D2, E1, E2, F1, F2)

    # strain = stress / E
    strain_a = stress_a / E
    strain_b = stress_b / E

    s1a, s2a, s3a, s4a = strain_a
    s1b, s2b, s3b, s4b = strain_b

    # axial strain from element forces
    axial = Fe[0] / (A * E)

    smaxa = max(s1a, s2a, s3a, s4a)
    smina = min(s1a, s2a, s3a, s4a)
    smaxb = max(s1b, s2b, s3b, s4b)
    sminb = min(s1b, s2b, s3b, s4b)
    MS_tension = np.nan
    MS_compression = np.nan

    out = (
        s1a, s2a, s3a, s4a, axial, smaxa, smina, MS_tension,
        s1b, s2b, s3b, s4b, smaxb, sminb, MS_compression,
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
    nu = E / (2 * G) - 1

    cdef = _get_bar_recovery_points(model, elem)

    stresses = np.full((neids, 15), np.nan, dtype=fdtype)

    for (ieid, Li, Ai, Ii, Ji, Ei, Gi, nui, Fei, cdefi,
         ) in zip(ieids, L, A, I, J, E, G, nu, Fe, cdef):
        stresses[ieid, :] = _recover_stressi_cbar(
            Li, Ai, Ii, Ji, Ei, Gi,
            Fei, cdefi, fdtype=fdtype,)

    data = stresses.reshape(1, *stresses.shape)
    table_name = "OES1"
    stress_obj = RealBarStressArray.add_static_case(
        table_name, "CBAR", eids, data, isubcase,
        is_sort1=True, is_random=False, is_msc=True,
        random_code=0, title=title, subtitle=subtitle, label=label,)

    stress = op2.op2_results.stress
    stress.cbar_stress[isubcase] = stress_obj

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

    C1, C2, D1, D2, E1, E2, F1, F2 = cdef
    stress_a, stress_b = beam_stress_at_points(
        Fe, A, I1, I2, J,
        C1, C2, D1, D2, E1, E2, F1, F2)

    s1a, s2a, s3a, s4a = stress_a
    s1b, s2b, s3b, s4b = stress_b
    axial = Fe[0] / A

    smaxa = max(s1a, s2a, s3a, s4a)
    smina = min(s1a, s2a, s3a, s4a)
    smaxb = max(s1b, s2b, s3b, s4b)
    sminb = min(s1b, s2b, s3b, s4b)
    MS_tension = np.nan
    MS_compression = np.nan

    out = (
        s1a, s2a, s3a, s4a, axial, smaxa, smina, MS_tension,
        s1b, s2b, s3b, s4b, smaxb, sminb, MS_compression,
    )
    return out
