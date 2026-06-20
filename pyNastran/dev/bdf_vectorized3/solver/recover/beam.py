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
from .utils import save_strain_energy, fix_xb_shape
from pyNastran.dev.bdf_vectorized3.solver.elements.beam import (
    timoshenko_stiffness,
    beam_transforms,
)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf_vectorized3.bdf import BDF, Subcase
    from pyNastran.bdf_vectorized3.cards import CBAR, PBAR, PBARL
DOF_MAP = dict[tuple[int, int], int]


def recover_beam_force(
    Ke: np.ndarray,
    Teb: np.ndarray,
    q_basic: np.ndarray,) -> np.ndarray:
    """Recover element internal forces from global displacements.

    Parameters
    ----------
    Ke : np.ndarray, shape (12, 12)
        Element stiffness in element local coordinates.
    Teb : np.ndarray, shape (12, 12)
        Transformation from element to basic coordinates.
    q_basic : np.ndarray, shape (12,)
        Nodal displacements in basic coordinates.

    Returns
    -------
    Fe : np.ndarray, shape (12,)
        Element forces in element local coordinates.
        [Fx1, Fy1, Fz1, Mx1, My1, Mz1, Fx2, Fy2, Fz2, Mx2, My2, Mz2]
    """
    #q_element = Teb @ q_basic
    #Fe = Ke @ q_element
    Fe = Ke @ Teb @ q_basic
    return Fe


def _get_beam_recovery_points(model: BDF, elem) -> np.ndarray:
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
        cdef[i_lookup, :] = prop.cdef()[i_all, :]
    return cdef


def _recover_force_cbeam(f06_file: TextIO, op2,
                         model: BDF, dof_map,
                         isubcase, xb, eids_str,
                         element_name: str, fdtype='float32',
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
    #L = LAIJEG[:, 0]
    #A = LAIJEG[:, 1]
    #I = LAIJEG[:, [2, 3, 4]]
    #J = LAIJEG[:, 5]
    #k1
    #k2
    #s1
    #s2
    #E = LAIJEG[:, 10]
    #G = LAIJEG[:, 11]

    Fe = get_beam_force_fe(
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

    data = force.reshape(nmode, neids, 8)
    table_name = 'OEF1'
    force_obj = RealCBarForceArray.add_static_case(
        table_name, 'CBAR', eids, data, isubcase,
        is_sort1=True, is_random=False, is_msc=True,
        random_code=0, title=title, subtitle=subtitle, label=label)

    force_dict = op2.op2_results.force
    force_dict.cbeam_force[isubcase] = force_obj

    force_obj.write_f06(f06_file, header=None, page_stamp=page_stamp,
                        page_num=page_num, is_mag_phase=False, is_sort1=True)
    return len(elem)


def get_beam_force_fe(model: BDF,
                      xb: np.ndarray,
                      dof_map: DOF_MAP,
                      elem: CBAR,
                      ieids, neids: int,
                      LAIJEG: np.ndarray, fdtype='float32'):

    xyz1 = model.grid.get_position_by_node_id(elem.nodes[:, 0])
    xyz2 = model.grid.get_position_by_node_id(elem.nodes[:, 1])

    v, ihat, yhat, zhat, wa, wb = elem.get_axes(xyz1, xyz2)
    Teb = beam_transforms(ihat, yhat, zhat)

    q_basic, Ke = get_beam_force_Ku(
        model, xb, dof_map,
        elem, ieids, neids,
        LAIJEG, fdtype=fdtype)

    #Fe = Ke @ Teb @ q_basic
    Fe = np.einsum('ntu,nuv,nv->nt', Ke, Teb, q_basic)
    return Fe


def get_beam_force_Ku(model: BDF,
                      xb: np.ndarray,
                      dof_map: DOF_MAP,
                      elem: CBAR,
                      ieids, neids: int,
                      LAIJEG: np.ndarray, fdtype='float32'):
    L = LAIJEG[:, 0]
    A = LAIJEG[:, 1]
    I = LAIJEG[:, [2, 3, 4]]
    J = LAIJEG[:, 5]

    k1 = LAIJEG[:, 6]
    k2 = LAIJEG[:, 7]
    #s1 = LAIJEG[:, 8]
    #s2 = LAIJEG[:, 9]

    E = LAIJEG[:, 10]
    G = LAIJEG[:, 11]

    q_basic = np.full((neids, 12), np.nan, dtype=fdtype)
    Ke = np.full((neids, 12, 12), np.nan, dtype=fdtype)

    for (ieid, nodes, Li, Ai, I1, Ji, Ei, Gi, k1i, k2i) in zip(
            ieids, elem.nodes, L, A, I, J, E, G,
            k1, k2):
        q_basici, Kei = _recover_Ku_cbeami(
            xb, dof_map, nodes,
            Li, Ai, I1, Ji, Ei, Gi,
            k1i, k2i)
        q_basic[ieid, :] = q_basici
        Ke[ieid, :] = Kei

        #(Fx1, Fy1, Fz1, Mx1, My1, Mz1,
        # Fx2, Fy2, Fz2, Mx2, My2, Mz2) = Fei
        #Fei = recover_beam_force(Ke, Teb, q_all)
        #Fei = Kei @ Tebi @ q_basici

        #axial = Fx1
        #torque = Mx1
        #shear1 = Fy1
        #shear2 = Fz1
        #bending_moment_a1 = My1
        #bending_moment_a2 = Mz1

        #bending_moment_b1 = My2
        #bending_moment_b2 = Mz2

        #force = (
        #    bending_moment_a1, bending_moment_a2,
        #    bending_moment_b1, bending_moment_b2,
        #    shear1, shear2,
        #    axial, torque)
    return q_basic, Ke


#def ke_cbeam(L, A, I, J, E, G,
#             Teb, wa, wb,
#             fdtype: str='float64'):
#    """Get the elemental stiffness matrix in the basic frame."""
#    I1, I2, I12 = I
#    k1 = k2 = 1e8
#    Ke = timoshenko_stiffness(A, E, G, L, I1, I2, J, k1, k2, pa=0, pb=0)
#    K = Teb.T @ Ke @ Teb
#    return True, K
#
#
#def ke_cbeam2(L, A, I, J, E, G, k1, k2,
#              fdtype: str='float64'):
#    """Get the elemental stiffness matrix in the basic frame."""
#    I1, I2, I12 = I
#    k1 = k2 = 1e8
#    Ke = timoshenko_stiffness(A, E, G, L, I1, I2, J, k1, k2, pa=0, pb=0)
#    K = Teb.T @ Ke @ Teb
#    return True, K


#def _recover_dxi_cbeam(xb: np.ndarray,
#                       dof_map: DOF_MAP,
#                       nodes: np.ndarray,
#                       fdtype: str='float64'):
#    """Get the static CBAR force."""
#    nid1, nid2 = nodes
#    i1 = dof_map[(nid1, 1)]
#    i2 = dof_map[(nid2, 1)]
#    q_all = np.hstack([xb[i1:i1 + 6], xb[i2:i2 + 6]])
#    return q_all


def _recover_Ku_cbeami(xb: np.ndarray,
                       dof_map: DOF_MAP,
                       nodes: np.ndarray,
                       L: np.ndarray,
                       A, I, J, E, G,
                       k1, k2,
                       fdtype: str='float64'):
    """Get the static CBAR force."""
    nid1, nid2 = nodes
    I1, I2, I12 = I

    i1 = dof_map[(nid1, 1)]
    i2 = dof_map[(nid2, 1)]

    q_basic = np.hstack([xb[i1:i1 + 6], xb[i2:i2 + 6]])

    Ke = timoshenko_stiffness(
        A, E, G, L, I1, I2, J,
        k1, k2, pa=0, pb=0)
    return q_basic, Ke


def _recover_strain_cbeam(
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
    Fe = get_beam_force_fe(
        model, xb, dof_map, elem, ieids, neids,
        LAIJEG, fdtype=fdtype)

    # columns: [length, area, I1, I2, I12, J, k1, k2, s1, s2, E, G]
    L = LAIJEG[:, 0]
    A = LAIJEG[:, 1]
    I = LAIJEG[:, [2, 3, 4]]
    J = LAIJEG[:, 5]

    #k1 = LAIJEG[:, 6]
    #k2 = LAIJEG[:, 7]
    #s1 = LAIJEG[:, 8]
    #s2 = LAIJEG[:, 9]

    E = LAIJEG[:, 10]
    G = LAIJEG[:, 11]

    cdef = _get_beam_recovery_points(model, elem)

    strain = np.full((neids, 15), np.nan, dtype=fdtype)
    for (ieid, eid, Ai, Ii, Ji, Ei, Gi, Fei, cdefi) in zip(
        ieids, eids, A, I, J, E, G, Fe, cdef,):
        strain[ieid, :] = _recover_straini_cbeam(
            Ai, Ii, Ji, Ei, Gi, Fei, cdefi, fdtype=fdtype,)

    _set_min_max_stress(strain)

    #strain = (
    #    s1a, s2a, s3a, s4a, axial, smaxa, smina, MS_tension,
    #    s1b, s2b, s3b, s4b, smaxb, sminb, MS_compression,
    #)
    data = strain.reshape(1, *strain.shape)
    table_name = "OSTR1"
    strain_obj = RealBarStrainArray.add_static_case(
        table_name, "CBAR", eids, data, isubcase,
        is_sort1=True, is_random=False, is_msc=True,
        random_code=0, title=title, subtitle=subtitle, label=label,
    )

    strain_dict = op2.op2_results.strain
    strain_dict.cbeam_strain[isubcase] = strain_obj

    strain_obj.write_f06(
        f06_file, header=None, page_stamp=page_stamp,
        page_num=page_num, is_mag_phase=False, is_sort1=True,
    )
    return neids


def _recover_straini_cbeam(A: float, I: np.ndarray, J: float,
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

    # smaxa = max(s1a, s2a, s3a, s4a)
    # smina = min(s1a, s2a, s3a, s4a)
    # smaxb = max(s1b, s2b, s3b, s4b)
    # sminb = min(s1b, s2b, s3b, s4b)
    # MS_tension = np.nan
    # MS_compression = np.nan

    #out = (
    #    s1a, s2a, s3a, s4a, axial, smaxa, smina, MS_tension,
    #    s1b, s2b, s3b, s4b, smaxb, sminb, MS_compression,
    #)
    out = (
        s1a, s2a, s3a, s4a, axial, np.nan, np.nan, np.nan,
        s1b, s2b, s3b, s4b, np.nan, np.nan, np.nan,
    )
    return out


def _recover_stress_cbeam(
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
    Fe = get_beam_force_fe(
        model, xb, dof_map, elem, ieids, neids,
        LAIJEG, fdtype=fdtype)

    # columns: [length, area, I1, I2, I12, J, k1, k2, s1, s2, E, G]
    L = LAIJEG[:, 0]
    A = LAIJEG[:, 1]
    I = LAIJEG[:, [2, 3, 4]]
    J = LAIJEG[:, 5]

    #k1 = LAIJEG[:, 6]
    #k2 = LAIJEG[:, 7]
    #s1 = LAIJEG[:, 8]
    #s2 = LAIJEG[:, 9]

    E = LAIJEG[:, 10]
    G = LAIJEG[:, 11]
    nu = E / (2 * G) - 1

    cdef = _get_beam_recovery_points(model, elem)

    stress = np.full((neids, 15), np.nan, dtype=fdtype)
    for (ieid, Li, Ai, Ii, Ji, Ei, Gi, nui,
         Fei, cdefi) in zip(
         ieids, L, A, I, J, E, G, nu, Fe, cdef):
        stress[ieid, :] = _recover_stressi_cbeam(
            Li, Ai, Ii, Ji, Ei, Gi,
            Fei, cdefi, fdtype=fdtype,)

    _set_min_max_stress(stress)
    #stress = (
    #    s1a, s2a, s3a, s4a, axial, smaxa, smina, MS_tension,
    #    s1b, s2b, s3b, s4b, smaxb, sminb, MS_compression,
    #)

    data = stress.reshape(1, *stress.shape)
    table_name = "OES1"
    stress_obj = RealBarStressArray.add_static_case(
        table_name, "CBAR", eids, data, isubcase,
        is_sort1=True, is_random=False, is_msc=True,
        random_code=0, title=title, subtitle=subtitle, label=label,)

    stress_dict = op2.op2_results.stress
    stress_dict.cbeam_stress[isubcase] = stress_obj

    stress_obj.write_f06(
        f06_file, header=None, page_stamp=page_stamp,
        page_num=page_num, is_mag_phase=False, is_sort1=True,
    )
    return neids



def _set_min_max_stress(stress):
    #out = (
    #    s1a, s2a, s3a, s4a, axial, smaxa, smina, MS_tension,
    #    s1b, s2b, s3b, s4b, smaxb, sminb, MS_compression,
    #)
    s1a, s2a, s3a, s4a = stress[:, 0], stress[:, 1], stress[:, 2], stress[:, 3]
    s1b, s2b, s3b, s4b = stress[:, 8], stress[:, 9], stress[:, 10], stress[:, 11]
    smaxa = np.max([s1a, s2a, s3a, s4a], axis=0)  # 5
    smina = np.min([s1a, s2a, s3a, s4a], axis=0)  # 6
    smaxb = np.max([s1b, s2b, s3b, s4b], axis=0)  # 12
    sminb = np.min([s1b, s2b, s3b, s4b], axis=0)  # 13
    assert smaxa.shape == s1a.shape
    stress[:, 5] = smaxa
    stress[:, 6] = smina
    stress[:, 12] = smaxb
    stress[:, 13] = sminb
    #MS_tension = np.nan
    #MS_compression = np.nan


def _recover_stressi_cbeam(
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


def beam_stress_at_points(
    Fe: np.ndarray,
    A: float,
    I1: float,
    I2: float,
    J: float,
    cdef: np.ndarray,) -> tuple[np.ndarray, np.ndarray]:
    """Compute beam longitudinal stress at C,D,E,F recovery points.

    Parameters
    ----------
    Fe : np.ndarray, shape (12,)
        Element forces in local coordinates.
    A : float
        Cross-sectional area.
    I1 : float
        Moment of inertia about axis 1 (for bending in plane 1).
    I2 : float
        Moment of inertia about axis 2 (for bending in plane 2).
    J : float
        Torsional constant.
    cdef : np.ndarray, shape (8,)
        Stress recovery point coordinates (y, z offsets).

    Returns
    -------
    stress_a : np.ndarray, shape (4,)
        Longitudinal stress at end A [C, D, E, F].
    stress_b : np.ndarray, shape (4,)
        Longitudinal stress at end B [C, D, E, F].

    Notes
    -----
    sigma = P/A + My*z/Iy - Mz*y/Iz
    The sign convention follows Nastran: point (c1,c2) means y=c1, z=c2.
    sigma = Fx/A + Mz*c1/Iz + My*c2/Iy  (Nastran convention)
    """
    c1, c2, d1, d2, e1, e2, f1, f2 = cdef

    # Safe division
    inv_I1 = 1.0 / I1 if I1 > 0.0 else 0.0
    inv_I2 = 1.0 / I2 if I2 > 0.0 else 0.0

    # force2stress maps [Fx, Fy, Fz, Mx, My, Mz] -> stress at each point
    # sigma = Fx/A + My*z/Iy + Mz*y/Iz  (Nastran: My=Fe[4], Mz=Fe[5])
    # But Nastran convention: c1=y coord, c2=z coord of recovery point
    # sigma = Fx/A + Mz*(c1)/Iz + My*(c2)/Iy
    #       = Fx/A + 0*Fy + 0*Fz + 0*Mx + c2/I1*My + c1/I2*Mz
    # where I1=Iy (bending about y, uses z-coord) and I2=Iz (bending about z, uses y-coord)
    # NX convention: sigma = P/A - M1*c1/I1 + M2*c2/I2
    # where M1=Mz (plane 1), M2=My (plane 2)
    force2stress = np.array([
        [1.0 / A, 0.0, 0.0, 0.0, c2 * inv_I1, -c1 * inv_I2],
        [1.0 / A, 0.0, 0.0, 0.0, d2 * inv_I1, -d1 * inv_I2],
        [1.0 / A, 0.0, 0.0, 0.0, e2 * inv_I1, -e1 * inv_I2],
        [1.0 / A, 0.0, 0.0, 0.0, f2 * inv_I1, -f1 * inv_I2],
    ])

    stress_a = -force2stress @ Fe[:6]
    stress_b = force2stress @ Fe[6:]
    return stress_a, stress_b

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
    # columns: [length, area, I1, I2, I12, J, k1, k2,
    #           E, G]
    L = LAIJEG[:, 0]
    A = LAIJEG[:, 1]
    I = LAIJEG[:, [2, 3, 4]]
    J = LAIJEG[:, 5]
    k1 = LAIJEG[:, 6]
    k2 = LAIJEG[:, 7]
    if element_name == 'CBAR':
        E = LAIJEG[:, 8]
        G = LAIJEG[:, 9]
        assert LAIJEG.shape[1] == 10, LAIJEG.shape
    elif element_name == 'CBEAM':
        # s1 = LAIJEG[:, 8]
        # s2 = LAIJEG[:, 9]
        E = LAIJEG[:, 10]
        G = LAIJEG[:, 11]
        assert LAIJEG.shape[1] == 12, LAIJEG.shape
    else:  # pragma: no cover
        raise NotImplementedError(element_name)
        
    assert E.min() > 0, E
    assert G.min() > 0, G

    v, ihat, yhat, zhat, wa, wb = elem.get_axes(xyz1, xyz2)

    Teb = beam_transforms(ihat, yhat, zhat)
    Ke = np.zeros((neids, 12, 12), dtype=fdtype)
    q_basic = np.zeros((neids, 12), dtype=fdtype)
    q_elem2 = np.zeros((neids, 12), dtype=fdtype)
    for (ieid, (nid1, nid2), Li, Ai, Ii, Ji, Ei, Gi, k1i, k2i, Tebi,
        ) in zip(ieids, elem.nodes, L, A, I, J, E, G, k1, k2, Teb):
        I1, I2, I12 = Ii
   
        Kei = timoshenko_stiffness(Ai, Ei, Gi, Li, I1, I2, Ji,
                                   k1i, k2i)
        Ke[ieid, :, :] = Kei

        gi1 = dof_map[(nid1, 1)]
        gi2 = dof_map[(nid2, 1)]
        q_basici = np.hstack([xb[gi1:gi1 + 6], xb[gi2:gi2 + 6]])
        q_basic[ieid, :] = q_basici

    q_elem = np.einsum('ntu,nu->nt', Teb, q_basic)
    strain_energy = 0.5 * np.einsum('nt,ntu,nu->n', q_elem, Ke, q_elem)
    #assert np.allclose(q_elem, q_elem2)
    #assert np.allclose(strain_energy, strain_energy2)

    strain_energy = strain_energy.reshape(nmode, neids, 1)
    save_strain_energy(
        op2, f06_file, page_num, page_stamp, element_name,
        strain_energy, eids, isubcase, title, subtitle, label)
    return neids
