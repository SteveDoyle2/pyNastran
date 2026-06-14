from __future__ import annotations
from typing import TextIO, TYPE_CHECKING
import numpy as np

from pyNastran.dev.bdf_vectorized3.solver.utils import (
    get_ieids_eids, get_element, lambda1d)
from pyNastran.dev.bdf_vectorized3.solver.elements.beam import (
    timoshenko_stiffness,
    beam_transform,
    recover_beam_force,
    beam_stress_at_points,
)
from pyNastran.op2.op2_interface.op2_classes import (
    RealRodStrainArray,
    RealBarStrainArray,
)
from .utils import get_plot_request
from .static_force import ke_cbar
from .static_stress import _get_bar_recovery_points

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.bdf import BDF, Subcase

    DOF_MAP = dict[tuple[int, int], int]


def recover_strain_101(
    f06_file: TextIO,
    op2,
    model: BDF,
    dof_map: DOF_MAP,
    subcase: Subcase,
    xb: np.ndarray,
    fdtype: str = "float32",
    title: str = "",
    subtitle: str = "",
    label: str = "",
    page_num: int = 1,
    page_stamp: str = "PAGE %s",
):
    """Recovers element strains from STRAIN = ALL."""
    if "STRAIN" not in subcase:
        return
    eid_str = "ALL"
    unused_eids_write, write_f06, write_op2, quick_return = get_plot_request(subcase, "STRAIN")
    if quick_return:
        return page_num
    isubcase = subcase.id

    nelements = 0
    nelements += _recover_strain_rod(
        f06_file,
        op2,
        model,
        dof_map,
        isubcase,
        xb,
        eid_str,
        "CROD",
        fdtype=fdtype,
        title=title,
        subtitle=subtitle,
        label=label,
        page_num=page_num,
        page_stamp=page_stamp,
    )
    nelements += _recover_strain_rod(
        f06_file,
        op2,
        model,
        dof_map,
        isubcase,
        xb,
        eid_str,
        "CONROD",
        fdtype=fdtype,
        title=title,
        subtitle=subtitle,
        label=label,
        page_num=page_num,
        page_stamp=page_stamp,
    )
    nelements += _recover_strain_rod(
        f06_file,
        op2,
        model,
        dof_map,
        isubcase,
        xb,
        eid_str,
        "CTUBE",
        fdtype=fdtype,
        title=title,
        subtitle=subtitle,
        label=label,
        page_num=page_num,
        page_stamp=page_stamp,
    )
    nelements += _recover_strain_cbar(
        f06_file,
        op2,
        model,
        dof_map,
        isubcase,
        xb,
        eid_str,
        "CBAR",
        fdtype=fdtype,
        title=title,
        subtitle=subtitle,
        label=label,
        page_num=page_num,
        page_stamp=page_stamp,
    )
    nelements += _recover_strain_cbar(
        f06_file,
        op2,
        model,
        dof_map,
        isubcase,
        xb,
        eid_str,
        "CBEAM",
        fdtype=fdtype,
        title=title,
        subtitle=subtitle,
        label=label,
        page_num=page_num,
        page_stamp=page_stamp,
    )
    if nelements == 0:
        model.log.warning(f"no strain output...{model.card_count}; {model.bdf_filename}")


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

    strains = np.full((neids, 4), np.nan, dtype=fdtype)

    if element_name == "CONROD":
        prop = elem
    elif element_name == "CROD":
        pid = elem.property_id
        prop = model.prod.slice_card_by_property_id(pid)
    elif element_name == "CTUBE":
        pid = elem.property_id
        prop = model.ptube.slice_card_by_property_id(pid)
    else:  # pragma: no cover
        raise NotImplementedError(element_name)

    mat1 = model.mat1.slice_card_by_material_id(prop.material_id)
    G = mat1.G
    if prop.type == "PTUBE":
        J = prop.J()
    else:
        J = prop.J

    for ieid, eid, nodes, xyz1i, xyz2i, Gi, Ji in zip(ieids, eids, elem.nodes, xyz1, xyz2, G, J):
        strains[ieid, :] = _recover_straini_rod(xb, dof_map, nodes, xyz1i, xyz2i, Gi, Ji)

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
    if element_name == "CONROD":
        strain.conrod_strain[isubcase] = strain_obj
    elif element_name == "CROD":
        strain.crod_strain[isubcase] = strain_obj
    elif element_name == "CTUBE":
        strain.ctube_strain[isubcase] = strain_obj
    else:  # pragma: no cover
        raise NotImplementedError(element_name)

    strain_obj.write_f06(
        f06_file,
        header=None,
        page_stamp=page_stamp,
        page_num=page_num,
        is_mag_phase=False,
        is_sort1=True,
    )
    return neids


def _recover_straini_rod(
    xb: np.ndarray,
    dof_map: DOF_MAP,
    nodes: np.ndarray,
    xyz1: np.ndarray,
    xyz2: np.ndarray,
    G: float,
    J: float,
):
    """Get static rod strain: [axial, SMa, torsion, SMt]."""
    nid1, nid2 = nodes
    i1 = dof_map[(nid1, 1)]
    i2 = dof_map[(nid2, 1)]

    q_axial = np.array([xb[i1], xb[i1 + 1], xb[i1 + 2], xb[i2], xb[i2 + 1], xb[i2 + 2]])
    q_torsion = np.array([xb[i1 + 3], xb[i1 + 4], xb[i1 + 5], xb[i2 + 3], xb[i2 + 4], xb[i2 + 5]])
    dxyz12 = xyz1 - xyz2
    Lambda = lambda1d(dxyz12, debug=False)

    u_axial = Lambda @ q_axial
    u_torsion = Lambda @ q_torsion
    du_axial = u_axial[0] - u_axial[1]
    du_torsion = u_torsion[0] - u_torsion[1]

    L = np.linalg.norm(dxyz12)
    axial_strain = du_axial / L
    torsional_strain = du_torsion / L
    # SM fields are NaN (safety margin not computed)
    return axial_strain, np.nan, torsional_strain, np.nan


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
    page_stamp: str = "PAGE %s",
) -> int:
    """Recovers static CBAR strain."""
    neids, ieids, eids = get_ieids_eids(model, element_name, eids_str)
    if not neids:
        return neids

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

    cdef = _get_bar_recovery_points(model, elem)

    strains = np.full((neids, 15), np.nan, dtype=fdtype)

    v, ihat, yhat, zhat, wa, wb = elem.get_axes(xyz1, xyz2)
    for (ieid, eid, nodes, xyz1i, xyz2i, Ai, Ii, Ji, Ei, Gi,
         vi, ihati, yhati, zhati, wai, wbi, cdefi) in zip(
        ieids, eids, elem.nodes, xyz1, xyz2, Avec, Ivec, Jvec, Evec, Gvec,
        v, ihat, yhat, zhat, wa, wb, cdef,
    ):
        strains[ieid, :] = _recover_straini_cbar(
            model, xb, dof_map, nodes, xyz1i, xyz2i,
            Ai, Ii, Ji, Ei, Gi, vi, ihati, yhati, zhati, wai, wbi,
            cdefi, fdtype=fdtype,
        )

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


def _recover_straini_cbar(
    model: BDF,
    xb: np.ndarray,
    dof_map: DOF_MAP,
    nodes: np.ndarray,
    xyz1: np.ndarray,
    xyz2: np.ndarray,
    A: float,
    I: np.ndarray,
    J: float,
    E: float,
    G: float,
    v,
    ihat,
    jhat,
    khat,
    wa,
    wb,
    cdef: np.ndarray,
    fdtype: str = "float64",
):
    """Get static CBAR strain at points C, D, E, F for ends A and B."""
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

    C1, C2, D1, D2, E1, E2, F1, F2 = cdef
    stress_a, stress_b = beam_stress_at_points(
        Fe, A, I1, I2, J, C1, C2, D1, D2, E1, E2, F1, F2,
    )

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
