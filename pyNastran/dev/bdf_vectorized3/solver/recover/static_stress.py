r"""
Hooke's law:
   {\sigma^0} = [D][B]{u_0}

"""

from __future__ import annotations
from typing import TextIO, TYPE_CHECKING
import numpy as np

from pyNastran.dev.bdf_vectorized3.cards.base_card import searchsorted_filter
from pyNastran.dev.bdf_vectorized3.solver.utils import (
    get_ieids_eids, get_element, lambda1d)
from .rod import _recover_stress_rod
from pyNastran.dev.bdf_vectorized3.solver.elements.beam import (
    timoshenko_stiffness,
    beam_transform,
    recover_beam_force,
    beam_stress_at_points,
)
from pyNastran.op2.op2_interface.op2_classes import (
    RealBarStressArray,
    RealSpringStressArray,
    RealPlateStressArray,
    RealPlateStrainArray,
)
from .utils import get_plot_request

#from .static_force import ke_cbar
from .static_shell import recover_shell_stress_cquad4, recover_shell_stress_ctria3


if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.bdf import BDF, Subcase
    DOF_MAP = dict[tuple[int, int], int]


def recover_stress_101(
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
    """Recovers element stresses from STRESS = ALL."""
    if "STRESS" not in subcase:
        return
    eid_str = "ALL"
    unused_eids_write, write_f06, write_op2, quick_return = get_plot_request(subcase, "STRESS")
    if quick_return:
        return page_num
    isubcase = subcase.id

    nelements = 0
    nelements += _recover_stress_celas_v3(
        f06_file,
        op2,
        model,
        dof_map,
        isubcase,
        xb,
        eid_str,
        "CELAS1",
        fdtype=fdtype,
        title=title,
        subtitle=subtitle,
        label=label,
        page_num=page_num,
        page_stamp=page_stamp,
    )
    nelements += _recover_stress_celas_v3(
        f06_file,
        op2,
        model,
        dof_map,
        isubcase,
        xb,
        eid_str,
        "CELAS2",
        fdtype=fdtype,
        title=title,
        subtitle=subtitle,
        label=label,
        page_num=page_num,
        page_stamp=page_stamp,
    )
    nelements += _recover_stress_celas_v3(
        f06_file,
        op2,
        model,
        dof_map,
        isubcase,
        xb,
        eid_str,
        "CELAS3",
        fdtype=fdtype,
        title=title,
        subtitle=subtitle,
        label=label,
        page_num=page_num,
        page_stamp=page_stamp,
    )
    nelements += _recover_stress_celas_v3(
        f06_file,
        op2,
        model,
        dof_map,
        isubcase,
        xb,
        eid_str,
        "CELAS4",
        fdtype=fdtype,
        title=title,
        subtitle=subtitle,
        label=label,
        page_num=page_num,
        page_stamp=page_stamp,
    )

    nelements += _recover_stress_rod(
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
    nelements += _recover_stress_rod(
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
    nelements += _recover_stress_rod(
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
    nelements += _recover_stress_cbar(
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
    nelements += _recover_stress_cbar(
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
    nelements += _recover_stress_shell(
        f06_file,
        op2,
        model,
        dof_map,
        isubcase,
        xb,
        fdtype=fdtype,
        title=title,
        subtitle=subtitle,
        label=label,
        page_num=page_num,
        page_stamp=page_stamp,
    )

    if nelements == 0:
        model.log.warning(f"no stress output...{model.card_count}; {model.bdf_filename}")


def _recover_stress_celas_v3(
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
    """Recovers static spring stress for vectorized3 BDF."""
    neids, ieids, eids = get_ieids_eids(model, element_name, eids_str)
    if not neids:
        return neids

    elem = get_element(model, element_name, ieids, eids)

    if element_name in {"CELAS1", "CELAS3"}:
        pid = elem.property_id[ieids]
        pelas = model.pelas.slice_card_by_property_id(pid)
        k = pelas.k.ravel()
        s = pelas.s.ravel()
    elif element_name == "CELAS2":
        k = elem.k[ieids].ravel()
        s = elem.s[ieids].ravel()
    elif element_name == "CELAS4":
        k = elem.k[ieids].ravel()
        s = np.ones(neids, dtype="float64")
    else:  # pragma: no cover
        raise NotImplementedError(element_name)

    stress = np.full((neids, 1), np.nan, dtype=fdtype)
    if element_name in {"CELAS1", "CELAS2"}:
        for ieid, ki, si, nodes, comps in zip(ieids, k, s, elem.nodes, elem.components):
            nid1, nid2 = nodes
            c1, c2 = comps
            i = dof_map[(nid1, c1)]
            j = dof_map[(nid2, c2)]
            dx = xb[j] - xb[i]
            stress[ieid] = ki * si * dx
    elif element_name in {"CELAS3", "CELAS4"}:
        for ieid, ki, si, nodes in zip(ieids, k, s, elem.spoints):
            nid1, nid2 = nodes
            i = dof_map[(nid1, 0)]
            j = dof_map[(nid2, 0)]
            dx = xb[j] - xb[i]
            stress[ieid] = ki * si * dx
    else:  # pragma: no cover
        raise NotImplementedError(element_name)

    data = stress.reshape(1, *stress.shape)
    table_name = "OES1"
    stress_obj = RealSpringStressArray.add_static_case(
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

    stress_results = op2.op2_results.stress
    slot_name = f"{element_name.lower()}_stress"
    getattr(stress_results, slot_name)[isubcase] = stress_obj

    stress_obj.write_f06(
        f06_file,
        header=None,
        page_stamp=page_stamp,
        page_num=page_num,
        is_mag_phase=False,
        is_sort1=True,
    )
    return neids


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
    page_stamp: str = "PAGE %s",
) -> int:
    """Recovers static CBAR stress."""
    neids, ieids, eids = get_ieids_eids(model, element_name, eids_str)
    if not neids:
        return neids

    elem = get_element(model, element_name, ieids, eids)
    xyz1 = model.grid.get_position_by_node_id(elem.nodes[:, 0])
    xyz2 = model.grid.get_position_by_node_id(elem.nodes[:, 1])

    LAIJEG = elem.stiffness_info()
    # columns: [length, area, I1, I2, I12, J, E, G]
    A = LAIJEG[:, 0]
    A = LAIJEG[:, 1]
    I = LAIJEG[:, [2, 3, 4]]
    J = LAIJEG[:, 5]
    E = LAIJEG[:, 6]
    G = LAIJEG[:, 7]

    cdef = _get_bar_recovery_points(model, elem)

    stresses = np.full((neids, 15), np.nan, dtype=fdtype)

    v, ihat, yhat, zhat, wa, wb = elem.get_axes(xyz1, xyz2)
    for (ieid, eid, nodes, xyz1i, xyz2i, Li, Ai, Ii, Ji, Ei, Gi,
         vi, ihati, yhati, zhati, wai, wbi, cdefi) in zip(
        ieids, eids, elem.nodes, xyz1, xyz2, L, A, I, J, E, G,
        nu, ihat, yhat, zhat, wa, wb, cdef,
    ):
        stresses[ieid, :] = _recover_stressi_cbar(
            model, xb, dof_map, nodes, xyz1i, xyz2i,
            Li, Ai, Ii, Ji, Ei, Gi, vi,
            ihati, yhati, zhati, wai, wbi,
            cdefi, fdtype=fdtype,)

    data = stresses.reshape(1, *stresses.shape)
    table_name = "OES1"
    stress_obj = RealBarStressArray.add_static_case(
        table_name, "CBAR", eids, data, isubcase,
        is_sort1=True, is_random=False, is_msc=True,
        random_code=0, title=title, subtitle=subtitle, label=label,
    )

    stress = op2.op2_results.stress
    stress.cbar_stress[isubcase] = stress_obj

    stress_obj.write_f06(
        f06_file, header=None, page_stamp=page_stamp,
        page_num=page_num, is_mag_phase=False, is_sort1=True,
    )
    return neids


def _recover_stressi_cbar(
    model: BDF,
    xb: np.ndarray,
    dof_map: DOF_MAP,
    nodes: np.ndarray,
    xyz1: np.ndarray,
    xyz2: np.ndarray,
    L: float,
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
    """Get static CBAR stress at points C, D, E, F for ends A and B."""
    nid1, nid2 = nodes
    I1, I2, I12 = I

    i1 = dof_map[(nid1, 1)]
    i2 = dof_map[(nid2, 1)]

    q_all = np.hstack([xb[i1:i1 + 6], xb[i2:i2 + 6]])

    k1 = k2 = 1e8
    Ke = timoshenko_stiffness(A, E, G, L, I1, I2, J, k1, k2, pa=0, pb=0)
    Teb = beam_transform(ihat, jhat, khat)
    Fe = recover_beam_force(Ke, Teb, q_all)

    C1, C2, D1, D2, E1, E2, F1, F2 = cdef
    stress_a, stress_b = beam_stress_at_points(
        Fe, A, I1, I2, J, C1, C2, D1, D2, E1, E2, F1, F2,
    )

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


def _recover_stress_shell(
    f06_file: TextIO,
    op2,
    model: BDF,
    dof_map: DOF_MAP,
    isubcase: int,
    xb: np.ndarray,
    fdtype: str = "float32",
    title: str = "",
    subtitle: str = "",
    label: str = "",
    page_num: int = 1,
    page_stamp: str = "PAGE %s",
) -> int:
    """Recover CQUAD4/CTRIA3 shell stresses and store in op2."""
    nelements = 0
    quad_results = recover_shell_stress_cquad4(model, dof_map, xb)
    tri_results = recover_shell_stress_ctria3(model, dof_map, xb)
    nelements = len(quad_results) + len(tri_results)

    if nelements == 0:
        return 0

    # Build RealPlateStressArray for CQUAD4
    if quad_results:
        n_quad = len(quad_results)
        element_node = np.zeros((n_quad * 2, 2), dtype="int32")
        fiber = np.zeros(n_quad * 2, dtype=fdtype)
        data = np.zeros((1, n_quad * 2, 8), dtype=fdtype)
        for i, (eid, stress_data) in enumerate(sorted(quad_results.items())):
            for j in range(2):
                row = 2 * i + j
                element_node[row] = [eid, 0]
                fiber[row] = stress_data[j, 0]
                data[0, row] = stress_data[j]

        quad_stress = RealPlateStressArray.add_static_case(
            table_name="OES1",
            element_name="CQUAD4",
            nnodes=1,
            element_node=element_node,
            fiber=fiber,
            data=data,
            isubcase=isubcase,
            title=title,
            subtitle=subtitle,
            label=label,
        )
        op2.op2_results.stress.cquad4_stress[isubcase] = quad_stress

    # Build RealPlateStressArray for CTRIA3
    if tri_results:
        n_tri = len(tri_results)
        element_node = np.zeros((n_tri * 2, 2), dtype="int32")
        fiber = np.zeros(n_tri * 2, dtype=fdtype)
        data = np.zeros((1, n_tri * 2, 8), dtype=fdtype)
        for i, (eid, stress_data) in enumerate(sorted(tri_results.items())):
            for j in range(2):
                row = 2 * i + j
                element_node[row] = [eid, 0]
                fiber[row] = stress_data[j, 0]
                data[0, row] = stress_data[j]

        tri_stress = RealPlateStressArray.add_static_case(
            table_name="OES1",
            element_name="CTRIA3",
            nnodes=1,
            element_node=element_node,
            fiber=fiber,
            data=data,
            isubcase=isubcase,
            title=title,
            subtitle=subtitle,
            label=label,
        )
        op2.op2_results.stress.ctria3_stress[isubcase] = tri_stress

    return nelements


def recover_modal_stress_shell(
    op2,
    model: BDF,
    dof_map: DOF_MAP,
    isubcase: int,
    modes: np.ndarray,
    eigenvalues: np.ndarray,
    mode_shapes: np.ndarray,
    fdtype: str = "float32",
    title: str = "",
    subtitle: str = "",
    label: str = "",
) -> int:
    """Recover CQUAD4/CTRIA3 shell stresses for modal analysis (SOL 103).

    Parameters
    ----------
    modes : (nmodes,) int array — mode numbers [1, 2, ...]
    eigenvalues : (nmodes,) float — eigenvalues (omega^2)
    mode_shapes : (ndof, nmodes) — displacement per mode
    """
    nmodes = len(modes)
    cycles = np.sqrt(np.abs(eigenvalues)) / (2.0 * np.pi)

    # Compute stress for each mode
    quad_all = []
    tri_all = []
    for imode in range(nmodes):
        xb_mode = mode_shapes[:, imode]
        quad_all.append(recover_shell_stress_cquad4(model, dof_map, xb_mode))
        tri_all.append(recover_shell_stress_ctria3(model, dof_map, xb_mode))

    # CQUAD4
    if quad_all and quad_all[0]:
        eids_sorted = sorted(quad_all[0].keys())
        n_quad = len(eids_sorted)
        element_node = np.zeros((n_quad * 2, 2), dtype="int32")
        fiber = np.zeros(n_quad * 2, dtype=fdtype)
        data = np.zeros((nmodes, n_quad * 2, 8), dtype=fdtype)

        for i, eid in enumerate(eids_sorted):
            for j in range(2):
                row = 2 * i + j
                element_node[row] = [eid, 0]
                fiber[row] = quad_all[0][eid][j, 0]
            for imode in range(nmodes):
                stress_data = quad_all[imode].get(eid)
                if stress_data is not None:
                    data[imode, 2 * i] = stress_data[0]
                    data[imode, 2 * i + 1] = stress_data[1]

        quad_stress = RealPlateStressArray.add_modal_case(
            table_name="OES1",
            element_name="CQUAD4",
            nnodes=1,
            element_node=element_node,
            fiber=fiber,
            data=data,
            isubcase=isubcase,
            modes=modes,
            eigns=eigenvalues,
            cycles=cycles,
            title=title,
            subtitle=subtitle,
            label=label,
        )
        op2.op2_results.stress.cquad4_stress[isubcase] = quad_stress

    # CTRIA3
    if tri_all and tri_all[0]:
        eids_sorted = sorted(tri_all[0].keys())
        n_tri = len(eids_sorted)
        element_node = np.zeros((n_tri * 2, 2), dtype="int32")
        fiber = np.zeros(n_tri * 2, dtype=fdtype)
        data = np.zeros((nmodes, n_tri * 2, 8), dtype=fdtype)

        for i, eid in enumerate(eids_sorted):
            for j in range(2):
                row = 2 * i + j
                element_node[row] = [eid, 0]
                fiber[row] = tri_all[0][eid][j, 0]
            for imode in range(nmodes):
                stress_data = tri_all[imode].get(eid)
                if stress_data is not None:
                    data[imode, 2 * i] = stress_data[0]
                    data[imode, 2 * i + 1] = stress_data[1]

        tri_stress = RealPlateStressArray.add_modal_case(
            table_name="OES1",
            element_name="CTRIA3",
            nnodes=1,
            element_node=element_node,
            fiber=fiber,
            data=data,
            isubcase=isubcase,
            modes=modes,
            eigns=eigenvalues,
            cycles=cycles,
            title=title,
            subtitle=subtitle,
            label=label,
        )
        op2.op2_results.stress.ctria3_stress[isubcase] = tri_stress

    return len(quad_all[0]) + len(tri_all[0]) if quad_all else 0
