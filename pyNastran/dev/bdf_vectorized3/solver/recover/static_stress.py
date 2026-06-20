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
from pyNastran.op2.op2_interface.op2_classes import (
    RealBarStressArray,
    RealSpringStressArray,
    RealPlateStressArray,
    RealPlateStrainArray,
)
from .utils import get_plot_request

#from .static_force import ke_cbar
from .static_shell import recover_shell_stress_cquad4, recover_shell_stress_ctria3
from .bar import _recover_stress_cbar
from .beam import _recover_stress_cbeam
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
    page_stamp: str = "PAGE %s",):
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
        page_stamp=page_stamp,)
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
        page_stamp=page_stamp,)
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
        page_stamp=page_stamp,)
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
        page_stamp=page_stamp,)

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
        page_stamp=page_stamp,)
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
        page_stamp=page_stamp,)
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
        page_stamp=page_stamp,)
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
        page_stamp=page_stamp,)
    nelements += _recover_stress_cbeam(
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
        page_stamp=page_stamp,)
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
        page_stamp=page_stamp,)

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
    page_stamp: str = "PAGE %s",) -> int:
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
        label=label,)

    stress_results = op2.op2_results.stress
    slot_name = f"{element_name.lower()}_stress"
    getattr(stress_results, slot_name)[isubcase] = stress_obj

    stress_obj.write_f06(
        f06_file,
        header=None,
        page_stamp=page_stamp,
        page_num=page_num,
        is_mag_phase=False,
        is_sort1=True)
    return neids


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
    page_stamp: str = "PAGE %s",) -> int:
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
    label: str = "",) -> int:
    """Recover CQUAD4/CTRIA3 shell stresses for modal analysis (SOL 103).

    Parameters
    ----------
    modes : (nmodes,) int array — mode numbers [1, 2, ...]
    eigenvalues : (nmodes,) float — eigenvalues (omega^2)
    mode_shapes : (ndof, nmodes) — displacement per mode
    """
    nmodes = len(modes)
    cycles = np.sqrt(np.abs(eigenvalues)) / (2.0 * np.pi)

    case = {
        'type': 'modal',
        'eigens': eigenvalues,
        'cycles': cycles,
        'modes': modes,
    }
    del case['type']

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
            title=title,
            subtitle=subtitle,
            label=label, **case)
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
            title=title,
            subtitle=subtitle,
            label=label, **case)
        op2.op2_results.stress.ctria3_stress[isubcase] = tri_stress

    return len(quad_all[0]) + len(tri_all[0]) if quad_all else 0
