"""
This is *NOT* how the code works. I'm spitballing a better API.

Analysis is done in the material coordinate frame for composites
and in the element frame for shells/solids/bars

def theta_func_0_90(thetas: list[float]):
    thetas_out = []
    for theta in thetas:
        if np.isclose(abs(theta), 0.) or np.isclose(abs(theta), 90.):
            thetas_out.append(theta)
    return thetas_out

# all criterion must be satisfied
criteria = {
    'PSHELL': {
        1: {'eid': eids, 'pid': pids,   'result': {'vm': 20000,}},
        2: {             'mid1': mid1s, 'result': {'max': 35000., 'vm', 21000}},
    },
    'PCOMP': {
        1: {'eid': eids, 'ply': [1, 2, 3], 'result': {'e1': 2500}},
        2: {'eid': eids, 'theta': [0., 90.], 'result': {'e1': 2500}},
        3: {'eid': eids, 'theta': theta_func_0_90, 'result': {'e1': 2500}},
    },
    'PCOMPG': {
        2: {'glply: [4, 5, 6], 'result': {'e1': 2500}},
    },
}
PSHELL:
   - stress: quantity/allowable (vm, max, min)
  #- strain: quantity/allowable
PCOMP:
  - strain quantity/allowable (ply, e1, e2, e12)
  - strain quantity/allowable (core)
SOLID:
  - stress quantity/allowable (max/vm)
"""

from collections import defaultdict
import warnings
from typing import cast, Any, Optional, Callable
import numpy as np

from cpylog import SimpleLogger
from pyNastran.utils import PathLike
from pyNastran.bdf.bdf import BDF, read_bdf
from pyNastran.op2.op2 import OP2, read_op2
from pyNastran.op2.data_in_material_coord import data_in_material_coord

# from pyNastran.converters.nastran.gui.result_objects.solid_stress_results import SolidStrainStressResults2
# from pyNastran.converters.nastran.gui.result_objects.composite_stress_results import CompositeStrainStressResults2
# from pyNastran.converters.nastran.gui.result_objects.plate_stress_results import PlateStrainStressResults2

BASE_STRESS_STRAIN = [
    "max",
    "min",
    "abs_max",
    "von_mises",
    "max_shear",
]
SOLID_STRESS_STRAIN_KEYS = [
    "xx",
    "yy",
    "zz",
    "xy",
    "xz",
    "yz",
    "max",
    "min",
    "abs_max",
    "von_mises",
    "max_shear",
]
ROD_STRESS_STRAIN_KEYS = [
    "xx",
    "xy",
    "max",
    "min",
    "abs_max",
    "von_mises",
    "max_shear",
]
BAR_STRESS_STRAIN_KEYS = [
    "xx",
    "xy",
    "max",
    "min",
    "abs_max",
    "von_mises",
    "max_shear",
]
BEAM_STRESS_STRAIN_KEYS = BAR_STRESS_STRAIN_KEYS
BEND_STRESS_STRAIN_KEYS = BAR_STRESS_STRAIN_KEYS
PLATE_STRESS_STRAIN_KEYS = [
    "xx",
    "yy",
    "xy",
    "max",
    "min",
    "abs_max",
    "von_mises",
    "max_shear",
]
COMP_PLATE_STRESS_STRAIN_KEYS = BASE_STRESS_STRAIN


def envelope(
    bdf_filename: PathLike | BDF,
    op2_filename: PathLike | OP2,
    include_results: [list[str]] = None,
    exclude_results: [list[str]] = None,
    # subcases=None,
    rod_stress: str = "",
    rod_strain: str = "",
    # rod_force: str='',
    bar_stress: str = "",
    bar_strain: str = "",
    # bar_force: str='',
    beam_stress: str = "",
    beam_strain: str = "",
    # beam_force: str='',
    # cbush_stress: str='',
    # cbush_strain: str='',
    bush_force: str = "",
    plate_stress: str = "",
    plate_strain: str = "",
    comp_plate_stress: str = "",
    comp_plate_strain: str = "",
    solid_stress: str = "",
    solid_strain: str = "",
    displacement: str = "",
    spc_force: str = "",
    mpc_force: str = "",
    consider_plate_nodes: bool = True,
    consider_solid_nodes: bool = True,
    node_ids: list[int] | np.ndarray = None,
    element_ids: list[int] | np.ndarray = None,
    transform_to_global_coord: bool = False,
    transform_to_material_coord: bool = True,
    percent_nids_target: float = 1.00,
    percent_eids_target: float = 1.00,
) -> np.ndarray:
    """
    Pick one (i.e., solid_stress or solid_strain) per group.

    The advantages of this approach:
     - few inputs (it could be much worse)
    The disadvantages of this approach:
     - elements have multiple criterion (e.g., max_shear and max principal)
       -> TODO: make combined quantities max_shear_principal??? -> no
       -> cbush_force = 'fxy_rss' (this is one criterion)
     - different materials have different allowables
       -> run it multiple times

    Enveloping works on:
     - rod stress/strain
     - bar stress/strain (no von_mises or max_shear)
     - plate stress/strain
     - comp plate stress/strain
     - solid stress/strain
     - cbush_force?

    Everything else is a placeholder
     - element force
     - cbeam
     - displacements
     - spc_forces
     - mpc_forces

    Parameters
    ----------
    bdf_filename : PathLike or BDF
        the geometry model
    op2_filename : PathLike or OP2
        the results model
    node_ids : list[int] or np.ndarray; default=None
        the nodes to consider; default=all
        TODO: currently doesn't work
    element_ids : list[int] or np.ndarray; default=None
        the elements to consider; default=all
        TODO: currently doesn't work
              only disables elements if entire group is disabled
    include_results list[str] or None
        see read_op2
    exclude_results list[str] or None
        see read_op2
    bush_force : str; default=''
        group0: 'f', 'm'
        group1: 'x', 'y', 'z', 'xy', 'yz', 'xz', 'xyz',
        group2: 'max', 'min', 'rss', 'abs_max'
        cbush_force = {group0}{group1}_{group2} = 'fxy_max'
    rod_stress : str; default=''
        'max', 'min', 'abs_max', 'von_mises', 'max_shear'
    rod_strain : str; default=''
        'max', 'min', 'abs_max', 'von_mises', 'max_shear'
    bar_stress : str; default=''
        'max', 'min', 'abs_max', 'von_mises', 'max_shear'
    bar_strain : str; default=''
        'max', 'min', 'abs_max', 'von_mises', 'max_shear'
    beam_stress : str; default=''
        'max', 'min', 'abs_max', 'von_mises', 'max_shear'
    beam_strain : str; default=''
        'max', 'min', 'abs_max', 'von_mises', 'max_shear'
    plate_stress : str; default=''
        'max', 'min', 'abs_max', 'von_mises', 'max_shear'
    plate_strain : str; default=''
        'max', 'min', 'abs_max', 'von_mises', 'max_shear'
    comp_plate_stress : str; default=''
        'max', 'min', 'abs_max', 'von_mises', 'max_shear'
    comp_plate_strain : str; default=''
        'max', 'min', 'abs_max', 'von_mises', 'max_shear'
    solid_stress : str; default=''
        'max', 'min', 'abs_max', 'von_mises', 'max_shear'
    solid_strain : str; default=''
        'max', 'min', 'abs_max', 'von_mises', 'max_shear'
    displacement; default=''
        see cbush_force except f->t, m->r
        TODO: add coordinate system support?
    spc_force; default=''
        see cbush_force
        TODO: add coordinate system support?
    mpc_force; default=''
        see cbush_force
        TODO: add coordinate system support?
    consider_plate_nodes : bool=True
        if there are no plate nodes, False is automatically selected
        True:  take the max/min of all extrapolated plate nodes
        False: take the max/min of only plate centers
    consider_solid_nodes : bool=True
        if there are no solid nodes, False is automatically selected
        True:  take the max/min of all extrapolated solid nodes
        False: take the max/min of only solid centers
    transform_to_global : bool; default=False
        Transform displacement/spc_force/mpc_force to the
        global (cid=0) coordinate frame.
    transform_to_material_coord : bool; default=True
        Transform shell stress/strain/force into the
        material coordinate frame.
        doesn't transform composites or solids
    percent_eids_target : float; default=1.00 -> all
        defines the percentage of critical elements that are considered

    Returns
    -------
    subcases : np.ndarray
        The worst subcases
    """
    # rod_stress_strain_keys = [
    #     'max', 'min', 'abs_max',
    #     'von_mises', 'max_shear',
    # ]
    # is_cbush1d_min, cbush1d_force_tuple = _validate_cbush_force(cbush1d_force)
    is_cbush1d_min = False

    # is_ctube_min = is_conrod_min = is_crod_min
    displacement = _str_to_list(displacement)
    spc_force = _str_to_list(spc_force)
    mpc_force = _str_to_list(mpc_force)

    rod_stress = _str_to_list(rod_stress)
    rod_strain = _str_to_list(rod_strain)
    bar_stress = _str_to_list(bar_stress)
    bar_strain = _str_to_list(bar_strain)

    beam_stress = _str_to_list(beam_stress)
    beam_strain = _str_to_list(beam_strain)
    is_cbend_min = False
    # -----------------------------------------------------------------
    # verify requests
    plate_stress = _str_to_list(plate_stress)
    plate_strain = _str_to_list(plate_strain)

    comp_plate_stress = _str_to_list(comp_plate_stress)
    comp_plate_strain = _str_to_list(comp_plate_strain)

    solid_stress = _str_to_list(solid_stress)
    solid_strain = _str_to_list(solid_strain)

    bush_force = _str_to_list(bush_force)
    is_cbush = False
    for cbush_forcei in bush_force:
        is_cbush_min, cbush_force_tuple = _validate_force_moment(cbush_forcei)
        is_cbush = True
    # is_bush = _check_force('bush', bush_force, bush_force_keys)

    is_disp, is_spc_force, is_mpc_force = is_oug(displacement, spc_force, mpc_force)
    is_ougi = is_disp or is_spc_force or is_mpc_force

    is_rod = _check_stress_strain("rod", rod_stress, rod_strain, ROD_STRESS_STRAIN_KEYS)
    is_bar = _check_stress_strain("bar", bar_stress, bar_strain, BAR_STRESS_STRAIN_KEYS)
    is_beam = _check_stress_strain("beam", beam_stress, beam_strain, BEAM_STRESS_STRAIN_KEYS)
    is_plate = _check_stress_strain("plate", plate_stress, plate_strain, PLATE_STRESS_STRAIN_KEYS)
    is_comp_plate = _check_stress_strain(
        "comp_plate", comp_plate_stress, comp_plate_strain, COMP_PLATE_STRESS_STRAIN_KEYS
    )
    is_solid = _check_stress_strain("solid", solid_stress, solid_strain, SOLID_STRESS_STRAIN_KEYS)

    # pick at least one
    is_results = (
        is_rod
        or is_bar
        or is_beam
        or is_cbush
        or is_comp_plate
        or is_plate
        or is_solid
        or is_disp
        or is_spc_force
        or is_mpc_force
        or is_ougi
    )
    if not is_results:
        raise RuntimeError("no plate/solid results were selected")

    # -----------------------------------------------------------------
    idtype = "int32"
    fdtype = "float32"

    log = SimpleLogger(level="info")
    model, model_results, node_id, element_id = _setup_models(
        bdf_filename,
        op2_filename,
        exclude_results=exclude_results,
        include_results=include_results,
        transform_to_global_coord=transform_to_global_coord,
        transform_to_material_coord=transform_to_material_coord,
        log=log,
        idtype=idtype,
    )

    if node_ids is None:
        allow_missing_nid = False
    else:
        node_id = np.asarray(node_ids, dtype=idtype)
        allow_missing_nid = True
    node_id.sort()

    if element_ids is None:
        allow_missing_eid = False
    else:
        element_id = np.asarray(element_ids, dtype=idtype)
        allow_missing_eid = True
    element_id.sort()

    dtype = "float32"  # if model_results.size == 4 else 'float64'
    out_eids = _get_base_etypes(model)
    crod_eids, ctube_eids, conrod_eids = out_eids["rod"]
    cbar_eids = out_eids["cbar"][0]
    cbeam_eids = out_eids["cbeam"][0]
    cbend_eids = out_eids["cbend"][0]
    cbush_eids = out_eids["cbush"][0]
    cbush1d_eids = out_eids["cbush1d"][0]
    (
        ctria3_plate_eids,
        ctria6_plate_eids,
        ctriar_plate_eids,
        cquad4_plate_eids,
        cquad8_plate_eids,
        cquadr_plate_eids,
    ) = out_eids["plate"]
    (
        ctria3_comp_eids,
        ctria6_comp_eids,
        ctriar_comp_eids,
        cquad4_comp_eids,
        cquad8_comp_eids,
        cquadr_comp_eids,
    ) = out_eids["comp_plate"]
    (ctetra_eids, cpenta_eids, cpyram_eids, chexa_eids) = out_eids["solid"]

    # 'rod': (crod_eids, ctube_eids, conrod_eids),
    # 'cbar': (cbar_eids),
    # 'cbeam': (cbeam_eids),
    # 'cbend': (cbend_eids),
    # 'cbush': (cbush_eids),
    # 'cbush1d': (cbush1d_eids,),
    # 'plate': (ctria3_plate_eids, ctria6_plate_eids, ctriar_plate_eids,
    #           cquad4_plate_eids, cquad8_plate_eids, cquadr_plate_eids)

    if allow_missing_eid:
        crod_eids = np.intersect1d(crod_eids, element_id)
        conrod_eids = np.intersect1d(conrod_eids, element_id)
        ctube_eids = np.intersect1d(ctube_eids, element_id)

    if allow_missing_eid:
        cbar_eids = np.intersect1d(cbar_eids, element_id)
        cbeam_eids = np.intersect1d(cbeam_eids, element_id)
        cbend_eids = np.intersect1d(cbend_eids, element_id)
    # all_bar_eids = get_all_eids((cbar_eids, cbeam_eids, cbend_eids))

    if allow_missing_eid:
        cbush_eids = np.intersect1d(cbush_eids, element_id)
        cbush1d_eids = np.intersect1d(cbush1d_eids, element_id)
    # all_bush_eids = get_all_eids((cbush_eids, cbush1d_eids))

    if allow_missing_eid:
        cquad4_plate_eids = np.intersect1d(cquad4_plate_eids, element_id)
        ctria3_plate_eids = np.intersect1d(ctria3_plate_eids, element_id)
        cquad8_plate_eids = np.intersect1d(cquad8_plate_eids, element_id)
        ctria6_plate_eids = np.intersect1d(ctria6_plate_eids, element_id)
        cquadr_plate_eids = np.intersect1d(cquadr_plate_eids, element_id)
        ctriar_plate_eids = np.intersect1d(ctriar_plate_eids, element_id)
    # all_plate_eids = get_all_eids((ctria3_plate_eids, cquad4_plate_eids,
    #                                ctria6_plate_eids, cquad8_plate_eids,
    #                                ctriar_plate_eids, cquadr_plate_eids,))

    if allow_missing_eid:
        cquad4_comp_eids = np.intersect1d(cquad4_comp_eids, element_id)
        ctria3_comp_eids = np.intersect1d(ctria3_comp_eids, element_id)
        cquad8_comp_eids = np.intersect1d(cquad8_comp_eids, element_id)
        ctria6_comp_eids = np.intersect1d(ctria6_comp_eids, element_id)
        cquadr_comp_eids = np.intersect1d(cquadr_comp_eids, element_id)
        ctriar_comp_eids = np.intersect1d(ctriar_comp_eids, element_id)
    # all_comp_plate_eids = get_all_eids((ctria3_comp_eids, cquad4_comp_eids,
    #                                     ctria6_comp_eids, cquad8_comp_eids,
    #                                     ctriar_comp_eids, cquadr_comp_eids,))

    if allow_missing_eid:
        ctetra_eids = np.intersect1d(ctetra_eids, element_id)
        cpenta_eids = np.intersect1d(cpenta_eids, element_id)
        chexa_eids = np.intersect1d(chexa_eids, element_id)
        cpyram_eids = np.intersect1d(cpyram_eids, element_id)
    # all_solid_eids = get_all_eids((ctetra_eids, cpenta_eids, chexa_eids, cpyram_eids))
    # all_eids = get_all_eids((
    #     all_rod_eids, all_bar_eids, all_bush_eids,
    #     all_comp_plate_eids, all_plate_eids,
    #     all_solid_eids), sort=True)
    all_eids = element_id
    all_nids = node_id
    # log.info(f'all_eids = {all_eids.tolist()}')

    # print(''.join(model.get_bdf_stats()))
    # print(''.join(model_results.get_op2_stats()))
    assert len(all_eids), "no elements found"

    # data = np.zeros((nsubcase, neid), dtype='float32')
    all_eids_list = []
    all_subcases_list = []

    cbend_results = []
    cbush1d_results = []
    results_dict = {
        "cbend": (cbend_results, cbend_eids, is_cbend_min),
        "cbush1d": (cbush1d_results, cbush1d_eids, is_cbush1d_min),
    }

    node_results_dict = {}
    all_nids_list = []
    _fill_results_disp_spc_mpc_force(
        log,
        "displacements",
        node_results_dict,
        model_results.displacements,
        displacement,
        node_id,
        all_nids,
        all_nids_list,
        all_subcases_list,
        dtype=dtype,
    )

    _fill_results_disp_spc_mpc_force(
        log,
        "spc_forces",
        node_results_dict,
        model_results.spc_forces,
        spc_force,
        node_id,
        all_nids,
        all_nids_list,
        all_subcases_list,
        dtype=dtype,
    )

    _fill_results_disp_spc_mpc_force(
        log,
        "mpc_forces",
        node_results_dict,
        model_results.mpc_forces,
        mpc_force,
        node_id,
        all_nids,
        all_nids_list,
        all_subcases_list,
        dtype=dtype,
    )

    _fill_results_stress_strain(
        "rod",
        model_results,
        results_dict,
        (crod_eids, ctube_eids, conrod_eids),
        _fill_rod_list,
        rod_stress,
        rod_strain,
        all_eids,
        all_eids_list,
        all_subcases_list,
        dtype=dtype,
    )

    _fill_results_stress_strain(
        "bar",
        model_results,
        results_dict,
        (cbar_eids,),
        _fill_bar_list,
        bar_stress,
        bar_strain,
        all_eids,
        all_eids_list,
        all_subcases_list,
        dtype=dtype,
    )

    _fill_results_stress_strain(
        "beam",
        model_results,
        results_dict,
        (cbeam_eids,),
        _fill_beam_list,
        beam_stress,
        beam_strain,
        all_eids,
        all_eids_list,
        all_subcases_list,
        dtype=dtype,
    )

    is_cbush_force = len(bush_force) > 0
    if len(cbush_eids):
        if is_cbush_force:
            bush_force_list = _fill_bush_list(model_results, cbush_eids, is_force=True)
            assert len(bush_force_list), "No CBUSH force results found"
            for cbush_forcei in bush_force:
                # is_cbush_min, cbush_force_tuple = _validate_cbush_force(cbush_forcei)
                cbush_results = []
                _envelope_stress_strain(
                    "cbush force",
                    cbush_forcei,
                    bush_force_list,
                    all_eids,
                    all_eids_list,
                    all_subcases_list,
                    cbush_results,
                    dtype=dtype,
                )
                is_bush_min = "min" in cbush_forcei
                results_dict[("bush", cbush_forcei)] = (cbush_results, cbush_eids, is_bush_min)
                # log.info('_envelope_stress_strain_end')

    _fill_results_stress_strain(
        "plate",
        model_results,
        results_dict,
        (
            ctria3_plate_eids,
            cquad4_plate_eids,
            ctria6_plate_eids,
            cquad8_plate_eids,
            ctriar_plate_eids,
            cquadr_plate_eids,
        ),
        _fill_plate_list,
        plate_stress,
        plate_strain,
        all_eids,
        all_eids_list,
        all_subcases_list,
        consider_corner_nodes=consider_plate_nodes,
        dtype=dtype,
    )

    _fill_results_stress_strain(
        "comp_plate",
        model_results,
        results_dict,
        (
            ctria3_comp_eids,
            ctria6_comp_eids,
            ctriar_comp_eids,
            cquad4_comp_eids,
            cquad8_comp_eids,
            cquadr_comp_eids,
        ),
        _fill_comp_plate_list,
        comp_plate_stress,
        comp_plate_strain,
        all_eids,
        all_eids_list,
        all_subcases_list,
        dtype=dtype,
    )

    _fill_results_stress_strain(
        "solid",
        model_results,
        results_dict,
        (ctetra_eids, cpenta_eids, chexa_eids, cpyram_eids),
        _fill_solid_list,
        solid_stress,
        solid_strain,
        all_eids,
        all_eids_list,
        all_subcases_list,
        dtype=dtype,
        consider_corner_nodes=consider_solid_nodes,
    )

    all_subcases = _unique_subcases(all_subcases_list)
    out_subcases, eid_combined_data = _envelope_post(
        all_subcases, all_eids, results_dict, percent_eids_target, model, log
    )

    if len(node_results_dict):
        out_subcases_nid, nid_combined_data = _envelope_post(
            all_subcases, all_nids, node_results_dict, percent_nids_target, model, log
        )
    return out_subcases


def _fill_results_stress_strain(
    name: str,
    model_results: OP2,
    results_dict: dict,
    eids_in: tuple[np.ndarray, ...],
    fill_element_list: Callable,
    stress_request_list: list[str],
    strain_request_list: list[str],
    all_eids: np.ndarray,
    all_eids_list: list[int],
    all_subcases_list: list[int],
    dtype: str = "float64",
    consider_corner_nodes: Optional[bool] = None,
) -> None:
    is_stress = len(stress_request_list) > 0
    is_strain = len(strain_request_list) > 0
    if len(eids_in) > 1:
        eids = np.hstack(eids_in)
    else:
        eids = eids_in[0]
        assert isinstance(eids, np.ndarray)
    if not len(eids):
        return
    if is_stress:
        # log.debug(f'check_{name}_stress')
        group_name = f"{name} stress"
        stress_list = fill_element_list(model_results, *eids_in, is_stress=True)
        assert len(stress_list), f"No {name} stress results found"
        # log.info('_envelope_stress_strain')
        for stressi in stress_request_list:
            results_list = []
            if consider_corner_nodes:
                _envelope_corner_stress_strain(
                    group_name,
                    stressi,
                    stress_list,
                    all_eids,
                    all_eids_list,
                    all_subcases_list,
                    results_list,
                    consider_corner_nodes=consider_corner_nodes,
                    dtype=dtype,
                )
            else:
                _envelope_stress_strain(
                    group_name,
                    stressi,
                    stress_list,
                    all_eids,
                    all_eids_list,
                    all_subcases_list,
                    results_list,
                    dtype=dtype,
                )
            is_min = "min" in stressi
            results_dict[(name, stressi)] = (results_list, eids, is_min)
        # log.info('_envelope_stress_strain_end')

    if is_strain:
        # log.debug(f'check_{name}_strain')
        group_name = f"{name} strain"
        strain_list = fill_element_list(model_results, *eids_in, is_stress=False)
        assert len(strain_list), f"No {name} strain results found"
        for straini in strain_request_list:
            results_list = []
            if consider_corner_nodes:
                _envelope_corner_stress_strain(
                    group_name,
                    straini,
                    strain_list,
                    all_eids,
                    all_eids_list,
                    all_subcases_list,
                    results_list,
                    consider_corner_nodes=consider_corner_nodes,
                    dtype=dtype,
                )
            else:
                _envelope_stress_strain(
                    group_name,
                    straini,
                    strain_list,
                    all_eids,
                    all_eids_list,
                    all_subcases_list,
                    results_list,
                    dtype=dtype,
                )
            is_min = "min" in straini
            results_dict[(name, straini)] = (results_list, eids, is_min)


def _fill_results_disp_spc_mpc_force(
    log: SimpleLogger,
    group_name: str,
    results_dict: dict,
    obj_dict: dict,
    request_list: list[str],
    nids: np.ndarray,
    all_nids: np.ndarray,
    nids_list: list[int],
    all_subcases_list: list[int],
    dtype: str = "float64",
) -> None:
    is_request = len(request_list) > 0
    assert isinstance(nids, np.ndarray)
    is_not_nodes = len(nids) == 0
    if is_not_nodes or not is_request:
        log.debug(f"{group_name!r}: is_not_nodes={is_not_nodes} request_list={request_list}")
        return

    request_data_list = [(nids, group_name, obj_dict)]
    assert len(request_data_list), f"No {group_name} results found"
    # log.info('_fill_results_disp_spc_mpc_force')
    for requesti in request_list:
        results_list = []
        _envelope_stress_strain(
            group_name,
            requesti,
            request_data_list,
            nids,
            nids_list,
            all_subcases_list,
            results_list,
            dtype=dtype,
        )
        is_min = "min" in requesti
        print(requesti)
        results_dict[(group_name, requesti)] = (results_list, nids, is_min)
        # log.info('_envelope_stress_strain_end')
    assert len(results_dict) > 0, results_dict


def _validate_force_moment(force_moment: str) -> tuple[bool, tuple[str, str]]:
    """
    works on CBUSH and SPC/MPC force
    fx_min, fxy_abs, fxy_rss, mxyz_abs_max
    """
    is_min = False
    if force_moment == "":
        return is_min, ("", "")
    group1, group2 = force_moment.lower().split("_", 1)
    group0, group1 = group1[0], group1[1:]
    assert group0 in ["f", "m"], (group0, group1, group2)
    group1 = "".join(sorted(group1))
    assert group1 in ["x", "y", "z", "xy", "yz", "xz", "xyz"], (group1, group2)
    assert group2 in ["min", "max", "rss", "abs_max"], (group1, group2)
    is_min = group2 == "min"
    return is_min, (group0, group1, group2)


def is_oug(
    displacement: list[str], spc_force: list[str], mpc_force: list[str]
) -> tuple[bool, bool, bool]:
    is_disp = False
    is_spc_force = False
    is_mpc_force = False
    for displacementi in displacement:
        _validate_disp(displacementi)
        is_disp = True
    for spc_forcei in spc_force:
        _validate_force_moment(spc_forcei)
        is_spc_force = True
    for mpc_forcei in mpc_force:
        _validate_force_moment(mpc_forcei)
        is_mpc_force = True
    return is_disp, is_spc_force, is_mpc_force


def _validate_disp(displacement: str) -> tuple[bool, tuple[str, str]]:
    """
    tx_min, txy_abs, txy_rss, rxyz_abs_max
    """
    is_min = False
    if displacement == "":
        return is_min, ("", "")
    group1, group2 = displacement.lower().split("_", 1)
    group0, group1 = group1[0], group1[1:]
    assert group0 in ["t", "r"], (group0, group1, group2)
    group1 = "".join(sorted(group1))
    assert group1 in ["x", "y", "z", "xy", "yz", "xz", "xyz"], (group1, group2)
    assert group2 in ["min", "max", "rss", "abs_max"], (group1, group2)
    is_min = group2 == "min"
    return is_min, (group0, group1, group2)


def _str_to_list(solid_stress: list[str] | str) -> list[str]:
    if solid_stress:
        solid_stress = [solid_stress] if isinstance(solid_stress, str) else solid_stress
    else:
        solid_stress = []
    return solid_stress


def _check_stress_strain(
    name: str, rod_stress: list[str], rod_strain: list[str], rod_stress_strain_keys: list[str]
) -> bool:
    is_rod = False
    for rod_stressi in rod_stress:
        assert rod_stressi in rod_stress_strain_keys, rod_stressi
        is_rod = True
    for rod_straini in rod_strain:
        assert rod_straini in rod_stress_strain_keys, rod_straini
        is_rod = True
    return is_rod


def _check_force(
    name: str,
    # rod_stress: list[str],
    # rod_strain: list[str],
    rod_force: list[str],
    rod_force_keys: list[str],
) -> bool:
    is_rod = False
    # for rod_stressi in rod_stress:
    #     assert rod_stressi in rod_stress_strain_keys, rod_stressi
    #     is_rod = True
    # for rod_straini in rod_strain:
    #     assert rod_straini in rod_stress_strain_keys, rod_straini
    #     is_rod = True
    for rod_forcei in rod_force:
        assert rod_force in rod_force_keys, rod_forcei
        is_rod = True
    return is_rod


def _setup_models(
    bdf_filename: PathLike | BDF,
    op2_filename: PathLike | OP2,
    include_results: [list[str]] = None,
    exclude_results: [list[str]] = None,
    transform_to_global_coord: bool = False,
    transform_to_material_coord: bool = True,
    log: Optional[SimpleLogger] = None,
    idtype: str = "int32",
) -> tuple[BDF, OP2, np.ndarray, np.ndarray]:
    """
    Parameters
    ----------
    include_results list[str] or None
        see read_op2
    exclude_results list[str] or None
        see read_op2
    """
    # cards_to_skip = None
    # cards_to_skip = cards_to_skip
    if isinstance(bdf_filename, BDF):
        model = bdf_filename
    else:
        model = read_bdf(bdf_filename, log=log)
    log.level = "debug"

    node_id = np.array(list(model.nodes), dtype=idtype)
    assert node_id is not None, node_id

    element_id = np.array(list(model.elements), dtype=idtype)
    assert element_id is not None, element_id

    if isinstance(op2_filename, OP2):
        model_results = op2_filename
        # print(model_results.get_op2_stats())
    else:
        subcases = None
        model_results = read_op2(
            op2_filename,
            subcases=subcases,
            include_results=include_results,
            exclude_results=exclude_results,
            debug=None,
        )

    if transform_to_global_coord:
        raise NotImplementedError(transform_to_global_coord)
    if transform_to_material_coord:
        data_in_material_coord(model, model_results, in_place=True)
    return model, model_results, node_id, element_id


def _get_base_etypes(model: BDF) -> dict[str, np.ndarray]:
    etype_to_eids, etype_ptype_to_eids = get_elements_dict(model)
    crod_eids = np.array(etype_to_eids.get("CROD", []))
    ctube_eids = np.array(etype_to_eids.get("CTUBE", []))
    conrod_eids = np.array(etype_to_eids.get("CONROD", []))

    cbar_eids = np.array(etype_to_eids.get("CBAR", []))
    cbeam_eids = np.array(etype_to_eids.get("CBEAM", []))
    cbend_eids = np.array(etype_to_eids.get("CBEND", []))
    cbush_eids = np.array(etype_to_eids.get("CBUSH", []))
    cbush1d_eids = np.array(etype_to_eids.get("CBUSH1D", []))

    cquad4_plate_eids = np.array(etype_ptype_to_eids.get(("CQUAD4", "PSHELL"), []))
    ctria3_plate_eids = np.array(etype_ptype_to_eids.get(("CTRIA3", "PSHELL"), []))
    cquad8_plate_eids = np.array(etype_ptype_to_eids.get(("CQUAD8", "PSHELL"), []))
    ctria6_plate_eids = np.array(etype_ptype_to_eids.get(("CTRIA6", "PSHELL"), []))
    cquadr_plate_eids = np.array(etype_ptype_to_eids.get(("CQUADR", "PSHELL"), []))
    ctriar_plate_eids = np.array(etype_ptype_to_eids.get(("CTRIAR", "PSHELL"), []))

    cquad4_comp_eids = np.array(etype_ptype_to_eids.get(("CQUAD4", "PCOMP"), []))
    ctria3_comp_eids = np.array(etype_ptype_to_eids.get(("CTRIA3", "PCOMP"), []))
    cquad8_comp_eids = np.array(etype_ptype_to_eids.get(("CQUAD8", "PCOMP"), []))
    ctria6_comp_eids = np.array(etype_ptype_to_eids.get(("CTRIA6", "PCOMP"), []))
    cquadr_comp_eids = np.array(etype_ptype_to_eids.get(("CQUADR", "PCOMP"), []))
    ctriar_comp_eids = np.array(etype_ptype_to_eids.get(("CTRIAR", "PCOMP"), []))

    ctetra_eids = np.array(etype_to_eids.get("CTETRA", []))
    cpenta_eids = np.array(etype_to_eids.get("CPENTA", []))
    cpyram_eids = np.array(etype_to_eids.get("CPYRAM", []))
    chexa_eids = np.array(etype_to_eids.get("CHEXA", []))
    out = {
        "rod": (crod_eids, ctube_eids, conrod_eids),
        # 'crod': crod_eids,
        # 'ctube': ctube_eids,
        # 'conrod': conrod_eids,
        "cbar": (cbar_eids,),
        "cbeam": (cbeam_eids,),
        "cbend": (cbend_eids,),
        "cbush": (cbush_eids,),
        "cbush1d": (cbush1d_eids,),
        "plate": (
            ctria3_plate_eids,
            ctria6_plate_eids,
            ctriar_plate_eids,
            cquad4_plate_eids,
            cquad8_plate_eids,
            cquadr_plate_eids,
        ),
        "comp_plate": (
            ctria3_comp_eids,
            ctria6_comp_eids,
            ctriar_comp_eids,
            cquad4_comp_eids,
            cquad8_comp_eids,
            cquadr_comp_eids,
        ),
        "solid": (ctetra_eids, cpenta_eids, cpyram_eids, chexa_eids),
    }
    return out


def get_all_eids(eids_tuple: tuple[np.ndarray, ...], sort: bool = False) -> np.ndarray:
    eids_list = [eids for eids in eids_tuple if len(eids)]
    if len(eids_list) == 0:
        return np.array([], dtype="int32")
    all_eids = np.hstack(eids_list)
    if sort:
        all_eids.sort()
    return all_eids


def _envelope_post(
    all_subcases: np.ndarray,
    all_eids: np.ndarray,
    results: dict[str, tuple[list, np.ndarray, bool]],
    percent_eids_target: float,
    model: BDF,
    log: SimpleLogger,
    debug: bool = False,
) -> tuple[np.ndarray, Any]:
    nsubcase_all = len(all_subcases)
    neid_all = len(all_eids)
    assert nsubcase_all > 0 and neid_all > 0, (nsubcase_all, neid_all)

    comp_plate_combined_data = np.array([])
    solid_combined_data = np.array([])

    icase_critical_list = []
    icase_critical = np.full(neid_all, -1, dtype="int32")

    # result_groups = ['rod', 'cbush', 'plate', 'comp_plate', 'solid']
    # for result_group in result_groups:
    #     result_data = results[result_group]
    for result_group, result_data in results.items():
        elem_results, elem_eids, is_elem_min = result_data
        if len(elem_results):
            log.info(f"result_group={result_group!r}")
            assert "bush" not in result_group, result_group
            icase_criticali, data_critical, combined_data = _get_combined_data(
                all_subcases, result_group, elem_results, nsubcase_all, neid_all, is_elem_min
            )
            icase_critical_list.append(icase_criticali)
            icase = np.where(icase_criticali != -1)[0]
            icase_critical[icase] = icase_criticali[icase]
        else:
            log.debug(f"***NO result_group={result_group!r}?")
            print(elem_results)
        del elem_results, elem_eids, is_elem_min

    if debug:
        log.debug(f"icase_critical = {icase_critical}")
    ncritical_cases = len(icase_critical)
    assert ncritical_cases == neid_all, (ncritical_cases, neid_all)

    # potential for nan in icase_max if no result
    #   filter the -1 cases
    imissing_eids = icase_critical == -1
    if imissing_eids.sum():
        eids_missing = all_eids[imissing_eids]
        etypes_missing = np.array(list(set([model.elements[eid].type for eid in eids_missing])))
        etypes_missing.sort()
        log.warning(f"eids={eids_missing} not found; etypes_missing={etypes_missing}")

    icase_reduced = np.where(~imissing_eids)[0]
    icase_critical = icase_critical[icase_reduced]
    if debug:
        log.info(f"icase_critical (reduced) = {icase_critical}")
    assert len(icase_critical) > 0, "No results"

    # get the critical subcases
    subcase_critical = all_subcases[icase_critical]
    log.info(f"subcase_critical = {subcase_critical}")

    usubcase_critical, counts = np.unique(subcase_critical, return_counts=True)
    # print(f'usubcase_critical = {usubcase_critical}')
    # print(f'counts = {counts}')

    isort = np.argsort(counts)
    counts = counts[isort]
    usubcase_critical = usubcase_critical[isort]
    log.info(f"usubcase_critical = {usubcase_critical.tolist()}")
    log.info(f"           counts = {counts.tolist()}")

    counts_sum = counts.sum()
    counts_running_sum = np.cumsum(counts)
    counts_target = np.ceil(percent_eids_target * counts_sum)

    icount_slice = np.where(counts_running_sum > counts_target)[0]
    if len(icount_slice):
        # print(f'icount_slice = {icount_slice}')
        i = icount_slice[0]
        out_subcases = usubcase_critical[:i]
    else:
        out_subcases = usubcase_critical

    # out = model.get_xyz_in_coord_array(
    #     cid=0, fdtype='float64', idtype='int64')
    # nid_cp_cd, unused_xyz_cid, unused_xyz_cp, unused_icd_transform, unused_icp_transform = out
    # node_id = nid_cp_cd[:, 0]

    # is_stress = False
    # obj2 = PlateStrainStressResults2.load_from_code(
    #     subcase_id, model, model_results, element_id,
    #     is_stress, eid_to_nid_map,
    #     node_id=node_id,
    #     # is_variable_data_format=False,
    #     require_results=False,
    # )
    # plates = PlateStrainStressResults2(
    #     subcase_id, model,
    #     node_idy,
    #     element_id,
    #     cases, #: list[RealPlateArray],
    #     result, # str
    #     title='plate',
    #     is_fiber_distance,
    #     eid_to_nid_map,
    #     # dim_max: float=1.0,
    #     data_format='%g',
    #     is_variable_data_format=False,
    #     nlabels=None, labelsize=None, ncolors=None,
    #     colormap='', set_max_min=False,
    #     uname='PlateStressStrainResults2')
    return out_subcases, (comp_plate_combined_data, solid_combined_data)


def _unique_subcases(all_subcases_list: list) -> np.ndarray:
    """Get sorted unique subcases, handling both int and tuple keys."""
    if not all_subcases_list:
        return np.array([], dtype="int32")
    if isinstance(all_subcases_list[0], (int, np.integer)):
        return np.unique(all_subcases_list)
    unique = sorted(set(all_subcases_list))
    out = np.empty(len(unique), dtype=object)
    for i, key in enumerate(unique):
        out[i] = key
    return out


def _subcase_indices(all_subcases: np.ndarray, subcases: list) -> np.ndarray:
    """Map subcases to indices in all_subcases, handling tuple keys."""
    if all_subcases.dtype == object:
        lookup = {key: i for i, key in enumerate(all_subcases)}
        return np.array([lookup[s] for s in subcases], dtype="int32")
    return np.searchsorted(all_subcases, subcases)


def _get_combined_data(
    all_subcases: np.ndarray,
    result_group: str,
    results: list[tuple[str, np.ndarray, list[int], np.ndarray, np.ndarray]],
    nsubcase_all: int,
    neid_all: int,
    is_min: bool,
    dtype: str = "float32",
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    assert len(all_subcases) == len(set(all_subcases.tolist()))

    # print('----------------')
    # print(result_group)
    # print(f'nsubcase_all={nsubcase_all} neid_all={neid_all}')
    combined_data = np.full((nsubcase_all, neid_all), np.nan, dtype=dtype)
    if len(results) == 1:
        obj_name, subcases, eids, ieid, data = results[0]
        isubcase = _subcase_indices(all_subcases, subcases)
        combined_data[isubcase, ieid] = data
        # print('', obj_name, eids, ieid)
        # print(combined_data)
    else:
        assert len(results)
        for obj_name, subcases, eids, ieid, data in results:
            # print('', obj_name, eids, ieid)
            # print(data)
            isubcase = _subcase_indices(all_subcases, subcases)
            combined_data[isubcase, ieid] = data

    # find the driving case (and data) for each element
    icase_critical = np.full(neid_all, -1, dtype="int32")
    is_not_nan = ~np.all(np.isnan(combined_data), axis=0)
    # print('is_not_nan', is_not_nan, neid_all)

    # print('combined_data\n', combined_data)
    combined_datai = combined_data[:, is_not_nan]
    # print('combined_datai\n', combined_datai)
    # print('shape', combined_datai.shape)
    # print('is_not_nan', is_not_nan)
    # print('icase_critical', icase_critical)

    func = np.nanargmin if is_min else np.nanargmax
    icase_critical[is_not_nan] = func(combined_datai, axis=0)
    # print(f'{result_group}_icase_critical = {icase_critical}')

    data_critical = combined_data[icase_critical]
    del result_group
    return icase_critical, data_critical, combined_data


def get_elements_dict(model: BDF) -> tuple[dict[str, list[int]], dict[str, list[int]]]:
    # log = model.log
    eid_to_nid_map = {}
    etype_to_eids = defaultdict(list)
    etype_ptype_to_eids = defaultdict(list)
    for eid, elem in model.elements.items():
        eid_to_nid_map[eid] = elem.nodes
        prop = elem.pid_ref
        etype_to_eids[elem.type].append(eid)
        prop_type = prop.type
        if prop_type in {"PBARL", "PBEAML", "PCOMPG"}:
            prop_type = prop_type[:-1]
        etype_ptype_to_eids[(elem.type, prop_type)].append(eid)
    return dict(etype_to_eids), dict(etype_ptype_to_eids)


def _fill_rod_list(
    model_results: OP2,
    crod_eids: np.ndarray,
    ctube_eids: np.ndarray,
    conrod_eids: np.ndarray,
    is_stress: bool,
) -> list:
    if is_stress:
        word = "stress"
        stress = model_results.op2_results.stress
        crod_obj = stress.crod_stress
        ctube_obj = stress.ctube_stress
        conrod_obj = stress.conrod_stress
    else:
        word = "strain"
        strain = model_results.op2_results.strain
        crod_obj = strain.crod_strain
        ctube_obj = strain.ctube_strain
        conrod_obj = strain.conrod_strain

    rod_list = []
    # print('eids-rods', crod_eids, ctube_eids, conrod_eids)
    if len(crod_eids):
        rod_list.append((crod_eids, f"crod_{word}", crod_obj))
    if len(ctube_eids):
        rod_list.append((ctube_eids, f"ctube_{word}", ctube_obj))
    if len(conrod_eids):
        rod_list.append((conrod_eids, f"conrod_{word}", conrod_obj))
    return rod_list


def _fill_bush_list(
    model_results: OP2,
    cbush_eids: np.ndarray,
    is_stress: bool = False,
    is_strain: bool = False,
    is_force: bool = False,
) -> list:
    assert any([is_stress, is_strain, is_force]), (is_stress, is_strain, is_force)
    assert sum([is_stress, is_strain, is_force]) == 1, (is_stress, is_strain, is_force)
    if is_stress:
        word = "stress"
        stress = model_results.op2_results.stress
        cbush_obj = stress.cbush_stress
    if is_strain:
        word = "strain"
        strain = model_results.op2_results.strain
        cbush_obj = strain.cbush_strain
    if is_force:
        word = "force"
        force = model_results.op2_results.force
        cbush_obj = force.cbush_force

    bush_list = []
    if len(cbush_eids):
        bush_list.append((cbush_eids, f"cbush_{word}", cbush_obj))
    return bush_list


def _fill_bar_list(model_results: OP2, cbar_eids: np.ndarray, is_stress: bool) -> list:
    if is_stress:
        word = "stress"
        stress = model_results.op2_results.stress
        cbar_obj = stress.cbar_stress
    else:
        word = "strain"
        strain = model_results.op2_results.strain
        cbar_obj = strain.cbar_strain

    bar_list = []
    if len(cbar_eids):
        bar_list.append((cbar_eids, f"cbar_{word}", cbar_obj))
    return bar_list


def _fill_beam_list(model_results: OP2, cbeam_eids: np.ndarray, is_stress: bool) -> list:
    if is_stress:
        word = "stress"
        stress = model_results.op2_results.stress
        cbeam_obj = stress.cbeam_stress
    else:
        word = "strain"
        strain = model_results.op2_results.strain
        cbeam_obj = strain.cbeam_strain

    beam_list = []
    if len(cbeam_eids):
        beam_list.append((cbeam_eids, f"cbeam_{word}", cbeam_obj))
    return beam_list


def _fill_plate_list(
    model_results: OP2,
    ctria3_eids: np.ndarray,
    cquad4_eids: np.ndarray,
    ctria6_eids: np.ndarray,
    cquad8_eids: np.ndarray,
    ctriar_eids: np.ndarray,
    cquadr_eids: np.ndarray,
    is_stress: bool,
) -> list:
    if is_stress:
        word = "stress"
        stress = model_results.op2_results.stress
        ctria3_obj = stress.ctria3_stress
        cquad4_obj = stress.cquad4_stress
        ctria6_obj = stress.ctria6_stress
        cquad8_obj = stress.cquad8_stress
        ctriar_obj = stress.ctriar_stress
        cquadr_obj = stress.cquadr_stress
    else:
        word = "strain"
        strain = model_results.op2_results.strain
        ctria3_obj = strain.ctria3_strain
        cquad4_obj = strain.cquad4_strain
        ctria6_obj = strain.ctria6_strain
        cquad8_obj = strain.cquad8_strain
        ctriar_obj = strain.ctriar_strain
        cquadr_obj = strain.cquadr_strain

    plate_list = []
    if len(ctria3_eids):
        plate_list.append((ctria3_eids, f"ctria3_{word}", ctria3_obj))
    if len(cquad4_eids):
        plate_list.append((cquad4_eids, f"cquad4_{word}", cquad4_obj))
    if len(ctria6_eids):
        plate_list.append((ctria6_eids, f"ctria6_{word}", ctria6_obj))
    if len(cquad8_eids):
        plate_list.append((cquad8_eids, f"cquad8_{word}", cquad8_obj))
    if len(ctriar_eids):
        plate_list.append((ctriar_eids, f"ctriar_{word}", ctriar_obj))
    if len(cquadr_eids):
        plate_list.append((cquadr_eids, f"cquadr_{word}", cquadr_obj))
    return plate_list


def _fill_comp_plate_list(
    model_results: OP2,
    ctria3_eids: np.ndarray,
    cquad4_eids: np.ndarray,
    ctria6_eids: np.ndarray,
    cquad8_eids: np.ndarray,
    ctriar_eids: np.ndarray,
    cquadr_eids: np.ndarray,
    is_stress: bool,
) -> list:
    if is_stress:
        word = "stress"
        stress = model_results.op2_results.stress
        ctria3_obj = stress.ctria3_composite_stress
        cquad4_obj = stress.cquad4_composite_stress
        ctria6_obj = stress.ctria6_composite_stress
        cquad8_obj = stress.cquad8_composite_stress
        ctriar_obj = stress.ctriar_composite_stress
        cquadr_obj = stress.cquadr_composite_stress
    else:
        word = "strain"
        strain = model_results.op2_results.strain
        ctria3_obj = strain.ctria3_composite_strain
        cquad4_obj = strain.cquad4_composite_strain
        ctria6_obj = strain.ctria6_composite_strain
        cquad8_obj = strain.cquad8_composite_strain
        ctriar_obj = strain.ctriar_composite_strain
        cquadr_obj = strain.cquadr_composite_strain

    plate_list = []
    if len(ctria3_eids):
        plate_list.append((ctria3_eids, f"ctria3_{word}", ctria3_obj))
    if len(cquad4_eids):
        plate_list.append((cquad4_eids, f"cquad4_{word}", cquad4_obj))
    if len(ctria6_eids):
        plate_list.append((ctria6_eids, f"ctria6_{word}", ctria6_obj))
    if len(cquad8_eids):
        plate_list.append((cquad8_eids, f"cquad8_{word}", cquad8_obj))
    if len(ctriar_eids):
        plate_list.append((ctriar_eids, f"ctriar_{word}", ctriar_obj))
    if len(cquadr_eids):
        plate_list.append((cquadr_eids, f"cquadr_{word}", cquadr_obj))
    return plate_list


def _fill_solid_list(
    model_results: OP2,
    ctetra_eids: np.ndarray,
    cpenta_eids: np.ndarray,
    chexa_eids: np.ndarray,
    cpyram_eids: np.ndarray,
    is_stress: bool,
) -> list:
    if is_stress:
        word = "stress"
        stress = model_results.op2_results.stress
        ctetra_obj = stress.ctetra_stress
        cpenta_obj = stress.cpenta_stress
        chexa_obj = stress.chexa_stress
        cpyram_obj = stress.cpyram_stress
    else:
        word = "strain"
        strain = model_results.op2_results.strain
        ctetra_obj = strain.ctetra_strain
        cpenta_obj = strain.cpenta_strain
        chexa_obj = strain.chexa_strain
        cpyram_obj = strain.cpyram_strain

    solid_list = []
    if len(ctetra_eids):
        solid_list.append((ctetra_eids, f"ctetra_{word}", ctetra_obj))
    if len(chexa_eids):
        solid_list.append((chexa_eids, f"chexa_{word}", chexa_obj))
    if len(cpenta_eids):
        solid_list.append((cpenta_eids, f"cpenta_{word}", cpenta_obj))
    if len(cpyram_eids):
        solid_list.append((cpyram_eids, f"cpyram_{word}", cpyram_obj))
    return solid_list


def _envelope_corner_stress_strain(
    result_group: str,
    result_name: str,
    stress_list: list[tuple[np.ndarray, str, dict]],
    all_eids: np.ndarray,
    all_eids_list: list[list[int]],
    all_subcases_list: list[int],
    #  results.append((obj_name, subcases, eids, ieid, comp_data))
    results_list: list[tuple[str, list[int], list[int], np.ndarray, np.ndarray]],
    consider_corner_nodes: bool = True,
    dtype: str = "float32",
):
    for eids, obj_name, obj_dict in stress_list:
        subcases = list(obj_dict)
        nsubcase_obj = len(subcases)
        if nsubcase_obj == 0:
            warnings.warn(f"no results for {obj_name}")
            break

        neid_obj = len(eids)
        solid_data = np.full((nsubcase_obj, neid_obj), np.nan, dtype=dtype)
        for isubcase, subcase in enumerate(subcases):
            # print(f'{obj_name}: isubcase={isubcase} subcase={subcase}')
            obj = obj_dict[subcase]
            if isubcase == 0:
                ieid = np.searchsorted(all_eids, eids)
                # print(f'  solid_eids = {eids}')
                # print(f'  ieid = {ieid}')
                missing_eids = np.setdiff1d(eids, all_eids)
                if len(missing_eids):
                    raise RuntimeError(f"missing solid eids={missing_eids.tolist()}")
                assert np.array_equal(all_eids[ieid], eids)
            datai = obj.envelope(eids, result_name, consider_corner_nodes)
            solid_data[isubcase, :] = datai
        # print(f'subcases = {subcases}')
        obj_name = cast(str, obj_name)
        results_list.append((obj_name, subcases, eids, ieid, solid_data))
        all_eids_list.append(eids)
        all_subcases_list.extend(subcases)


def _envelope_stress_strain(
    group_name: str,
    result_name: str,
    comp_plate_list: list[tuple[np.ndarray, str, dict]],
    all_eids: np.ndarray,
    all_eids_list: list[list[int]],
    all_subcases_list: list[int],
    #  results.append((obj_name, subcases, eids, ieid, comp_data))
    results: list[tuple[str, np.ndarray, list[int], np.ndarray, np.ndarray]],
    dtype: str = "float32",
):
    """no consideration for corner nodes"""
    assert isinstance(result_name, str), result_name
    for eids, obj_name, obj_dict in comp_plate_list:
        subcases = list(obj_dict)
        nsubcase_obj = len(subcases)
        if nsubcase_obj == 0:
            warnings.warn(f"no results for {obj_name}")
            break

        neid_obj = len(eids)
        comp_data = np.full((nsubcase_obj, neid_obj), np.nan, dtype=dtype)
        for isubcase, subcase in enumerate(subcases):
            # print(f'{obj_name}: isubcase={isubcase} subcase={subcase}')
            obj = obj_dict[subcase]
            if isubcase == 0:
                ieid = np.searchsorted(all_eids, eids)
                # print(f'  eids = {eids}')
                # print(f'  ieid = {ieid}')
                missing_eids = np.setdiff1d(eids, all_eids)
                if len(missing_eids):
                    raise RuntimeError(f"missing {group_name!r} eids={missing_eids.tolist()}")
                assert np.array_equal(all_eids[ieid], eids)

            print(group_name)
            if group_name == "displacements":
                is_min, (group0, group1, group2) = _validate_disp(result_name)
                result_name = (group0, group1, group2)
            datai = obj.envelope(eids, result_name)
            comp_data[isubcase, :] = datai
        # print(f'subcases = {subcases}')

        # subcases, eids, ieid, combined_data
        results.append((obj_name, subcases, eids, ieid, comp_data))
        all_eids_list.append(eids)
        all_subcases_list.extend(subcases)
