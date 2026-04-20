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
from typing import Any
import numpy as np

from cpylog import SimpleLogger
from pyNastran.utils import PathLike
from pyNastran.bdf.bdf import BDF, read_bdf
from pyNastran.op2.op2 import OP2, read_op2
# from pyNastran.converters.nastran.gui.result_objects.solid_stress_results import SolidStrainStressResults2
# from pyNastran.converters.nastran.gui.result_objects.composite_stress_results import CompositeStrainStressResults2
# from pyNastran.converters.nastran.gui.result_objects.plate_stress_results import PlateStrainStressResults2

def envelope(
        bdf_filename: PathLike | BDF,
        op2_filename: OP2,
        include_results: [list[str]]=None,
        exclude_results: [list[str]]=None,
        # eids=None,
        # nids=None,
        # subcases=None,
        rod_stress: str='',
        rod_strain: str='',
        # rod_force: str='',
        bar_stress: str='',
        bar_strain: str='',
        # bar_force: str='',
        beam_stress: str='',
        beam_strain: str='',
        # beam_force: str='',
        # cbush_stress: str='',
        # cbush_strain: str='',
        # cbush_force: str='',
        plate_stress: str='',
        plate_strain: str='',
        comp_plate_stress: str='',
        comp_plate_strain: str='',
        solid_stress: str='',
        solid_strain: str='',
        consider_plate_nodes: bool=True,
        consider_solid_nodes: bool=True,
        element_ids: list[int] | np.ndarray=None,
        percent_eids_target: float=1.00,
    ) -> np.ndarray:
    """
    Pick one (i.e., solid_stress or solid_strain) per group.

    The advantages of this approach:
     - few inputs (it could be much worse)
    The disadvantages of this approach:
     - elements have multiple criterion (e.g., max_shear and max principal)
       -> TODO: make combined quantities max_shear_principal???
     - different materials have different allowables
       -> run it multiple times

    Enveloping works on:
     - rod stress/strain
     - bar stress/strain (no von_mises or max_shear)
     - plate stres/strain
     - comp plate stres/strain
     - solid stress/strain

    Everything else is a placeholder
     - element force
     - cbeam/cbush
     - displacements
     - spc_forces
     - mpc_forces

    Parameters
    ----------
    bdf_filename : PathLike or BDF
        the geometry model
    op2_filename : PathLike or OP2
        the results model
    element_ids : list[int] or np.ndarray; default=None
        the elements to consider; default=all
        TODO: currently doesn't work
              only disables elements if entire group is disabled
    include_results list[str] or None
        see read_op2
    exclude_results list[str] or None
        see read_op2
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
        CTRIA3, CQUAD4 only
        'max', 'min', 'abs_max', 'von_mises', 'max_shear'
    plate_strain : str; default=''
        CTRIA3, CQUAD4 only
        'max', 'min', 'abs_max', 'von_mises', 'max_shear'
    comp_plate_stress : str; default=''
        CTRIA3, CQUAD4 only
        'max', 'min', 'abs_max', 'von_mises', 'max_shear'
    comp_plate_strain : str; default=''
        CTRIA3, CQUAD4 only
        'max', 'min', 'abs_max', 'von_mises', 'max_shear'
    solid_stress : str; default=''
        'max', 'min', 'abs_max', 'von_mises', 'max_shear'
    solid_strain : str; default=''
        'max', 'min', 'abs_max', 'von_mises', 'max_shear'
    consider_plate_nodes : bool=True
        if there are no plate nodes, False is automatically selected
        True:  take the max/min of all extrapolated plate nodes
        False: take the max/min of only plate centers
    consider_solid_nodes : bool=True
        if there are no solid nodes, False is automatically selected
        True:  take the max/min of all extrapolated solid nodes
        False: take the max/min of only solid centers

    Returns
    -------
    subcases : np.ndarray
        The worst subcases
    """
    solid_stress_strain_keys = [
        'max', 'min', 'abs_max',
        'von_mises', 'max_shear',
    ]
    is_cbush_min = is_cbush1d_min = False

    assert rod_stress == '' or rod_stress in solid_stress_strain_keys, rod_stress
    assert rod_strain == '' or rod_strain in solid_stress_strain_keys, rod_strain
    is_crod_min = ((rod_stress == 'min') or (rod_strain == 'min'))
    # is_ctube_min = is_conrod_min = is_crod_min

    assert bar_stress == '' or bar_stress in solid_stress_strain_keys, bar_stress
    assert bar_strain == '' or bar_strain in solid_stress_strain_keys, bar_strain
    is_cbar_min = ((bar_stress == 'min') or (bar_strain == 'min'))

    assert beam_stress == '' or beam_stress in solid_stress_strain_keys, beam_stress
    assert beam_strain == '' or beam_strain in solid_stress_strain_keys, beam_strain
    is_cbeam_min = ((beam_stress == 'min') or (beam_strain == 'min'))
    #-----------------------------------------------------------------
    # verify requests
    assert plate_stress == '' or plate_stress in solid_stress_strain_keys, plate_stress
    assert plate_strain == '' or plate_strain in solid_stress_strain_keys, plate_strain
    is_plate_min = ((plate_stress == 'min') or (plate_strain == 'min'))

    assert comp_plate_stress == '' or comp_plate_stress in solid_stress_strain_keys, comp_plate_stress
    assert comp_plate_strain == '' or comp_plate_strain in solid_stress_strain_keys, comp_plate_strain
    is_comp_plate_min = ((comp_plate_stress == 'min') or (comp_plate_strain == 'min'))

    assert solid_stress == '' or solid_stress in solid_stress_strain_keys, solid_stress
    assert solid_strain == '' or solid_strain in solid_stress_strain_keys, solid_strain
    is_solid_min = ((solid_stress == 'min') or (solid_strain == 'min'))
    check_solid_stress = solid_stress != ''
    check_solid_strain = solid_strain != ''
    check_plate_stress = plate_stress != ''
    check_plate_strain = plate_strain != ''
    check_comp_plate_stress = comp_plate_stress != ''
    check_comp_plate_strain = comp_plate_strain != ''

    check_rod_stress = rod_stress != ''
    check_rod_strain = rod_strain != ''

    check_bar_stress = bar_stress != ''
    check_bar_strain = bar_strain != ''

    check_beam_stress = beam_stress != ''
    check_beam_strain = beam_strain != ''

    # pick one
    is_results = (
            # check_rod_stress or check_rod_strain or
            # check_bar_stress or check_bar_strain or
            check_comp_plate_stress or check_comp_plate_strain or
            check_plate_stress or check_plate_strain or
            check_solid_stress or check_solid_strain)
    if not is_results:
        raise RuntimeError('no plate/solid results were selected')

    if check_rod_stress and check_rod_strain:
        raise RuntimeError('Cannot downselect rod stress and strain')
    if check_bar_stress and check_bar_strain:
        raise RuntimeError('Cannot downselect bar stress and strain')
    # if check_bush_stress and check_bush_strain:
    #     raise RuntimeError('Cannot downselect bush stress and strain')

    if check_comp_plate_stress and check_comp_plate_strain:
        raise RuntimeError('Cannot downselect comp plate stress and strain')
    if check_plate_stress and check_plate_strain:
        raise RuntimeError('Cannot downselect plate stress and strain')
    if check_solid_stress and check_solid_strain:
        raise RuntimeError('Cannot downselect solid stress and strain')
    #-----------------------------------------------------------------

    # cards_to_skip = None
    # cards_to_skip = cards_to_skip
    log = SimpleLogger(level='info')
    if isinstance(bdf_filename, BDF):
        model = bdf_filename
    else:
        model = read_bdf(bdf_filename, log=log)
    log.level = 'debug'

    idtype = 'int32'
    if element_ids is None:
        element_id = np.array(list(model.elements), dtype=idtype)
        allow_missing = False
    else:
        element_id = np.asarray(element_ids, dtype=idtype)
        allow_missing = True
    element_id.sort()

    if isinstance(op2_filename, OP2):
        model_results = op2_filename
        # print(model_results.get_op2_stats())
    else:
        subcases = None
        model_results = read_op2(
            op2_filename, subcases=subcases,
            include_results=include_results,
            exclude_results=exclude_results,
            debug=None)
    dtype = 'float32'
    etype_to_eids, etype_ptype_to_eids = get_elements_dict(model)

    crod_eids = np.array(etype_to_eids.get('CROD', []))
    ctube_eids = np.array(etype_to_eids.get('CTUBE', []))
    conrod_eids = np.array(etype_to_eids.get('CONROD', []))
    if allow_missing:
        crod_eids = np.intersect1d(crod_eids, element_id)
        conrod_eids = np.intersect1d(conrod_eids, element_id)
        ctube_eids = np.intersect1d(ctube_eids, element_id)
    all_rod_eids = get_all_eids((crod_eids, ctube_eids, conrod_eids))

    cbar_eids = np.array(etype_to_eids.get('CBAR', []))
    cbeam_eids = np.array(etype_to_eids.get('CBEAM', []))
    cbend_eids = np.array(etype_to_eids.get('CBEND', []))
    if allow_missing:
        cbar_eids = np.intersect1d(cbar_eids, element_id)
        cbeam_eids = np.intersect1d(cbeam_eids, element_id)
        cbend_eids = np.intersect1d(cbend_eids, element_id)
    all_bar_eids = get_all_eids((cbar_eids, cbeam_eids, cbend_eids))

    cbush_eids = np.array(etype_to_eids.get('CBUSH', []))
    cbush1d_eids = np.array(etype_to_eids.get('CBUSH1D', []))
    if allow_missing:
        cbush_eids = np.intersect1d(cbush_eids, element_id)
        cbush1d_eids = np.intersect1d(cbush1d_eids, element_id)
    all_bush_eids = get_all_eids((cbush_eids, cbush1d_eids))

    cquad4_plate_eids = np.array(etype_ptype_to_eids.get(('CQUAD4', 'PSHELL'), []))
    ctria3_plate_eids = np.array(etype_ptype_to_eids.get(('CTRIA3', 'PSHELL'), []))
    cquad8_plate_eids = np.array(etype_ptype_to_eids.get(('CQUAD8', 'PSHELL'), []))
    ctria6_plate_eids = np.array(etype_ptype_to_eids.get(('CTRIA6', 'PSHELL'), []))
    cquadr_plate_eids = np.array(etype_ptype_to_eids.get(('CQUADR', 'PSHELL'), []))
    ctriar_plate_eids = np.array(etype_ptype_to_eids.get(('CTRIAR', 'PSHELL'), []))
    if allow_missing:
        cquad4_plate_eids = np.intersect1d(cquad4_plate_eids, element_id)
        ctria3_plate_eids = np.intersect1d(ctria3_plate_eids, element_id)
        cquad8_plate_eids = np.intersect1d(cquad8_plate_eids, element_id)
        ctria6_plate_eids = np.intersect1d(ctria6_plate_eids, element_id)
        cquadr_plate_eids = np.intersect1d(cquadr_plate_eids, element_id)
        ctriar_plate_eids = np.intersect1d(ctriar_plate_eids, element_id)
    all_plate_eids = get_all_eids((ctria3_plate_eids, cquad4_plate_eids,
                                   ctria6_plate_eids, cquad8_plate_eids,
                                   ctriar_plate_eids, cquadr_plate_eids,))

    cquad4_comp_eids = np.array(etype_ptype_to_eids.get(('CQUAD4', 'PCOMP'), []))
    ctria3_comp_eids = np.array(etype_ptype_to_eids.get(('CTRIA3', 'PCOMP'), []))
    cquad8_comp_eids = np.array(etype_ptype_to_eids.get(('CQUAD8', 'PCOMP'), []))
    ctria6_comp_eids = np.array(etype_ptype_to_eids.get(('CTRIA6', 'PCOMP'), []))
    cquadr_comp_eids = np.array(etype_ptype_to_eids.get(('CQUADR', 'PCOMP'), []))
    ctriar_comp_eids = np.array(etype_ptype_to_eids.get(('CTRIAR', 'PCOMP'), []))
    if allow_missing:
        cquad4_comp_eids = np.intersect1d(cquad4_comp_eids, element_id)
        ctria3_comp_eids = np.intersect1d(ctria3_comp_eids, element_id)
        cquad8_comp_eids = np.intersect1d(cquad8_comp_eids, element_id)
        ctria6_comp_eids = np.intersect1d(ctria6_comp_eids, element_id)
        cquadr_comp_eids = np.intersect1d(cquadr_comp_eids, element_id)
        ctriar_comp_eids = np.intersect1d(ctriar_comp_eids, element_id)
    all_comp_plate_eids = get_all_eids((ctria3_comp_eids, cquad4_comp_eids,
                                        ctria6_comp_eids, cquad8_comp_eids,
                                        ctriar_comp_eids, cquadr_comp_eids,))

    ctetra_eids = np.array(etype_to_eids.get('CTETRA', []))
    cpenta_eids = np.array(etype_to_eids.get('CPENTA', []))
    chexa_eids = np.array(etype_to_eids.get('CHEXA', []))
    cpyram_eids = np.array(etype_to_eids.get('CPYRAM', []))
    if allow_missing:
        ctetra_eids = np.intersect1d(ctetra_eids, element_id)
        cpenta_eids = np.intersect1d(cpenta_eids, element_id)
        chexa_eids = np.intersect1d(chexa_eids, element_id)
        cpyram_eids = np.intersect1d(cpyram_eids, element_id)
    all_solid_eids = get_all_eids((ctetra_eids, cpenta_eids, chexa_eids, cpyram_eids))
    # all_eids = get_all_eids((
    #     all_rod_eids, all_bar_eids, all_bush_eids,
    #     all_comp_plate_eids, all_plate_eids,
    #     all_solid_eids), sort=True)
    all_eids = element_id
    # log.info(f'all_eids = {all_eids.tolist()}')

    # print(''.join(model.get_bdf_stats()))
    # print(''.join(model_results.get_op2_stats()))
    assert len(all_eids), 'no elements found'

    # neid_solid = sum([len(mylist) for mylist in solid_eids_tuple])
    # neid_plate = sum([len(mylist) for mylist in
    #                  (ctria3_eids, cquad4_eids)])
    # neid_rod = sum([len(mylist) for mylist in
    #                (crod_eids, ctube_eids, conrod_eids)])
    # data = np.zeros((nsubcase, neid), dtype='float32')

    all_eids_list = []
    all_subcases_list = []

    rod_results = []
    cbar_results = []
    cbeam_results = []
    cbush_results = []
    cbush1d_results = []
    plate_results = []
    comp_results = []
    solid_results = []

    if check_rod_stress and len(crod_eids):
        # log.debug('check_rod_stress')
        rod_stress_list = _fill_rod_list(
            model_results,
            crod_eids, ctube_eids, conrod_eids,
            is_stress=True)
        assert len(rod_stress_list), 'No rod stress results found'
        # log.info('_envelope_stress_strain')
        _envelope_stress_strain(
            rod_stress, rod_stress_list,
            all_eids,
            all_eids_list, all_subcases_list, rod_results, dtype=dtype)
        # log.info('_envelope_stress_strain_end')

    if check_rod_strain and len(crod_eids):
        # log.debug('check_rod_strain')
        rod_strain_list = _fill_rod_list(
            model_results,
            crod_eids, ctube_eids, conrod_eids,
            is_stress=True)
        assert len(rod_strain_list), 'No rod strain results found'
        _envelope_stress_strain(
            rod_strain, rod_strain_list,
            all_eids,
            all_eids_list, all_subcases_list, rod_results, dtype=dtype)

    if check_bar_stress and len(cbar_eids):
        bar_stress_list = _fill_bar_list(
            model_results, cbar_eids, is_stress=True)
        assert len(bar_stress_list), 'No bar stress results found'
        _envelope_stress_strain(
            bar_stress, bar_stress_list,
            all_eids,
            all_eids_list, all_subcases_list, cbar_results, dtype=dtype)

    if check_bar_strain and len(cbar_eids):
        bar_strain_list = _fill_bar_list(
            model_results, cbar_eids, is_stress=True)
        assert len(bar_strain_list), 'No bar strain results found'
        _envelope_stress_strain(
            bar_strain, bar_strain_list,
            all_eids,
            all_eids_list, all_subcases_list, cbar_results, dtype=dtype)

    if check_beam_stress and len(cbeam_eids):
        beam_stress_list = _fill_beam_list(
            model_results, cbeam_eids, is_stress=True)
        assert len(beam_stress_list), 'No beam stress results found'
        _envelope_stress_strain(
            beam_stress, beam_stress_list,
            all_eids,
            all_eids_list, all_subcases_list, cbeam_results, dtype=dtype)

    if check_beam_strain and len(cbeam_eids):
        beam_strain_list = _fill_beam_list(
            model_results, cbar_eids, is_stress=True)
        assert len(beam_strain_list), 'No beam strain results found'
        _envelope_stress_strain(
            beam_strain, beam_strain_list,
            all_eids,
            all_eids_list, all_subcases_list, cbeam_results, dtype=dtype)

    if check_comp_plate_stress and len(all_comp_plate_eids):
        # log.debug('check_comp_plate_stress')
        comp_plate_stress_list = _fill_comp_plate_list(
            model_results,
            ctria3_comp_eids, cquad4_comp_eids,
            ctria6_comp_eids, cquad8_comp_eids,
            ctriar_comp_eids, cquadr_comp_eids,
            is_stress=True)
        assert len(comp_plate_stress_list), 'No comp plate stress results found'
        _envelope_stress_strain(
            comp_plate_stress,
            comp_plate_stress_list,
            all_eids,
            all_eids_list, all_subcases_list, comp_results, dtype=dtype)

    if check_comp_plate_strain and len(all_comp_plate_eids):
        # log.debug('check_comp_plate_strain')
        comp_plate_strain_list = _fill_comp_plate_list(
            model_results,
            ctria3_comp_eids, cquad4_comp_eids,
            ctria6_comp_eids, cquad8_comp_eids,
            ctriar_comp_eids, cquadr_comp_eids,
            is_stress=False)
        assert len(comp_plate_strain_list), 'No comp plate strain results found'
        _envelope_stress_strain(
            comp_plate_strain,
            comp_plate_strain_list,
            all_eids,
            all_eids_list, all_subcases_list, comp_results, dtype=dtype)

    if check_plate_stress and len(all_plate_eids):
        # log.debug('check_plate_stress')
        plate_stress_list = _fill_plate_list(
            model_results,
            ctria3_plate_eids, cquad4_plate_eids,
            ctria6_plate_eids, cquad8_plate_eids,
            ctriar_plate_eids, cquadr_plate_eids,
            is_stress=True)
        assert len(plate_stress_list), 'No plate stress results found'
        _envelope_solid_stress_strain(
            plate_stress,
            plate_stress_list,
            all_eids,
            all_eids_list, all_subcases_list, plate_results,
            consider_corner_nodes=consider_plate_nodes, dtype=dtype)

    if check_plate_strain and len(all_plate_eids):
        # log.debug('check_plate_strain')
        plate_strain_list = _fill_plate_list(
            model_results,
            ctria3_plate_eids, cquad4_plate_eids,
            ctria6_plate_eids, cquad8_plate_eids,
            ctriar_plate_eids, cquadr_plate_eids,
            is_stress=False)
        assert len(plate_strain_list), 'No plate strain results found'
        _envelope_solid_stress_strain(
            plate_strain,
            plate_strain_list,
            all_eids,
            all_eids_list, all_subcases_list, plate_results,
            consider_corner_nodes=consider_plate_nodes, dtype=dtype)

    if check_solid_stress and len(all_solid_eids):
        # log.debug('check_solid_stress')
        solid_stress_list = _fill_solid_list(
            model_results,
            ctetra_eids, cpenta_eids, chexa_eids, cpyram_eids,
            is_stress=True)
        assert len(solid_stress_list), 'No solid stress results found'
        _envelope_solid_stress_strain(
            solid_stress,
            solid_stress_list,
            all_eids,
            all_eids_list, all_subcases_list, solid_results,
            consider_corner_nodes=consider_solid_nodes)

    if check_solid_strain and len(all_solid_eids):
        # log.debug('check_solid_strain')
        solid_strain_list = _fill_solid_list(
            model_results,
            ctetra_eids, cpenta_eids, chexa_eids, cpyram_eids,
            is_stress=False)
        assert len(solid_strain_list), 'No solid strain results found'
        _envelope_solid_stress_strain(
            solid_strain,
            solid_strain_list,
            all_eids,
            all_eids_list, all_subcases_list, solid_results,
            consider_corner_nodes=consider_solid_nodes)

    results = {
        'rod': (rod_results, all_rod_eids, is_crod_min),
        # 'crod': (crod_results, crod_eids, is_crod_min),
        # 'ctube': (ctube_results, ctube_eids, is_ctube_min),
        # 'conrod': (conrod_results, conrod_eids, is_conrod_min),
        'cbar': (cbar_results, cbar_eids, is_cbar_min),
        'cbeam': (cbeam_results, cbeam_eids, is_cbeam_min),
        'cbush': (cbush_results, cbush_eids, is_cbush_min),
        'cbush1d': (cbush1d_results, cbush1d_eids, is_cbush1d_min),

        'plate': (plate_results, all_plate_eids, is_plate_min),
        'comp_plate': (comp_results, all_comp_plate_eids, is_comp_plate_min),
        'solid': (solid_results, all_solid_eids, is_solid_min),
    }
    out_subcases, combined_data = _envelope_post(
        all_subcases_list, all_eids,
        results,
        percent_eids_target, model, log)
    return out_subcases

def get_all_eids(eids_tuple: tuple[np.ndarray, ...],
                 sort: bool=False) -> np.ndarray:
    eids_list = [eids for eids in eids_tuple if len(eids)]
    if len(eids_list) == 0:
        return np.array([], dtype='int32')
    all_eids = np.hstack(eids_list)
    if sort:
        all_eids.sort()
    return all_eids

def _envelope_post(all_subcases_list: list[int],
                   all_eids: np.ndarray,
                   results: dict[str, tuple[np.ndarray, np.ndarray, bool]],
                   percent_eids_target: float,
                   model: BDF,
                   log: SimpleLogger) -> tuple[np.ndarray, Any]:
    all_subcases = np.unique(all_subcases_list)
    nsubcase_all = len(all_subcases)
    neid_all = len(all_eids)
    assert nsubcase_all > 0 and neid_all > 0, (nsubcase_all, neid_all)

    comp_plate_combined_data = np.array([])
    solid_combined_data = np.array([])

    icase_critical_list = []
    icase_critical = np.full(neid_all, -1, dtype='int32')

    rod_results, rod_eid, is_rod_min = results['rod']
    if len(rod_results):
        icase_criticali, data_critical, comp_plate_combined_data = _get_combined_data(
            all_subcases, 'rods', rod_results,
            nsubcase_all, neid_all, is_rod_min)
        icase_critical_list.append(icase_criticali)
        icase = np.where(icase_criticali != -1)[0]
        icase_critical[icase] = icase_criticali[icase]
    del rod_results, rod_eid, is_rod_min

    plate_results, plate_eid, is_plate_min = results['plate']
    if len(plate_results):
        icase_criticali, data_critical, comp_plate_combined_data = _get_combined_data(
            all_subcases, 'plates', plate_results,
            nsubcase_all, neid_all, is_plate_min)
        icase_critical_list.append(icase_criticali)
        icase = np.where(icase_criticali != -1)[0]
        icase_critical[icase] = icase_criticali[icase]
    del plate_results, plate_eid, is_plate_min

    comp_results, comp_plate_eid, is_comp_plate_min = results['comp_plate']
    if len(comp_results):
        icase_criticali, data_critical, comp_plate_combined_data = _get_combined_data(
            all_subcases, 'comp_plates', comp_results,
            nsubcase_all, neid_all, is_comp_plate_min)
        icase_critical_list.append(icase_criticali)
        icase = np.where(icase_criticali != -1)[0]
        icase_critical[icase] = icase_criticali[icase]
    del comp_results, comp_plate_eid, is_comp_plate_min

    solid_results, solid_eids, is_solid_min = results['solid']
    if len(solid_results):
        icase_criticali, data_critical, solid_combined_data = _get_combined_data(
            all_subcases, 'solids', solid_results,
            nsubcase_all, neid_all, is_solid_min)
        icase_critical_list.append(icase_criticali)
        icase = np.where(icase_criticali != -1)[0]
        icase_critical[icase] = icase_criticali[icase]
    del solid_results, solid_eids, is_solid_min

    # icase_critical = np.vstack(icase_critical_list)
    log.debug(f'icase_critical = {icase_critical}')
    ncritical_cases = len(icase_critical)
    assert ncritical_cases == neid_all, (ncritical_cases, neid_all)

    # potential for nan in icase_max if no result
    #   filter the -1 cases
    imissing_eids = (icase_critical == -1)
    if imissing_eids.sum():
        eids_missing = all_eids[imissing_eids]
        etypes_missing = np.array(list(set([model.elements[eid].type for eid in eids_missing])))
        etypes_missing.sort()
        log.warning(f'eids={eids_missing} not found; etypes_missing={etypes_missing}')

    icase_reduced = np.where(~imissing_eids)[0]
    icase_critical = icase_critical[icase_reduced]
    log.info(f'icase_critical (reduced) = {icase_critical}')
    assert len(icase_critical) > 0, 'No results'

    # get the critical subcases
    subcase_critical = all_subcases[icase_critical]
    log.info(f'subcase_critical = {subcase_critical}')

    usubcase_critical, counts = np.unique(subcase_critical, return_counts=True)
    # print(f'usubcase_critical = {usubcase_critical}')
    # print(f'counts = {counts}')

    isort = np.argsort(counts)
    counts = counts[isort]
    usubcase_critical = usubcase_critical[isort]
    log.info(f'usubcase_critical = {usubcase_critical.tolist()}')
    log.info(f'           counts = {counts.tolist()}')

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

def _get_combined_data(all_subcases: np.ndarray,
                       result_group: str, results: list[np.ndarray],
                       nsubcase_all: int, neid_all: int,
                       is_min: bool,
                       dtype: str='float32') -> np.ndarray:
    assert len(all_subcases) == len(np.unique(all_subcases))

    # print('----------------')
    # print(result_group)
    # print(f'nsubcase_all={nsubcase_all} neid_all={neid_all}')
    combined_data = np.full((nsubcase_all, neid_all), np.nan, dtype=dtype)
    if len(results) == 1:
        obj_name, subcases, eids, ieid, data = results[0]
        isubcase = np.searchsorted(all_subcases, subcases)
        combined_data[isubcase, ieid] = data
        # print('', obj_name, eids, ieid)
        # print(combined_data)
    else:
        assert len(results)
        for obj_name, subcases, eids, ieid, data in results:
            # print('', obj_name, eids, ieid)
            # print(data)
            isubcase = np.searchsorted(all_subcases, subcases)
            combined_data[isubcase, ieid] = data

    # find the driving case (and data) for each element
    icase_critical = np.full(neid_all, -1, dtype='int32')
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
    return icase_critical, data_critical, combined_data

def get_elements_dict(model: BDF) -> tuple[dict[str, list[int]],
                                           dict[str, list[int]]]:
    log = model.log
    eid_to_nid_map = {}
    etype_to_eids = defaultdict(list)
    etype_ptype_to_eids = defaultdict(list)
    for eid, elem in model.elements.items():
        eid_to_nid_map[eid] = elem.nodes
        prop = elem.pid_ref
        etype_to_eids[elem.type].append(eid)
        prop_type = prop.type
        if prop_type in {'PBARL', 'PBEAML', 'PCOMPG'}:
            prop_type = prop_type[:-1]
        etype_ptype_to_eids[(elem.type, prop_type)].append(eid)
    return dict(etype_to_eids), dict(etype_ptype_to_eids)

def _fill_rod_list(model_results: OP2,
                   crod_eids: np.ndarray, ctube_eids: np.ndarray, conrod_eids: np.ndarray,
                   is_stress: bool) -> list:
    if is_stress:
        word = 'stress'
        stress = model_results.op2_results.stress
        crod_obj = stress.crod_stress
        ctube_obj = stress.ctube_stress
        conrod_obj = stress.conrod_stress
    else:
        word = 'strain'
        strain = model_results.op2_results.strain
        crod_obj = strain.crod_strain
        ctube_obj = strain.ctube_strain
        conrod_obj = strain.conrod_strain

    rod_list = []
    # print('eids-rods', crod_eids, ctube_eids, conrod_eids)
    if len(crod_eids):
        rod_list.append((crod_eids, f'crod_{word}', crod_obj))
    if len(ctube_eids):
        rod_list.append((ctube_eids, f'ctube_{word}', ctube_obj))
    if len(conrod_eids):
        rod_list.append((conrod_eids, f'conrod_{word}', conrod_obj))
    return rod_list

def _fill_bar_list(model_results: OP2,
                   cbar_eids: np.ndarray,
                   is_stress: bool) -> list:
    if is_stress:
        word = 'stress'
        stress = model_results.op2_results.stress
        cbar_obj = stress.cbar_stress
    else:
        word = 'strain'
        strain = model_results.op2_results.strain
        cbar_obj = strain.cbar_strain

    bar_list = []
    if len(cbar_eids):
        bar_list.append((cbar_eids, f'cbar_{word}', cbar_obj))
    return bar_list

def _fill_beam_list(model_results: OP2,
                    cbeam_eids: np.ndarray,
                    is_stress: bool) -> list:
    if is_stress:
        word = 'stress'
        stress = model_results.op2_results.stress
        cbeam_obj = stress.cbeam_stress
    else:
        word = 'strain'
        strain = model_results.op2_results.strain
        cbeam_obj = strain.cbeam_strain

    beam_list = []
    if len(cbeam_eids):
        beam_list.append((cbeam_eids, f'cbeam_{word}', cbeam_obj))
    return beam_list

def _fill_plate_list(model_results: OP2,
                     ctria3_eids: np.ndarray, cquad4_eids: np.ndarray,
                     ctria6_eids: np.ndarray, cquad8_eids: np.ndarray,
                     ctriar_eids: np.ndarray, cquadr_eids: np.ndarray,
                     is_stress: bool) -> list:
    if is_stress:
        word = 'stress'
        stress = model_results.op2_results.stress
        ctria3_obj = stress.ctria3_stress
        cquad4_obj = stress.cquad4_stress
        ctria6_obj = stress.ctria6_stress
        cquad8_obj = stress.cquad8_stress
        ctriar_obj = stress.ctriar_stress
        cquadr_obj = stress.cquadr_stress
    else:
        word = 'strain'
        strain = model_results.op2_results.strain
        ctria3_obj = strain.ctria3_strain
        cquad4_obj = strain.cquad4_strain
        ctria6_obj = strain.ctria6_strain
        cquad8_obj = strain.cquad8_strain
        ctriar_obj = strain.ctriar_strain
        cquadr_obj = strain.cquadr_strain

    plate_list = []
    if len(ctria3_eids):
        plate_list.append((ctria3_eids, f'ctria3_{word}', ctria3_obj))
    if len(cquad4_eids):
        plate_list.append((cquad4_eids, f'cquad4_{word}', cquad4_obj))
    if len(ctria6_eids):
        plate_list.append((ctria6_eids, f'ctria6_{word}', ctria6_obj))
    if len(cquad8_eids):
        plate_list.append((cquad8_eids, f'cquad8_{word}', cquad8_obj))
    if len(ctriar_eids):
        plate_list.append((ctriar_eids, f'ctriar_{word}', ctriar_obj))
    if len(cquadr_eids):
        plate_list.append((cquadr_eids, f'cquadr_{word}', cquadr_obj))
    return plate_list

def _fill_comp_plate_list(model_results: OP2,
                          ctria3_eids: np.ndarray, cquad4_eids: np.ndarray,
                          ctria6_eids: np.ndarray, cquad8_eids: np.ndarray,
                          ctriar_eids: np.ndarray, cquadr_eids: np.ndarray,
                          is_stress: bool) -> list:
    if is_stress:
        word = 'stress'
        stress = model_results.op2_results.stress
        ctria3_obj = stress.ctria3_composite_stress
        cquad4_obj = stress.cquad4_composite_stress
        ctria6_obj = stress.ctria6_composite_stress
        cquad8_obj = stress.cquad8_composite_stress
        ctriar_obj = stress.ctriar_composite_stress
        cquadr_obj = stress.cquadr_composite_stress
    else:
        word = 'strain'
        strain = model_results.op2_results.strain
        ctria3_obj = strain.ctria3_composite_strain
        cquad4_obj = strain.cquad4_composite_strain
        ctria6_obj = strain.ctria6_composite_strain
        cquad8_obj = strain.cquad8_composite_strain
        ctriar_obj = strain.ctriar_composite_strain
        cquadr_obj = strain.cquadr_composite_strain

    plate_list = []
    if len(ctria3_eids):
        plate_list.append((ctria3_eids, f'ctria3_{word}', ctria3_obj))
    if len(cquad4_eids):
        plate_list.append((cquad4_eids, f'cquad4_{word}', cquad4_obj))
    if len(ctria6_eids):
        plate_list.append((ctria6_eids, f'ctria6_{word}', ctria6_obj))
    if len(cquad8_eids):
        plate_list.append((cquad8_eids, f'cquad8_{word}', cquad8_obj))
    if len(ctriar_eids):
        plate_list.append((ctriar_eids, f'ctriar_{word}', ctriar_obj))
    if len(cquadr_eids):
        plate_list.append((cquadr_eids, f'cquadr_{word}', cquadr_obj))
    return plate_list

def _fill_solid_list(model_results: OP2,
                     ctetra_eids: np.ndarray, cpenta_eids: np.ndarray,
                     chexa_eids: np.ndarray, cpyram_eids: np.ndarray,
                     is_stress: bool) -> list:
    if is_stress:
        word = 'stress'
        stress = model_results.op2_results.stress
        ctetra_obj = stress.ctetra_stress
        cpenta_obj = stress.cpenta_stress
        chexa_obj = stress.chexa_stress
        cpyram_obj = stress.cpyram_stress
    else:
        word = 'strain'
        strain = model_results.op2_results.strain
        ctetra_obj = strain.ctetra_strain
        cpenta_obj = strain.cpenta_strain
        chexa_obj = strain.chexa_strain
        cpyram_obj = strain.cpyram_strain

    solid_list = []
    if len(ctetra_eids):
        solid_list.append((ctetra_eids, f'ctetra_{word}', ctetra_obj))
    if len(chexa_eids):
        solid_list.append((chexa_eids, f'chexa_{word}', chexa_obj))
    if len(cpenta_eids):
        solid_list.append((cpenta_eids, f'cpenta_{word}', cpenta_obj))
    if len(cpyram_eids):
        solid_list.append((cpyram_eids, f'cpyram_{word}', cpyram_obj))
    return solid_list

def _envelope_solid_stress_strain(result_name: str,
                                  solid_stress_list: list[tuple[np.ndarray, Any]],
                                  all_eids: np.ndarray,
                                  all_eids_list: list[list[int]],
                                  all_subcases_list: list[int],
                                  results: list[np.ndarray],
                                  consider_corner_nodes: bool=True,
                                  dtype: str='float32'):
    for (eids, obj_name, obj_dict) in solid_stress_list:
        subcases = list(obj_dict)
        nsubcase_obj = len(subcases)
        if nsubcase_obj == 0:
            warnings.warn(f'no results for {obj_name}')
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
                    raise RuntimeError(f'missing solid eids={missing_eids.tolist()}')
                assert np.array_equal(all_eids[ieid], eids)
            datai = obj.envelope(eids, result_name, consider_corner_nodes)
            solid_data[isubcase, :] = datai
        # print(f'subcases = {subcases}')
        results.append((obj_name, subcases, eids, ieid, solid_data))
        all_eids_list.append(eids)
        all_subcases_list.extend(subcases)

def _envelope_stress_strain(result_name: str,
                            comp_plate_list: list[tuple[np.ndarray, Any]],
                            all_eids: np.ndarray,
                            all_eids_list: list[list[int]],
                            all_subcases_list: list[int],
                            results: list[np.ndarray],
                            dtype: str='float32'):
    """no consideration for corner nodes"""
    for (eids, obj_name, obj_dict) in comp_plate_list:
        subcases = list(obj_dict)
        nsubcase_obj = len(subcases)
        if nsubcase_obj == 0:
            warnings.warn(f'no results for {obj_name}')
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
                    raise RuntimeError(f'missing composite plate eids={missing_eids.tolist()}')
                assert np.array_equal(all_eids[ieid], eids)
            datai = obj.envelope(eids, result_name)
            comp_data[isubcase, :] = datai
        # print(f'subcases = {subcases}')

        # subcases, eids, ieid, combined_data
        results.append((obj_name, subcases, eids, ieid, comp_data))
        all_eids_list.append(eids)
        all_subcases_list.extend(subcases)
