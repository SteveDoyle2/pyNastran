from copy import deepcopy
import numpy as np
from pyNastran.dev.bdf_vectorized3.bdf import BDF


def remove_unused(model: BDF, inplace: bool=False) -> BDF:
    if not inplace:
        model = deepcopy(model)
    log = model.log
    cards = [card for card in model._cards_to_setup if card.n > 0]
    used_dict = {
        'coord_id': [0],
        'spoint_id': [],
        'point_id': [],
        'node_id': [],
        'element_id': [],

        'property_id': [],
        'pelast_id': [],
        'pdampt_id': [],
        'pbusht_id': [],

        'material_id': [],
        'spc_id': [],          # SPCADD -> SPC
        'mpc_id': [],          # MPCADD -> MPC
        'load_id': [],         # LOAD -> FORCE, FORCE1
        'suport_id': [0],      # SUPORT, SUPORT1
        'delay_id': [],        # TLOAD1/2 -> DELAY
        'dphase_id': [],       # RLOAD1/2 -> DPHASE
        'tf_id': [],       # TF

        'desvar_id': [],       # DVPREL1/DVPREL2 -> DESVAR
        'dconstr_id': [],      # DCONADD -> DCONSTR
        'dvprel_id': [],       # DVPREL2 -> DVPREL1
        'dvmrel_id': [],       # DVMREL2 -> DVMREL1
        'dvcrel_id': [],       # DVCREL2 -> DVCREL1
        'dresp_id': [],        # DVPREL1/DVPREL2 -> DRESP1/DRESP2
        'ddval_id': [],
        'deqatn_id': [],

        'tablem_id': [],       # MATTx  -> TABLEMx
        'tabled_id': [0],      # TLOADx -> TABLEDx
        'contact_set_id': [],  # BCTADD -> BCTSET
        'glue_id': [],         # BGADD  -> BGSET
        'contact_id': [],      # BGSET/BCTSET -> BSURF/BSURFS
        'bfric_id': [],        # BCONP -> BFRIC
        'bolt_id': [],
    }
    #for card in model._cards_to_setup:
        #print(card)


    for subcase_id, subcase in model.subcases.items():
        if 'SPC' in subcase:
            spc_id, _options = subcase['SPC']
            used_dict['spc_id'].append(spc_id)
            del spc_id, _options
        if 'MPC' in subcase:
            mpc_id, _options = subcase['MPC']
            used_dict['mpc_id'].append(mpc_id)
            del mpc_id, _options
        if 'SUPORT' in subcase:
            suport_id, _options = subcase['SUPORT']
            used_dict['suport_id'].append(suport_id)
            del suport_id, _options
        if 'TFL' in subcase:
            tf_id, _options = subcase['TFL']
            used_dict['tf_id'].append(tf_id)
            del tf_id, _options

    set_cards_to_skip = {
        'NSM1', 'NSML1'}
    for card in cards:
        if card.type in set_cards_to_skip:
            continue

        #print(card.type)
        if hasattr(card, 'set_used'):
            card.set_used(used_dict)
        else:
            #card.set_used(used_dict)
            log.warning(f'{card.type} does not support set_used for remove_unused')
            #raise RuntimeError(card.type)
        used_arrays = to_dict_array(card, used_dict)
        assert -1 not in used_arrays['coord_id'], card.type
        #for key, datai in used_arrays.items():
            #if len(datai):
                #print(card.type, key, datai)
        #print('')
    del used_dict

    if len(cards) == 0:
        return model

    coord_id = used_arrays['coord_id']
    #property_id = used_arrays['property_id']
    assert len(coord_id) > 0, coord_id
    not_supported_cards = {'DLINK', 'POINT', 'BEDGE', 'BCONP', 'BFRIC', 'BGSET', 'BGADD'}
    for card in cards:
        if hasattr(card, 'remove_unused'):
            card.remove_unused(used_arrays)
        elif card.type in not_supported_cards:
            log.error(f'{card.type} does not support remove_unused')
        else:
            log.error(f'{card.type} does not support remove_unused')
            #raise RuntimeError(card.type)
    return model

def to_dict_array(card, used_dict: dict[str, list[np.ndarray]]) -> dict[str, np.ndarray]:
    assert isinstance(used_dict['coord_id'], list), card.type
    used_arrays = {}
    for key, used_list in used_dict.items():
        if len(used_list) == 0:
            used_arrays[key] = np.array([], dtype='int32')
            continue
        try:
            values_stacked = np.hstack(used_list)
        except ValueError:  # pragma: no cover
            raise RuntimeError((card.type, key))

        try:
            values = np.unique(values_stacked)
        except ValueError:  # pragma: no cover
            raise RuntimeError((card.type, key))

        min_value = 1
        if key == 'coord_id':
            min_value = -1
        elif key == 'suport_id':
            min_value = 0
        elif key == 'node_id':
            min_value = 0
        elif key == 'tabled_id':
            min_value = 0
        elif key == 'tablem_id':
            min_value = 0

        if len(values):
            assert values.min() >= min_value, (card.type, key, values.min())
        used_arrays[key] = values

    coord_id = used_arrays['coord_id']
    assert coord_id.min() >= -1, (card.type, coord_id.min())
    return used_arrays
