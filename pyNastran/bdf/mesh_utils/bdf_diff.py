import os
from collections import defaultdict
import numpy as np
from cpylog import SimpleLogger
from pyNastran.utils import PathLike
from pyNastran.bdf.bdf import read_bdf, BDF, CaseControlDeck, Subcase

scalar_obj_keys: list[str] = [
    'aero', 'aeros',
    #'axic', 'axif', # removed
    'cyax', 'baror', 'beamor',
    'acmodl', 'modtrak',
    'doptprm',
    'dtable', 'grdset', 'radset', 'seqgp',
    'sol',
    # 'zona',
]
dict_cards = [
    'nodes', 'elements', 'properties',
    'pbusht', 'pdampt', 'pelast', 'properties_mass',
    'methods', 'cMethods',
    'materials', 'thermal_materials', 'creep_materials', 'hyperelastic_materials',

    # aero
    'caeros', 'paeros', 'splines', 'sets', 'flfacts',
    'trims', 'csschds', 'gusts',
    'aecomps', 'aefacts', 'aelists', 'aeparams',
    'aestats', 'aesurf', 'aesurfs',
    'divergs', 'dlinks',

    # contact
    'bconp', 'bcrparas', 'bctadds',
    'bctparas', 'bctsets', 'blseg', 'bsurf', 'bsurfs',
    'bgadds', 'bgsets', 'bctparms',
    'bfric', 'bcbodys',

    # optimization
    'desvars',

    # other
    'ao_element_flags',
    'dareas',
    'dmig', 'dmiji', 'dmij', 'dmik', 'dmi', 'dmiax', 'dti',
    'dphases', 'delays',
    'epoints',
    #'gridb',  # removed
    'nlparms', 'nlpcis',
    'normals',
    'nxstrats',
    'phbdys', 'points',
    'radcavs', 'radmtx',
    'rotors',
    'spcoffs', 'spoints', 'suport1',
    'tables_d', 'tables_m', 'tables_sdamping', 'random_tables',
    'tempds', 'tics',
    'tstepnls', 'tsteps',
    'view3ds', 'views', 'convection_properties',

    'MATT1', 'MATT2', 'MATT3',
    'MATT4', 'MATT5', 'MATT8', 'MATT9', 'MATT11',
    'MATS1', 'MATS3', 'MATS8',
]
seen = set()
dupes = [x for x in dict_cards if x in seen or seen.add(x)]
assert len(dupes) == 0, dupes
del seen, dupes

dict_int_list_obj_attrs: list[str] = [
    'spcs', 'spcadds',
    'mpcs', 'mpcadds',
    'loads', 'load_combinations',
    'dloads', 'dload_entries',
    # 'usets', # has string keys
    'nsms', 'nsmadds',
    'frequencies',
    'bcs', 'transfer_functions',
    'dvgrids',

    # parametric
    'pval',
]


def get_diff_bdfs(bdf_filename1: PathLike, bdf_filename2: PathLike,
                  added_bdf_filename: PathLike='',
                  removed_bdf_filename: PathLike='',
                  log=None) -> tuple[bool, bool, BDF, BDF]:
    """
    diffs two bdfs

    doesn't consider:
     - constraints
     - other things that are lists (e.g., monitor_points) or unsortable
    """
    base1, ext1 = os.path.splitext(bdf_filename1)
    base2, ext2 = os.path.splitext(bdf_filename2)
    if added_bdf_filename == '':
        added_bdf_filename = base1 + '.added' + ext1
    if removed_bdf_filename == '':
        removed_bdf_filename = base2 + '.removed' + ext2

    old_model = read_bdf(bdf_filename1, xref=False, log=log)
    new_model = read_bdf(bdf_filename2, xref=False, log=log)

    removed_model = BDF(log=log)
    added_model = BDF(log=log)

    # dict_int_list_obj_attrs: list[str] = [
    #     'spcs', 'spcadds',
    #     'mpcs', 'mpcadds',
    #     'loads', 'load_combinations',
    #     'dloads', 'dload_entries',
    #     # 'usets', # has string keys
    #     'nsms', 'nsmadds',
    #     'frequencies',
    #     'bcs', 'transfer_functions',
    #     'dvgrids',
    #
    #     # parametric
    #     'pval',
    # ]

    added_cards = False
    removed_cards = False
    added_cards, removed_cards = _save_case_control(
        old_model, new_model,
        added_model, removed_model,
        added_cards, removed_cards,
    )
    added_cards, removed_cards = _save_scalars(
        old_model, new_model,
        added_model, removed_model,
        scalar_obj_keys,
        added_cards, removed_cards,
    )
    added_cards, removed_cards = _save_dict_list_cards(
        old_model, new_model,
        added_model, removed_model,
        dict_int_list_obj_attrs,
        added_cards, removed_cards,
    )
    added_cards, removed_cards = _save_dict_cards(
        old_model, new_model,
        added_model, removed_model,
        dict_cards,
        added_cards, removed_cards,
    )
    if not added_cards:
        log.warning(f'added_bdf={added_bdf_filename} is empty...')
    added_model.write_bdf(added_bdf_filename)
    if not removed_cards:
        log.warning(f'removed_bdf={removed_bdf_filename} is empty...')
    removed_model.write_bdf(removed_bdf_filename)
    # bdf_merge(bdf_filenames, bdf_filename_out, renumber=True,
    #           encoding=None, size=size, is_double=False, cards_to_skip=cards_to_skip,
    #           log=log)
    return added_cards, removed_cards, added_model, removed_model

def _save_case_control(old_model: BDF, new_model: BDF,
                       added_model: BDF, removed_model: BDF,
                       added: bool, removed: bool) -> tuple[bool, bool]:
    log = old_model.log
    old_subcases = old_model.case_control_deck.subcases
    new_subcases = new_model.case_control_deck.subcases
    isubcases = np.unique(np.array(list(old_subcases) + list(new_subcases)))
    added_subcases = {}
    removed_subcases = {}
    for isubcase in isubcases:
        old_subcase = old_subcases[isubcase]
        new_subcase = new_subcases[isubcase]
        added_subcase = Subcase(id=isubcase)
        removed_subcase = Subcase(id=isubcase)
        if old_subcase == new_subcase:
            #log.info(f'same subcase; isubcase={isubcase}')
            #print(old_subcase)
            #print(new_subcase)
            #print('-'*40)
            continue
        #print(old_subcase.params)
        old_params = old_subcase.params
        new_params = new_subcase.params
        keys = np.unique(np.array(list(old_params) + list(new_params)))
        added_params = {}
        removed_params = {}
        for key in keys:
            if key not in old_params:
                #log.debug(f'{isubcase}: added   {key} param: {new_params[key]}')
                added_params[key] = new_params[key]
                continue
            if key not in new_params:
                #log.debug(f'{isubcase}: removed {key} param: {old_params[key]}')
                removed_params[key] = old_params[key]
                continue
            old_param = old_params[key]
            new_param = new_params[key]
            if old_param == new_param:
                continue
            if str(new_param) == str(old_param):
                #print(new_param)
                continue

            added_params[key] = new_param
            removed_params[key] = old_param
            #print(type(old_param), type(new_param))
            #log.debug(f'{isubcase}: changed {key} param: {old_param!r} {new_param!r}')
            #log.debug(f'  old={old_param!r} -> {new_param!r}')

        if added_params:
            #print(isubcase, type(isubcase), type(added_subcases))
            added_subcase.params = added_params
            added_subcases[isubcase] = added_subcase
        if removed_params:
            #print('removed_params =', removed_params)
            removed_subcase.params = removed_params
            removed_subcases[isubcase] = removed_subcase
            #log.debug(f'{isubcase}: changed isubcase={isubcase}:\n{str(removed_subcase)}')

    if len(added_subcases):
        added = True
        #print(f'added_subcases = {list(added_subcases)}')
        added_case_control = CaseControlDeck([])
        if 0 not in added_subcases:
            added_subcases[0] = Subcase(id=0)
        added_case_control.subcases = added_subcases
        added_model.case_control_deck = added_case_control
    if len(removed_subcases):
        removed = True
        #print(f'removed_subcases = {list(removed_subcases)}')
        removed_case_control = CaseControlDeck([])
        if 0 not in removed_subcases:
            removed_subcases[0] = Subcase(id=0)
        removed_case_control.subcases = removed_subcases
        removed_model.case_control_deck = removed_case_control

    return added, removed

def _save_scalars(old_model: BDF, new_model: BDF,
                  added_model: BDF, removed_model: BDF,
                  scalars: list[str],
                  added_cards: bool, removed_cards: bool) -> tuple[bool, bool]:
    for card in scalars:
        old_value = getattr(old_model, card)
        new_value = getattr(new_model, card)
        #added_group = getattr(added_model, card)
        #removed_group = getattr(removed_model, card)
        if old_value is None and new_value is None:
            continue
        elif old_value is None:
            #print(f'{card} was added')
            added_cards = True
            setattr(added_model, card, new_value)
        elif new_value is None:
            #print(f'{card} was removed')
            removed_cards = True
            setattr(removed_model, card, old_value)
        elif old_value == new_value:
            continue
        else:
            #print(f'{card} was added/removed; {old_value} to {new_value}')
            #print()
            setattr(added_model, card, new_value)
            setattr(removed_model, card, old_value)
            added_cards = True
            removed_cards = True
    return added_cards, removed_cards

def _save_dict_cards(old_model: BDF, new_model: BDF,
                     added_model: BDF, removed_model: BDF,
                     dict_cards: list[str],
                     added_cards: bool, removed_cards: bool) -> tuple[bool, bool]:

    for card in dict_cards:
        old_group = getattr(old_model, card)
        new_group = getattr(new_model, card)
        added_group = getattr(added_model, card)
        removed_group = getattr(removed_model, card)

        all_ids = np.unique(list(old_group.keys()) +
                            list(new_group.keys()))
        for idi in all_ids:
            if idi not in old_group:
                #print(f'{card} id={idi} not in old_group; added')
                added_group[idi] = new_group[idi]
            elif idi not in new_group:
                #print(f'{card} id={idi} not in new_group; removed')
                removed_group[idi] = old_group[idi]
            else:
                added_card = new_group[idi]
                removed_card = old_group[idi]
                if added_card == removed_card:
                    #print(f'{card} id={idi} in both')
                    continue
                #print(f'{card} id={idi} in both, but different')
                added_group[idi] = added_card
                removed_group[idi] = removed_card
        if len(added_group):
            added_cards = True
        if len(removed_group):
            removed_cards = True
    return added_cards, removed_cards

def _save_dict_list_cards(old_model: BDF, new_model: BDF,
                          added_model: BDF, removed_model: BDF,
                          dict_list_cards: dict[int, list[str]],
                          added_cards: bool, removed_cards: bool) -> tuple[bool, bool]:

    for card in dict_list_cards:
        old_groups = getattr(old_model, card)
        new_groups = getattr(new_model, card)
        added_groups = getattr(added_model, card)
        removed_groups = getattr(removed_model, card)

        all_ids = np.unique(list(old_groups.keys()) +
                            list(new_groups.keys()))
        for idi in all_ids:
            if idi not in old_groups:
                #print(f'{card} id={idi} not in old_group; added')
                added_groups[idi] = new_groups[idi]  # model.loads
                continue
            if idi not in new_groups:
                #print(f'{card} id={idi} not in new_group; removed')
                removed_groups[idi] = old_groups[idi]
                continue

            added_group = []
            removed_group = []

            # make tuples for duplication, but we're going to add the objects
            old_cards = defaultdict(list)
            new_cards = defaultdict(list)
            for old_card in old_groups[idi]:
                old_cards[tuple(old_card.raw_fields())].append(old_card)
            for new_card in old_groups[idi]:
                new_cards[tuple(new_card.raw_fields())].append(new_card)

            keys = list(set(list(old_cards.keys())).union(set(list(new_cards.keys()))))
            for key in keys:
                if key not in old_cards:
                    for cardi in new_cards[key]:
                        added_group.append(cardi)
                    continue
                if key not in new_cards:
                    for cardi in old_cards[key]:
                        removed_group.append(cardi)
                    continue
                ## TODO: haven't checked for inconsistent lengths
                ##       (e.g., duplicate FORCE cards)
                assert len(old_cards[key]) == len(new_cards[key])

            # key is in both
            if len(added_group):
                added_groups[idi] = added_group
            if len(removed_group):
                removed_groups[idi] = removed_group

        if len(added_groups):
            added_cards = True
        if len(removed_groups):
            removed_cards = True
    return added_cards, removed_cards
