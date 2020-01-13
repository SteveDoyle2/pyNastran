from __future__ import annotations
from typing import List, Set, Dict, Any, Union, TYPE_CHECKING
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


def get_bdf_stats(model: BDF, return_type: str='string',
                  word: str='') -> Union[str, List[str]]:
    """
    Print statistics for the BDF

    Parameters
    ----------
    return_type : str (default='string')
        the output type ('list', 'string')
            'list' : list of strings
            'string' : single, joined string
    word : str; default=''
        model flag

    Returns
    -------
    return_data : str, optional
        the output data

    .. note:: if a card is not supported and not added to the proper
              lists, this method will fail

    .. todo:: RBE3s from OP2s can show up as ???s

    """
    card_dict_groups = [
        'params', 'nodes', 'spoints', 'epoints', 'points', 'gridb',

        'elements', 'ao_element_flags', 'normals', 'rigid_elements', 'plotels',

        'properties', 'pbusht', 'pdampt', 'pelast',
        'properties_mass', 'masses',

        'materials', 'creep_materials', 'hyperelastic_materials',
        'MATT1', 'MATT2', 'MATT3', 'MATT4', 'MATT5', 'MATT8', 'MATT9',
        'MATS1', 'MATS3', 'MATS8', 'MATT8',
        'coords', 'mpcs',

        # axisysmmetric

        # dynamic cards
        'dareas', 'delays', 'dphases', 'nlparms', 'nlpcis',
        'tsteps', 'tstepnls',
        'rotors',

        # direct matrix input - DMIG - dict
        'dmi', 'dmig', 'dmij', 'dmiji', 'dmik', 'dmiax',
        'dequations',
        'transfer_functions',
        'tics',

        # frequencies - dict[List[FREQ]]
        'frequencies',

        # optimization - dict
        'dconadds', 'dconstrs', 'desvars', 'topvar', 'ddvals', 'dlinks', 'dresps',
        'dvcrels', 'dvmrels', 'dvprels', 'dvgrids',

        # SESETx - dict
        'suport1',

        # tables
        'tables', 'tables_d', 'tables_m', 'random_tables', 'tables_sdamping',

        # methods
        'methods', 'cMethods',

        # aero
        'caeros', 'paeros', 'aecomps', 'aefacts', 'aelinks',
        'aelists', 'aeparams', 'aesurf', 'aesurfs', 'aestats', 'gusts', 'flfacts',
        'flutters', 'splines', 'trims', 'divergs', 'csschds',

        # thermal
        'bcs', 'thermal_materials', 'phbdys', 'views', 'view3ds',
        'convection_properties',

        # contact
        'bsurf', 'bsurfs', 'blseg', 'bfric',
        'bconp', 'bcrparas', 'bctadds', 'bctparas', 'bctsets',

        # sets
        'sets', 'usets',

        # superelements
        'csuper', 'csupext',
        'sebulk', 'sebndry', 'seconct', 'seelt', 'seexcld',
        'selabel', 'seloc', 'seload', 'sempln', 'senqset',
        'setree',
        'se_sets', 'se_usets', 'release',

        # ???
        'dscreen', 'dti', 'nxstrats', 'radcavs', 'radmtx', 'ringaxs', 'ringfl',
        'tempds', 'spcoffs',

        # cyclic
        'cyjoin',

        # parametric
        'feedge', 'feface', 'gmcurv', 'gmsurf', 'pset', 'pval',
    ]
    scalar_attrs = [
        'aero', 'aeros', 'grdset', # handled below
        'axic', 'axif', 'cyax', 'modtrak',

        # not handled
        'acmodl',
        'baror', 'beamor', 'doptprm', 'dtable',
        'zona',
    ]

    list_attrs = [
        'asets', 'bsets', 'csets', 'omits', 'qsets',
        'se_bsets', 'se_csets', 'se_qsets',
        'suport', 'se_suport',
        'monitor_points',
    ]
    skip_attrs = [
        'active_filename', 'active_filenames', 'debug', # 'log',
        'reject_lines',
        'is_nx', 'is_msc', 'is_bdf_vectorized', 'dumplines', 'values_to_skip',
        'system_command_lines', 'executive_control_lines', 'case_control_lines',
        'case_control_deck',
        'is_superelements', 'special_cards', 'units',
        'sol', 'sol_iline', 'sol_method', 'cards_to_read', 'card_count',
        'superelement_models', 'wtmass', 'echo', 'force_echo_off',
        'read_includes', 'reject_cards', 'reject_count', 'punch',
        'include_dir', 'include_filenames', 'save_file_structure',
        'rsolmap_to_str', 'nastran_format', 'nid_map', 'bdf_filename',
        'initial_superelement_models',
        'is_zona', 'is_nasa95', 'type_slot_str', 'dict_of_vars', 'code_block',

        # handled below
        'mpcadds', 'mpcs', 'spcadds', 'spcs',
        'loads', 'load_combinations',
        'dloads', 'dload_entries',
        'aero', 'aeros', 'mkaeros',
        'nsmadds', 'nsms',
        'seqgp',

        # unhandled
        'radset',
        'dmigs', 'dmijis', 'dmijs', 'dmiks', 'dmis',


        # vector
        'cbar', 'cbeam', 'cbush',
        'conm2',
        'cshear', 'cvisc',
        'conrod', 'ctube', 'crod',
        'cdamp1', 'cdamp2', 'cdamp3', 'cdamp4',
        'celas1', 'celas2', 'celas3', 'celas4',
        'chexa20', 'chexa8', 'cpenta15', 'cpenta6',
        'cpyram13', 'cpyram5', 'ctetra10', 'ctetra4',
        'cquad', 'cquad4', 'cquad8', 'cquadr',
        'ctria3', 'ctria6', 'ctriar',
        'plotel',

        'dampers', 'elements2',
        'bars', 'beams', 'bushes', 'rods', 'shears',
        'shells', 'solids', 'springs',

        'force', 'force1', 'force2', 'grav', 'grid', 'lseq', 'masses2',
        'moment', 'moment1', 'moment2', 'pload', 'pload1', 'pload2', 'pload4',
        'sload', 'spcd',
        'temp', 'tempd',

    ] + list_attrs + card_dict_groups + scalar_attrs
    missed_attrs = []
    for attr in model.object_attributes(filter_properties=True,
                                        keys_to_skip=skip_attrs):
        #if attr in skip_attrs:
            #continue
        missed_attrs.append(attr)
    if model.__class__.__name__ == 'BDF':
        assert missed_attrs == [], missed_attrs


    # These are ignored because they're lists
    #ignored_types = set([
        #'spoints', 'spointi',  # singleton
        #'grdset',  # singleton

        #'spcs',

        #'suport', 'se_suport', # suport, suport1 - list
        #'doptprm',  # singleton

        ## SETx - list
        #'sets', 'asets', 'bsets', 'csets', 'qsets',
        #'se_bsets', 'se_csets', 'se_qsets',
    #])

    ## TODO: why are some of these ignored?
    #ignored_types2 = set([
        #'case_control_deck', 'caseControlDeck',

        ## done
        #'sol', 'loads', 'mkaeros',
        #'reject_lines', 'reject_cards',

        ## not cards
        #'debug', 'executive_control_lines',
        #'case_control_lines', 'cards_to_read', 'card_count',
        #'is_structured', 'uniqueBulkDataCards',
        #'model_type', 'include_dir',
        #'sol_method', 'log',
        #'sol_iline',
        #'reject_count', '_relpath',
        #'special_cards',])

    #unsupported_types = ignored_types.union(ignored_types2)
    #all_params = object_attributes(model, keys_to_skip=unsupported_types)

    msg = ['---BDF Statistics%s---' % word]
    # sol
    if 'Superelement' not in word:
        msg.append('SOL %s\n' % model.sol)
    msg.extend(_get_bdf_stats_loads(model))

    # load_combinations / loads: handled below
    #handled_explicitly = [
        #'dloads', 'dload_entries',
        #'spcadds', 'spcs',
        #'mpcadds', 'mpcs',
        #'nsmadds', 'nsms',
        #'aero', 'aeros', 'mkaeros', 'radset',
        #'seqgp',
    #]

    # dloads
    for (lid, loads) in sorted(model.dloads.items()):
        msg.append('bdf.dloads[%s]' % lid)
        groups_dict = {}  # type: Dict[str, Any]
        for loadi in loads:
            groups_dict[loadi.type] = groups_dict.get(loadi.type, 0) + 1
        for name, count_name in sorted(groups_dict.items()):
            msg.append('  %-8s %s' % (name + ':', count_name))
        msg.append('')

    for (lid, loads) in sorted(model.dload_entries.items()):
        msg.append('bdf.dload_entries[%s]' % lid)
        groups_dict = {}
        for loadi in loads:
            groups_dict[loadi.type] = groups_dict.get(loadi.type, 0) + 1
        for name, count_name in sorted(groups_dict.items()):
            msg.append('  %-8s %s' % (name + ':', count_name))
        msg.append('')

    _constraint_stats(model, msg)
    _nsm_stats(model, msg)

    _cyclic_stats(model, msg)
    _aero_stats(model, msg)

    # radset
    if model.radset:
        msg.append('bdf:radset')
        msg.append('  %-8s 1' % ('RADSET:'))

    #seqgp
    if model.seqgp:
        msg.append('bdf:seqgp')
        msg.append('  %-8s 1' % ('SEQGP:'))

    for card_group_name in card_dict_groups:
        try:
            card_group = getattr(model, card_group_name)
        except AttributeError:
            msgi = 'cant find card_group_name=%r' % card_group_name
            raise AttributeError(msgi)

        groups = set() # type: Set[str]

        if not isinstance(card_group, dict):
            msgi = '%s is a %s; not dictionary, which is required by get_bdf_stats()' % (
                card_group_name, type(card_group))
            model.log.error(msgi)
            continue
            #raise RuntimeError(msg)

        for card in card_group.values():
            if isinstance(card, list):
                for card2 in card:
                    groups.add(card2.type)
            else:
                groups.add(card.type)

        group_msg = []
        for card_name in sorted(groups):
            try:
                ncards = model.card_count[card_name]
                group_msg.append('  %-8s : %s' % (card_name, ncards))
            except KeyError:
                # we get in here because we used add_grid or similar method, which
                # doesn't increase the card_count, so instead we'll use _type_to_id_map
                counter = '???'
                if card_name in model._type_to_id_map:
                    counter = len(model._type_to_id_map[card_name])
                if card_name == 'CORD2R' and counter == '???':
                    # there is always 1 CORD2R that isn't added to card_count/_type_to_id_map
                    continue
                group_msg.append('  %-8s : %s' % (card_name, counter))
                #assert card_name == 'CORD2R', model.card_count
        if group_msg:
            msg.append('bdf.%s' % card_group_name)
            msg.append('\n'.join(group_msg))
            msg.append('')

    if model.reject_lines:  # List[card]; card = List[str]
        msg.append('Rejected Cards')
        for name, counter in sorted(model.card_count.items()):
            if name not in model.cards_to_read:
                msg.append('  %-8s %s' % (name + ':', counter))
    msg.append('')

    for super_id, superelement in model.superelement_models.items():
        msg += get_bdf_stats(superelement, return_type='list', word=' (Superelement %i)' % super_id)

    if return_type == 'string':
        return '\n'.join(msg)
    return msg

def _constraint_stats(model: BDF, msg: List[str]) -> None:
    """helper for ``get_bdf_stats(...)``"""
    # spcs
    for (spc_id, spcadds) in sorted(model.spcadds.items()):
        msg.append('bdf.spcadds[%s]' % spc_id)
        groups_dict = {}
        for spcadd in spcadds:
            groups_dict[spcadd.type] = groups_dict.get(spcadd.type, 0) + 1
        for name, count_name in sorted(groups_dict.items()):
            msg.append('  %-8s %s' % (name + ':', count_name))
        msg.append('')

    for (spc_id, spcs) in sorted(model.spcs.items()):
        msg.append('bdf.spcs[%s]' % spc_id)
        groups_dict = {}
        for spc in spcs:
            groups_dict[spc.type] = groups_dict.get(spc.type, 0) + 1
        for name, count_name in sorted(groups_dict.items()):
            msg.append('  %-8s %s' % (name + ':', count_name))
        msg.append('')

    # mpcs
    for (mpc_id, mpcadds) in sorted(model.mpcadds.items()):
        msg.append('bdf.mpcadds[%s]' % mpc_id)
        groups_dict = {}
        for mpcadd in mpcadds:
            groups_dict[mpcadd.type] = groups_dict.get(mpcadd.type, 0) + 1
        for name, count_name in sorted(groups_dict.items()):
            msg.append('  %-8s %s' % (name + ':', count_name))
        msg.append('')

    for (mpc_id, mpcs) in sorted(model.mpcs.items()):
        msg.append('bdf.mpcs[%s]' % mpc_id)
        groups_dict = {}
        for mpc in mpcs:
            groups_dict[mpc.type] = groups_dict.get(mpc.type, 0) + 1
        for name, count_name in sorted(groups_dict.items()):
            msg.append('  %-8s %s' % (name + ':', count_name))
        msg.append('')

def _cyclic_stats(model: BDF, msg: List[str]) -> None:
    """helper for ``get_bdf_stats(...)``"""
    if model.cyax:
        msg.append('bdf.cyax')
        msg.append('  %-8s 1' % ('CYAX:'))
    if model.cyjoin:
        msg.append('bdf:cyjoin')
        msg.append('  %-8s %s' % ('CYJOIN:', len(model.cyjoin)))

def _aero_stats(model: BDF, msg: List[str]) -> None:
    """helper for ``get_bdf_stats(...)``"""
    if model.aero:
        msg.append('bdf.aero')
        msg.append('  %-8s 1' % ('AERO:'))

    # aero
    if model.aero:
        msg.append('bdf.aero')
        msg.append('  %-8s 1' % ('AERO:'))

    # aeros
    if model.aeros:
        msg.append('bdf:aeros')
        msg.append('  %-8s 1' % ('AEROS:'))

    #mkaeros
    if model.mkaeros:
        msg.append('bdf:mkaeros')
        msg.append('  %-8s %s' % ('MKAERO:', len(model.mkaeros)))

def _nsm_stats(model: BDF, msg: List[str]) -> None:
    """helper for ``get_bdf_stats(...)``"""
    # nsms
    for (nsm_id, nsmadds) in sorted(model.nsmadds.items()):
        msg.append('bdf.nsmadds[%s]' % nsm_id)
        groups_dict = {}
        for nsmadd in nsmadds:
            groups_dict[nsmadd.type] = groups_dict.get(nsmadd.type, 0) + 1
        for name, count_name in sorted(groups_dict.items()):
            msg.append('  %-8s %s' % (name + ':', count_name))
        msg.append('')

    for (mpc_id, nsms) in sorted(model.nsms.items()):
        msg.append('bdf.nsms[%s]' % mpc_id)
        groups_dict = {}
        for nsm in nsms:
            groups_dict[nsm.type] = groups_dict.get(nsm.type, 0) + 1
        for name, count_name in sorted(groups_dict.items()):
            msg.append('  %-8s %s' % (name + ':', count_name))
        msg.append('')

def _get_bdf_stats_loads(model: BDF) -> List[str]:
    """helper for ``get_bdf_stats(...)``"""
    # loads
    msg = []
    if model.is_bdf_vectorized:
        ## kind of hackish
        for (lid, load_combination) in sorted(model.load_combinations.items()):
            msg.append('bdf.load_combinations[%s]' % lid)
            msg.append('')
            if len(model.loads):
                msg.append('bdf.loads[%s] : ???')

    else:
        for (lid, load_combinations) in sorted(model.load_combinations.items()):
            msg.append('bdf.load_combinations[%s]' % lid)
            groups_dict = {}  # type: Dict[str, int]
            for load_combination in load_combinations:
                groups_dict[load_combination.type] = groups_dict.get(load_combination.type, 0) + 1
            for name, count_name in sorted(groups_dict.items()):
                msg.append('  %-8s %s' % (name + ':', count_name))
            msg.append('')

        for (lid, loads) in sorted(model.loads.items()):
            msg.append('bdf.loads[%s]' % lid)
            groups_dict = {}
            for loadi in loads:
                groups_dict[loadi.type] = groups_dict.get(loadi.type, 0) + 1
            for name, count_name in sorted(groups_dict.items()):
                msg.append('  %-8s %s' % (name + ':', count_name))
            msg.append('')
    return msg
