def get_bdf_stats(model, return_type='string', word=''):
    # type: (str) -> Union[str, List[str]]
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
        'params', 'nodes', 'points', 'elements', 'normals', 'rigid_elements',
        'properties', 'materials', 'creep_materials',
        'MATT1', 'MATT2', 'MATT3', 'MATT4', 'MATT5', 'MATT8', 'MATT9',
        'MATS1', 'MATS3', 'MATT8',
        'coords', 'mpcs',

        # dynamic cards
        'dareas', 'dphases', 'nlparms', 'nlpcis', 'tsteps', 'tstepnls',
        'rotors',

        # direct matrix input - DMIG - dict
        'dmis', 'dmigs', 'dmijs', 'dmijis', 'dmiks',
        'dequations',

        # frequencies - dict[List[FREQ]]
        'frequencies',

        # optimization - dict
        'dconadds', 'dconstrs', 'desvars', 'ddvals', 'dlinks', 'dresps',
        'dvcrels', 'dvmrels', 'dvprels', 'dvgrids',

        # SESETx - dict
        'suport1',
        'se_sets',
        'se_usets',

        # tables
        'tables', 'tables_d', 'tables_m', 'random_tables',

        # methods
        'methods', 'cMethods',

        # aero
        'caeros', 'paeros', 'aecomps', 'aefacts', 'aelinks',
        'aelists', 'aeparams', 'aesurfs', 'aestats', 'gusts', 'flfacts',
        'flutters', 'splines', 'trims',

        # thermal
        'bcs', 'thermal_materials', 'phbdys', 'views', 'view3ds',
        'convection_properties', ]

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
    msg.append('SOL %s\n' % model.sol)
    msg.extend(_get_bdf_stats_loads(model))

    # dloads
    for (lid, loads) in sorted(model.dloads.items()):
        msg.append('bdf.dloads[%s]' % lid)
        groups_dict = {}
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

    # aero
    if model.aero:
        msg.append('bdf:aero')
        msg.append('  %-8s 1' % ('AERO:'))

    # aeros
    if model.aeros:
        msg.append('bdf:aeros')
        msg.append('  %-8s 1' % ('AEROS:'))

    #mkaeros
    if model.mkaeros:
        msg.append('bdf:mkaeros')
        msg.append('  %-8s %s' % ('MKAERO:', len(model.mkaeros)))

    # radset
    if model.radset:
        msg.append('bdf:radset')
        msg.append('  %-8s 1' % ('RADSET:'))

    #mkaeros
    if model.seqgp:
        msg.append('bdf:seqgp')
        msg.append('  %-8s 1' % ('SEQGP:'))

    for card_group_name in card_dict_groups:
        try:
            card_group = getattr(model, card_group_name)
        except AttributeError:
            msgi = 'cant find card_group_name=%r' % card_group_name
            raise AttributeError(msgi)

        groups = set([]) # type: Set[str]

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
                if card_name == 'CORD2R':
                    continue
                group_msg.append('  %-8s : ???' % card_name)
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

def _get_bdf_stats_loads(model):
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
            groups_dict = {}  # type: Dict[str, int]
            for loadi in loads:
                groups_dict[loadi.type] = groups_dict.get(loadi.type, 0) + 1
            for name, count_name in sorted(groups_dict.items()):
                msg.append('  %-8s %s' % (name + ':', count_name))
            msg.append('')
    return msg
