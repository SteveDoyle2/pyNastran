"""
defines:
 - bdf_merge(bdf_filenames, bdf_filename_out=None, renumber=True, encoding=None, size=8,
             is_double=False, cards_to_skip=None, log=None, skip_case_control_deck=False)

"""
from io import StringIO
from pyNastran.bdf.bdf import BDF, read_bdf
from pyNastran.bdf.case_control_deck import CaseControlDeck
from pyNastran.bdf.mesh_utils.bdf_renumber import bdf_renumber, get_renumber_starting_ids_from_model

def bdf_merge(bdf_filenames, bdf_filename_out=None, renumber=True, encoding=None, size=8,
              is_double=False, cards_to_skip=None, log=None, skip_case_control_deck=False):
    """
    Merges multiple BDF into one file

    Parameters
    ----------
    bdf_filenames : List[str]
        list of bdf filenames
    bdf_filename_out : str / None
        the output bdf filename (default=None; None -> no writing)
    renumber : bool
        should the bdf be renumbered (default=True)
    encoding : str
        the unicode encoding (default=None; system default)
    size : int; {8, 16}; default=8
        the bdf write precision
    is_double : bool; default=False
        the field precision to write
    cards_to_skip : List[str]; (default=None -> don't skip any cards)
        There are edge cases (e.g. FLUTTER analysis) where things can break due to
        uncross-referenced cards.  You need to disable entire classes of cards in
        that case (e.g. all aero cards).
    skip_case_control_deck : bool, optional, default : False
        If true, don't consider the case control deck while merging.

    Returns
    -------
    model : BDF
        Merged model.
    mappers_all : List[mapper]
        mapper : Dict[bdf_attribute] : old_id_to_new_id_dict
            List of mapper dictionaries of original ids to merged

            bdf_attribute : str
                a BDF attribute (e.g., 'nodes', 'elements')
            old_id_to_new_id_dict : dict[id_old] : id_new
                a sub dictionary that is used to map the node/element/etc. ids
            mapper = {
                'elements' : eid_map,
                'nodes' : nid_map,
                'coords' : cid_map,
                ...
            }

    Supports
      nodes:      GRID
      coords:     CORDx
      elements:   CQUAD4, CTRIA3, CTETRA, CPENTA, CHEXA, CELASx, CBAR, CBEAM
                  CONM1, CONM2, CMASS
      properties: PSHELL, PCOMP, PSOLID, PMASS
      materials:  MAT1, MAT8

    .. todo:: doesn't support SPOINTs/EPOINTs
    .. warning:: still very preliminary

    """
    if not isinstance(bdf_filenames, (list, tuple)):
        raise TypeError('bdf_filenames is not a list/tuple...%s' % str(bdf_filenames))

    if not len(bdf_filenames) > 1:
        raise RuntimeError("You can't merge one BDF...bdf_filenames=%s" % str(bdf_filenames))
    for bdf_filename in bdf_filenames:
        if not isinstance(bdf_filename, (str, BDF, StringIO)):
            raise TypeError('bdf_filenames is not a string/BDF...%s' % bdf_filename)

        #bdf_filenames = [bdf_filenames]

    #starting_id_dict_default = {
        #'cid' : max(model.coords.keys()),
        #'nid' : max(model.nodes.keys()),
        #'eid' : max([
            #max(model.elements.keys()),
            #max(model.masses.keys()),
        #]),
        #'pid' : max([
            #max(model.properties.keys()),
            #max(model.properties_mass.keys()),
        #]),
        #'mid' : max(model.material_ids),
    #}
    bdf_filename0 = bdf_filenames[0]
    if isinstance(bdf_filename0, BDF):
        model = bdf_filename0
    else:
        model = BDF(debug=False, log=log)
        model.disable_cards(cards_to_skip)
        model.read_bdf(bdf_filename0, encoding=encoding, validate=False)

    if skip_case_control_deck:
        model.case_control_deck = CaseControlDeck([], log=None)
    model.log.info('primary=%s' % bdf_filename0)

    _mapper_0 = _get_mapper_0(model) # mapper for first model

    data_members = [
        'coords', 'nodes', 'elements', 'masses', 'properties', 'properties_mass',
        'materials', 'sets', 'rigid_elements', 'mpcs', 'caeros', 'splines',
    ]
    mappers = []
    for bdf_filename in bdf_filenames[1:]:
        starting_id_dict = get_renumber_starting_ids_from_model(model)
        #for param, val in sorted(starting_id_dict.items()):
            #print('  %-3s %s' % (param, val))

        model.log.info('secondary=%s' % bdf_filename)
        if isinstance(bdf_filename, BDF):
            model2_renumber = bdf_filename
        else:
            model2_renumber = BDF(debug=False, log=log)
            model2_renumber.disable_cards(cards_to_skip)
            model2_renumber.read_bdf(bdf_filename)

        _apply_scalar_cards(model, model2_renumber)

        bdf_dump = StringIO() # 'bdf_merge_temp.bdf'
        _, mapperi = bdf_renumber(model2_renumber, bdf_dump, starting_id_dict=starting_id_dict,
                                  size=size, is_double=is_double, cards_to_skip=cards_to_skip)
        bdf_dump.seek(0)

        mappers.append(mapperi)
        model2 = BDF(debug=False, log=log)
        model2.disable_cards(cards_to_skip)
        model2.read_bdf(bdf_dump)

        #model.log.info('model2.node_ids = %s' % np.array(model2.node_ids))
        for data_member in data_members:
            data1 = getattr(model, data_member)
            data2 = getattr(model2, data_member)
            if isinstance(data1, dict):
                #model.log.info('  working on %s' % (data_member))
                for key, value in data2.items():
                    if data_member in 'coords' and key == 0:
                        continue

                    if isinstance(value, list):
                        assert key not in data1, key
                        data1[key] = value
                    else:
                        assert key not in data1, f'{data_member} key={key}\n{data1}'
                        data1[key] = value
                        #print('   %s' % key)
            else:  # pragma: no cover
                raise NotImplementedError(type(data1))
    #if bdf_filenames_out:
        #model.write_bdf(bdf_filenames_out, size=size)

    mapper_renumber = None
    if renumber:
        model.log.info('final renumber...')

        starting_id_dict = {
            'cid' : 1,
            'nid' : 1,
            'eid' : 1,
            'pid' : 1,
            'mid' : 1,
        }
        _, mapper_renumber = bdf_renumber(model, bdf_filename_out,
                                          starting_id_dict=starting_id_dict, size=size,
                                          is_double=is_double, cards_to_skip=cards_to_skip)
        bdf_filename_temp = StringIO()
        model.write_bdf(bdf_filename_temp, size=size, is_double=False, interspersed=False,
                        enddata=None, close=False)
        bdf_filename_temp.seek(0)
        model = read_bdf(bdf_filename_temp, validate=False, xref=model._xref, punch=False,
                         log=model.log, debug=True, mode=model._nastran_format)

    elif bdf_filename_out:
        model.write_bdf(out_filename=bdf_filename_out, encoding=None,
                        size=size, is_double=is_double,
                        interspersed=True,
                        enddata=None)

    mappers_final = _assemble_mapper(mappers, _mapper_0, data_members,
                                     mapper_renumber=mapper_renumber)
    return model, mappers_final

def _apply_scalar_cards(model, model2_renumber):
    """apply cards from model2 to model if they don't exist in model"""
    if model.aero is None and model2_renumber.aero:
        model.aero = model2_renumber.aero
    if model.aeros is None and model2_renumber.aeros:
        model.aeros = model2_renumber.aeros
    model.mkaeros += model2_renumber.mkaeros
    for key, param in model2_renumber.params.items():
        if key not in model.params:
            model.params[key] = param

def _assemble_mapper(mappers, mapper_0, data_members, mapper_renumber=None):
    """
    Assemble final mappings from all original ids to the ids in the merged and possibly
    renumbered model.

    Parameters
    ----------
    mappers : List[mapper]
       mapper : dict[key] : value
            key : ???
                ???
            value : ???
                ???
    mapper_0 : mapper
        key : ???
            ???
        value : ???
            ???
    data_members : List[str]
       list of things to include in the mappers?
        data_members = [
            'coords', 'nodes', 'elements', 'masses', 'properties', 'properties_mass',
            'materials', 'sets', 'rigid_elements', 'mpcs',
        ]
    mapper_renumber : dict[key]: value; default=None
        key : str
            a BDF attribute
        value : dict[id_old] : id_new
            a sub dictionary that is used to map the node/element/etc. ids
        mapper = {
            'elements' : eid_map,
            'nodes' : nid_map,
            'coords' : cid_map,
            ...
        }

    Returns
    -------
    mappers_all : List[mappers]
        One mapper for each bdf_filename

    """
    if mapper_renumber is not None:
        mappers_all = [_renumber_mapper(mapper_0, mapper_renumber)]

        for mapper in mappers:
            mapper_temp = {}
            for map_type in data_members:
            #for map_type, sub_mappper in mapper.items():
                sub_mappper = mapper[map_type]
                mapper_temp[map_type] = {}
                for id_orig, id_merge in sub_mappper.items():
                    # map from original to renumbered
                    mapper_temp[map_type][id_orig] = mapper_renumber[map_type][id_merge]
            mappers_all.append(mapper_temp)
    else:
        # the first model nids are unchanged
        mappers_all = [mapper_0] + mappers

    return mappers_all

def _get_mapper_0(model):
    """
    Get the mapper for the first model.

    Parameters
    ----------
    model : BDF()
        the bdf model object

    Returns
    -------
    mapper : dict[key]: value
        key : str
            a BDF attribute
        value : dict[id_old] : id_new
            a sub dictionary that is used to map the node/element/etc. ids
        mapper = {
            'elements' : eid_map,
            'nodes' : nid_map,
            'coords' : cid_map,
            ...
        }

    """
    # build the maps
    eids_all = (
        list(model.elements.keys()) +
        list(model.masses.keys()) +
        list(model.rigid_elements.keys())
    )
    eid_map = {eid : eid for eid in eids_all}
    nid_map = {nid : nid for nid in model.point_ids}
    cid_map = {cid : cid for cid in model.coord_ids}
    mid_map = {mid : mid for mid in model.material_ids}
    spc_map = _dicts_key_to_key((model.spcs, model.spcadds))
    mpc_map = _dicts_key_to_key((model.mpcs, model.mpcadds))
    method_map = _dict_key_to_key(model.methods)
    properties_map = _dict_key_to_key(model.properties)
    rigid_elements_map = _dict_key_to_key(model.rigid_elements)
    cmethod_map = _dict_key_to_key(model.cMethods)
    flfact_map = _dict_key_to_key(model.flfacts)
    flutter_map = _dict_key_to_key(model.flutters)
    caero_map = _dict_key_to_key(model.caeros)
    freq_map = _dict_key_to_key(model.frequencies)

    dload_map = _dicts_key_to_key((model.dload_entries, model.dloads))
    load_map = _dicts_key_to_key((model.loads, model.load_combinations))
    lseq_map = load_map # wrong???
    temp_map = load_map # wrong???

    tstep_map = _dict_key_to_key(model.tsteps)
    tstepnl_map = _dict_key_to_key(model.tstepnls)
    suport1_map = _dict_key_to_key(model.suport1)
    #suport_map = {}

    nlparm_map = _dict_key_to_key(model.nlparms)
    #nlpci_map = _dict_key_to_key(model.nlpcis)
    table_sdamping_map = _dict_key_to_key(model.tables_sdamping)
    dconadd_map = _dict_key_to_key(model.dconadds)
    dconstr_map = _dict_key_to_key(model.dconstrs)
    dessub_map = dconadd_map
    for key, value in dconstr_map.items():
        if key in dessub_map:
            raise NotImplementedError()
        dessub_map[key] = value
    dresp_map = _dict_key_to_key(model.dresps)
    gust_map = _dict_key_to_key(model.gusts)
    trim_map = _dict_key_to_key(model.trims)
    tic_map = _dict_key_to_key(model.tics)
    csschd_map = _dict_key_to_key(model.csschds)
    tranfer_function_map = _dict_key_to_key(model.transfer_functions)

    mapper = {
        'elements' : eid_map,
        'nodes' : nid_map,
        'coords' : cid_map,
        'materials' : mid_map,
        'properties' : properties_map,
        'rigid_elements': rigid_elements_map,
        'spcs' : spc_map,
        'mpcs' : mpc_map,
        'METHOD' : method_map,
        'CMETHOD' : cmethod_map,
        'FLFACT' : flfact_map,
        'FMETHOD' : flutter_map,
        'caeros' : caero_map,
        'FREQUENCY' : freq_map,

        'DLOAD' : dload_map,
        'LOAD' : load_map,
        'LOADSET' : lseq_map,
        'TSTEP' : tstep_map,
        'TSTEPNL' : tstepnl_map,
        'SUPORT1' : suport1_map,
        'NLPARM' : nlparm_map,
        'SDAMPING' : table_sdamping_map,
        'DESSUB' : dessub_map,
        'DESOBJ' : dresp_map,
        'GUST' : gust_map,
        'TRIM' : trim_map,
        'IC' : tic_map,
        'CSSCHD' : csschd_map,
        'TFL' : tranfer_function_map,
        #'DESSUB' : dessub_map,
        # bad...
        'TEMPERATURE(LOAD)' : temp_map,
        'TEMPERATURE(INITIAL)' : temp_map,
        #'DATAREC' : datarec_map,
        #'ADAPT' : adapt_map,
        #'SUPER' : super_map,
        #'BOUTPUT' : boutput_map,
        #'OUTRCV' : outrcv_map,
    }

    return mapper

def _renumber_mapper(mapper_0, mapper_renumber):
    """
    Renumbers a mapper

    Parameters
    ----------
    mapper_0 : dict[key]: value
        key : str
            a BDF attribute
        value : dict[id_old] : id_new
            a sub dictionary that is used to map the node/element/etc. ids
        mapper = {
            'elements' : eid_map,
            'nodes' : nid_map,
            'coords' : cid_map,
            ...
        }
    mapper_renumber : ???
        ???

    Returns
    -------
    mapper : dict[map_type] : sub_mapper
        map_type : ???
            ???
        sub_mapper : dict[key] : value
            key : ???
                ???
            value : ???
                ???

    """
    mapper = mapper_0.copy()
    # apply any renumbering
    for map_type, sub_mapper in mapper.items():
        for id_ in sub_mapper.keys():
            if sub_mapper[id_] == mapper_renumber[map_type][id_]:
                continue
            sub_mapper[id_] = mapper_renumber[map_type][id_]
    return mapper

def _dict_key_to_key(dictionary):
    """creates a dummy map from the nominal key to the nominal key"""
    return {key : key for key in dictionary.keys()}

def _dicts_key_to_key(dictionaries):
    """
    creates a dummy map from the nominal key to the nominal key for
    multiple input dictionaries

    """
    out = {}
    for dicti in dictionaries:
        for key in dicti:
            out[key] = key
    return out
