"""
defines:
    bdf_renumber(bdf_filename, bdf_filename_out, size=8, is_double=False,
                 starting_id_dict=None, round_ids=False, cards_to_skip=None,
                 log=None, debug=False)
    superelement_renumber(bdf_filename, bdf_filename_out=None, size=8, is_double=False,
                          starting_id_dict=None, cards_to_skip=None,
                          log=None, debug=False)

"""
from itertools import chain
from io import StringIO, IOBase
from typing import List, Dict, Optional, Union

import numpy as np

from pyNastran.bdf.bdf import BDF
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.utils.mathematics import roundup


def bdf_renumber(bdf_filename: Union[str, BDF, StringIO],
                 bdf_filename_out: str,
                 size=8, is_double=False,
                 starting_id_dict=None, round_ids: bool=False,
                 cards_to_skip: Optional[List[str]]=None,
                 log=None, debug=False) -> BDF:
    """
    Renumbers a BDF

    Parameters
    ----------
    bdf_filename : str / BDF
        str : a bdf_filename (string; supported)
        BDF : a BDF model that has been cross referenced and is
        fully valid (an equivalenced deck is not valid)
    bdf_filename_out : str / None
        str : a bdf_filename to write
        None : don't write the BDF
    size : int; {8, 16}; default=8
        the bdf write precision
    is_double : bool; default=False
        the field precision to write
    starting_id_dict : dict, None (default=None)
        None : renumber everything starting from 1
        dict : {key : starting_id}
            key : str
                the key (e.g. eid, nid, cid, ...)
            starting_id : int, None
                int : the value to start from
                None : don't renumber this key
    round_ids : bool; default=False
        Should a rounding up be applied for each variable?
        This makes it easier to read a deck and verify that it's been
        renumbered properly.
        This only really applies when starting_id_dict is None
    cards_to_skip : List[str]; (default=None -> don't skip any cards)
        There are edge cases (e.g. FLUTTER analysis) where things can
        break due to uncross-referenced cards.  You need to disable
        entire classes of cards in that case (e.g. all aero cards).

    Returns
    -------
    model : BDF()
        a renumbered BDF object corresponding to bdf_filename_out
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


    .. todo:: bdf_model option for bdf_filename hasn't been tested
    .. todo:: add support for subsets (e.g. renumber only a subset of nodes/elements)
    .. todo:: doesn't support partial renumbering
    .. todo:: doesn't support element material coordinate systems

    ..warning :: spoints might be problematic...check
    ..warning :: still in development, but it usually brutally crashes
                 if it's not supported
    ..warning :: be careful of card unsupported cards (e.g. ones not read in)

    Supports
     - GRIDs
       - no superelements
     - COORDx
     - elements
        - CELASx/CONROD/CBAR/CBEAM/CQUAD4/CTRIA3/CTETRA/CPENTA/CHEXA
        - RBAR/RBAR1/RBE1/RBE2/RBE3/RSPLINE/RSSCON
     - properties
        - PSHELL/PCOMP/PCOMPG/PSOLID/PSHEAR/PBAR/PBARL
          PROD/PTUBE/PBEAM
     - mass
        - CMASSx/CONMx/PMASS
     - aero
       - FLFACT
       - SPLINEx
       - FLUTTER
     - partial case control
       - METHOD/CMETHOD/FREQENCY
       - LOAD/DLOAD/LSEQ/LOADSET...LOADSET/LSEQ is iffy
       - SET cards
         - nodes
         - elements
       - SPC/MPC/FLUTTER/FLFACT
    - constraints
       - SPC/SPCADD/SPCAX/SPCD
       - MPC/MPCADD
       - SUPORT/SUPORT1
    - solution control/methods
       - TSTEP/TSTEPNL
       - NLPARM
       - EIGB/EIGC/EIGRL/EIGR
    - sets
       - USET
    - other
      - tables
      - materials
      - loads/dloads

    Not Done

     - SPOINT
     - any cards with SPOINTs
       - DMIG/DMI/DMIJ/DMIJI/DMIK/etc.
       - CELASx
       - CDAMPx
     - superelements
     - aero cards
       - CAEROx
       - PAEROx
     - thermal cards?
     - optimization cards
     - SETx
     - PARAM,GRDPNT,x; where x>0
     - GRID SEID
     - case control
       - STATSUB
       - SUBCASE
       - global SET cards won't be renumbered properly

    Examples
    --------
    **Renumber Everything; Start from 1**

    >>> bdf_renumber(bdf_filename, bdf_filename_out, size=8, is_double=False,
                     round_ids=False)

    **Renumber Everything; Start Material IDs from 100**

    >>> starting_id_dict = {
        'mid' : 100,
    }
    >>> bdf_renumber(bdf_filename, bdf_filename_out, size=8, is_double=False,
                     starting_ids_dict=starting_ids_dict, round_ids=False)

    **Only Renumber Material IDs**

    >>> starting_id_dict = {
        'cid' : None,
        'nid' : None,
        'eid' : None,
        'pid' : None,
        'mid' : 1,
        'spc_id' : None,
        'mpc_id' : None,
        'load_id' : None,
        'dload_id' : None,

        'method_id' : None,
        'cmethod_id' : None,
        'spline_id' : None,
        'table_id' : None,
        'flfact_id' : None,
        'flutter_id' : None,
        'freq_id' : None,
        'tstep_id' : None,
        'tstepnl_id' : None,
        'suport_id' : None,
        'suport1_id' : None,
        'tf_id' : None,
        'set_id' : None,
    }
    >>> bdf_renumber(bdf_filename, bdf_filename_out, size=8, is_double=False,
                     starting_ids_dict=starting_ids_dict, round_ids=False)

    """
    assert size in [8, 16], size
    assert isinstance(is_double, bool), is_double
    starting_id_dict_default = {
        'cid' : 1,
        'nid' : 1,
        'eid' : 1,
        'pid' : 1,
        'mid' : 1,
        'spc_id' : 1,
        'mpc_id' : 1,
        'load_id' : 1,
        'set_id' : 1,
        'dload_id' : 1,

        'method_id' : 1,
        'cmethod_id' : 1,
        'spline_id' : 1,
        'caero_id' : 1,
        'table_id' : 1,
        'flfact_id' : 1,
        'flutter_id' : 1,
        'freq_id' : 1,
        'tstep_id' : 1,
        'tstepnl_id' : 1,
        'suport_id' : 1,
        'suport1_id' : 1,
        'tf_id' : 1,
    }
    # fill up starting_id_dict
    if starting_id_dict is None:
        starting_id_dict = starting_id_dict_default
    else:
        for key, value in starting_id_dict_default.items():
            if key not in starting_id_dict:
                starting_id_dict[key] = value

    nid = None
    cid = None
    eid = None
    pid = None
    mid = None
    nsm_id = None
    spc_id = None
    mpc_id = None
    load_id = None
    dload_id = None
    method_id = None
    cmethod_id = None
    spline_id = None
    table_id = None
    flfact_id = None
    flutter_id = None
    freq_id = None
    ##tstep_id = None
    tstepnl_id = None
    set_id = None
    suport_id = None
    suport1_id = None
    tf_id = None

    # turn them into variables
    for key, value in sorted(starting_id_dict.items()):
        #assert isinstance(key, str), key
        assert key in starting_id_dict_default, 'key=%r is invalid' % (key)
        #assert isidentifier(key), 'key=%s is invalid' % key
        if value is None:
            pass
        else:
            if not isinstance(value, integer_types):
                msg = 'key=%r value=%r must be an integer; type(value)=%s' % (
                    key, value, type(value))
                raise TypeError(msg)

        if key == 'nid':
            nid = int(value)
        elif key == 'cid':
            if value is None:
                cid = None
            else:
                cid = int(value)
        elif key == 'set_id':
            if value is None:
                set_id = None
            else:
                set_id = int(value)
        elif key == 'eid':
            eid = int(value)
        elif key == 'pid':
            if value is None:
                pid = None
            else:
                pid = int(value)
        elif key == 'mid':
            if value is None:
                mid = None
            else:
                mid = int(value)
        elif key == 'spc_id':
            spc_id = int(value)
        elif key == 'mpc_id':
            mpc_id = int(value)

        elif key == 'load_id':
            load_id = int(value)
        elif key == 'dload_id':
            dload_id = int(value)
        elif key == 'method_id':
            method_id = int(value)
        elif key == 'cmethod_id':
            cmethod_id = int(value)
        elif key == 'spline_id':
            spline_id = int(value)
        elif key == 'caero_id':
            caero_id = int(value)
        elif key == 'table_id':
            table_id = int(value)
        elif key == 'flfact_id':
            flfact_id = int(value)
        elif key == 'flutter_id':
            flutter_id = int(value)
        elif key == 'freq_id':
            freq_id = int(value)
        elif key == 'tstep_id':
            tstep_id = int(value)
        elif key == 'tstepnl_id':
            tstepnl_id = int(value)
        elif key == 'set_id':
            set_id = int(value)
        elif key == 'suport_id':
            suport_id = int(value)
        elif key == 'suport1_id':
            suport1_id = int(value)
        elif key == 'tf_id':
            tf_id = int(value)
        else:
            raise NotImplementedError('key=%r' % key)

    # just to make pylint be quiet
    str(
        (nsm_id, method_id, cmethod_id, flfact_id, flutter_id,
         tstep_id, tstepnl_id,
         suport_id, suport1_id)
    )

    # build the maps
    mass_id_map = {}
    properties_map = {}
    properties_mass_map = {}
    eid_map = {}
    rigid_elements_map = {}
    #nsm_map = {}
    cid_map = {}
    mpc_map = {}
    spc_map = {}
    dload_map = {}
    load_map = {}

    cmethod_map = {}
    method_map = {}
    flfact_map = {}
    flutter_map = {}
    freq_map = {}
    tstep_map = {}
    tstepnl_map = {}
    #suport_map = {}
    suport1_map = {}

    model = _get_bdf_model(bdf_filename, cards_to_skip=cards_to_skip, log=log, debug=debug)

    nid_map, unused_reverse_nid_map = _create_nid_maps(model, starting_id_dict, nid)
    mid_map, all_materials = _create_mid_map(model, mid)

    _update_nodes(
        model, starting_id_dict, nid,
        nid_map)

    _update_properties(
        model, starting_id_dict, pid,
        properties_map, properties_mass_map)

    _update_elements(
        model, starting_id_dict, eid,
        eid_map, mass_id_map, rigid_elements_map)

    _update_materials(
        model, starting_id_dict, mid,
        mid_map, all_materials)

    _update_spcs(
        model, starting_id_dict, spc_id,
        spc_map)

    _update_mpcs(
        model, starting_id_dict, mpc_id,
        mpc_map)

    _update_coords(
        model, starting_id_dict, cid,
        cid_map)

    if 'freq_id' in starting_id_dict and freq_id is not None:
        # frequencies
        for freqi, freqs in sorted(model.frequencies.items()):
            freq_map[freqi] = freq_id
            for freq in freqs:
                freq.sid = freqi
            freq_id += 1

    set_map = {}
    if 'set_id' in starting_id_dict and set_id is not None:
        # sets
        for sidi, set_ in sorted(model.sets.items()):
            set_.sid = set_id
            set_map[sidi] = set_id
            set_id += 1

    if 'spline_id' in starting_id_dict and spline_id is not None:
        # set up spline1 box mapping
        delta_box1_map = {}
        delta_box2_map = {}
        for sidi, spline in sorted(model.splines.items()):
            if spline.type in ['SPLINE1', 'SPLINE2']:
                delta_box1_map[sidi] = spline.box1 - spline.caero
                delta_box2_map[sidi] = spline.box2 - spline.caero
            else:
                # should be handled by the xref?
                pass
            #else:
                #raise NotImplementedError(spline)

    caero_id_map = {}
    if 'caero_id' in starting_id_dict and caero_id is not None:
        # caeros
        for caero_idi, caero in sorted(model.caeros.items()):
            if caero.type in ['CAERO1', 'CAERO3', 'CAERO4']: # not CAERO5
                caero.eid = caero_id
                caero_id_map[caero_idi] = caero_id
                caero_id += caero.shape[0] * caero.shape[1]
            elif caero.type == 'CAERO2':
                caero.eid = caero_id
                caero_id_map[caero_idi] = caero_id
                caero_id += caero.nboxes
            else:
                raise NotImplementedError(caero)

    spline_id_map = {}
    if 'spline_id' in starting_id_dict and spline_id is not None:
        # splines
        for sidi, spline in sorted(model.splines.items()):
            spline.eid = spline_id
            #spline.cross_reference(model)
            if spline.type in ['SPLINE1', 'SPLINE2']:
                spline.box1 = caero_id_map[spline.caero] + delta_box1_map[sidi]
                spline.box2 = caero_id_map[spline.caero] + delta_box2_map[sidi]
            else:
                # should be handled by the xref?
                pass
            #else:
                #raise NotImplementedError(spline)
            spline_id_map[sidi] = spline_id
            spline_id += 1

    nlparm_map = {}
    nlpci_map = {}
    table_sdamping_map = {}
    dconstr_map = {}
    dconadd_map = {}
    dresp_map = {}
    gust_map = {}
    trim_map = {}
    tic_map = {}
    csschd_map = {}
    tranfer_function_map = {}
    data = (
        (model.methods, 'sid', method_map),
        (model.cMethods, 'sid', cmethod_map),
        (model.flfacts, 'sid', flfact_map),
        (model.flutters, 'sid', flutter_map),
        (model.tsteps, 'sid', tstep_map),
        (model.tstepnls, 'sid', tstepnl_map),
        (model.suport1, 'conid', suport1_map),
        (model.nlparms, 'nlparm_id', nlparm_map),
        (model.nlpcis, 'nlpci_id', nlpci_map),
        (model.tables_sdamping, 'tid', table_sdamping_map),
        (model.dconadds, 'dcid', dconadd_map),
        #(model.dconstrs, 'oid', dconstr_map),
        (model.dresps, 'dresp_id', dresp_map),
        (model.gusts, 'sid', gust_map),
        (model.trims, 'sid', trim_map),
        (model.tics, 'sid', tic_map),
        (model.csschds, 'sid', csschd_map),
        (model.aefacts, 'sid', None),
        (model.aelinks, 'sid', None),
        (model.aelists, 'sid', None),
        (model.paeros, 'pid', None),

        #(model.asets, 'sid', None),
        (model.dareas, 'sid', None),
        (model.transfer_functions, 'sid', tranfer_function_map),
        #(model.bsets, 'sid', None),
        #(model.csets, 'sid', None),
        #(model.qsets, 'sid', None),
        #(model.usets, 'sid', None),

        #(model.se_sets, 'sid', None),
        #(model.se_asets, 'sid', None),
        #(model.se_bsets, 'sid', None),
        #(model.se_csets, 'sid', None),
        #(model.se_qsets, 'sid', None),
        #(model.se_usets, 'sid', None),
    )

    # apply the simple to update parameters
    param_id = 9999
    for (dict_obj, param_name, mmap) in data:
        if round_ids:
            param_id = roundup(param_id, 1000) + 1
        else:
            param_id = 1
        for idi, param in sorted(dict_obj.items()):
            #print('working on id=%s param=%s' % (str(idi), str(param)))
            try:
                msg = '%s has no %r; use %s' % (param.type, param_name, param.object_attributes())
            except AttributeError:
                model.log.error('param = %r' % param)
                raise
            assert hasattr(param, param_name), msg
            setattr(param, param_name, param_id)
            if mmap is not None:
                mmap[idi] = param_id
            param_id += 1

    # start the complicated set
    # dconstr
    dessub_map = dconadd_map
    for key, value in dconstr_map.items():
        if key in dessub_map:
            raise NotImplementedError()
        dessub_map[key] = value

    # tables
    for unused_table_idi, table in sorted(model.tables.items()):
        assert hasattr(table, 'tid')
        table.tid = table_id
        table_id += 1
    for unused_table_idi, table in sorted(model.random_tables.items()):
        assert hasattr(table, 'tid')
        table.tid = table_id
        table_id += 1

    # dloads
    for dload_idi, dloads in sorted(model.dloads.items()):
        for dload in dloads:
            assert hasattr(dload, 'sid')
            dload.sid = dload_id
        dload_map[dload_idi] = dload_id
        dload_id += 1
    for dload_idi, dloads in sorted(model.dload_entries.items()):
        for dload in dloads:
            assert hasattr(dload, 'sid')
            dload.sid = dload_id
        dload_map[dload_idi] = dload_id
        dload_id += 1

    # loads
    for load_idi, load_combinations in sorted(model.load_combinations.items()):
        for load_combination in load_combinations:
            assert hasattr(load_combination, 'sid')
            load_combination.sid = load_id
        load_map[load_idi] = load_id
        load_id += 1
    for load_idi, loads in sorted(model.loads.items()):
        for load in loads:
            assert hasattr(load, 'sid')
            load.sid = load_id
        load_map[load_idi] = load_id
        load_id += 1

    # transfer_functions
    for tf_idi, transfer_functions in sorted(model.transfer_functions.items()):
        for transfer_function in transfer_functions:
            assert hasattr(transfer_function, 'sid')
            transfer_function.sid = tf_id
        tranfer_function_map[tf_idi] = tf_id
        load_id += 1

    lseq_map = load_map # wrong???
    temp_map = load_map # wrong???

    mapper = {
        'elements' : eid_map,
        'masses' : mass_id_map,
        'rigid_elements' : rigid_elements_map,
        'nodes' : nid_map,
        'coords' : cid_map,
        'materials' : mid_map,
        'properties' : properties_map,
        'properties_mass' : properties_mass_map,
        'spcs' : spc_map,
        'mpcs' : mpc_map,
        'METHOD' : method_map,
        'CMETHOD' : cmethod_map,
        'FLFACT' : flfact_map,
        'FMETHOD' : flutter_map,
        'FREQUENCY' : freq_map,
        'sets' : set_map,
        'splines' : spline_id_map,
        'caeros' : caero_id_map,

        'DLOAD' : dload_map,
        'RANDOM' : dload_map, # RANDPS, RANDT1
        'LOAD' : load_map,
        'LOADSET' : lseq_map,
        'CLOAD' : lseq_map,
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
    #print('****suport1_map', suport1_map)
    #print('****dessub_map', dessub_map)
    #print('****dresp_map', dresp_map)
    _update_case_control(model, mapper)

    _write_bdf(model, bdf_filename_out, size=size, is_double=is_double)
    return model, mapper


def _write_bdf(model, bdf_filename_out, size=8, is_double=False):
    """helper method"""
    if bdf_filename_out is not None:
        close = True
        if isinstance(bdf_filename_out, IOBase):
            close = False
        model.write_bdf(bdf_filename_out, size=size, is_double=is_double,
                        interspersed=False, close=close)


def get_renumber_starting_ids_from_model(model: BDF) -> Dict[str, int]:
    """
    Get the starting ids dictionary used for renumbering with ids greater
    than those in model.

    Parameters
    ----------
    model : BDF
        BDF object to get maximum ids from.

    Returns
    -------
    starting_id_dict : dict {str : int, ...}
        Dictionary from id type to starting id.
    """
    eid_max = max([
        max(model.elements),
        max(model.masses) if model.masses else 0,
        max(model.rigid_elements) if model.rigid_elements else 0,
    ])
    pid_max = max([
        max(model.properties) if model.properties else 0,
        max(model.properties_mass) if model.properties_mass else 0,
    ])
    mids = model.material_ids
    starting_id_dict = {
        'cid' : max(model.coords) + 1,
        'nid' : max(model.point_ids) + 1,
        'eid' : eid_max + 1,
        'pid' : pid_max + 1,
        'mid' : max(mids) + 1 if mids else 1,
        'set_id' : max(model.sets) + 1 if model.sets else 1,
        'spline_id' : max(model.splines) + 1 if model.splines else 1,
        'caero_id' : max(caero.box_ids[-1, -1]
                         for caero in model.caeros.values()) + 1 if model.caeros else 1,
    }
    return starting_id_dict


def get_starting_ids_dict_from_mapper(model, mapper):
    starting_id_dict2 = {}
    missed_keys = []
    name_map = {
        'nodes' : 'nid',
        'elements' : 'eid',
        'properties' : 'pid',
        'materials' : 'mid',
        'coords' : 'cid',
        'TFL' : 'tf_id',
        'FREQUENCY' : 'freq_id',
        'splines' : 'spline_id',
        'METHOD' : 'method_id',
        'CMETHOD' : 'cmethod_id',
        'caeros' : 'caero_id',
        'TSTEP' : 'tstep_id',
        'TSTEPNL' : 'tstepnl_id',
        'FLFACT' : 'flfact_id',
        'FMETHOD' : 'flutter_id',
        'spcs' : 'spc_id',
        'mpcs' : 'mpc_id',
        'LOAD' : 'load_id',
        'DLOAD' : 'dload_id', 'RANDOM' : 'dload_id',
        'SUPORT' : 'suport_id',
        'SUPORT1' : 'suport1_id',
        'sets' : 'set_id',

        # other valid names in starting_id_dict
        #'table_id' : 1,
    }
    #print('------------------------------------')
    for key, old_new_map in mapper.items():
        if key in name_map:
            key2 = name_map[key]
            if not old_new_map:
                continue
            #print("%s map = %s" % (key, old_new_map))
            #nold_mapper = 0
            #if key2 in old_mapper:
                #nold_mapper = len(old_mapper[key2].keys())

            #offset = 1
            #if key2 in starting_id_dict:
                #offset = max(old_new_map.keys())
                #offset = starting_id_dict[key2] + nold_mapper
                #offset = 1

            #if key2 == 'nid':
                #print("  nid_offset = %s" % offset)
                #print("  nid map = %s" % old_new_map)
            starting_id_dict2[key2] = max(old_new_map.values()) + 1 # offset
        else:
            if len(old_new_map):
                missed_keys.append(key)
    if missed_keys:
        model.log.warning('bdf_renumber: missed_keys = %s' % missed_keys)
    return starting_id_dict2


def superelement_renumber(bdf_filename, bdf_filename_out=None, size=8, is_double=False,
                          starting_id_dict=None, cards_to_skip=None,
                          log=None, debug=False):
    """
    Renumbers a superelement

    Parameters
    ----------
    bdf_filename : str / BDF
        str : a bdf_filename (string; supported)
        BDF : a BDF model that has been cross referenced and is
        fully valid (an equivalenced deck is not valid)
    bdf_filename_out : str / None
        str : a bdf_filename to write
        None : don't write the BDF
    size : int; {8, 16}; default=8
        the bdf write precision
    is_double : bool; default=False
        the field precision to write
    starting_id_dict : dict, None (default=None)
        None : renumber everything starting from 1
        dict : {key : starting_id}
            key : str
                the key (e.g. eid, nid, cid, ...)
            starting_id : int, None
                int : the value to start from
                None : don't renumber this key
    cards_to_skip : List[str]; (default=None -> don't skip any cards)
        There are edge cases (e.g. FLUTTER analysis) where things can
        break due to uncross-referenced cards.  You need to disable
        entire classes of cards in that case (e.g. all aero cards).

    Returns
    -------
    model : BDF()
        a renumbered BDF object corresponding to bdf_filename_out
    """
    if starting_id_dict is None:
        starting_id_dict = {
            'cid' : 1,
            'nid' : 1,
            'eid' : 1,
            'pid' : 1,
            'mid' : 1,
        }

    model = _get_bdf_model(bdf_filename, cards_to_skip=cards_to_skip, log=log, debug=debug)

    _bdf_filename_out = None
    _model, mapper = bdf_renumber(
        model, _bdf_filename_out,
        size=8, is_double=False, starting_id_dict=starting_id_dict, round_ids=False,
        cards_to_skip=None, log=None, debug=False)

    starting_id_dict_new = get_starting_ids_dict_from_mapper(model, mapper)
    #mapper_short = {key : value for key, value in mapper.items() if len(value)}
    #print('mapper_short =', mapper_short)

    for unused_super_id, superelement in sorted(model.superelement_models.items()):
        _smodel, superelement_mapper = bdf_renumber(
            superelement, _bdf_filename_out,
            size=8, is_double=False, starting_id_dict=starting_id_dict_new, round_ids=False,
            cards_to_skip=None, log=None, debug=False)

        #mapper2 = {key : value for key, value in superelement_mapper.items() if len(value)}
        starting_id_dict_new = get_starting_ids_dict_from_mapper(
            model, superelement_mapper)

    if bdf_filename_out is not None:
        model.write_bdf(bdf_filename_out)

    _write_bdf(model, bdf_filename_out, size=size, is_double=is_double)
    return model #, mapper


def _create_nid_maps(model, starting_id_dict, nid):
    """builds the nid_maps"""
    nid_map = {}
    reverse_nid_map = {}

    spoints = list(model.spoints.keys())
    epoints = list(model.epoints.keys())
    nids = model.nodes.keys()

    nids_spoints_epoints = sorted(chain(nids, spoints, epoints))

    #spoints_nids.sort()
    i = 1
    #nnodes = len(spoints_nids)

    i = 1
    #j = 1
    #print(spoints_nids)
    #k = 0
    #model.log.debug(starting_id_dict)
    if 'nid' in starting_id_dict and nid is not None:
        i = nid
        #banned_nodes = spoints
        for nidi in nids_spoints_epoints:
            if nidi in spoints or nidi in epoints:
                pass
                #print('sid=%s -> %s' % (nid, i))
                #i += 1
            else:
                while i in spoints or i in epoints:
                    #print('*bump')
                    i += 1
                #print('nid=%s -> %s' % (nid, i))
                nid_map[nidi] = i
                reverse_nid_map[i] = nidi
                i += 1
        #for nid in sorted(nids):
            #nid_map[nid] = i
            #reverse_nid_map[i] = nid
            #i += 1
        #print(nid_map)
        #print(reverse_nid_map)
    else:
        for nid in nids_spoints_epoints:
            nid_map[nid] = nid
        reverse_nid_map = nid_map
    return nid_map, reverse_nid_map


def _create_mid_map(model, mid):
    """builds the mid_map"""
    mid_map = {}
    all_materials = (
        model.materials,
        model.creep_materials,
        model.thermal_materials,
        model.hyperelastic_materials,
        model.MATT1,
        model.MATT2,
        model.MATT3,
        model.MATT4,
        model.MATT5,
        #model.MATT6,
        #model.MATT7,
        model.MATT8,
        model.MATT9,
        model.MATS1,
        model.MATS3,
        model.MATS8,
    )

    if mid is not None:
        mids = []
        for materials in all_materials:
            mids += materials.keys()
        mids = np.unique(mids)
        mids.sort()
        nmaterials = len(mids)

        for i in range(nmaterials):
            midi = mids[i]
            mid_map[midi] = mid + i
    return mid_map, all_materials


def _get_bdf_model(bdf_filename, cards_to_skip=None, log=None, debug=False):
    """helper method"""
    if isinstance(bdf_filename, BDF):
        model = bdf_filename
    else:
        model = BDF(log=log, debug=debug)
        model.disable_cards(cards_to_skip)
        model.read_bdf(bdf_filename)
    return model


def _update_nodes(model, starting_id_dict, nid, nid_map):
    """updates the nodes"""
    if 'nid' in starting_id_dict and nid is not None:
        #spoints2 = arange(1, len(spoints) + 1)
        #nid = _create_dict_mapper(model.nodes, nid_map, 'nid', nid)

        for nid, node in sorted(model.nodes.items()):
            nid_new = nid_map[nid]
            #print('nid=%s -> %s' % (nid, nid_new))
            node.nid = nid_new


def _update_properties(model, starting_id_dict, pid,
                       properties_map, properties_mass_map):
    """updates the properties"""
    if 'pid' in starting_id_dict and pid is not None:
        # properties
        #pid = _create_dict_mapper(model.properties, properties_map, 'pid', pid)
        #pid = _create_dict_mapper(model.properties_mass, properties_mass_map, 'pid', pid)
        #pid = _update(model.convection_properties, properties_mass_map, pid)

        for pidi, prop in sorted(model.properties.items()):
            prop.pid = pid
            properties_map[pidi] = pid
            pid += 1
        for pidi, prop in sorted(model.properties_mass.items()):
            # PMASS
            prop.pid = pid
            properties_mass_map[pidi] = pid
            pid += 1
        for pidi, prop in sorted(model.convection_properties.items()):
            # PCONV
            prop.pid = pid
            pid += 1
        for pidi, prop in sorted(model.phbdys.items()):
            # PHBDY
            prop.pid = pid
            pid += 1


def _update_elements(model, starting_id_dict, eid,
                     eid_map, mass_id_map, rigid_elements_map):
    """updates the elements"""
    if 'eid' in starting_id_dict and eid is not None:
        # elements
        #eid = _create_dict_mapper(model.elements, eid_map, 'eid', eid)

        for eidi, element in sorted(model.elements.items()):
            element.eid = eid
            eid_map[eidi] = eid
            eid += 1
        for eidi, element in sorted(model.masses.items()):
            # CONM1, CONM2, CMASSx
            element.eid = eid
            eid_map[eidi] = eid
            mass_id_map[eidi] = eid
            eid += 1
        for eidi, elem in sorted(model.rigid_elements.items()):
            # RBAR/RBAR1/RBE1/RBE2/RBE3/RSPLINE/RSSCON
            elem.eid = eid
            eid_map[eidi] = eid
            rigid_elements_map[eidi] = eid
            eid += 1
        #for eidi, elem in model.caeros.items():
            #pass


def _update_materials(unused_model, starting_id_dict, mid,
                      mid_map, all_materials):
    if 'mid' in starting_id_dict and mid is not None:
        #mid = 1
        for materials in all_materials:
            for midi, material in materials.items():
                mid = mid_map[midi]
                assert hasattr(material, 'mid')
                material.mid = mid


def _update_spcs(model, starting_id_dict, spc_id,
                 spc_map):
    """updates the spcs"""
    if 'spc_id' in starting_id_dict and spc_id is not None:
        # spc
        for spc_idi, spc_group in sorted(model.spcadds.items()):
            for unused_i, spc in enumerate(spc_group):
                assert hasattr(spc, 'conid')
                spc.conid = spc_id
            spc_map[spc_idi] = spc_id
            spc_id += 1

        for spc_idi, spc_group in sorted(model.spcs.items()):
            for unused_i, spc in enumerate(spc_group):
                assert hasattr(spc, 'conid')
                spc.conid = spc_id
            spc_map[spc_idi] = spc_id
            spc_id += 1
    else:
        # TODO: why are we doing this?
        for spc_id in model.spcadds:
            spc_map[spc_id] = spc_id
        for spc_id in model.spcs:
            spc_map[spc_id] = spc_id


def _update_mpcs(model, starting_id_dict, mpc_id,
                 mpc_map):
    """updates the mpcs"""
    if 'mpc_id' in starting_id_dict and mpc_id is not None:
        # mpc
        for mpc_idi, mpc_group in sorted(model.mpcadds.items()):
            for unused_i, mpc in enumerate(mpc_group):
                assert hasattr(mpc, 'conid')
                mpc.conid = mpc_id
            mpc_map[mpc_idi] = mpc_id
            mpc_id += 1

        for mpc_idi, mpc_group in sorted(model.mpcs.items()):
            for unused_i, mpc in enumerate(mpc_group):
                assert hasattr(mpc, 'conid')
                mpc.conid = mpc_id
            mpc_map[mpc_idi] = mpc_id
            mpc_id += 1
    else:
        # TODO: why are we doing this?
        for mpc_id in model.mpcadds:
            mpc_map[mpc_id] = mpc_id
        for mpc_id in model.mpcs:
            mpc_map[mpc_id] = mpc_id


def _update_coords(model, starting_id_dict, cid,
                   cid_map):
    """updates the coords"""
    if 'cid' in starting_id_dict and cid is not None:
        # coords
        for cidi, coord in sorted(model.coords.items()):
            if cidi == 0:
                cid_map[0] = 0
                continue
            coord.cid = cid
            cid_map[cidi] = cid
            cid += 1


def _update_case_control(model, mapper):
    """
    Updates the case control deck; helper method for ``bdf_renumber``.

    Parameters
    ----------
    model : BDF()
        the BDF object
    mapper : dict[str] = List[int]
        Defines the possible case control header values for each entry (e.g. `LOAD`)

    """
    # unified system that uses the same key
    # for bdf_merge and _update_case_control
    case_control_card_to_pynastran_key = {
        'MPC' : 'mpcs',
        'SPC' : 'spcs',
    }

    elemental_quantities = ['STRESS', 'STRAIN', 'FORCE', 'ESE', 'EKE']
    nodal_quantities = [
        'DISPLACEMENT', 'VELOCITY', 'ACCELERATION', 'SPCFORCES', 'MPCFORCES',
        'GPFORCES', 'SDISPLACEMENT', 'OLOAD',
    ]
    mapper_quantities = [
        'FREQUENCY', 'DLOAD', 'LOAD', 'LOADSET', 'SPC', 'MPC', 'METHOD', 'CMETHOD',
        'TSTEP', 'TSTEPNL', 'NLPARM', 'SDAMPING', 'DESSUB', 'DESOBJ', 'GUST',
        'SUPORT1', 'TRIM', 'BOUTPUT', 'IC', 'CSSCHD', 'FMETHOD', 'TFL', 'CLOAD',
    ]

    # TODO: remove this...
    skip_keys_temp = [
        'DESSUB', 'ADAPT', 'DATAREC', 'DSAPRT(END=SENS)=ALL', 'TEMPERATURE(LOAD)',
        'CURVELINESYMBOL', 'DSAPRT=(END=SENS)', 'SUPER', 'BOUTPUT', 'IC',
        'OUTRCV', 'TEMPERATURE(INITIAL)',
    ]

    nid_map = mapper['nodes']
    eid_map = mapper['elements']
    skip_keys = [
        'TITLE', 'ECHO', 'ANALYSIS', 'SUBTITLE', 'LABEL', 'SUBSEQ', 'OUTPUT',
        'TCURVE', 'XTITLE', 'YTITLE', 'AECONFIG', 'AESYMXZ', 'MAXLINES', 'PARAM', 'CONTOUR',
        'PTITLE', 'PLOTTER', 'K2PP', 'CSCALE', 'XGRID LINES', 'YGRID LINES', 'YMIN', 'YMAX',
        'LINE', 'XAXIS', 'YAXIS',
        ] + skip_keys_temp
    warn_keys = ['RANDOM']

    sets_analyzed = set()
    # sets in the global don't get updated....
    # so we're going to find all the sets and
    # map them
    # TODO: renumber the sets
    set_locations = {}
    case_control = model.case_control_deck
    if case_control is None:
        return

    for isubcase, subcase in sorted(case_control.subcases.items()):
        for key, values in sorted(subcase.params.items()):
            value, options, param_type = values
            if 'SET ' in key:
                #print(isubcase, key, value, options, param_type)
                #assert isubcase != 0, isubcase
                if isubcase not in set_locations:
                    set_locations[isubcase] = [key]
                else:
                    set_locations[isubcase].append(key)

    #print('set_locations =', set_locations)
    #iset = 1
    global_subcase = case_control.subcases[0]
    for isubcase, subcase in sorted(case_control.subcases.items()):
        #print('-----------------------')
        for key, values in sorted(subcase.params.items()):
            mapper_key = case_control_card_to_pynastran_key.get(key, key)

            if key in skip_keys:
                pass
            elif key in warn_keys:
                model.log.warning('skipping %s=%s' % (key, values))
            elif 'SET ' in key:
                # does this need to be updated...I don't think so...
                continue
            elif 'SET ' not in key:
                value, options, param_type = values
                if isinstance(value, integer_types):
                    seti = 'SET %i' % value
                    msg = ', which is needed by %s=%s' % (key, value)

                    if key in mapper_quantities:
                        #print('mapper = %s, value=%s' % (key, value))
                        kmap = mapper[mapper_key]
                        try:
                            value2 = kmap[value]
                        except KeyError:
                            msg = 'Could not find id=%s in %s dictionary\n' % (value, key)
                            msg += str(kmap)
                            raise KeyError(msg)
                        subcase.update(
                            key, value2, options, param_type)

                    elif key in elemental_quantities + nodal_quantities:
                        #msg += '  allowed_keys=%s' % sorted(kmap.keys())
                        if seti in subcase:
                            seti2, seti_key = subcase.get_parameter(seti, msg=msg)
                            if seti_key in sets_analyzed:
                                continue
                            sets_analyzed.add(seti_key)
                            msgi = 'seti_key=%s must be an integer; type(seti_key)=%s\n'  % (
                                seti_key, type(seti_key))
                            msgi += '  key=%r value=%r options=%r param_type=%r\n' % (
                                key, value, options, param_type)
                            msgi += '  seti=%r\n' % seti
                            #print(msgi)
                            assert isinstance(seti_key, int), msgi

                            #print('key=%s options=%s param_type=%s value=%s' % (
                                #key, options, param_type, value))
                            #print(seti2)
                        else:
                            seti2 = [value]
                            model.log.error('key=%s options=%s param_type=%s value=%s' % (
                                key, options, param_type, value))
                            raise NotImplementedError(key)

                        eids_missing, nids_missing, values2 = _update_case_key(
                            key, elemental_quantities, seti2, eid_map, nid_map)
                        if eids_missing:
                            model.log.warning(f"  couldn't find eids={eids_missing}...dropping")
                        if nids_missing:
                            model.log.warning(f"  couldn't find nids={nids_missing}...dropping")


                        param_type = 'SET-type'
                        #print('adding seti=%r values2=%r seti_key=%r param_type=%r'  % (
                            #seti, values2, seti_key, param_type))
                        assert len(values2) > 0, 'key=%r values2=%s' % (key, values2)

                        if isubcase in set_locations and key in set_locations[isubcase]:
                            # or not key in global_subcase
                            gset = subcase.get_parameter(seti)
                            lset = subcase.get_parameter(seti)
                            #print('gset', gset)
                            #print('lset', lset)
                            if gset != lset:
                                msg = 'gset=%s lset=%s' % (str(gset), str(lset))
                                raise NotImplementedError(msg)
                                #subcase.update(
                                    #seti, values2, seti_key, param_type)
                            else:
                                global_subcase.update(
                                    seti, values2, seti_key, param_type)
                            #subcase.update(
                                #seti, values2, seti_key, param_type)
                        elif not key not in global_subcase:
                            subcase.update(
                                seti, values2, seti_key, param_type)
                        else:
                            global_subcase.update(
                                seti, values2, seti_key, param_type)
                    else:
                        #print('key=%s seti2=%s' % (key, seti2))
                        model.log.error('key=%r options=%r param_type=%r value=%r' % (
                            key, options, param_type, value))
                        raise RuntimeError(key)
                elif value in ['NONE', 'ALL']:
                    #print('*ALL -> key=%s options=%s param_type=%s value=%s' % (
                        #key, options, param_type, value))
                    #print('*all')
                    pass
                elif key == '':
                    pass
                else:
                    msg = 'key=%s options=%s param_type=%s value=%s' % (
                        key, options, param_type, value)
                    raise RuntimeError(msg)

            else:
                raise RuntimeError(key)
                    #if value ==
        #print()

def _update_case_key(key, elemental_quantities, seti2, eid_map, nid_map):
    """Updates a Case Control SET card.  A set may have an elemental result
    or a nodal result."""
    values2 = []
    eids_missing = []
    nids_missing = []
    if key in elemental_quantities:
        # renumber eids
        for eid in seti2:
            if eid not in eid_map:
                eids_missing.append(eid)
                continue
            eid_new = eid_map[eid]
            values2.append(eid_new)
        #print('updating element SET %r' % options)
    else:
        # renumber nids
        for nid in seti2:
            if nid not in nid_map:
                nids_missing.append(nid)
                continue
            nid_new = nid_map[nid]
            values2.append(nid_new)
        #print('updating node SET %r' % options)
    return eids_missing, nids_missing, values2

#def _create_dict_mapper(properties, properties_map, pid_name, pid):
    #for pidi, prop in sorted(mydict.items()):
        #setattr(prop, pid_name, pid)
        #properties_map[pidi] = pid
        #pid += 1
    #return pid
