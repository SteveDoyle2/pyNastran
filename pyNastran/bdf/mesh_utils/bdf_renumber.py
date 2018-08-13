"""
defines:
    bdf_renumber(bdf_filename, bdf_filename_out, size=8, is_double=False,
                 starting_id_dict=None, round_ids=False, cards_to_skip=None)
"""
from __future__ import print_function
from itertools import chain

import io
from six import PY2, PY3, iteritems, itervalues, StringIO
import numpy as np

from pyNastran.bdf.bdf import BDF
from pyNastran.utils import integer_types, object_attributes
from pyNastran.utils.mathematics import roundup


def bdf_renumber(bdf_filename, bdf_filename_out, size=8, is_double=False,
                 starting_id_dict=None, round_ids=False, cards_to_skip=None,
                 log=None, debug=False):
    """
    Renumbers a BDF

    Parameters
    ----------
    bdf_filename : str / BDF
        str : a bdf_filename (string; supported)
        BDF : a BDF model that has been cross referenced and is
        fully valid (an equivalenced deck is not valid)
    bdf_filename_out : str
        a bdf_filename to write
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

    .. todo:: bdf_model option for bdf_filename hasn't been tested
    .. todo:: add support for subsets (e.g. renumber only a subset of nodes/elements)
    .. todo:: doesn't support partial renumbering
    .. todo:: doesn't support element material coordinate systems

    ..warning :: spoints might be problematic...check
    ..warning :: still in development, but it usually brutally crashes
                 if it's not supported
    ..warning :: be careful of card unsupported cards (e.g. ones not read in)

    Supports
    ========
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
    ========
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
        for key, value in iteritems(starting_id_dict_default):
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
    tstep_id = None
    tstepnl_id = None
    set_id = None
    suport_id = None
    suport1_id = None
    tf_id = None

    # turn them into variables
    for key, value in sorted(iteritems(starting_id_dict)):
        #assert isinstance(key, string_types), key
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

    # build the maps
    mass_id_map = {}
    nid_map = {}
    properties_map = {}
    properties_mass_map = {}
    reverse_nid_map = {}
    eid_map = {}
    rigid_elements_map = {}
    nsm_map = {}
    mid_map = {}
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

    if isinstance(bdf_filename, BDF):
        model = bdf_filename
    else:
        model = BDF(log=log, debug=debug)
        model.disable_cards(cards_to_skip)
        model.read_bdf(bdf_filename)

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
            reverse_nid_map[nid] = nid

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

    if 'nid' in starting_id_dict and nid is not None:
        #spoints2 = arange(1, len(spoints) + 1)
        #nid = _create_dict_mapper(model.nodes, nid_map, 'nid', nid)

        for nid, node in sorted(iteritems(model.nodes)):
            nid_new = nid_map[nid]
            #print('nid=%s -> %s' % (nid, nid_new))
            node.nid = nid_new

    if 'pid' in starting_id_dict and pid is not None:
        # properties
        #pid = _create_dict_mapper(model.properties, properties_map, 'pid', pid)
        #pid = _create_dict_mapper(model.properties_mass, properties_mass_map, 'pid', pid)
        #pid = _update(model.convection_properties, properties_mass_map, pid)

        for pidi, prop in sorted(iteritems(model.properties)):
            prop.pid = pid
            properties_map[pidi] = pid
            pid += 1
        for pidi, prop in sorted(iteritems(model.properties_mass)):
            # PMASS
            prop.pid = pid
            properties_mass_map[pidi] = pid
            pid += 1
        for pidi, prop in sorted(iteritems(model.convection_properties)):
            # PCONV
            prop.pid = pid
            pid += 1
        for pidi, prop in sorted(iteritems(model.phbdys)):
            # PHBDY
            prop.pid = pid
            pid += 1

    if 'eid' in starting_id_dict and eid is not None:
        # elements
        #eid = _create_dict_mapper(model.elements, eid_map, 'eid', eid)

        for eidi, element in sorted(iteritems(model.elements)):
            element.eid = eid
            eid_map[eidi] = eid
            eid += 1
        for eidi, element in sorted(iteritems(model.masses)):
            # CONM1, CONM2, CMASSx
            element.eid = eid
            eid_map[eidi] = eid
            mass_id_map[eidi] = eid
            eid += 1
        for eidi, elem in sorted(iteritems(model.rigid_elements)):
            # RBAR/RBAR1/RBE1/RBE2/RBE3/RSPLINE/RSSCON
            elem.eid = eid
            eid_map[eidi] = eid
            rigid_elements_map[eidi] = eid
            eid += 1
        #for eidi, elem in iteritems(model.caeros):
            #pass

    if 'mid' in starting_id_dict and mid is not None:
        #mid = 1
        for materials in all_materials:
            for midi, material in iteritems(materials):
                mid = mid_map[midi]
                assert hasattr(material, 'mid')
                material.mid = mid

    if 'spc_id' in starting_id_dict and spc_id is not None:
        # spc
        for spc_idi, spc_group in sorted(iteritems(model.spcadds)):
            for i, spc in enumerate(spc_group):
                assert hasattr(spc, 'conid')
                spc.conid = spc_id
            spc_map[spc_idi] = spc_id
            spc_id += 1
        for spc_idi, spc_group in sorted(iteritems(model.spcs)):
            for i, spc in enumerate(spc_group):
                assert hasattr(spc, 'conid')
                spc.conid = spc_id
            spc_map[spc_idi] = spc_id
            spc_id += 1
    else:
        for spc_id in model.spcadds:
            spc_map[spc_id] = spc_id
        for spc_id in model.spcs:
            spc_map[spc_id] = spc_id

    if 'mpc_id' in starting_id_dict and mpc_id is not None:
        # mpc
        for mpc_idi, mpc_group in sorted(iteritems(model.mpcadds)):
            for i, mpc in enumerate(mpc_group):
                assert hasattr(mpc, 'conid')
                mpc.conid = mpc_id
            mpc_map[mpc_idi] = mpc_id
            mpc_id += 1
        for mpc_idi, mpc_group in sorted(iteritems(model.mpcs)):
            for i, mpc in enumerate(mpc_group):
                assert hasattr(mpc, 'conid')
                mpc.conid = mpc_id
            mpc_map[mpc_idi] = mpc_id
            mpc_id += 1
    else:
        for mpc_id in model.mpcadds:
            mpc_map[mpc_id] = mpc_id
        for mpc_id in model.mpcs:
            mpc_map[mpc_id] = mpc_id

    if 'cid' in starting_id_dict and cid is not None:
        # coords
        for cidi, coord in sorted(iteritems(model.coords)):
            if cidi == 0:
                cid_map[0] = 0
                continue
            coord.cid = cid
            cid_map[cidi] = cid
            cid += 1

    if 'freq_id' in starting_id_dict and freq_id is not None:
        # frequencies
        for freqi, freqs in sorted(iteritems(model.frequencies)):
            freq_map[freqi] = freq_id
            for freq in freqs:
                freq.sid = freqi
            freq_id += 1

    set_map = {}
    if 'set_id' in starting_id_dict and set_id is not None:
        # sets
        for sidi, set_ in sorted(iteritems(model.sets)):
            set_.sid = set_id
            set_map[sidi] = set_id
            set_id += 1

    if 'spline_id' in starting_id_dict and spline_id is not None:
        # set up spline1 box mapping
        delta_box1_map = {}
        delta_box2_map = {}
        for sidi, spline in sorted(iteritems(model.splines)):
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
        for caero_idi, caero in sorted(iteritems(model.caeros)):
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
        for sidi, spline in sorted(iteritems(model.splines)):
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
        (model.transfer_functions, 'sid', tranfer_function_map)
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
        for idi, param in sorted(iteritems(dict_obj)):
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
    for key, value in iteritems(dconstr_map):
        if key in dessub_map:
            raise NotImplementedError()
        dessub_map[key] = value

    # tables
    for table_idi, table in sorted(sorted(iteritems(model.tables))):
        assert hasattr(table, 'tid')
        table.tid = table_id
        table_id += 1
    for table_idi, table in sorted(sorted(iteritems(model.random_tables))):
        assert hasattr(table, 'tid')
        table.tid = table_id
        table_id += 1

    # dloads
    for dload_idi, dloads in sorted(iteritems(model.dloads)):
        for dload in dloads:
            assert hasattr(dload, 'sid')
            dload.sid = dload_id
        dload_map[dload_idi] = dload_id
        dload_id += 1
    for dload_idi, dloads in sorted(iteritems(model.dload_entries)):
        for dload in dloads:
            assert hasattr(dload, 'sid')
            dload.sid = dload_id
        dload_map[dload_idi] = dload_id
        dload_id += 1

    # loads
    for load_idi, load_combinations in sorted(iteritems(model.load_combinations)):
        for load_combination in load_combinations:
            assert hasattr(load_combination, 'sid')
            load_combination.sid = load_id
        load_map[load_idi] = load_id
        load_id += 1
    for load_idi, loads in sorted(iteritems(model.loads)):
        for load in loads:
            assert hasattr(load, 'sid')
            load.sid = load_id
        load_map[load_idi] = load_id
        load_id += 1

    # transfer_functions
    for tf_idi, tfs in sorted(iteritems(model.transfer_functions)):
        for tf in tfs:
            assert hasattr(tf, 'sid')
            tf.sid = tf_id
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
        'SPC' : spc_map,
        'MPC' : mpc_map,   # TODO: come up with unified system that uses the same key
        'mpcs' : mpc_map,  #       for bdf_merge and _update_case_control
        'METHOD' : method_map,
        'CMETHOD' : cmethod_map,
        'FLFACT' : flfact_map,
        'FMETHOD' : flutter_map,
        'FREQUENCY' : freq_map,
        'sets' : set_map,
        'splines' : spline_id_map,
        'caeros' : caero_id_map,

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
    #print('****suport1_map', suport1_map)
    #print('****dessub_map', dessub_map)
    #print('****dresp_map', dresp_map)
    _update_case_control(model, mapper)
    if bdf_filename_out is not None:
        close = True
        if PY2 and isinstance(bdf_filename_out, (file, StringIO)):
            close = False
        elif PY3 and isinstance(bdf_filename_out, io.IOBase):
            close = False
        model.write_bdf(bdf_filename_out, size=size, is_double=is_double,
                        interspersed=False, close=close)
    return model, mapper

def get_renumber_starting_ids_from_model(model):
    """
    Get the starting ids dictionary used for renumbering with ids greater than those in model.
    Parameters
    -----------
    model : BDF
        BDF object to get maximum ids from.
    Returns
    --------
    starting_id_dict : dict {str : int, ...}
        Dictionary from id type to starting id.
    """
    starting_id_dict = {
        'cid' : max(model.coords.keys()) + 1,
        'nid' : max(model.point_ids) + 1,
        'eid' : max([max(model.elements.keys()),
                     max(model.masses.keys()) if model.masses else 0,
                     max(model.rigid_elements.keys()) if model.rigid_elements else 0,
                     ]) + 1,
        'pid' : max([max(model.properties.keys()),
                     0 if len(model.properties_mass) == 0 else max(model.properties_mass.keys()),
                     ]) + 1,
        'mid' : max(model.material_ids) + 1,
        'set_id' : max(model.sets.keys()) + 1 if model.sets else 1,
        'spline_id' : max(model.splines.keys()) + 1 if model.splines else 1,
        'caero_id' : max(caero.box_ids[-1, -1]
                         for caero in itervalues(model.caeros)) + 1 if model.caeros else 1,
    }
    return starting_id_dict

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
    elemental_quantities = ['STRESS', 'STRAIN', 'FORCE', 'ESE', 'EKE']
    nodal_quantities = [
        'DISPLACEMENT', 'VELOCITY', 'ACCELERATION', 'SPCFORCES', 'MPCFORCES',
        'GPFORCES', 'SDISPLACEMENT', 'OLOAD',
    ]
    mapper_quantities = [
        'FREQUENCY', 'DLOAD', 'LOAD', 'LOADSET', 'SPC', 'MPC', 'METHOD', 'CMETHOD',
        'TSTEP', 'TSTEPNL', 'NLPARM', 'SDAMPING', 'DESSUB', 'DESOBJ', 'GUST',
        'SUPORT1', 'TRIM', 'BOUTPUT', 'IC', 'CSSCHD', 'FMETHOD', 'TFL',
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
        'LINE',
        ] + skip_keys_temp

    sets_analyzed = set([])
    # sets in the global don't get updated....
    # so we're going to find all the sets and
    # map them
    # TODO: renumber the sets
    set_locations = {}
    case_control = model.case_control_deck
    if case_control is None:
        return

    for isubcase, subcase in sorted(iteritems(case_control.subcases)):
        for key, values in sorted(iteritems(subcase.params)):
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
    for isubcase, subcase in sorted(iteritems(case_control.subcases)):
        #print('-----------------------')
        #print(subcase)
        for key, values in sorted(iteritems(subcase.params)):
            if key in skip_keys:
                pass
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
                        kmap = mapper[key]
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

                        values2 = []
                        if key in elemental_quantities:
                            # renumber eids
                            for eid in seti2:
                                if eid not in eid_map:
                                    model.log.warning("  couldn't find eid=%s...dropping" % eid)
                                    continue
                                eid_new = eid_map[eid]
                                values2.append(eid_new)
                            #print('updating element SET %r' % options)
                        else:
                            # renumber nids
                            for nid in seti2:
                                if nid not in nid_map:
                                    model.log.warning("  couldn't find nid=%s...dropping" % nid)
                                    continue
                                nid_new = nid_map[nid]
                                values2.append(nid_new)
                            #print('updating node SET %r' % options)

                        param_type = 'SET-type'
                        #print('adding seti=%r values2=%r seti_key=%r param_type=%r'  % (
                            #seti, values2, seti_key, param_type))
                        assert len(values2) > 0, values2
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

#def _create_dict_mapper(properties, properties_map, pid_name, pid):
    #for pidi, prop in sorted(iteritems(mydict)):
        #setattr(prop, pid_name, pid)
        #properties_map[pidi] = pid
        #pid += 1
    #return pid
