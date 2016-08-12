from __future__ import print_function
from math import ceil

from six import iteritems, PY2, string_types
import numpy as np

from pyNastran.bdf.bdf import BDF
from pyNastran.utils import integer_types, object_attributes


if PY2:
    import re
    _name_re = re.compile(r"[a-zA-Z_][a-zA-Z0-9_]*$")
    def isidentifier(string, dotted=False):
        """Is the string a valid Python variable name?"""
        assert dotted is False, dotted
        return bool(_name_re.match(string))
else:
    def isidentifier(string, dotted=False):
        """Is the string a valid Python variable name?"""
        assert dotted is False, dotted
        return string.isidentifier()


def bdf_renumber(bdf_filename, bdf_filename_out, size=8, is_double=False,
                 starting_id_dict=None, round_ids=False, cards_to_skip=None):
    """
    Renumbers a BDF

    Parameters
    ----------
    bdf_filename : str
        a bdf_filename (string; supported) or a BDF model (BDF)
        that has been cross referenced and is fully valid (a equivalenced deck is not valid)
    bdf_filename_out : str
        a bdf_filename to write
    size : int; {8, 16}; default=8
        the bdf write precision
    is_double : bool; default=False
        the field precision to write

    .. todo:: bdf_model option for bdf_filename hasn't been tested
    .. todo:: add support for subsets (e.g. renumber only a subset of nodes/elements)
    ..warning :: spoints might be problematic...check
    ..warning :: still in development, but it usually brutally crashes if it's not supported
    ..warning :: be careful of unsupported cards

    Supports
    ========
     - GRIDs
       - no superelements
     - COORDx

     - elements
        - CELASx/CONROD/CBAR/CBEAM/CQUAD4/CTRIA3/CTETRA/CPENTA/CHEXA
        - RBAR/RBAR1/RBE1/RBE2/RBE3/RSPLINE

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

    Example 1 - Renumber Everything; Start from 1
    ---------------------------------------------
    bdf_renumber(bdf_filename, bdf_filename_out, size=8, is_double=False,
                 round_ids=False)

    Example 2 - Renumber Everything; Start Material IDs from 100
    ------------------------------------------------------------
    starting_id_dict = {
        'mid' : 100,
    }
    bdf_renumber(bdf_filename, bdf_filename_out, size=8, is_double=False,
                 starting_ids_dict=starting_ids_dict, round_ids=False)

    Example 3 - Only Renumber Material IDs
    --------------------------------------
    starting_id_dict = {
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
    }
    bdf_renumber(bdf_filename, bdf_filename_out, size=8, is_double=False,
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
        'dload_id' : 1,

        'method_id' : 1,
        'cmethod_id' : 1,
        'spline_id' : 1,
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
    if starting_id_dict is None:
        starting_id_dict = starting_id_dict_default
    else:
        for key, value in iteritems(starting_id_dict_default):
            if key not in starting_id_dict:
                starting_id_dict[key] = value

    for key, value in sorted(iteritems(starting_id_dict)):
        assert isinstance(key, string_types), key
        assert key in starting_id_dict_default, 'key=%s is invalid' % (key)
        assert isidentifier(key), 'key=%s is invalid' % key
        if value is None:
            pass
        else:
            if not isinstance(value, integer_types):
                msg = 'value=%s must be an integer; type(value)=%s' % (value, type(value))
                raise TypeError(msg)
        call = '%s = %s' % (key, value)

        # this exec is safe because we checked the identifier
        exec(call)


    eid_map = {}
    nid_map = {}
    reverse_nid_map = {}
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
    suport_map = {}
    suport1_map = {}

    if isinstance(bdf_filename, string_types):
        model = BDF(debug=False)
        model.disable_cards(cards_to_skip)
        model.read_bdf(bdf_filename)
    else:
        model = bdf_filename

    if model.spoints is None:
        spoints = []
    else:
        spoints = list(model.spoints.points)
    if model.epoints is None:
        epoints = []
    else:
        epoints = list(model.epoints.points)

    nids = model.nodes.keys()

    spoints_nids = spoints + nids
    spoints_nids.sort()
    i = 1
    #nnodes = len(spoints_nids)

    i = 1
    #j = 1
    #print(spoints_nids)
    #k = 0

    if 'nid' in starting_id_dict and nid is not None:
        i = nid
        #banned_nodes = spoints
        for nidi in spoints_nids:
            if nidi in spoints:
                pass
                #print('sid=%s -> %s' % (nid, i))
                #i += 1
            else:
                while i in spoints:
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
        for nid in spoints_nids:
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
        for nid, node in sorted(iteritems(model.nodes)):
            nid_new = nid_map[nid]
            node.nid = nid_new

    if 'pid' in starting_id_dict and pid is not None:
        # properties
        for pidi, prop in sorted(iteritems(model.properties)):
            prop.pid = pid
            pid += 1
        for pidi, prop in sorted(iteritems(model.properties_mass)):
            # PMASS
            prop.pid = pid
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
        for eidi, element in sorted(iteritems(model.elements)):
            element.eid = eid
            eid_map[eidi] = eid
            eid += 1
        for eidi, element in sorted(iteritems(model.masses)):
            # CONM1, CONM2, CMASSx
            element.eid = eid
            eid_map[eidi] = eid
            eid += 1
        for eidi, elem in sorted(iteritems(model.rigid_elements)):
            # RBAR/RBAR1/RBE1/RBE2/RBE3/RSPLINE
            elem.eid = eid
            eid_map[eidi] = eid
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
        for spc_idi, spc_group in sorted(iteritems(model.spcs)):
            for i, spc in enumerate(spc_group):
                assert hasattr(spc, 'conid')
                spc.conid = spc_id
            spc_map[spc_idi] = spc_id
            spc_id += 1
        for spc_idi, spcadd in sorted(iteritems(model.spcadds)):
            assert hasattr(spcadd, 'conid')
            spcadd.conid = spc_id
            spc_map[spc_idi] = spc_id
            spc_id += 1
    else:
        for spc_id in model.spcs:
            spc_map[spc_id] = spc_id
        for spc_id in model.spcadds:
            spc_map[spc_id] = spc_id

    if 'mpc_id' in starting_id_dict and mpc_id is not None:
        # mpc
        for mpc_idi, mpc_group in sorted(iteritems(model.mpcs)):
            for i, mpc in enumerate(mpc_group):
                assert hasattr(mpc, 'conid')
                mpc.conid = mpc_id
            mpc_map[mpc_idi] = mpc_id
            mpc_id += 1
        for mpc_idi, mpcadd in sorted(iteritems(model.mpcadds)):
            assert hasattr(mpcadd, 'conid')
            mpcadd.conid = mpc_id
            mpc_map[mpc_idi] = mpc_id
            mpc_id += 1
    else:
        for mpc_id in model.mpcs:
            mpc_map[mpc_id] = mpc_id
        for mpc_id in model.mpcadds:
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
        (model.frequencies, 'sid', freq_map),
        (model.tsteps, 'sid', tstep_map),
        (model.tstepnls, 'sid', tstepnl_map),
        (model.splines, 'eid', None),
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

        (model.sets, 'sid', None),
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
    param_id = 9999
    for (dict_obj, param_name, mmap) in sorted(data):
        if round_ids:
            param_id = _roundup(param_id, 1000) + 1
        else:
            param_id = 1
        for idi, param in sorted(iteritems(dict_obj)):
            try:
                msg = '%s has no %r; use %s' % (param.type, param_name, object_attributes(param))
            except AttributeError:
                print('param = %r' % param)
                raise
            assert hasattr(param, param_name), msg
            setattr(param, param_name, param_id)
            if mmap is not None:
                mmap[idi] = param_id
            param_id += 1

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
        'nodes' : nid_map,
        'coords' : cid_map,
        'materials' : mid_map,
        'SPC' : spc_map,
        'MPC' : mpc_map,
        'METHOD' : method_map,
        'CMETHOD' : cmethod_map,
        'FLFACT' : flfact_map,
        'FMETHOD' : flutter_map,
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
    #print('****suport1_map', suport1_map)
    #print('****dessub_map', dessub_map)
    #print('****dresp_map', dresp_map)
    _update_case_control(model, mapper)
    if bdf_filename_out is not None:
        model.write_bdf(bdf_filename_out, size=size, is_double=is_double,
                        interspersed=False)
    return model


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
                        subcase.update_parameter_in_subcase(
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
                            print('key=%s options=%s param_type=%s value=%s' % (
                                key, options, param_type, value))
                            raise NotImplementedError(key)

                        values2 = []
                        if key in elemental_quantities:
                            # renumber eids
                            for eid in seti2:
                                if eid not in eid_map:
                                    print("  couldn't find eid=%s...dropping" % eid)
                                    continue
                                eid_new = eid_map[eid]
                                values2.append(eid_new)
                            #print('updating element SET %r' % options)
                        else:
                            # renumber nids
                            for nid in seti2:
                                if nid not in nid_map:
                                    print("  couldn't find nid=%s...dropping" % nid)
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
                                #subcase.update_parameter_in_subcase(
                                    #seti, values2, seti_key, param_type)
                            else:
                                global_subcase.update_parameter_in_subcase(
                                    seti, values2, seti_key, param_type)
                            #subcase.update_parameter_in_subcase(
                                #seti, values2, seti_key, param_type)
                        elif not key not in global_subcase:
                            subcase.update_parameter_in_subcase(
                                seti, values2, seti_key, param_type)
                        else:
                            global_subcase.update_parameter_in_subcase(
                                seti, values2, seti_key, param_type)
                    else:
                        #pass
                        #print('key=%s seti2=%s' % (key, seti2))
                        print('key=%r options=%r param_type=%r value=%r' % (
                            key, options, param_type, value))
                        raise RuntimeError(key)
                elif value in ['NONE', 'ALL']:
                    # print('*ALL -> key=%s options=%s param_type=%s value=%s' % (
                        # key, options, param_type, value))
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

def _roundup(value, round_increment=100):
    """
    Rounds up to the next N.

    Parameters
    ----------
    value : int
        the value to round up
    round_increment : int
        the increment to round by

    .. python

      >>> 100 = _roundup(10)
      >>> 200 = _roundup(105)
      >>> 300 = _roundup(200)
      >>> 1000 = _roundup(200, 1000)
      >>> 2000 = _roundup(1000, 1000)
      >>> 2000 = _roundup(1001, 1000)

    .. note :: this function is used to ensure that renumbering is more
               obvious when testing
    """
    return int(ceil(value / float(round_increment))) * round_increment
