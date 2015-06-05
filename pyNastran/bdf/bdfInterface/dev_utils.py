from __future__ import print_function
from six import iteritems, itervalues, integer_types
from six.moves import zip, range
import scipy
from math import ceil
from pyNastran.bdf.bdf import BDF
from numpy import array, unique, where, arange, hstack, vstack, searchsorted, unique, log10
from pyNastran.bdf.cards.baseCard import expand_thru
from pyNastran.utils import object_attributes


def remove_unassociated_nodes(bdf_filename, bdf_filename_out, renumber=False):
    """
    Removes nodes from a model that are not referenced.

    .. warning only considers elements
    .. renumber=False is not supported
    """
    model = BDF()
    model.read_bdf(bdf_filename, xref=True)

    nids_used = set([])
    for element in itervalues(model.elements):
        nids_used.update(element.node_ids)
    all_nids = set(model.nodes.keys())

    nodes_to_remove = all_nids - nids_used
    #print('nodes_to_remove = %s' % nodes_to_remove)
    for nid in nodes_to_remove:
        del model.nodes[nid]
    model.write_bdf(bdf_filename_out)


def bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol):
    """
    Equivalences nodes; keeps the lower node id; creates two nodes with the same

    .. warning:: only handles CQUAD4, CTRIA3
    .. warning:: assumes cid=0
    """
    model = BDF()
    model.read_bdf(bdf_filename, xref=True)

    # quads / tris
    quads = []
    quadmap = []
    tris = []
    trimap = []

    inode = 1
    nid_map = {}
    for nid, node in sorted(iteritems(model.nodes)):
        node.nid = inode
        nid_map[inode - 1] = nid
        inode += 1
    #model.write_bdf('A_' + bdf_filename_out)

    nids = array([node.nid for nid, node in sorted(iteritems(model.nodes))], dtype='int32')
    #nids_used = set([])
    #for element in itervalues(model.elements):
    #    nids_used.update(element.node_ids)
    #nids_used = array(list(nids_used), dtype='int32')
    #nids_used.sort()
    #nids = nids_used

    nnodes = len(nids)
    i = arange(nnodes, dtype='int32')
    nids2 = vstack([i, nids]).T

    nodes_xyz = array([node.xyz for nid, node in sorted(iteritems(model.nodes))], dtype='float32')

    nids_new = set([])
    for eid, element in sorted(iteritems(model.elements)):
        emap = []

        if element.type == 'CQUAD4':
            quads.append(element.node_ids)
            quadmap.append(element.eid)
        elif element.type == 'CTRIA3':
            tris.append(element.node_ids)
            trimap.append(element.eid)
        else:
            raise NotImplementedError(element.type)

    quads = array(quads, dtype='int32') - 1
    quadmap = array(quadmap, dtype='int32')
    tris = array(tris, dtype='int32') - 1
    trimap = array(trimap, dtype='int32')

    # build the kdtree
    try:
        kdt = scipy.spatial.cKDTree(nodes_xyz)
    except RuntimeError:
        print(nodes_xyz)
        raise RuntimeError(nodes_xyz)

    # find the node ids of interest
    nids_new = hstack([unique(quads), unique(tris)])
    nids_new.sort()
    inew = searchsorted(nids, nids_new, side='left')

    # check the closest 10 nodes for equality
    deq, ieq = kdt.query(nodes_xyz[inew, :], k=10, distance_upper_bound=tol)

    # get the ids of the duplicate nodes
    slots = where(ieq[:, 1:] < nnodes)
    irows, icols = slots
    #replacer = unique(ieq[slots])

    skip_nodes = []
    for (irow, icol) in zip(irows, icols):
        nid1 = nid_map[irow]
        nid2 = nid_map[icol]

        node1 = model.nodes[nid1]
        node2 = model.nodes[nid2]

        node2.nid = node1.nid
        node2.xyz = node1.xyz
        assert node2.cp == node1.cp
        assert node2.cd == node1.cd
        assert node2.ps == node1.ps
        assert node2.seid == node1.seid
        skip_nodes.append(nid2)

    model.remove_nodes = skip_nodes
    #model._write_nodes = _write_nodes
    model.write_bdf(bdf_filename_out)


def cut_model(model, axis='-y'):
    """
    Removes the elements on one side of a model.

    Typically aircraft are defined as x-aft, y-right, z-up, so
    all node xyz locations are positive.  We then have a xz plane
    of symmetry with the axis of symmetry being y.

    Considers
    =========
      - nodes
      - elements
      - loads (LOAD/PLOAD4)

    ..note doesn't support +cuts (e.g. +x, +y, +z)
    """
    #model.read_bdf(bdf_filename, xref=False)
    #model.cross_reference(xref_loads=xref_loads)

    if axis == '-x':
        iaxis = 0
    elif axis == '-y':
        iaxis = 1
    elif axis == '-z':
        iaxis = 2
    else:
        raise NotImplementedError(axis)

    remove_nids = []
    for nid, node in iteritems(model.nodes):
        c = node.Position()
        if c[iaxis] < 0.0:
            remove_nids.append(nid)

    remove_eids = []
    for eid, element in iteritems(model.elements):
        c = element.Centroid()
        if c[iaxis] < 0.0:
            remove_eids.append(eid)

    for nid in remove_nids:
        del model.nodes[nid]
    for eid in remove_eids:
        del model.elements[eid]

    loads2 = {}
    for load_id, loadcase in iteritems(model.loads):
        loadcase2 = []
        for load in loadcase:
            if load.type == 'LOAD':
                loadcase2.append(load)
            elif load.type == 'PLOAD4':
                if load.eid not in remove_eids:
                    loadcase2.append(load)
            else:
                loadcase2.append(load)
        loads2[load_id] = loadcase2
    model.loads = loads2

def _write_nodes(self, outfile, size, is_double):
    """
    Writes the NODE-type cards

    :param self: the BDF object
    """
    if self.spoints:
        msg = []
        msg.append('$SPOINTS\n')
        msg.append(self.spoints.write_card(size, is_double))
        outfile.write(''.join(msg))

    if self.nodes:
        msg = []
        msg.append('$NODES\n')
        if self.gridSet:
            msg.append(self.gridSet.print_card(size))
        for (nid, node) in sorted(iteritems(self.nodes)):
            if nid not in self.remove_nodes:
                msg.append(node.write_card(size, is_double))
        outfile.write(''.join(msg))


def _roundup(x, n=100):
    return int(ceil(x / float(n))) * n


def bdf_renumber(bdf_filename, bdf_filename_out, size=8, is_double=False):
    """
    Supports
    ========
     - GRIDs
       - no superelements
     - COORDx

     - elements
        - CELASx/CONROD/CBAR/CBEAM/CQUAD4/CTRIA3/CTETRA/CPENTA/CHEXA
        - RBAR/RBAR1/RBE1/RBE2/RBE3

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


    ..warning:: spoints might be problematic...check
    ..warning:: still in development, but it usually brutally crashes if it's not supported
    ..warning:: be careful of unsupported cards
    """
    cid = 50
    nid = 101
    eid = 301
    pid = 401
    mid = 501
    spc_id = 501
    mpc_id = 601
    load_id = 701
    dload_id = 801

    method_id = 901
    cmethod_id = 1001
    spline_id = 1101
    table_id = 1201
    flfact_id = 1301
    flutter_id = 1401
    freq_id = 1501
    tstep_id = 1601
    tstepnl_id = 1701
    suport_id = 1801
    suport1_id = 1901
    tf_id = 2001

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

    model = BDF()
    model.read_bdf(bdf_filename)

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
    nnodes = len(spoints_nids)

    i = 1
    j = 1
    #print(spoints_nids)
    k = 0

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

    all_materials = (
        model.materials,
        model.creepMaterials,
        model.thermalMaterials,
        model.hyperelasticMaterials,
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
    mids = []
    for materials in all_materials:
        mids += materials.keys()
    mids = unique(mids)
    mids.sort()
    nmaterials = len(mids)

    for i in range(nmaterials):
        midi = mids[i]
        mid_map[midi] = mid + i

    #spoints2 = arange(1, len(spoints) + 1)
    for nid, node in sorted(iteritems(model.nodes)):
        nid_new = nid_map[nid]
        node.nid = nid_new

    # properties
    for pidi, prop in sorted(iteritems(model.properties)):
        prop.pid = pid
        pid += 1
    for pidi, prop in sorted(iteritems(model.properties_mass)):
        # PMASS
        prop.pid = pid
        pid += 1
    for pidi, prop in sorted(iteritems(model.convectionProperties)):
        # PCONV
        prop.pid = pid
        pid += 1
    for pidi, prop in sorted(iteritems(model.phbdys)):
        # PHBDY
        prop.pid = pid
        pid += 1

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
    for eidi, elem in sorted(iteritems(model.rigidElements)):
        # RBAR/RBAR1/RBE1/RBE2/RBE3
        elem.eid = eid
        eid_map[eidi] = eid
        eid += 1
    #for eidi, elem in iteritems(model.caeros):
        #pass

    #mid = 1
    for materials in all_materials:
        for midi, material in iteritems(materials):
            mid = mid_map[midi]
            #if midi == 15:
                #print(material)
            assert hasattr(material, 'mid')
            material.mid = mid

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
        (model.dconstrs, 'oid', dconstr_map),
        (model.dresps, 'oid', dresp_map),
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
        #(model.transfer_functions, 'sid', tranfer_function_map)
        #(model.bsets, 'sid', None),
        #(model.csets, 'sid', None),
        #(model.qsets, 'sid', None),
        #(model.usets, 'sid', None),

        (model.se_sets, 'sid', None),
        #(model.se_asets, 'sid', None),
        #(model.se_bsets, 'sid', None),
        #(model.se_csets, 'sid', None),
        #(model.se_qsets, 'sid', None),
        #(model.se_usets, 'sid', None),
    )
    param_id = 9999
    for (dict_obj, param_name, mmap) in sorted(data):
        param_id = _roundup(param_id, 1000) + 1
        for idi, param in sorted(iteritems(dict_obj)):
            msg = '%s has no %r; use %s' % (param.type, param_name, object_attributes(param))
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
    for table_idi, table in sorted(sorted(iteritems(model.randomTables))):
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
    model.write_bdf(bdf_filename_out, size=8, is_double=False,
                    interspersed=False)

def _update_case_control(model, mapper):
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
        ] + skip_keys_temp

    sets_analyzed = set([])
    # sets in the global don't get updated....
    # so we're going to find all the sets and
    # map them
    # TODO: renumber the sets
    set_locations = {}
    for isubcase, subcase in sorted(iteritems(model.caseControlDeck.subcases)):
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
    case_control = model.caseControlDeck
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
                        subcase.update_parameter_in_subcase(key, value2, options, param_type)

                    elif key in elemental_quantities + nodal_quantities:
                        #msg += '  allowed_keys=%s' % sorted(kmap.keys())
                        if seti in subcase:
                            seti2, seti_key = subcase.get_parameter(seti, msg=msg)
                            if seti_key in sets_analyzed:
                                continue
                            sets_analyzed.add(seti_key)
                            msgi = 'seti_key=%s must be an integer; type(seti_key)=%s\n'  %(seti_key, type(seti_key))
                            msgi += '  key=%r value=%r options=%r param_type=%r\n' % (key, value, options, param_type)
                            msgi += '  seti=%r\n' % seti
                            #print(msgi)
                            assert isinstance(seti_key, int), msgi

                            #print('key=%s options=%s param_type=%s value=%s' % (key, options, param_type, value))
                            #print(seti2)
                        else:
                            seti2 = [value]
                            print('key=%s options=%s param_type=%s value=%s' % (key, options, param_type, value))
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
                            # or not global_subcase.has_parameter(key)
                            gset = subcase.get_parameter(seti)
                            lset = subcase.get_parameter(seti)
                            #print('gset', gset)
                            #print('lset', lset)
                            if gset != lset:
                                asdf
                                subcase.update_parameter_in_subcase(seti, values2, seti_key, param_type)
                            else:
                                global_subcase.update_parameter_in_subcase(seti, values2, seti_key, param_type)
                            #subcase.update_parameter_in_subcase(seti, values2, seti_key, param_type)
                        elif not global_subcase.has_parameter(key):
                            subcase.update_parameter_in_subcase(seti, values2, seti_key, param_type)
                        else:
                            global_subcase.update_parameter_in_subcase(seti, values2, seti_key, param_type)
                    else:
                        pass
                        #print('key=%s seti2=%s' % (key, seti2))
                        print('key=%r options=%r param_type=%r value=%r' % (key, options, param_type, value))
                        raise RuntimeError(key)
                elif value in ['NONE', 'ALL']:
                    #print('*ALL -> key=%s options=%s param_type=%s value=%s' % (key, options, param_type, value))
                    #print('*all')
                    pass
                elif key == '':
                    pass
                else:
                    print('key=%s options=%s param_type=%s value=%s' % (key, options, param_type, value))
                    raise RuntimeError(key)

            else:
                raise RuntimeError(key)

                    #if value ==
        #print()


def main():
    msg = 'CEND\n'
    msg += 'BEGIN BULK\n'
    msg += 'GRID,1,,0.,0.,0.\n'
    msg += 'GRID,2,,0.,0.,0.5\n'
    msg += 'GRID,3,,0.,0.,0.51\n'
    msg += 'GRID,10,,0.,0.,1.\n'
    msg += 'GRID,11,,0.,0.,1.\n'
    msg += 'CTRIA3,1,1,1,2,11\n'
    msg += 'CTRIA3,2,1,1,2,11\n'
    msg += 'CTRIA3,3,1,2,3,11\n'
    msg += 'CTRIA3,4,1,1,2,10\n'
    msg += 'PSHELL,1,1,0.1\n'
    msg += 'MAT1,1,3.0,, 0.3\n'
    msg += 'ENDDATA'


    bdf_filename = 'nonunique.bdf'

    f = open(bdf_filename, 'wb')
    f.write(msg)
    f.close()

    bdf_filename_out = 'unique.bdf'
    tol = 0.2
    bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol)

if __name__ == '__main__':
    main()
