from __future__ import print_function
from six import iteritems, itervalues, integer_types, PY2
from six.moves import zip, range

import os
from itertools import count
from math import ceil

from numpy import array, unique, where, arange, hstack, vstack, searchsorted, unique, log10, array_equal
from numpy.linalg import norm
import scipy

from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.cards.baseCard import expand_thru
from pyNastran.utils import object_attributes


if PY2:
    import re
    _name_re = re.compile(r"[a-zA-Z_][a-zA-Z0-9_]*$")
    def isidentifier(s, dotted=False):
        return bool(_name_re.match(s))
else:
    def isidentifier(s, dotted=False):
        return s.isidentifier()


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
    #for element in itervalues(model.masses):
        #nids_used.update(element.node_ids)
    all_nids = set(model.nodes.keys())

    nodes_to_remove = all_nids - nids_used
    for nid in nodes_to_remove:
        del model.nodes[nid]

    if renumber:
        starting_id_dict = {
            'nid' : 1,
            'eid' : 1,
            'pid' : 1,
            'mid' : 1,
        }
        bdf_renumber(model, bdf_filename_out, size=8, is_double=False,
                     starting_id_dict=starting_id_dict)
    else:
        model.write_bdf(bdf_filename_out)


def bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                          renumber_nodes=False, neq_max=4, xref=True,
                          node_set=None, crash_on_collapse=False):
    """
    Equivalences nodes; keeps the lower node id; creates two nodes with the same

    :param bdf_filename: a bdf_filename (string) or a BDF model (BDF)
        that is fully valid (see xref)
    :param bdf_filename_out: a bdf_filename to write
    :param tol:              the spherical tolerance (float)
    :param renumber_nodes:   should the nodes be renumbered (default=False; not supported)
    :param neq_max:          the number of "close" points (default=4)
    :param xref:             does the model need to be cross_referenced
                             (default=True; only applies to model option)
    :param nodes_set:        the list/array of nodes to consider (not supported)
    :param crash_on_collapse: stop if nodes have been collapsed
                              (default=False;
                               False: blindly move on
                               True: rereads the BDF which catches doubled nodes (temporary);
                                     in the future collapse=True won't need to double read;
                                     an alternative is to do Patran's method of avoiding collapse)

    .. warning:: I doubt SPOINTs/EPOINTs work correctly
    .. warning:: xref not fully implemented (assumes cid=0)
    """
    assert isinstance(tol, float), tol
    if isinstance(bdf_filename, str):
        xref = True
        model = BDF()
        model.read_bdf(bdf_filename, xref=True)
    else:
        model = bdf_filename
        model.cross_reference(xref=xref)

    coord_ids = model.coord_ids
    needs_get_position = True if coord_ids == [0] else False

    # quads / tris
    nids_quads = []
    eids_quads = []
    nids_tris = []
    eids_tris = []

    # map the node ids to the slot in the nids array
    renumber_nodes = False

    inode = 0
    nid_map = {}
    if renumber_nodes:
        for nid, node in sorted(iteritems(model.nodes)):
            node.nid = inode + 1
            nid_map[inode] = nid
            inode += 1
        nnodes = len(model.nodes)
        nids = arange(1, inode + 1, dtype='int32')
        assert nids[-1] == nnodes

    else:
        for nid, node in sorted(iteritems(model.nodes)):
            nid_map[inode] = nid
            inode += 1
        nids = array([node.nid for nid, node in sorted(iteritems(model.nodes))], dtype='int32')
    #model.write_bdf('A_' + bdf_filename_out)

    #nids_used = set([])
    #for element in itervalues(model.elements):
    #    nids_used.update(element.node_ids)
    #nids_used = array(list(nids_used), dtype='int32')
    #nids_used.sort()
    #nids = nids_used

    nnodes = len(nids)
    i = arange(nnodes, dtype='int32')
    #nids2 = vstack([i, nids]).T

    if needs_get_position:
        nodes_xyz = array([node.get_position()
                           for nid, node in sorted(iteritems(model.nodes))], dtype='float32')
    else:
        nodes_xyz = array([node.xyz
                           for nid, node in sorted(iteritems(model.nodes))], dtype='float32')

    if 0:
        i = 0
        for nid, xyz in zip(nids, nodes_xyz):
            #print('%i %-5s %s' % (i, nid, list_print(xyz)))
            i += 1

    # there is some set of points that are used on the elements that
    # will be considered.
    #
    # Presumably this is enough to capture all the node ids and NOT
    # spoints, but I doubt it...
    spoint_epoint_nid_set = set([])
    for eid, element in sorted(iteritems(model.elements)):
        spoint_epoint_nid_set.update(element.node_ids)
    for eid, element in sorted(iteritems(model.masses)):
        spoint_epoint_nid_set.update(element.node_ids)

    if model.spoints and model.epoints:
        nids_new = spoint_epoint_nid_set - model.spoints.points - model.epoints.points
    elif model.spoints:
        nids_new = spoint_epoint_nid_set - model.spoints.points
    elif model.epoints:
        nids_new = spoint_epoint_nid_set - model.epoints.points
    else:
        nids_new = spoint_epoint_nid_set

    if None in nids_new:
        nids_new.remove(None)

    # autosorts the data
    nids_new = unique(list(nids_new))
    assert isinstance(nids_new[0], integer_types), type(nids_new[0])

    missing_nids = list(set(nids_new) - set(nids))
    if missing_nids:
        missing_nids.sort()
        msg = 'There are missing nodes...\n'
        msg = 'missing nids=%s' % str(missing_nids)
        raise RuntimeError(msg)

    # get the node_id mapping for the kdtree
    inew = searchsorted(nids, nids_new, side='left')
    #assert array_equal(nids[inew], nids_new), 'some nodes are not defined'

    # build the kdtree
    try:
        kdt = scipy.spatial.cKDTree(nodes_xyz)
    except RuntimeError:
        print(nodes_xyz)
        raise RuntimeError(nodes_xyz)

    # check the closest 10 nodes for equality
    deq, ieq = kdt.query(nodes_xyz[inew, :], k=neq_max, distance_upper_bound=tol)

    # get the ids of the duplicate nodes
    slots = where(ieq[:, :] < nnodes)
    irows, icols = slots
    #replacer = unique(ieq[slots])  ## TODO: turn this back on?

    skip_nodes = []
    for (islot, irow, icol) in zip(count(), irows, icols):
        inid2 = ieq[irow, icol]
        nid1 = nids[irow]
        nid2 = nids[inid2]

        #if nid1 == nid2:
            #continue

        node1 = model.nodes[nid1]
        node2 = model.nodes[nid2]

        # TODO: doesn't use get position...
        R = norm(node1.xyz - node2.xyz)

        #print('  irow=%s->n1=%s icol=%s->n2=%s' % (irow, nid1, icol, nid2))
        if R > tol:
            #print('  *n1=%-4s xyz=%s\n  *n2=%-4s xyz=%s\n  *R=%s\n' % (
            #    nid1, list_print(node1.xyz),
            #    nid2, list_print(node2.xyz),
            #    R))
            continue
        #print('  n1=%-4s xyz=%s\n  n2=%-4s xyz=%s\n  R=%s\n' % (
        #    nid1, list_print(node1.xyz),
        #    nid2, list_print(node2.xyz),
        #    R))

        node2.nid = node1.nid
        node2.xyz = node1.xyz
        node2.cp = node1.cp
        assert node2.cd == node1.cd
        assert node2.ps == node1.ps
        assert node2.seid == node1.seid
        skip_nodes.append(nid2)

    if bdf_filename_out is not None:
        model.write_bdf(bdf_filename_out)
    if crash_on_collapse:
        # lazy way to make sure there aren't any collapsed nodes
        model2 = BDF()
        model.read_bdf(bdf_filename_out)
    return model


def cut_model(model, axis='-y'):
    """
    Removes the elements on one side of a model.

    Typically aircraft are defined as x-aft, y-right, z-up, so
    all node xyz locations are positive.  We then have a xz plane
    of symmetry with the axis of symmetry being y, and typically
    save the +y elements.

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
        if c[iaxis] <= 0.0:
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


def bdf_merge(bdf_filenames, bdf_filenames_out=None, renumber=True):
    """
    .. todo :: doesn't support SPOINTs/EPOINTs
    .. warning :: still very preliminary
    """
    if isinstance(bdf_filenames, str):
        bdf_filenames = [bdf_filenames]
    elif not (isinstance(bdf_filenames, list) or isinstance(bdf_filenames, tuple)):
        raise TypeError(bdf_filenames)

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
    model = BDF()
    bdf_filename0 = bdf_filenames[0]
    model.read_bdf(bdf_filename0)
    print('primary=%s' % bdf_filename0)

    data_members = [
        'coords', 'nodes', 'elements', 'masses', 'properties',  'properties_mass',
        'materials',
    ]
    for bdf_filename in bdf_filenames[1:]:
        print('model.masses = %s' % model.masses)
        starting_id_dict = {
            'cid' : max(model.coords.keys()) + 1,
            'nid' : max(model.nodes.keys()) + 1,
            'eid' : max([
                max(model.elements.keys()),
                0 if len(model.masses) == 0 else max(model.masses.keys()),
            ]) + 1,
            'pid' : max([
                max(model.properties.keys()),
                0 if len(model.properties_mass) == 0 else max(model.properties_mass.keys()),
            ]) + 1,
            'mid' : max(model.material_ids) + 1,
        }
        #for param, val in sorted(iteritems(starting_id_dict)):
            #print('  %-3s %s' % (param, val))

        print('secondary=%s' % bdf_filename)
        model2 = BDF()
        bdf_dump = 'bdf_merge_temp.bdf'
        #model2.read_bdf(bdf_filename, xref=False)

        bdf_renumber(bdf_filename, bdf_dump, starting_id_dict=starting_id_dict)
        model2 = BDF()
        model2.read_bdf(bdf_dump)
        os.remove(bdf_dump)

        print('model2.node_ids = %s' % array(model2.node_ids))
        for data_member in data_members:
            data1 = getattr(model, data_member)
            data2 = getattr(model2, data_member)
            if isinstance(data1, dict):
                print('  working on %s' % (data_member))
                for key, value in iteritems(data2):
                    if data_member in 'coords' and key == 0:
                        continue
                    if isinstance(value, list):
                        raise NotImplementedError(type(value))
                    else:
                        assert key not in data1, key
                        data1[key] = value
                        #print('   %s' % key)
            else:
                raise NotImplementedError(type(data1))
    #if bdf_filenames_out:
        #model.write_bdf(bdf_filenames_out)

    if renumber:
        print('final renumber...')
        starting_id_dict = {
            'cid' : 1,
            'nid' : 1,
            'eid' : 1,
            'pid' : 1,
            'mid' : 1,
        }
        bdf_renumber(model, bdf_filenames_out, starting_id_dict=starting_id_dict)

    return model


def bdf_renumber(bdf_filename, bdf_filename_out, size=8, is_double=False,
                 starting_id_dict=None):
    """
    Renumbers a BDF

    :param bdf_filename: a bdf_filename (string; supported) or a BDF model (BDF)
        that has been cross referenced and is fully valid (a equivalenced deck is not valid)
    :param bdf_filename_out: a bdf_filename to write
    :param size:       the field size to write (default=8; 8 or 16)
    :param is_double:  the field precision to write (default=True)

    ..todo :: bdf_model option for bdf_filename hasn't been tested
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
    """
    starting_id_dict_default = {
        'cid' : 50,
        'nid' : 101,
        'eid' : 301,
        'pid' : 401,
        'mid' : 501,
        'spc_id' : 501,
        'mpc_id' : 601,
        'load_id' : 701,
        'dload_id' : 801,

        'method_id' : 901,
        'cmethod_id' : 1001,
        'spline_id' : 1101,
        'table_id' : 1201,
        'flfact_id' : 1301,
        'flutter_id' : 1401,
        'freq_id' : 1501,
        'tstep_id' : 1601,
        'tstepnl_id' : 1701,
        'suport_id' : 1801,
        'suport1_id' : 1901,
        'tf_id' : 2001,
    }
    if starting_id_dict is None:
        starting_id_dict = starting_id_dict_default
    else:
        for key, value in iteritems(starting_id_dict_default):
            if key not in starting_id_dict:
                starting_id_dict[key] = value

    for key, value in sorted(iteritems(starting_id_dict)):
        assert isinstance(key, str), key
        assert key in starting_id_dict_default, 'key=%s is invalid' % (key)
        assert isidentifier(key), 'key=%s is invalid' % key
        assert isinstance(value, integer_types), 'value=%s must be an integer; type(value)' % (value, type(value))
        call = '%s = %s' % (key, value)
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

    if isinstance(bdf_filename, str):
        model = BDF()
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
    if bdf_filename_out is not None:
        model.write_bdf(bdf_filename_out, size=8, is_double=False,
                        interspersed=False)
    return model


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
    case_control = model.caseControlDeck
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


def eq2():
    lines = [
        '$pyNastran: version=msc',
        '$pyNastran: punch=True',
        '$pyNastran: encoding=ascii',
        '$NODES',
        '$ Nodes to merge:',
        '$ 5987 10478',
        '$   GRID        5987           35.46     -6.      0.',
        '$   GRID       10478           35.46     -6.      0.',
        '$ 5971 10479',
        '$   GRID        5971           34.92     -6.      0.',
        '$   GRID       10479           34.92     -6.      0.',
        '$ 6003 10477',
        '$   GRID        6003             36.     -6.      0.',
        '$   GRID       10477             36.     -6.      0.',
        'GRID        5971           34.92     -6.      0.',
        'GRID        5972           34.92-5.73333      0.',
        'GRID        5973           34.92-5.46667      0.',
        'GRID        5987           35.46     -6.      0.',
        'GRID        5988           35.46-5.73333      0.',
        'GRID        5989           35.46-5.46667      0.',
        'GRID        6003             36.     -6.      0.',
        'GRID        6004             36.-5.73333      0.',
        'GRID        6005             36.-5.46667      0.',
        'GRID       10476             36.     -6.    -1.5',
        'GRID       10477             36.     -6.      0.',
        'GRID       10478           35.46     -6.      0.',
        'GRID       10479           34.92     -6.      0.',
        'GRID       10561           34.92     -6.    -.54',
        '$ELEMENTS_WITH_PROPERTIES',
        'PSHELL         1       1      .1',
        'CQUAD4      5471       1    5971    5987    5988    5972',
        'CQUAD4      5472       1    5972    5988    5989    5973',
        'CQUAD4      5486       1    5987    6003    6004    5988',
        'CQUAD4      5487       1    5988    6004    6005    5989',
        'PSHELL        11       1      .1',
        'CTRIA3      9429      11   10561   10476   10478',
        'CTRIA3      9439      11   10478   10479   10561',
        'CTRIA3      9466      11   10476   10477   10478',
        '$MATERIALS',
        'MAT1           1      3.              .3',
    ]
    bdf_filename = 'nonunique2.bdf'
    f = open(bdf_filename, 'wb')
    f.write('\n'.join(lines))
    f.close()
    bdf_equivalence_nodes('nonunique2.bdf', 'unique2.bdf', 0.01)

def eq1():
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
    eq1()
    eq2()
