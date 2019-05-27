"""
defines:
  - extract_bodies(bdf_filename)

"""
from collections import defaultdict
import numpy as np
from pyNastran.bdf.bdf import BDF, read_bdf, print_card_16

def extract_bodies(bdf_filename, mpc_id=0):
    """
    Finds the isolated bodies

    Parameters
    ----------
    bdf_filename : str/BDF
        str : the path the the *.bdf file
        BDF : a BDF() boject
    mpc_id : int; default=0
        0 : consider all MPCs
        >0 : use this MPC set
        not supported

    Considers:
     - elements
     - rigid_elements

    Doesn't consider:
      - elements_mass
      - MPC
      - MPCADD
      - DMIx

    Doesn't support:
      - xref
      - duplicate element ids
      - large values

    """
    if isinstance(bdf_filename, BDF):
        model = bdf_filename
    else:
        model = read_bdf(bdf_filename, xref=False)

    debug = False
    nnodes = len(model.nodes)
    nspoints = 0
    nepoints = 0
    if model.spoints:
        nspoints = len(model.spoints)
    if model.epoints:
        nepoints = len(model.epoints)
    npoints = nnodes + nspoints + nepoints
    nelements = len(model.elements) + len(model.rigid_elements) + len(model.masses)
    if npoints == 0 or nelements == 0:
        return {}

    nid_to_eid_map = defaultdict(list)
    eid_to_nid_map = defaultdict(list)
    for eid, elem in model.elements.items():
        if debug:  # pragma: no cover
            print(print_card_16(elem.repr_fields()))
        node_ids = elem.node_ids
        eid_to_nid_map[eid] = node_ids
        for nid in node_ids:
            if nid is None:
                continue
            nid_to_eid_map[nid].append(eid)

    rigid_offset = 0
    if len(model.elements):
        rigid_offset = max(model.elements)
    for eid, elem in model.rigid_elements.items():
        if debug:  # pragma: no cover
            print(print_card_16(elem.repr_fields()))
        eid += rigid_offset
        assert eid not in model.elements, 'eid=%s cannot be used twice' % eid
        #node_ids = elem.node_ids
        node_ids = elem.independent_nodes + elem.dependent_nodes
        eid_to_nid_map[eid] = node_ids
        for nid in node_ids:
            if nid is None:
                raise RuntimeError(elem)
            nid_to_eid_map[nid].append(eid)
    #mpc_offset = rigid_offset + max(model.rigid_elements)

    if len(nid_to_eid_map) == 0:
        raise RuntimeError(model.get_bdf_stats())
        #return {}

    nids_used = set()
    eids_used = set()
    ibody = 0
    keys = list(nid_to_eid_map.keys())
    all_nids_to_check = set(keys)
    key0 = keys[0]
    nids_to_check = set([key0])
    if debug:  # pragma: no cover
        print('all_nids_to_check= ', all_nids_to_check)
    body_eids = {ibody : set()}
    nbodies_max = 3
    while all_nids_to_check:
        if debug:  # pragma: no cover
            print(all_nids_to_check)
        while nids_to_check:
            #if len(nids_to_check) < 10:
                #print('nids_to_check =', nids_to_check)
            nid = nids_to_check.pop()
            if debug:  # pragma: no cover
                print('nid = ', nid)
                print('all_nids_to_check =', all_nids_to_check)
            if nid not in all_nids_to_check:
                msg = 'nids_to_check = %s' % nids_to_check
                raise RuntimeError(msg)
            if debug:  # pragma: no cover
                print('nids_used         =', nids_used)
                print('all_nids_to_check =', all_nids_to_check)

            nids_used.add(nid)
            if debug:  # pragma: no cover
                print('  adding nid=%s' % nid)
            all_nids_to_check.remove(nid)
            if debug:  # pragma: no cover
                print('  nids_used         =', nids_used)
                print('  all_nids_to_check =', all_nids_to_check)

            eids = nid_to_eid_map[nid]
            if debug:  # pragma: no cover
                print('  eids = %s' % eids)
            for eidi in eids:
                if eidi in eids_used:
                    continue
                eids_used.add(eidi)
                try:
                    elem = model.elements[eidi]
                except KeyError:
                    elem = model.rigid_elements[eidi - rigid_offset]
                if debug:  # pragma: no cover
                    print(print_card_16(elem.repr_fields()))
                    print('adding eidi=%s' % eidi)
                body_eids[ibody].add(eidi)
                nidsi = eid_to_nid_map[eidi]
                for nidi in nidsi:
                    if nidi not in nids_used:
                        nids_to_check.add(nidi)
                        if debug:  # pragma: no cover
                            print('  adding nid=%s' % nidi)
                del nidsi # , elem
            if debug:  # pragma: no cover
                print('next eid\n')
                print('---------------------------')
            del nid, eids

        if len(body_eids[ibody]) == 0:
            return body_eids
            # unused node
            #try:
                #nid = nids_to_check.pop()
            #except KeyError:
                #break
            #continue
            #msg = 'cannot find a new body...nbodies=%s\nelements:' % ibody
            #for nid, eids in sorted(nid_to_eid_map.items()):
                #msg += '  nid=%r eids=%s\n' % (nid, eids)
                #for eid in eids:
                    #try:
                        #elem = model.elements[eid]
                    #except KeyError:
                        #elem = model.rigid_elements[eid - rigid_offset]
                    #msg += print_card_16(elem.repr_fields())

                    #msg += 'eid=%s used=%s\n\n' % (eid, eid in eids_used)
                #msg += '----------------\n'
            #msg += ''
            #raise RuntimeError(msg + model.get_bdf_stats())
        ibody += 1
        if ibody > nbodies_max:
            raise RuntimeError('Too many bodies...\n' + model.get_bdf_stats())
        body_eids[ibody] = set()
        #print('--------------------------------------')
    if len(body_eids[ibody]) == 0:
        del body_eids[ibody]

    body_eids2 = {}
    for ibody, body in sorted(body_eids.items()):
        abody = np.unique(np.array(list(body), dtype='int64'))
        ielem = np.where(abody <= rigid_offset)
        irigid = np.where(abody > rigid_offset)
        #ielem = np.asarray(ielem, dtype='int32')
        #ielem = np.asarray(ielem, dtype='int32')
        body_eids2[ibody] = [
            np.asarray(abody[ielem], dtype='int32'),
            np.asarray(abody[irigid], dtype='int32')
        ]
    #print('body_eids = %s' % body_eids2)
    nbodies = len(body_eids2)
    if nbodies > 1:
        print('nbodies = %i' % nbodies)
    return body_eids2
