"""
defines:
  - extract_bodies(bdf_filename)
"""
from __future__ import print_function
from collections import defaultdict
from six import iteritems, iterkeys
import numpy as np
from pyNastran.bdf.bdf import BDF, read_bdf

def extract_bodies(bdf_filename, mpc_id=0):
    """
    Finds the isolated bodies

    Paramters
    ---------
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
    assert len(model.nodes) > 0

    nid_to_eid_map = defaultdict(list)
    eid_to_nid_map = defaultdict(list)
    for eid, _elem in iteritems(model.elements):
        #print(elem)
        eid_to_nid_map[eid] = model.nodes
        for nid in model.nodes:
            nid_to_eid_map[nid].append(eid)

    rigid_offset = max(model.elements)
    for eid, _elem in iteritems(model.rigid_elements):
        #print(elem)
        eid += rigid_offset
        assert eid not in model.rigid_elements, 'eid=%s cannot be used twice' % eid
        eid_to_nid_map[eid] = model.nodes
        for nid in model.nodes:
            nid_to_eid_map[nid].append(eid)
    #mpc_offset = rigid_offset + max(model.rigid_elements)

    nids_used = set([])
    eids_used = set([])
    ibody = 0
    nids_to_check = set([next(iterkeys(model.nodes))])
    all_nids_to_check = set(list(model.nodes.keys()))

    body_eids = {ibody : set([])}
    while all_nids_to_check:
        #print(all_nids_to_check)
        while nids_to_check:
            #if len(nids_to_check) < 10:
                #print('nids_to_check =', nids_to_check)
            nid = nids_to_check.pop()
            #print('nid = ', nid)
            #print('all_nids_to_check =', all_nids_to_check)
            if nid not in all_nids_to_check:
                msg = 'nids_to_check = %s' % nids_to_check
                raise RuntimeError(msg)
            #print('nids_used         =', nids_used)
            #print('all_nids_to_check =', all_nids_to_check)

            nids_used.add(nid)
            #print('  adding nid=%s' % nid)
            all_nids_to_check.remove(nid)
            #print('  nids_used         =', nids_used)
            #print('  all_nids_to_check =', all_nids_to_check)

            eids = nid_to_eid_map[nid]
            #print('  eids = %s' % eids)
            for eidi in eids:
                if eidi in eids_used:
                    continue
                eids_used.add(eidi)
                #try:
                    #elem = model.elements[eidi]
                #except KeyError:
                    #elem = model.rigid_elements[eidi - rigid_offset]
                #print(elem)
                #print('adding eidi=%s' % eidi)
                body_eids[ibody].add(eidi)
                nidsi = eid_to_nid_map[eid]
                for nidi in nidsi:
                    if nidi not in nids_used:
                        nids_to_check.add(nidi)
                        #print('  adding nid=%s' % nidi)
                del nidsi # , elem
            #print('next eid\n')
            del nid, eids
        ibody += 1
        body_eids[ibody] = set([])
        print('--------------------------------------')
    if len(body_eids[ibody]) == 0:
        del body_eids[ibody]

    body_eids2 = {}
    for ibody, body in sorted(iteritems(body_eids)):
        abody = np.unique(np.array(list(body), dtype='int32'))
        ielem = np.where(abody <= rigid_offset)
        irigid = np.where(abody > rigid_offset)
        body_eids2[ibody] = [abody[ielem], abody[irigid]]
    #print('body_eids = %s' % body_eids2)
    nbodies = len(body_eids2)
    if nbodies > 1:
        print('nbodies = %s' % nbodies)
    return body_eids2


