"""
defines:
 - convert_bad_quads_to_tris(model, eids_to_check=None, xyz_cid0=None,
                             min_edge_length=0.0)

"""
import numpy as np
from pyNastran.bdf.cards.elements.shell import CTRIA3


def convert_bad_quads_to_tris(model, eids_to_check=None, xyz_cid0=None, min_edge_length=0.0):
    """
    A standard quad is a nice rectangle.  If an edge is collapsed, it's a triangle.
    Change the element type.

    Parameters
    ----------
    model : BDF()
        a BDF model that has not had it's properties/load xref'd, but is valid
        such that it could
    eids : list; (default=None -> all CQUAD4s)
        the subset of element ids to check
    xyz_cid0 : (n, 3) ndarray
        nodes in cid=0
    min_edge_length : float; default=0.0
        what is classified as "short"

    .. warning::  Don't cross reference properties/loads

    .. todo::  check for bad xref

    """
    out = model.get_card_ids_by_card_types('CQUAD4')
    cquad4s = out['CQUAD4']
    if eids_to_check is None:
        cquad4s_to_check = cquad4s
    else:
        cquad4s_to_check = list(set(eids_to_check).intersection(set(cquad4s)))

    elements = model.elements
    unused_eids_to_remove = []

    if xyz_cid0 is None:
        xyz_cid0 = model.get_xyz_in_coord(cid=0)
    nid_cd = np.array([[nid, node.Cd()] for nid, node in sorted(model.nodes.items())])
    all_nids = nid_cd[:, 0]

    if len(cquad4s_to_check) == 0:
        model.log.warning('no quads in the model...')
        return

    for eid in sorted(cquad4s_to_check):
        elem = elements[eid]
        nids = elem.node_ids
        nids_long = nids + nids

        nids_to_remove = []
        for inode in range(4):
            nid0 = nids_long[inode]
            nid1 = nids_long[inode + 1]
            edge = [nid0, nid1]
            i = np.searchsorted(all_nids, edge)
            xyz = xyz_cid0[i, :]
            edge_length = np.linalg.norm(xyz[0, :] - xyz[1, :])
            if edge_length <= min_edge_length:
                nids_to_remove.append(nid1)

        if len(nids_to_remove) == 0:
            continue

        for nid in nids_to_remove:
            nids.remove(nid)

        if len(nids) < 3:
            model.log.debug('found eid=%s is a line/point...removing' % eid)
            etype = elem.type
            model.card_count[etype] -= 1
            del model.elements[eid]
            continue

        nids_to_remove = []
        nids_long = nids + nids
        for inode in range(3):
            nid0 = nids_long[inode]
            nid1 = nids_long[inode + 1]
            edge = [nid0, nid1]
            i = np.searchsorted(all_nids, edge)
            xyz = xyz_cid0[i, :]
            edge_length = np.linalg.norm(xyz[0, :] - xyz[1, :])
            if edge_length <= min_edge_length:
                nids_to_remove.append(nid1)

        for nid in nids_to_remove:
            nids.remove(nid)

        if len(nids) < 3:
            model.log.debug('found eid=%s is a line/point...removing' % eid)
            etype = elem.type
            model.card_count[etype] -= 1
            del model.elements[eid]
            continue

        model.log.debug('found eid=%s is a triangle...replacing CQUAD4 with CTRIA3' % eid)
        #print('  nids2=%s' % (nids))
        assert elem.tflag == 0, elem
        assert elem.T1 is None, elem.T1
        assert elem.T2 is None, elem.T2
        assert elem.T3 is None, elem.T3
        assert elem.T4 is None, elem.T4
        elem2 = CTRIA3(eid, elem.Pid(), nids, elem.zoffset,
                       theta_mcid=elem.theta_mcid, tflag=0, T1=None, T2=None, T3=None,
                       comment='$ was a CQUAD4\n')
        model.increase_card_count('CTRIA3')
        del model.elements[eid]
        model.elements[eid] = elem2
