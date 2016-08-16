from __future__ import print_function
#from collections import defaultdict
#from functools import reduce

from six import iteritems, string_types, PY2
#from six.moves import zip, range


import numpy as np
from numpy import (array, unique, arange, searchsorted,
                   setdiff1d, intersect1d, asarray)

from numpy.linalg import norm
import scipy

from pyNastran.utils import integer_types
from pyNastran.bdf.bdf import BDF


def bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                          renumber_nodes=False, neq_max=4, xref=True,
                          node_set=None,
                          size=8, is_double=False,
                          remove_collapsed_elements=False,
                          avoid_collapsed_elements=False,
                          crash_on_collapse=False, debug=True):
    """
    Equivalences nodes; keeps the lower node id; creates two nodes with the same

    Parameters
    ----------
    bdf_filename : str / BDF
        str : bdf file path
        BDF : a BDF model that is fully valid (see xref)
    bdf_filename_out : str
        a bdf_filename to write
    tol : float
        the spherical tolerance
    renumber_nodes : bool
        should the nodes be renumbered (default=False)
    neq_max : int
        the number of "close" points (default=4)
    xref bool: bool
        does the model need to be cross_referenced
        (default=True; only applies to model option)
    node_set : List[int] / (n, ) ndarray
        the list/array of nodes to consider (not supported with renumber_nodes=True)
    size : int; {8, 16}; default=8
        the bdf write precision
    is_double : bool; default=False
        the field precision to write
    crash_on_collapse : bool; default=False
        stop if nodes have been collapsed
           False: blindly move on
           True: rereads the BDF which catches doubled nodes (temporary);
                 in the future collapse=True won't need to double read;
                 an alternative is to do Patran's method of avoiding collapse)
    remove_collapsed_elements : bool; default=False (unsupported)
        True  : 1D/2D/3D elements will not be collapsed;
                CELASx/CDAMP/MPC/etc. are not considered
        False : no elements will be removed
    avoid_collapsed_elements : bool; default=False (unsupported)
        True  : only collapses that don't break 1D/2D/3D elements will be considered;
                CELASx/CDAMP/MPC/etc. are considered
        False : element can be collapsed

    Returns
    -------
    model : BDF()
        The BDF model corresponding to bdf_filename_out

    .. warning:: I doubt SPOINTs/EPOINTs work correctly
    .. warning:: xref not fully implemented (assumes cid=0)

    .. todo:: node_set stil does work on the all the nodes in the big
               kdtree loop, which is very inefficient
    .. todo:: remove_collapsed_elements is not supported
    .. todo:: avoid_collapsed_elements is not supported
    """
    nodes_xyz, model, nids, inew = _eq_nodes_setup(
        bdf_filename, tol, renumber_nodes=renumber_nodes,
        xref=xref, node_set=node_set, debug=debug)
    ieq, slots = _eq_nodes_build_tree(nodes_xyz, nids, tol,
                                      inew=inew, node_set=node_set,
                                      neq_max=neq_max)[1:]

    nid_pairs = _eq_nodes_find_pairs(nids, slots, ieq, node_set=node_set)
    _eq_nodes_final(nid_pairs, model, tol, node_set=node_set)

    if bdf_filename_out is not None:
        model.write_bdf(bdf_filename_out, size=size, is_double=is_double)
    if crash_on_collapse:
        # lazy way to make sure there aren't any collapsed nodes
        model2 = BDF(debug=debug)
        model2.read_bdf(bdf_filename_out)
    return model


def _eq_nodes_setup(bdf_filename, tol,
                    renumber_nodes=False, xref=True,
                    node_set=None, debug=True):
    """helper function for `bdf_equivalence_nodes`"""
    assert isinstance(tol, float), tol
    if node_set is not None:
        if renumber_nodes:
            raise NotImplementedError('node_set is not None & renumber_nodes=True')

        #print(type(node_set))
        #print('*node_set', node_set)
        assert len(node_set) > 0, node_set
        if isinstance(node_set, set):
            node_set = asarray(list(node_set), dtype='int32')
        else:
            node_set = asarray(node_set, dtype='int32')

    if isinstance(bdf_filename, string_types):
        xref = True
        model = BDF(debug=debug)
        model.read_bdf(bdf_filename, xref=True)
    else:
        model = bdf_filename
        model.cross_reference(xref=xref)

    coord_ids = model.coord_ids
    needs_get_position = True if coord_ids == [0] else False

    # quads / tris
    #nids_quads = []
    #eids_quads = []
    #nids_tris = []
    #eids_tris = []

    # map the node ids to the slot in the nids array
    renumber_nodes = False

    inode = 0
    nid_map = {}
    if node_set is not None:
        if PY2:
            all_nids = array(model.nodes.keys(), dtype='int32')
        else:
            all_nids = array(list(model.nodes.keys()), dtype='int32')

        # B - A
        # these are all the nodes that are requested from node_set that are missing
        #   thus len(diff_nodes) == 0
        diff_nodes = setdiff1d(node_set, all_nids)
        if len(diff_nodes) != 0:
            msg = ('The following nodes cannot be found, but are included'
                   ' in the reduced set; nids=%s' % diff_nodes)
            raise RuntimeError(msg)

        # A & B
        # the nodes to analyze are the union of all the nodes and the desired set
        # which is basically the same as:
        #   nids = unique(node_set)
        nids = intersect1d(all_nids, node_set, assume_unique=True)  # the new values

        if renumber_nodes:
            raise NotImplementedError('node_set is not None & renumber_nodes=True')
        else:
            for nid in all_nids:
                nid_map[inode] = nid
                inode += 1
        #nids = array([node.nid for nid, node in sorted(iteritems(model.nodes))
                        #if nid in node_set], dtype='int32')

    else:
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
        all_nids = nids

    if needs_get_position:
        nodes_xyz = array([model.nodes[nid].get_position()
                           for nid in nids], dtype='float32')
    else:
        nodes_xyz = array([model.nodes[nid].xyz
                           for nid in nids], dtype='float32')

    if node_set is not None:
        assert nodes_xyz.shape[0] == len(nids)

    if 0:
        # I forget entirely what this block of code is for, but my general
        # recollection was that it checked that all the nodes that were
        # referenced were included in the nids list.  I'd rather break that
        # check in order to support nodes_set.
        #
        # It's also possible that it's here, so you only consider nodes that
        # are associated...

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

        missing_nids = list(set(nids_new) - set(all_nids))
        if missing_nids:
            missing_nids.sort()
            msg = 'There are missing nodes...\n'  # TODO: in what???
            msg = 'missing nids=%s' % str(missing_nids)
            raise RuntimeError(msg)

        # get the node_id mapping for the kdtree
        inew = searchsorted(nids, nids_new, side='left')
        # print('nids_new =', nids_new)
    else:
        inew = slice(None)
    #assert np.array_equal(nids[inew], nids_new), 'some nodes are not defined'
    return nodes_xyz, model, nids, inew


def _eq_nodes_find_pairs(nids, slots, ieq, node_set=None):
    """helper function for `bdf_equivalence_nodes`"""
    irows, icols = slots
    #replacer = unique(ieq[slots])  ## TODO: turn this back on?

    #skip_nodes = []
    nid_pairs = []
    for (irow, icol) in zip(irows, icols):
        inid2 = ieq[irow, icol]
        nid1 = nids[irow]
        nid2 = nids[inid2]
        if nid1 == nid2:
            continue
        if node_set is not None:
            if nid1 not in node_set and nid2 not in node_set:
                continue
        nid_pairs.append((nid1, nid2))
    return nid_pairs


def _eq_nodes_final(nid_pairs, model, tol, node_set=None):
    """apply nodal equivalencing to model"""
    for (nid1, nid2) in nid_pairs:
        node1 = model.nodes[nid1]
        node2 = model.nodes[nid2]

        # TODO: doesn't use get position...
        distance = norm(node1.xyz - node2.xyz)

        #print('  irow=%s->n1=%s icol=%s->n2=%s' % (irow, nid1, icol, nid2))
        if distance > tol:
            #print('  *n1=%-4s xyz=%s\n  *n2=%-4s xyz=%s\n  *distance=%s\n' % (
            #    nid1, list_print(node1.xyz),
            #    nid2, list_print(node2.xyz),
            #    distance))
            continue

        if node_set is not None:
            assert nid1 in node_set, 'nid1=%s node_set=%s' % (nid1, node_set)
            assert nid2 in node_set, 'nid2=%s node_set=%s' % (nid2, node_set)
            #print('  n1=%-4s xyz=%s\n  n2=%-4s xyz=%s\n  distance=%s\n' % (
                #nid1, str(node1.xyz),
                #nid2, str(node2.xyz),
                #distance))

        # if hasattr(node2, 'new_node_id'):

        # else:
        node2.nid = node1.nid
        node2.xyz = node1.xyz
        node2.cp = node1.cp
        assert node2.cd == node1.cd
        assert node2.ps == node1.ps
        assert node2.seid == node1.seid
        # node2.new_nid = node1.nid
        #skip_nodes.append(nid2)
    return

def _eq_nodes_build_tree(nodes_xyz, nids, tol, inew=None, node_set=None, neq_max=4):
    """helper function for `bdf_equivalence_nodes`"""
    kdt = _get_tree(nodes_xyz)

    # check the closest 10 nodes for equality
    deq, ieq = kdt.query(nodes_xyz[inew, :], k=neq_max, distance_upper_bound=tol)

    if node_set is not None:
        assert len(deq) == len(nids)
    nnodes = len(nids)

    # get the ids of the duplicate nodes
    slots = np.where(ieq[:, :] < nnodes)
    return kdt, ieq, slots


def _get_tree(nodes_xyz):
    """gets the kdtree"""
    assert isinstance(nodes_xyz, np.ndarray), type(nodes_xyz)

    # build the kdtree
    if scipy.__version__ < '0.18.1':
        try:
            kdt = scipy.spatial.KDTree(nodes_xyz)
        except RuntimeError:
            print(nodes_xyz)
            raise RuntimeError(nodes_xyz)
    else:
        try:
            kdt = scipy.spatial.cKDTree(nodes_xyz)
        except RuntimeError:
            print(nodes_xyz)
            raise RuntimeError(nodes_xyz)
    return kdt
