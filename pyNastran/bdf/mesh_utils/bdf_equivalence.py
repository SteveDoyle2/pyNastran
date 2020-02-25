"""
defines:
    model = bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                                  renumber_nodes=False, neq_max=4, xref=True,
                                  node_set=None,
                                  size=8, is_double=False,
                                  remove_collapsed_elements=False,
                                  avoid_collapsed_elements=False,
                                  crash_on_collapse=False, log=None, debug=True)

"""
from itertools import combinations
from typing import List, Tuple, Dict, Union, Optional, Any
import numpy as np
from numpy import (array, unique, arange, searchsorted,
                   setdiff1d, intersect1d, asarray)
from numpy.linalg import norm  # type: ignore
import scipy
import scipy.spatial

from pyNastran.nptyping import NDArrayNint, NDArrayN3float
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.mesh_utils.internal_utils import get_bdf_model


def bdf_equivalence_nodes(bdf_filename: str, bdf_filename_out: str, tol: float,
                          renumber_nodes: bool=False, neq_max: int=4, xref: bool=True,
                          node_set: Union[List[int], NDArrayNint]=None,
                          size: int=8, is_double: bool=False,
                          remove_collapsed_elements: bool=False,
                          avoid_collapsed_elements: bool=False,
                          crash_on_collapse: bool=False,
                          log=None, debug: bool=True, method: str='new'):
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
    xref : bool
        does the model need to be cross_referenced
        (default=True; only applies to model option)
    node_set : List[int] / (n, ) ndarray; default=None
        the list/array of nodes to consider
        (not supported with renumber_nodes=True)
    size : int; {8, 16}; default=8
        the bdf write precision
    is_double : bool; default=False
        the field precision to write
    remove_collapsed_elements : bool; default=False (unsupported)
        True  : 1D/2D/3D elements will not be collapsed;
                CELASx/CDAMP/MPC/etc. are not considered
        False : no elements will be removed
    avoid_collapsed_elements : bool; default=False (unsupported)
        True  : only collapses that don't break 1D/2D/3D elements will be considered;
                CELASx/CDAMP/MPC/etc. are considered
        False : element can be collapsed
    crash_on_collapse : bool; default=False
        stop if nodes have been collapsed
           False: blindly move on
           True: rereads the BDF which catches doubled nodes (temporary);
                 in the future collapse=True won't need to double read;
                 an alternative is to do Patran's method of avoiding collapse)
    debug : bool
        bdf debugging
    method: str; default='new'
        'new': doesn't require neq_max; new in v1.3
        'old': use neq_max; used in v1.2
    log : logger(); default=None
        bdf logging

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
    if not isinstance(tol, float):
        tol = float(tol)
    nodes_xyz, model, nids, inew = _eq_nodes_setup(
        bdf_filename, tol, renumber_nodes=renumber_nodes,
        xref=xref, node_set=node_set, log=log, debug=debug)

    nid_pairs = _nodes_xyz_nids_to_nid_pairs(
        nodes_xyz, nids, tol, log, inew,
        node_set=node_set, neq_max=neq_max, method=method, debug=debug)
    _eq_nodes_final(nid_pairs, model, tol, node_set=node_set, debug=debug)

    if bdf_filename_out is not None:
        model.write_bdf(bdf_filename_out, size=size, is_double=is_double)
    if crash_on_collapse:
        # lazy way to make sure there aren't any collapsed nodes
        model2 = BDF(log=log, debug=debug)
        model2.read_bdf(bdf_filename_out)
    return model


def _eq_nodes_setup(bdf_filename, unused_tol,
                    renumber_nodes=False, xref=True,
                    node_set=None, log=None, debug=True):
    """helper function for ``bdf_equivalence_nodes``"""
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

    model = get_bdf_model(bdf_filename, xref=xref, log=log, debug=debug)

    # quads / tris
    #nids_quads = []
    #eids_quads = []
    #nids_tris = []
    #eids_tris = []

    # map the node ids to the slot in the nids array
    renumber_nodes = False
    if node_set is not None:
        nids, all_nids, unused_nid_map = _eq_nodes_setup_node_set(
            model, node_set, renumber_nodes=renumber_nodes)
    else:
        nids, all_nids, unused_nid_map = _eq_nodes_setup_node(
            model, renumber_nodes=renumber_nodes)

    nodes_xyz = _get_xyz_cid0(model, nids)
    inew = _check_for_referenced_nodes(model, node_set, nids, all_nids, nodes_xyz)

    #assert np.array_equal(nids[inew], nids_new), 'some nodes are not defined'
    return nodes_xyz, model, nids, inew

def _get_xyz_cid0(model, nids):
    """gets xyz_cid0"""
    coord_ids = model.coord_ids
    needs_get_position = (coord_ids == [0])

    if needs_get_position:
        nodes_xyz = array([model.nodes[nid].get_position()
                           for nid in nids], dtype='float32')
    else:
        nodes_xyz = array([model.nodes[nid].xyz
                           for nid in nids], dtype='float32')
    return nodes_xyz

def _eq_nodes_setup_node_set(model: BDF, node_set, renumber_nodes: bool=False,
                             ) -> Tuple[NDArrayNint, NDArrayNint, Dict[int, int]]:
    """helper function for ``_eq_nodes_setup``"""
    inode = 0
    nid_map = {}
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
    #nids = array([node.nid for nid, node in sorted(model.nodes.items())
                    #if nid in node_set], dtype='int32')
    return nids, all_nids, nid_map

def _eq_nodes_setup_node(model: BDF, renumber_nodes: bool=False,
                         ) -> Tuple[NDArrayNint, NDArrayNint, Dict[int, int]]:
    """helper function for ``_eq_nodes_setup``"""
    inode = 0
    nid_map = {}
    if renumber_nodes:
        for nid, node in sorted(model.nodes.items()):
            node.nid = inode + 1
            nid_map[inode] = nid
            inode += 1
        nnodes = len(model.nodes)
        nids = arange(1, inode + 1, dtype='int32')
        assert nids[-1] == nnodes
    else:
        for nid, node in sorted(model.nodes.items()):
            nid_map[inode] = nid
            inode += 1
        nids = array([node.nid for nid, node in sorted(model.nodes.items())], dtype='int32')
    all_nids = nids
    return nids, all_nids, nid_map

def _check_for_referenced_nodes(model: BDF,
                                node_set: NDArrayNint,
                                nids: NDArrayNint,
                                all_nids: NDArrayNint,
                                nodes_xyz: NDArrayN3float) -> Optional[NDArrayNint]:
    """helper function for ``_eq_nodes_setup``"""
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
        spoint_epoint_nid_set = set()
        for unused_eid, element in sorted(model.elements.items()):
            spoint_epoint_nid_set.update(element.node_ids)
        for unused_eid, element in sorted(model.masses.items()):
            spoint_epoint_nid_set.update(element.node_ids)

        nids_new = spoint_epoint_nid_set - set(model.spoints) - set(model.epoints)

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
        inew = None
    #assert np.array_equal(nids[inew], nids_new), 'some nodes are not defined'
    return inew

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


def _eq_nodes_final(nid_pairs, model: BDF, tol: float,
                    node_set=None, debug: bool=False) -> None:
    """apply nodal equivalencing to model"""
    for (nid1, nid2) in nid_pairs:
        node1 = model.nodes[nid1]
        node2 = model.nodes[nid2]

        # TODO: doesn't use get position...
        distance = norm(node1.get_position() - node2.get_position())

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

def _nodes_xyz_nids_to_nid_pairs(nodes_xyz: NDArrayN3float,
                                 nids: NDArrayNint,
                                 tol: float,
                                 log,
                                 inew: NDArrayNint,
                                 node_set=None, neq_max: int=4, method: str='new',
                                 debug: bool=False) -> None:
    """helper for equivalencing"""
    unused_kdt, nid_pairs = _eq_nodes_build_tree(
        nodes_xyz, nids, tol, log,
        inew=inew, node_set=node_set,
        neq_max=neq_max, method=method, debug=debug)
    return nid_pairs

def _nodes_xyz_nids_to_nid_pairs_new(kdt, nids, node_set, tol: float):
    """
    helper function for `bdf_equivalence_nodes`
    """
    ieq3 = kdt.query_ball_tree(kdt, tol)
    nid_pairs = []
    for pair in ieq3:
        if len(pair) == 1:
            continue
        # the combinations should be paired with 2 nodes in each group
        for inid1, inid2 in combinations(pair, 2):
            nid1 = nids[inid1]
            nid2 = nids[inid2]
            #if nid1 == nid2:
                #continue
            if node_set is not None:
                if nid1 not in node_set and nid2 not in node_set:
                    continue
            if (nid1, nid2) in nid_pairs:
                continue
            nid_pairs.append((nid1, nid2))
    return nid_pairs

def _eq_nodes_build_tree(nodes_xyz, nids, tol, log,
                         inew=None, node_set=None, neq_max=4, method='new', msg='',
                         debug=False) -> Tuple[Any, Any, Any]:
    """
    helper function for `bdf_equivalence_nodes`

    Parameters
    ----------
    nodes_xyz : (nnodes, 3) float ndarray
        the xyzs to equivalence
    nids : (nnodes,) int ndarray
        the node ids
    tol : float
        the spherical equivalence tolerance
    inew : int ndarray; default=None -> slice(None)
        a slice on nodes_xyz to exclude some nodes from the equivalencing
    node_set : List[int] / (n, ) ndarray; default=None
        the list/array of nodes to consider
    neq_max : int; default=4
        the number of nodes to consider for equivalencing
    msg : str; default=''
        custom message used for errors

    Returns
    -------
    kdt : cKDTree()
        the kdtree object
    ieq : (nnodes, neq) int ndarray
        ???
    slots : (nnodes, neq) int ndarray
        ???

    """
    nnodes = len(nids)
    if inew is None:
        inew = slice(None)

    assert isinstance(tol, float), 'tol=%r' % tol
    kdt = _get_tree(nodes_xyz, msg=msg)

    is_not_node_set = inew is None or inew == slice(None)

    # check the closest 10 nodes for equality
    if method == 'new' and is_not_node_set:
        kdt, nid_pairs = _eq_nodes_build_tree_new(
            kdt, nodes_xyz,
            nids, nnodes, is_not_node_set,
            tol, log,
            inew=inew, node_set=node_set, neq_max=neq_max, msg=msg,
            debug=debug)
    else:
        if method == 'new':
            log.warning(f'setting method to "old" because node_set is specified')
        deq, ieq = kdt.query(nodes_xyz[inew, :], k=neq_max, distance_upper_bound=tol)
        if node_set is not None:
            assert len(deq) == len(nids)

        # get the ids of the duplicate nodes
        slots = np.where(ieq[:, :] < nnodes)
        nid_pairs = _eq_nodes_find_pairs(nids, slots, ieq, node_set=node_set)
    assert isinstance(nid_pairs, list), nid_pairs
    return kdt, nid_pairs

def _eq_nodes_build_tree_new(kdt,
                             nodes_xyz: NDArrayN3float,
                             nids, nnodes, is_not_node_set: bool,
                             tol: float,
                             log,
                             inew=None, node_set=None, neq_max: int=4, msg: str='',
                             debug: float=False) -> Tuple[Any, List[Tuple[int, int]]]:
    deq, ieq = kdt.query(nodes_xyz[inew, :], k=neq_max, distance_upper_bound=tol)
    slots = np.where(ieq[:, :] < nnodes)
    nid_pairs_expected = _eq_nodes_find_pairs(nids, slots, ieq, node_set=node_set)
    if is_not_node_set:
        nid_pairs = _nodes_xyz_nids_to_nid_pairs_new(kdt, nids, node_set, tol)
    else:
        raise NotImplementedError(f'node_set = {node_set}')

    snid_pairs = set(nid_pairs)
    snid_pairs_expected = set(nid_pairs_expected)
    diff_bad = snid_pairs - snid_pairs_expected
    diff_missed = snid_pairs - snid_pairs_expected
    if debug and len(diff_bad) or len(diff_missed):  # pragma: no cover
        #log.warning(f'nid_pairs          = {nid_pairs}')
        #log.warning(f'nid_pairs_expected = {nid_pairs_expected}')
        log.warning(f'diff_bad = {diff_bad}')
        log.warning(f'diff_missed = {diff_missed}')

    return kdt, nid_pairs

def _get_tree(nodes_xyz, msg=''):
    """gets the kdtree"""
    assert isinstance(nodes_xyz, np.ndarray), type(nodes_xyz)
    assert nodes_xyz.shape[0] > 0, 'nnodes=0%s' % msg

    # build the kdtree
    try:
        kdt = scipy.spatial.cKDTree(nodes_xyz)
    except RuntimeError:
        print(nodes_xyz)
        raise RuntimeError(nodes_xyz)
    return kdt
