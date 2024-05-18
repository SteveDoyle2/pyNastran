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
from __future__ import annotations
from itertools import combinations
from typing import Union, Optional, TYPE_CHECKING
import numpy as np
from numpy.linalg import norm  # type: ignore
import scipy
from scipy.spatial import KDTree

from pyNastran.nptyping_interface import NDArrayNint, NDArrayN3float
from pyNastran.femutils.utils import unique2d
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.mesh_utils.internal_utils import get_bdf_model
if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger
    from pyNastran.bdf.bdf import GRID


def bdf_equivalence_nodes(bdf_filename: str,
                          bdf_filename_out: Optional[str],
                          tol: float,
                          renumber_nodes: bool=False, neq_max: int=4, xref: bool=True,
                          node_set: Optional[Union[list[int], NDArrayNint]]=None,
                          size: int=8, is_double: bool=False,
                          remove_collapsed_elements: bool=False,
                          avoid_collapsed_elements: bool=False,
                          crash_on_collapse: bool=False,
                          log: Optional[SimpleLogger]=None,
                          debug: bool=True, method: str='new') -> BDF:
    """
    Equivalences nodes; keeps the lower node id; creates two nodes with the same

    Parameters
    ----------
    bdf_filename : str / BDF
        str : bdf file path
        BDF : a BDF model that is fully valid (see xref)
    bdf_filename_out : str
        str: a bdf_filename to write
        None: don't write the deck
    tol : float
        the spherical tolerance
    renumber_nodes : bool
        should the nodes be renumbered (default=False)
    neq_max : int
        the number of "close" points (default=4)
    xref : bool
        does the model need to be cross_referenced
        (default=True; only applies to model option)
    node_set : list[int] / (n, ) ndarray; default=None
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

    .. todo:: node_set still does work on the all the nodes in the big
               kdtree loop, which is very inefficient
    .. todo:: remove_collapsed_elements is not supported
    .. todo:: avoid_collapsed_elements is not supported

    """
    if not isinstance(tol, float):
        tol = float(tol)
    #neq_max = 15
    #method = 'old'

    node_set = _simplify_node_set(node_set)
    model, nid_pairs = _bdf_equivalence_nodes(
        bdf_filename, tol,
        renumber_nodes=renumber_nodes, neq_max=neq_max,
        xref=xref, node_set=node_set, log=log, debug=debug,
        method=method,
        idtype='int32', fdtype='float64')
    model.log.debug(f'equivalence {len(nid_pairs):d} nodes')

    if bdf_filename_out is not None:
        model.write_bdf(bdf_filename_out, size=size, is_double=is_double)
    else:
        model.log.debug('skipping equivalence write')

    if crash_on_collapse:
        # lazy way to make sure there aren't any collapsed nodes
        model2 = BDF(log=log, debug=debug)
        model2.read_bdf(bdf_filename_out)
    return model

#def _simplify_node_set_old(node_set: Optional[Union[list[int], list[NDArrayNint]]],
                           #idtype: str='int32') -> Optional[list[NDArrayNint]]:
    #if node_set is None:
        #return
    #if isinstance(node_set, np.ndarray):
        #return node_set
    ## list
    #assert isinstance(node_set, list), type(node_set)
    #node = node_set[0]
    #if isinstance(node, integer_types):
        #return np.array(node_set, dtype=idtype)
    #raise NotImplementedError(type(node))
    # list of ndarrays
    #return node_set

def _simplify_node_set(node_set: Optional[Union[list[int],
                                                set[int],
                                                list[NDArrayNint]]],
                       idtype: str='int32') -> Optional[list[NDArrayNint]]:  # pragma: no cover
    """
    accepts multiple forms of the node_set parameter
     - list[int]
     - set[int]
     - int ndarray
     - list[int ndarray]

     """
    if node_set is None:
        return
    if isinstance(node_set, np.ndarray):
        return [node_set]
    elif isinstance(node_set, set):
        node_set_array = np.asarray(list(node_set), dtype=idtype)
        node_set_array.sort()
        return [node_set_array]

    assert isinstance(node_set, list), type(node_set)
    node = node_set[0]
    if isinstance(node, integer_types):
        # list
        node_set_array = np.array(node_set, dtype=idtype)
        return [node_set_array]
    else:
        # list of lists
        # list of numpy arrays
        node_set_list = []
        for node_seti in node_set:
            if isinstance(node_seti, list):
                node_set_array = np.array(node_seti, dtype=idtype)
            elif isinstance(node_seti, set):
                node_set_array = np.array(list(node_seti), dtype=idtype)
            else:
                assert isinstance(node_seti, np.ndarray), type(node_seti)
                node_set_array = node_seti
            node_set_array.sort()
            node_set_list.append(node_set_array)
    # list of ndarrays
    return node_set_list

def _bdf_equivalence_nodes(bdf_filename: str, tol: float,
                           renumber_nodes: bool=False, neq_max: int=4, xref: bool=True,
                           node_set: Optional[list[NDArrayNint]]=None,
                           log: Optional[SimpleLogger]=None,
                           debug: bool=True,
                           method: str='new',
                           idtype: str='int32',
                           fdtype: str='float64') -> tuple[BDF,
                                                           list[tuple[int, int]]]:
    """helper for bdf_equivalence_nodes"""
    all_node_set = get_all_node_set(node_set)
    nodes_xyz, model, nids, inew = _eq_nodes_setup(
        bdf_filename, renumber_nodes=renumber_nodes,
        xref=xref, node_set=node_set, log=log, debug=debug,
        idtype=idtype, fdtype=fdtype)

    log = model.log
    log.debug(f'bdf_equivalence_nodes; tol={tol}')

    nid_pairs = _nodes_xyz_nids_to_nid_pairs(
        nodes_xyz, nids, all_node_set,
        tol, log, inew,
        node_set=node_set, neq_max=neq_max, method=method, debug=debug)
    _eq_nodes_final(nid_pairs, model, tol, all_node_set, debug=debug)
    return model, nid_pairs

def _eq_nodes_setup(bdf_filename,
                    renumber_nodes=False, xref=True,
                    node_set: Optional[list[NDArrayNint]]=None,
                    log: Optional[SimpleLogger]=None,
                    debug: bool=True,
                    idtype: str='int32',
                    fdtype: str='float64') -> tuple[NDArrayN3float, BDF,
                                                    NDArrayNint, NDArrayNint]:
    """helper function for ``bdf_equivalence_nodes``"""
    if node_set is not None:
        if renumber_nodes:
            raise NotImplementedError('node_set is not None & renumber_nodes=True')
        #print('*node_set', node_set)
        assert len(node_set) > 0, node_set
        assert isinstance(node_set, list), type(node_set)
    all_node_set = get_all_node_set(node_set)

    model = get_bdf_model(bdf_filename, xref=xref, log=log, debug=debug)

    # map the node ids to the slot in the nids array
    renumber_nodes = False
    if node_set is not None:
        nids, all_nids = _eq_nodes_setup_node_set(
            model, node_set, all_node_set,
            renumber_nodes=renumber_nodes, idtype=idtype)
    else:
        nids, all_nids = _eq_nodes_setup_node(
            model, renumber_nodes=renumber_nodes, idtype=idtype)

    nodes_xyz = _get_xyz_cid0(model, nids, fdtype=fdtype)
    inew = _check_for_referenced_nodes(model, node_set, nids, all_nids, nodes_xyz)

    #assert np.array_equal(nids[inew], nids_new), 'some nodes are not defined'
    return nodes_xyz, model, nids, inew

def _get_xyz_cid0(model: BDF, nids: NDArrayNint,
                  fdtype: str='float32') -> NDArrayN3float:
    """gets xyz_cid0"""
    coord_ids = model.coord_ids
    needs_get_position = (coord_ids == [0])

    if needs_get_position:
        nodes_xyz = np.array([model.nodes[nid].get_position()
                              for nid in nids], dtype=fdtype)
    else:
        nodes_xyz = np.array([model.nodes[nid].xyz
                              for nid in nids], dtype=fdtype)
    return nodes_xyz

def _eq_nodes_setup_node_set(model: BDF,
                             node_set: list[NDArrayNint],
                             all_node_set: NDArrayNint,
                             renumber_nodes: bool=False,
                             idtype:str='int32') -> tuple[NDArrayNint, NDArrayNint]:
    """helper function for ``_eq_nodes_setup`` that handles node_sets"""
    if len(node_set) > 1:
        model.log.warning(f'multi node_sets; n={len(node_set)}')

    node_list = list(model.nodes.keys())
    all_nids = np.array(node_list, dtype=idtype)
    #all_nids.sort()

    # B - A
    # these are all the nodes that are requested from all_node_set that are missing
    #   thus len(diff_nodes) == 0
    diff_nodes = np.setdiff1d(all_node_set, all_nids)
    if len(diff_nodes) != 0:
        msg = ('The following nodes cannot be found, but are included'
               ' in the reduced set; nids=%s' % diff_nodes)
        raise RuntimeError(msg)

    # A & B
    # the nodes to analyze are the union of all the nodes and the desired set
    # which is basically the same as:
    #   nids = unique(all_node_set)
    nids = np.intersect1d(all_nids, all_node_set, assume_unique=True)  # the new values

    if renumber_nodes:
        raise NotImplementedError('node_set is not None & renumber_nodes=True')

    #nids = array([node.nid for nid, node in sorted(model.nodes.items())
                    #if nid in node_set], dtype='int32')
    return nids, all_nids

def _eq_nodes_setup_node(model: BDF, renumber_nodes: bool=False,
                         idtype: str='int32') -> tuple[NDArrayNint, NDArrayNint]:
    """helper function for ``_eq_nodes_setup`` that doesn't handle node sets

    Returns
    -------
    nids : (nnode,) int array
        ???
    all_nids : (nnode,) int array
        all the GRID ids

    """
    inode = 0
    if renumber_nodes:
        model.log.info('renumbering nodes')
        for nid, node in sorted(model.nodes.items()):
            node.nid = inode + 1
            inode += 1
        nnodes = len(model.nodes)
        nids = np.arange(1, inode + 1, dtype=idtype)
        assert nids[-1] == nnodes
    else:
        nodes_list = list(model.nodes.keys())
        nids = np.array(nodes_list, dtype=idtype)
        nids.sort()
    all_nids = nids
    return nids, all_nids

def _check_for_referenced_nodes(model: BDF,
                                node_set: Optional[NDArrayNint],
                                nids: NDArrayNint,
                                all_nids: NDArrayNint,
                                nodes_xyz: NDArrayN3float) -> Optional[NDArrayNint]:
    """helper function for ``_eq_nodes_setup``

    Parameters
    ----------
    nodes_xyz : (nnode, 3) float array
        the xyz for nids

    """
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

        nids_new: set[Optional[int]] = (
            spoint_epoint_nid_set - set(model.spoints) - set(model.epoints)
        )
        if None in nids_new:
            nids_new.remove(None)

        # autosorts the data
        nids_new = np.unique(list(nids_new))
        assert isinstance(nids_new[0], integer_types), type(nids_new[0])

        missing_nids = list(set(nids_new) - set(all_nids))
        if missing_nids:
            missing_nids.sort()
            msg = 'There are missing nodes...\n'  # TODO: in what???
            msg = 'missing nids=%s' % str(missing_nids)
            raise RuntimeError(msg)

        # get the node_id mapping for the kdtree
        inew = np.searchsorted(nids, nids_new, side='left')
        # print('nids_new =', nids_new)
    else:
        inew = None
    #assert np.array_equal(nids[inew], nids_new), 'some nodes are not defined'
    return inew

def get_all_node_set(node_set: Optional[list[NDArrayNint]]) -> NDArrayNint:
    if node_set is None:
        all_node_set = np.array([])
    else:
        all_node_set = np.hstack(node_set)
    return all_node_set

def _eq_nodes_find_pairs(nids: NDArrayNint,
                         slots, ieq,
                         log: SimpleLogger,
                         all_node_set: NDArrayNint,
                         node_set: Optional[list[NDArrayNint]]=None) -> list[tuple[int, int]]:
    """helper function for `bdf_equivalence_nodes`"""
    irows, icols = slots
    if len(irows) == 0:
        return []
        #return np.array([], dtype='int32')

    all_node_set = get_all_node_set(node_set)
    if node_set is not None and len(node_set) > 1:
        log.warning(f'multi node_sets; n={len(node_set)}')

    #replacer = unique(ieq[slots])  ## TODO: turn this back on?

    #skip_nodes = []
    nid_pairs = []
    if node_set is None:
        inids2 = ieq[irows, icols]
        nids1 = nids[irows]
        nids2 = nids[inids2]
        is_unique = (nids1 != nids2)
        nid1_nid2 = np.column_stack([
            nids1[is_unique],
            nids2[is_unique],
        ])
        nid_pairs.append(nid1_nid2)
    else:
        #node_set is not None
        inids2 = ieq[irows, icols]
        nids1 = nids[irows]
        nids2 = nids[inids2]
        is_unique = (nids1 != nids2)
        nid1_nid2 = np.column_stack([
            nids1[is_unique],
            nids2[is_unique],
        ])
        for nid1, nid2 in nid1_nid2:
            if nid1 not in all_node_set and nid2 not in all_node_set:
                continue
            for seti in node_set:
                if nid1 in seti and nid2 in seti:
                    nid_pairs.append((nid1, nid2))
                    #print(f'({nid1}, {nid2})')
                    break

    if len(nid_pairs) == 0:
        return []
        #return np.array([], dtype='int32')

    nid_pairs_array = np.vstack(nid_pairs)
    nid_pairs_array.sort(axis=1) #  inplace sort [min, max]
    #assert nid_pairs_array.shape == nid_pairs_array2.shape

    nid_pairs_array2 = unique2d(nid_pairs_array).tolist()
    nid_pairs_3 = [tuple(nid_pair) for nid_pair in nid_pairs_array2]
    #if len(nid_pairs) > 2_000:
        #asdf
    #return nid_pairs_array3
    return nid_pairs_3


def _eq_nodes_final(nid_pairs: list[tuple[int, int]],
                    model: BDF,
                    tol: float,
                    all_node_set: NDArrayNint,
                    debug: bool=False) -> None:
    """apply nodal equivalencing to model"""
    #log = model.log
    #log.info('_eq_nodes_final')
    if len(all_node_set) == 0:
        # node_sets is None
        for (nid1, nid2) in nid_pairs:
            node1 = model.nodes[nid1]
            node2 = model.nodes[nid2]

            xyz1 = node1.get_position()
            xyz2 = node2.get_position()
            distance = norm(xyz1 - xyz2)
            if distance > tol:
                continue
            _update_grid(node1, node2)
        return

    for (nid1, nid2) in nid_pairs:
        node1 = model.nodes[nid1]
        node2 = model.nodes[nid2]

        xyz1 = node1.get_position()
        xyz2 = node2.get_position()
        distance = norm(xyz1 - xyz2)

        #print('  irow=%s->n1=%s icol=%s->n2=%s' % (irow, nid1, icol, nid2))
        if distance > tol:
            #print('  *n1=%-4s xyz=%s\n  *n2=%-4s xyz=%s\n  *distance=%s\n' % (
               #nid1, xyz1.tolist(),
               #nid2, xyz2.tolist(),
               #distance))
            continue

        assert nid1 in all_node_set, 'nid1=%s all_node_set=%s' % (nid1, all_node_set)
        assert nid2 in all_node_set, 'nid2=%s all_node_set=%s' % (nid2, all_node_set)
        #print('  n1=%-4s xyz=%s\n  n2=%-4s xyz=%s\n  distance=%s\n' % (
            #nid1, str(node1.xyz),
            #nid2, str(node2.xyz),
            #distance))

        # if hasattr(node2, 'new_node_id'):
        # else:
        #if xyz1[1] < 1.0:
            #print('  *n1=%-4s xyz=%s\n  *n2=%-4s xyz=%s\n  *distance=%s\n' % (
               #nid1, xyz1.tolist(),
               #nid2, xyz2.tolist(),
               #distance))
        _update_grid(node1, node2)

        # node2.new_nid = node1.nid
        #skip_nodes.append(nid2)
    return

def _update_grid(node1: GRID, node2: GRID) -> None:
    """helper method for _eq_nodes_final"""
    node2.nid = node1.nid
    node2.xyz = node1.xyz
    node2.cp = node1.cp
    assert node2.cd == node1.cd
    assert node2.ps == node1.ps
    assert node2.seid == node1.seid
    #node1.cp_ref = None
    #node2.cp_ref = None

def _nodes_xyz_nids_to_nid_pairs(nodes_xyz: NDArrayN3float,
                                 nids: NDArrayNint,
                                 all_node_set: NDArrayNint,
                                 tol: float,
                                 log: SimpleLogger,
                                 inew: NDArrayNint,
                                 node_set: Optional[NDArrayNint]=None,
                                 neq_max: int=4,
                                 method: str='new',
                                 debug: bool=False) -> list[tuple[int, int]]:
    """
    Helper for equivalencing

    Returns
    -------
    nid_pairs : list[tuple[int, int]]
        a series of (nid1, nid2) pairs

    """
    if tol < 0.0:
        return []
    unused_kdt, nid_pairs = _eq_nodes_build_tree(
        nodes_xyz, nids, all_node_set,
        tol, log,
        inew=inew, node_set=node_set,
        neq_max=neq_max, method=method, debug=debug)
    return nid_pairs

def _nodes_xyz_nids_to_nid_pairs_new(kdt: KDTree,
                                     nids: NDArrayNint,
                                     all_node_set: NDArrayNint,
                                     node_set: Optional[NDArrayNint],
                                     tol: float,
                                     log: SimpleLogger) -> list[tuple[int, int]]:
    """
    helper function for `bdf_equivalence_nodes`
    """
    ieq3 = kdt.query_ball_tree(kdt, tol)
    #log.info('made ieq3')
    nid_pairs = set()

    #print(ieq3)

    if node_set is None:
        for pair in ieq3:
            if len(pair) == 1:
                continue
            # the combinations should be paired with 2 nodes in each group
            for inid1, inid2 in combinations(pair, 2):
                nid1 = nids[inid1]
                nid2 = nids[inid2]
                pair = (nid1, nid2)
                nid_pairs.add(pair)
        return list(nid_pairs)

    nsets = len(node_set)
    #if nsets > 1:
        #warnings.warn('multiple node_sets not handled in _nodes_xyz_nids_to_nid_pairs_new')

    for pair in ieq3:
        if len(pair) == 1:
            continue
        # the combinations should be paired with 2 nodes in each group
        for inid1, inid2 in combinations(pair, 2):
            nid1 = nids[inid1]
            nid2 = nids[inid2]
            #if nid1 == nid2:
                #continue
            if nid1 not in all_node_set and nid2 not in all_node_set:
                continue
            pair = (nid1, nid2)
            nid_pairs.add(pair)

    if nsets > 1:
        # nid_pairs was simply the set of all potential connections
        # now we filter connections that aren't part of an explicit set
        nid_pairs2 = set([])
        for pair in nid_pairs:
            for seti in node_set:
                nid1, nid2 = pair
                if nid1 in seti and nid2 in seti:
                    nid_pairs2.add(pair)
                    break
        #return nid_pairs2
        nid_pairs = nid_pairs2
    nid_pairs3 = np.array(list(nid_pairs))
    nid_pairs4 = unique2d(nid_pairs3).tolist()
    nid_pairs5 = [tuple(nid_pair) for nid_pair in nid_pairs4]
    return nid_pairs5

def _eq_nodes_build_tree(nodes_xyz: NDArrayN3float,
                         nids: NDArrayNint,
                         all_node_set: NDArrayNint,
                         tol: float,
                         log: SimpleLogger,
                         inew=None,
                         node_set: Optional[NDArrayNint]=None,
                         neq_max: int=4,
                         method: str='new',
                         msg: str='',
                         debug: bool=False) -> tuple[KDTree,
                                                     list[tuple[int, int]]]:
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
    node_set : list[int] / (n, ) ndarray; default=None
        the list/array of nodes to consider
    neq_max : int; default=4
        the number of nodes to consider for equivalencing
    msg : str; default=''
        custom message used for errors

    Returns
    -------
    kdt : KDTree()
        the kdtree object
    nid_pairs : list[tuple[int, int]]
        a series of (nid1, nid2) pairs

    """
    #log.info('_eq_nodes_build_tree')
    nnodes = len(nids)
    if inew is None:
        inew = slice(None)

    assert isinstance(tol, float), 'tol=%r' % tol
    kdt = _get_tree(nodes_xyz, msg=msg)
    #log.info('made kdt')

    is_not_node_set = inew is None or inew == slice(None)

    # check the closest 10 nodes for equality
    if method == 'new' and is_not_node_set:
        #log.info('method=new is_not_node_set')
        nid_pairs = _eq_nodes_build_tree_new(
            kdt, nodes_xyz,
            nids, all_node_set,
            nnodes, is_not_node_set,
            tol, log,
            inew=inew, node_set=node_set, neq_max=neq_max, msg=msg,
            debug=debug)
    else:
        #log.info('method=old')
        if method == 'new':
            log.warning(f'setting method to "old" because node_set is specified')

        #ieq : (nnodes, neq) int ndarray
        #    the node indices that are close
        #slots : (nnodes, neq) int ndarray
        #    the location of where
        deq, ieq = kdt.query(nodes_xyz[inew, :], k=neq_max, distance_upper_bound=tol)
        if node_set is not None:
            assert len(deq) == len(nids)

        # get the ids of the duplicate nodes
        slots = np.where(ieq[:, :] < nnodes)

        if ieq[:, -1].max() == nnodes:
            log.warning(f'neq_max={neq_max} and should be increased')
        nid_pairs = _eq_nodes_find_pairs(nids, slots, ieq, log,
                                         all_node_set, node_set=node_set)
    assert isinstance(nid_pairs, list), nid_pairs
    return kdt, nid_pairs

def _eq_nodes_build_tree_new(kdt: KDTree,
                             nodes_xyz: NDArrayN3float,
                             nids: NDArrayNint,
                             all_node_set: NDArrayNint,
                             nnodes: int,
                             is_not_node_set: bool,
                             tol: float,
                             log: SimpleLogger,
                             inew=None, node_set=None, neq_max: int=4, msg: str='',
                             debug: float=False) -> list[tuple[int, int]]:
    assert isinstance(nnodes, int), nnodes
    #log.info('calling query')
    deq, ieq = kdt.query(nodes_xyz[inew, :], k=neq_max, distance_upper_bound=tol)
    #log.info('building slots')
    slots = np.where(ieq[:, :] < nnodes)
    #log.info('building pairs')
    nid_pairs_expected = _eq_nodes_find_pairs(
        nids, slots, ieq, log, all_node_set, node_set=node_set)
    if is_not_node_set:
        #log.info('creating new pairs')
        nid_pairs = _nodes_xyz_nids_to_nid_pairs_new(
            kdt, nids, all_node_set,
            node_set, tol, log)
    else:  # pragma: no cover
        raise NotImplementedError(f'node_set = {node_set}')

    #log.info('final checks...')
    snid_pairs = set(nid_pairs)
    snid_pairs_expected = set(nid_pairs_expected)
    diff_bad = snid_pairs - snid_pairs_expected
    diff_missed = snid_pairs - snid_pairs_expected
    if debug and len(diff_bad) or len(diff_missed):  # pragma: no cover
        #log.warning(f'nid_pairs          = {nid_pairs}')
        #log.warning(f'nid_pairs_expected = {nid_pairs_expected}')
        if len(diff_bad) + len(diff_missed) < 20:
            log.warning(f'diff_bad = {diff_bad}')
            log.warning(f'diff_missed = {diff_missed}')
        else:
            log.warning(f'diff_bad diff_missed are large...cause?')

    return nid_pairs

def _get_tree(nodes_xyz: NDArrayN3float, msg: str='') -> KDTree:
    """gets the kdtree"""
    assert isinstance(nodes_xyz, np.ndarray), type(nodes_xyz)
    assert nodes_xyz.shape[0] > 0, 'nnodes=0%s' % msg

    # build the kdtree
    try:
        kdt = KDTree(nodes_xyz)
    except RuntimeError:
        print(nodes_xyz)
        raise RuntimeError(nodes_xyz)
    return kdt
