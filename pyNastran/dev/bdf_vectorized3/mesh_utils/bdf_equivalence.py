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
from pathlib import PurePath
from io import StringIO
from typing import Union, Optional, Any, TYPE_CHECKING
import numpy as np
#from numpy import array, arange, setdiff1d, intersect1d
from numpy.linalg import norm  # type: ignore
#import scipy
#import scipy.spatial

from pyNastran.nptyping_interface import NDArrayNint, NDArrayN3float
#from pyNastran.utils import int_version
#from pyNastran.utils.numpy_utils import integer_types
#from pyNastran.bdf.mesh_utils.internal_utils import get_bdf_model
from pyNastran.bdf.mesh_utils.bdf_equivalence import (
    _simplify_node_set, _check_for_referenced_nodes,
    get_all_node_set, KDTree)

from pyNastran.dev.bdf_vectorized3.bdf import BDF
if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger
    from pyNastran.bdf.bdf import GRID

#SCIPY_VERSION = int_version('scipy', scipy.__version__)
#import scipy.spatial
#if SCIPY_VERSION > [1, 6, 0]:
    #KDTree = scipy.spatial.KDTree
#else:
    #KDTree = scipy.spatial.cKDTree


BDF_FILETYPE = Union[BDF, str, StringIO, PurePath]
def get_bdf_model(bdf_filename: BDF_FILETYPE,
                  xref: bool=True,
                  cards_to_skip=None,
                  validate: bool=True,
                  log=None, debug: bool=False) -> BDF:
    if isinstance(bdf_filename, (str, StringIO, PurePath)):
        #model = read_bdf(bdf_filename=bdf_filename, validate=True, xref=True,
                        #punch=False, skip_cards=None,
                        #read_cards=None,
                        #encoding=None, log=None,
                        #debug=True, mode='msc')

        xref = True
        model = BDF(log=log, debug=debug)
        model.disable_cards(cards_to_skip)
        model.read_bdf(bdf_filename, validate=validate, xref=xref,
                       punch=False, read_includes=True,
                       save_file_structure=False, encoding=None)
    elif isinstance(bdf_filename, BDF):
        model = bdf_filename
        if xref:
            model.cross_reference(run_geom_check=xref)
    #elif isinstance(bdf_filename, StringIO):
        #model = BDF(log=log, debug=debug)
        #model.read_bdf(bdf_filename, xref=xref)
    else:  # pragma: no cover
        raise NotImplementedError(bdf_filename)
    return model

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

def _bdf_equivalence_nodes(bdf_filename: str, tol: float,
                           renumber_nodes: bool=False, neq_max: int=4, xref: bool=True,
                           node_set: Optional[list[NDArrayNint]]=None,
                           log: Optional[SimpleLogger]=None,
                           debug: bool=True, method: str='new',
                           idtype: str='int32', fdtype: str='float64') -> tuple[BDF,
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

        #print(type(node_set))
        #print('*node_set', node_set)
        assert len(node_set) > 0, node_set
        assert isinstance(node_set, list), type(node_set)
    all_node_set = get_all_node_set(node_set)

    model = get_bdf_model(bdf_filename, xref=xref, log=log, debug=debug)

    # quads / tris
    #nids_quads = []
    #eids_quads = []
    #nids_tris = []
    #eids_tris = []

    # map the node ids to the slot in the nids array
    renumber_nodes = False
    if node_set is not None:
        nids, all_nids = _eq_nodes_setup_node_set(
            model, node_set, all_node_set,
            renumber_nodes=renumber_nodes, idtype=idtype)
        #print(f'sub-set:\n{model.grid.write()}')
    else:
        nids, all_nids = _eq_nodes_setup_node(
            model, renumber_nodes=renumber_nodes, idtype=idtype)
        #print(f'all-set:\n{model.grid.write()}')

    nodes_xyz = model.grid.slice_card_by_id(nids).xyz_cid0()
    inew = _check_for_referenced_nodes(model, node_set, nids, all_nids, nodes_xyz)

    #assert np.array_equal(nids[inew], nids_new), 'some nodes are not defined'
    return nodes_xyz, model, nids, inew

def _eq_nodes_setup_node_set(model: BDF,
                             node_set: list[NDArrayNint],
                             all_node_set: NDArrayNint,
                             renumber_nodes: bool=False,
                             idtype:str='int32') -> tuple[NDArrayNint, NDArrayNint]:
    """helper function for ``_eq_nodes_setup`` that handles node_sets"""
    if len(node_set) > 1:
        model.log.warning(f'multi node_sets; n={len(node_set)}')

    all_nids = model.grid.node_id
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
    assert renumber_nodes is False, renumber_nodes
    if renumber_nodes:
        model.log.info('renumbering nodes')
        for nid, node in sorted(model.nodes.items()):
            node.nid = inode + 1
            inode += 1
        nnodes = len(model.nodes)
        nids = np.arange(1, inode + 1, dtype=idtype)
        assert nids[-1] == nnodes
    else:
        nids = model.grid.node_id.copy()
        unids = np.unique(nids)
        assert len(nids) == len(unids), 'cant equivalence when there are duplicate node ids'
        nids = unids
        del unids
    all_nids = nids
    return nids, all_nids

def _eq_nodes_find_pairs(nids: NDArrayNint,
                         slots, ieq,
                         log: SimpleLogger,
                         all_node_set: NDArrayNint,
                         node_set: Optional[list[NDArrayNint]]=None) -> list[tuple[int, int]]:
    """helper function for `bdf_equivalence_nodes`"""
    irows, icols = slots
    all_node_set = get_all_node_set(node_set)
    if node_set is not None and len(node_set) > 1:
        log.warning(f'multi node_sets; n={len(node_set)}')

    #replacer = unique(ieq[slots])  ## TODO: turn this back on?

    #skip_nodes = []
    nid_pairs = []
    if node_set is None:
        for (irow, icol) in zip(irows, icols):
            inid2 = ieq[irow, icol]
            nid1 = nids[irow]
            nid2 = nids[inid2]
            if nid1 == nid2:
                continue
            nid_pairs.append((nid1, nid2))
        return nid_pairs

    for (irow, icol) in zip(irows, icols):
        inid2 = ieq[irow, icol]
        nid1 = nids[irow]
        nid2 = nids[inid2]
        if nid1 == nid2:
            continue
        if node_set is not None:
            if nid1 not in all_node_set and nid2 not in all_node_set:
                continue
            for seti in node_set:
                if nid1 in seti and nid2 in seti:
                    nid_pairs.append((nid1, nid2))
                    #print(f'({nid1}, {nid2})')
                    break
    return nid_pairs


def _eq_nodes_final(nid_pairs: list[tuple[int, int]],
                    model: BDF,
                    tol: float,
                    all_node_set: NDArrayNint,
                    debug: bool=False) -> None:
    """apply nodal equivalencing to model"""
    if len(nid_pairs) == 0:
        return
    #log = model.log
    #log.info('_eq_nodes_final')
    grid = model.grid

    node_id = grid.node_id
    xyz_cid0 = grid.xyz_cid0()
    idtype = model.idtype
    nid_pairs2i = [(nid1, nid2) if nid1 < nid2 else (nid2, nid1)
                   for (nid1, nid2) in nid_pairs]
    nid_pairs_array = np.array(nid_pairs2i, dtype=idtype)
    inode = np.searchsorted(node_id, nid_pairs_array)

    # update the nodes
    nid_old_to_new = {}
    for i1, i2 in inode:
        xyz1 = xyz_cid0[i1, :]
        xyz2 = xyz_cid0[i2, :]
        dxyz = xyz2 - xyz1
        normi = np.linalg.norm(dxyz)
        #assert len(normi) == len(nid_pairs_array)
        if normi <= tol:
            #iupdate = (normi <= tol)
            for nid1, nid2 in nid_pairs_array:
                nid_new = grid.node_id[i1]
                nid_old = grid.node_id[i2]
                grid.node_id[i2] = nid_new
                grid.xyz[i2, :] = grid.xyz[i1, :]
                grid.cp[i2] = grid.cp[i1]
                grid.cd[i2] = grid.cd[i1]
                grid.seid[i2] = grid.seid[i1]
                grid.ps[i2] = grid.ps[i1]
                nid_old_to_new[nid2] = nid1
    #print(f'reduced-set:\n{model.grid.write()}')

    # =0: node_sets is None -> no sets are used
    if len(all_node_set) != 0:
        for nid1, nid2 in nid_old_to_new.items():
            assert nid1 in all_node_set, 'nid1=%s all_node_set=%s' % (nid1, all_node_set)
            assert nid2 in all_node_set, 'nid2=%s all_node_set=%s' % (nid2, all_node_set)

    #if np.any(iupdate):
        #print(grid.write())
        #inode_update = inode[iupdate, :]
        #for i1, i2 in inode_update:
            #grid.node_id[i2] = grid.node_id[i1]
            #grid.xyz[i2, :] = grid.xyz[i1, :]
            #grid.cp[i2] = grid.cp[i1]
            #grid.cd[i2] = grid.cd[i1]
            #grid.seid[i2] = grid.seid[i1]
            #grid.ps[i2] = grid.ps[i1]
        #print('---------------\nnew:')
        #print(grid.write())

    update_cards(model, nid_old_to_new)
    return


def update_cards(model: BDF,
                 nid_old_to_new: dict[int, int]) -> None:
    skip_cards = {
        'GRID', 'SPOINT', 'EPOINT',
        'MAT1', 'MAT2', 'MAT4', 'MAT5', 'MAT8', 'MAT9', 'MAT10',
        'MATT1', 'MATT2', 'MATT8',
        'MATS1',
        'PELAS', 'PELAST',
        'PDAMP', 'PDAMPT',
        'PBUSH', 'PBUSHT',
        'PMASS', 'PGAP', 'PVISC',
        'PROD', 'PTUBE',
        'PBAR', 'PBARL',
        'PBEAM', 'PBEAML', 'PBCOMP', 'PBEND', 'PBEAM3',
        'PSHELL', 'PCOMP', 'PSHEAR',
        'PSOLID', 'PLSOLID',
        'BCONP', 'BFRIC',
        'SPCADD', # 'MPCADD',
        'DESVAR', 'DVPREL1', 'DVMREL1', 'DVPREL2', 'DVMREL2',
        'DCONSTR', 'DCONADD', 'DSCREEN',
        'GRAV', 'PLOAD1', 'PLOAD2', 'SLOAD',
        # time/freq/random loads
        'TLOAD1', 'TLOAD2', 'RLOAD1', 'RLOAD2', 'RANDPS', 'DLOAD',
    }
    grid = model.grid
    ids = np.unique(grid.node_id)

    nids_to_remove_list = list(nid_old_to_new.keys())
    nids_to_remove = np.unique(nids_to_remove_list)

    nids_to_keep = np.setdiff1d(ids, nids_to_remove)
    reverse_index = grid.index(nids_to_keep, assume_sorted=False, check_index=True)
    grid.__apply_slice__(grid, reverse_index)
    #grid.slice_card_by_id(ids, assume_sorted=True, sort_ids=False)

    #print(f'final-set:\n{model.grid.write()}')

    unid = np.unique(model.grid.node_id)
    assert len(unid) == len(model.grid)

    supported_cards = {}
    cards = [card for card in model._cards_to_setup if card.n > 0]
    for card in cards:
        if card.type in skip_cards:
            continue
        try:
            card.equivalence_nodes(nid_old_to_new)
        except AttributeError:
            print(card.get_stats())
            raise

#def _update_grid(node1: GRID, node2: GRID):
    #"""helper method for _eq_nodes_final"""
    #node2.nid = node1.nid
    #node2.xyz = node1.xyz
    #node2.cp = node1.cp
    #assert node2.cd == node1.cd
    #assert node2.ps == node1.ps
    #assert node2.seid == node1.seid
    #node1.cp_ref = None
    #node2.cp_ref = None

def _nodes_xyz_nids_to_nid_pairs(nodes_xyz: NDArrayN3float,
                                 nids: NDArrayNint,
                                 all_node_set: NDArrayNint,
                                 tol: float,
                                 log: SimpleLogger,
                                 inew: NDArrayNint,
                                 node_set: Optional[NDArrayNint]=None,
                                 neq_max: int=4, method: str='new',
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
                                     node_set: Optional[NDArrayNint], tol: float):
    """
    helper function for `bdf_equivalence_nodes`
    """
    ieq3 = kdt.query_ball_tree(kdt, tol)
    nid_pairs = []

    if node_set is None:
        for pair in ieq3:
            if len(pair) == 1:
                continue
            # the combinations should be paired with 2 nodes in each group
            for inid1, inid2 in combinations(pair, 2):
                nid1 = nids[inid1]
                nid2 = nids[inid2]
                pair = (nid1, nid2)
                if pair in nid_pairs:
                    continue
                nid_pairs.append(pair)
        return nid_pairs

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
            if pair in nid_pairs:
                continue
            nid_pairs.append(pair)

    if nsets > 1:
        # nid_pairs was simply the set of all potential connections
        # now we filter connections that aren't part of an explicit set
        nid_pairs2 = []
        for pair in nid_pairs:
            for seti in node_set:
                nid1, nid2 = pair
                if nid1 in seti and nid2 in seti:
                    nid_pairs2.append(pair)
                    break
        return nid_pairs2
    return nid_pairs

def _eq_nodes_build_tree(nodes_xyz: NDArrayN3float,
                         nids: NDArrayNint,
                         all_node_set: NDArrayNint,
                         tol: float,
                         log: SimpleLogger,
                         inew=None,
                         node_set: Optional[NDArrayNint]=None,
                         neq_max: int=4, method: str='new', msg: str='',
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
            nids, all_node_set,
            nnodes, is_not_node_set,
            tol, log,
            inew=inew, node_set=node_set, neq_max=neq_max, msg=msg,
            debug=debug)
    else:
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
                             debug: float=False) -> tuple[Any, list[tuple[int, int]]]:
    assert isinstance(nnodes, int), nnodes
    deq, ieq = kdt.query(nodes_xyz[inew, :], k=neq_max, distance_upper_bound=tol)
    slots = np.where(ieq[:, :] < nnodes)
    nid_pairs_expected = _eq_nodes_find_pairs(nids, slots, ieq, log, all_node_set, node_set=node_set)
    if is_not_node_set:
        nid_pairs = _nodes_xyz_nids_to_nid_pairs_new(kdt, nids, all_node_set,
                                                     node_set, tol)
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
