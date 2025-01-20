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
from pathlib import PurePath
from io import StringIO
from typing import Optional, TYPE_CHECKING
import numpy as np

from pyNastran.nptyping_interface import NDArrayNint, NDArrayN3float
#from pyNastran.bdf.mesh_utils.internal_utils import get_bdf_model
from pyNastran.bdf.mesh_utils.bdf_equivalence import (
    get_all_node_set,
    _simplify_node_set, _check_for_referenced_nodes,
    _nodes_xyz_nids_to_nid_pairs,
)

from pyNastran.dev.bdf_vectorized3.bdf import BDF
if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger


BDF_FILETYPE = BDF | str | StringIO | PurePath
def get_bdf_model(bdf_filename: BDF_FILETYPE,
                  xref: bool=True,
                  cards_to_skip=None,
                  validate: bool=True,
                  log=None, debug: bool=False) -> BDF:
    if isinstance(bdf_filename, (str, StringIO, PurePath)):
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
    else:  # pragma: no cover
        raise NotImplementedError(bdf_filename)
    return model

def bdf_equivalence_nodes(bdf_filename: str,
                          bdf_filename_out: Optional[str],
                          tol: float,
                          renumber_nodes: bool=False, neq_max: int=4, xref: bool=True,
                          node_set: Optional[list[int] | NDArrayNint]=None,
                          size: int=8, is_double: bool=False,
                          remove_collapsed_elements: bool=False,
                          avoid_collapsed_elements: bool=False,
                          crash_on_collapse: bool=False,
                          log: Optional[SimpleLogger]=None,
                          debug: bool=True, method: str='old') -> BDF:
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


def _eq_nodes_final(nid_pairs: list[tuple[int, int]],
                    model: BDF,
                    tol: float,
                    all_node_set: NDArrayNint,
                    debug: bool=False) -> None:
    """apply nodal equivalencing to model"""
    if len(nid_pairs) == 0:
        return
    #log = model.log
    #log.warning('_eq_nodes_final')
    #log.info('_eq_nodes_final')
    grid = model.grid

    node_id = grid.node_id
    all_node_ids = node_id.copy()
    new_node_ids = node_id.copy()

    xyz_cid0 = grid.xyz_cid0()
    idtype = model.idtype
    nid_pairs2i = [(nid1, nid2) if nid1 < nid2 else (nid2, nid1)
                   for (nid1, nid2) in nid_pairs]
    nid_pairs_array = np.array(nid_pairs2i, dtype=idtype)
    inode = np.searchsorted(node_id, nid_pairs_array)

    use_old_method = False
    # update the nodes
    nid_old_to_new = {}
    for i1, i2 in inode:
        xyz1 = xyz_cid0[i1, :]
        xyz2 = xyz_cid0[i2, :]
        dxyz = xyz2 - xyz1
        normi = np.linalg.norm(dxyz)
        if normi <= tol:
            #iupdate = (normi <= tol)
            for nid1, nid2 in nid_pairs_array:
                nid_new0 = new_node_ids[i1]
                new_node_ids[i2] = nid_new0
                if use_old_method:
                    nid_new = grid.node_id[i1]
                    #nid_old = grid.node_id[i2]
                    grid.node_id[i2] = nid_new
                    grid.xyz[i2, :] = grid.xyz[i1, :]
                    grid.cp[i2] = grid.cp[i1]
                    grid.cd[i2] = grid.cd[i1]
                    grid.seid[i2] = grid.seid[i1]
                    grid.ps[i2] = grid.ps[i1]
                    nid_old_to_new[nid2] = nid1

    #print(f'reduced-set:\n{model.grid.write()}')

    if use_old_method:
        # =0: node_sets is None -> no sets are used
        if len(all_node_set) != 0:
            for nid1, nid2 in nid_old_to_new.items():
                assert nid1 in all_node_set, 'nid1=%s all_node_set=%s' % (nid1, all_node_set)
                assert nid2 in all_node_set, 'nid2=%s all_node_set=%s' % (nid2, all_node_set)

    update_cards(model, all_node_ids, new_node_ids, nid_old_to_new)
    return


def update_cards(model: BDF,
                 all_node_ids: np.ndarray,
                 new_node_ids: np.ndarray,
                 nid_old_to_new: dict[int, int]) -> None:
    """
    The ids have been updated, but now we need to slice out
    invalid ids
    """
    log = model.log
    log.warning('update_cards')
    #print(f'new_node_ids = {new_node_ids}')
    unique_new_node_ids = np.unique(new_node_ids)
    #print(f'unique_new_node_ids = {unique_new_node_ids}')

    no_equiv_cards = {
        'MAT1', 'MAT2', 'MAT3', 'MAT4', 'MAT5', 'MAT8', 'MAT9', 'MAT10',
        'MATT1', 'MATT2', 'MATT8',
        'MATS1',
        'PELAS', 'PELAST',
        'PDAMP', 'PDAMPT',
        'PBUSH', 'PBUSHT', 'PBUSH1D',
        'PMASS', 'PGAP', 'PVISC', 'PFAST',
        'PROD', 'PTUBE',
        'PBAR', 'PBARL',
        'PBEAM', 'PBEAML', 'PBCOMP', 'PBEND', 'PBEAM3',
        'PSHELL', 'PCOMP', 'PCOMPG', 'PSHEAR',
        'PSOLID', 'PLSOLID',
        # aero
        'PAERO1',
        'SPLINE1', 'SPLINE2', 'SPLINE3', 'SPLINE4', 'SPLINE5', 'SPLINE6',
        'AESURF', 'AESURFS', 'AESTAT',
        'AEFACT',
        # loads
        'GRAV',
        # acoustic
        'PAABSF',
        # contact
        'BSURF', 'BSURFS', 'BCPROP', 'BCPROPS', 'BCTSET',
        'BGADD', 'BCTADD', 'BFRIC', 'BCONP',
        # bolt
        'BOLTLD',
    }
    skip_cards = {
        'GRID', 'SPOINT', 'EPOINT',
        'CBARAO',
        'CELAS3', 'CELAS4',
        'CDAMP3', 'CDAMP4',
        'CMASS3', 'CMASS4',
        'SPCADD', 'MPCADD', 'SET1',
        # aero
        'AELIST',
        'CAERO1', 'CAERO2', 'CAERO3', 'CAERO4', 'CAERO5', 'CAERO7',
        'PAERO1', 'PAERO2', 'PAERO3', 'PAERO4', 'PAERO5',
        #  optimization
        'DESVAR', 'DVPREL1', 'DVMREL1', 'DVPREL2', 'DVMREL2',
        'DCONSTR', 'DCONADD', 'DSCREEN',
        'PLOAD1', 'PLOAD2', 'SLOAD', 'LOAD', 'DEFORM',
        # thermal
        'RADM', 'TEMPD', 'RADSET',
        # time/freq/random loads
        'TLOAD1', 'TLOAD2', 'RLOAD1', 'RLOAD2', 'RANDPS', 'ACSRCE',
        'DLOAD', 'LSEQ',
        # not supported
        'SEBSET', 'SECSET', 'SEQSET',
        'BOLT', 'BOLTFOR', 'BOLTSEQ',
    } | no_equiv_cards
    grid = model.grid
    ids = np.unique(grid.node_id)

    if 0:  # pragma: no cover
        nids_to_remove_list = list(nid_old_to_new.keys())
        nids_to_remove = np.unique(nids_to_remove_list)
        nids_to_keep = np.setdiff1d(ids, nids_to_remove)
        reverse_index = grid.index(nids_to_keep, assume_sorted=False, check_index=True)
    else:
        #reverse_index = np.searchsorted(all_node_ids, ids)
        reverse_index = np.searchsorted(all_node_ids, unique_new_node_ids)
    #print(f'reverse_index = {reverse_index}')
    # inplace operation to remove duplicate ids
    grid.__apply_slice__(grid, reverse_index)
    grid.sort()
    #grid.slice_card_by_id(ids, assume_sorted=True, sort_ids=False)

    #print(f'final-set:\n{model.grid.write()}')

    unid = np.unique(model.grid.node_id)
    assert len(unid) == len(model.grid)

    if len(model.set1):
        log.error('SET1 requires a node/element flag from another card')

    #supported_cards = {}
    cards = [card for card in model._cards_to_setup if card.n > 0]
    for card in cards:
        if card.type in skip_cards:
            continue
        try:
            card.equivalence_nodes(nid_old_to_new)
        except AttributeError:
            print(card.get_stats())
            raise
