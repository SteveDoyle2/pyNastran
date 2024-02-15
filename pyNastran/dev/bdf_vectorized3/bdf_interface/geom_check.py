from __future__ import annotations
import numpy as np
from typing import TYPE_CHECKING
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.cards.base_card import VectorizedBaseCard


def geom_check(self: VectorizedBaseCard,
               missing: dict[str, np.ndarray],
               node=None, filter_node0: bool=False,
               coord=None,
               property_id=None, material_id=None,
               element_id=None,
               mpc=None, spc=None,
               spoint=None,
               aecomp=None,
               set1=None,
               caero=None,
               delay=None,
               dphase=None,
               dconstr_id=None) -> None:

    _geom_check_node(missing, node, filter_node0=filter_node0)
    if node is not None:
        nid, node_id = node
        if len(nid) == 0:
            unid = np.unique(node_id)
            #self.log.error(f'node_id={nid} missing nodes={unid}')
            missing['node_id'] = unid
            return
        assert len(nid) > 0, nid
        inode, umissing = find_missing(nid, node_id, 'nodes')
        if filter_node0:
            umissing = np.setdiff1d(umissing, [0])
        if len(umissing):
            missing['node_id'] = umissing
        del inode, umissing, nid, node_id

    _geom_check(missing, coord, 'coord_id', 'coords')
    _geom_check(missing, property_id, 'property_id', 'properties')
    _geom_check(missing, element_id, 'element_id', 'elements')
    _geom_check(missing, material_id, 'material_id', 'materials')
    _geom_check(missing, aecomp, 'aecomp_id', 'aecomps')
    _geom_check(missing, set1, 'set1', 'set1s')
    _geom_check(missing, delay, 'delay', 'delays', filter0=True)
    _geom_check(missing, dphase, 'dphase', 'dphases', filter0=True)
    _geom_check(missing, dconstr_id, 'dconstr_id', 'dconstr')

    if caero is not None:
        all_caero_ids, used_caero_ids = caero
        assert len(all_caero_ids) > 0, all_caero_ids
        if len(used_caero_ids) > 0:
            #assert len(used_caero_ids) > 0, used_caero_ids
            icaero, umissing = find_missing(all_caero_ids, used_caero_ids, 'caeros')
            if len(umissing):
                missing['caero_id'] = umissing
            del icaero, umissing
        del all_caero_ids, used_caero_ids

def _geom_check_node(missing: dict[str, np.ndarray],
                     node, filter_node0: bool):
    if node is None:
        return
    nid, node_id = node
    if len(nid) == 0:
        unid = np.unique(node_id)
        #self.log.error(f'node_id={nid} missing nodes={unid}')
        missing['node_id'] = unid
        return
    assert len(nid) > 0, nid
    inode, umissing = find_missing(nid, node_id, 'nodes')
    if filter_node0:
        umissing = np.setdiff1d(umissing, [0])
    if len(umissing):
        missing['node_id'] = umissing
    #del inode, umissing, nid, node_id

def _geom_check(missing: dict[str, np.ndarray],
                group: list[np.ndarray, np.ndarray],
                name: str,
                names: str, filter0: bool=True):
    """
    _geom_check(missing, coord, 'coord_id', 'coords')
    """
    if group is not None:
        all_ids, used_ids = group
        if filter0:
            # you don't need a DELAY if the only used ID=0
            used_ids = np.setdiff1d(used_ids, [0])
            if len(used_ids) == 0:
                return

        if len(all_ids) == 0:
            #log.warning(f'{name}: no ids: {all_ids}')
            umissing = np.unique(used_ids)
        else:
            iset, umissing = find_missing(all_ids, used_ids, names)
        if len(umissing):
            missing[name] = umissing

def find_missing(all_nodes_sorted: np.ndarray, my_nodes: np.ndarray,
                 name: str) -> tuple[np.ndarray, np.ndarray]:
    assert isinstance(all_nodes_sorted, np.ndarray), all_nodes_sorted
    assert isinstance(my_nodes, (np.ndarray, list)), my_nodes

    unids = np.unique(all_nodes_sorted)
    if not np.array_equal(unids, all_nodes_sorted):
        raise RuntimeError(f'all_{name}_sorted={all_nodes_sorted} is not unique')
    inode = np.searchsorted(all_nodes_sorted, my_nodes)
    if my_nodes.ndim == 1:
        nnodes = len(all_nodes_sorted)
        ibad = np.where(inode == nnodes)
        #if len(ibad[0]):
            #raise NotImplementedError('check this case...probably fine')
        inode[ibad] = -1
    elif my_nodes.ndim == 2:
        nnodes = len(all_nodes_sorted)
        ibad = np.where(inode == nnodes)
        inode[ibad] = -1
        #if len(ibad[0]):
            #print('all_nodes_sorted', all_nodes_sorted)
            #print('my_nodes\n', my_nodes)
            #print(inode)
            #raise NotImplementedError('check this case...probably fine')
    else:
        raise RuntimeError(my_nodes.shape)

    try:
        actual = all_nodes_sorted[inode]
    except IndexError:
        unique_my_nodes_sorted = np.unique(my_nodes.ravel())
        print(f'all_{name}_sorted={all_nodes_sorted}')
        print(f'unique_my_{name}_sorted={unique_my_nodes_sorted}')
        raise
    expected = my_nodes
    is_missing = (actual != expected)
    imissing = expected[is_missing]
    umissing = np.unique(imissing)
    return inode, umissing
