"""
defines:
    edges = free_edges(model, eids=None)
    edges = non_paired_edges(model, eids=None)

"""
from __future__ import annotations
from collections import defaultdict
from typing import Tuple, List, Optional, TYPE_CHECKING
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF

def free_edges(model: BDF, eids: Optional[List[int]]=None, maps=None) -> List[Tuple[int, int]]:
    """
    Gets the free edges for shell elements.
    A free edge is an edge that is only connected to 1 shell element.

    Parameters
    ----------
    model : BDF()
        the BDF model
    eids : List[int]; default=None
        a subset of elements to consider
    maps : List[...] (default=None -> calculate)
        the output from _get_maps(eids, map_names=None,
                                  consider_0d=False, consider_0d_rigid=False,
                                  consider_1d=False, consider_2d=True, consider_3d=False)

    Returns
    -------
    edges: List[Tuple[int,int]]
        list of node ids of each edges

    """
    if maps is not None:
        edge_to_eid_map = maps['edge_to_eid_map']
    else:
        edge_to_eid_map = _get_edge_to_eids_map(model, eids=eids)

    edges = []
    for edge, eids in edge_to_eid_map.items():
        if len(eids) == 1:
            edges.append(edge)
    return edges

def non_paired_edges(model: BDF, eids: List[int]=None, maps=None) -> List[Tuple[int, int]]:
    """
    Gets the edges not shared by exactly 2 elements.
    This is useful for identifying rib/spar intersections.

    Parameters
    ----------
    model : BDF()
        the BDF model
    eids : List[int]; default=None
        a subset of elements to consider
    maps : List[...] (default=None -> calculate)
        the output from _get_maps(eids, map_names=None,
                                  consider_0d=False, consider_0d_rigid=False,
                                  consider_1d=False, consider_2d=True, consider_3d=False)

    Returns
    -------
    non_paired_edges : List[(int nid1, int nid2), ...]
        the non-paired edges

    """
    if maps is not None:
        edge_to_eid_map = maps['edge_to_eid_map']
    else:
        edge_to_eid_map = _get_edge_to_eids_map(model, eids=eids)

    edges = []
    for edge, eids in edge_to_eid_map.items():
        if len(eids) != 2:
            edges.append(edge)
    return edges

def _get_edge_to_eids_map(model, eids=None):
    """helper method"""
    edge_to_eids = defaultdict(set)
    shell_elements = ['CTRIA3', 'CTRIAX', 'CTRIA6', 'CTRIAX6',
                      'CQUAD4', 'CQUAD', 'CQUAD8', 'CQUADR', 'CQUADX', 'CQUADX8',
                      'CSHEAR']
    if eids is None:
        for eid, elem in model.elements.items():
            if elem.type not in shell_elements:
                continue
            edges = elem.get_edge_ids()
            for edge in edges:
                edge_to_eids[edge].add(eid)
    else:
        if isinstance(eids, int):
            eids = [eids]
        for eid in eids:
            elem = model.elements[eid]
            if elem.type not in shell_elements:
                continue
            edges = elem.get_edge_ids()
            for edge in edges:
                edge_to_eids[edge].add(eid)

    return edge_to_eids
