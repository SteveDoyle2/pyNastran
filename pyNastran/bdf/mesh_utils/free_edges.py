"""
defines:
    edges = free_edges(model)
    edges = non_paired_edges(model)
"""
from __future__ import print_function
from collections import defaultdict


def free_edges(model):
    """
    Gets the free edges for shell elements.
    A free edge is an edge that is only connected to 1 shell element.

    Parameters
    ----------
    model : BDF()
        the BDF model

    """
    edge_to_eids = _get_edge_to_eids_map(model)

    free_edges = []
    for edge, eids in edge_to_eids.items():
        if len(eids) == 1:
            free_edges.append(edge)
    return free_edges

def non_paired_edges(model):
    """
    Gets the edges not shared by exactly 2 elements.
    This is useful for identifying rib/spar intersections.

    Parameters
    ----------
    model : BDF()
        the BDF model

    """
    edge_to_eids = _get_edge_to_eids_map(model)

    non_paired_edges = []
    for edge, eids in edge_to_eids.items():
        if len(eids) != 2:
            non_paired_edges.append(edge)
    return non_paired_edges

def _get_edge_to_eids_map(model):
    """helper method"""
    edge_to_eids = defaultdict(set)
    shell_elements = ['CTRIA3', 'CTRIAX', 'CTRIA6', 'CTRIAX6',
                      'CQUAD4', 'CQUAD', 'CQUAD8', 'CQUADR', 'CQUADX', 'CQUADX8',
                      'CSHEAR']
    for eid, elem in model.elements.items():
        if elem.type not in shell_elements:
            continue
        edges = elem.get_edge_ids()
        for edge in edges:
            edge_to_eids[edge].add(eid)
    return edge_to_eids
