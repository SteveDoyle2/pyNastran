"""
defines:
    free_edges(model)
"""
from __future__ import print_function
from collections import defaultdict
from six import iteritems


def free_edges(model):
    """gets the free edges for shell elements"""
    edge_to_eids = defaultdict(set)
    shell_elements = ['CTRIA3', 'CTRIAX', 'CTRIA6', 'CTRIAX6',
                      'CQUAD4', 'CQUAD', 'CQUAD8', 'CQUADR', 'CQUADX', 'CQUADX8',
                      'CSHEAR']
    for eid, elem in iteritems(model.elements):
        if elem.type not in shell_elements:
            continue
        edges = elem.get_edge_ids()
        for edge in edges:
            edge_to_eids[edge].add(eid)

    free_edges = []
    for edge, eids in iteritems(edge_to_eids):
        if len(eids) == 1:
            free_edges.append(edge)
    return free_edges
