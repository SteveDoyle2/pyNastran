from __future__ import annotations
from collections import defaultdict
from typing import TYPE_CHECKING

import numpy as np
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


def get_normals_at_elements(model: BDF) -> dict[int, np.ndarray]:
    normals = {}
    for eid, elem in model.elements.items():
        try:
            normal = elem.Normal()
        except AttributeError:
            continue
        normals[eid] = normal
    return normals

def get_normals_at_nodes(model: BDF) -> dict[int, np.ndarray]:
    normals_at_nodes = {}
    #node_count = defaultdict(int)
    nid_to_eid_map = defaultdict(list)

    normals = {}
    for eid, elem in model.elements.items():
        try:
            normal = elem.Normal()
        except AttributeError:
            continue
        normals[eid] = normal

        for nid in elem.nodes:
            #node_count[nid] += 1
            nid_to_eid_map[nid].append(eid)
    del normal, nid, eid

    for nid, eids in nid_to_eid_map.items():
        normal = np.zeros(3, dtype='float64')
        for eid in eids:
            normal += normals[eid]
        normal /= np.linalg.norm(normal)
        normals_at_nodes[nid] = normal

    return normals_at_nodes
