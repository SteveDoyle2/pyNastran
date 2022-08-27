import copy
from itertools import combinations
from typing import List # , TYPE_CHECKING

import numpy as np
#if TYPE_CHECKING:
from pyNastran.bdf.bdf import BDF, read_bdf

def break_elements(model: BDF, groups, idtype: str='int32') -> BDF:
    """
    when elements overlap, we can split them to create new nodes

    1     2      3       1     2  7      3
    +-----+------+       +-----+  +------+
    |  1  |   2  |   ->  |  1  |  |   2  |
    +-----+------+       +-----+  +------+
    4     5      6       4     5  8      6

    1     2      3       1    10  11     3
    +-----+------+       +-----+  +------+
    |  1  |   2  |   ->  |  1  |  |   2  |
    4-----5------6       +-----+  +------+
    |  3  |   3  |       12    13 14     15
    +-----+------+
    7     8      9       16     17      18
                         +------+--------+
                         |   3  |   3    |
                         +------+--------+
                         4      8        9
    group_1 = [1]
    group_2 = [2]
    groups = [group_1, group_2]
    break_elements(model, groups)

    group_1 = [1]
    group_2 = [2]
    group_3 = [3, 4]
    groups = [group_1, group_2, group_3]
    break_elements(model, groups)

    # or more efficiently???
    groups = [group_3, group_1, group_2]
    break_elements(model, groups)

    """
    ngroups = len(groups)
    groupsi = range(ngroups)
    combos = combinations(groupsi, 2)

    if isinstance(model, str):
        model = read_bdf(model, xref=False)
    else:
        assert isinstance(model, BDF)
    node_id = max(model.nodes) + 1

    #nodes_by_group = {}
    #for igroup, group in enumerate(groups):

    mapped_nodes = {}
    for pair in combos:
        igroup0, igroup1 = pair
        nids_group0 = get_nodes_by_group(model, groups[igroup0], idtype=idtype)
        nids_group1 = get_nodes_by_group(model, groups[igroup1], idtype=idtype)
        overlapping_nids = np.intersect1d(nids_group0, nids_group1)
        if len(overlapping_nids) == 0:
            print('overlapping_nids', overlapping_nids)
            continue

        print(f'noverlapping_nids = {len(overlapping_nids)}')
        for nid in overlapping_nids:
            mapped_nodes[nid] = node_id
            node_id += 1

        # we will update group2
        for eid in groups[igroup1]:
            try:
                elem = model.elements[eid]
            except KeyError:
                continue
            nodes = elem.nodes
            for i, nid in enumerate(nodes):
                mapped_nid = mapped_nodes.get(nid, nid)
                nodes[i] = mapped_nid

    for nid, mapped_nid in mapped_nodes.items():
        node = model.nodes[nid]
        node2 = copy.deepcopy(node)
        node2.nid = mapped_nid
        model.nodes[mapped_nid] = node2
    return model

def get_nodes_by_group(model: BDF, group: list[int], idtype='int32'):
    nodes_by_group = []
    missing_eids = []
    for eid in group:
        try:
            elem = model.elements[eid]
        except KeyError:
            missing_eids.append(eid)
            continue

        nids = elem.nodes
        nodes_by_group.extend(nids)
    if missing_eids:
        model.log.warning(f'missing elements to break...eids={missing_eids}')
    return np.unique(nodes_by_group)


