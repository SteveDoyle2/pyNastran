from collections import defaultdict
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.cards.elements.rigid import RBE2

def merge_rbe2s(model: BDF, ):
    eid_to_dependent_nodes = {}
    # eid_to_eids = {}
    dependent_nid_to_eids_map = defaultdict(set)
    nrigid = 0
    for eid, elem in model.rigid_elements.items():
        if elem.type not in {'RBE2'}:
            continue
        dependent_nodes = elem.dependent_nodes
        eid_to_dependent_nodes[eid] = dependent_nodes
        for nid in dependent_nodes:
            dependent_nid_to_eids_map[nid].add(eid)
        nrigid += 1
    dependent_nid_to_eids_map = dict(dependent_nid_to_eids_map)
    print(f'dependent_nid_to_eids_map = {dependent_nid_to_eids_map}')

    eid_to_nids_map = get_eid_to_nid_map(eid_to_dependent_nodes, dependent_nid_to_eids_map)
    eid_to_nids_map = get_eid_to_nid_map(eid_to_nids_map, dependent_nid_to_eids_map)
    eid_to_nids_map = get_eid_to_nid_map(eid_to_nids_map, dependent_nid_to_eids_map)
    print(f'eid_to_nids_map = {eid_to_nids_map}')

    nids_to_eids_map = defaultdict(list)
    for eid, nids in eid_to_nids_map.items():
        nids_list = list(nids)
        nids_list.sort()
        nids_to_eids_map[tuple(nids_list)].append(eid)
    print(f'nids_to_eids_map = {set(nids_to_eids_map)}')

    for dep_nids_tuple, eids in nids_to_eids_map.items():
        eid0 = eids[0]
        elem: RBE2 = model.rigid_elements[eid0]
        for eid in eids[1:]:
            del model.rigid_elements[eid]

        # TODO: update the CONM2 mass

        elem.Gmi = list(dep_nids_tuple)
        print(elem)
    return

def get_eid_to_nid_map(
        eid_to_nids_map: dict[int, set[int]],
        nid_to_eids_map: dict[int, set[int]],):
    eid_to_nids_map2 = defaultdict(set)
    for eid, nids in eid_to_nids_map.items():
        for nid in nids:
            eid_to_nids_map2[eid].add(nid)
            eids2 = nid_to_eids_map[nid]
            for eid2 in eids2:
                nids2 = eid_to_nids_map[eid2]
                eid_to_nids_map2[eid].update(nids2)
    return dict(eid_to_nids_map2)

def main():
    model = BDF()
    model.add_grid(1, [0., 0., 0.])
    model.add_grid(2, [0., 0., 0.])
    model.add_grid(3, [0., 0., 0.])
    model.add_grid(4, [0., 0., 0.])
    model.add_grid(5, [0., 0., 0.])

    model.add_grid(10, [0., 0., 0.])
    model.add_grid(11, [0., 0., 0.])
    model.add_grid(12, [0., 0., 0.])
    model.add_rbe2(1, 10, '123456', [1, 2])
    model.add_rbe2(2, 11, '123456', [2, 3])
    model.add_rbe2(3, 12, '123456', [4, 5])

    merge_rbe2s(model)

if __name__ == '__main__':
    main()
