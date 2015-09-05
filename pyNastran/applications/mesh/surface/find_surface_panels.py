from pyNastran.bdf.bdf import BDF


def find_ribs_and_spars(model, edge_to_eid_map, eid_to_edge_map):
    # find all the edges of any patch
    patch_edges = []
    for edge, eids in iteritems(edge_to_eid_map):
        if len(eids) != 2:
            patch_edges.append(edge)

    #--------------------------
    # we'll not get the patches
    patches = []
    used_eids = set([])
    patch_edge_count = len(patch_edges)
    assert patch_edge_count > 2, patch_edge_count

    # dummy value
    new_patch_edge_count = patch_edge_count - 1

    # loop until we're grouped all the patch edges into a patch
    # patches have one or more PSHELL/PCOMP properties
    while new_patch_edge_count != patch_edge_count:  # prevent infinite looping
        edge0 = patch_edges.pop()
        eids = edge_to_eid_map[edge0]
        eid0 = eids.pop()
        patch = set([eid0])
        used_eids.add(eid0)
        len_patch = 1
        old_patch_len = 0

        edges_to_check = [edge0]
        while len(edges_to_check):
            for edge in edges_to_check:
                eids_to_check = [eid for eid in edge_to_eid_map[edge]
                                 if eid not in used_eids]

                for eid in eids_to_check:
                    for edge in eid_to_edge_map[eid0]:
                        if edge in used_edges:
                            continue
                        if edge not in patch_edges:
                            eid1 = eids.pop()
                            edges_to_check.append(edge)
                            used_eids.add(eid1)
                            patch.add(eid1)
                            patch_len += 1
        patches.append(patch)
    return patches

def main():
    model = BDF()
    model.read_bdf(bdf_filename)
    out = model._get_maps()

    (edge_to_eid_map, eid_to_edge_map, nid_to_edge_map) = out
    find_ribs_and_spars(model, edge_to_eid_map, eid_to_edge_map)


if __name__ == '__main__':
    main()