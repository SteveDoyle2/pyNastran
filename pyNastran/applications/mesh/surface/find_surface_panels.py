from __future__ import print_function
import os
from copy import deepcopy
from six import iteritems
from pyNastran.bdf.bdf import BDF
from numpy import hstack, unique, allclose, savetxt, vstack, array


def find_ribs_and_spars(model, xyz_cid0, edge_to_eid_map, eid_to_edge_map):

    # find all the edges of any patch
    patch_edges = []
    eids_on_edge = set([])
    for edge, eids in iteritems(edge_to_eid_map):
        yedge = [True for nid in edge
                 if allclose(xyz_cid0[nid][0], 0.0)]

        if len(eids) != 2 or len(yedge) > 1:
            patch_edges.append(edge)
            eids_on_edge.update(eids)
    patch_edges_array = array(patch_edges, dtype='int32')
    #--------------------------
    # we'll not get the patches
    patches = []
    patch_edge_count = len(patch_edges)
    assert patch_edge_count > 2, patch_edge_count

    # dummy value
    new_patch_edge_count = patch_edge_count - 1
    patch_count = -1
    new_patch_count = 0
    # loop until we're grouped all the patch edges into a patch
    # patches have one or more PSHELL/PCOMP properties
    while new_patch_count != patch_count:  # prevent infinite looping
        # get elements on other patches
        used_edges = set([])
        if patches:
            used_eids = set(unique(hstack(patches)))
        else:
            used_eids = set([])
        print('new_patch_edge_count=%s patch_edge_count=%s' % (new_patch_edge_count, patch_edge_count))
        patch_edge_count = new_patch_edge_count
        if len(patch_edges) == 0:
            break

        # We'll start from an edge and move inwards by finding an element on
        # the edge to start from.  We'll assume all edges from an element
        # associated with a bounding edge are starting edges.
        #
        # Then we'll use a different edge from some an element with that edge
        # and check for a non-patch-edge edge.
        #
        # We'll scoop up new elements and make sure we haven't checked elements
        # twice using used_eids
        edge0 = patch_edges.pop()
        eids = edge_to_eid_map[edge0]
        for eid0 in eids:
            if eid0 not in used_eids:
                break
        else:
            new_patch_edge_count += 1
            continue
        print('*eid0 =', eid0)
        patch = set([eid0])
        used_eids.add(eid0)
        patch_len = 1
        old_patch_len = 0

        edges_to_check = set(eid_to_edge_map[eid0])
        i = 0
        while len(edges_to_check):
            print('len(edges_to_check) = %i' % len(edges_to_check))
            if len(edges_to_check) < 5:
                print('  edges =', edges_to_check)

            edges_to_check_next = deepcopy(edges_to_check)
            for edge in edges_to_check:
                print('  edge =', edge)
                edges_to_check_next.remove(edge)
                if edge in patch_edges:
                    print('    continuing on patch edge')
                    continue

                eids_to_check = [eid for eid in edge_to_eid_map[edge]
                                 if eid not in used_eids]
                print('  map = ', edge_to_eid_map[edge])
                print('  eids_to_check[%s] = %s' % (edge, eids_to_check))
                #print('  used_eids =', used_eids)
                for eid in eids_to_check:
                    if eid in used_eids:
                        print('  eid=%s is used' % eid)
                        continue
                    edgesi_to_check = eid_to_edge_map[eid]
                    print('  edgesi_to_check[%i] = %s' % (eid, edgesi_to_check))
                    for edgei in edgesi_to_check:
                        #print('  edgei =', edgei)
                        if edgei in used_edges:
                            print('    edge=%s is used' % str(edgei))
                            continue
                        if edgei in patch_edges:
                            print('    continuing on patch edge=%s' % str(edgei))
                            continue
                        used_edges.add(edgei)
                        #edges_to_check_next.add(edgei)
                        #edges_to_check_next.remove(edge)
                        eidsi = edge_to_eid_map[edgei]
                        for eidii in eidsi:
                            edges_to_check_next.update(eid_to_edge_map[eidii])
                        used_eids.update(eidsi)
                        patch.update(eidsi)
                        print('*    adding eids=%s' % eidsi)
                        patch_len += 1
                used_edges.add(edge)
                print()

            edges_to_check = edges_to_check_next
            i += 1
            if i == 1000:
                raise RuntimeError('too many iterations')
            print('len(edges_to_check) = %s' % (len(edges_to_check)))
        if len(patch):
            print('patch=%s; len(patch)=%s' % (patch, len(patch)))
            new_patch_count += 1
            patches.append(list(patch))
        else:
            break
        print('*' * 80)
    return patch_edges_array, eids_on_edge, patches

def main():
    model = BDF()
    if not os.path.exists('model.obj') or 0:
        bdf_filename = r'F:\work\pyNastran\set1_test\BWB_afl_static_analysis_short.bdf'
        #bdf_filename = r'F:\work\pyNastran\pyNastran\master2\models\solid_bending\solid_bending.bdf'
        model.read_bdf(bdf_filename)
        model.save_object('model.obj')
    else:
        model.load_object('model.obj')
        model.cross_reference()

    out = model._get_maps(consider_1d=False, consider_2d=True, consider_3d=False)
    (edge_to_eid_map, eid_to_edge_map, nid_to_edge_map) = out
    xyz_cid0 = {}
    for nid, node in iteritems(model.nodes):
        xyz_cid0[nid] = node.get_position()
    patch_edges, eids_on_edge, patches = find_ribs_and_spars(model, xyz_cid0, edge_to_eid_map, eid_to_edge_map)
    unique_nids = unique(array(patch_edges, dtype='int32').ravel())
    print('unique_nids =', unique_nids)

    savetxt('nodal_edges.txt', [xyz_cid0[nid] for nid in unique_nids], delimiter=',')

    eids_all = model.element_ids
    pids_all = [-10] * len(eids_all)
    npatches = len(patches)

    ipatch = 1
    for eids in patches:
        if len(eids) == 1:
            continue
        for eid in eids:
            i = eids_all.index(eid)
            pids_all[i] = ipatch
        ipatch += 1

    g = open('element_edges.txt', 'wb')
    g.write('# is_edge\n')
    for eid in eids_all:
        if eid in eids_on_edge:
            g.write('1\n')
        else:
            g.write('0\n')
    g.close()

    f = open('element_patches.txt', 'wb')
    f.write('# patch\n')
    for pid in pids_all:
        f.write('%s\n' % pid)
    f.close()

   # these should be the same if we did this right
    counter = [1 for patch in patches
               if len(patch) > 1]
    print('npatches=%s npids=%s' % (npatches, sum(counter)))

if __name__ == '__main__':
    main()