"""
Takes a 2D mesh of a structure (e.g. an aircraft) and finds the 2D skin panels
in order to do a buckling analysis.
"""
from __future__ import print_function
import os
from copy import deepcopy
from six import iteritems

from numpy import hstack, unique, allclose, savetxt, array
from numpy import zeros, abs
from numpy import searchsorted, vstack

from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import OP2
from pyNastran.bdf.fieldWriter import print_card_8


def find_ribs_and_spars(xyz_cid0, edge_to_eid_map, eid_to_edge_map):
    """
    Parameters
    ----------
    xyz_cid0 : (n, 3) ndarray
        nodes in cid=0
    edge_to_eid_map : dict[(int, int)] = [int, int, ...]
        maps edges to element ids
    eid_to_edge_map : dict[int] = [(int, int), (int, int), ...]
        maps element ids to edges

    Returns
    -------
    patch_edges_array : (n, 3) ndarray
        all the edges on the geometry
    eids_on_edge : Set[int]
        all the elements on the edges
    patches : List[List[int]]
        the patches
    """

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
        # old_patch_len = 0

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

def save_patch_info(model, xyz_cid0, patch_edges, eids_on_edge, patches):
    """saves the patch data"""
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

    element_edges_file = open('element_edges.txt', 'wb')
    element_edges_file.write('# is_edge\n')
    for eid in eids_all:
        if eid in eids_on_edge:
            element_edges_file.write('1\n')
        else:
            element_edges_file.write('0\n')
    element_edges_file.close()

    element_patches_file = open('element_patches.txt', 'wb')
    element_patches_file.write('# patch\n')
    for pid in pids_all:
        element_patches_file.write('%s\n' % pid)
    element_patches_file.close()

   # these should be the same if we did this right
    counter = [1 for patch in patches
               if len(patch) > 1]
    print('npatches=%s npids=%s' % (npatches, sum(counter)))


def run(model, op2_filename, mode):
    """
    create the input decks for buckling
    """
    out = model._get_maps(consider_1d=False, consider_2d=True, consider_3d=False)
    (edge_to_eid_map, eid_to_edge_map, nid_to_edge_map) = out
    xyz_cid0 = {}
    for nid, node in iteritems(model.nodes):
        xyz_cid0[nid] = node.get_position()
    patch_edges, eids_on_edge, patches = find_ribs_and_spars(xyz_cid0, edge_to_eid_map, eid_to_edge_map)

    save_patch_info(model, xyz_cid0, patch_edges, eids_on_edge, patches)

    header = ''
    header += 'SOL 105\n'
    header += 'CEND\n'
    header += 'TITLE = BUCKLING\n'
    header += 'ECHO = NONE\n'
    header += 'SUBCASE 1\n'
    # header += '  SUPORT = 1\n'
    if mode == 'load':
        header += '  LOAD = 55\n'
    header += '  METHOD = 42\n'
    header += '  STRESS(PLOT,PRINT,VONMISES,CENTER) = ALL\n'
    header += '  SPCFORCES(PLOT,PRINT) = ALL\n'
    header += '  STRAIN(PLOT,PRINT,VONMISES,FIBER,CENTER) = ALL\n'
    header += '  DISPLACEMENT(PLOT,PRINT) = ALL\n'
    header += '  SPC = 100\n'
    #header += '  MPC = 1\n'
    # header += '  TRIM = 1\n'
    header += 'BEGIN BULK\n'
    header += 'PARAM,GRDPNT,0\n'
    header += 'PARAM,COUPMASS,1\n'
    header += 'PARAM,AUNITS,0.00259\n'
    header += 'PARAM,WTMASS,0.00259\n'
    header += 'PARAM,BAILOUT,-1\n'
    header += 'PARAM,PRTMAXIM,YES\n'
    header += 'PARAM,POST,-1    \n'

    eig1 = 0.0
    eig2 = 100.
    nroots = 20
    method = 42
    spc_id = 100
    load_id = 55
    header += 'EIGB,%s,INV,%s,%s,%s\n' % (method, eig1, eig2, nroots)

    isubcase = 1
    out_model = OP2()
    out_model.read_op2(op2_filename)
    if mode == 'displacement':
        displacements = out_model.displacements[isubcase]
        node_ids_full_model = displacements.node_gridtype[:, 0]
    elif mode == 'load':
        nodal_forces = out_model.grid_point_forces[isubcase]
    else:
        raise RuntimeError(mode)

    for ipatch, patch in enumerate(patches):
        patch = array(patch, dtype='int32')
        eids = patch
        patch_file = open('patch_%i.bdf' % ipatch, 'w')
        edge_file = open('edge_%i.bdf' % ipatch, 'w')
        patch_file.write(header)

        # get list of all nodes on patch; write elements
        all_nids = []
        for eid in eids:
            elem = model.elements[eid]
            node_ids = elem.node_ids
            all_nids.append(node_ids)
            patch_file.write(str(elem))
        all_nids = unique(hstack(all_nids))
        nnodes = len(all_nids)

        # get nodal cd; write nodes
        cd = zeros(nnodes, dtype='int32')
        for i, nid in enumerate(all_nids):
            node = model.nodes[nid]
            cdi = node.Cd()
            patch_file.write(str(node))
            cd[i] = cdi

        # model._write_common(patch_file, size=8, is_double=False)
        size = 8
        is_double = False
        model._write_coords(patch_file, size, is_double)
        model._write_materials(patch_file, size, is_double)
        model._write_properties(patch_file, size, is_double)

        # msg = ['                                          G R I D   P O I N T   F O R C E   B A L A N C E\n',
            # '  POINT-ID    ELEMENT-ID    SOURCE         T1             T2             T3             R1             R2             R3\n', ]
            #'0     13683          3736    TRIAX6         4.996584E+00   0.0            1.203093E+02   0.0            0.0            0.0'
            #'      13683          3737    TRIAX6        -4.996584E+00   0.0           -1.203093E+02   0.0            0.0            0.0'
            #'      13683                  *TOTALS*       6.366463E-12   0.0           -1.364242E-12   0.0            0.0            0.0'
        zero = ' '
        if 0:
            msg = []
            self = nodal_forces
            for ekey, Force in sorted(iteritems(self.force_moment)):
                for iLoad, force in enumerate(Force):
                    (f1, f2, f3, m1, m2, m3) = force
                    elem_name = self.elemName[ekey][iLoad]
                    eid = self.eids[ekey][iLoad]
                    # vals = [f1, f2, f3, m1, m2, m3]
                    # (vals2, is_all_zeros) = writeFloats13E(vals)
                    # [f1, f2, f3, m1, m2, m3] = vals2
                    if eid == 0:
                        eid = ''
                    msg.append('%s  %8s    %10s    %-8s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (zero, ekey, eid, elem_name,
                                                                                                         f1, f2, f3, m1, m2, m3))
                    zero = ' '
                zero = '0'

        if mode == 'displacement':
            patch_edge_nids = []
            for edge in patch_edges:
                n1, n2 = edge
                if n1 in all_nids and n2 in all_nids:
                    patch_edge_nids.append(edge)
            patch_edge_nids = unique(vstack(patch_edge_nids))
            #print('patch_edge_nids = %s' % patch_edge_nids)
            #print('node_ids_full_model = %s' % node_ids_full_model)
            idisp_full_model = searchsorted(node_ids_full_model, patch_edge_nids)
            nids2 = node_ids_full_model[idisp_full_model]
            #print('idisp_full_model = %s' % idisp_full_model)
            #print('nids2 = %s' % nids2)
            #print('delta = %s' % (node_ids_full_model[idisp_full_model] - patch_edge_nids))
            disp = displacements.data[0, idisp_full_model, :]
            #disp = displacements.data[:, 3]

            ipack = 0
            spc_pack = []

            #$SPCD         100   20028       1.0009906   20028       24.3233-4
            #$SPCD         100   20028       3.3169889   20028       41.6345-4
            #$SPCD         100   20028       5.0004592   20028       6-2.092-3
            #FORCE,1,20036,0,0.000001,1.0,1.0,1.0
            #SPC1,100,123456,20028
            n1 = patch_edge_nids[0]
            n2 = patch_edge_nids[1]
            patch_file.write('SPC1,%i,123456,%i\n' % (spc_id, n1))

            # dummy load
            patch_file.write('FORCE,%i,%i,0,0.000001,1.0,1.0,1.0\n' % (load_id, n2))
            for i, nid in enumerate(patch_edge_nids):
                if i == 0:
                    continue
                xyz = xyz_cid0[nid]
                edge_file.write('%f, %f, %f\n' % (xyz[0], xyz[1], xyz[2]))
                assert nids2[i] == nid, 'nid=%s nid2=%s' % (nid, node_ids_full_model[i])
                for j in range(6):
                    dispi = disp[i, j]
                    if abs(dispi) > 0.0:
                        # SPCD, sid, g1, c1, d1
                        #patch_file.write(print_card_8(['SPCD', spc_id, nid, j + 1, dispi]))
                        if ipack == 0:
                            spc_pack = ['SPCD', spc_id, nid, j + 1, dispi]
                            ipack += 1
                        else:
                            # ipack = 1
                            spc_pack += [nid, j + 1, dispi]
                            patch_file.write(print_card_8(spc_pack))
                            ipack = 0
                            spc_pack = []
            if len(spc_pack):
                patch_file.write(print_card_8(spc_pack))

        elif mode == 'load':
            for node_id, cdi in zip(all_nids[1:], cd[1:]):
                Force = nodal_forces.force_moment[node_id]
                force_moment_sum = zeros(6, dtype='float32')
                for iload, force in enumerate(Force):
                    eid = nodal_forces.eids[node_id][iload]
                    element_name = nodal_forces.elemName[node_id][iload]
                    if eid not in eids: # neighboring element
                        force_moment_sum += Force[iload]
                    elif eid == 0 and element_name != '*TOTALS*':
                        print(element_name)
                        force_moment_sum += Force[iload]

                abs_force_moment_sum = abs(force_moment_sum)
                forcei = abs_force_moment_sum[:3]
                momenti = abs_force_moment_sum[3:]
                if forcei.sum() > 0.0:
                    # write force
                    card = ['FORCE', load_id, node_id, cdi, 1.0, ] + list(forcei)
                    patch_file.write(print_card_8(card))
                if momenti.sum() > 0.0:
                    # write moment
                    card = ['MOMENT', load_id, node_id, cdi, 1.0, ] + list(momenti)
                    patch_file.write(print_card_8(card))
        else:
            raise RuntimeError(mode)
        patch_file.write('ENDDATA\n')
        patch_file.close()
        edge_file.close()
    return

def main():
    """prevents bleedover of data"""
    model = BDF()
    bdf_filename = r'F:\work\pyNastran\pyNastran\master2\pyNastran\applications\mesh\surface\BWB_afl_static_analysis_short.bdf'
    op2_filename = r'F:\work\pyNastran\pyNastran\master2\pyNastran\applications\mesh\surface\BWB_afl_static_analysis_short.op2'
    #bdf_filename = r'F:\work\pyNastran\pyNastran\master2\models\solid_bending\solid_bending.bdf'
    if not os.path.exists(bdf_filename):
        bdf_filename = os.path.join(os.path.expanduser('~'),
                                    'Desktop', 'move', '3_LoadCases_Final',
                                    'BWB_afl_static_analysis_short.bdf')
        op2_filename = os.path.join(os.path.expanduser('~'),
                                    'Desktop', 'move', '3_LoadCases_Final',
                                    'BWB_afl_static_analysis_short.op2')
        assert os.path.exists(bdf_filename), bdf_filename

    if not os.path.exists('model.obj') or 1:
        model.read_bdf(bdf_filename)
        # model.save_object('model.obj')
    else:
        model.load_object('model.obj')
        model.cross_reference()
    #run(model, op2_filename, 'load')
    run(model, op2_filename, 'displacement')


if __name__ == '__main__':
    main()
