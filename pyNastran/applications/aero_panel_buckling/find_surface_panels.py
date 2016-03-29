"""
Takes a 2D mesh of a structure (e.g. an aircraft) and finds the 2D skin panels
in order to do a buckling analysis.
"""
# pylint: disable=E1101
from __future__ import print_function
import os
import sys
import glob
from copy import deepcopy
from itertools import count
from six import iteritems, string_types

import numpy as np

from pyNastran.bdf.case_control_deck import CaseControlDeck
from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import OP2
from pyNastran.bdf.field_writer import print_card_8

def get_patch_edges(edge_to_eid_map, xyz_cid0, is_symmetric=True,
                    model=None, consider_pids=False):
    """
    Find all the edges of any patch

    Parameters
    ----------
    xyz_cid0 : (n, 3) ndarray
        nodes in cid=0
    eid_to_edge_map : dict[int] = [(int, int), (int, int), ...]
        maps element ids to edges
    is_symmetric : bool
        enables an xz symmetry (+/-y is a mirror)
        (should this actually be yz symmetry????)
        True: half model or model separated by a central gap
    model : BDF(); default=None
        required for consider_pids=True
    consider_pids : bool; default=False
        split model by properties

    Returns
    -------
    patch_edges : List[(int, int), (int, int), ...]
        the edges as sorted (int, int) pairs
    eids_on_edge : Set([int, int, int])
        the elements on the edges of patches
    free_edges : List[(int, int), (int, int), ...]
        the free edges of the model (will probably be removed)
    free_eids : Set([int, int, int])
        the free elements of the model (will probably be removed)
    """
    patch_edges = []
    eids_on_edge = set([])

    free_edges = []
    free_eids = set([])

    if not consider_pids:
        pids = []
    if is_symmetric:

        key0 = list(xyz_cid0.keys())[0]
        ymin = xyz_cid0[key0][0]
        for nid, xyz in iteritems(xyz_cid0):
            ymin = min(ymin, xyz[1])

        # xz symmetry  (should this actually be yz symmetry????)
        for edge, eids in iteritems(edge_to_eid_map):
            yedge = [True for nid in edge
                     if np.allclose(xyz_cid0[nid][1], ymin)]

            if consider_pids:
                pids = np.unique([model.elements[eid].Pid() for eid in eids])

            if len(eids) != 2 or len(yedge) > 1 or len(pids) > 1:
                patch_edges.append(edge)
                eids_on_edge.update(eids)
                if len(eids) == 1:
                    free_edges.append(edge)
                    free_eids.update(eids)
    else:
        for edge, eids in iteritems(edge_to_eid_map):
            if consider_pids:
                pids = np.unique([model.elements[eid].Pid() for eid in eids])

            if len(eids) != 2 or len(pids):
                patch_edges.append(edge)
                eids_on_edge.update(eids)
    return patch_edges, eids_on_edge, free_edges, free_eids

def find_ribs_and_spars(xyz_cid0, edge_to_eid_map, eid_to_edge_map,
                        workpath='results', is_symmetric=True,
                        model=None, consider_pids=False):
    """
    Extracts rib/spar/skin patches based on geometry.  You can have a
    single property id and this will still work.

    Parameters
    ----------
    xyz_cid0 : (n, 3) ndarray
        nodes in cid=0
    edge_to_eid_map : dict[(int, int)] = [int, int, ...]
        maps edges to element ids
    eid_to_edge_map : dict[int] = [(int, int), (int, int), ...]
        maps element ids to edges
    workpath : str; default='results'
        the results directory
    is_symmetric : bool; default=True
        enables yz symmetry
    model : BDF(); default=None
        required for consider_pids=True
    consider_pids : bool; default=False
        split model by properties

    Returns
    -------
    patch_edges_array : (n, 2) int ndarray
        all the edges on the geometry as (n1, n2) integer pairs
        where n1 < n2
    eids_on_edge : Set[int]
        all the elements on the edges
    patches : List[List[int]]
        the patches

    .. note :: We're only considering shell elements, which probably is a poor
               assumption given the stiffness of ring frames/longerons, which are
               typically modeled as CBAR/CBEAMs.

    Explanation of Approach
    =======================
    In the following case, we'll assume that line 6-7-8-9-10-11 is
    a hard spar (=) and 10-15 is a rib (#).  '|' is an element bounadary.

    1    2    3    4    5
    +----+----+----+----+
    |    |    |    |    |
    | A1 | B1 | C1 | D1 |    <------ patch 1
    |6   |7   |8   |9   |10  11
    +====+====+====+====+====+
    |    |    |    |    #    |
    | A2 | B2 | C2 | D2 # E3 |   <------ patch 2, 3
    |    |    |    |    #    |
    +----+----+----+----+----+
    11   12   13   14   15   16

    We'll start with the patch edges as:
        patch_edges = [
            # free edges
            (1, 2), (2, 3), (3, 4), (4, 5), # top
            (5, 10), (10, 15), (10, 11), (11, 16), # corner right
            (1, 6), (6, 11), # left
            (11, 12), (12, 13), (13, 14), (14, 15), (15, 16), # bottom

            # patch boundaries
            (6, 7), (7, 8), (8, 9), (9, 10), # spar
            (10, 15), # rib
        ]
    Note that node IDs tuples (n1, n2) have the property of n1 < n2.

    Given, an arbitrary starting edge, we can search for the neighboring
    elements and find all elements with a neighbor.  We have two cases:

    Free Edge
    =========
    1.  If we start from edge (1, 2), so we'll start at element A1 and
        identify the neighbors with no problems.

    Common Edge
    ===========
    1.  If we start from edge (9, 10), we'll find D1 and D2, which jumps
        patches.

    2.  However, had we started from edge (1, 2), removed patch 1, and
        then analyzed edge (9, 10), we'd wouldn't have jumped edges.

    Solution
    ========
    There are two ways to handle this problem (Common Edge #1):

    1.  For an asymmetrical model, you can calculate the free edges,
        start from a free edge, and update the free edges after
        finishing each patch.

    2.  For a symmetrical model, you have to handle patch jumping by
        only checking one of the possible edges for neighboring
        elements.  This gets tricky because some edges lead to isolated
        elements.  As such, you need to check to see if a valid patch was
        created.

    Note that option #2 is the general case and will probably be faster
    than option #1 because you don't need to caculate/recalculate the
    free edges ater each loop.
    """
    if not os.path.exists(workpath):
        os.makedirs(workpath)
    patch_edges, eids_on_edge, free_edges, free_eids = get_patch_edges(
        edge_to_eid_map, xyz_cid0, is_symmetric=is_symmetric,
        model=model, consider_pids=consider_pids)

    _free_edges = np.array(free_edges, dtype='int32')
    _free_edge_nids = np.unique(_free_edges.ravel())

    free_edge_nodes_filename = os.path.join(workpath, 'free_edge_nodes.csv')
    with open(free_edge_nodes_filename, 'w') as free_edge_nodes_file:
        free_edge_nodes_file.write('# X, Y, Z\n')
        free_edge_nodes_file.write('# Node IDs corresponding to free_edge_nodes_xyz.csv\n')
        np.savetxt(free_edge_nodes_file, _free_edge_nids, fmt='%i', delimiter=',')

    free_edge_nodes_xyz_filename = os.path.join(workpath, 'free_edge_nodes_xyz.csv')
    with open(free_edge_nodes_xyz_filename, 'w') as free_edge_nodes_xyz_file:
        for nid in _free_edge_nids:
            free_edge_nodes_xyz_file.write('%s, %s, %s\n' % tuple(xyz_cid0[nid]))

    #with open(free_edge_nodes_filename, 'w') as free_edge_nodes_file:
        #_free_edges = array(free_edges, dtype='int32')
        #for nid in unique(_free_edges.ravel()):
            #free_edge_nodes_file.write('%s, %s, %s\n' % tuple(xyz_cid0[nid]))

    patch_edges_filename = os.path.join(workpath, 'patch_edges.csv')
    patch_edges_array_filename = os.path.join(workpath, 'patch_edges_array.csv')
    _patch_edges = np.array(patch_edges, dtype='int32')
    header = 'Nid1, Nid2\n'
    header += 'NodeIDs corresponding to patch_edges_array.csv'
    np.savetxt(patch_edges_array_filename, _patch_edges, delimiter=',', header=header, fmt='%i')
    with open(patch_edges_filename, 'w') as patch_edges_file:
        patch_edges_file.write('# X, Y, Z\n')
        patch_edges_file.write('# xyz corresponding to patch_edges.csv\n')
        for nid in np.unique(_patch_edges.ravel()):
            patch_edges_file.write('%s, %s, %s\n' % tuple(xyz_cid0[nid]))

    patch_edges_array = np.array(patch_edges, dtype='int32')
    all_patch_edges = deepcopy(patch_edges)
    #--------------------------
    # we'll now get the patches
    patches = []
    patch_edge_count = len(patch_edges)
    assert patch_edge_count > 2, patch_edge_count

    # dummy value
    max_patches = 500
    icheck_max = 1000
    new_patch_edge_count = patch_edge_count - 1
    patch_count = -1
    new_patch_count = 0
    # loop until we're grouped all the patch edges into a patch
    # patches have one or more PSHELL/PCOMP properties
    while new_patch_count != patch_count:  # prevent infinite looping
        # get elements on other patches
        used_edges = set([])
        if patches:
            used_eids = set(np.unique(np.hstack(patches)))
        else:
            used_eids = set([])
        print('new_patch_edge_count=%s patch_edge_count=%s' % (
            new_patch_edge_count, patch_edge_count))
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

        # TODO: on only the first loop of the upper level while loop,
        #       downselect the edges to make sure that there are no
        #       edges that cross boundary lines.
        #eid0 = 97112
        edges_to_check_all = set(eid_to_edge_map[eid0])
        i = 0
        # see docstring for an explanation of this "big giant for loop"
        #for edge_to_check in edges_to_check_all:
            #edges_to_check = [edge_to_check]
        if 1:
            # this sometimes leads to patch joining

            if len(used_edges) == 0 and 0:
                # clean out patch edges if we're starting a new patch
                edges_to_check = []
                for edge in edges_to_check_all:
                    if edge in free_edges: # patch_edges?
                        continue
                    if edge in all_patch_edges: # patch_edges?
                        continue
                    if edge in used_edges: # patch_edges?
                        continue
                    edges_to_check.append(edge)
                edges_to_check = set(edges_to_check)
            else:
                edges_to_check = edges_to_check_all

            while len(edges_to_check):
                #print('len(edges_to_check) = %i' % len(edges_to_check))
                #if len(edges_to_check) < 5:
                    #print('  edges =', edges_to_check)

                edges_to_check_next, patch_len = get_next_edges(
                    eid_to_edge_map, edge_to_eid_map,
                    patch, all_patch_edges, patch_edges,
                    used_eids, used_edges,
                    edges_to_check, patch_len, free_edges)
                edges_to_check = edges_to_check_next
                i += 1
                if i == 1000:
                    raise RuntimeError('too many iterations')
                #print('len(edges_to_check) = %s' % (len(edges_to_check)))

            if len(patch):
                print('ipatch=%s len(patch)=%s; patch=%s' % (new_patch_count, len(patch), patch))
                sys.stdout.flush()
                new_patch_count += 1
                if len(patch) == 1:
                    break
                patches.append(list(patch))
                assert len(patch) > 1, len(patch)
            else:
                break
            print('*' * 80)

    return patch_edges_array, eids_on_edge, patches

def get_next_edges(eid_to_edge_map, edge_to_eid_map,
                   patch, all_patch_edges, patch_edges,
                   used_eids, used_edges,
                   edges_to_check, patch_len, free_edges):
    edges_to_check_next = deepcopy(edges_to_check)
    for edge in edges_to_check:
        #print('  edge =', edge)
        edges_to_check_next.remove(edge)
        #if edge in used_edges:
            #print('    continuing on used edge')
            #continue
        if edge in all_patch_edges:
            #print('    continuing on patch edge=%s' % str(edge))
            continue
        #if edge in free_edges:
            #print('    continuing on free edge=%s' % str(edge))
            #continue

        eids_to_check = [eid for eid in edge_to_eid_map[edge]
                         if eid not in used_eids]
        #print('  map = ', edge_to_eid_map[edge])
        #print('  eids_to_check[%s] = %s' % (edge, eids_to_check))
        #print('  used_eids =', used_eids)
        for eid in eids_to_check:
            if eid in used_eids:
                #print('  eid=%s is used' % eid)
                continue
            edgesi_to_check = eid_to_edge_map[eid]
            #print('  edgesi_to_check[%i] = %s' % (eid, edgesi_to_check))
            for edgei in edgesi_to_check:
                #print('  edgei =', edgei)
                if edgei in used_edges:
                    #print('    continuing on used edge=%s' % str(edgei))
                    continue
                if edgei in all_patch_edges:
                    #print('    continuing on patch edge=%s' % str(edgei))
                    continue
                #if edgei in free_edges:
                    #print('    continuing on free edge=%s' % str(edgei))
                    #continue
                used_edges.add(edgei)
                #edges_to_check_next.add(edgei)
                #edges_to_check_next.remove(edge)
                eidsi = edge_to_eid_map[edgei]

                eids2 = set([])
                for eidii in eidsi:
                    edges_temp = eid_to_edge_map[eidii]
                    #edges_temp2 = set([])
                    iedge_count = 0
                    for edge_temp in edges_temp:
                        #if edge_temp in used_edges:  # new
                            #continue
                        if edge_temp in all_patch_edges:  # new
                            continue
                        #if edge_temp in free_edges:  new
                            #continue
                        #eid0 = 97112
                        #assert (53480, 53481) != edge_temp, edges_temp
                        edges_to_check_next.add(edge_temp)
                        iedge_count += 1
                    if iedge_count:
                        eids2.add(eidii)
                used_eids.update(eids2)
                patch.update(eids2)

                #print('*    adding eids2=%s' % eids2)
                #eid0 = 97112
                #for eid_badi in [97063, 98749, 97062, 98748, 97061, 98747, 97060, 98746, 97059, 98745, 97058, 97090,
                                 #126619, 126620, 126621, 126622, 126623, 126624]:
                    #assert eid_badi not in patch, 'eid_badi=%s patch=%s' % (eid_badi, patch)

                patch_len += 1
        used_edges.add(edge)
        #print()
    return edges_to_check_next, patch_len

def save_patch_info(eids, xyz_cid0, patch_edges_array, eids_on_edge, patches,
                    all_eids=None,
                    workpath='results'):
    """
    Saves the patch data

    Parameters
    ----------
    eids : (nelements, ) int ndarray / List[int]
        list of element ids
    xyz_cid0 : dict[nid] = xyz
        key : int (node_id)
        value : (3, ) float ndarray of the grid point in cid=0
    eids_on_edge : Set[int]
        all the elements on the edges
    patches : List[List[int]]
        the patches
    patch_edges_array : (n, 2) int ndarray
        all the edges on the geometry as (n1, n2) integer pairs
        where n1 < n2
    all_eids : (mall_elements, ) int ndarray; default=None
        all the element ids to merge back in
        only necessary if we only analyzed a subset of elements
    workpath : str; default='results'
        the working directory
    """
    if isinstance(eids, np.ndarray):
        pass
    elif isinstance(eids, list):
        eids = np.asarray(eids)
    else:
        raise NotImplementedError('eids must be a list/ndarray; type=%s' % (type(eids)))

    #savetxt(patch_edges_array_filename, patch_edges_array, fmt='i', delimiter=',',
            #header='# patch edges array\n # n1, n2')
    #patch_edges_array_filename = os.path.join(workpath, 'patch_edges.csv')

    # is the element flagged as an edge
    element_edges_filename = os.path.join(workpath, 'element_edges.csv')

    # what patch is each element part of
    element_patches_filename = os.path.join(workpath, 'element_patches.csv')

    # what are the element ids in each patch
    patches_filename = os.path.join(workpath, 'patches.csv')

    # what elements are on the edge of a patch
    eids_on_edge_filename = os.path.join(workpath, 'eids_on_edge.txt')


    #with open(patch_edges_array_filename, 'w') as patch_edges_array_file:
        #for patch_edge in patch_edges_array:
            #patch_edges_array_file.write(str(patch_edge)[1:-1] + '\n')

    # xyz points of patch
    nodal_edges_filename = os.path.join(workpath, 'nodal_edges.csv')
    unique_nids = np.unique(np.array(patch_edges_array, dtype='int32').ravel())
    np.savetxt(nodal_edges_filename, [xyz_cid0[nid] for nid in unique_nids],
               header='xyz points for all patches', delimiter=',')

    patch_ids_all = np.ones(len(eids), dtype='int32') * -10
    npatches = len(patches)

    ipatch = 0
    for patch_eids in patches:
        if len(patch_eids) == 1:
            continue

        #for eid in patch_eids:
            #i = eids.index(eid)
            #patch_ids_all[i] = ipatch
        i = np.searchsorted(eids, patch_eids)
        patch_ids_all[i] = ipatch
        ipatch += 1

    with open(element_edges_filename, 'wb') as element_edges_file:
        element_edges_file.write('# eid(%i), is_edge(%i)\n')
        for eid in eids:
            if eid in eids_on_edge:
                element_edges_file.write('%s, 1\n' % eid)
            else:
                element_edges_file.write('%s, 0\n' % eid)

    with open(eids_on_edge_filename, 'wb') as eids_on_edge_file:
        for eid in eids_on_edge:
            eids_on_edge_file.write('%s\n' % eid)

    if len(eids) != len(patch_ids_all):
        msg = 'len(eids)=%i len(patch_ids_all)=%i' % (len(eids), len(patch_ids_all))
        raise ValueError(msg)
    with open(element_patches_filename, 'wb') as element_patches_file:
        element_patches_file.write('# eid(%i), patch(%i), element_count(%i)\n')
        for (i, eid, patch_id) in zip(count(), eids, patch_ids_all):
            neids_on_patch = len(patches[patch_id])
            element_patches_file.write('%s, %s, %s\n' % (eid, patch_id, neids_on_patch))

    with open(patches_filename, 'wb') as patches_file:
        patches_file.write('# patch_id(%i), eids(%i)\n')
        for patch_id, patch in enumerate(patches):
            patches_file.write('%s, ' % patch_id)
            patches_file.write(str(patch)[1:-1])
            patches_file.write('\n')
   # these should be the same if we did this right
    counter = [1 for patch in patches
               if len(patch) > 1]
    print('npatches=%s npids=%s' % (npatches, sum(counter)))


def load_patch_info(workpath='results'):
    """
    Loads the patch data

    Parameters
    ----------
    workpath : str; default='results'
        the working directory

    Returns
    -------
    eids_on_edge : Set[int]
        all the elements on the edges
    patches : List[List[int]]
        the patches
    patch_edges_array : (n, 2) int ndarray
        all the edges on the geometry as (n1, n2) integer pairs
        where n1 < n2
    """
    patch_edges_array_filename = os.path.join(workpath, 'patch_edges_array.csv')
    patch_edges_array = np.loadtxt(patch_edges_array_filename, delimiter=',', dtype='int32')


    eids_on_edge_filename = os.path.join(workpath, 'eids_on_edge.txt')
    eids_on_edge = set(np.loadtxt(eids_on_edge_filename, dtype='int32'))
    #nodal_edges_filename = os.path.join(workpath, 'nodal_edges.csv')
    #nodal_edges = load_nodal_edges_file(nodal_edges_filename)

    #element_edges_filename = os.path.join(workpath, 'element_edges.csv')
    #element_edges = load_element_edges(element_edges_filename)

    #element_patches_filename = os.path.join(workpath, 'element_patches.csv')
    #element_patches, _patch_counter = load_element_patches(element_patches_filename)
    #print('patches/element_patches =', element_patches)

    # this is patches
    patches_filename = os.path.join(workpath, 'patches.csv')
    patches = load_patches(patches_filename)
    return patch_edges_array, eids_on_edge, patches

#def load_nodal_edges_file(nodal_edges_filename):

#def load_element_edges(element_edges_filename):

def load_element_patches(element_patches_filename):
    A = np.loadtxt(element_patches_filename, delimiter=',', type='int32')
    element_patches = A[:, 1]
    counter = A[:, 2]
    return element_patches, counter

def load_patches(patches_filename):
    with open(patches_filename, 'r') as patches_file:
        lines = patches_file.readlines()

    patches = []
    for line in lines[1:]:
        sline = line.strip().split(',')
        vals = [int(val) for val in sline]
        patch_id = vals[0]
        patch = vals[1:]
        patches.append(patch)
    return patches

def create_plate_buckling_models(model, op2_filename, mode, isubcase=1,
                                 workpath='results', is_symmetric=True,
                                 rebuild_patches=True, consider_pids=False):
    """
    Create the input decks for buckling

    Parameters
    ----------
    model : BDF()
        a bdf object
    op2_filename : str
        the path to an OP2 file
    mode : str
        'load' : extract loads from the op2 and apply interface loads to
            each patch and constrain a boundary node
        'displacement' : extract displacements from the op2 and apply
            interface displacements for each patch
    is_symmetric : bool; default=True
        is the model symmetric (???)
        full model : set to False (???)
        half model : set to True (???)
    """
    # cleanup
    if workpath != '':
        print('workpath =', workpath)
        if not os.path.exists(workpath):
            os.makedirs(workpath)
        os.chdir(workpath)

    patch_filenames = glob.glob('patch_*.bdf')
    edge_filenames = glob.glob('edge_*.csv')
    for fname in patch_filenames + edge_filenames:
        os.remove(fname)

    #-------------------------
    assert mode in ['load', 'displacement'], 'mode=%r' % mode
    assert isinstance(op2_filename, string_types), 'op2_filename=%r' % op2_filename
    out = model._get_maps(consider_1d=False, consider_2d=True, consider_3d=False)
    (edge_to_eid_map, eid_to_edge_map, nid_to_edge_map) = out


    xyz_cid0 = {}
    #cds = {}
    for nid, node in iteritems(model.nodes):
        xyz_cid0[nid] = node.get_position()
        #cds[nid] = node.Cd()

    #rebuild_patches = False
    if rebuild_patches:
        patch_edges_array, eids_on_edge, patches = find_ribs_and_spars(
            xyz_cid0, edge_to_eid_map, eid_to_edge_map,
            workpath=workpath, is_symmetric=is_symmetric,
            model=model, consider_pids=consider_pids)

        eids = model.element_ids
        save_patch_info(eids, xyz_cid0, patch_edges_array, eids_on_edge, patches,
                        all_eids=None,
                        workpath=workpath)
    else:
        #msg = 'add ability to start from existing patch/edge files; rebuild_patches=%s' % rebuild_patches
        #raise NotImplementedError(msg)

        #workpath = 'results'
        patch_edges_array, eids_on_edge, patches = load_patch_info(workpath=workpath)
    write_buckling_bdfs(model, op2_filename, xyz_cid0,
                        patches, patch_edges_array,
                        subcase_id=isubcase, mode=mode, workpath=workpath)


def write_buckling_bdfs(model, op2_filename, xyz_cid0, patches, patch_edges_array,
                        subcase_id=1,
                        eig_min=-1.0, eig_max=3., nroots=20,
                        mode='displacement', workpath='results'):
    """
    Creates a series of buckling BDFs from a static analysis

    Parameters
    ----------
    model : BDF()
        the BDF object
    op2_filename : str
        the name of the OP2 to analyze
    xyz_cid0 : (nnodes, 3) ndarray
        the nodes in the global CID=0 frame
    patches : ???
        ???
    patch_edges_array : (n, 2) int ndarray
        all the edges on the geometry as (n1, n2) integer pairs
        where n1 < n2
    subcase_id : int; default=1
        the subcase to analyze
    mode : str; default='displacement'
        'displacement' : map displacements from the OP2 to the new BDFs
        'load' : map the grid point forces from the OP2 to the new BDFs
    workpath : str; default='results'
        the location of where the output files should be placed
    """
    # TODO: Create 'BDF' of deflected shape that doesn't include the wing up bend from the CAERO1
    #       boxes.  This will allow you to visually see where the wing undulates relative to the
    #       center plane of the wing
    assert mode in ['load', 'displacement'], 'mode=%r' % mode

    subcase = model.subcases[subcase_id]

    case_control_lines, bulk_data_cards, spc_id, mpc_id, load_id = create_buckling_header(subcase, eig_min, eig_max, nroots)

    out_model = OP2()
    print('**** workpath', workpath)
    out_model.read_op2(op2_filename)
    if mode == 'displacement':
        #print('out_model.displacements =', out_model.displacements)
        displacements = out_model.displacements[subcase_id]
        node_ids_full_model = displacements.node_gridtype[:, 0]
        ## TODO: check for cd != 0
    elif mode == 'load':
        nodal_forces = out_model.grid_point_forces[subcase_id]
    else:
        raise RuntimeError(mode)

    edge_nids = np.unique(patch_edges_array.ravel())
    #free_nodes = setdiff1d(node_ids_full_model, edge_nids) # A-B

    # TODO: support excluded pids
    excluded_pids = []

    patch_dir = os.path.join(workpath, 'patches')
    edge_dir = os.path.join(workpath, 'edges')
    if not os.path.exists(patch_dir):
        os.makedirs(patch_dir)
    if not os.path.exists(edge_dir):
        os.makedirs(edge_dir)
    for ipatch, patch_eids in enumerate(patches):
        patch_eids_array = np.array(patch_eids, dtype='int32')

        patch_filename = os.path.join(patch_dir, 'patch_%i.bdf' % ipatch)
        edge_filename = os.path.join(edge_dir, 'edge_%i.csv' % ipatch)
        patch_file = open(patch_filename, 'w')
        edge_file = open(edge_filename, 'w')

        patch_file.write('SOL 105\n')
        patch_file.write('CEND\n')
        patch_file.write(''.join(case_control_lines))
        patch_file.write('BEGIN BULK\n')
        for card in bulk_data_cards:
            patch_file.write(print_card_8(card))

        # get list of all nodes on patch_eids_array; write elements
        all_nids_on_patch = []
        pids_on_patch = set([])
        for eid in patch_eids_array:
            elem = model.elements[eid]
            node_ids = elem.node_ids

            all_nids_on_patch.append(node_ids)
            patch_file.write(str(elem))
            pids_on_patch.add(elem.Pid())

        # TODO: assumes a patch has only one property region
        pid = elem.Pid()

        if pid in excluded_pids:
            print('pid=%s is excluded; patch_filename=%s' % patch_filename)
            patch_file.close()
            edge_file.close()
            os.remove(patch_filename)
            os.remove(edge_filename)
            continue

        all_nids_on_patch = np.unique(np.hstack(all_nids_on_patch))

        nnodes = len(all_nids_on_patch)

        # get nodal cd; write nodes
        if mode == 'load':
            cd = np.zeros(nnodes, dtype='int32')
            for i, nid in enumerate(all_nids_on_patch):
                node = model.nodes[nid]
                cdi = node.Cd()
                patch_file.write(str(node))
                cd[i] = cdi
        elif mode == 'displacement':
            for i, nid in enumerate(all_nids_on_patch):
                node = model.nodes[nid]
                patch_file.write(str(node))
        else:
            raise NotImplementedError(mode)

        # model._write_common(patch_file, size=8, is_double=False)
        size = 8
        is_double = False
        model._write_coords(patch_file, size, is_double)
        model._write_materials(patch_file, size, is_double)

        for pid in pids_on_patch:
            prop = model.properties[pid]
            patch_file.write(prop.write_card(size=size, is_double=is_double))
            #model._write_properties(patch_file, size, is_double)

        # msg = ['                                          G R I D   P O I N T   F O R C E   B A L A N C E\n',
            # '  POINT-ID    ELEMENT-ID    SOURCE         T1             T2             T3             R1             R2             R3\n', ]
            #'0     13683          3736    TRIAX6         4.996584E+00   0.0            1.203093E+02   0.0            0.0            0.0'
            #'      13683          3737    TRIAX6        -4.996584E+00   0.0           -1.203093E+02   0.0            0.0            0.0'
            #'      13683                  *TOTALS*       6.366463E-12   0.0           -1.364242E-12   0.0            0.0            0.0'
        zero = ' '
        if 0:
            msg = []
            self = nodal_forces
            for ekey, force in sorted(iteritems(self.force_moment)):
                for iload, forcei in enumerate(force):
                    (f1, f2, f3, m1, m2, m3) = forcei
                    elem_name = self.elemName[ekey][iload]
                    eid = self.eids[ekey][iload]
                    # vals = [f1, f2, f3, m1, m2, m3]
                    # vals2 = write_floats_13e(vals)
                    # [f1, f2, f3, m1, m2, m3] = vals2
                    if eid == 0:
                        eid = ''
                    msg.append('%s  %8s    %10s    %-8s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                        zero, ekey, eid, elem_name, f1, f2, f3, m1, m2, m3))
                    zero = ' '
                zero = '0'

        if mode == 'displacement':
            # the displacement are in the Cd coordinate frame
            # because we're just carrying along the GRID cards
            # we shouldn't need to do anything
            patch_edge_nids = []
            for edge in patch_edges_array:
                n1, n2 = edge
                if n1 in all_nids_on_patch and n2 in all_nids_on_patch:
                    patch_edge_nids.append(edge)
            patch_edge_nids = np.unique(np.vstack(patch_edge_nids))
            #print('patch_edge_nids = %s' % patch_edge_nids)
            #print('node_ids_full_model = %s' % node_ids_full_model)
            idisp_full_model = np.searchsorted(node_ids_full_model, patch_edge_nids)
            nids2 = node_ids_full_model[idisp_full_model]

            #print('idisp_full_model = %s' % idisp_full_model)
            #print('nids2 = %s' % nids2)
            #print('delta = %s' % (node_ids_full_model[idisp_full_model] - patch_edge_nids))


            # make sure we're in Cd == 0
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

            spc1 = ['SPC1', spc_id, 123] + list(patch_edge_nids)
            #patch_file.write('SPC1,%i,123456,%i\n' % (spc_id, n1))
            patch_file.write(print_card_8(spc1))

            # dummy load
            # TODO: figure out free node
            #edge_nids = unique(patch_edges_array.ravel())
            #free_nodes = setdiff1d(node_ids_full_model, edge_nids) # A-B
            #free_nodes_on_patch = setdiff1d(free_nodes, all_nids_on_patch) # A-B
            #free_nodes_on_patch = in1d(free_nodes, all_nids_on_patch) # A and B

            free_nodes_on_patch = np.setdiff1d(all_nids_on_patch, edge_nids)
            #print('free_nodes_on_patch =', free_nodes_on_patch)
            if len(free_nodes_on_patch) == 0:
                print('couldnt find free node for patch_filename=%s' % patch_filename)
                patch_file.close()
                edge_file.close()
                #os.remove(patch_filename)
                #os.remove(edge_filename)
                continue
            free_node = free_nodes_on_patch[0]
            if free_node not in all_nids_on_patch:
                msg = 'free_node=%s is not patch_filename=in %s' % (free_node, patch_filename)
                raise RuntimeError(msg)
            patch_file.write('FORCE,%i,%i,0,0.000001,1.0,1.0,1.0\n' % (load_id, free_node))

            for i, nid in enumerate(patch_edge_nids):
                #if i == 0:
                    #continue
                xyz = xyz_cid0[nid]
                edge_file.write('%f, %f, %f\n' % (xyz[0], xyz[1], xyz[2]))
                assert nids2[i] == nid, 'nid=%s nid2=%s' % (nid, node_ids_full_model[i])
                for j in range(3):
                #for j in range(6):
                    dispi = disp[i, j]
                    if np.abs(dispi) > 0.0:
                        # SPCD, sid, g1, c1, d1
                        #patch_file.write(print_card_8(['SPCD', spc_id, nid, j + 1, dispi]))
                        if ipack == 0:
                            # SPCDs are loads, not SPCs.  Don't change this!!!
                            spc_pack = ['SPCD', load_id, nid, j + 1, dispi]
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
            # TODO: does this work for vectorized classes?
            for node_id, cdi in zip(all_nids_on_patch[1:], cd[1:]):
                force = nodal_forces.force_moment[node_id]
                force_moment_sum = np.zeros(6, dtype='float32')
                for iload, forcei in enumerate(force):
                    eid = nodal_forces.eids[node_id][iload]
                    element_name = nodal_forces.elemName[node_id][iload]
                    if eid not in patch_eids_array: # neighboring element
                        force_moment_sum += force[iload]
                    elif eid == 0 and element_name != '*TOTALS*':
                        print(element_name)
                        force_moment_sum += force[iload]

                abs_force_moment_sum = np.abs(force_moment_sum)
                forcei = abs_force_moment_sum[:3]
                momenti = abs_force_moment_sum[3:]
                if forcei.sum() > 0.0:
                    card = ['FORCE', load_id, node_id, cdi, 1.0, ] + list(forcei)
                    patch_file.write(print_card_8(card))
                if momenti.sum() > 0.0:
                    card = ['MOMENT', load_id, node_id, cdi, 1.0, ] + list(momenti)
                    patch_file.write(print_card_8(card))
        else:
            raise RuntimeError(mode)
        patch_file.write('ENDDATA\n')
        patch_file.close()
        edge_file.close()
    return


def write_buckling_bdf(model, op2_filename, nids_to_constrain,
                       bdf_filenames=None,
                       subcase_ids=None,
                       eig_min=-1.0, eig_max=1.1, nroots=100,
                       mode='displacement', size=8, is_double=False,
                       debug=False,
                       workpath='results', base_buckling='buckling_%i.bdf'):
    """
    Creates a single buckling BDF from a static analysis

    Parameters
    ----------
    model : BDF()
        the BDF object
        WARNING: the BDF object will be modified
    op2_filename : str
        the name of the OP2 to analyze
    nids_to_constrain : (nnodes, ) int ndarray
        the list of nodes IDs to constrain

    bdf_filenames : str; default=None -> ['buckling_%i.bdf']
        the name of the filename to write
        %i in default is the subcase_id
    subcase_ids : List[int]; default=None -> all subcases
        the subcase_ids to analyze
    mode : str; default='displacement'
        'displacement' : map displacements from the OP2 to the new BDFs
        'load' : map the grid point forces from the OP2 to the new BDFs (TODO: removed for now)
    workpath : str; default='results'
        the location of where the output files should be placed
    size : int; default=8
        the field width for the output BDF
    is_double : bool; default=False
        the large/double precision flag
    debug : bool; default=False
        the BDF/OP2 debug flag

    Returns
    -------
    bdf_filenames : List[str]
        the bdf filenames that were created
    subcase_ids : List[int]
        the subcase ids that were analyzed from the original BDF (not the new ones)
    """
    # TODO: Create 'BDF' of deflected shape that doesn't include the wing up bend from the CAERO1
    #       boxes.  This will allow you to visually see where the wing undulates relative to the
    #       center plane of the wing
    assert mode in ['displacement'], 'mode=%r' % mode
    if subcase_ids is None:
        subcase_ids = list(model.subcases.keys())
    elif isinstance(subcase_ids, int):
        subcase_ids = [subcase_ids]
    elif isinstance(subcase_ids, (list, tuple)):
        pass
    else:
        msg = 'subcase_ids=%s; expected List[int]; type=%s' % (subcase_ids, type(subcase_ids))
        raise TypeError(msg)
    subcase_ids = np.unique([subcase_id if subcase_id > 0 else 1 for subcase_id in subcase_ids])

    if bdf_filenames is None:
        bdf_filenames = []
        for isubcase in subcase_ids:
            bdf_filenames.append(base_buckling % isubcase)
    elif isinstance(bdf_filenames, int):
        bdf_filenames = [bdf_filenames]
    elif isinstance(bdf_filenames, (list, tuple)):
        pass
    else:
        msg = 'bdf_filenames=%s; expected List[str]; type=%s' % (bdf_filenames, type(bdf_filenames))
        raise TypeError(msg)

    if len(subcase_ids) != len(bdf_filenames):
        msg = 'nsubcases != nbdfs; nsubcases=%s nbdfs=%s; subcase_ids=%s bdf_filenames=%s' % (
            len(subcase_ids), len(bdf_filenames), subcase_ids, bdf_filenames)

    if not os.path.exists(workpath):
        os.makedirs(workpath)

    out_model = OP2(debug=debug)
    out_model.read_op2(op2_filename)

    case_control_deck = deepcopy(model.case_control_deck)
    print('**** workpath', workpath)

    if model.suport:
        # we don't need rigid body modes in buckling
        model.suport = []

    for subcase_id, bdf_filename in zip(subcase_ids, bdf_filenames):
        subcase = case_control_deck.subcases[subcase_id]

        model.sol = 105
        case_control_lines, bulk_data_cards, spc_id, mpc_id, load_id = create_buckling_header(subcase, eig_min, eig_max, nroots)
        if load_id in model.loads:
            del model.loads[load_id]
        if spc_id in model.spcs:
            del model.spcs[spc_id]


        dependent_nid_to_components = model.get_dependent_nid_to_components(mpc_id)

        if mode == 'displacement':
            displacements = out_model.displacements[subcase_id]
            node_ids_full_model = displacements.node_gridtype[:, 0]
        else:
            raise RuntimeError(mode)

        #edge_nids = np.unique(patch_edges_array.ravel())
        #free_nodes = setdiff1d(node_ids_full_model, edge_nids) # A-B


        model.case_control_deck = CaseControlDeck(case_control_lines)
        for card in bulk_data_cards:
            model.add_card(card, card[0])
        #nnodes = len(all_nids_on_patch)

        if mode == 'displacement':
            # the displacement are in the Cd coordinate frame
            # because we're just carrying along the GRID cards
            # we shouldn't need to do anything

            nids_extra = np.setdiff1d(nids_to_constrain, node_ids_full_model)
            assert len(nids_extra) == 0, 'Some nodes dont exist; nids_extra=%s' % str(nids_extra)
            idisp_full_model = np.searchsorted(node_ids_full_model, nids_to_constrain)
            #nids2 = node_ids_full_model[idisp_full_model]
            disp = displacements.data[0, idisp_full_model, :]

            ipack = 0
            spc_pack = []
            #$SPCD         100   20028       1.0009906   20028       24.3233-4
            #$SPCD         100   20028       3.3169889   20028       41.6345-4
            #$SPCD         100   20028       5.0004592   20028       6-2.092-3
            #FORCE,1,20036,0,0.000001,1.0,1.0,1.0

            free_nodes = np.setdiff1d(node_ids_full_model, nids_to_constrain)
            if len(free_nodes) == 0:
                raise RuntimeError('couldnt find any free_nodes')

            free_node = free_nodes[0]
            cid = 0
            force = ['FORCE', load_id, free_node, cid, 0.000001, 1.0, 1.0, 1.0]
            model.add_card(force, 'FORCE')

            spc1_1 = []
            spc1_2 = []
            spc1_3 = []
            for i, nid in enumerate(nids_to_constrain):
                dofs = ''
                if nid in dependent_nid_to_components:
                    dofs = dependent_nid_to_components[nid]
                for j in range(3):
                    dof = str(j + 1)
                    if dof in dofs:
                        print('skipping nid=%s dof=%s' % (nid, dof))
                        continue
                    if j == 0:
                        spc1_1.append(nid)
                    elif j == 1:
                        spc1_2.append(nid)
                    elif j == 2:
                        spc1_3.append(nid)
                    else:
                        raise NotImplementedError(j)

                    dispi = disp[i, j]
                    if np.abs(dispi) > 0.0:
                        # SPCD, sid, g1, c1, d1
                        if ipack == 0:
                            # SPCDs are loads, not SPCs.  Don't change this!!!
                            spc_pack = ['SPCD', load_id, nid, j + 1, dispi]
                            ipack += 1
                        else:
                            spc_pack += [nid, j + 1, dispi]
                            model.add_card(spc_pack, 'SPCD')
                            ipack = 0
                            spc_pack = []
            if len(spc_pack):
                model.add_card(spc_pack, 'SPCD')

            if len(spc1_1):
                spc1 = ['SPC1', spc_id, 1] + spc1_1
                model.add_card(spc1, 'SPC1')
            if len(spc1_2):
                spc1 = ['SPC1', spc_id, 2] + spc1_2
                model.add_card(spc1, 'SPC1')
            if len(spc1_3):
                spc1 = ['SPC1', spc_id, 3] + spc1_3
                model.add_card(spc1, 'SPC1')
        else:
            raise RuntimeError(mode)

        bdf_filename2 = os.path.join(workpath, bdf_filename)
        model.log.info('writing %s'% bdf_filename2)
        model.write_bdf(bdf_filename2, size=size, is_double=is_double)
    return bdf_filenames, subcase_ids


def create_buckling_header(subcase, eig_min=0., eig_max=100., nroots=20):
    # TODO: add title, subtitle from the actual load case for cross validation
    case_control = []
    #header += 'SOL 105\n'
    #header += 'CEND\n'
    case_control.append('ECHO = NONE\n')

    if 'TITLE' in subcase:
        case_control.append('  TITLE = %s\n' % subcase.get_parameter('TITLE')[0])
        #print(subcase.get_parameter('TITLE'))
    else:
        case_control.append('TITLE = BUCKLING\n')
    case_control.append('SUBCASE 1\n')

    keys_to_copy = ['SUBTITLE', 'LINE', 'MAXLINES', 'ECHO']
    for key in keys_to_copy:
        if key in subcase:
            line = '  %s = %s\n' % (key, subcase.get_parameter(key)[0])
            case_control.append(line)

    spc_id = 100
    if 'SPC' in subcase:
        spc_id = subcase.get_parameter('SPC')[0]
    case_control.append('  SPC = %i\n' % spc_id)

    mpc_id = None
    if 'MPC' in subcase:
        mpc_id = subcase.get_parameter('SPC')[0]
        case_control.append('  MPC = %i\n' % mpc_id)

    #case_control.append('  SUPORT = 1\n')  # no rigid body modes
    #if mode == 'load':
    case_control.append('  LOAD = 55\n')
    case_control.append('  METHOD = 42\n')
    #case_control.append('  STRESS(PLOT,PRINT,VONMISES,CENTER) = ALL\n')
    #case_control.append('  SPCFORCES(PLOT,PRINT) = ALL\n')
    #case_control.append('  STRAIN(PLOT,PRINT,VONMISES,FIBER,CENTER) = ALL\n')
    case_control.append('  DISPLACEMENT(PLOT,PRINT) = ALL\n')
    #header += '  TRIM = 1\n'

    bulk_data_cards = []
    bulk_data_cards.append(['PARAM', 'GRDPNT', 0])
    bulk_data_cards.append(['PARAM', 'COUPMASS', 1])
    #bulk_data_cards.append(['PARAM', 'AUNITS', 0.00259])
    #bulk_data_cards.append(['PARAM', 'WTMASS', 0.00259])
    #bulk_data_cards.append(['PARAM', 'BAILOUT', -1])
    bulk_data_cards.append(['PARAM', 'PRTMAXIM', 'YES'])
    bulk_data_cards.append(['PARAM', 'POST', -1])

    method = 42
    load_id = 55
    # EIGB has severe performance issue
    #bulk_data_cards.append(['EIGB', method, 'INV', eig_min, eig_max, nroots])
    bulk_data_cards.append(['EIGRL', eig_min, eig_max, nroots, ])
    return case_control, bulk_data_cards, spc_id, mpc_id, load_id


def find_surface_panels(bdf_filename=None, op2_filename=None, isubcase=1,
                        consider_pids=False, rebuild_patches=True,
                        workpath='results'):
    """prevents bleedover of data"""
    model = BDF()
    #bdf_filename = r'F:\work\pyNastran\pyNastran\master2\pyNastran\applications\mesh\surface\BWB_afl_static_analysis_short.bdf'
    #op2_filename = r'F:\work\pyNastran\pyNastran\master2\pyNastran\applications\mesh\surface\BWB_afl_static_analysis_short.op2'

    if bdf_filename is None:
        bdf_filename = 'model_144.bdf'
    if op2_filename is None:
        op2_filename = 'model_144.op2'

    #bdf_filename = r'F:\work\pyNastran\pyNastran\master2\models\solid_bending\solid_bending.bdf'
    if 0:
        if not os.path.exists(bdf_filename):
            bdf_filename = os.path.join(os.path.expanduser('~'),
                                        'Desktop', 'move', '3_LoadCases_Final',
                                        'BWB_afl_static_analysis_short.bdf')
            op2_filename = os.path.join(os.path.expanduser('~'),
                                        'Desktop', 'move', '3_LoadCases_Final',
                                        'BWB_afl_static_analysis_short.op2')
            assert os.path.exists(bdf_filename), bdf_filename

    if not os.path.exists('model.obj') or 1:
        model.read_bdf(bdf_filename, xref=False)
        #model.save('model.obj')
    else:
        model.load('model.obj')
    model.cross_reference()

    #create_plate_buckling_models(model, op2_filename, 'load')
    create_plate_buckling_models(model, op2_filename, 'displacement',
                                 workpath=workpath, isubcase=isubcase,
                                 consider_pids=consider_pids,
                                 rebuild_patches=rebuild_patches)


if __name__ == '__main__':
    find_surface_panels()
