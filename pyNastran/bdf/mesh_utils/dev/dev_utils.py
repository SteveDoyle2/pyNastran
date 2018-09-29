"""
defines some half-implemented functions with significant restrictions
 - find_closest_nodes(...)
 - create_spar_cap(...)
 - split_model_by_material_id(...)
 - extract_surface_patches(...)
 - get_joints(...)
 - cut_model(...)
 - get_free_edges(...)
"""
from __future__ import print_function
from collections import defaultdict
from functools import reduce

from six import string_types


import numpy as np
from numpy import (array, where, hstack, searchsorted, float32, arccos, dot, degrees)

from numpy.linalg import norm  # type: ignore

from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.cards.nodes import GRID
from pyNastran.bdf.cards.elements.shell import CQUAD4
from pyNastran.bdf.cards.elements.rigid import RBE3
from pyNastran.bdf.mesh_utils.bdf_equivalence import (
    _eq_nodes_setup, _eq_nodes_build_tree, _eq_nodes_find_pairs)
from pyNastran.bdf.mesh_utils.find_closest_nodes import (
    find_closest_nodes, find_closest_nodes_index)


def create_rbe3s_between_close_nodes(bdf_filename, bdf_filename_out, tol,
                                     renumber_nodes=False, neq_max=4, xref=True,
                                     node_set=None, size=8, is_double=False,
                                     debug=True):
    """
    Creates semi-rigid RBE3 elements between two nodes within tolerance.

    Parameters
    ----------
    bdf_filename : str / BDF
        str : bdf file path
        BDF : a BDF model that is fully valid (see xref)
    bdf_filename_out : str
        a bdf_filename to write
    tol : float
        the spherical tolerance
    renumber_nodes : bool
        should the nodes be renumbered (default=False)
    neq_max : int
        the number of "close" points (default=4)
    xref bool: bool
        does the model need to be cross_referenced
        (default=True; only applies to model option)
    node_set : List[int] / (n, ) ndarray
        the list/array of nodes to consider (not supported with renumber_nodes=True)
    size : int; {8, 16}; default=8
        the bdf write precision
    is_double : bool; default=False
        the field precision to write

    Returns
    -------
    model : BDF()
        The BDF model corresponding to bdf_filename_out

    .. warning:: I doubt SPOINTs/EPOINTs work correctly
    .. warning:: xref not fully implemented (assumes cid=0)

    .. todo:: node_set stil does work on the all the nodes in the big
               kdtree loop, which is very inefficient
    """
    nodes_xyz, model, nids, inew = _eq_nodes_setup(
        bdf_filename, tol, renumber_nodes=renumber_nodes,
        xref=xref, node_set=node_set, debug=debug)

    eid = 1
    if len(model.rigid_elements):
        eid = max(model.rigid_elements) + 1

    ieq, slots = _eq_nodes_build_tree(nodes_xyz, nids, tol,
                                      inew=inew, node_set=node_set,
                                      neq_max=neq_max)[1:]

    nid_pairs = _eq_nodes_find_pairs(nids, slots, ieq, node_set=node_set)

    for (nid1, nid2) in nid_pairs:
        node1 = model.nodes[nid1]
        node2 = model.nodes[nid2]

        # TODO: doesn't use get position...
        distance = norm(node1.xyz - node2.xyz)
        if distance > tol:
            continue
        refgrid = nid1
        refc = '123456'
        weights = [1.0]
        comps = ['123456']
        gijs = [[nid2]]
        rbe3 = RBE3(eid, refgrid, refc, weights, comps, gijs,
                    comment='')
        model.rigid_elements[eid] = rbe3
        eid += 1
    #_eq_nodes_final(nid_pairs, model, tol, node_set=node_set)

    if bdf_filename_out is not None:
        model.write_bdf(bdf_filename_out, size=size, is_double=is_double)
    return model


def cut_model(model, axis='-y'):
    """
    Removes the elements on one side of a model.

    Typically aircraft are defined as x-aft, y-right, z-up, so
    all node xyz locations are positive.  We then have a xz plane
    of symmetry with the axis of symmetry being y, and typically
    save the +y elements.

    Parameters
    ----------
    model : BDF()
        the BDF object.
    axis : {'-x', '-y', '-z'}
        What direction should be removed?

    Considers:
      - nodes
      - elements
      - loads (LOAD/PLOAD4)

    ..note doesn't support +cuts (e.g. +x, +y, +z)
    """
    #model.read_bdf(bdf_filename, xref=False)
    #model.cross_reference(xref_loads=xref_loads)

    if axis == '-x':
        iaxis = 0
    elif axis == '-y':
        iaxis = 1
    elif axis == '-z':
        iaxis = 2
    else:
        raise NotImplementedError(axis)

    remove_nids = []
    for nid, node in model.nodes.items():
        xyz = node.get_position()
        if xyz[iaxis] < 0.0:
            remove_nids.append(nid)

    remove_eids = []
    for eid, element in model.elements.items():
        centroid = element.Centroid()
        if centroid[iaxis] <= 0.0:
            remove_eids.append(eid)

    #print('remove_nids =', remove_nids)
    #print('remove_eids =', remove_eids)
    for nid in remove_nids:
        del model.nodes[nid]
    for eid in remove_eids:
        del model.elements[eid]

    loads2 = {}
    for load_id, loadcase in model.loads.items():
        loadcase2 = []
        for load in loadcase:
            if load.type == 'LOAD':
                loadcase2.append(load)
            elif load.type == 'PLOAD4':
                if load.eid not in remove_eids:
                    loadcase2.append(load)
            else:
                loadcase2.append(load)
        loads2[load_id] = loadcase2
    model.loads = loads2


def _write_nodes(self, outfile, size, is_double):
    """
    Writes the NODE-type cards
    """
    if self.spoints:
        msg = []
        msg.append('$SPOINTS\n')
        msg.append(self.spoints.write_card(size, is_double))
        outfile.write(''.join(msg))

    if self.nodes:
        msg = []
        msg.append('$NODES\n')
        if self.grdset:
            msg.append(self.grdset.print_card(size))
        for (nid, node) in sorted(self.nodes.items()):
            if nid not in self.remove_nodes:
                msg.append(node.write_card(size, is_double))
        outfile.write(''.join(msg))


def get_free_edges(bdf_filename, eids=None, maps=None):
    """
    assumes no solids/bars

    Parameters
    ----------
    bdf_filename : str, BDF()
        a BDF object or filename
    maps : List[...] (default=None -> calculate)
        the output from _get_maps(eids, map_names=None,
                                  consider_0d=False, consider_0d_rigid=False,
                                  consider_1d=False, consider_2d=True, consider_3d=False)

    Returns
    -------
    free_edges : List[(int nid1, int nid2), ...]
        the free edges
    """
    if isinstance(bdf_filename, string_types):
        model = BDF(debug=False)
        model.read_bdf(bdf_filename)
    else:
        model = bdf_filename

    free_edges = []
    if maps is None:
        maps = model._get_maps(eids, map_names='edge_to_eid_map',
                               consider_0d=False, consider_0d_rigid=False,
                               consider_1d=False, consider_2d=True, consider_3d=False)
    edge_to_eid_map = maps['edge_to_eid_map']
    for edge, eids in edge_to_eid_map.items():
        if len(eids) == 2:
            continue
        free_edges.append(edge)
    return free_edges


def get_joints(model, pid_sets):
    """
    Gets the nodes at the joints of multiple property id sets

    Parameters
    ----------
    model : BDF()
        a BDF object
    pid_sets : List[pid_set, pid_set]
        set of properties IDs to boolean
        pid_set : List[int, int, int]
            set of property ids to boolean

    Returns
    -------
    joints : (N, ) int nodes
        the list of joint node ids

    Examples
    --------
    For a set of ribs, spars, and skins:
     - ribs intersect spars
     - ribs and spars intersect skins

    We want the nodes on the skin at the intersections.
    >>> pid_sets = [
        [rib1_pid, rib2_pid, rib3_pid, ...],
        [spar1_pid, spar2_pid, spar3_pid, ...],
        [skin],
    ]
    """
    nid_sets = defaultdict(set)
    for eid, elem in model.elements.items():
        pid = elem.Pid()
        nid_sets[pid].update(elem.node_ids)

    nid_array_sets = []
    for pid_set in pid_sets:
        inside = np.hstack([list(nid_sets[pid]) for pid in pid_set])
        nid_array_set = np.unique(
            inside
        )
        nid_array_sets.append(nid_array_set)
        #print('nid_array_set =', nid_array_set)
    joint_nids = reduce(np.intersect1d, nid_array_sets)
    #joint_nids = intersect1d_multi(nid_array_sets, assume_unique=True)
    return joint_nids


def extract_surface_patches(bdf_filename, starting_eids, theta_tols=40.):
    """
    Extracts the unique patches of a model based on a list of starting
    element ids and surface curvature.

    Parameters
    ----------
    bdf_filename : str
        the bdf_filename
    starting_eids : List[int]
        a list of starting element ids
    theta_tols : List[float]
        a list of tolerances for each element id
        (e.g. the nose has a different tolerance than the base)

    Returns
    -------
    model : BDF()
        the BDF object
    groups : List[Set[int]]
        the list of element ids in each group

    .. warning:: only supports CTRIA3 & CQUAD4
    """
    if isinstance(theta_tols, (float, float32)):
        theta_tols = [theta_tols] * len(starting_eids)

    model = BDF(debug=False)
    model.read_bdf(bdf_filename)

    # get data for all shell elemenst
    card_types = ['CTRIA3', 'CQUAD4']
    out = model.get_card_ids_by_card_types(card_types,
                                           reset_type_to_slot_map=False,
                                           stop_on_missing_card=False)
    shell_eids = hstack([out[card_type] for card_type in card_types])
    shell_eids.sort()

    # set_shell_eids = set(shell_eids)
    neids = len(shell_eids)
    assert neids > 0, neids
    normals = np.zeros((neids, 3), dtype='float32')
    for i, eid in enumerate(shell_eids):
        element = model.elements[eid]
        normal = element.Normal()
        normals[i, :] = normal

    #print('shell_eids = %s' % shell_eids)

    # get neighboring shell elements
    out_map = model._get_maps(eids=shell_eids, map_names='edge_to_eid_map')
    edge_to_eid_map = out_map['edge_to_eid_map']

    eid_to_eid_map = defaultdict(set)
    #if 1:
    for edge, eids in edge_to_eid_map.items():
        for eid_a in eids:
            for eid_b in eids:
                eid_to_eid_map[eid_a].add(eid_b)
    # else:
        # for edge, eids in edge_to_eid_map.items():
            # for eid_a in eids:
                # for eid_b in eids:
                    # if eid_a < eid_b:
                        # eid_to_eid_map[eid_a].add(eid_b)
                        # eid_to_eid_map[eid_b].add(eid_a)

    #print('\neid_to_eid_map:')
    #for eid, eids in eid_to_eid_map.items():
        #print('%-6s %s' % (eid, eids))

    groups = []
    # now trace the model
    for starting_eid, theta_tol in zip(starting_eids, theta_tols):
        print('starting_eid = %s' % starting_eid)
        group = set([])
        check = set([starting_eid])
        while check:
            eid = next(iter(check))
            #print('  eid = %s' % eid)
            neighboring_eids = eid_to_eid_map[eid]
            #print('    neighbors = %s' % neighboring_eids)

            # don't double check eids
            neigboring_eids_to_check = array(list(neighboring_eids - group), dtype='int32')
            nneighbors = len(neigboring_eids_to_check)
            if nneighbors == 0:
                check.remove(eid)
                continue
            assert nneighbors > 0, neigboring_eids_to_check

            i = searchsorted(shell_eids, eid)
            #assert len(i) > 0, 'eid=%s cant be found' % eid
            local_normal = normals[i, :]

            # find the normals
            i = searchsorted(shell_eids, neigboring_eids_to_check)
            assert len(i) > 0, 'eid=%s cant be found' % eid
            normal = normals[i, :]

            # a * b = |a| |b| cos(theta)
            # |a| = |b| = 1
            # cos^-1(a*b) = theta

            # we flip the dimensions because matrix shapes are stupid
            #theta = degrees(arccos(dot(local_normal, normal)))
            theta = degrees(arccos(dot(local_normal, normal.T).T))

            #print('    theta[eid=%s] = %s; n=%s' % (eid, theta, nneighbors))
            assert len(theta) == nneighbors, len(theta)

            itol = where(theta <= theta_tol)[0]
            eids_within_tol = neigboring_eids_to_check[itol]
            group.update(eids_within_tol)
            check.update(eids_within_tol)
            check.remove(eid)
        groups.append(group)
    return model, groups


def split_model_by_material_id(bdf_filename, bdf_filename_base,
                               encoding=None, size=8, is_double=False):
    """
    Splits a model based on the material ID

    Parameters
    ----------
    bdf_filename : str
        the filename to read (e.g. fem.bdf)
    bdf_filename_base : str
        the prefix to the output file (e.g. fem)
    encoding : str
        the unicode encoding
    size : int; default=8
        the field width to write
    is__double : bool; default=False
        should double precision be used

    .. warning:: only considers elements with materials (so no CBUSH, but yes to CONROD)
    .. warning:: PSHELL only considers mid1 (not mid2, mid3, mid4)
    .. warning:: doesn't consider PCOMPs
    .. warning:: doesn't consider SPCs/loads/etc.

    .. warning:: hasn't really been tested
    """
    model = BDF(debug=False)
    model.read_bdf(bdf_filename, xref=True, encoding=encoding)

    mid_to_eids_map = defaultdict(set)
    #nids_to_write = set([])
    elements_with_properties = [
        'CQUAD4', 'CQUAD', 'CQUAD8',
        'TRIA3', 'CTRIA6', 'CTRIAX', 'CTRIAX6',
        'CTETRA', 'CPENTA', 'CHEXA', 'CPYRAM',
        'CROD', 'CONROD', 'CBEAM', 'CBEND',
    ]

    elements_without_properties_with_materials = [
        'CONROD',
    ]
    #invalid_properties = ['PCOMP', 'PCOMPG']
    for eid, elem in model.elements.items():
        etype = elem.type
        if etype in elements_with_properties:
            pid = elem.pid
            ptype = pid.type
            if ptype == 'PSHELL':
                mid1 = pid.Mid1()
                mid_to_eids_map[mid1].add(eid)
            else:
                msg = 'skipping eid=%s elem.type=%s pid=%s prop.type=%s' % (
                    eid, etype, pid.pid, ptype)
                model.log.warning(msg)
        elif etype in elements_without_properties_with_materials:
            mid = elem.Mid()
            mid_to_eids_map[mid].add(eid)
        else:
            model.log.warning('skipping eid=%s elem.type=%s' % (eid, etype))

    for mid, eids in mid_to_eids_map.items():
        if not eids:
            continue
        bdf_filename_out = bdf_filename_base + '_mid=%s.bdf' % mid

        with open(bdf_filename_out, 'w') as bdf_file:
            bdf_file.write('$ pyNastran : punch=True\n')
            mat = model.materials[mid]
            bdf_file.write(mat.write_card(size=size, is_double=is_double))

            nids_to_write = set([])
            for eid in eids:
                elem = model.elements[eid]
                nids_to_write.update(elem.node_ids)
                bdf_file.write(elem.write_card(size=size, is_double=is_double))

            for nid in nids_to_write:
                node = model.nodes[nid]
                bdf_file.write(node.write_card(size=size, is_double=is_double))
            bdf_file.write('ENDDATA\n')


def create_spar_cap(model, eids, nids, width, nelements=1, symmetric=True, xyz_cid0=None,
                    vector1=None, vector2=None, idir=None, eid_start=1):
    """
    Builds elements along a line of nodes that are normal to the element face.

    .. code-block:: console

        1    2     3      4     5  D  4  E  6
        +----+-----+------+     *-----+-----*
        | A  |   B |   C  |           |
        +----+-----+------+           |
             (side view)       (along the axis)

    Nodes at the upper corner of Element A, B, C.

    Parameters
    ----------
    xyz_cid0 : (n, 3) ndarray
        nodes in cid=0
    width : float
        the width of the spar cap
    nelements : int; default=1
        the number of nodes to create that are normal to the plane
    vector1 : (float, float, float)
        overwrite the normal of the normals (not supported)
    vector2 : (float, float, float)
        overwrite the normal of the normals (not supported)

    idir : int; default=1
        A primary direction to which the normals should be consistent (not supported).
        0 : x
        1 : y
        2 : z
        For an aircraft rib, y would be the preferred direction for the normal.
        For an aircraft spar, x/z would work.

    Examples
    --------
    >>> eids = [A, B, C]
    >>> nids = [1, 2, 3, 4]
    >>> width = 3.0
    >>> nodes, elements = create_spar_cap(model, eids, nids, width)
    >>> nodes
    [5, 6]
    >>> elements
    [D, E]
    """
    assert vector1 is None, vector1
    assert vector2 is None, vector2
    assert idir is None, idir
    assert nelements == 1, idir

    if xyz_cid0 is None:
        xyz_cid0 = model.get_xyz_in_coord(cid=0)
    #print(xyz_cid0)
    nid_cd = np.array([[nid, node.Cd()] for nid, node in sorted(model.nodes.items())])
    all_nids = nid_cd[:, 0]

    width_array = np.linspace(0., width, num=nelements, endpoint=True)[1:]

    nodes = []
    elements = []
    nid_start = all_nids.max() + 1
    eids = np.asarray(eids, dtype='int32')
    nids = np.asarray(nids, dtype='int32')

    # Create the CQUAD4s
    map_nid = {}
    map_nid2 = {}
    neids = eids.size
    all_common_nids = set([])
    normals = np.zeros((neids, 3), dtype='float64')
    for eid in eids:
        elem = model.elements[eid]
        pid = elem.Pid()
        enids = elem.node_ids
        common_nids = np.where(np.in1d(enids, nids)) # A in B
        all_common_nids.update(common_nids)
        i = searchsorted(all_nids, enids)
        xyz = xyz_cid0[i, :]
        etype = elem.type
        if etype == 'CQUAD4':
            # 3---2
            # |   |
            # 0---1
            normal = np.cross(
                xyz[2, :] - xyz[0, :],
                xyz[3, :] - xyz[1, :],
            )
        elif etype == 'CTRIA3':
            normal = np.cross(
                xyz[1, :] - xyz[0, :],
                xyz[2, :] - xyz[0, :],
            )
        else:
            raise NotImplementedError(etype)
        normi = np.linalg.norm(normal)
        assert np.allclose(normi, 1.0), normi
        normal /= normi

        # add some new nodes
        nid0, nid1 = common_nids
        if nid0 not in map_nid:
            map_nid[nid0] = np.arange(nid_start, nid_start + nelements)
            nid_start += nelements
        if nid1 not in map_nid:
            map_nid[nid1] = np.arange(nid_start, nid_start + nelements)
            nid_start += nelements

        nid3 = map_nid[nid1][0]
        nid4 = map_nid[nid0][0]
        nid0i = nid0
        nid1i = nid1
        for ielement in range(nelements):
            nid3 = map_nid[nid1][ielement]
            nid4 = map_nid[nid0][ielement]
            new_nids = [nid0i, nid1i, nid3, nid4]
            cquad4 = CQUAD4(eid_start, pid, new_nids)
            elements.append(cquad4)
            eid_start += 1
            nid0i = nid3
            nid1i = nid4

        if symmetric:
            if nid0 not in map_nid2:
                map_nid2[nid0] = np.arange(nid_start, nid_start + nelements)
                nid_start += nelements
            if nid1 not in map_nid2:
                map_nid2[nid1] = np.arange(nid_start, nid_start + nelements)
                nid_start += nelements

            nid3 = map_nid2[nid1][0]
            nid4 = map_nid2[nid0][0]
            nid0i = nid0
            nid1i = nid1
            for ielement in range(nelements):
                nid3 = map_nid2[nid1][ielement]
                nid4 = map_nid2[nid0][ielement]
                new_nids = [nid0i, nid1i, nid3, nid4]
                cquad4 = CQUAD4(eid_start, pid, new_nids)
                elements.append(cquad4)
                eid_start += 1
                nid0i = nid3
                nid1i = nid4

    # Create the GRIDs
    inids = searchsorted(all_nids, nids)
    xyz = xyz_cid0[inids, :]
    for nid in all_common_nids:
        mapped_nids = map_nid[nid]
        avg_normal_at_node = np.zeros(3, dtype='float64')
        node = model.nodes[nid]
        node_elems = nid.elements
        nelems = len(node_elems)
        for elem in node_elems:
            eid = elem.eid
            j = searchsorted(eids, eid)
            avg_normal_at_node += normals[j, :]
        avg_normal_at_node /= nelems

        for mapped_nid in mapped_nids:
            xyzi = xyz[i, :] + avg_normal_at_node * width_array[i]
            node1 = GRID(mapped_nid, xyz=xyzi)
            nodes.append(node1)

        if symmetric:
            mapped_nids2 = map_nid2[nid]
            for imap, mapped_nid2 in enumerate(mapped_nids2):
                xyzi = xyz[i, :] + avg_normal_at_node * width_array[i]
                node2 = GRID(mapped_nid2, xyz=xyzi)
                nodes.append(node2)
    return nodes, elements
