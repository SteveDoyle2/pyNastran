from itertools import count
#from copy import deepcopy
from typing import Tuple, Dict, Union

import numpy as np
from pyNastran.utils.numpy_utils import zip_strict
from pyNastran.bdf.bdf import BDF, CTRIA3, CQUAD4, GRID, CBAR, CHEXA8 #CTRIA6, CQUAD8,
from pyNastran.bdf.mesh_utils.internal_utils import get_bdf_model, BDF_FILETYPE


def refine_model(bdf_filename: BDF_FILETYPE, refinement_ratio: int=2) -> BDF:
    """xref should be turned off"""
    model = get_bdf_model(bdf_filename, xref=False, cards_to_skip=None,
                          validate=True, log=None, debug=False)
    model.cross_reference(
        xref=True, xref_nodes=True, xref_elements=True,
        xref_nodes_with_elements=False, xref_properties=False, xref_masses=False,
        xref_materials=False, xref_loads=False, xref_constraints=False,
        xref_aero=False, xref_sets=False, xref_optimization=False, word='')

    out = model.get_displacement_index_xyz_cp_cd(fdtype='float64', idtype='int32', sort_ids=True)
    icd_transform, icp_transform, xyz_cp, nid_cp_cd = out

    xyz_cid0 = xyz_cp
    all_nodes = nid_cp_cd[:, 0]

    #nodes_new = []
    #xyz_new = []
    #elems = []
    nid0 = max(model.nodes) + 1
    eid0 = max(model.elements) + 1
    nodes = model.nodes
    nelements = 0
    elements = {}

    debug = False
    nnodes_to_add = 1
    nnodes_to_add_with_ends = nnodes_to_add + 2
    nid0, edges_to_center, faces_to_center = _get_edge_to_center_nodes(
        model,
        all_nodes, xyz_cid0,
        nid0, nnodes_to_add)

    if debug:
        for edge, cen in edges_to_center.items():
            print(f'e={edge} cen={cen}')

    for eid, elem in model.elements.items():
        if elem.type == 'CTRIA3':
            pass
            nid0, eid0, nelements = _refine_tri(
                model, all_nodes, xyz_cid0,
                nodes, elements,
                edges_to_center, faces_to_center,
                eid, elem,
                nid0, eid0,
                nelements, nnodes_to_add_with_ends)
        elif elem.type == 'CQUAD4':
            pass
            nid0, eid0, nelements = _refine_quad(
                model, all_nodes, xyz_cid0,
                nodes, elements,
                edges_to_center, faces_to_center,
                eid, elem,
                nid0, eid0,
                nelements, nnodes_to_add_with_ends)
        elif elem.type == 'CHEXA':
            nid0, eid0, nelements = _refine_hexa(
                model, all_nodes, xyz_cid0,
                nodes, elements,
                edges_to_center, faces_to_center,
                eid, elem,
                nid0, eid0,
                nelements, nnodes_to_add_with_ends)
        elif elem.type == 'CBAR':
            n1, n2 = elem.nodes
            n3 = nid0
            (in1, in2) = np.searchsorted(all_nodes, elem.nodes)
            xyz1 = xyz_cid0[in1, :]
            xyz2 = xyz_cid0[in2, :]
            centroid = (xyz1 + xyz2) / 2.
            nodes.update({
                n3: GRID(n3, centroid, cp=0, cd=0, ps='', seid=0, comment=''),
            })
            elem1 = CBAR(eid, elem.pid, [n1, n3], elem.x, elem.g0,
                         pa=elem.pa, wa=elem.wa, wb=None, offt=elem.offt, comment=elem.comment)
            elem2 = CBAR(eid0, elem.pid, [n3, n2], elem.x, elem.g0,
                         pb=elem.pb, wa=None, wb=elem.wb, offt=elem.offt)
            elements.update({
                elem1.eid: elem1,
                elem2.eid: elem2})
            nid0 += 1
            eid0 += 1
        #elif elem.type in {'CROD', 'CONROD', 'CTUBE', 'CBAR', 'CBEAM'}:
            #continue
        #elif elem.type == 'CHEXA':
            #continue
        else:
            print(elem)
    #model.nodes = nodes
    model.elements = elements
    model.loads = {}
    model.load_combinations = {}
    return model

def _refine_tri(model: BDF,
                all_nodes, xyz_cid0,
                nodes, elements,
                edges_to_center, faces_to_center,
                eid, elem,
                nid0, eid0,
                nelements, nnodes_to_add_with_ends, debug=False):
    #continue
    n1i, n2i, n3i = elem.nodes
    #n4, n5, n6 = nid0, nid0 + 1, nid0 + 2
    (in1, in2, in3) = np.searchsorted(all_nodes, elem.nodes)
    #nodes = [n1, n2, n3, n4, n5, n6]

    xyz1 = xyz_cid0[in1, :]
    xyz2 = xyz_cid0[in2, :]
    xyz3 = xyz_cid0[in3, :]

    edges = elem.get_edge_ids()
    forward_edges = [(n1i, n2i), (n2i, n3i), (n3i, n1i)]

    nids_array = np.zeros((nnodes_to_add_with_ends, nnodes_to_add_with_ends), dtype='int32')
    nids_array[0, 0] = n1i
    nids_array[0, -1] = n2i
    nids_array[-1, -1] = n3i
    assert nids_array[0, 0] == n1i
    assert nids_array[0, -1] == n2i
    assert nids_array[-1, -1] == n3i
    #print('nids_array0')
    #print(nids_array)
    nid0 = _insert_tri_nodes(
        nodes, nids_array, nid0,
        edges, forward_edges, edges_to_center,
        nnodes_to_add_with_ends,
        xyz1, xyz2, xyz3, debug=debug,
    )
    assert nids_array[0, 0] == n1i
    assert nids_array[0, -1] == n2i
    assert nids_array[-1, -1] == n3i

    n1i, n4i, n2i = nids_array[0, :]
    n6i, n5i = nids_array[1, 1:]
    n3i = nids_array[2, 2]

    args = {'theta_mcid': elem.theta_mcid,
            'zoffset': elem.zoffset,
            'tflag': elem.tflag,
            'T1': None, 'T2': None, 'T3': None, }

    elem1 = CTRIA3(eid,    elem.pid, [n1i, n4i, n6i], comment=elem.comment, **args)
    elem2 = CTRIA3(eid0,   elem.pid, [n4i, n2i, n5i], **args)
    elem3 = CTRIA3(eid0+1, elem.pid, [n6i, n5i, n3i], **args)
    elem4 = CTRIA3(eid0+2, elem.pid, [n4i, n5i, n6i], **args)
    elements.update({
        elem1.eid: elem1,
        elem2.eid: elem2,
        elem3.eid: elem3,
        elem4.eid: elem4})
    nelements += 4

    #centroid = (xyz1 + xyz2 + xyz3) / 3.
    #xyz4 = (xyz1 + xyz2) / 2.
    #xyz5 = (xyz2 + xyz3) / 2.
    #xyz6 = (xyz3 + xyz1) / 2.
    #nodes.update({
        #n4: GRID(n4, xyz4, cp=0, cd=0, ps='', seid=0, comment=''),
        #n5: GRID(n5, xyz5, cp=0, cd=0, ps='', seid=0, comment=''),
        #n6: GRID(n6, xyz6, cp=0, cd=0, ps='', seid=0, comment=''),
    #})
    elem1.validate()
    elem2.validate()
    elem3.validate()
    elem4.validate()
    elem1.cross_reference(model)
    elem2.cross_reference(model)
    elem3.cross_reference(model)
    elem4.cross_reference(model)
    elem1.Normal()
    elem2.Normal()
    elem3.Normal()
    elem4.Normal()
    #xyz_new.extend([xyz4, xyz5, xyz6])
    #nodes_new.extend([n4, n5, n6])
    #nid0 += 3
    eid0 += 3
    del n1i, n2i, n3i, n4i, n5i, n6i, xyz1, xyz2, xyz3 #, centroid
    assert len(elements) == nelements
    return nid0, eid0, nelements

def _refine_quad(model: BDF,
                all_nodes, xyz_cid0,
                nodes, elements,
                edges_to_center, faces_to_center,
                eid, elem,
                nid0, eid0,
                nelements, nnodes_to_add_with_ends):
    n1i, n2i, n3i, n4i = elem.nodes
    #n5, n6, n7, n8, n9 = nid0, nid0 + 1, nid0 + 2, nid0 + 3, nid0 + 4
    (in1, in2, in3, in4) = np.searchsorted(all_nodes, elem.nodes)
    xyz1 = xyz_cid0[in1, :]
    xyz2 = xyz_cid0[in2, :]
    xyz3 = xyz_cid0[in3, :]
    xyz4 = xyz_cid0[in4, :]
    #centroid = (xyz1 + xyz2 + xyz3 + xyz4) / 4.

    edges = elem.get_edge_ids()
    forward_edges = [(n1i, n2i), (n2i, n3i), (n3i, n4i), (n4i, n1i)]

    nids_array = np.zeros((nnodes_to_add_with_ends, nnodes_to_add_with_ends), dtype='int32')
    nids_array[0, 0] = n1i
    nids_array[0, -1] = n2i
    nids_array[-1, -1] = n3i
    nids_array[-1, 0] = n4i
    assert nids_array[0, 0] == n1i
    assert nids_array[0, -1] == n2i
    assert nids_array[-1, -1] == n3i
    assert nids_array[-1, 0] == n4i

    #if eid == 28286:
        #debug = True
        #x = 1

    debug = False
    if debug:
        print('nids_array0:\n', nids_array)
    nid0 = _insert_quad_nodes(
        nodes, nids_array, nid0,
        edges, forward_edges,
        edges_to_center,
        faces_to_center,
        nnodes_to_add_with_ends,
        [n1i, n2i, n3i, n4i],
        xyz1, xyz2, xyz3, xyz4, debug=debug,
    )
    unids1 = np.unique([n1i, n2i, n3i, n4i])
    assert len(unids1) == 4, unids1
    assert nids_array[0, 0] == n1i
    assert nids_array[0, -1] == n2i
    assert nids_array[-1, -1] == n3i
    assert nids_array[-1, 0] == n4i
    #assert unids1 == unids2

    n1, n2, n3, n4 = _quad_nids_to_node_ids(nids_array)
    nelementsi = len(n1) - 1
    eids = [eid] + [eid0 + i for i in range(nelementsi)]
    args = {'theta_mcid': elem.theta_mcid,
            'zoffset': elem.zoffset,
            'tflag': elem.tflag,
            'T1': None, 'T2': None, 'T3': None, 'T4': None}

    pid = elem.pid
    for eidi, n1i, n2i, n3i, n4i in zip_strict(eids, n1, n2, n3, n4):
        nidsi = [n1i, n2i, n3i, n4i]
        if debug:
            print('quad', eidi, nidsi)
        elem1 = CQUAD4(eidi, pid, nidsi, **args)
        #print(elem1)
        elem1.validate()
        elem1.cross_reference(model)
        try:
            elem1.Normal()
        except:
            print(elem1.nodes)
            raise
        elements[eidi] = elem1
        nelements += 1
    eid0 += nelementsi

    #xyz5 = (xyz1 + xyz2) / 2.
    #xyz6 = (xyz2 + xyz3) / 2.
    #xyz7 = (xyz3 + xyz4) / 2.
    #xyz8 = (xyz4 + xyz1) / 2.
    #nodes.update({
        #n5: GRID(n5, xyz5, cp=0, cd=0, ps='', seid=0, comment=''),
        #n6: GRID(n6, xyz6, cp=0, cd=0, ps='', seid=0, comment=''),
        #n7: GRID(n7, xyz7, cp=0, cd=0, ps='', seid=0, comment=''),
        #n8: GRID(n8, xyz8, cp=0, cd=0, ps='', seid=0, comment=''),
        #n9: GRID(n9, centroid, cp=0, cd=0, ps='', seid=0, comment=''),
    #})
    #nodes_new.extend([n5, n6, n7, n8, n9])
    #xyz_new.extend([xyz5, xyz6, xyz7, xyz8, centroid])

    #elem1.cross_reference(model)
    #elem2.cross_reference(model)
    #elem3.cross_reference(model)
    #elem4.cross_reference(model)
    #elem1.Normal()
    #elem2.Normal()
    #elem3.Normal()
    #elem4.Normal()
    #nid0 += 5
    #del n1, n2, n3, n4, xyz1, xyz2, xyz3, xyz4
    #print(list(elements.keys()))
    assert len(elements) == nelements
    #break
    return nid0, eid0, nelements

def _refine_hexa(model: BDF,
                all_nodes, xyz_cid0,
                nodes, elements,
                edges_to_center, faces_to_center,
                eid, elem,
                nid0, eid0,
                nelements, nnodes_to_add_with_ends):
    n1i, n2i, n3i, n4i, n5i, n6i, n7i, n8i = elem.nodes
    (in1, in2, in3, in4, in5, in6, in7, in8) = np.searchsorted(all_nodes, elem.nodes)
    xyz1 = xyz_cid0[in1, :]
    xyz2 = xyz_cid0[in2, :]
    xyz3 = xyz_cid0[in3, :]
    xyz4 = xyz_cid0[in4, :]

    xyz5 = xyz_cid0[in5, :]
    xyz6 = xyz_cid0[in6, :]
    xyz7 = xyz_cid0[in7, :]
    xyz8 = xyz_cid0[in8, :]

    #centroid = (xyz1 + xyz2 + xyz3 + xyz4) / 4.

    edges = elem.get_edge_ids()
    forward_edges = [
        (n1i, n2i), (n2i, n3i), (n3i, n4i), (n4i, n1i),
        (n5i, n6i), (n6i, n7i), (n7i, n8i), (n8i, n5i),
        (n1i, n5i), (n2i, n6i), (n3i, n7i), (n4i, n8i),
    ]
    assert len(edges) == len(forward_edges), len(edges)
    nids_array = np.zeros(
        (nnodes_to_add_with_ends, nnodes_to_add_with_ends, nnodes_to_add_with_ends),
        dtype='int32')
    nids_array[0, 0, 0] = n1i
    nids_array[-1, 0, 0] = n2i
    nids_array[-1, -1, 0] = n3i
    nids_array[0, -1, 0] = n4i

    nids_array[0, 0, -1] = n5i
    nids_array[-1, 0, -1] = n6i
    nids_array[-1, -1, -1] = n7i
    nids_array[0, -1, -1] = n8i

    assert nids_array[0, 0, 0] == n1i
    assert nids_array[-1, 0, 0] == n2i
    assert nids_array[-1, -1, 0] == n3i
    assert nids_array[0, -1, 0] == n4i
    assert nids_array[0, 0, -1] == n5i
    assert nids_array[-1, 0, -1] == n6i
    assert nids_array[-1, -1, -1] == n7i
    assert nids_array[0, -1, -1] == n8i

    #if eid == 28286:
        #debug = True
        #x = 1

    debug = False
    if debug:
        print('nids_array0:\n', nids_array)
    faces = hexa_get_sorted_faces(elem)
    nid0 = _insert_hexa_nodes(
        nodes, nids_array, nid0,
        edges, forward_edges,
        edges_to_center, faces_to_center, faces,
        nnodes_to_add_with_ends,
        xyz1, xyz2, xyz3, xyz4,
        xyz5, xyz6, xyz7, xyz8,
        debug=debug,
    )
    unids1 = np.unique([n1i, n2i, n3i, n4i, n5i, n6i, n7i, n8i])
    assert len(unids1) == 8, unids1
    assert nids_array[0, 0, 0] == n1i
    assert nids_array[-1, 0, 0] == n2i
    assert nids_array[-1, -1, 0] == n3i
    assert nids_array[0, -1, 0] == n4i
    assert nids_array[0, 0, -1] == n5i
    assert nids_array[-1, 0, -1] == n6i
    assert nids_array[-1, -1, -1] == n7i
    assert nids_array[0, -1, -1] == n8i
    #assert unids1 == unids2

    nodes_hexas = _hexa_nids_to_node_ids(nids_array)

    nelementsi = nodes_hexas.shape[0] - 1
    eids = [eid] + [eid0 + i for i in range(nelementsi)]

    debug = True
    pid = elem.pid
    for eidi, nidsi in zip_strict(eids, nodes_hexas):
        if debug:
            print('hexa', eidi, nidsi)
        #elem1 = CHEXA(eidi, pid, nidsi, **args)
        elem1 = CHEXA8(eidi, pid, nidsi, comment='')
        #print(elem1)
        elem1.validate()
        elem1.cross_reference(model)
        try:
            elem1.Volume()
        except:
            print(elem1.nodes)
            raise
        elements[eidi] = elem1
        nelements += 1
    eid0 += nelementsi

    #xyz5 = (xyz1 + xyz2) / 2.
    #xyz6 = (xyz2 + xyz3) / 2.
    #xyz7 = (xyz3 + xyz4) / 2.
    #xyz8 = (xyz4 + xyz1) / 2.
    #nodes.update({
        #n5: GRID(n5, xyz5, cp=0, cd=0, ps='', seid=0, comment=''),
        #n6: GRID(n6, xyz6, cp=0, cd=0, ps='', seid=0, comment=''),
        #n7: GRID(n7, xyz7, cp=0, cd=0, ps='', seid=0, comment=''),
        #n8: GRID(n8, xyz8, cp=0, cd=0, ps='', seid=0, comment=''),
        #n9: GRID(n9, centroid, cp=0, cd=0, ps='', seid=0, comment=''),
    #})
    #nodes_new.extend([n5, n6, n7, n8, n9])
    #xyz_new.extend([xyz5, xyz6, xyz7, xyz8, centroid])

    #elem1.cross_reference(model)
    #elem2.cross_reference(model)
    #elem3.cross_reference(model)
    #elem4.cross_reference(model)
    #elem1.Normal()
    #elem2.Normal()
    #elem3.Normal()
    #elem4.Normal()
    #nid0 += 5
    #del n1, n2, n3, n4, xyz1, xyz2, xyz3, xyz4
    #print(list(elements.keys()))
    assert len(elements) == nelements
    #break
    return nid0, eid0, nelements

def _extract_edges(edges_to_center,
                   nodes: Dict[int, GRID],
                   all_nodes, xyz_cid0,
                   nnodes_to_add: int,
                   nnodes_to_add_with_ends: int,
                   nid0: int,
                   elem: Union[CTRIA3, CQUAD4, CBAR, CHEXA8],
                   debug: bool=False) -> int:
    edges = elem.get_edge_ids()
    for edge in edges:
        if edge in edges_to_center:
            continue

        n1, n2 = edge
        nodes_to_add = [nid0 + i for i in range(nnodes_to_add)]
        in1, in2 = np.searchsorted(all_nodes, edge)
        xyz1 = xyz_cid0[in1, :]
        xyz2 = xyz_cid0[in2, :]
        if debug:
            print('edge =', edge)
        edges_to_center[edge] = (n1, *nodes_to_add, n2)
        xi = np.linspace(0., 1., num=nnodes_to_add_with_ends, endpoint=True)[1:-1]
        for xii in xi:
            xyz = xyz1 * (1. - xii) + xyz2 * xii
            nodes[nid0] = GRID(nid0, xyz)
            nid0 += 1
        x = 1
    return nid0

def _get_edge_to_center_nodes(model: BDF,
                              all_nodes: np.ndarray, xyz_cid0: np.ndarray,
                              nid0: int, nnodes_to_add: int,
                              debug=False):
    nnodes_to_add_with_ends = nnodes_to_add + 2
    nodes = model.nodes
    edges_to_center = {}
    for elem in model.elements.values():
        if elem.type in {'CTRIA3', 'CQUAD4', 'CBAR', 'CHEXA'}:
            nid0 = _extract_edges(
                edges_to_center,
                nodes, all_nodes, xyz_cid0,
                nnodes_to_add, nnodes_to_add_with_ends,
                nid0, elem)
        else:
            print(elem)

    elements_1d = {
        'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
        'CONROD', 'CROD', 'CTUBE',
        'CBAR', 'CBEAM', 'CBUSH', 'CBUSH1D', 'CBUSH2D',
        '',
    }
    nodes = model.nodes
    faces_to_center = {}
    for elem in model.elements.values():
        if elem.type in {'CQUAD4'}:
            xi = np.linspace(0., 1., num=nnodes_to_add_with_ends, endpoint=True)[1:-1]
            xj = np.linspace(0., 1., num=nnodes_to_add_with_ends, endpoint=True)[1:-1]
            face0 = elem.nodes
            face = _sort_face(face0)
            if face in faces_to_center:
                continue

            nodes_list = []
            (in1, in2, in3, in4) = np.searchsorted(all_nodes, elem.nodes)
            xyz1 = xyz_cid0[in1, :]
            xyz2 = xyz_cid0[in2, :]
            xyz3 = xyz_cid0[in3, :]
            xyz4 = xyz_cid0[in4, :]
            for xii, xjj in zip(xi, xj):
                xyz12 = xyz1 * (1. - xii) + xyz2 * xii
                xyz34 = xyz3 * (1. - xii) + xyz4 * xii # TODO: probably backwards
                xyzc = xyz12 * (1. - xjj) + xyz34 * xjj
                nodes[nid0] = GRID(nid0, xyzc)
                nodes_list.append(nid0)
                nid0 += 1
            faces_to_center[face] = nodes_list
            print(f'*adding face={face} nodes={nodes_list}')
            x = 1

        elif elem.type in {'CHEXA'}:
            #faces = elem.get_sorted_faces()
            xarray = np.linspace(0., 1., num=nnodes_to_add_with_ends, endpoint=True)
            faces = hexa_get_sorted_faces(elem)

            #(in1, in2, in3, in4,
             #in5, in6, in7, in8) = np.searchsorted(all_nodes, elem.nodes)
            #xyz1 = xyz_cid0[in1, :]
            #xyz2 = xyz_cid0[in2, :]
            #xyz3 = xyz_cid0[in3, :]
            #xyz4 = xyz_cid0[in4, :]

            #xyz5 = xyz_cid0[in5, :]
            #xyz6 = xyz_cid0[in6, :]
            #xyz7 = xyz_cid0[in7, :]
            #xyz8 = xyz_cid0[in8, :]

            for face in faces:
                if face in faces_to_center:
                    continue
                nodes_list = []
                (in1, in2, in3, in4) = np.searchsorted(all_nodes, face)
                xyz1 = xyz_cid0[in1, :]
                xyz2 = xyz_cid0[in2, :]
                xyz3 = xyz_cid0[in3, :]
                xyz4 = xyz_cid0[in4, :]

                for xi in xarray:
                    xyz12 = xyz1 * (1. - xi) + xyz2 * xi
                    xyz43 = xyz4 * (1. - xi) + xyz3 * xi
                    for xj in xarray:
                        if xi in {0., 1.} or xj in {0., 1.}:
                            # or because it's 2d
                            # and edges have been handled
                            continue
                        xyz = xyz12 * (1. - xj) + xyz43 * xj

                        print(xi, xj, xyz)
                        nodes[nid0] = GRID(nid0, xyz)
                        nodes_list.append(nid0)
                        nid0 += 1
                        y = 1
                print(f'face={face} nodes_list={nodes_list}')
                faces_to_center[face] = nodes_list
                # end of face
            # end of faces

            x = 1
        elif elem.type in elements_1d:
            continue
        else:
            print(elem)
    return nid0, edges_to_center, faces_to_center

def hexa_get_sorted_faces(elem: CHEXA8):
    (n1, n2, n3, n4, n5, n6, n7, n8) = elem.nodes
    faces = [
        # top/btm
        [n1, n2, n3, n4],
        [n5, n6, n7, n8],

        #fore/aft
        [n1, n2, n6, n5],
        [n4, n3, n7, n8],

        # left/right
        [n1, n4, n8, n5],
        [n2, n3, n7, n6],
    ]
    sorted_faces = []
    for face in faces:
        face2 = _sort_face(face)
        sorted_faces.append(face2)
    return sorted_faces

def _sort_face(face):
    """
    A sorted face starts at the minimum id
    and steps up to the next lowest. Then it just
    continues and order doesn't matter
    """
    iface = face.index(min(face))
    face2 = face[iface:] + face[:iface]
    assert len(face) == len(face2), face2
    # flip face because:
    #  [n1, n4, n3, n2]
    # is not in simplest form, so we change it to:
    #  [n1, n2, n3, n4]
    #
    # by reversing it:
    # [n2, n3, n4, n1]
    # and:
    # slicing it
    #print('-----')
    #print(f'face = {face}')
    #print(f'face2 = {face2}')

    # reverse if n2 < n4
    if face2[1] > face2[-1]:
        face2 = [face2[0]] + face2[1:][::-1]
        #print(f'*face3 = {face3}')
        #x = 1
    assert len(face2) == len(face)
    return tuple(face2)

def _insert_tri_nodes(nodes: Dict[int, GRID],
                       nids_array, nid0: int,
                       edges, forward_edges, edges_to_center, nnodes_to_add_with_ends: int,
                       xyz1, xyz2, xyz3, debug=False):
    for i, edge, fwd_edge in zip(count(), edges, forward_edges):
        nids_center = edges_to_center[edge]
        if debug:
            print('nids_center =', nids_center)
        flag = ' '
        if edge == fwd_edge:
            nids_set = nids_center
        else:
            flag = '*'
            nids_set = list(nids_center)
            nids_set.reverse()
            if debug:
                print('nids_set =', nids_set)
        if i == 0:
            nids_array[0, :] = nids_set
        elif i == 1:
            nids_array[:, -1] = nids_set
        else:
            # this diagonal is added in reverse
            nids_set2 = nids_set[::-1]
            for i in range(nnodes_to_add_with_ends):
                nids_array[i, i] = nids_set2[i]
            #nids_array[, ::-1] = nids_set
            y = 1
        if debug:
            print(f'{flag}i={i} nids={nids_set} edge={edge} fwd_edge={fwd_edge}')
            print(nids_array)
            print('---------')


    if nnodes_to_add_with_ends == 3:
        return nid0
    raise NotImplementedError(nnodes_to_add_with_ends)
    if nids_array.min() == 0:
        i, j = np.where(nids_array == 0)
        xi = np.linspace(0., 1., num=nnodes_to_add_with_ends, endpoint=True)#[1:-1]
        xj = np.linspace(0., 1., num=nnodes_to_add_with_ends, endpoint=True)#[1:-1]

        # bilinear interpolation
        for ii, jj in zip(i, j):
            xii = xi[ii]
            xjj = xj[jj]
            xyz12 = xyz1 * (1. - xii) + xyz2 * xii
            xyz34 = xyz3 * (1. - xii)
            xyzc = xyz12 * (1. - xjj) + xyz34 * xjj
            nodes[nid0] = GRID(nid0, xyzc)
            nids_array[i, j] = nid0
            nid0 += 1
    x = 1

    #print('nids_array1:\n', nids_array)
    new_corner_nids = [
        nids_array[0, 0],
        nids_array[0, -1],
        nids_array[-1, -1],
        nids_array[-1, 0],
    ]
    unids2 = np.unique(new_corner_nids)
    assert len(unids2) == 3, unids2

    return nid0

def _insert_quad_nodes(nodes: Dict[int, GRID],
                       nids_array, nid0: int,
                       edges, forward_edges,
                       edges_to_center,
                       faces_to_center,
                       nnodes_to_add_with_ends: int,
                       face,
                       xyz1, xyz2, xyz3, xyz4,
                       debug=False):
    for i, edge, fwd_edge in zip(count(), edges, forward_edges):
        nids_center = edges_to_center[edge]
        if debug:
            print('nids_center =', nids_center)
        flag = ' '
        if edge == fwd_edge:
            nids_set = nids_center
        else:
            #flag = '*'
            nids_set = list(nids_center)
            nids_set.reverse()
            if debug:
                print('nids_set =', nids_set)
        if i == 0:
            nids_array[0, :] = nids_set
        elif i == 1:
            nids_array[:, -1] = nids_set
        elif i == 2:
            nids_array[-1, ::-1] = nids_set
        else: #i=3
            # this column is added in reverse
            #nids_set = list(nids_set)
            #nids_set.reverse()
            nids_array[::-1, 0] = nids_set
        if debug:
            print(f'{flag}i={i} nids={nids_set} edge={edge} fwd_edge={fwd_edge}')
            print(nids_array)
            print('---------')

    sorted_face = _sort_face(face)
    ids = faces_to_center[sorted_face]
    if nnodes_to_add_with_ends == 3 and 1:
        nids_array[1, 1] = ids[0]
        #xyzc = (xyz1 + xyz2 + xyz3 + xyz4) / 4.
        #nodes[nid0] = GRID(nid0, xyzc)
        #nid0 += 1
    else:
        raise RuntimeError(ids)
        if nids_array.min() == 0:
            i, j = np.where(nids_array == 0)
            xi = np.linspace(0., 1., num=nnodes_to_add_with_ends, endpoint=True)#[1:-1]
            xj = np.linspace(0., 1., num=nnodes_to_add_with_ends, endpoint=True)#[1:-1]

            xii_ = xi[i]
            xjj_ = xj[j]
            # probably wrong...
            # bilinear interpolation
            for ii, jj, xii, xjj in zip(i, j, xii_, xjj_):
                xyz12 = xyz1 * (1. - xii) + xyz2 * xii
                xyz34 = xyz3 * (1. - xii) + xyz4 * xii # TODO: probably backwards
                xyzc = xyz12 * (1. - xjj) + xyz34 * xjj
                nodes[nid0] = GRID(nid0, xyzc)
                nids_array[ii, jj] = nid0
                nid0 += 1
        x = 1

    if debug:
        print('nids_array1:\n', nids_array)
    new_corner_nids = [
        nids_array[0, 0],
        nids_array[0, -1],
        nids_array[-1, -1],
        nids_array[-1, 0],
    ]
    unids2 = np.unique(new_corner_nids)
    assert len(unids2) == 4, unids2
    #print('nid0*** =', nid0)
    assert len(np.unique(nids_array)) == nids_array.size
    return nid0

def _insert_hexa_nodes(nodes: Dict[int, GRID],
                       nids_array, nid0: int,
                       edges, forward_edges,
                       edges_to_center,
                       faces_to_center, faces,
                       nnodes_to_add_with_ends: int,
                       xyz1, xyz2, xyz3, xyz4,
                       xyz5, xyz6, xyz7, xyz8,
                       debug=False):

    # apply edges_to_center
    for i, edge, fwd_edge in zip(count(), edges, forward_edges):
        nids_center = edges_to_center[edge]
        if debug:
            print('nids_center =', nids_center)
        flag = ' '
        if edge == fwd_edge:
            nids_set = nids_center
        else:
            #flag = '*'
            nids_set = list(nids_center)
            nids_set.reverse()
            if debug:
                print('nids_set =', nids_set)
        if i == 0:
            nids_array[:, 0, 0] = nids_set # [1,2]
        elif i == 1:
            nids_array[-1, :, 0] = nids_set #[2,3]
        elif i == 2:
            nids_array[::-1, -1, 0] = nids_set #[3,4]
        elif i == 3:
            # this column is added in reverse
            #nids_set = list(nids_set)
            #nids_set.reverse()
            nids_array[0, ::-1, 0] = nids_set #[4, 1]

        elif i == 4:
            nids_array[:, 0, -1] = nids_set
        elif i == 5:
            nids_array[-1, :, -1] = nids_set
        elif i == 6:
            nids_array[::-1, -1, -1] = nids_set
        elif i == 7:
            # this column is added in reverse
            #nids_set = list(nids_set)
            #nids_set.reverse()
            nids_array[0, :, -1] = nids_set

        # verticals
        elif i == 8:
            nids_array[0, 0, :] = nids_set
        elif i == 9:
            nids_array[-1, 0, :] = nids_set
        elif i == 10:
            nids_array[-1, -1, :] = nids_set
        elif i == 11:
            nids_array[0, -1, :] = nids_set
        else:
            raise NotImplementedError(i)

        if debug:
            print(f'{flag}i={i} nids={nids_set} edge={edge} fwd_edge={fwd_edge}')
            print(nids_array)
            print('---------')


    if nnodes_to_add_with_ends == 3:
        # apply faces_to_center
        for i, face in enumerate(faces):
            # bottom/top
            # left/right
            # front/back
            center_nids = faces_to_center[face] # length=1
            assert len(center_nids) == 1, center_nids
            center_nid = center_nids[0]
            if i == 0:
                plane = nids_array[:, :, 0]
            elif i == 1:
                plane = nids_array[:, :, -1]
            elif i == 2:
                plane = nids_array[:, 0, :]
            elif i == 3:
                plane = nids_array[:, -1, :]
            elif i == 4:
                plane = nids_array[0, :, :]
            elif i == 5:
                plane = nids_array[-1, :, :]
            else:
                asdf
            plane[1, 1] = center_nid

    else:
        raise NotImplementedError(nnodes_to_add_with_ends)

    #if nnodes_to_add_with_ends == 3 and 0:
        ##nids_array[1, 1] = nid0
        #xyzc = (xyz1 + xyz2 + xyz3 + xyz4 + xyz5 + xyz6 + xyz7 + xyz8) / 8.
        #nodes[nid0] = GRID(nid0, xyzc)
        #nid0 += 1
    #else:
    if nids_array.min() == 0:
        i, j, k = np.where(nids_array == 0)
        xarray = np.linspace(0., 1., num=nnodes_to_add_with_ends, endpoint=True)

        xi = xarray[i]
        xj = xarray[j]
        xk = xarray[k]
        # trilinear interpolation
        for ii, jj, kk, xii, xjj, xkk in zip(i, j, k, xi, xj, xk):
            xyz12 = xyz1 * (1. - xii) + xyz2 * xii
            xyz43 = xyz4 * (1. - xii) + xyz3 * xii
            xyz56 = xyz5 * (1. - xii) + xyz6 * xii
            xyz87 = xyz8 * (1. - xii) + xyz7 * xii

            xyz1234 = xyz12 * (1. - xjj) + xyz43 * xjj
            xyz5678 = xyz56 * (1. - xjj) + xyz87 * xjj

            xyz = xyz1234 * (1. - xkk) + xyz5678 * xkk
            nodes[nid0] = GRID(nid0, xyz)
            nids_array[ii, jj, kk] = nid0
            print(f'ijk=({ii},{jj},{kk}) nid={nid0} xyz={xyz}')
            nid0 += 1
        x = 1

    if debug:
        print('nids_array1:\n', nids_array)
    new_corner_nids = [
        nids_array[0, 0, 0],
        nids_array[0, -1, 0],
        nids_array[-1, -1, 0],
        nids_array[-1, 0, 0],

        nids_array[0, 0, -1],
        nids_array[0, -1, -1],
        nids_array[-1, -1, -1],
        nids_array[-1, 0, -1],
    ]
    unids2 = np.unique(new_corner_nids)
    assert len(unids2) == 8, unids2
    #print('nid0*** =', nid0)
    assert nids_array.min() >= 1, nids_array
    assert len(np.unique(nids_array)) == nids_array.size
    return nid0

def _quad_nids_to_node_ids(nids_array):
    assert nids_array.min() == 1, nids_array
    n1 = nids_array[:-1, :-1].ravel()
    n2 = nids_array[:-1, 1:].ravel()
    n3 = nids_array[1:, 1:].ravel()
    n4 = nids_array[1:, :-1].ravel()

    nnodes = len(n1)
    nodes = np.stack([n1, n2, n3, n4], axis=1)
    assert nodes.shape == (nnodes, 4), nodes.shape

    nodes = np.zeros((nnodes, 4), dtype=nids_array.dtype)
    nodes[:, 0] = n1
    nodes[:, 1] = n2
    nodes[:, 2] = n3
    nodes[:, 3] = n4
    return n1, n2, n3, n4

def _hexa_nids_to_node_ids(nids_array):
    assert nids_array.min() >= 1, nids_array
    n1 = nids_array[:-1, :-1, :-1].ravel()
    n2 = nids_array[1:, :-1,  :-1].ravel()
    n3 = nids_array[1:, 1:,   :-1].ravel()
    n4 = nids_array[:-1, 1:,  :-1].ravel()

    n5 = nids_array[:-1, :-1, 1:].ravel() # good
    n6 = nids_array[1:, :-1,  1:].ravel()
    n7 = nids_array[1:, 1:,   1:].ravel() # good
    n8 = nids_array[:-1, 1:,  1:].ravel()

    nodes = np.stack([n1, n2, n3, n4, n5, n6, n7, n8], axis=1)
    nnodes = len(n1)
    assert nodes.shape == (nnodes, 8), nodes.shape
    nodes = np.zeros((nnodes, 8), dtype=nids_array.dtype)
    nodes[:, 0] = n1
    nodes[:, 1] = n2
    nodes[:, 2] = n3
    nodes[:, 3] = n4

    nodes[:, 4] = n5
    nodes[:, 5] = n6
    nodes[:, 6] = n7
    nodes[:, 7] = n8
    return nodes

def test_insert_tri_nodes():
    #[[  10196 1206947   10184]
     #[      0 1206949 1206937]
     #[      0       0   10195]]
    nodes = np.array([
        [1, 4, 2],
        [0, 6, 5],
        [0, 0, 3]], dtype='int32')
    n1, n4, n2 = nodes[0, :]
    n6, n5 = nodes[1, 1:]
    n3 = nodes[2, 2]
    unids = np.unique([n1, n2, n3, n4, n5, n6])
    assert len(unids) == 6

def test_insert_quad_nodes():
    n1 = 1
    n2 = 2
    n3 = 3
    n4 = 4
    n5 = 5
    n6 = 6
    n7 = 7
    n8 = 8

    nodes = {
        1: np.array([0., 0., 0.]),
        2: np.array([1., 0., 0.]),
        3: np.array([1., 1., 0.]),
        4: np.array([0., 1., 0.]),
    }
    nid0 = max(nodes) + 1
    xyz1 = nodes[n1]
    xyz2 = nodes[n2]
    xyz3 = nodes[n3]
    xyz4 = nodes[n4]

    nnodes_to_add_with_ends = 3
    nids_array = np.zeros((nnodes_to_add_with_ends, nnodes_to_add_with_ends), dtype='int32')
    nids_array[0, 0] = n1
    nids_array[0, -1] = n2
    nids_array[-1, -1] = n3
    nids_array[-1, 0] = n4
    forward_edges = [(n1, n2), (n2, n3), (n3, n4), (n4, n1)]
    edges = forward_edges

    edges_to_center = {
        (n1, n2) : [n1, n5, n2],
        (n2, n3) : [n2, n6, n3],
        (n3, n4) : [n3, n7, n4],
        (n4, n1) : [n4, n8, n1],
    }
    nid0 = _insert_quad_nodes(
        nodes, nids_array, nid0,
        edges, forward_edges,
        edges_to_center, faces_to_center,
        nnodes_to_add_with_ends,
        n1, n2, n3, n4,
        xyz1, xyz2, xyz3, xyz4,
    )

def test_quad_nids_to_node_ids():
    nids = np.array([
        [1, 5, 2],
        [8, 9, 6],
        [4, 7, 3],
    ])
    n1, n2, n3, n4 = _quad_nids_to_node_ids(nids)
    assert np.array_equal(n1, [1, 5, 8, 9]), n1
    assert np.array_equal(n2, [5, 2, 9, 6]), n2
    assert np.array_equal(n3, [9, 6, 7, 3]), n3
    assert np.array_equal(n4, [8, 9, 4, 7]), n4

def test_hexa_nids_to_node_ids():
    nids = np.array([
        [
            [ 1, 17,  5],
            [ 9, 21, 13],
            [ 2, 18,  6]],

        [
            [12, 22, 16],
            [23, 24, 25],
            [10, 26, 14]],

        [
            [ 4, 20,  8],
            [11, 27, 15],
            [ 3, 19,  7]]]
    )
    nodes = _hexa_nids_to_node_ids(nids)
    n1 = nodes[:, 0]
    n2 = nodes[:, 1]
    n3 = nodes[:, 2]
    n4 = nodes[:, 3]
    n5 = nodes[:, 4]
    n6 = nodes[:, 5]
    n7 = nodes[:, 6]
    n8 = nodes[:, 7]
    un1 = np.unique(n1)
    un2 = np.unique(n2)
    un3 = np.unique(n3)
    un4 = np.unique(n4)
    un5 = np.unique(n5)
    un6 = np.unique(n6)
    un7 = np.unique(n7)
    un8 = np.unique(n8)
    assert np.array_equal(un1, np.unique([1, 9, 23, 12,
                                          17, 22, 24, 21])), un1
    assert np.array_equal(un2, np.unique([9, 2, 10, 23,
                                          22, 18, 27, 24])), un2
    assert np.array_equal(un3, np.unique([23, 10, 3, 11,
                                          24, 27, 19, 26])), un3
    assert np.array_equal(un4, np.unique([12, 23, 11, 4, 21, 24, 26, 20])), un4
    #assert np.array_equal(n1, [1, 5, 8, 9]), n1
    #assert np.array_equal(n2, [5, 2, 9, 6]), n2
    #assert np.array_equal(n3, [9, 6, 7, 3]), n3
    #assert np.array_equal(n4, [8, 9, 4, 7]), n4
    x = 1

def test_tri():
    import os
    import pyNastran
    pkg_path = pyNastran.__path__[0]
    model_path = os.path.join(pkg_path, '..', 'models')
    bwb_path = os.path.join(model_path, 'bwb')
    #bdf_filename = os.path.join(bwb_path, 'bwb_saero.bdf')
    bdf_filename_out = os.path.join(bwb_path, 'tri.bdf')

    model = BDF()
    model.add_grid(1, [0., 0., 0.])
    model.add_grid(2, [1., 0., 0.])
    model.add_grid(3, [1., 1., 0.])
    model.add_grid(4, [0., 1., 0.])
    model.add_ctria3(1, 1, [1, 2, 3])
    model.add_ctria3(2, 1, [1, 3, 4])
    model.add_pshell(1, mid1=1, t=0.1)
    model.add_mat1(1, 3.0e7, None, 0.3)
    ntris = len(model.elements)

    model = refine_model(model, refinement_ratio=2)
    model.write_bdf(bdf_filename_out)
    model.validate()
    assert len(model.elements) == ntris * 4
    model.cross_reference()
    x = 1

def test_quad():
    import os
    import pyNastran
    pkg_path = pyNastran.__path__[0]
    model_path = os.path.join(pkg_path, '..', 'models')
    bwb_path = os.path.join(model_path, 'bwb')
    bdf_filename_out = os.path.join(bwb_path, 'quad.bdf')

    model = BDF()
    model.add_grid(1, [0., 0., 0.])
    model.add_grid(2, [1., 0., 0.])
    model.add_grid(3, [1., 1., 0.])
    model.add_grid(4, [0., 1., 0.])
    model.add_cquad4(1, 1, [1, 2, 3, 4])
    model.add_pshell(1, mid1=1, t=0.1)
    model.add_mat1(1, 3.0e7, None, 0.3)
    nquads = len(model.elements)

    model = refine_model(model, refinement_ratio=2)
    model.write_bdf(bdf_filename_out)
    model.validate()
    assert len(model.elements) == nquads * 4
    model.cross_reference()
    x = 1

def test_hexa():
    import os
    import pyNastran
    pkg_path = pyNastran.__path__[0]
    model_path = os.path.join(pkg_path, '..', 'models')
    bwb_path = os.path.join(model_path, 'bwb')
    bdf_filename_out = os.path.join(bwb_path, 'hexa.bdf')

    model = BDF()
    model.add_grid(1, [0., 0., 0.])
    model.add_grid(2, [1., 0., 0.])
    model.add_grid(3, [1., 1., 0.])
    model.add_grid(4, [0., 1., 0.])

    model.add_grid(5, [0., 0., 1.])
    model.add_grid(6, [1., 0., 1.])
    model.add_grid(7, [1., 1., 1.])
    model.add_grid(8, [0., 1., 1.])

    model.add_chexa(1, 1, [1, 2, 3, 4, 5, 6, 7, 8])
    model.add_psolid(1, 1)
    model.add_mat1(1, 3.0e7, None, 0.3)
    nhexas = len(model.elements)

    model = refine_model(model, refinement_ratio=2)
    model.write_bdf(bdf_filename_out)
    model.validate()
    assert len(model.elements) == nhexas * 8
    model.cross_reference()
    x = 1

def test_hexa2():
    import os
    import pyNastran
    pkg_path = pyNastran.__path__[0]
    model_path = os.path.join(pkg_path, '..', 'models')
    bwb_path = os.path.join(model_path, 'bwb')
    bdf_filename_out = os.path.join(bwb_path, 'hexa2.bdf')

    model = BDF()
    model.add_grid(1, [0., 0., 0.])
    model.add_grid(2, [1., 0., 0.])
    model.add_grid(3, [1., 1., 0.])
    model.add_grid(4, [0., 1., 0.])

    model.add_grid(5, [0., 0., 1.])
    model.add_grid(6, [1., 0., 1.])
    model.add_grid(7, [1., 1., 1.])
    model.add_grid(8, [0., 1., 1.])

    model.add_grid(15, [0., 0., 2.])
    model.add_grid(16, [1., 0., 2.])
    model.add_grid(17, [1., 1., 2.])
    model.add_grid(18, [0., 1., 2.])

    model.add_chexa(1, 1, [1, 2, 3, 4,
                           5, 6, 7, 8])
    model.add_chexa(2, 1, [5, 6, 7, 8,
                           15, 16, 17, 18])
    model.add_psolid(1, 1)
    model.add_mat1(1, 3.0e7, None, 0.3)
    nhexas = len(model.elements)

    model = refine_model(model, refinement_ratio=2)
    model.write_bdf(bdf_filename_out)
    model.validate()
    assert len(model.elements) == nhexas * 8
    model.cross_reference()
    x = 1

def test_quad2():
    import os
    import pyNastran
    pkg_path = pyNastran.__path__[0]
    model_path = os.path.join(pkg_path, '..', 'models')
    bwb_path = os.path.join(model_path, 'bwb')
    bdf_filename_out = os.path.join(bwb_path, 'quad2.bdf')

    model = BDF()
    model.add_grid(1, [0., 0., 0.])
    model.add_grid(2, [1., 0., 0.])
    model.add_grid(3, [1., 1., 0.])
    model.add_grid(4, [0., 1., 0.])
    model.add_grid(5, [1., 2., 0.])
    model.add_grid(6, [0., 2., 0.])
    model.add_cquad4(1, 1, [1, 2, 3, 4])
    model.add_cquad4(2, 1, [4, 3, 5, 6])
    model.add_pshell(1, mid1=1, t=0.1)
    model.add_mat1(1, 3.0e7, None, 0.3)
    nquads = len(model.elements)

    model = refine_model(model, refinement_ratio=2)
    model.write_bdf(bdf_filename_out)
    model.validate()
    assert len(model.elements) == nquads * 4
    model.cross_reference()
    x = 1

def test_split_quad2():
    nids_array =np.array([
        [1, 7, 2],
        [10, 14, 8],
        [4, 9, 3], ], dtype='int32')
    n1, n2, n3, n4 = _quad_nids_to_node_ids(nids_array)
    x = 1

def test_refine():
    import os
    import pyNastran
    pkg_path = pyNastran.__path__[0]
    model_path = os.path.join(pkg_path, '..', 'models')
    bwb_path = os.path.join(model_path, 'bwb')
    bdf_filename = os.path.join(bwb_path, 'bwb_saero.bdf')
    bdf_filename_out = os.path.join(bwb_path, 'bwb_saero_fine.bdf')

    #bdf_filename = os.path.join(model_path, 'plate', 'plate.bdf')
    #model = BDF()
    #model.add_grid(1, [0., 0., 0.])
    #model.add_grid(2, [1., 0., 0.])
    #model.add_grid(3, [1., 1., 0.])
    #model.add_grid(4, [0., 1., 0.])
    #model.add_cquad4(1, 1, [1, 2, 3, 4])
    #model.add_pshell(1, mid1=1, t=0.1)
    #model.add_mat1(1, 3.0e7, None, 0.3)
    model = refine_model(bdf_filename, refinement_ratio=2)
    model.write_bdf(bdf_filename_out)
    model.validate()
    model.cross_reference()
    x = 1

if __name__ == '__main__':
    #test_insert_quad_nodes()
    #test_quad_nids_to_node_ids()
    #test_insert_tri_nodes()
    #test_tri()
    test_quad()
    #test_split_quad2()
    #test_quad2()
    #test_refine()
    #test_hexa_nids_to_node_ids()
    test_hexa()
    test_hexa2()
