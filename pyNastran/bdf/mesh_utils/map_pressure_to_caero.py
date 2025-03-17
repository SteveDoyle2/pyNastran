import os
import numpy as np
from pyNastran.utils import PathLike
from pyNastran.bdf.bdf import BDF, read_bdf

def map_caero(bdf_filename: PathLike):
    model = read_bdf(bdf_filename, xref=False)
    model.xref_obj.cross_reference_coordinates()
    model.xref_obj.cross_reference_nodes()
    neid = len(model.elements)
    nodes, xyz = model_to_node_xyz(model)
    out = model_to_tri_quads(model, nodes, xyz)
    quad_centroid = out['CQUAD4']['centroid']
    x = quad_centroid[:, 0]
    z = quad_centroid[:, 2]
    quad_pressure = np.sin(x/1000)*10 + z * 100.

    #-----------------------------------------------
    base, ext = os.path.splitext(bdf_filename)
    caero_bdf_filename = f'{base}.caero{ext}'
    from pyNastran.bdf.mesh_utils.export_caero_mesh import export_caero_mesh
    export_caero_mesh(model, caero_bdf_filename, is_subpanel_model=True,
                      write_panel_xyz=False)
    caero_model = read_bdf(caero_bdf_filename, xref=False)
    #-----------------------------------------------
    caero_nodes, caero_xyz = model_to_node_xyz(caero_model)
    quad_to_tri_model(caero_model)
    out_caero = model_to_tri_quads(caero_model, caero_nodes, caero_xyz)
    tri = out_caero['CTRIA3']
    tri_nodes = tri['nodes']
    assert len(tri_nodes) > 0
    itri_nodes = np.searchsorted(caero_nodes, tri_nodes)
    tri_centroid = tri['centroid']
    n1 = caero_xyz[itri_nodes[:, 0], :]
    n2 = caero_xyz[itri_nodes[:, 1], :]
    n3 = caero_xyz[itri_nodes[:, 2], :]
    assert len(tri_centroid) > 0

def quad_to_tri_model(model: BDF) -> None:
    eid = 1
    elements = model.elements
    model.elements = {}
    for _eid, elem in elements.items():
        if elem.type == 'CQUAD4':
            pid = elem.pid
            (n1, n2, n3, n4) = elem.nodes
            model.add_ctria3(eid, pid, [n1, n2, n3])
            model.add_ctria3(eid+1, pid, [n1, n3, n4])
            eid += 2
        elif elem.type == 'CTRIA3':
            pid = elem.pid
            nodes = elem.nodes
            model.add_ctria3(eid, pid, nodes)
            eid += 1
    return

def model_to_node_xyz(model: BDF):
    xyz_list = []
    node_list = []
    for nid, node in sorted(model.nodes.items()):
        node_list.append(nid)
        xyz_list.append(node.get_position())
    nodes = np.array(node_list, dtype='int32')
    xyz = np.array(xyz_list, dtype='float64')
    return nodes, xyz

def model_to_tri_quads(model: BDF,
                       nodes: np.ndarray,
                       xyz: np.ndarray,
                       ) -> dict[str, dict]:
    tris_list = []
    quads_list = []
    tri_nodes_list = []
    quad_nodes_list = []
    for i, (eid, elem) in enumerate(sorted(model.elements.items())):
        if elem.type == 'CTRIA3':
            tri_nodes_list.append(elem.nodes)
            tris_list.append(eid)
        elif elem.type == 'CQUAD4':
            quad_nodes_list.append(elem.nodes)
            quads_list.append(eid)
    tri_nodes = np.array(tri_nodes_list, dtype='int32')
    quad_nodes = np.array(quad_nodes_list, dtype='int32')
    tris = np.array(tris_list, dtype='int32')
    quads = np.array(quads_list, dtype='int32')

    nquad = len(quads)
    ntri = len(tris)
    out = {}
    if ntri:
        itri_nodes = np.searchsorted(nodes, tri_nodes)
        tri_xyz = xyz[itri_nodes]
        n1 = tri_xyz[:, 0, :]
        n2 = tri_xyz[:, 1, :]
        n3 = tri_xyz[:, 2, :]
        tri_cross = np.cross(n2 - n1, n3 - n1, axis=1)
        normi = np.linalg.norm(tri_cross, axis=1)
        tri_area = 0.5 * normi
        tri_normal = tri_cross / normi[:, np.newaxis]
        tri_centroid = tri_xyz.mean(axis=1)
        assert tri_centroid.shape == (ntri, 3), (ntri, tri_xyz.shape)
        out['CTRIA3'] = {
            'nodes': tri_nodes,
            'centroid': tri_centroid,
            'area': tri_area,
            'normal': tri_normal,
        }

    if nquad:
        iquad_nodes = np.searchsorted(nodes, quad_nodes)
        quad_xyz = xyz[iquad_nodes]
        n1 = quad_xyz[:, 0, :]
        n2 = quad_xyz[:, 1, :]
        n3 = quad_xyz[:, 2, :]
        n4 = quad_xyz[:, 3, :]
        quad_cross = np.cross(n3 - n1, n4 - n2, axis=1)
        normi = np.linalg.norm(quad_cross, axis=1)
        quad_area = 0.5 * normi
        quad_normal = quad_cross / normi[:, np.newaxis]

        # quads: (9236, 4, 3) -> (9236, 3)
        quad_centroid = quad_xyz.mean(axis=1)
        assert quad_centroid.shape == (nquad, 3), (nquad, quad_xyz.shape)
        out['CQUAD4'] = {
            'nodes': quad_nodes,
            'centroid': quad_centroid,
            'area': quad_area,
            'normal': quad_normal,
        }
    return out

def point_in_trangle(xyz, xyz1, xyz2, xyz3):
    x = xyz[:, 0]
    y = xyz[:, 1]
    x1 = xyz1[:, 0]
    x2 = xyz2[:, 0]
    x3 = xyz3[:, 0]
    y1 = xyz1[:, 1]
    y2 = xyz2[:, 1]
    y3 = xyz3[:, 1]

    #Calculate u, v using the following formulas:
    u = ((x2-x3)(y-y1) - (x-x1)(y2-y1)) / ((x2-x3)(y2-y1) - (x2-x1)(y3-y1))
    v = ((x3-x1)(y-y2) - (x-x2)(y3-y2)) / ((x2-x3)(y2-y1) - (x2-x1)(y3-y1))
    #If u >= 0, v >= 0, and 1-u-v >= 0, then P is inside the triangle.
