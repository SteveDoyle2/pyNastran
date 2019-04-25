# pylint: disable=C0103
"""
defines NastranGuiAttributes, which defines
GUI specific geometry functions that don't involve PyQt/VTK
this is no longer true...but should be
"""
from __future__ import print_function
import sys
from collections import defaultdict
import numpy as np
from numpy.linalg import norm
import vtk
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.cards.elements.beam_connectivity import (
    rod_faces, tube_faces, chan1_faces,
    bar_faces, box_faces, i_faces, t_faces, t1_faces, t2_faces,
    h_faces, i1_faces, chan_faces, l_faces, z_faces,
    hexa_faces, hat_faces,
)
from pyNastran.bdf.cards.elements.bars import rotate_v_wa_wb
from pyNastran.gui.utils.vtk.vtk_utils import numpy_to_vtk_points, numpy_to_vtk


from pyNastran.gui.qt_files.colors import BLUE_FLOAT
BEAM_GEOM_TYPES = [
    'BAR', 'BOX', 'BOX1', 'CHAN', 'CHAN1', 'CHAN2', 'CROSS', 'DBOX',
    'H', 'HAT', 'HAT1', 'HEXA', 'I', 'I1', 'L', 'ROD',
    'T', 'T1', 'T2', 'TUBE', 'TUBE2', 'Z',
]


class NastranGuiAttributes(object):
    """GUI specific geometry functions that don't involve PyQt/VTK"""
    def __init__(self):
        # new options, no way to access them through the gui
        # they control results generation
        self.make_xyz = False
        self.make_offset_normals_dim = True # creates normals, which is required for loads
        self.make_nnodes_result = False  # make_offset_normals_dim must be True for this to work
        self.make_released_dofs1 = False
        self.make_released_dofs2 = False
        self.plot_applied_loads = True
        self.plot_pressures = True  # make_offset_normals_dim must be True for this to work

        # these are local variables
        self.model = None
        self.model_results = None
        self.bar_lines = None
        self.bar_eids = None
        self.xyz_cid0 = None
        self.normals = None
        self.save_data = True # was False
        #--------------------------------------

        #: flips the nastran CAERO subpaneling
        #:   False -> borders of CAEROs can be seen
        #:   True  -> individual subpanels can be seen
        self.show_caero_sub_panels = False

        #: coordinate systems can be messy, so this is the
        #: list of coords to show
        self.show_cids = []

        self.show_caero_actor = True  # show the caero mesh
        self.show_control_surfaces = True
        self.show_conm = True

        self.element_ids = None
        self.node_ids = None
        self.nid_map = None
        self.eid_map = None
        self.nnodes = None
        self.nelements = None
        self.model_type = None
        self.isubcase_name_map = None
        self.has_caero = False
        self.dependents_nodes = set()
        self.icd_transform = {}


class NastranGeometryHelper(NastranGuiAttributes):
    """
    Defines VTK/PyQt-less methods used by NastranIO-Geometry
    """
    def __init__(self):
        super(NastranGeometryHelper, self).__init__()

    def _get_bar_yz_arrays(self, model, bar_beam_eids, scale, debug):
        lines_bar_y = []
        lines_bar_z = []
        points_list = []

        bar_types = {
            # PBAR
            'bar' : [],

            # PBEAML/PBARL
            "ROD": [],
            "TUBE": [],
            "TUBE2" : [],
            "I": [],
            "CHAN": [],
            "T": [],
            "BOX": [],
            "BAR": [],
            "CROSS": [],
            "H": [],
            "T1": [],
            "I1": [],
            "CHAN1": [],
            "Z": [],
            "CHAN2": [],
            "T2": [],
            "BOX1": [],
            "HEXA": [],
            "HAT": [],
            "HAT1": [],
            "DBOX": [],  # was 12

            # PBEAM
            'beam' : [],

            # PBEAML specfic
            "L" : [],
        }  # for GROUP="MSCBML0"
        allowed_types = [
            'BAR', 'BOX', 'BOX1', 'CHAN', 'CHAN1', 'CHAN2', 'CROSS', 'DBOX',
            'H', 'HAT', 'HAT1', 'HEXA', 'I', 'I1', 'L', 'ROD',
            'T', 'T1', 'T2', 'TUBE', 'TUBE2', 'Z', 'bar', 'beam',
        ]

        # bar_types['bar'] = [ [...], [...], [...] ]
        #bar_types = defaultdict(lambda : defaultdict(list))

        found_bar_types = set()
        #neids = len(self.element_ids)
        for bar_type, data in bar_types.items():
            eids = []
            lines_bar_y = []
            lines_bar_z = []
            bar_types[bar_type] = (eids, lines_bar_y, lines_bar_z)
            #bar_types[bar_type] = [eids, lines_bar_y, lines_bar_z]

        ugrid = vtk.vtkUnstructuredGrid()
        node0 = 0

        nid_release_map = defaultdict(list)

        #debug = True
        bar_nids = set()
        #print('bar_beam_eids = %s' % bar_beam_eids)
        for eid in bar_beam_eids:
            if eid not in self.eid_map:
                self.log.error('eid=%s is not a valid bar/beam element...' % eid)
                if debug:  # pragma: no cover
                    print('eid=%s is not a valid bar/beam element...' % eid)
                continue
            #unused_ieid = self.eid_map[eid]
            elem = model.elements[eid]
            pid_ref = elem.pid_ref
            if pid_ref is None:
                pid_ref = model.Property(elem.pid)
            assert not isinstance(pid_ref, integer_types), elem

            ptype = pid_ref.type
            bar_type = _get_bar_type(ptype, pid_ref)

            if debug:  # pragma: no cover
                print('%s' % elem)
                print('  bar_type = %s' % bar_type)
            found_bar_types.add(bar_type)

            (nid1, nid2) = elem.node_ids
            bar_nids.update([nid1, nid2])
            node1 = model.nodes[nid1]
            node2 = model.nodes[nid2]
            n1 = node1.get_position()
            n2 = node2.get_position()

            # wa/wb are not considered in i_offset
            # they are considered in ihat
            i = n2 - n1
            Li = norm(i)
            ihat = i / Li

            if elem.pa != 0:
                nid_release_map[nid1].append((eid, elem.pa))
            if elem.pb != 0:
                nid_release_map[nid2].append((eid, elem.pb))

            unused_v, wa, wb, xform = rotate_v_wa_wb(
                model, elem,
                n1, n2, node1, node2,
                ihat, i, eid, Li, model.log)
            if wb is None:
                # one or more of v, wa, wb are bad
                continue

            yhat = xform[1, :]
            zhat = xform[2, :]

            ## concept has a GOO

            #if debug:  # pragma: no cover
                #print('  centroid = %s' % centroid)
                #print('  ihat = %s' % ihat)
                #print('  yhat = %s' % yhat)
                #print('  zhat = %s' % zhat)
                #print('  scale = %s' % scale)
            #if eid == 616211:
                #print('  check - eid=%s yhat=%s zhat=%s v=%s i=%s n%s=%s n%s=%s' % (
                      #eid, yhat, zhat, v, i, nid1, n1, nid2, n2))

                #print('adding bar %s' % bar_type)
                #print('   centroid=%s' % centroid)
                #print('   yhat=%s len=%s' % (yhat, np.linalg.norm(yhat)))
                #print('   zhat=%s len=%s' % (zhat, np.linalg.norm(zhat)))
                #print('   Li=%s scale=%s' % (Li, scale))
            if bar_type not in allowed_types:
                msg = 'bar_type=%r allowed=[%s]' % (bar_type, ', '.join(allowed_types))
                raise RuntimeError(msg)

            if bar_type in BEAM_GEOM_TYPES:
                node0 = add_3d_bar_element(
                    bar_type, ptype, pid_ref,
                    n1+wa, n2+wb, xform,
                    ugrid, node0, points_list)

            centroid = (n1 + n2) / 2.
            bar_types[bar_type][0].append(eid)
            bar_types[bar_type][1].append((centroid, centroid + yhat * Li * scale))
            bar_types[bar_type][2].append((centroid, centroid + zhat * Li * scale))

        if node0: # and '3d_bars' not in self.alt_grids:
            def update_grid_function(unused_nid_map, ugrid, points, nodes):  # pragma: no cover
                """custom function to update the 3d bars"""
                points_list = []
                node0b = 0
                for eid in bar_beam_eids:
                    elem = self.model.elements[eid]
                    pid_ref = elem.pid_ref
                    if pid_ref is None:
                        pid_ref = self.model.Property(elem.pid)
                    assert not isinstance(pid_ref, integer_types), elem

                    ptype = pid_ref.type
                    bar_type = _get_bar_type(ptype, pid_ref)

                    #nids = elem.nodes
                    (nid1, nid2) = elem.node_ids
                    node1 = model.nodes[nid1]
                    node2 = model.nodes[nid2]

                    i1, i2 = np.searchsorted(self.node_ids, [nid1, nid2])
                    n1 = nodes[i1, :]
                    n2 = nodes[i2, :]
                    #centroid = (n1 + n2) / 2.

                    i = n2 - n1
                    Li = norm(i)
                    ihat = i / Li

                    unused_v, wa, wb, xform = rotate_v_wa_wb(
                        model, elem,
                        n1, n2, node1, node2,
                        ihat, i, eid, Li, model.log)
                    if wb is None:
                        # one or more of v, wa, wb are bad
                        continue

                    ugridi = None
                    node0b = add_3d_bar_element(
                        bar_type, ptype, pid_ref,
                        n1+wa, n2+wb, xform,
                        ugridi, node0b, points_list, add_to_ugrid=False)

                points_array = _make_points_array(points_list)

                points_array2 = numpy_to_vtk(
                    num_array=points_array,
                    deep=1,
                    array_type=vtk.VTK_FLOAT,
                )
                points.SetData(points_array2)

                ugrid.SetPoints(points)
                points.Modified()
                ugrid.Modified()
                return

            if points_list:
                if not sys.argv[0].startswith('test_'):
                    update_grid_function = None
                self.gui.create_alternate_vtk_grid(
                    '3d_bars', color=BLUE_FLOAT, opacity=0.2,
                    representation='surface', is_visible=True,
                    follower_function=update_grid_function,
                    ugrid=ugrid,
                )
                points_array = _make_points_array(points_list)
                points = numpy_to_vtk_points(points_array)
                ugrid.SetPoints(points)

        #print('bar_types =', bar_types)
        for bar_type in list(bar_types):
            bars = bar_types[bar_type]
            if len(bars[0]) == 0:
                del bar_types[bar_type]
                continue
            #bar_types[bar_type][1] = np.array(bars[1], dtype='float32')  # lines_bar_y
            #bar_types[bar_type][2] = np.array(bars[2], dtype='float32')  # lines_bar_z

        debug = False
        if debug:  # pragma: no cover
            #np.set_printoptions(formatter={'float': '{: 0.3f}'.format})
            for bar_type, data in sorted(bar_types.items()):
                eids, lines_bar_y, lines_bar_z = data
                if len(eids):
                    #print('barsi =', barsi)
                    #print('bar_type = %r' % bar_type)
                    for eid, line_y, line_z  in zip(eids, lines_bar_y, lines_bar_z):
                        print('eid=%s centroid=%s cy=%s cz=%s' % (
                            eid, line_y[0], line_y[1], line_z[1]))

        #print('found_bar_types =', found_bar_types)
        #no_axial_torsion = (no_axial, no_torsion)
        #no_shear_bending = (no_shear_y, no_shear_z, no_bending_y, no_bending_z)
        #no_dofs = (no_bending, no_bending_bad, no_6_16, no_0_456,
                   #no_0_56, no_56_456, no_0_6, no_0_16)
        return bar_nids, bar_types, nid_release_map

def _make_points_array(points_list):
    if len(points_list) == 1:
        points_array = points_list[0]
    else:
        points_array = np.vstack(points_list)
    return points_array

def _apply_points_list(points_list, ugrid):
    if points_list:
        points_array = _make_points_array(points_list)
        points = numpy_to_vtk_points(points_array)
        ugrid.SetPoints(points)

def _get_bar_type(ptype, pid_ref):
    """helper method for _get_bar_yz_arrays"""
    if ptype in ['PBAR', 'PBEAM']:
        bar_type = 'bar'
    elif ptype in ['PBEAM']:
        bar_type = 'beam'
    elif ptype in ['PBARL', 'PBEAML']:
        bar_type = pid_ref.Type
    else:
        raise NotImplementedError(pid_ref)
    return bar_type

def get_bar_yz_transform(v, ihat, eid, n1, n2, nid1, nid2, i, Li):
    """helper method for _get_bar_yz_arrays"""
    vhat = v / norm(v) # j
    try:
        z = np.cross(ihat, vhat) # k
    except ValueError:
        msg = 'Invalid vector length\n'
        msg += 'n1  =%s\n' % str(n1)
        msg += 'n2  =%s\n' % str(n2)
        msg += 'nid1=%s\n' % str(nid1)
        msg += 'nid2=%s\n' % str(nid2)
        msg += 'i   =%s\n' % str(i)
        msg += 'Li  =%s\n' % str(Li)
        msg += 'ihat=%s\n' % str(ihat)
        msg += 'v   =%s\n' % str(v)
        msg += 'vhat=%s\n' % str(vhat)
        msg += 'z=cross(ihat, vhat)'
        print(msg)
        raise ValueError(msg)

    zhat = z / norm(z)
    yhat = np.cross(zhat, ihat) # j

    if norm(ihat) == 0.0 or norm(yhat) == 0.0 or norm(z) == 0.0:
        print('  invalid_orientation - eid=%s yhat=%s zhat=%s v=%s i=%s n%s=%s n%s=%s' % (
            eid, yhat, zhat, v, i, nid1, n1, nid2, n2))
    elif not np.allclose(norm(yhat), 1.0) or not np.allclose(norm(zhat), 1.0) or Li == 0.0:
        print('  length_error        - eid=%s Li=%s Lyhat=%s Lzhat=%s'
              ' v=%s i=%s n%s=%s n%s=%s' % (
                  eid, Li, norm(yhat), norm(zhat), v, i, nid1, n1, nid2, n2))
    return yhat, zhat

def get_suport_node_ids(model, suport_id):
    """gets the nodes where SUPORTs and SUPORT1s are defined"""
    node_ids = []
    # list
    #for suport in model.suport:
        #node_ids += suport.IDs

    # dict
    if suport_id in model.suport1:
        suport1 = model.suport1[suport_id]
        node_ids += suport1.nodes
    else:
        # TODO: shouldn't this block always be included?
        for suport in model.suport:
            if suport_id in suport.nodes:
                node_ids.append(suport_id)
    return np.unique(node_ids)

def get_material_arrays(model, mids):
    """gets e11, e22, e33"""
    #e11 = np.zeros(mids.shape, dtype='float32')
    #e22 = np.zeros(mids.shape, dtype='float32')
    #e33 = np.zeros(mids.shape, dtype='float32')
    #rho = np.zeros(mids.shape, dtype='float32')
    #bulk = np.zeros(mids.shape, dtype='float32')
    #speed_of_sound = np.zeros(mids.shape, dtype='float32')

    e11 = np.full(mids.shape, np.nan, dtype='float32')
    e22 = np.full(mids.shape, np.nan, dtype='float32')
    e33 = np.full(mids.shape, np.nan, dtype='float32')
    rho = np.full(mids.shape, np.nan, dtype='float32')
    bulk = np.full(mids.shape, np.nan, dtype='float32')
    speed_of_sound = np.full(mids.shape, np.nan, dtype='float32')

    has_mat8 = False
    #has_mat10 = False
    has_mat11 = False
    materials = model.materials
    for superelement in model.superelement_models.values():
        materials.update(superelement.materials)
    for umid in np.unique(mids):
        if umid == 0:
            continue
        e11i = e22i = e33i = 0.
        rhoi = 0.
        bulki = 0.
        speed_of_soundi = 0.
        try:
            mat = materials[umid]
        except KeyError:
            print("can't find mid=%s" % umid)
            print('  mids = %s' % mids)
            print('  mids = %s' % materials.keys())
            continue
            #raise
        if mat.type == 'MAT1':
            e11i = e22i = e33i = mat.e
            rhoi = mat.rho
        elif mat.type == 'MAT8':
            e11i = e33i = mat.e11
            e22i = mat.e22
            has_mat8 = True
            rhoi = mat.rho
        #elif mat.type == 'MAT9':
            # Defines the material properties for linear, temperature-independent,
            #anisotropic materials for solid isoparametric elements (PSOLID)
            #g11i = mat.G11
        elif mat.type in ['MAT11', 'MAT3D']:
            e11i = mat.e1
            e22i = mat.e2
            e33i = mat.e3
            has_mat11 = True
            rhoi = mat.rho
        elif mat.type == 'MAT10':
            bulki = mat.bulk
            rhoi = mat.rho
            speed_of_soundi = mat.c
            #has_mat10 = True
            #self.log.info('skipping\n%s' % mat)
            #continue
        else:
            msg = 'skipping\n%s' % mat
            #print(msg.encode('utf16'))
            print(msg.encode(model._encoding))
            continue
            #raise NotImplementedError(mat)
        #print('mid=%s e11=%e e22=%e' % (umid, e11i, e22i))
        i = np.where(umid == mids)[0]
        e11[i] = e11i
        e22[i] = e22i
        e33[i] = e33i
        rho[i] = rhoi
        bulk[i] = bulki
        speed_of_sound[i] = speed_of_soundi
    return has_mat8, has_mat11, e11, e22, e33


def faces_to_element_facelist(faces, node0):
    """creates a series of faces for the custom elements"""
    face_idlist = vtk.vtkIdList()

    nfaces = len(faces)
    face_idlist.InsertNextId(nfaces) # Number faces that make up the cell.
    for face in faces: # Loop over all the faces
        #print(face)
        face_idlist.InsertNextId(len(face)) # Number of points in face

        # Insert the pointIds for the face
        #for i in face:
            #face_idlist.InsertNextId(i + node0)
        [face_idlist.InsertNextId(i + node0) for i in face]

    return face_idlist

def add_3d_bar_element(bar_type, ptype, pid_ref,
                       n1, n2, xform,
                       ugrid, node0, points_list, add_to_ugrid=True):
    """adds a 3d bar element to the unstructured grid"""
    if ptype in ['PBARL']:
        dim1 = dim2 = pid_ref.dim
        #bar_type = pid_ref.Type
    elif ptype in ['PBEAML']:
        dim1 = pid_ref.dim[0, :]
        dim2 = pid_ref.dim[-1, :]
    else:
        dim1 = dim2 = None
        return node0

    if bar_type == 'BAR':
        pointsi = bar_faces(n1, n2, xform, dim1, dim2)
        elem = vtk.vtkHexahedron()
        point_ids = elem.GetPointIds()
        point_ids.SetId(0, node0 + 0)
        point_ids.SetId(1, node0 + 1)
        point_ids.SetId(2, node0 + 2)
        point_ids.SetId(3, node0 + 3)
        point_ids.SetId(4, node0 + 4)
        point_ids.SetId(5, node0 + 5)
        point_ids.SetId(6, node0 + 6)
        point_ids.SetId(7, node0 + 7)
        if add_to_ugrid:
            ugrid.InsertNextCell(12, point_ids)
        points_list.append(pointsi)
        node0 += 8
        return node0
    elif bar_type == 'ROD':
        faces, pointsi, dnode = rod_faces(n1, n2, xform, dim1, dim2)
        face_idlist = faces_to_element_facelist(faces, node0)
        node0 += dnode
    elif bar_type == 'TUBE':
        faces, pointsi, dnode = tube_faces(n1, n2, xform, dim1, dim2)
        face_idlist = faces_to_element_facelist(faces, node0)
        node0 += dnode
    elif bar_type == 'BOX':
        faces, pointsi = box_faces(n1, n2, xform, dim1, dim2)
        face_idlist = faces_to_element_facelist(faces, node0)
        node0 += 16
    elif bar_type == 'L':
        faces, pointsi = l_faces(n1, n2, xform, dim1, dim2)
        face_idlist = faces_to_element_facelist(faces, node0)
        node0 += 12
    elif bar_type == 'CHAN':
        faces, pointsi = chan_faces(n1, n2, xform, dim1, dim2)
        face_idlist = faces_to_element_facelist(faces, node0)
        node0 += 16
    elif bar_type == 'CHAN1':
        faces, pointsi = chan1_faces(n1, n2, xform, dim1, dim2)
        face_idlist = faces_to_element_facelist(faces, node0)
        node0 += 16
    elif bar_type == 'T':
        faces, pointsi = t_faces(n1, n2, xform, dim1, dim2)
        face_idlist = faces_to_element_facelist(faces, node0)
        node0 += 16
    elif bar_type == 'T1':
        faces, pointsi = t1_faces(n1, n2, xform, dim1, dim2)
        face_idlist = faces_to_element_facelist(faces, node0)
        node0 += 16
    elif bar_type == 'T2':
        faces, pointsi = t2_faces(n1, n2, xform, dim1, dim2)
        face_idlist = faces_to_element_facelist(faces, node0)
        node0 += 16
    elif bar_type == 'I':
        faces, pointsi = i_faces(n1, n2, xform, dim1, dim2)
        assert pointsi.shape[0] == 24, pointsi.shape
        face_idlist = faces_to_element_facelist(faces, node0)
        node0 += 24
    elif bar_type == 'I1':
        faces, pointsi = i1_faces(n1, n2, xform, dim1, dim2)
        assert pointsi.shape[0] == 24, pointsi.shape
        face_idlist = faces_to_element_facelist(faces, node0)
        node0 += 24
    elif bar_type == 'H':
        faces, pointsi = h_faces(n1, n2, xform, dim1, dim2)
        assert pointsi.shape[0] == 24, pointsi.shape
        face_idlist = faces_to_element_facelist(faces, node0)
        node0 += 24
    elif bar_type == 'Z':
        faces, pointsi = z_faces(n1, n2, xform, dim1, dim2)
        face_idlist = faces_to_element_facelist(faces, node0)
        node0 += 16
    elif bar_type == 'HEXA':
        faces, pointsi = hexa_faces(n1, n2, xform, dim1, dim2)
        face_idlist = faces_to_element_facelist(faces, node0)
        node0 += 12
    elif bar_type == 'HAT':
        faces, pointsi = hat_faces(n1, n2, xform, dim1, dim2)
        face_idlist = faces_to_element_facelist(faces, node0)
        node0 += 24
    else:
        #print('skipping 3d bar_type = %r' % bar_type)
        return node0
    if add_to_ugrid:
        ugrid.InsertNextCell(vtk.VTK_POLYHEDRON, face_idlist)
    points_list.append(pointsi)
    return node0
