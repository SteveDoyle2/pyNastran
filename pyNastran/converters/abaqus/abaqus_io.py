"""
Defines how the GUI reads Abaqus files
"""
from __future__ import print_function
from collections import OrderedDict

import numpy as np

import vtk
from vtk import vtkLine, vtkTriangle, vtkQuad, vtkTetra, vtkHexahedron
from pyNastran.gui.utils.vtk.vtk_utils import numpy_to_vtk

from pyNastran.gui.gui_objects.gui_result import GuiResult, NormalResult
#from pyNastran.gui.qt_files.result import Result
from pyNastran.converters.abaqus.abaqus import Abaqus


class AbaqusIO(object):
    """Defines the GUI class for Abaqus."""
    def __init__(self, gui):
        self.gui = gui

    def get_abaqus_wildcard_geometry_results_functions(self):
        """dynamic named method for loading abaqus input files"""
        data = (
            'Abaqus',
            'Abaqus (*.inp)', self.load_abaqus_geometry,
            None, None
        )
        return data

    def load_abaqus_geometry(self, abaqus_filename, name='main', plot=True):
        """loads abaqus input files into the gui"""
        model_name = name
        skip_reading = self.gui._remove_old_geometry(abaqus_filename)
        if skip_reading:
            return

        self.gui.eid_maps[name] = {}
        self.gui.nid_maps[name] = {}
        model = Abaqus(log=self.gui.log, debug=False)
        self.gui.model_type = 'abaqus'
        #self.model_type = model.model_type
        model.read_abaqus_inp(abaqus_filename)

        self.gui.nid_map = {}
        nnodes, all_nodes, nelements = get_nodes_nnodes_nelements(model)
        self.gui.log.info('nnodes=%s nelements=%s' % (nnodes, nelements))
        assert nelements > 0, nelements
        #nodes = model.nodes
        #elements = model.elements


        self.gui.nnodes = nnodes
        self.gui.nelements = nelements

        grid = self.gui.grid
        grid.Allocate(self.gui.nelements, 1000)

        assert len(all_nodes) > 0 is not None, len(all_nodes)
        if len(all_nodes) == 1:
            nodes = all_nodes[0]
        else:
            nodes = np.vstack(all_nodes)

        mmax = np.amax(nodes, axis=0)
        mmin = np.amin(nodes, axis=0)
        dim_max = (mmax - mmin).max()
        self.gui.create_global_axes(dim_max)

        points = vtk.vtkPoints()
        points.SetNumberOfPoints(self.gui.nnodes)

        data_type = vtk.VTK_FLOAT
        points_array = numpy_to_vtk(
            num_array=nodes,
            deep=True,
            array_type=data_type
        )
        points.SetData(points_array)

        nid_offset = -1
        nids = []
        for unused_part_name, part in model.parts.items():
            self.gui.log.info('part_name = %r' % unused_part_name)
            nnodesi = part.nodes.shape[0]
            nidsi = part.nids
            nids.append(nidsi)

            add_lines(grid, nidsi, part.r2d2, nid_offset)

            add_tris(grid, nidsi, part.cps3, nid_offset)
            add_tris(grid, nidsi, part.cpe3, nid_offset)

            add_quads(grid, nidsi, part.cpe4, nid_offset)
            add_quads(grid, nidsi, part.cpe4r, nid_offset)
            add_quads(grid, nidsi, part.coh2d4, nid_offset)
            add_quads(grid, nidsi, part.cohax4, nid_offset)

            add_quads(grid, nidsi, part.cax4r, nid_offset)
            add_tris(grid, nidsi, part.cax3, nid_offset)

            # solids
            add_tetras(grid, nidsi, part.c3d10h, nid_offset)
            add_hexas(grid, nidsi, part.c3d8r, nid_offset)

            nid_offset += nnodesi
        nids = np.hstack(nids)
        grid.SetPoints(points)
        grid.Modified()

        # loadCart3dResults - regions/loads
        self.gui.scalar_bar_actor.VisibilityOn()
        self.gui.scalar_bar_actor.Modified()

        note = ''
        self.gui.isubcase_name_map = {1: ['Abaqus%s' % note, '']}
        #form = []
        cases = OrderedDict()
        ID = 1
        form, cases, unused_icase, node_ids, element_ids = self._fill_abaqus_case(
            cases, ID, nids, nodes, nelements, model)
        #self._fill_cart3d_results(cases, form, icase, ID, model)

        self.gui.node_ids = node_ids
        self.gui.element_ids = element_ids
        self.gui._finish_results_io2(model_name, form, cases)

    def clear_abaqus(self):
        """does nothing"""
        pass

    def load_abaqus_results(self, abaqusd_filename):
        """does nothing"""
        raise NotImplementedError()

    def _fill_abaqus_case(self, cases, ID, node_ids, nodes, nelements, unused_model):
        """creates the result objects for abaqus"""
        #return [], {}, 0
        #nelements = elements.shape[0]
        nnodes = nodes.shape[0]

        element_ids = np.arange(1, nelements + 1)
        #print(nodes)
        #node_ids = np.arange(1, nnodes + 1)
        #cnormals = model.get_normals(shift_nodes=False)
        #cnnodes = cnormals.shape[0]
        #assert cnnodes == nelements, len(cnnodes)

        #print('nnodes =', nnodes)
        #print('nelements =', nelements)
        #print('regions.shape =', regions.shape)
        #subcase_id = 0
        #labels = ['NodeID', 'ElementID']
        #cart3d_geo = Cart3dGeometry(subcase_id, labels,
                                    #nids, eids, regions, cnormals,
                                    #uname='Cart3dGeometry')
        colormap = 'jet'
        nid_res = GuiResult(ID, header='NodeID', title='NodeID',
                            location='node', scalar=node_ids)
        eid_res = GuiResult(ID, header='ElementID', title='ElementID',
                            location='centroid', scalar=element_ids)
        nxyz_res = NormalResult(0, 'Normals', 'Normals',
                                nlabels=2, labelsize=5, ncolors=2,
                                colormap=colormap, data_format='%.1f',
                                uname='NormalResult')

        cases[0] = (nid_res, (0, 'NodeID'))
        cases[1] = (eid_res, (0, 'ElementID'))
        cases[2] = (nxyz_res, (0, 'Normal'))

        geometry_form = [
            ('NodeID', 0, []),
            ('ElementID', 1, []),
            ('Normal', 2, []),
        ]
        form = [
            ('Geometry', None, geometry_form),
        ]
        icase = 2
        return form, cases, icase, node_ids, element_ids


def get_nodes_nnodes_nelements(model):
    """helper method"""
    n_r2d2 = 0

    n_cps3 = 0
    n_cpe3 = 0
    n_cpe4 = 0
    n_cpe4r = 0
    n_coh2d4 = 0
    n_c3d10h = 0

    n_cohax4 = 0
    n_cax3 = 0
    n_cax4r = 0
    n_c3d8r = 0

    nnodes = 0
    nelements = 0
    all_nodes = []
    for unused_part_name, part in model.parts.items():
        #unused_nids = part.nids - 1
        nodes = part.nodes

        nnodes += nodes.shape[0]
        if part.r2d2 is not None:
            n_r2d2 += part.r2d2.shape[0]

        # shells
        if part.cps3 is not None:
            n_cps3 += part.cps3.shape[0]
        if part.cpe3 is not None:
            n_cpe3 += part.cpe3.shape[0]
        if part.cpe4 is not None:
            n_cpe4 += part.cpe4.shape[0]
        if part.cpe4r is not None:
            n_cpe4r += part.cpe4r.shape[0]
        if part.coh2d4 is not None:
            n_coh2d4 += part.coh2d4.shape[0]

        if part.cohax4 is not None:
            n_cohax4 += part.cohax4.shape[0]
        if part.cax3 is not None:
            n_cax3 += part.cax3.shape[0]
        if part.cax4r is not None:
            n_cax4r += part.cax4r.shape[0]

        if part.c3d10h is not None:
            n_c3d10h += part.c3d10h.shape[0]
        if part.c3d8r is not None:
            n_c3d8r += part.c3d8r.shape[0]

        all_nodes.append(nodes)
    nelements += (
        n_r2d2 +
        n_cps3 + n_cpe3 + n_cpe4 + n_cpe4r +
        n_coh2d4 + n_c3d10h + n_cohax4 + n_cax3 + n_cax4r + n_c3d8r
    )
    return nnodes, all_nodes, nelements

def add_lines(grid, nids, eids_lines, nid_offset):
    """adds line elements to the vtkUnstructuredGrid"""
    nelements = 0
    if eids_lines is not None:
        nelements = eids_lines.shape[0]
        eids = eids_lines[:, 0]
        elem_nids = eids_lines[:, 1:]
        #inids = np.searchsorted(nids, elem_nids)
        node_ids = elem_nids + nid_offset
        for unused_eid, node_idsi in zip(eids, node_ids):
            elem = vtkLine()
            elem.GetPointIds().SetId(0, node_idsi[0])
            elem.GetPointIds().SetId(1, node_idsi[1])
            grid.InsertNextCell(3, elem.GetPointIds())
    return nelements


def add_tris(grid, nids, eids_tris, nid_offset):
    """adds tri elements to the vtkUnstructuredGrid"""
    nelements = 0
    if eids_tris is not None:
        nelements = eids_tris.shape[0]
        eids = eids_tris[:, 0]
        elem_nids = eids_tris[:, 1:]
        #inids = np.searchsorted(nids, elem_nids)
        node_ids = elem_nids + nid_offset
        for unused_eid, node_idsi in zip(eids, node_ids):
            elem = vtkTriangle()
            elem.GetPointIds().SetId(0, node_idsi[0])
            elem.GetPointIds().SetId(1, node_idsi[1])
            elem.GetPointIds().SetId(2, node_idsi[2])
            grid.InsertNextCell(5, elem.GetPointIds())
    return nelements


def add_quads(grid, nids, eids_quads, nid_offset):
    """adds quad elements to the vtkUnstructuredGrid"""
    nelements = 0
    if eids_quads is not None:
        nelements = eids_quads.shape[0]
        eids = eids_quads[:, 0]
        elem_nids = eids_quads[:, 1:]
        #inids = np.searchsorted(nids, elem_nids)
        #print(inids)
        node_ids = elem_nids + nid_offset
        #node_ids = inids # + nid_offset + 1
        for unused_eid, node_idsi in zip(eids, node_ids):
            elem = vtkQuad()
            elem.GetPointIds().SetId(0, node_idsi[0])
            elem.GetPointIds().SetId(1, node_idsi[1])
            elem.GetPointIds().SetId(2, node_idsi[2])
            elem.GetPointIds().SetId(3, node_idsi[3])
            grid.InsertNextCell(9, elem.GetPointIds())
    return nelements


def add_tetras(grid, nids, eids_tetras, nid_offset):
    """adds tet elements to the vtkUnstructuredGrid"""
    nelements = 0
    if eids_tetras is not None:
        nelements = eids_tetras.shape[0]
        eids = eids_tetras[:, 0]
        elem_nids = eids_tetras[:, 1:]
        #inids = np.searchsorted(nids, elem_nids)
        node_ids = elem_nids + nid_offset
        for unused_eid, node_idsi in zip(eids, node_ids):
            elem = vtkTetra()
            elem.GetPointIds().SetId(0, node_idsi[0])
            elem.GetPointIds().SetId(1, node_idsi[1])
            elem.GetPointIds().SetId(2, node_idsi[2])
            elem.GetPointIds().SetId(3, node_idsi[3])
            grid.InsertNextCell(10, elem.GetPointIds())
    return nelements


def add_hexas(grid, nids, eids_hexas, nid_offset):
    """adds hex elements to the vtkUnstructuredGrid"""
    nelements = 0
    if eids_hexas is not None:
        nelements = eids_hexas.shape[0]
        eids = eids_hexas[:, 0]
        elem_nids = eids_hexas[:, 1:]
        #inids = np.searchsorted(nids, elem_nids)
        node_ids = elem_nids + nid_offset
        for unused_eid, node_idsi in zip(eids, node_ids):
            elem = vtkHexahedron()
            elem.GetPointIds().SetId(0, node_idsi[0])
            elem.GetPointIds().SetId(1, node_idsi[1])
            elem.GetPointIds().SetId(2, node_idsi[2])
            elem.GetPointIds().SetId(3, node_idsi[3])
            elem.GetPointIds().SetId(4, node_idsi[4])
            elem.GetPointIds().SetId(5, node_idsi[5])
            elem.GetPointIds().SetId(6, node_idsi[6])
            elem.GetPointIds().SetId(7, node_idsi[7])
            grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
    return nelements
