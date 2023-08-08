"""Defines how the GUI reads Abaqus files"""
from typing import Any

import numpy as np

from pyNastran.gui.vtk_common_core import vtkPoints, VTK_FLOAT
from pyNastran.gui.vtk_interface import (
    vtkLine, vtkTriangle, vtkQuad, vtkTetra, vtkHexahedron,
    vtkUnstructuredGrid)
from pyNastran.gui.utils.vtk.vtk_utils import numpy_to_vtk

from pyNastran.gui.gui_objects.gui_result import GuiResult, NormalResult
#from pyNastran.gui.qt_files.result import Result
from pyNastran.converters.abaqus.abaqus import Abaqus, get_nodes_nnodes_nelements


class AbaqusIO:
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

    def load_abaqus_geometry(self, abaqus_filename: str,
                             name: str='main', plot: bool=True) -> None:
        """loads abaqus input files into the gui"""
        model_name = name
        skip_reading = self.gui._remove_old_geometry(abaqus_filename)
        if skip_reading:
            return

        self.gui.eid_maps[name] = {}
        self.gui.nid_maps[name] = {}
        model = Abaqus(log=self.gui.log, debug=False)
        self.gui.model_type = 'abaqus'
        model.read_abaqus_inp(abaqus_filename)

        self.gui.nid_map = {}
        nnodes, nids, nodes, nelements = get_nodes_nnodes_nelements(model)
        self.gui.log.info('nnodes=%s nelements=%s' % (nnodes, nelements))
        assert nelements > 0, nelements
        #nodes = model.nodes
        #elements = model.elements


        self.gui.nnodes = nnodes
        self.gui.nelements = nelements

        grid = self.gui.grid
        grid.Allocate(self.gui.nelements, 1000)

        assert len(nodes) > 0, len(nodes)

        mmax = np.amax(nodes, axis=0)
        mmin = np.amin(nodes, axis=0)
        dim_max = (mmax - mmin).max()
        self.gui.create_global_axes(dim_max)

        points = vtkPoints()
        points.SetNumberOfPoints(self.gui.nnodes)

        data_type = VTK_FLOAT
        points_array = numpy_to_vtk(
            num_array=nodes,
            deep=True,
            array_type=data_type,
        )
        points.SetData(points_array)

        nid_offset = -1
        nids = []
        if model.nodes is not None:
            nid_offset = add_part(grid, model, nids, nid_offset)

        for part_name, part in model.parts.items():
            self.gui.log.info(f'part_name = {part_name!r}')
            nid_offset = add_part(grid, part, nids, nid_offset)
        assert len(nids) > 0, nids
        nids = np.hstack(nids)
        grid.SetPoints(points)
        grid.Modified()

        self.gui.scalar_bar_actor.VisibilityOn()
        self.gui.scalar_bar_actor.Modified()

        note = ''
        self.gui.isubcase_name_map = {1: ['Abaqus%s' % note, '']}
        #form = []
        cases = {}
        ID = 1
        form, cases, unused_icase, node_ids, element_ids = self._fill_abaqus_case(
            cases, ID, nids, nodes, nelements, model)
        #self._fill_cart3d_results(cases, form, icase, ID, model)

        self.gui.node_ids = node_ids
        self.gui.element_ids = element_ids
        self.gui._finish_results_io2(model_name, form, cases)

    def clear_abaqus(self) -> None:
        """does nothing"""
        pass

    def load_abaqus_results(self, abaqus_filename: str):
        """does nothing"""
        raise NotImplementedError()

    def _fill_abaqus_case(self, cases, ID: int, node_ids, nodes, nelements: int,
                          unused_model: Abaqus) -> tuple[Any, Any, int, np.ndarray, np.ndarray]:
        """creates the result objects for abaqus"""
        #return [], {}, 0
        #nelements = elements.shape[0]
        #nnodes = nodes.shape[0]

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

def add_part(grid: vtkUnstructuredGrid,
             part,
             nids: list[np.ndarray],
             nid_offset: int) -> int:
    nnodesi = part.nodes.shape[0]
    nidsi = part.nids
    nids.append(nidsi)
    elements = part.elements
    add_lines(grid, nidsi, elements.r2d2_eids, elements.r2d2, nid_offset)
    add_lines(grid, nidsi, elements.b31h_eids, elements.b31h, nid_offset)

    add_tris(grid, nidsi, elements.cps3_eids, elements.cps3, nid_offset)
    add_tris(grid, nidsi, elements.cpe3_eids, elements.cpe3, nid_offset)

    add_quads(grid, nidsi, elements.cpe4_eids, elements.cpe4, nid_offset)
    add_quads(grid, nidsi, elements.cpe4r_eids, elements.cpe4r, nid_offset)

    add_quads(grid, nidsi, elements.cps4_eids, elements.cps4, nid_offset)
    add_quads(grid, nidsi, elements.cps4r_eids, elements.cps4r, nid_offset)

    add_quads(grid, nidsi, elements.coh2d4_eids, elements.coh2d4, nid_offset)
    add_quads(grid, nidsi, elements.cohax4_eids, elements.cohax4, nid_offset)

    add_tris(grid, nidsi, elements.cax3_eids, elements.cax3, nid_offset)
    #add_quads(grid, nidsi, elements.cax4_eids, elements.cax4, nid_offset)
    add_quads(grid, nidsi, elements.cax4r_eids, elements.cax4r, nid_offset)

    # solids
    add_tetras(grid, nidsi, elements.c3d10h_eids, elements.c3d10h, nid_offset)
    add_hexas(grid, nidsi, elements.c3d8r_eids, elements.c3d8r, nid_offset)

    nid_offset += nnodesi
    return nid_offset

def add_lines(grid: vtkUnstructuredGrid,
              nids, eids, elem_nids: np.ndarray, nid_offset: int) -> int:
    """adds line elements to the vtkUnstructuredGrid"""
    nelements = 0
    if elem_nids is not None:
        nelements = len(eids)
        #eids = eids_lines[:, 0]
        #elem_nids = eids_lines[:, 1:]
        #inids = np.searchsorted(nids, elem_nids)
        node_ids = elem_nids + nid_offset
        for node_idsi in node_ids:
            elem = vtkLine()
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, node_idsi[0])
            point_ids.SetId(1, node_idsi[1])
            grid.InsertNextCell(3, point_ids)
    return nelements


def add_tris(grid: vtkUnstructuredGrid,
             nids, eids, elem_nids: np.ndarray, nid_offset: int) -> int:
    """adds tri elements to the vtkUnstructuredGrid"""
    nelements = 0
    if eids is not None:
        nelements = len(elem_nids)
        #inids = np.searchsorted(nids, elem_nids)
        node_ids = elem_nids + nid_offset
        for node_idsi in node_ids:
            elem = vtkTriangle()
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, node_idsi[0])
            point_ids.SetId(1, node_idsi[1])
            point_ids.SetId(2, node_idsi[2])
            grid.InsertNextCell(5, point_ids)
    return nelements


def add_quads(grid: vtkUnstructuredGrid,
              nids, eids: np.ndarray, elem_nids: np.ndarray, nid_offset: int) -> int:
    """adds quad elements to the vtkUnstructuredGrid"""
    nelements = 0
    if elem_nids is not None:
        nelements = len(eids)
        #inids = np.searchsorted(nids, elem_nids)
        #print(inids)
        node_ids = elem_nids + nid_offset
        #node_ids = inids # + nid_offset + 1
        for unused_eid, node_idsi in zip(eids, node_ids):
            elem = vtkQuad()
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, node_idsi[0])
            point_ids.SetId(1, node_idsi[1])
            point_ids.SetId(2, node_idsi[2])
            point_ids.SetId(3, node_idsi[3])
            grid.InsertNextCell(9, point_ids)
    return nelements


def add_tetras(grid: vtkUnstructuredGrid,
               nids, eids: np.ndarray, elem_nids: np.ndarray, nid_offset: int) -> int:
    """adds tet elements to the vtkUnstructuredGrid"""
    nelements = 0
    if elem_nids is not None:
        nelements = len(elem_nids)
        #inids = np.searchsorted(nids, elem_nids)
        node_ids = elem_nids + nid_offset
        for node_idsi in node_ids:
            elem = vtkTetra()
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, node_idsi[0])
            point_ids.SetId(1, node_idsi[1])
            point_ids.SetId(2, node_idsi[2])
            point_ids.SetId(3, node_idsi[3])
            grid.InsertNextCell(10, point_ids)
    return nelements


def add_hexas(grid: vtkUnstructuredGrid,
              nids, eids: np.ndarray, elem_nids: np.ndarray, nid_offset: int) -> int:
    """adds hex elements to the vtkUnstructuredGrid"""
    nelements = 0
    if elem_nids is not None:
        nelements = len(elem_nids)
        #inids = np.searchsorted(nids, elem_nids)
        node_ids = elem_nids + nid_offset
        for node_idsi in node_ids:
            elem = vtkHexahedron()
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, node_idsi[0])
            point_ids.SetId(1, node_idsi[1])
            point_ids.SetId(2, node_idsi[2])
            point_ids.SetId(3, node_idsi[3])
            point_ids.SetId(4, node_idsi[4])
            point_ids.SetId(5, node_idsi[5])
            point_ids.SetId(6, node_idsi[6])
            point_ids.SetId(7, node_idsi[7])
            grid.InsertNextCell(elem.GetCellType(), point_ids)
    return nelements
