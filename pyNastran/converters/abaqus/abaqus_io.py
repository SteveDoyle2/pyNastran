"""Defines how the GUI reads Abaqus files"""
from collections import OrderedDict
from typing import Tuple, Any, TYPE_CHECKING

import numpy as np

import vtk
from vtk import vtkLine, vtkTriangle, vtkQuad, vtkTetra, vtkHexahedron
from pyNastran.gui.utils.vtk.vtk_utils import numpy_to_vtk

from pyNastran.gui.gui_objects.gui_result import GuiResult, NormalResult
#from pyNastran.gui.qt_files.result import Result
from pyNastran.converters.abaqus.abaqus import Abaqus
if TYPE_CHECKING:  # pragma: no cover
    from vtk import vtkUnstructuredGrid


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

        assert len(all_nodes) > 0, len(all_nodes)
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
            array_type=data_type,
        )
        points.SetData(points_array)

        nid_offset = -1
        nids = []
        for unused_part_name, part in model.parts.items():
            self.gui.log.info('part_name = %r' % unused_part_name)
            nnodesi = part.nodes.shape[0]
            nidsi = part.nids
            nids.append(nidsi)

            elements = part
            add_lines(grid, nidsi, elements.r2d2, nid_offset)

            add_tris(grid, nidsi, elements.cps3, nid_offset)
            add_tris(grid, nidsi, elements.cpe3, nid_offset)

            add_quads(grid, nidsi, elements.cpe4, nid_offset)
            add_quads(grid, nidsi, elements.cpe4r, nid_offset)

            add_quads(grid, nidsi, elements.cps4, nid_offset)
            add_quads(grid, nidsi, elements.cps4r, nid_offset)

            add_quads(grid, nidsi, elements.coh2d4, nid_offset)
            add_quads(grid, nidsi, elements.cohax4, nid_offset)

            add_tris(grid, nidsi, elements.cax3, nid_offset)
            #add_quads(grid, nidsi, elements.cax4, nid_offset)
            add_quads(grid, nidsi, elements.cax4r, nid_offset)

            # solids
            add_tetras(grid, nidsi, elements.c3d10h, nid_offset)
            add_hexas(grid, nidsi, elements.c3d8r, nid_offset)

            nid_offset += nnodesi
        nids = np.hstack(nids)
        grid.SetPoints(points)
        grid.Modified()

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

    def clear_abaqus(self) -> None:
        """does nothing"""
        pass

    def load_abaqus_results(self, abaqus_filename: str):
        """does nothing"""
        raise NotImplementedError()

    def _fill_abaqus_case(self, cases, ID: int, node_ids, nodes, nelements: int,
                          unused_model: Abaqus) -> Tuple[Any, Any, int, np.ndarray, np.ndarray]:
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


def get_nodes_nnodes_nelements(model, stop_for_no_elements=True):
    """helper method"""
    nnodes = 0
    nelements = 0
    all_nodes = []
    for unused_part_name, part in model.parts.items():
        #unused_nids = part.nids - 1
        nodes = part.nodes

        nnodes += nodes.shape[0]
        nelements += part.nelements

        all_nodes.append(nodes)
    if nelements == 0 and stop_for_no_elements:
        raise RuntimeError('nelements=0')
    return nnodes, all_nodes, nelements

def add_lines(grid: vtk.vtkUnstructuredGrid,
              nids, eids_lines: np.ndarray, nid_offset: int) -> int:
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
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, node_idsi[0])
            point_ids.SetId(1, node_idsi[1])
            grid.InsertNextCell(3, point_ids)
    return nelements


def add_tris(grid: vtk.vtkUnstructuredGrid,
             nids, eids_tris: np.ndarray, nid_offset: int) -> int:
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
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, node_idsi[0])
            point_ids.SetId(1, node_idsi[1])
            point_ids.SetId(2, node_idsi[2])
            grid.InsertNextCell(5, point_ids)
    return nelements


def add_quads(grid: vtk.vtkUnstructuredGrid,
              nids, eids_quads: np.ndarray, nid_offset: int) -> int:
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
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, node_idsi[0])
            point_ids.SetId(1, node_idsi[1])
            point_ids.SetId(2, node_idsi[2])
            point_ids.SetId(3, node_idsi[3])
            grid.InsertNextCell(9, point_ids)
    return nelements


def add_tetras(grid: vtk.vtkUnstructuredGrid,
               nids, eids_tetras: np.ndarray, nid_offset: int) -> int:
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
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, node_idsi[0])
            point_ids.SetId(1, node_idsi[1])
            point_ids.SetId(2, node_idsi[2])
            point_ids.SetId(3, node_idsi[3])
            grid.InsertNextCell(10, point_ids)
    return nelements


def add_hexas(grid: vtk.vtkUnstructuredGrid,
              nids, eids_hexas: np.ndarray, nid_offset: int) -> int:
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
