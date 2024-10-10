"""Defines the GUI IO file for STL."""
from __future__ import annotations
from typing import Any, TYPE_CHECKING
import numpy as np

import vtkmodules
from pyNastran.gui.vtk_interface import vtkTriangle, vtkQuad

from pyNastran.converters.fluent.fluent import read_fluent
from pyNastran.gui.gui_objects.gui_result import GuiResult, NormalResult
from pyNastran.gui.utils.vtk.vtk_utils import (
    create_vtk_cells_of_constant_element_type, numpy_to_vtk_points)
from pyNastran.gui.utils.vtk.vectorized_geometry import (
    create_offset_arrays, build_vtk_geometry)
if TYPE_CHECKING:
    from pyNastran.gui.main_window import MainWindow
    from vtk import vtkUnstructuredGrid

class FluentIO:
    def __init__(self, gui: MainWindow):
        self.gui = gui

    def get_fluent_wildcard_geometry_results_functions(self):
        data = ('Fluent',
                'Fluent Field (*.cel)', self.load_fluent_geometry,
                None, None)
        return data

    def load_fluent_geometry(self, fld_filename: str,
                             name: str='main', plot: bool=True):
        model_name = name
        skip_reading = self.gui._remove_old_geometry(fld_filename)
        if skip_reading:
            return

        log = self.gui.log
        model = read_fluent(fld_filename, log=log, debug=False)
        #self.model_type = model.model_type
        nodes = model.xyz

        #self.node = node
        #self.xyz = xyz
        #self.element_id = element_id
        #self.results = results
        #self.quads = quads
        #self.tris = tris

        # support multiple results
        titles = model.titles
        results = model.results

        if 0:
            element_ids = model.element_ids
            region = model.region
            elements_list = model.elements_list
            nelement = len(element_ids)
            is_list = True
        else:
            is_list = False
            tris = model.tris
            quads = model.quads
            element_ids = model.element_id #np.arange(1, nelement+1)
            assert np.array_equal(element_ids, np.unique(element_ids))

            # we reordered the tris/quads to be continuous to make them easier to add
            itri = np.searchsorted(element_ids, tris[:, 0])
            iquad = np.searchsorted(element_ids, quads[:, 0])

            tri_results = results[itri, :]
            quad_results = results[iquad, :]

            region = np.hstack([tris[:, 1], quads[:, 1],])
            results = np.vstack([tri_results, quad_results])
            ntri = len(tris)
            nquad = len(quads)
            # else:
            #     ntri, ntri_col = tris.shape
            #     nquad, nquad_col = quads.shape
            #     ntri = 0
            #     region = quads[:, 1]
            #     tris = np.zeros((ntri, ntri_col))
            #     results = quad_results
            nelement = ntri + nquad
            log.debug(f'results.shape = {results.shape}')

        node_ids = model.node_id
        nnodes = len(nodes)
        #log.debug(f'nelement={nelement} nregion={len(region)}')
        #log.debug('nnodes={nnodes} nnodes2={len(node_ids)}')
        #log.debug(f'nquad={nquad} ntri={ntri}')
        assert len(element_ids) == len(region)

        self.gui.nnodes = nnodes
        self.gui.nelements = nelement

        self.gui.log.info('nnodes=%s nelements=%s' % (self.gui.nnodes, self.gui.nelements))
        grid = self.gui.grid
        grid.Allocate(self.gui.nelements, 1000)

        points = numpy_to_vtk_points(nodes)
        if is_list:
            _create_elements_list(grid, node_ids, elements_list)
        else:
            _create_elements(grid, node_ids, model.tris, model.quads)

        assert len(element_ids) == len(region)
        log.info(f'created vtk points')
        self.gui.nid_map = {}
        #elem.SetNumberOfPoints(nnodes)

        assert nodes is not None
        xmax, ymax, zmax = nodes.max(axis=0)
        xmin, ymin, zmin = nodes.min(axis=0)
        self.gui.log.info('xmax=%s xmin=%s' % (xmax, xmin))
        self.gui.log.info('ymax=%s ymin=%s' % (ymax, ymin))
        self.gui.log.info('zmax=%s zmin=%s' % (zmax, zmin))
        dim_max = max(xmax-xmin, ymax-ymin, zmax-zmin)
        self.gui.create_global_axes(dim_max)

        log.info(f'created vtk elements')

        grid.SetPoints(points)
        grid.Modified()

        # loadSTLResults - regions/loads
        self.gui.scalar_bar_actor.VisibilityOff()
        self.gui.scalar_bar_actor.Modified()

        cases = {}
        self.gui.isubcase_name_map = {}
        ID = 1
        if is_list:
            form, cases = self._fill_fluent_case_list(
                cases, ID, node_ids, element_ids,
                region, results, titles)
        else:
            form, cases = self._fill_fluent_case(
                cases, ID, node_ids, element_ids,
                #tris, quads,
                region, results, titles)
        self.gui.node_ids = node_ids
        self.gui.element_ids = element_ids
        log.debug(f'running _finish_results_io2')
        self.gui._finish_results_io2(model_name, form, cases)
        log.info(f'finished')

    def _fill_fluent_case(self,
                          cases: dict[int, Any],
                          ID: int,
                          node_ids: np.ndarray,
                          element_ids: np.ndarray,
                          region: np.ndarray,
                          #tris: np.ndarray,
                          #quads: np.ndarray,
                          results: np.ndarray,
                          titles: np.ndarray) -> None:
        """adds the sidebar results"""
        self.gui.isubcase_name_map[ID] = ('Fluent', '')
        # reorg the ids
        #element_ids = np.unique(np.hstack([tris[:, 0], quads[:, 0]]))
        #colormap = 'jet'
        icase = 0
        itime = 0
        nid_res = GuiResult(ID, header='NodeID', title='NodeID',
                            location='node', scalar=node_ids)
        eid_res = GuiResult(ID, header='ElementID', title='ElementID',
                            location='centroid', scalar=element_ids)
        region_res = GuiResult(ID, header='Region', title='Region',
                               location='centroid', scalar=region)

        assert len(element_ids) == len(region), f'neids={len(element_ids)} nregion={len(region)}'
        cases[icase] = (nid_res, (itime, 'NodeID'))
        cases[icase + 1] = (eid_res, (itime, 'ElementID'))
        cases[icase + 2] = (region_res, (itime, 'Region'))

        form = [
            ('NodeID', icase, []),
            ('ElementID', icase + 1, []),
            ('Region', icase + 2, []),
        ]
        icase += 3

        for i, title in enumerate(titles[1:]):
            result = results[:, i]
            assert len(element_ids) == len(result)
            pressure_res = GuiResult(ID, header=title, title=title,
                                     location='centroid', scalar=result)
            cases[icase] = (pressure_res, (itime, title))
            formi = (title, icase, [])
            form.append(formi)
            icase += 1

        return form, cases

    def _fill_fluent_case_list(self,
                               cases: dict[int, Any],
                               ID: int,
                               node_ids: np.ndarray,
                               element_ids: np.ndarray,
                               region: np.ndarray,
                               results: np.ndarray,
                               titles: np.ndarray) -> None:
        """adds the sidebar results"""
        self.gui.isubcase_name_map[ID] = ('Fluent', '')
        #colormap = 'jet'
        icase = 0
        itime = 0
        nid_res = GuiResult(ID, header='NodeID', title='NodeID',
                            location='node', scalar=node_ids)
        eid_res = GuiResult(ID, header='ElementID', title='ElementID',
                            location='centroid', scalar=element_ids)
        region_res = GuiResult(ID, header='Region', title='Region',
                               location='centroid', scalar=region)
        cases[icase] = (nid_res, (itime, 'NodeID'))
        cases[icase + 1] = (eid_res, (itime, 'ElementID'))
        cases[icase + 2] = (region_res, (itime, 'Region'))

        form = [
            ('NodeID', icase, []),
            ('ElementID', icase + 1, []),
            ('Region', icase + 2, []),
        ]
        icase += 3

        for i, title in enumerate(titles[1:]):
            result = results[:, i]
            pressure_res = GuiResult(ID, header=title, title=title,
                                     location='centroid', scalar=result)
            cases[icase] = (pressure_res, (itime, title))
            formi = (title, icase, [])
            form.append(formi)
            icase += 1

        return form, cases

def _create_elements_list(grid: vtkUnstructuredGrid,
                          node_ids: np.ndarray,
                          elements_list: list[list[int]]):
    assert np.array_equal(node_ids, np.unique(node_ids))
    nid_to_index = {nid : i for i, nid in enumerate(node_ids)}

    for facei in elements_list:
        face = np.array(facei, dtype='int32')
        if len(face) == 3:
            elem = vtkTriangle()
            epoints = elem.GetPointIds()
            epoints.SetId(0, nid_to_index[face[0]])
            epoints.SetId(1, nid_to_index[face[1]])
            epoints.SetId(2, nid_to_index[face[2]])
            grid.InsertNextCell(5, epoints)
        elif len(face) == 4:
            elem = vtkQuad()
            epoints = elem.GetPointIds()
            epoints.SetId(0, nid_to_index[face[0]])
            epoints.SetId(1, nid_to_index[face[1]])
            epoints.SetId(2, nid_to_index[face[2]])
            epoints.SetId(3, nid_to_index[face[3]])
            grid.InsertNextCell(9, epoints)
        else:  # pragma: no cover
            raise RuntimeError(face)

def _create_elements(ugrid: vtkUnstructuredGrid,
                     node_ids: np.ndarray,
                     tris: np.ndarray,
                     quads: np.ndarray) -> None:
    #element_ids_list = []
    #region_list = []
    cell_type_list = []
    cell_offset_list = []

    # cell_offset0, n_nodesi, cell_typei, cell_offseti = _create_solid_vtk_arrays(
    #     element, grid_id, cell_offset0)
    # n_nodes_.append(n_nodesi.ravel())
    # element_ids.append(element.element_id)
    # property_ids.append(element.property_id)
    # cell_type_.append(cell_typei)
    # cell_offset_.append(cell_offseti)
    # del n_nodesi, cell_typei, cell_offseti

    #quad_nodes = quads[:, 2:] - 1
    #tris_nodes = tris[:, 2:] - 1
    #elements_list = []
    nquad = len(quads)
    ntri = len(tris)
    nelement_total = nquad + ntri
    cell_offset0 = 0
    n_nodes_list = []
    if nquad:
        cell_type = 9
        dnode = 4
        cell_offset0, n_nodesi, cell_typei, cell_offseti = create_offset_arrays(
            node_ids, quads[:, 2:],
            nquad, cell_type, cell_offset0, dnode)

        #cell_typei = np.ones(nquad) * etype
        #cell_offseti = np.ones(nquad) * offset
        #elements_list.append()
        #elem.GetCellType() = 9  # vtkQuad

        #element_ids.append(element.element_id)
        #property_ids.append(element.property_id)
        n_nodes_list.append(n_nodesi.ravel())
        cell_type_list.append(cell_typei)
        cell_offset_list.append(cell_offseti)

        if 0:
            quad_nodes = np.searchsorted(node_ids, quads[:, 2:])
            for face in quad_nodes:
                assert len(face) == 4, face
                elem = vtkQuad()
                epoints = elem.GetPointIds()
                epoints.SetId(0, face[0])
                epoints.SetId(1, face[1])
                epoints.SetId(2, face[2])
                epoints.SetId(3, face[3])
                ugrid.InsertNextCell(9, epoints)

    if ntri:
        cell_type = 5
        dnode = 3
        #cell_typei = np.ones(ntri) * cell_type
        #cell_offseti = np.ones(ntri) * offset
        #cell_type_list.append(cell_typei)
        #cell_offset_list.append(cell_offseti)

        cell_offset0, n_nodesi, cell_typei, cell_offseti = create_offset_arrays(
            node_ids, tris[:, 2:],
            ntri, cell_type, cell_offset0, dnode)
        n_nodes_list.append(n_nodesi.ravel())
        cell_type_list.append(cell_typei)
        cell_offset_list.append(cell_offseti)

        #elem.GetCellType() = 5  # vtkTriangle
        if 0:
            tris_nodes = np.searchsorted(node_ids, tris[:, 2:])
            for face in tris_nodes:
                assert len(face) == 3, face
                elem = vtkTriangle()
                epoints = elem.GetPointIds()
                epoints.SetId(0, face[0])
                epoints.SetId(1, face[1])
                epoints.SetId(2, face[2])
                ugrid.InsertNextCell(5, epoints)

    n_nodes = np.hstack(n_nodes_list)
    cell_type = np.hstack(cell_type_list)
    cell_offset = np.hstack(cell_offset_list)
    build_vtk_geometry(
        nelement_total, ugrid,
        n_nodes, cell_type, cell_offset)

    # TODO: make this faster...
    #etype = 5  # vtkTriangle().GetCellType()
    #create_vtk_cells_of_constant_element_type(grid, tris, etype)
