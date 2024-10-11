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

        #self.node = node
        nodes = model.xyz
        #self.element_id = element_id
        tris = model.tris
        quads = model.quads

        # support multiple results
        titles = model.titles
        results = model.results

        element_id = model.element_id #np.arange(1, nelement+1)
        #assert np.array_equal(element_id, np.unique(element_id))

        # we reordered the tris/quads to be continuous to make them easier to add
        iquad = np.searchsorted(element_id, quads[:, 0])
        itri = np.searchsorted(element_id, tris[:, 0])

        quad_results = results[iquad, :]
        tri_results = results[itri, :]

        region = np.hstack([quads[:, 1], tris[:, 1]])
        results = np.vstack([quad_results, tri_results])
        nquad = len(quads)
        ntri = len(tris)
        nelement = nquad + ntri
        log.debug(f'results.shape = {results.shape}')

        node_id = model.node_id
        nnodes = len(nodes)
        assert len(element_id) == len(region)

        self.gui.nnodes = nnodes
        self.gui.nelements = nelement

        self.gui.log.info('nnodes=%s nelements=%s' % (self.gui.nnodes, self.gui.nelements))
        grid = self.gui.grid
        grid.Allocate(self.gui.nelements, 1000)

        points = numpy_to_vtk_points(nodes)
        _create_elements(grid, node_id, model.quads, model.tris)

        assert len(element_id) == len(region)
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

        self.gui.isubcase_name_map[ID] = ('Fluent', '')
        form, cases = _fill_fluent_case(
            cases, ID, node_id, element_id,
            region, results, titles)

        self.gui.node_ids = node_id
        self.gui.element_ids = element_id
        log.debug(f'running _finish_results_io2')
        self.gui._finish_results_io2(model_name, form, cases)
        log.info(f'finished')

def _create_elements(ugrid: vtkUnstructuredGrid,
                     node_ids: np.ndarray,
                     tris: np.ndarray,
                     quads: np.ndarray) -> None:
    cell_type_list = []
    cell_offset_list = []

    nquad = len(quads)
    ntri = len(tris)
    nelement_total = nquad + ntri
    cell_offset0 = 0
    n_nodes_list = []
    if nquad:
        #elem.GetCellType() = 9  # vtkQuad
        cell_type = 9
        dnode = 4
        cell_offset0, n_nodesi, cell_typei, cell_offseti = create_offset_arrays(
            node_ids, quads[:, 2:],
            nquad, cell_type, cell_offset0, dnode)
        n_nodes_list.append(n_nodesi.ravel())
        cell_type_list.append(cell_typei)
        cell_offset_list.append(cell_offseti)

    if ntri:
        #elem.GetCellType() = 5  # vtkTriangle
        cell_type = 5
        dnode = 3
        cell_offset0, n_nodesi, cell_typei, cell_offseti = create_offset_arrays(
            node_ids, tris[:, 2:],
            ntri, cell_type, cell_offset0, dnode)
        n_nodes_list.append(n_nodesi.ravel())
        cell_type_list.append(cell_typei)
        cell_offset_list.append(cell_offseti)

    n_nodes = np.hstack(n_nodes_list)
    cell_type = np.hstack(cell_type_list)
    cell_offset = np.hstack(cell_offset_list)
    build_vtk_geometry(
        nelement_total, ugrid,
        n_nodes, cell_type, cell_offset)
    return

def _fill_fluent_case(cases: dict[int, Any],
                      ID: int,
                      node_ids: np.ndarray,
                      element_ids: np.ndarray,
                      region: np.ndarray,
                      results: np.ndarray,
                      titles: np.ndarray) -> None:
    """adds the sidebar results"""
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
