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
        model = read_fluent(
            fld_filename, #auto_read_write_h5=False,
            log=log, debug=False)
        #self.model_type = model.model_type

        node_id = model.node_id
        nodes = model.xyz
        nnodes = len(nodes)

        # support multiple results
        titles = model.titles
        results = model.results
        if 0:  # pragma: no cover
            element_id = model.element_ids
            region = model.region
            elements_list = model.elements_list
            nelement = len(element_id)
            is_list = True
        else:
            is_list = False
            tris = model.tris
            quads = model.quads

            tri_regions = tris[:, 1]
            quad_regions = quads[:, 1]

            # we reordered the tris/quads to be continuous to make them easier to add
            iquad = np.searchsorted(model.element_id, quads[:, 0])
            itri = np.searchsorted(model.element_id, tris[:, 0])

            region_split = False
            if region_split:
                #regions_to_remove = [3]
                regions_to_remove = []
                regions_to_include = [7, 4]
                is_remove = (len(regions_to_remove) == 0)
                is_include = (len(regions_to_include) == 0)
                assert (is_remove and not is_include) or (not is_remove and is_include)
                if regions_to_remove:
                    itri_regions = np.logical_and.reduce([(tri_regions != regioni) for regioni in regions_to_remove])
                    iquad_regions = np.logical_and.reduce([(quad_regions != regioni) for regioni in regions_to_remove])
                else:
                    itri_regions = np.logical_or.reduce([(tri_regions == regioni) for regioni in regions_to_include])
                    iquad_regions = np.logical_or.reduce([(quad_regions == regioni) for regioni in regions_to_include])

                quad_results = results[iquad, :][iquad_regions, :]
                tri_results = results[itri, :][itri_regions, :]

                tris = tris[itri_regions, :]
                quads = quads[iquad_regions, :]
            else:
                quad_results = results[iquad, :]
                tri_results = results[itri, :]

        self.gui.nnodes = nnodes
        self.gui.nelements = nelement

        self.gui.log.info('nnodes=%s nelements=%s' % (self.gui.nnodes, self.gui.nelements))
        ugrid = self.gui.grid
        ugrid.Allocate(self.gui.nelements, 1000)

        assert nodes is not None
        points = numpy_to_vtk_points(nodes)
        ugrid.SetPoints(points)
        log.info(f'created vtk points')

        xmax, ymax, zmax = nodes.max(axis=0)
        xmin, ymin, zmin = nodes.min(axis=0)
        log.info('xmax=%s xmin=%s' % (xmax, xmin))
        log.info('ymax=%s ymin=%s' % (ymax, ymin))
        log.info('zmax=%s zmin=%s' % (zmax, zmin))
        dim_max = max(xmax-xmin, ymax-ymin, zmax-zmin)

        assert len(node_id) == len(np.unique(node_id))
        if is_list:
            _create_elements_list(ugrid, node_id, elements_list)
        else:
            _create_elements(ugrid, node_id, tris, quads)
        log.info(f'created vtk elements')

        self.gui.nid_map = {}
        self.gui.create_global_axes(dim_max)

        ugrid.Modified()

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

def _create_elements_list(ugrid: vtkUnstructuredGrid,
                          node_id: np.ndarray,
                          elements_list: list[list[int]]):  # pragma: no cover5
    assert np.array_equal(node_id, np.unique(node_id))
    assert node_id.min() >= 0, node_id.min()
    nid_to_index = {nid : i for i, nid in enumerate(node_id)}
    #print(min(node_id), max(node_id))

    for facei in elements_list:
        face = np.array(facei, dtype='int32')  # cast str to int
        if len(face) == 3:
            elem = vtkTriangle()
            epoints = elem.GetPointIds()
            epoints.SetId(0, nid_to_index[face[0]])
            epoints.SetId(1, nid_to_index[face[1]])
            epoints.SetId(2, nid_to_index[face[2]])
            ugrid.InsertNextCell(5, epoints)
        elif len(face) == 4:
            elem = vtkQuad()
            epoints = elem.GetPointIds()
            epoints.SetId(0, nid_to_index[face[0]])
            epoints.SetId(1, nid_to_index[face[1]])
            epoints.SetId(2, nid_to_index[face[2]])
            epoints.SetId(3, nid_to_index[face[3]])
            ugrid.InsertNextCell(9, epoints)
        else:  # pragma: no cover
            raise RuntimeError(face)

def _create_elements(ugrid: vtkUnstructuredGrid,
                     node_id: np.ndarray,
                     tris: np.ndarray,
                     quads: np.ndarray) -> None:
    cell_type_list = []
    cell_offset_list = []

    nquad = len(quads)
    ntri = len(tris)
    nelement_total = nquad + ntri
    cell_offset0 = 0
    n_nodes_list = []

    assert tris.shape[1] == 5, tris.shape
    assert quads.shape[1] == 6, quads.shape
    tri_nodes = tris[:, 2:]
    quad_nodes = quads[:, 2:]
    assert tri_nodes.shape[1] == 3, tri_nodes.shape
    assert quad_nodes.shape[1] == 4, quad_nodes.shape
    assert node_id.min() >= 1, node_id.min()
    all_nodes = np.unique(np.hstack([tri_nodes.ravel(), quad_nodes.ravel()]))
    assert all_nodes.min() >= 1, all_nodes.min()
    missing_nodes = np.setdiff1d(all_nodes, node_id)
    assert len(missing_nodes) == 0, missing_nodes
    if nquad:
        #elem.GetCellType() = 9  # vtkQuad
        cell_type = 9
        dnode = 4
        cell_offset0, n_nodesi, cell_typei, cell_offseti = create_offset_arrays(
            node_id, quad_nodes,
            nquad, cell_type, cell_offset0, dnode)
        n_nodes_list.append(n_nodesi.ravel())
        cell_type_list.append(cell_typei)
        cell_offset_list.append(cell_offseti)

    if ntri:
        #elem.GetCellType() = 5  # vtkTriangle
        cell_type = 5
        dnode = 3
        cell_offset0, n_nodesi, cell_typei, cell_offseti = create_offset_arrays(
            node_id, tri_nodes,
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
