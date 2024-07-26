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
if TYPE_CHECKING:
    from pyNastran.gui.main_window import MainWindow

class FluentIO:
    def __init__(self, gui: MainWindow):
        self.gui = gui

    def get_fluent_wildcard_geometry_results_functions(self):
        data = ('Fluent',
                'Fluent Field (*.cel)', self.load_fluent_geometry,
                None, None)
        return data

    def load_fluent_geometry(self, fld_filename: str, name: str='main', plot: bool=True):
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
        results = model.results
        log.info(f'results.shape = {results.shape}')
        #pressure = results[:, -1]
        results = results.ravel()

        if 1:
            element_ids = model.element_ids
            region = model.region
            elements_list = model.elements_list
            pressure = results
            nelement = len(element_ids)
            is_list = True
        else:
            is_list = False
            tris = model.tris
            quads = model.quads
            element_ids = model.element_id #np.arange(1, nelement+1)

            # we reordered the tris/quads to be continuous to make them easier to add
            itri = np.searchsorted(element_ids, tris[:, 0])
            iquad = np.searchsorted(element_ids, quads[:, 0])

            tri_pressure = results[itri]
            quad_pressure = results[iquad]
            region = np.hstack([tris[:, 1], quads[:, 1],])
            pressure = np.hstack([tri_pressure, quad_pressure])
        
            ntri = len(tris)
            nquad = len(quads)
            nelement = ntri + nquad

        nnodes = len(nodes)
        node_ids = model.node_id

        self.gui.nnodes = nnodes
        self.gui.nelements = nelement

        self.gui.log.info('nnodes=%s nelements=%s' % (self.gui.nnodes, self.gui.nelements))
        grid = self.gui.grid
        grid.Allocate(self.gui.nelements, 1000)

        points = numpy_to_vtk_points(nodes)
        if is_list:
            _create_elements_list(grid, elements_list)
        else:
            _create_elements(grid, model.tris[:, 2:]-1, model.quads[:, 2:]-1)

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

        form, cases = self._fill_fluent_case(
            cases, ID, node_ids, element_ids, region, pressure)
        self.gui.node_ids = node_ids
        self.gui.element_ids = element_ids
        log.info(f'running _finish_results_io2')
        self.gui._finish_results_io2(model_name, form, cases)
        log.info(f'finished')

    def _fill_fluent_case(self, cases: dict[int, Any], ID: int,
                          node_ids: np.ndarray, 
                          element_ids: np.ndarray,
                          region: np.ndarray,
                          pressure: np.ndarray) -> None:
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
        pressure_res = GuiResult(ID, header='Pressure', title='Pressure',
                                 location='centroid', scalar=pressure)
        cases[icase] = (nid_res, (itime, 'NodeID'))
        cases[icase + 1] = (eid_res, (itime, 'ElementID'))
        cases[icase + 2] = (region_res, (itime, 'Region'))
        cases[icase + 3] = (pressure_res, (itime, 'Pressure'))

        form = [
            ('NodeID', icase, []),
            ('ElementID', icase + 1, []),
            ('Region', icase + 2, []),
            ('Pressure', icase + 3, []),
        ]
        return form, cases

def _create_elements_list(grid, elements_list: list[list[int]]):
        for facei in elements_list:
            face = np.array(facei, dtype='int32')
            if len(face) == 3:
                elem = vtkTriangle()
                epoints = elem.GetPointIds()
                epoints.SetId(0, face[0]-1)
                epoints.SetId(1, face[1]-1)
                epoints.SetId(2, face[2]-1)
                grid.InsertNextCell(5, epoints)
            elif len(face) == 4:
                elem = vtkQuad()
                epoints = elem.GetPointIds()
                epoints.SetId(0, face[0]-1)
                epoints.SetId(1, face[1]-1)
                epoints.SetId(2, face[2]-1)
                epoints.SetId(3, face[3]-1)
                grid.InsertNextCell(9, epoints)
            else:
                raise RuntimeError(face)

def _create_elements(grid, tris: np.ndarray, quads: np.ndarray) -> None:
    if len(quads):
        #elem.GetCellType() = 9  # vtkQuad
        for face in quads:
            assert len(face) == 4, face
            elem = vtkQuad()
            epoints = elem.GetPointIds()
            epoints.SetId(0, face[0])
            epoints.SetId(1, face[1])
            epoints.SetId(2, face[2])
            epoints.SetId(3, face[3])
            grid.InsertNextCell(9, epoints)

    if len(tris):
        #elem.GetCellType() = 5  # vtkTriangle
        for face in tris:
            assert len(face) == 3, face
            elem = vtkTriangle()
            epoints = elem.GetPointIds()
            epoints.SetId(0, face[0])
            epoints.SetId(1, face[1])
            epoints.SetId(2, face[2])
            grid.InsertNextCell(5, epoints)

    # TODO: make this faster...
    #etype = 5  # vtkTriangle().GetCellType()
    #create_vtk_cells_of_constant_element_type(grid, tris, etype)
