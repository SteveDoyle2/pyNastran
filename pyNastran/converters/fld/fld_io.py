"""Defines the GUI IO file for STL."""
from __future__ import annotations
from typing import Any, TYPE_CHECKING
import numpy as np

import vtkmodules

from pyNastran.converters.fld.fld import read_fld
from pyNastran.gui.gui_objects.gui_result import GuiResult, NormalResult
from pyNastran.gui.utils.vtk.vtk_utils import (
    create_vtk_cells_of_constant_element_type, numpy_to_vtk_points)
if TYPE_CHECKING:
    from pyNastran.gui.main_window import MainWindow

class FLD_IO:
    def __init__(self, gui: MainWindow):
        self.gui = gui

    def get_fld_wildcard_geometry_results_functions(self):
        data = ('FLD',
                'Fluent Field (*.fld)', self.load_fld_geometry,
                None, None)
        return data

    def load_fld_geometry(self, fld_filename: str, name: str='main', plot: bool=True):
        model_name = name
        skip_reading = self.gui._remove_old_geometry(fld_filename)
        if skip_reading:
            return

        model = read_fld(fld_filename, log=self.gui.log, debug=False)
        #self.model_type = model.model_type
        nodes = model.xyzp[:, :3]
        pressure = model.xyzp[:, 3:]

        nnodes = len(nodes)
        nelements = nnodes

        element_ids = np.arange(nnodes)
        node_ids = element_ids
        elements = element_ids.reshape(nelements, 1)

        self.gui.nnodes = nnodes
        self.gui.nelements = len(element_ids)

        self.gui.log.info('nnodes=%s nelements=%s' % (self.gui.nnodes, self.gui.nelements))
        grid = self.gui.grid
        grid.Allocate(self.gui.nelements, 1000)

        points = numpy_to_vtk_points(nodes)
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


        etype = 1  # vtkVertex().GetCellType()
        create_vtk_cells_of_constant_element_type(grid, elements, etype)

        grid.SetPoints(points)
        grid.Modified()

        # loadSTLResults - regions/loads
        self.gui.scalar_bar_actor.VisibilityOff()
        self.gui.scalar_bar_actor.Modified()

        cases = {}
        self.gui.isubcase_name_map = {}
        ID = 1

        form, cases = self._fill_fld_case(
            cases, ID, node_ids, pressure)
        self.gui.node_ids = node_ids
        self.gui.element_ids = element_ids
        self.gui._finish_results_io2(model_name, form, cases)

    def _fill_fld_case(self, cases: dict[int, Any], ID: int,
                       node_ids: np.ndarray, pressure: np.ndarray) -> None:
        """adds the sidebar results"""
        self.gui.isubcase_name_map[ID] = ('FLD', '')
        colormap = 'jet'
        icase = 0
        itime = 0
        nid_res = GuiResult(ID, header='NodeID', title='NodeID',
                            location='node', scalar=node_ids)
        pressure_res = GuiResult(ID, header='Pressure', title='Pressure',
                            location='node', scalar=pressure)
        cases[icase] = (nid_res, (itime, 'NodeID'))
        cases[icase + 1] = (pressure_res, (itime, 'Pressure'))

        form = [
            ('NodeID', icase, []),
            ('Pressure', icase + 1, []),
        ]
        return form, cases
