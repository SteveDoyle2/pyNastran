"""Defines the GUI IO file for Tetegen."""
import os
from collections import OrderedDict

import numpy as np

from pyNastran.converters.tetgen.tetgen import Tetgen
from pyNastran.gui.gui_objects.gui_result import GuiResult
from pyNastran.gui.utils.vtk.vtk_utils import (
    create_vtk_cells_of_constant_element_type, numpy_to_vtk_points)


class TetgenIO:
    def __init__(self, gui):
        self.gui = gui

    def get_tetgen_wildcard_geometry_results_functions(self):
        data = ('Tetgen',
                'Tetgen (*.smesh)', self.load_tetgen_geometry,
                None, None)
        return data

    def load_tetgen_geometry(self, smesh_filename, name='main', plot=True):
        model_name = name
        skip_reading = self.gui._remove_old_geometry(smesh_filename)
        if skip_reading:
            return

        model = Tetgen(log=self.gui.log, debug=False)

        base_filename, ext = os.path.splitext(smesh_filename)
        ext = ext.lower()
        node_filename = base_filename + '.node'
        ele_filename = base_filename + '.ele'
        if ext == '.smesh':
            dimension_flag = 2
        elif ext == '.ele':
            dimension_flag = 3
        else:
            raise RuntimeError('unsupported extension.  Use "smesh" or "ele".')
        model.read_tetgen(node_filename, smesh_filename, ele_filename, dimension_flag)
        nodes = model.nodes
        tris = model.tris
        tets = model.tets
        nnodes = nodes.shape[0]

        self.gui.nnodes = nodes.shape[0]
        ntris = 0
        ntets = 0
        if dimension_flag == 2:
            ntris = tris.shape[0]
        elif dimension_flag == 3:
            ntets = tets.shape[0]
        else:
            raise RuntimeError()
        nelements = ntris + ntets
        self.gui.nelements = nelements

        #print("nnodes = ",self.gui.nnodes)
        #print("nelements = ", self.gui.nelements)

        grid = self.gui.grid
        grid.Allocate(self.gui.nelements, 1000)

        assert nodes is not None
        points = numpy_to_vtk_points(nodes)
        self.gui.nid_map = {}

        #elements -= 1
        if dimension_flag == 2:
            etype = 5  # vtkTriangle().GetCellType()
            create_vtk_cells_of_constant_element_type(grid, tris, etype)
        elif dimension_flag == 3:
            etype = 10  # vtkTetra().GetCellType()
            create_vtk_cells_of_constant_element_type(grid, tets, etype)
        else:
            raise RuntimeError('dimension_flag=%r; expected=[2, 3]' % dimension_flag)

        grid.SetPoints(points)
        grid.Modified()

        # loadTetgenResults - regions/loads
        self.gui.scalar_bar_actor.VisibilityOff()
        self.gui.scalar_bar_actor.Modified()


        form, cases, node_ids, element_ids = self._fill_tetgen_case(nnodes, nelements)
        self.gui.node_ids = node_ids
        self.gui.element_ids = element_ids
        self.gui._finish_results_io2(model_name, form, cases, reset_labels=True)

    def _fill_tetgen_case(self, nnodes, nelements):
        subcase_id = 0
        self.gui.isubcase_name_map = {subcase_id : ('Tetgen', '')}

        icase = 0
        form = [
            ('NodeID', 0, []),
            ('ElementID', 1, []),
        ]
        nids = np.arange(1, nnodes+1)
        nid_res = GuiResult(subcase_id, 'NodeID', 'NodeID', 'node', nids,
                            mask_value=None, nlabels=None, labelsize=None, ncolors=None,
                            colormap='jet', data_format=None, uname='GuiResult')

        eids = np.arange(1, nelements+1)
        eid_res = GuiResult(subcase_id, 'ElementID', 'ElementID', 'centroid', eids,
                            mask_value=None, nlabels=None, labelsize=None, ncolors=None,
                            colormap='jet', data_format=None, uname='GuiResult')

        cases = OrderedDict()
        cases[icase] = (nid_res, (0, 'NodeID'))
        cases[icase + 1] = (eid_res, (0, 'ElementID'))
        return form, cases, nids, eids
