"""Defines the GUI IO file for STL."""
from collections import OrderedDict
from numpy import arange

import vtk

from pyNastran.converters.stl.stl import read_stl
from pyNastran.gui.gui_objects.gui_result import GuiResult, NormalResult
from pyNastran.gui.utils.vtk.vtk_utils import (
    create_vtk_cells_of_constant_element_type, numpy_to_vtk_points)


class STL_IO:
    def __init__(self, gui):
        self.gui = gui

    def get_stl_wildcard_geometry_results_functions(self):
        data = ('STL',
                'STereoLithography (*.STL)', self.load_stl_geometry,
                None, None)
        return data

    def load_stl_geometry(self, stl_filename, name='main', plot=True):
        model_name = name
        skip_reading = self.gui._remove_old_geometry(stl_filename)
        if skip_reading:
            return

        model = read_stl(stl_filename, remove_elements_with_bad_normals=True,
                         log=self.gui.log, debug=False)
        #self.model_type = model.model_type
        nodes = model.nodes
        elements = model.elements

        normals = model.get_normals(elements, stop_on_failure=False)
        areas = model.get_area(elements, stop_on_failure=False)
        #nnormals = model.get_normals_at_nodes(elements)
        self.gui.nnodes = nodes.shape[0]
        self.gui.nelements = elements.shape[0]

        self.gui.log.info('nnodes=%s nelements=%s' % (self.gui.nnodes, self.gui.nelements))
        grid = self.gui.grid
        grid.Allocate(self.gui.nelements, 1000)

        points = numpy_to_vtk_points(nodes)
        self.gui.nid_map = {}
        #elem.SetNumberOfPoints(nnodes)

        assert nodes is not None
        unused_nnodes = nodes.shape[0]
        xmax, ymax, zmax = nodes.max(axis=0)
        xmin, ymin, zmin = nodes.min(axis=0)
        self.gui.log.info('xmax=%s xmin=%s' % (xmax, xmin))
        self.gui.log.info('ymax=%s ymin=%s' % (ymax, ymin))
        self.gui.log.info('zmax=%s zmin=%s' % (zmax, zmin))
        dim_max = max(xmax-xmin, ymax-ymin, zmax-zmin)
        self.gui.create_global_axes(dim_max)


        etype = 5  # vtkTriangle().GetCellType()
        create_vtk_cells_of_constant_element_type(grid, elements, etype)

        grid.SetPoints(points)
        grid.Modified()

        # loadSTLResults - regions/loads
        self.gui.scalar_bar_actor.VisibilityOff()
        self.gui.scalar_bar_actor.Modified()

        cases = OrderedDict()
        self.gui.isubcase_name_map = {}
        ID = 1

        form, cases, node_ids, element_ids = self._fill_stl_case(
            cases, ID, elements, nodes, normals, areas)
        self.gui.node_ids = node_ids
        self.gui.element_ids = element_ids
        self.gui._finish_results_io2(model_name, form, cases)

    def _fill_stl_case(self, cases, ID, elements, nodes, normals, areas):
        """adds the sidebar results"""
        self.gui.isubcase_name_map[ID] = ('STL', '')
        colormap = 'jet'
        nelements = elements.shape[0]
        nnodes = nodes.shape[0]
        icase = 0
        itime = 0
        eids = arange(1, nelements + 1, dtype='int32')
        nids = arange(1, nnodes + 1, dtype='int32')
        eid_res = GuiResult(ID, header='ElementID', title='ElementID',
                            location='centroid', scalar=eids)
        nid_res = GuiResult(ID, header='NodeID', title='NodeID',
                            location='node', scalar=nids)
        area_res = GuiResult(ID, header='Area', title='Area',
                             location='centroid', scalar=areas)
        nx_res = GuiResult(ID, header='NormalX', title='NormalX',
                           location='centroid', scalar=normals[:, 0])
        ny_res = GuiResult(ID, header='NormalY', title='NormalY',
                           location='centroid', scalar=normals[:, 1])
        nz_res = GuiResult(ID, header='NormalZ', title='NormalZ',
                           location='centroid', scalar=normals[:, 2])
        nxyz_res = NormalResult(0, 'Normals', 'Normals',
                                nlabels=2, labelsize=5, ncolors=2,
                                colormap=colormap, data_format='%.1f',
                                uname='NormalResult')
        cases[icase] = (eid_res, (itime, 'ElementID'))
        cases[icase + 1] = (nid_res, (itime, 'NodeID'))
        cases[icase + 2] = (area_res, (itime, 'Area'))
        cases[icase + 3] = (nx_res, (itime, 'NormalX'))
        cases[icase + 4] = (ny_res, (itime, 'NormalY'))
        cases[icase + 5] = (nz_res, (itime, 'NormalZ'))
        cases[icase + 6] = (nxyz_res, (itime, 'Normal'))

        form = [
            ('ElementID', icase, []),
            ('NodeID', icase + 1, []),
            ('Area', icase + 2, []),
            ('NormalX', icase + 3, []),
            ('NormalY', icase + 4, []),
            ('NormalZ', icase + 5, []),
            ('Normal', icase + 6, []),
        ]
        return form, cases, nids, eids
