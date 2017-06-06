"""
Defines the GUI IO file for STL.
"""
from __future__ import print_function
from six import iteritems
from six.moves import range
from numpy import arange

import vtk
from vtk import vtkTriangle

from pyNastran.converters.stl.stl import read_stl
from pyNastran.gui.gui_objects.gui_result import GuiResult


class STL_IO(object):
    def __init__(self):
        pass

    def get_stl_wildcard_geometry_results_functions(self):
        data = ('STL',
                'STereoLithography (*.STL)', self.load_stl_geometry,
                None, None)
        return data

    def load_stl_geometry(self, stl_filename, dirname, name='main', plot=True):
        #print("load_stl_geometry...")
        skip_reading = self._remove_old_geometry(stl_filename)
        if skip_reading:
            return

        model = read_stl(stl_filename, log=self.log, debug=False)
        #self.model_type = model.model_type
        nodes = model.nodes
        elements = model.elements

        normals = model.get_normals(elements, stop_on_failure=False)
        areas = model.get_area(elements, stop_on_failure=False)
        #nnormals = model.get_normals_at_nodes(elements)
        self.nNodes = nodes.shape[0]
        self.nElements = elements.shape[0]

        self.log.info('nnodes=%s nelements=%s' % (self.nNodes, self.nElements))
        grid = self.grid
        grid.Allocate(self.nElements, 1000)
        #self.gridResult.SetNumberOfComponents(self.nElements)

        points = self.numpy_to_vtk_points(nodes)
        self.nid_map = {}
        #elem.SetNumberOfPoints(nnodes)
        if 0:
            fraction = 1. / self.nNodes  # so you can color the nodes by ID
            for nid, node in sorted(iteritems(nodes)):
                self.gridResult.InsertNextValue(nid * fraction)
                #print str(element)

                #elem = vtk.vtkVertex()
                #elem.GetPointIds().SetId(0, i)
                #self.aQuadGrid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                #vectorResult.InsertTuple3(0, 0.0, 0.0, 1.0)

        assert nodes is not None
        nnodes = nodes.shape[0]
        xmax, ymax, zmax = nodes.max(axis=0)
        xmin, ymin, zmin = nodes.min(axis=0)
        self.log.info('xmax=%s xmin=%s' % (xmax, xmin))
        self.log.info('ymax=%s ymin=%s' % (ymax, ymin))
        self.log.info('zmax=%s zmin=%s' % (zmax, zmin))
        dim_max = max(xmax-xmin, ymax-ymin, zmax-zmin)
        self.create_global_axes(dim_max)


        etype = 5  # vtkTriangle().GetCellType()
        self.create_vtk_cells_of_constant_element_type(grid, elements, etype)

        grid.SetPoints(points)
        grid.Modified()
        if hasattr(grid, 'Update'):
            grid.Update()
            self.log_info("updated grid")

        # loadSTLResults - regions/loads
        self.scalarBar.VisibilityOff()
        self.scalarBar.Modified()

        cases = {}
        self.iSubcaseNameMap = {}
        ID = 1

        form, cases = self._fill_stl_case(cases, ID, elements, nodes, normals, areas)
        self._finish_results_io2(form, cases)

    def _fill_stl_case(self, cases, ID, elements, nodes, normals, areas):
        """adds the sidebar results"""
        self.iSubcaseNameMap[ID] = ('STL', '')

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
        cases[icase] = (eid_res, (itime, 'ElementID'))
        cases[icase + 1] = (nid_res, (itime, 'NodeID'))
        cases[icase + 2] = (area_res, (itime, 'Area'))
        cases[icase + 3] = (nx_res, (itime, 'NormalX'))
        cases[icase + 4] = (ny_res, (itime, 'NormalY'))
        cases[icase + 5] = (nz_res, (itime, 'NormalZ'))

        form = [
            ('ElementID', icase, []),
            ('NodeID', icase + 1, []),
            ('Area', icase + 2, []),
            ('NormalX', icase + 3, []),
            ('NormalY', icase + 4, []),
            ('NormalZ', icase + 5, []),
        ]
        return form, cases
