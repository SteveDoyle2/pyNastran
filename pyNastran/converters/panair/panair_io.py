from __future__ import print_function
from six import iteritems
import os
import numpy as np
from numpy import zeros, array, cross, ravel, amax, amin, arange
from numpy.linalg import norm

import vtk
from vtk import vtkQuad

from pyNastran.converters.panair.panair_grid import PanairGrid
from pyNastran.converters.panair.agps import AGPS
from pyNastran.gui.gui_objects.gui_result import GuiResult


class PanairIO(object):
    def __init__(self):
        pass

    def get_panair_wildcard_geometry_results_functions(self):
        data = ('Panair',
                'Panair (*.inp)', self.load_panair_geometry,
                #'Panair (*.agps);;Panair (*.out)',  self.load_panair_results)
                'Panair (*agps)', self.load_panair_results)
        return data

    def load_panair_geometry(self, panair_filename, dirname, name='main', plot=True):
        self.nid_map = {}

        #key = self.case_keys[self.icase]
        #case = self.result_cases[key]

        skip_reading = self._remove_old_geometry(panair_filename)
        if skip_reading:
            return

        model = PanairGrid(log=self.log, debug=self.debug)
        self.model_type = model.model_type
        model.read_panair(panair_filename)

        nodes, elements, regions, kt, cp_norm = model.get_points_elements_regions()
        #for nid, node in enumerate(nodes):
            #print "node[%s] = %s" % (nid, str(node))

        self.nNodes = len(nodes)
        self.nElements = len(elements)

        #print("nNodes = ",self.nNodes)
        #print("nElements = ", self.nElements)

        self.grid.Allocate(self.nElements, 1000)
        #self.gridResult.SetNumberOfComponents(self.nElements)

        points = vtk.vtkPoints()
        points.SetNumberOfPoints(self.nNodes)
        #self.gridResult.Allocate(self.nNodes, 1000)
        #vectorReselt.SetNumberOfComponents(3)
        #elem.SetNumberOfPoints(nNodes)
        if 0:
            fraction = 1. / nnodes  # so you can color the nodes by ID
            for nid, node in sorted(iteritems(nodes)):
                points.InsertPoint(nid - 1, *point)
                self.gridResult.InsertNextValue(nid * fraction)
                #print str(element)

                #elem = vtk.vtkVertex()
                #elem.GetPointIds().SetId(0, i)
                #self.aQuadGrid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                #vectorResult.InsertTuple3(0, 0.0, 0.0, 1.0)

        assert len(nodes) > 0
        mmax = amax(nodes, axis=0)
        mmin = amin(nodes, axis=0)
        dim_max = (mmax - mmin).max()
        self.create_global_axes(dim_max)
        for nid, node in enumerate(nodes):
            points.InsertPoint(nid, *node)

        assert len(elements) > 0
        elem = vtkQuad()
        quad_type = elem.GetCellType()
        for eid, element in enumerate(elements):
            (p1, p2, p3, p4) = element
            elem = vtkQuad()
            elem.GetPointIds().SetId(0, p1)
            elem.GetPointIds().SetId(1, p2)
            elem.GetPointIds().SetId(2, p3)
            elem.GetPointIds().SetId(3, p4)
            self.grid.InsertNextCell(quad_type, elem.GetPointIds())

        self.grid.SetPoints(points)
        #self.grid.GetPointData().SetScalars(self.gridResult)
        #print dir(self.grid) #.SetNumberOfComponents(0)
        #self.grid.GetCellData().SetNumberOfTuples(1);
        #self.grid.GetCellData().SetScalars(self.gridResult)
        self.grid.Modified()
        if hasattr(self.grid, 'Update'):
            self.grid.Update()

        # loadPanairResults - regions/loads
        self. turn_text_on()
        if plot:
            self.scalarBar.VisibilityOn()
            self.scalarBar.Modified()

        self.iSubcaseNameMap = {1: ['Panair', '']}
        cases = {}
        ID = 1

        #print "nElements = ", nElements
        loads = []
        form, cases = self._fill_panair_geometry_case(cases, ID, nodes, elements, regions, kt, cp_norm, loads)
        #if plot:
        self._finish_results_io2(form, cases)
        #else:
        #self._set_results([form], cases)

    def clear_panair(self):
        del self.elements

    def _fill_panair_geometry_case(self, cases, ID, nodes, elements, regions, kt, cp_norm, loads):
        self.elements = elements
        nnids = nodes.shape[0]
        neids = elements.shape[0]
        nids = arange(0., nnids, dtype='int32') + 1
        eids = arange(0., neids, dtype='int32') + 1

        icase = 0
        location_form = [
            ('normal_x', icase + 4, []),
            ('normal_y', icase + 5, []),
            ('normal_z', icase + 6, []),

            ('centroid_x', icase + 7, []),
            ('centroid_y', icase + 8, []),
            ('centroid_z', icase + 9, []),

            ('node_x', icase + 10, []),
            ('node_y', icase + 11, []),
            ('node_z', icase + 12, []),
        ]

        geometry_form = [
            ('Patch', icase, []),
            ('ElementID', icase + 1, []),
            ('NodeID', icase + 2, []),
            ('Area', icase + 3, []),
            ('Location', None, location_form),
            ('Kt', icase + 13, []),
            ('CpNorm', icase + 14, []),
        ]
        form = [
            ('Geometry', None, geometry_form),
        ]

        # centroidal
        nelements = len(elements)
        #print('nelements = ', nelements)
        #print('nnodes = ', nodes.shape)

        p1 = nodes[elements[:, 0]]
        p2 = nodes[elements[:, 1]]
        p3 = nodes[elements[:, 2]]
        p4 = nodes[elements[:, 3]]
        xyz_centroid = (p1 + p2 + p3 + p4) / 4.

        n = np.cross(p3 - p1, p4 - p2)
        n_norm = np.linalg.norm(n, axis=1)
        area = n_norm / 2.
        normal = n / n_norm[:, np.newaxis]

        itime = 0
        # ID, header, title, location, values, format, uname
        region_res = GuiResult(ID, 'Region', 'Region', 'centroid', regions, data_format=None, uname='Region')
        eid_res = GuiResult(ID, 'ElementID', 'ElementID', 'centroid', eids, data_format=None, uname='ElementID')
        nid_res = GuiResult(ID, 'NodeID', 'NodeID', 'node', nids, data_format=None, uname='NodeID')
        area_res = GuiResult(ID, 'Area', 'Area', 'centroid', area, data_format=None, uname='Area')

        nx_res = GuiResult(ID, 'normal_x', 'NormalX', 'centroid', normal[:, 0], data_format='%.3f', uname='NormalX')
        ny_res = GuiResult(ID, 'normal_y', 'NormalY', 'centroid', normal[:, 1], data_format='%.3f', uname='NormalY')
        nz_res = GuiResult(ID, 'normal_z', 'NormalZ', 'centroid', normal[:, 2], data_format='%.3f', uname='NormalZ')
        cenx_res = GuiResult(ID, 'centroid_x', 'CentroidX', 'centroid', xyz_centroid[:, 0], data_format=None, uname='CentroidX')
        ceny_res = GuiResult(ID, 'centroid_y', 'CentroidY', 'centroid', xyz_centroid[:, 1], data_format=None, uname='CentroidY')
        cenz_res = GuiResult(ID, 'centroid_z', 'CentroidZ', 'centroid', xyz_centroid[:, 2], data_format=None, uname='CentroidZ')

        kt_res = GuiResult(ID, 'Kt', 'Kt', 'centroid', kt, data_format=None, uname='Kt')
        cp_res = GuiResult(ID, 'CpNorm', 'CpNorm', 'centroid', cp_norm, data_format=None, uname='CpNorm')

        cases[icase] = (region_res, (itime, 'Region'))
        cases[icase + 1] = (eid_res, (itime, 'ElementID'))
        cases[icase + 2] = (nid_res, (itime, 'NodeID'))
        cases[icase + 3] = (area_res, (itime, 'Area'))

        cases[icase + 4] = (nx_res, (itime, 'NormalX'))
        cases[icase + 5] = (ny_res, (itime, 'NormalY'))
        cases[icase + 6] = (nz_res, (itime, 'NormalZ'))
        cases[icase + 7] = (cenx_res, (itime, 'CentroidX'))
        cases[icase + 8] = (ceny_res, (itime, 'CentroidY'))
        cases[icase + 9] = (cenz_res, (itime, 'CentroidZ'))

        cases[icase + 13] = (kt_res, (itime, 'Kt'))
        cases[icase + 14] = (cp_res, (itime, 'CpNorm'))

        # nodal
        cenx_res = GuiResult(ID, 'node_x', 'NodeX', 'node', nodes[:, 0], data_format='%.2f', uname='node_x')
        ceny_res = GuiResult(ID, 'node_y', 'NodeY', 'node', nodes[:, 1], data_format='%.2f', uname='node_y')
        cenz_res = GuiResult(ID, 'node_z', 'NodeZ', 'node', nodes[:, 2], data_format='%.2f', uname='node_z')
        cases[icase + 10] = (cenx_res, (itime, 'node_x'))
        cases[icase + 11] = (ceny_res, (itime, 'node_y'))
        cases[icase + 12] = (cenz_res, (itime, 'node_z'))

        return form, cases

    def load_panair_results(self, panair_filename, dirname):
        cases = self.result_cases
        if os.path.basename(panair_filename) == 'agps':
            model = AGPS(log=self.log, debug=self.debug)
            model.read_agps(panair_filename)
        else:
            raise RuntimeError('only files named "agps" files are supported')

        # get the Cp on the nodes
        Cp_array = zeros(self.nNodes, dtype='float32')
        imin = 0
        for ipatch, Cp in sorted(iteritems(model.pressures)):
            Cpv = ravel(Cp)
            nCp = len(Cpv)
            try:
                Cp_array[imin:imin + nCp] = Cpv
            except ValueError:
                # agps stores implicit and explicit wakes
                # we're skipping all wakes
                pass
            imin += nCp

        Cp_array2 = (Cp_array[self.elements[:, 0]] +
                     Cp_array[self.elements[:, 1]] +
                     Cp_array[self.elements[:, 2]] +
                     Cp_array[self.elements[:, 3]]) / 4.

        icase = len(self.result_cases)

        form = self.get_form()
        results_form = [
            ('Cp Nodal', icase, []),
            ('Cp Centroidal', icase + 1, [],),
        ]
        form.append(('Results', None, results_form))

        ID = 1
        Cpn_res = GuiResult(ID, 'Cp Nodal', 'Cp', 'node', Cp_array, data_format='%.3f', uname='Cp_nodal')
        Cpc_res = GuiResult(ID, 'Cp Centroidal', 'Cp', 'centroid', Cp_array2, data_format='%.3f', uname='Cp_centroidal')
        cases[icase] = (Cpn_res, (0, 'Cp_nodal'))
        cases[icase + 1] = (Cpc_res, (0, 'Cp_centroidal'))

        self._finish_results_io2(form, cases)
