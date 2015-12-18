from __future__ import print_function
from six import iteritems
import os
from numpy import zeros, array, cross, dot, ravel, amax, amin, arange
from numpy.linalg import det, norm

import vtk
from vtk import vtkQuad

from pyNastran.converters.panair.panairGrid import PanairGrid
from pyNastran.converters.panair.agps import AGPS

class PanairIO(object):
    def __init__(self):
        pass

    def get_panair_wildcard_geometry_results_functions(self):
        data = ('Panair',
                'Panair (*.inp)', self.load_panair_geometry,
                #'Panair (*.agps);;Panair (*.out)',  self.load_panair_results)
                'Panair (*agps)', self.load_panair_results)
        return data

    def load_panair_geometry(self, panairFileName, dirname, plot=True):
        self.nidMap = {}

        #key = self.caseKeys[self.iCase]
        #case = self.resultCases[key]

        skip_reading = self.removeOldGeometry(panairFileName)
        if skip_reading:
            return

        model = PanairGrid(log=self.log, debug=self.debug)
        self.modelType = model.modelType
        model.read_panair(panairFileName)

        nodes, elements, regions = model.getPointsElementsRegions()
        #for nid,node in enumerate(nodes):
            #print "node[%s] = %s" %(nid,str(node))

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
            fraction = 1. / nNodes  # so you can color the nodes by ID
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
        self.update_axes_length(dim_max)
        for nid, node in enumerate(nodes):
            points.InsertPoint(nid, *node)

        assert len(elements) > 0
        for eid, element in enumerate(elements):
            (p1, p2, p3, p4) = element
            elem = vtkQuad()
            elem.GetPointIds().SetId(0, p1)
            elem.GetPointIds().SetId(1, p2)
            elem.GetPointIds().SetId(2, p3)
            elem.GetPointIds().SetId(3, p4)
            self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())

        self.grid.SetPoints(points)
        #self.grid.GetPointData().SetScalars(self.gridResult)
        #print dir(self.grid) #.SetNumberOfComponents(0)
        #self.grid.GetCellData().SetNumberOfTuples(1);
        #self.grid.GetCellData().SetScalars(self.gridResult)
        self.grid.Modified()
        if hasattr(self.grid, 'Update'):
            self.grid.Update()

        # loadPanairResults - regions/loads
        self.TurnTextOn()
        self.scalarBar.VisibilityOn()
        self.scalarBar.Modified()

        self.iSubcaseNameMap = {1: ['Panair', '']}
        cases = {}
        ID = 1

        #print "nElements = ",nElements
        loads = []
        form, cases = self._fill_panair_geometry_case(cases, ID, nodes, elements, regions, loads)
        self._finish_results_io2(form, cases)
        #self._finish_results_io(cases)

    def clear_panair(self):
        del self.elements

    def _fill_panair_geometry_case(self, cases, ID, nodes, elements, regions, loads):
        self.elements = elements
        nnids = nodes.shape[0]
        neids = elements.shape[0]
        nids = arange(0., nnids, dtype='int32') + 1
        eids = arange(0., neids, dtype='int32') + 1

        icase = 0
        location_form = [
            ('centroid_x', icase + 4, []),
            ('centroid_y', icase + 5, []),
            ('centroid_z', icase + 6, []),

            ('node_x', icase + 7, []),
            ('node_y', icase + 8, []),
            ('node_z', icase + 9, []),
        ]

        geometry_form = [
            ('Patch', icase, []),
            ('ElementID', icase + 1, []),
            ('NodeID', icase + 2, []),
            ('Area', icase + 3, []),
            ('Location', None, location_form),
        ]
        form = [
            ('Geometry', None, geometry_form),
        ]

        # centroidal
        Xc = zeros(len(elements), dtype='float32')
        Yc = zeros(len(elements), dtype='float32')
        Zc = zeros(len(elements), dtype='float32')
        area = zeros(len(elements), dtype='float32')
        for i, element in enumerate(elements):
            p1, p2, p3, p4 = element
            P1 = array(nodes[p1])
            P2 = array(nodes[p2])
            P3 = array(nodes[p3])
            P4 = array(nodes[p4])
            a = P3 - P1
            b = P4 - P2
            A = 0.5 * norm(cross(a, b))
            x, y, z = (P1 + P2 + P3 + P4) / 4.0
            Xc[i] = x
            Yc[i] = y
            Zc[i] = z
            area[i] = A

        cases[(ID, icase, 'Region', 1, 'centroid', '%i', '')] = regions
        cases[(ID, icase + 1, 'ElementID', 1, 'centroid', '%i', '')] = eids
        cases[(ID, icase + 2, 'NodeID', 1, 'node', '%i', '')] = nids
        cases[(ID, icase + 3, 'Area', 1, 'centroid', '%i', '')] = area
        cases[(ID, icase + 4, 'centroid_x', 1, 'centroid', '%.2f', '')] = Xc
        cases[(ID, icase + 5, 'centroid_y', 1, 'centroid', '%.2f', '')] = Yc
        cases[(ID, icase + 6, 'centroid_z', 1, 'centroid', '%.2f', '')] = Zc

        # nodal
        Xn = zeros(len(nodes), dtype='float32')
        Yn = zeros(len(nodes), dtype='float32')
        Zn = zeros(len(nodes), dtype='float32')
        for i, node in enumerate(nodes):
            Xn[i] = node[0]
            Yn[i] = node[1]
            Zn[i] = node[2]
        cases[(ID, icase + 7, 'node_x', 1, 'node', '%.2f', '')] = Xn
        cases[(ID, icase + 8, 'node_y', 1, 'node', '%.2f', '')] = Yn
        cases[(ID, icase + 9, 'node_z', 1, 'node', '%.2f', '')] = Zn

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

        #key = (1, 'Cp', 1, 'centroid', '%.3f')
        #self.result_cases[key] = Cp_array2
        icase = len(self.result_cases)

        form = self.get_form()
        results_form = [
            ('Cp', icase, []),
            ('Cp_centroidal', icase + 1, [],),
        ]
        form.append(('Results', None, results_form))
        key = (1, icase, 'Cp', 1, 'node', '%.3f', '')
        self.result_cases[key] = Cp_array

        key = (1, icase + 1, 'Cp_centroidal', 1, 'centroid', '%.3f', '')
        self.result_cases[key] = Cp_array2

        self._finish_results_io2(form, cases)
