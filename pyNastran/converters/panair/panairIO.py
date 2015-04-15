from __future__ import print_function
from six import iteritems
import os
from numpy import zeros, array, cross, dot, ravel, amax, amin
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
                'Panair (*.agps);;Panair (*.out)',  self.load_panair_results)
        return data

    def load_panair_geometry(self, panairFileName, dirname, plot=True):
        self.nidMap = {}

        #key = self.caseKeys[self.iCase]
        #case = self.resultCases[key]

        skipReading = self.removeOldGeometry(panairFileName)
        if skipReading:
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
        self.grid2.Allocate(1, 1000)

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
            #print "element = ",element
            elem = vtkQuad()
            elem.GetPointIds().SetId(0, p1)
            elem.GetPointIds().SetId(1, p2)
            elem.GetPointIds().SetId(2, p3)
            elem.GetPointIds().SetId(3, p4)
            self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())

        #print("eid = ", eid)
        self.grid.SetPoints(points)
        #self.grid2.SetPoints(points2)
        #self.grid.GetPointData().SetScalars(self.gridResult)
        #print dir(self.grid) #.SetNumberOfComponents(0)
        #self.grid.GetCellData().SetNumberOfTuples(1);
        #self.grid.GetCellData().SetScalars(self.gridResult)
        self.grid.Modified()
        self.grid2.Modified()
        if hasattr(self.grid, 'Update'):
            self.grid.Update()
            self.grid2.Update()
            print("updated grid")

        #return

        # loadCart3dResults - regions/loads
        self.TurnTextOn()
        self.scalarBar.VisibilityOn()
        self.scalarBar.Modified()

        self.iSubcaseNameMap = {1: ['Panair', '']}
        cases = {}
        ID = 1

        #print "nElements = ",nElements
        loads = []
        cases = self.fillPanairGeometryCase(cases, ID, nodes, elements, regions, loads)
        self._finish_results_io(cases)

    def clear_panair(self):
        del self.elements

    def fillPanairGeometryCase(self, cases, ID, nodes, elements, regions, loads):
        assert self.is_centroidal != self.is_nodal

        self.elements = elements
        cases[(ID, 'Region', 1, 'centroid', '%i')] = regions
        if self.is_centroidal:
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
            cases[(ID, 'centroid_x', 1, 'centroid', '%.2f')] = Xc
            cases[(ID, 'centroid_y', 1, 'centroid', '%.2f')] = Yc
            cases[(ID, 'centroid_z', 1, 'centroid', '%.2f')] = Zc
            cases[(ID, 'Area', 1, 'centroid', '%.2f')] = area
        elif self.is_nodal:
            Xn = zeros(len(nodes), dtype='float32')
            Yn = zeros(len(nodes), dtype='float32')
            Zn = zeros(len(nodes), dtype='float32')
            for i, node in enumerate(nodes):
                Xn[i] = node[0]
                Yn[i] = node[1]
                Zn[i] = node[2]
            cases[(ID, 'node_x', 1, 'node', '%.2f')] = Xn
            cases[(ID, 'node_y', 1, 'node', '%.2f')] = Yn
            cases[(ID, 'node_z', 1, 'node', '%.2f')] = Zn
        return cases

    def load_panair_results(self, panairFileName, dirname):
        if os.path.basename(panairFileName) == 'agps':
            model = AGPS(log=self.log, debug=self.debug)
            model.read_agps(panairFileName)
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

        if self.is_centroidal:
            Cp_array2 = (Cp_array[self.elements[:, 0]] +
                         Cp_array[self.elements[:, 1]] +
                         Cp_array[self.elements[:, 2]] +
                         Cp_array[self.elements[:, 3]]) / 4.
            key = (1, 'Cp', 1, 'centroid', '%.3f')
            self.resultCases[key] = Cp_array2

        elif self.is_nodal:
            key = (1, 'Cp', 1, 'node', '%.3f')
            self.resultCases[key] = Cp_array

        #self.resultCases = cases
        self.caseKeys = sorted(self.resultCases.keys())
        self.iCase = -1
        self.nCases = len(self.resultCases) - 1  # number of keys in dictionary
        #self.nCases = 1
        self.cycleResults()  # start at nCase=0

def main():
    def removeOldGeometry(self):
        pass
    def cycleResults(self):
        pass

    test = PanairIO()
    test.is_nodal = True
    test.is_centroidal = False
    test.removeOldGeometry = removeOldGeometry
    test.cycleResults = cycleResults

    #test.load_panair_geometry('SWB.INP','')
    test.load_panair_geometry('models/NAC6.INP', '')

if __name__ == '__main__':  # pragma: no cover
    main()


#if __name__=='__main__':
#    lawgs = LaWGS('tmx1242.wgs')
#    lawgs.run()
