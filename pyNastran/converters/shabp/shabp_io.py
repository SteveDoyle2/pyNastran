from six import iteritems
import os
from numpy import zeros, array, cross, dot, ravel, amax, amin
from numpy.linalg import det, norm

import vtk
from vtk import vtkQuad

from pyNastran.converters.shabp.shabp import SHABP

is_shabp = True


class ShabpIO(object):
    def __init__(self):
        pass

    def get_shabp_wildcard_geometry_results_functions(self):
        data = ('S/HABP',
                'Shabp (*.geo; *.mk5; *.inp)', self.load_shabp_geometry,
                'Shabp (*.out)', self.load_shabp_results)
        return data

    def load_shabp_geometry(self, shabpFilename, dirname, plot=True):
        self.nidMap = {}

        #key = self.caseKeys[self.iCase]
        #case = self.resultCases[key]

        skipReading = self.removeOldGeometry(shabpFilename)
        if skipReading:
            return

        self.model = SHABP(log=self.log, debug=self.debug)
        self.modelType = 'shabp' # model.modelType
        self.model.read_shabp(shabpFilename)

        nodes, elements, patches, components, impact, shadow = self.model.getPointsElementsRegions()
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
                points.InsertPoint(nid - 1, *node)
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
            #print("updated grid")

        #return

        # loadCart3dResults - regions/loads
        self.TurnTextOn()
        self.scalarBar.VisibilityOn()
        self.scalarBar.Modified()

        self.iSubcaseNameMap = {1: ['S/HABP', '']}
        cases = {}
        ID = 1

        self.log.debug("nNodes=%i nElements=%i" % (self.nNodes, self.nElements))
        cases = self.fillShabpGeometryCase(cases, ID, nodes, elements, patches, components, impact, shadow)
        self._finish_results_io(cases)

    def clear_shabp(self):
        del self.elements
        del self.model

    def fillShabpGeometryCase(self, cases, ID, nodes, elements, patches, components, impact, shadow):
        assert self.is_centroidal != self.is_nodal

        self.elements = elements
        if self.is_centroidal:
            cases[(ID, 'Component', 1, 'centroid', '%i')] = components
            cases[(ID, 'PatchID', 1, 'centroid', '%i')] = patches
            cases[(ID, 'Impact', 1, 'centroid', '%i')] = impact
            cases[(ID, 'Shadow', 1, 'centroid', '%i')] = shadow

            XYZc = zeros((len(elements),3), dtype='float32')
            #Normal = zeros((len(elements),3), dtype='float32')
            area = zeros(len(elements), dtype='float32')

            for i,element in enumerate(elements):
                p1, p2, p3, p4 = element
                P1 = array(nodes[p1])
                P2 = array(nodes[p2])
                P3 = array(nodes[p3])
                P4 = array(nodes[p4])
                a = P3 - P1
                b = P4 - P2
                n = cross(a, b)
                nnorm = norm(n)
                #normal = n / nnorm
                A = 0.5 * nnorm

                XYZc[i,:] = (P1 + P2 + P3 + P4) / 4.0
                #Normal[i, :] = normal
                area[i] = A
            cases[(ID, 'centroid_x', 1, 'centroid', '%.2f')] = XYZc[:,0]
            cases[(ID, 'centroid_y', 1, 'centroid', '%.2f')] = XYZc[:,1]
            cases[(ID, 'centroid_z', 1, 'centroid', '%.2f')] = XYZc[:,2]

            #cases[(ID, 'normal_x', 1, 'centroid', '%.2f')] = Normal[:,0]
            #cases[(ID, 'normal_y', 1, 'centroid', '%.2f')] = Normal[:,1]
            #cases[(ID, 'normal_z', 1, 'centroid', '%.2f')] = Normal[:,2]
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

    def load_shabp_results(self, shabp_filename, dirname):
        Cpd, deltad = self.model.read_shabp_out(shabp_filename)

        if self.is_centroidal:
            self.resultCases = {}
            for case_id, Cp in sorted(iteritems(Cpd)):
                Cp = Cpd[case_id]
                #delta = deltad[case_id]

                mach, alpha, beta = self.model.shabp_cases[case_id]
                #name = 'Mach=%g Alpha=%g' % (mach, alpha)
                name = 'Mach=%g Alpha=%g' % (mach, alpha)
                self.resultCases[(name, 'Cp', 1, 'centroid', '%.3f')] = Cp
                #self.resultCases[(name, 'delta', 1, 'centroid', '%.3f')] = delta
        elif self.is_nodal:
            #key = (1, 'Cp', 1, 'node', '%.3f')
            #self.resultCases[key] = Cp_array
            pass
        self._finish_results_io(self.resultCases)

def main():
    def removeOldGeometry(self):
        pass
    def cycleResults(self):
        pass

    test = ShabpIO()
    test.is_nodal = True
    test.is_centroidal = False
    test.removeOldGeometry = removeOldGeometry
    test.cycleResults = cycleResults

    #test.load_shabp_geometry('SWB.INP','')
    test.load_shabp_geometry('models/NAC6.INP', '')

if __name__ == '__main__':  # pragma: no cover
    main()


#if __name__=='__main__':
#    lawgs = LaWGS('tmx1242.wgs')
#    lawgs.run()
