#VTK_TRIANGLE = 5

from numpy import zeros, arange

import vtk
from vtk import vtkTriangle

from pyNastran.converters.cart3d.cart3d_reader import generic_cart3d_reader


class Cart3dIO(object):
    def __init__(self):
        pass

    def removeOldGeometry(self, fileName):
        if fileName is None:
            self.grid = vtk.vtkUnstructuredGrid()
            self.gridResult = vtk.vtkFloatArray()
            #self.emptyResult = vtk.vtkFloatArray()
            #self.vectorResult = vtk.vtkFloatArray()
            self.grid2 = vtk.vtkUnstructuredGrid()
            self.scalarBar.VisibilityOff()
            skipReading = True
        else:
            self.TurnTextOff()
            self.grid.Reset()
            self.grid2.Reset()
            self.gridResult = vtk.vtkFloatArray()
            self.gridResult.Reset()
            self.gridResult.Modified()
            self.eidMap = {}
            self.nidMap = {}

            self.resultCases = {}
            self.nCases = 0
            try:
                del self.caseKeys
                del self.iCase
                del self.iSubcaseNameMap
            except:
                print("cant delete geo")
                pass

            #print(dir(self))
            skipReading = False
        self.scalarBar.VisibilityOff()
        self.scalarBar.Modified()
        return skipReading

    def load_cart3d_geometry(self, cart3dFileName, dirname, isNodal, isCentroidal):
        self.isNodal = isNodal
        self.isCentroidal = isCentroidal
        #key = self.caseKeys[self.iCase]
        #case = self.resultCases[key]

        skipReading = self.removeOldGeometry(cart3dFileName)
        if skipReading:
            return

        model = generic_cart3d_reader(cart3dFileName)
        self.modelType = model.modelType
        (nodes, elements, regions, loads) = model.read_cart3d(cart3dFileName)
        self.nNodes = model.nPoints
        self.nElements = model.nElementsRead

        #print("nNodes = ",self.nNodes)
        print("nElements = ", self.nElements)

        self.grid.Allocate(self.nElements, 1000)
        self.gridResult.SetNumberOfComponents(self.nElements)
        self.grid2.Allocate(1, 1000)

        points = vtk.vtkPoints()
        points.SetNumberOfPoints(self.nNodes)
        self.gridResult.Allocate(self.nNodes, 1000)
        #vectorReselt.SetNumberOfComponents(3)
        self.nidMap = {}
        #elem.SetNumberOfPoints(nNodes)
        if 0:
            fraction = 1. / self.nNodes  # so you can color the nodes by ID
            for nid, node in sorted(nodes.iteritems()):
                points.InsertPoint(nid - 1, *node)
                self.gridResult.InsertNextValue(nid * fraction)
                #print str(element)

                #elem = vtk.vtkVertex()
                #elem.GetPointIds().SetId(0, i)
                #self.aQuadGrid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                #vectorResult.InsertTuple3(0, 0.0, 0.0, 1.0)

        assert nodes is not None
        nnodes, three = nodes.shape
        
        nid = 0
        print "nnodes=%s" % nnodes
        for i in xrange(nnodes):
            points.InsertPoint(nid, nodes[i, :])
            nid += 1

        nelements, three = elements.shape
        elements -= 1
        for eid in xrange(nelements):
            elem = vtkTriangle()
            node_ids = elements[eid, :]
            elem.GetPointIds().SetId(0, node_ids[0])
            elem.GetPointIds().SetId(1, node_ids[1])
            elem.GetPointIds().SetId(2, node_ids[2])
            self.grid.InsertNextCell(5, elem.GetPointIds())  #elem.GetCellType() = 5  # vtkTriangle

        self.grid.SetPoints(points)
        #self.grid2.SetPoints(points2)
        #self.grid.GetPointData().SetScalars(self.gridResult)
        #print dir(self.grid) #.SetNumberOfComponents(0)
        #self.grid.GetCellData().SetNumberOfTuples(1);
        #self.grid.GetCellData().SetScalars(self.gridResult)
        self.grid.Modified()
        #self.grid2.Modified()
        self.grid.Update()
        #self.grid2.Update()
        print("updated grid")

        # loadCart3dResults - regions/loads
        self.TurnTextOn()
        self.scalarBar.VisibilityOn()
        self.scalarBar.Modified()

        from numpy import mean
        avgMach = mean(loads['Mach'])
        self.iSubcaseNameMap = {1: ['Cart3d:  avg(Mach)=%g' % avgMach, '']}
        cases = {}
        ID = 1

        #print "nElements = ",nElements
        cases = self.fillCart3dCase(cases, ID, elements, regions, loads)

        self.resultCases = cases
        self.caseKeys = sorted(cases.keys())
        #print "caseKeys = ",self.caseKeys
        #print "type(caseKeys) = ",type(self.caseKeys)
        self.iCase = -1
        self.nCases = len(self.resultCases) - 1  # number of keys in dictionary
        self.cycleResults()  # start at nCase=0

    def fillCart3dCase(self, cases, ID, elements, regions, loads):
        #print "regions**** = ",regions
        isNodal = self.isNodal
        isCentroidal = self.isCentroidal

        #nNodes = self.nNodes
        #nElements = self.nElements

        print "isCentroidal=%s isNodal=%s" % (isCentroidal, isNodal)
        assert isCentroidal != isNodal
        
        result_names = ['Cp', 'Mach', 'U', 'V', 'W', 'E', 'rho',
                                      'rhoU', 'rhoV', 'rhoW', 'rhoE']
        if isCentroidal:
            nelements, three = elements.shape
            cases[(ID, 'Region', 1, 'centroid', '%.0f')] = regions
            cases[(ID, 'Eids', 1, 'centroid', '%.0f')] = arange(1, nelements+1)
            
            for key in result_names:
                if key in loads:
                    nodal_data = loads[key]
                    n1 = elements[:, 0]
                    n2 = elements[:, 1]
                    n3 = elements[:, 2]
                    elemental_result = (nodal_data[n1] + nodal_data[n2] + nodal_data[n3])/3.0
                    cases[(ID, key, 1, 'centroid', '%.3f')] = elemental_result

        elif isNodal:
            print("load.keys() = ", loads.keys())
            for key in result_names:
                if key in loads:
                    nodal_data = loads[key]
                    cases[(ID, key, 1, 'nodal', '%.3f')] = nodal_data
        return cases
