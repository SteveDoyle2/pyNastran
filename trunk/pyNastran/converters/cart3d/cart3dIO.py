#VTK_TRIANGLE = 5

from numpy import zeros

import vtk
from vtk import vtkTriangle

from pyNastran.converters.cart3d.cart3d_reader import genericCart3DReader


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
                print "cant delete geo"
                pass
            ###
            #print dir(self)
            skipReading = False
        self.scalarBar.VisibilityOff()
        self.scalarBar.Modified()
        return skipReading

    def loadCart3dGeometry(self, cart3dFileName, dirname, isNodal, isCentroidal):
        self.isNodal = isNodal
        self.isCentroidal = isCentroidal
        #key = self.caseKeys[self.iCase]
        #case = self.resultCases[key]

        skipReading = self.removeOldGeometry(cart3dFileName)
        if skipReading:
            return

        model = genericCart3DReader(cart3dFileName)
        self.modelType = model.modelType
        (nodes, elements, regions, loads) = model.readCart3d(cart3dFileName)
        self.nNodes = model.nPoints
        self.nElements = model.nElementsRead

        #print "nNodes = ",self.nNodes
        print "nElements = ", self.nElements

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
        for nid, node in sorted(nodes.iteritems()):
            points.InsertPoint(nid - 1, *node)

        for eid, nodeIDs in sorted(elements.iteritems()):
            #print "ctria3"
            elem = vtkTriangle()
            elem.GetPointIds().SetId(0, nodeIDs[0] - 1)
            elem.GetPointIds().SetId(1, nodeIDs[1] - 1)
            elem.GetPointIds().SetId(2, nodeIDs[2] - 1)
            self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
        ###
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
        print "updated grid"

        # loadCart3dResults - regions/loads
        self.TurnTextOn()
        self.scalarBar.VisibilityOn()
        self.scalarBar.Modified()

        self.iSubcaseNameMap = {1: ['Cart3d', '']}
        cases = {}
        ID = 1

        #print "nElements = ",nElements
        cases = self.fillCart3dCase(cases, ID, regions, loads)

        self.resultCases = cases
        self.caseKeys = sorted(cases.keys())
        #print "caseKeys = ",self.caseKeys
        #print "type(caseKeys) = ",type(self.caseKeys)
        self.iCase = -1
        self.nCases = len(self.resultCases) - 1  # number of keys in dictionary
        self.cycleResults()  # start at nCase=0

    def fillCart3dCase(self, cases, ID, regions, loads):
        #print "regions**** = ",regions
        plotNodal = self.isNodal
        plotCentroidal = self.isCentroidal

        nNodes = self.nNodes
        nElements = self.nElements
        Regions = zeros(nElements)
        Eids = zeros(nElements)

        #u    = zeros(nNodes)
        #v    = zeros(nNodes)
        #w    = zeros(nNodes)
        #Mach = zeros(nNodes)
        #rhoU = zeros(nNodes)
        #rhoV = zeros(nNodes)
        #rhoW = zeros(nNodes)

        if plotCentroidal:
            for eid, regioni in sorted(regions.iteritems()):
                Eids[eid - 1] = eid
                Regions[eid - 1] = regioni
            ###
            del regions
            cases[(ID, 'Region', 1, 'centroid', '%.0f')] = Regions
            cases[(ID, 'Eids', 1, 'centroid', '%.0f')] = Eids

        print "load.keys() = ", loads.keys()

        if 'Cp' in loads and plotNodal:
            cp = loads['Cp']
            Cp = zeros(nNodes)
            for nid, cpi in sorted(cp.iteritems()):
                Cp[nid - 1] = cpi
            ###
            #del loads['Cp']
            #del cp
            cases[(ID, 'Cp', 1, 'nodal', '%.3f')] = Cp

        if 'Mach' in loads and plotNodal:
            mach = loads['Mach']
            Mach = zeros(nNodes)
            for nid, machi in sorted(mach.iteritems()):
                Mach[nid - 1] = machi
            ###
            cases[(ID, 'Mach', 1, 'nodal', '%.3f')] = Mach

        if 'U' in loads and plotNodal:
            u = loads['U']
            U = zeros(nNodes)
            for nid, ui in sorted(u.iteritems()):
                U[nid - 1] = ui
            ###
            cases[(ID, 'U', 1, 'nodal', '%.3f')] = U

        if 'V' in loads and plotNodal:
            v = loads['V']
            V = zeros(nNodes)
            for nid, vi in sorted(v.iteritems()):
                V[nid - 1] = vi
            ###
            cases[(ID, 'V', 1, 'nodal', '%.3f')] = V

        if 'W' in loads and plotNodal:
            w = loads['W']
            W = zeros(nNodes)
            for nid, wi in sorted(w.iteritems()):
                W[nid - 1] = wi
            ###
            cases[(ID, 'W', 1, 'nodal', '%.3f')] = W

        #cases[(ID,'rhoU',   1,'nodal','%.3f')] = rhoU
        #cases[(ID,'rhoV',   1,'nodal','%.3f')] = rhoV
        #cases[(ID,'rhoW',   1,'nodal','%.3f')] = rhoW
        return cases
