import vtk
from vtk import vtkQuad
from wgsReader import LaWGS


class LaWGS_IO(object):
    def __init__(self):
        pass
#if __name__=='__main__':
#    lawgs = LaWGS('tmx1242.wgs')
#    lawgs.run()

    def loadLaWGSGeometry(self, lawgsFileName, dirname, isNodal, isCentroidal):
        self.isNodal = isNodal
        self.isCentroidal = isCentroidal
        #key = self.caseKeys[self.iCase]
        #case = self.resultCases[key]

        skipReading = self.removeOldGeometry(lawgsFileName)
        if skipReading:
            return

        model = LaWGS(lawgsFileName)
        self.modelType = model.modelType
        model.readLaWGS()

        nodes, elements = model.getPointsElements()
        self.nNodes = len(nodes)
        self.nElements = len(elements)

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

        for nid, node in enumerate(nodes):
            points.InsertPoint(nid, *node)
        print "nid = ", nid

        for eid, element in enumerate(elements):
            (p1, p2, p3, p4) = element
            #print "element = ",element
            elem = vtkQuad()
            elem.GetPointIds().SetId(0, p1)
            elem.GetPointIds().SetId(1, p2)
            elem.GetPointIds().SetId(2, p3)
            elem.GetPointIds().SetId(3, p4)
            self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
        ###
        print "eid = ", eid
        self.grid.SetPoints(points)
        #self.grid2.SetPoints(points2)
        #self.grid.GetPointData().SetScalars(self.gridResult)
        #print dir(self.grid) #.SetNumberOfComponents(0)
        #self.grid.GetCellData().SetNumberOfTuples(1);
        #self.grid.GetCellData().SetScalars(self.gridResult)
        self.grid.Modified()
        self.grid2.Modified()
        self.grid.Update()
        self.grid2.Update()
        print "updated grid"

        return

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

if __name__ == '__main__':
    print ''

    def removeOldGeometry(self):
        pass
    test = LaWGS_IO()
    test.removeOldGeometry = removeOldGeometry

    test.loadLaWGSGeometry('tmx1242.wgs', '', True, True)
