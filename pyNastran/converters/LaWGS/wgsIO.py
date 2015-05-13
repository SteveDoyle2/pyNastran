from six import iteritems
import vtk
from vtk import vtkQuad
from pyNastran.converters.LaWGS.wgsReader import LaWGS


class LaWGS_IO(object):
    def __init__(self):
        pass
#if __name__=='__main__':
#    lawgs = LaWGS('tmx1242.wgs')
#    lawgs.run()

    def get_lawgs_wildcard_geometry_results_functions(self):
        data = ('LaWGS',
                'LaWGS (*.inp; *.wgs)', self.load_lawgs_geometry,
                None, None)
        return data

    def load_lawgs_geometry(self, lawgsFileName, dirname, plot=True):
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

        #print("nNodes = ",self.nNodes)
        #print("nElements = ", self.nElements)

        self.grid.Allocate(self.nElements, 1000)
        #self.gridResult.SetNumberOfComponents(self.nElements)
        self.grid2.Allocate(1, 1000)

        points = vtk.vtkPoints()
        points.SetNumberOfPoints(self.nNodes)
        #self.gridResult.Allocate(self.nNodes, 1000)
        #vectorReselt.SetNumberOfComponents(3)
        self.nidMap = {}
        #elem.SetNumberOfPoints(nNodes)
        if 0:
            fraction = 1. / self.nNodes  # so you can color the nodes by ID
            for nid, node in sorted(iteritems(nodes)):
                points.InsertPoint(nid - 1, *node)
                self.gridResult.InsertNextValue(nid * fraction)
                #print(str(element))

                #elem = vtk.vtkVertex()
                #elem.GetPointIds().SetId(0, i)
                #self.aQuadGrid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                #vectorResult.InsertTuple3(0, 0.0, 0.0, 1.0)

        for nid, node in enumerate(nodes):
            points.InsertPoint(nid, *node)
        print("nid = ", nid)

        for eid, element in enumerate(elements):
            (p1, p2, p3, p4) = element
            #print("element = ",element)
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
        #print(dir(self.grid) #.SetNumberOfComponents(0))
        #self.grid.GetCellData().SetNumberOfTuples(1);
        #self.grid.GetCellData().SetScalars(self.gridResult)
        self.grid.Modified()
        self.grid2.Modified()
        if hasattr(self.grid, 'Update'):
            self.grid.Update()
            self.grid2.Update()
            print("updated grid")

        return

        # loadCart3dResults - regions/loads
        self.TurnTextOn()
        self.scalarBar.VisibilityOn()
        self.scalarBar.Modified()

        self.iSubcaseNameMap = {1: ['LaWGS', '']}
        cases = {}
        ID = 1

        #print("nElements = %s" % nElements)
        #cases = self._fill_lawgs_case(cases, ID, regions, loads)
        self._finish_results_io(cases)

    def _fill_lawgs_case(self, cases, ID, regions, loads):
        pass

if __name__ == '__main__':  # pragma: no cover
    print('')

    def removeOldGeometry(self):
        pass
    test = LaWGS_IO()
    test.removeOldGeometry = removeOldGeometry

    test.load_lawgs_geometry('tmx1242.wgs', '', True, True)
