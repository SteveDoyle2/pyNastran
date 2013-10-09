import vtk
from vtk import vtkQuad
from panairGrid import PanairGrid


class PanairIO(object):
    def __init__(self):
        pass
#if __name__=='__main__':
#    lawgs = LaWGS('tmx1242.wgs')
#    lawgs.run()

    def load_panair_geometry(self, panairFileName, dirname):
        self.nidMap = {}

        #key = self.caseKeys[self.iCase]
        #case = self.resultCases[key]

        skipReading = self.removeOldGeometry(panairFileName)
        if skipReading:
            return

        model = PanairGrid(panairFileName, log=self.log, debug=self.debug)
        self.modelType = model.modelType
        model.readGrid()

        nodes, elements, regions = model.getPointsElementsRegions()
        #for nid,node in enumerate(nodes):
            #print "node[%s] = %s" %(nid,str(node))

        self.nNodes = len(nodes)
        self.nElements = len(elements)

        #print("nNodes = ",self.nNodes)
        print("nElements = ", self.nElements)

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
            for nid, node in sorted(nodes.iteritems()):
                points.InsertPoint(nid - 1, *point)
                self.gridResult.InsertNextValue(nid * fraction)
                #print str(element)

                #elem = vtk.vtkVertex()
                #elem.GetPointIds().SetId(0, i)
                #self.aQuadGrid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                #vectorResult.InsertTuple3(0, 0.0, 0.0, 1.0)

        for nid, node in enumerate(nodes):
            points.InsertPoint(nid, *node)
        #print "nid = ",nid

        for eid, element in enumerate(elements):
            (p1, p2, p3, p4) = element
            #print "element = ",element
            elem = vtkQuad()
            elem.GetPointIds().SetId(0, p1)
            elem.GetPointIds().SetId(1, p2)
            elem.GetPointIds().SetId(2, p3)
            elem.GetPointIds().SetId(3, p4)
            self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())

        print("eid = ", eid)
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

        self.resultCases = cases
        self.caseKeys = sorted(cases.keys())
        print "caseKeys = ",self.caseKeys
        #print "type(caseKeys) = ",type(self.caseKeys)
        self.iCase = -1
        self.nCases = len(self.resultCases) #- 1  # number of keys in dictionary
        if self.nCases > 1:
            self.nCases -= 1
        self.cycleResults()  # start at nCase=0

    def fillPanairGeometryCase(self, cases, ID, nodes, elements, regions, loads):
        #print "regions**** = ",regions
        #nNodes = self.nNodes
        #nElements = self.nElements

        print "is_centroidal=%s isNodal=%s" % (self.is_centroidal, self.is_nodal)
        assert self.is_centroidal!= self.is_nodal

        #result_names = ['Cp', 'Mach', 'U', 'V', 'W', 'E', 'rho',
                                      #'rhoU', 'rhoV', 'rhoW', 'rhoE']
        #if self.is_centroidal:
        #nelements, three = elements.shape
        #print regions
        cases[(ID, 'Region', 1, 'centroid', '%.0f')] = regions

        from numpy import zeros, array, cross, dot
        from numpy.linalg import det, norm
        # centroidal
        if self.is_centroidal:
            Xc = zeros(len(elements), 'float64')
            Yc = zeros(len(elements), 'float64')
            Zc = zeros(len(elements), 'float64')
            area = zeros(len(elements), 'float64')
            for i,element in enumerate(elements):
                p1, p2, p3, p4 = element
                P1 = array(nodes[p1])
                P2 = array(nodes[p2])
                P3 = array(nodes[p3])
                P4 = array(nodes[p4])
                a = P3 - P1
                b = P4 - P2
                A = 0.5 * norm(cross(a, b))
                #assert -1 > 0, 'A =%s' % str(A)
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
            # nodal
            Xn = zeros(len(nodes), 'float64')
            Yn = zeros(len(nodes), 'float64')
            Zn = zeros(len(nodes), 'float64')
            for i, node in enumerate(nodes):
                Xn[i] = node[0]
                Yn[i] = node[1]
                Zn[i] = node[2]
            cases[(ID, 'node_x', 1, 'nodal', '%.2f')] = Xn
            cases[(ID, 'node_y', 1, 'nodal', '%.2f')] = Yn
            cases[(ID, 'node_z', 1, 'nodal', '%.2f')] = Zn


        #elif self.is_nodal:
            #pass
            #print("load.keys() = ", loads.keys())
            #break
            #for key in result_names:
                #if key in loads:
                    #nodal_data = loads[key]
                    #cases[(ID, key, 1, 'nodal', '%.3f')] = nodal_data
        return cases

    def load_panair_results(self, panairFileName, dirname):
        #self.resultCases = {}
        pass

if __name__ == '__main__':
    print('')

    def removeOldGeometry(self):
        pass

    test = PanairIO()
    test.removeOldGeometry = removeOldGeometry

    #test.load_panair_geometry('SWB.INP','',True,True)
    test.load_panair_geometry('models/NAC6.INP', '', True, True)
