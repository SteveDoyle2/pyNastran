#VTK_TRIANGLE = 5

from numpy import zeros, arange, mean

import vtk
from vtk import vtkTriangle

from pyNastran.converters.cart3d.cart3d_reader import Cart3DReader


class Cart3dIO(object):
    def __init__(self):
        pass

    def removeOldGeometry(self, fileName):
        self.eidMap = {}
        self.nidMap = {}
        if fileName is None:
            #self.emptyResult = vtk.vtkFloatArray()
            #self.vectorResult = vtk.vtkFloatArray()
            self.scalarBar.VisibilityOff()
            skipReading = True
        else:
            self.TurnTextOff()
            self.grid.Reset()
            self.grid2.Reset()
            #print(dir(self.grid2))
            #self.grid2.VisibilityOff()
            self.gridResult.Reset()
            self.gridResult.Modified()

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
        #self.scalarBar.VisibilityOff()
        self.scalarBar.Modified()
        return skipReading

    def load_cart3d_geometry(self, cart3dFileName, dirname):
        #key = self.caseKeys[self.iCase]
        #case = self.resultCases[key]

        skipReading = self.removeOldGeometry(cart3dFileName)
        if skipReading:
            return

        model = Cart3DReader(log=self.log, debug=False)
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

        assert loads is not None
        if 'Mach' in loads:
            avgMach = mean(loads['Mach'])
            note = ':  avg(Mach)=%g' % avgMach
        else:
            note = ''
        self.iSubcaseNameMap = {1: ['Cart3d%s' % note, '']}
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
        #nNodes = self.nNodes
        #nElements = self.nElements

        print "is_centroidal=%s isNodal=%s" % (self.is_centroidal, self.is_nodal)
        assert self.is_centroidal!= self.is_nodal
        
        result_names = ['Cp', 'Mach', 'U', 'V', 'W', 'E', 'rho',
                                      'rhoU', 'rhoV', 'rhoW', 'rhoE']
        if self.is_centroidal:
            nelements, three = elements.shape
            cases[(ID, 'Region', 1, 'centroid', '%.0f')] = regions
            cases[(ID, 'Eids', 1, 'centroid', '%.0f')] = arange(1, nelements+1)
            
            #print("load.keys() = ", loads.keys())
            #print("type(loads)", type(loads))
            for key in result_names:
                if key in loads:
                    nodal_data = loads[key]
                    n1 = elements[:, 0]
                    n2 = elements[:, 1]
                    n3 = elements[:, 2]
                    elemental_result = (nodal_data[n1] + nodal_data[n2] + nodal_data[n3])/3.0
                    cases[(ID, key, 1, 'centroid', '%.3f')] = elemental_result

        elif self.is_nodal:
            #print("load.keys() = ", loads.keys())
            for key in result_names:
                if key in loads:
                    nodal_data = loads[key]
                    cases[(ID, key, 1, 'nodal', '%.3f')] = nodal_data
        return cases
