#VTK_TRIANGLE = 5

from numpy import zeros, arange, mean

import vtk
from vtk import vtkTriangle

from pyNastran.converters.stl.stl_reader import STLReader


class STL_IO(object):
    def __init__(self):
        pass

    def _removeOldGeometry(self, fileName):
        # unused...
        self.eidMap = {}
        self.nidMap = {}
        if fileName is None:
            self.scalarBar.VisibilityOff()
            skipReading = True
        else:
            self.TurnTextOff()
            self.grid.Reset()
            self.grid2.Reset()

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

    def load_stl_geometry(self, stl_filename, dirname):
        print "load_stl_geometry..."
        skipReading = self.removeOldGeometry(stl_filename)
        if skipReading:
            return

        model = STLReader(log=self.log, debug=False)
        #self.modelType = model.modelType
        model.read_stl(stl_filename)
        nodes = model.nodes
        elements = model.elements
        self.nNodes, three = nodes.shape
        self.nElements, three = elements.shape

        print("nNodes = ",self.nNodes)
        print("nElements = ", self.nElements)

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
        #elements -= 1
        for eid in xrange(nelements):
            elem = vtkTriangle()
            node_ids = elements[eid, :]
            elem.GetPointIds().SetId(0, node_ids[0])
            elem.GetPointIds().SetId(1, node_ids[1])
            elem.GetPointIds().SetId(2, node_ids[2])
            self.grid.InsertNextCell(5, elem.GetPointIds())  #elem.GetCellType() = 5  # vtkTriangle

        self.grid.SetPoints(points)
        self.grid.Modified()
        self.grid.Update()
        print("updated grid")

        # loadSTLResults - regions/loads
        self.TurnTextOn()
        self.scalarBar.VisibilityOff()
        self.scalarBar.Modified()

        cases = {}
        ID = 1

        cases = self._fill_stl_case(cases, ID, elements)

        self.resultCases = cases
        self.caseKeys = sorted(cases.keys())
        #print "caseKeys = ",self.caseKeys
        #print "type(caseKeys) = ",type(self.caseKeys)
        self.nCases = min(0, len(self.resultCases) - 1)  # number of keys in dictionary
        self.iCase = 0 if self.nCases == 0 else -1
        self.cycleResults()  # start at nCase=0

    def _fill_stl_case(self, cases, ID, elements):
        return cases
        #print "regions**** = ",regions
        #nNodes = self.nNodes
        #nElements = self.nElements

        print "is_centroidal=%s isNodal=%s" % (self.is_centroidal, self.is_nodal)
        assert self.is_centroidal!= self.is_nodal

        if self.is_centroidal:
            nelements, three = elements.shape
            cases[(ID, 'Region', 1, 'centroid', '%.0f')] = regions
            cases[(ID, 'Eids', 1, 'centroid', '%.0f')] = arange(1, nelements+1)

        elif self.is_nodal:
            #print("load.keys() = ", loads.keys())
            for key in result_names:
                if key in loads:
                    nodal_data = loads[key]
                    cases[(ID, key, 1, 'nodal', '%.3f')] = nodal_data
        return cases
