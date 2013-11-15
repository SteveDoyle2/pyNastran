import os
import vtk
from vtk import vtkTriangle
from pyNastran.converters.tetgen.tetgen_reader import TetgenReader


class TetgenIO(object):
    def __init__(self):
        pass

    def load_tetgen_geometry(self, smesh_filename, dirname):
        print "load_tetgen_geometry..."
        skipReading = self.removeOldGeometry(smesh_filename)
        if skipReading:
            return

        model = TetgenReader(log=self.log, debug=False)

        base_filename, ext = os.path.splitext(smesh_filename)
        node_filename = base_filename + '.node'
        model.read_tetgen(node_filename, smesh_filename)
        nodes = model.nodes
        tris = model.tri

        self.nNodes, three = nodes.shape
        self.nElements, three = tris.shape

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

        #elements -= 1
        for (n0, n1, n2) in tris:
            elem = vtkTriangle()
            #node_ids = elements[eid, :]
            elem.GetPointIds().SetId(0, n0)
            elem.GetPointIds().SetId(1, n1)
            elem.GetPointIds().SetId(2, n2)
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

        #cases = self._fill_stl_case(cases, ID, elements)

        self.resultCases = cases
        self.caseKeys = sorted(cases.keys())
        #print "caseKeys = ",self.caseKeys
        #print "type(caseKeys) = ",type(self.caseKeys)
        self.nCases = min(0, len(self.resultCases) - 1)  # number of keys in dictionary
        self.iCase = 0 if self.nCases == 0 else -1
        self.cycleResults()  # start at nCase=0

    def _fill_stl_case(self, cases, ID, elements):
        return cases
    def read_tetgen():
        mesh = Tetgen_Reader()
        mesh.read_tetgen(node_filename, smesh_filename)


