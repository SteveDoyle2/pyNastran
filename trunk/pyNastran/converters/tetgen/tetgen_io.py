from six.moves import range
import os

import vtk
from vtk import vtkTriangle, vtkTetra

from pyNastran.converters.tetgen.tetgen_reader import TetgenReader


class TetgenIO(object):
    def __init__(self):
        pass

    def load_tetgen_geometry(self, smesh_filename, dirname):
        print("load_tetgen_geometry...")
        skipReading = self.removeOldGeometry(smesh_filename)
        if skipReading:
            return

        model = TetgenReader(log=self.log, debug=False)

        base_filename, ext = os.path.splitext(smesh_filename)
        node_filename = base_filename + '.node'
        ele_filename = base_filename + '.ele'
        if '.smesh' == ext:
            dimension_flag = 2
        elif '.ele' == ext:
            dimension_flag = 3
        else:
            raise RuntimeError('unsupported extension.  Use "smesh" or "ele".')
        model.read_tetgen(node_filename, smesh_filename, ele_filename, dimension_flag)
        nodes = model.nodes
        tris = model.tri
        tets = model.tet

        self.nNodes, three = nodes.shape
        ntris = 0
        ntets = 0
        if dimension_flag == 2:
            ntris, three = tris.shape
        elif dimension_flag == 3:
            ntets, four = tets.shape
        else:
            raise RuntimeError()
        self.nElements = ntris + ntets

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
                #print(str(element))

                #elem = vtk.vtkVertex()
                #elem.GetPointIds().SetId(0, i)
                #self.aQuadGrid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                #vectorResult.InsertTuple3(0, 0.0, 0.0, 1.0)

        assert nodes is not None
        nnodes, three = nodes.shape

        nid = 0
        print("nnodes=%s" % nnodes)
        for i in range(nnodes):
            points.InsertPoint(nid, nodes[i, :])
            nid += 1

        #elements -= 1
        if dimension_flag == 2:
            for (n0, n1, n2) in tris:
                elem = vtkTriangle()
                #node_ids = elements[eid, :]
                elem.GetPointIds().SetId(0, n0)
                elem.GetPointIds().SetId(1, n1)
                elem.GetPointIds().SetId(2, n2)
                self.grid.InsertNextCell(5, elem.GetPointIds())  #elem.GetCellType() = 5  # vtkTriangle
        elif dimension_flag == 3:
            for (n0, n1, n2, n3) in tets:
                elem = vtkTetra()
                assert elem.GetCellType() == 10, elem.GetCellType()
                elem.GetPointIds().SetId(0, n0)
                elem.GetPointIds().SetId(1, n1)
                elem.GetPointIds().SetId(2, n2)
                elem.GetPointIds().SetId(3, n3)
                self.grid.InsertNextCell(10, elem.GetPointIds())  #elem.GetCellType() = 5  # vtkTriangle
        else:
            raise RuntimeError()

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

        #cases = self._fill_tetgen_case(cases, ID, elements)
        self.finish_nastran_io(cases)

    def _fill_tetgen_case(self, cases, ID, elements):
        return cases

