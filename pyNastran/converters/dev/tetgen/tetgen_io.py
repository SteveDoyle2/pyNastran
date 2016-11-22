from __future__ import print_function
from six import iteritems
from six.moves import range
import os

import vtk
from vtk import vtkTriangle, vtkTetra

from pyNastran.converters.dev.tetgen.tetgen import Tetgen
from pyNastran.gui.gui_objects.gui_result import GuiResult


class TetgenIO(object):
    def __init__(self):
        pass

    def get_tetgen_wildcard_geometry_results_functions(self):
        data = ('Tetgen',
                'Tetgen (*.smesh)', self.load_tetgen_geometry,
                None, None)
        return data

    def load_tetgen_geometry(self, smesh_filename, dirname, name='main', plot=True):
        print("load_tetgen_geometry...")
        skip_reading = self.removeOldGeometry(smesh_filename)
        if skip_reading:
            return

        model = Tetgen(log=self.log, debug=False)

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

        self.nNodes = nodes.shape[0]
        ntris = 0
        ntets = 0
        if dimension_flag == 2:
            ntris = tris.shape[0]
        elif dimension_flag == 3:
            ntets = tets.shape[0]
        else:
            raise RuntimeError()
        self.nElements = ntris + ntets

        #print("nNodes = ",self.nNodes)
        #print("nElements = ", self.nElements)

        self.grid.Allocate(self.nElements, 1000)
        #self.gridResult.SetNumberOfComponents(self.nElements)

        points = vtk.vtkPoints()
        points.SetNumberOfPoints(self.nNodes)
        #self.gridResult.Allocate(self.nNodes, 1000)
        #vectorReselt.SetNumberOfComponents(3)
        self.nid_map = {}
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

        assert nodes is not None
        nnodes = nodes.shape[0]

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
        if hasattr(self.grid, 'Update'):
            self.grid.Update()

        # loadSTLResults - regions/loads
        self. turn_text_on()
        self.scalarBar.VisibilityOff()
        self.scalarBar.Modified()

        cases = {}
        ID = 1

        #cases = self._fill_tetgen_case(cases, ID, elements)
        self._finish_results_io(cases)

    def _fill_tetgen_case(self, cases, ID, elements):
        return cases

