from __future__ import print_function
from six import iteritems
from six.moves import range
from numpy import arange

import vtk
from vtk import vtkTriangle

from pyNastran.converters.stl.stl_reader import STLReader


class STL_IO(object):
    def __init__(self):
        pass

    def get_stl_wildcard_geometry_results_functions(self):
        data = ('STL',
                'STereoLithography (*.STL)', self.load_stl_geometry,
                None, None)
        return data


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

    def load_stl_geometry(self, stl_filename, dirname, plot=True):
        print("load_stl_geometry...")
        skipReading = self.removeOldGeometry(stl_filename)
        if skipReading:
            return

        model = STLReader(log=self.log, debug=False)
        #self.modelType = model.modelType
        model.read_stl(stl_filename)
        nodes = model.nodes
        elements = model.elements

        if self.is_centroidal:
            normals = model.get_normals(elements)
            areas = model.get_area(elements)
        elif self.is_nodal:
            normals = None
            areas = None
            #normals = model.get_normals_at_nodes(elements)

        self.nNodes, three = nodes.shape
        self.nElements, three = elements.shape

        print("nNodes = %s" % self.nNodes)
        print("nElements = %s" % self.nElements)

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
                #print str(element)

                #elem = vtk.vtkVertex()
                #elem.GetPointIds().SetId(0, i)
                #self.aQuadGrid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                #vectorResult.InsertTuple3(0, 0.0, 0.0, 1.0)

        assert nodes is not None
        nnodes = nodes.shape[0]
        xmax, ymax, zmax = nodes.max(axis=0)
        xmin, ymin, zmin = nodes.min(axis=0)
        self.log.info('xmax=%s xmin=%s' % (xmax, xmin))
        self.log.info('ymax=%s ymin=%s' % (ymax, ymin))
        self.log.info('zmax=%s zmin=%s' % (zmax, zmin))
        dim_max = max(xmax-xmin, ymax-ymin, zmax-zmin)
        self.update_axes_length(dim_max)


        nid = 0
        print("nnodes=%s" % nnodes)
        for i in range(nnodes):
            points.InsertPoint(nid, nodes[i, :])
            nid += 1

        nelements = elements.shape[0]
        #elements -= 1
        for eid in range(nelements):
            elem = vtkTriangle()
            node_ids = elements[eid, :]
            elem.GetPointIds().SetId(0, node_ids[0])
            elem.GetPointIds().SetId(1, node_ids[1])
            elem.GetPointIds().SetId(2, node_ids[2])
            self.grid.InsertNextCell(5, elem.GetPointIds())  #elem.GetCellType() = 5  # vtkTriangle

        self.grid.SetPoints(points)
        self.grid.Modified()
        if hasattr(self.grid, 'Update'):
            self.grid.Update()
            print("updated grid")

        # loadSTLResults - regions/loads
        self.TurnTextOn()
        self.scalarBar.VisibilityOff()
        self.scalarBar.Modified()

        cases = {}
        self.iSubcaseNameMap = {}
        ID = 1

        cases = self._fill_stl_case(cases, ID, elements, nodes, normals, areas)
        self._finish_results_io(cases)

    def _fill_stl_case(self, cases, ID, elements, nodes, normals, areas):
        self.iSubcaseNameMap[ID] = ('STL', '')
        #print("is_centroidal=%s isNodal=%s" % (self.is_centroidal, self.is_nodal))
        assert self.is_centroidal != self.is_nodal

        if self.is_centroidal:
            nelements, three = elements.shape
            #cases[(ID, 'Region', 1, 'centroid', '%i')] = regions
            cases[(ID, 'ElementID', 1, 'centroid', '%i')] = arange(1, nelements+1)
            cases[(ID, 'NormalX', 1, 'centroid', '%.3f')] = normals[:, 0]
            cases[(ID, 'NormalY', 1, 'centroid', '%.3f')] = normals[:, 1]
            cases[(ID, 'NormalZ', 1, 'centroid', '%.3f')] = normals[:, 2]
            cases[(ID, 'Area', 1, 'centroid', '%.4e')] = areas

        elif self.is_nodal:
            nnodes, three = nodes.shape
            #cases[(ID, 'Region', 1, 'centroid', '%i')] = regions
            cases[(ID, 'NodeID', 1, 'node', '%i')] = arange(1, nnodes+1)
            #cases[(ID, 'NormalX', 1, 'node', '%.3f')] = normals[:, 0]
            #cases[(ID, 'NormalY', 1, 'node', '%.3f')] = normals[:, 1]
            #cases[(ID, 'NormalZ', 1, 'node', '%.3f')] = normals[:, 2]
        return cases
