from __future__ import print_function
from six import iteritems
import vtk
from vtk import vtkQuad
from numpy import array, arange, cross
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

        skip_reading = self.removeOldGeometry(lawgsFileName)
        if skip_reading:
            return

        model = LaWGS(lawgsFileName)
        self.modelType = model.modelType
        model.readLaWGS()

        nodes, elements, regions = model.get_points_elements_regions()
        self.nNodes = len(nodes)
        self.nElements = len(elements)

        nodes = array(nodes, dtype='float32')
        elements = array(elements, dtype='int32')

        #print("nNodes = ",self.nNodes)
        #print("nElements = ", self.nElements)

        self.grid.Allocate(self.nElements, 1000)
        #self.gridResult.SetNumberOfComponents(self.nElements)

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

        assert len(nodes) > 0, len(nodes)
        assert len(elements) > 0, len(elements)
        for nid, node in enumerate(nodes):
            points.InsertPoint(nid, *node)

        elem = vtkQuad()
        etype = elem.GetCellType()
        for eid, element in enumerate(elements):
            (p1, p2, p3, p4) = element
            elem = vtkQuad()
            pts = elem.GetPointIds()
            pts.SetId(0, p1)
            pts.SetId(1, p2)
            pts.SetId(2, p3)
            pts.SetId(3, p4)
            self.grid.InsertNextCell(etype, elem.GetPointIds())

        self.grid.SetPoints(points)
        #self.grid.GetPointData().SetScalars(self.gridResult)
        #print(dir(self.grid) #.SetNumberOfComponents(0))
        #self.grid.GetCellData().SetNumberOfTuples(1);
        #self.grid.GetCellData().SetScalars(self.gridResult)
        self.grid.Modified()
        if hasattr(self.grid, 'Update'):
            self.grid.Update()

        # loadCart3dResults - regions/loads
        #self.TurnTextOn()
        #self.scalarBar.VisibilityOn()
        #self.scalarBar.Modified()

        self.iSubcaseNameMap = {1: ['LaWGS', '']}
        cases = {}
        ID = 1

        #print("nElements = %s" % nElements)
        form, cases = self._fill_lawgs_case(cases, ID, nodes, elements, regions)
        self._finish_results_io2(form, cases)

    def _fill_lawgs_case(self, cases, ID, nodes, elements, regions):
        eids = arange(1, len(elements) + 1, dtype='int32')
        nids = arange(1, len(nodes) + 1, dtype='int32')
        regions = array(regions, dtype='int32')

        icase = 0
        geometry_form = [
            ('Region', icase, []),
            ('ElementID', icase + 1, []),
            ('NodeID', icase + 2, []),
            ('X', icase + 3, []),
            ('Y', icase + 4, []),
            ('Z', icase + 5, []),
            ('NormalX', icase + 6, []),
            ('NormalY', icase + 7, []),
            ('NormalZ', icase + 8, []),
        ]
        cases[(ID, icase, 'Region', 1, 'centroid', '%i', '')] = regions
        cases[(ID, icase + 1, 'ElementID', 1, 'centroid', '%i', '')] = eids
        cases[(ID, icase + 2, 'NodeID', 1, 'node', '%i', '')] = nids

        #nnids = len(nids)
        neids = len(elements)

        a = nodes[elements[:, 2], :] - nodes[elements[:, 0], :]
        b = nodes[elements[:, 3], :] - nodes[elements[:, 1], :]
        normals = cross(a, b, axis=1)

        assert normals.shape[0] == neids, normals.shape
        assert normals.shape[1] == 3, normals.shape

        cases[(ID, icase + 3, 'X', 1, 'node', '%.3f', '')] = nodes[:, 0]
        cases[(ID, icase + 4, 'Y', 1, 'node', '%.3f', '')] = nodes[:, 1]
        cases[(ID, icase + 5, 'Z', 1, 'node', '%.3f', '')] = nodes[:, 2]

        cases[(ID, icase + 6, 'NormalX', 1, 'centroid', '%.3f', '')] = normals[:, 0]
        cases[(ID, icase + 7, 'NormalY', 1, 'centroid', '%.3f', '')] = normals[:, 1]
        cases[(ID, icase + 8, 'NormalZ', 1, 'centroid', '%.3f', '')] = normals[:, 2]
        return geometry_form, cases

