from __future__ import print_function
from six import iteritems
import vtk
from vtk import vtkQuad
from numpy import array, arange, cross
from pyNastran.converters.LaWGS.wgs_reader import LaWGS
from pyNastran.gui.gui_objects.gui_result import GuiResult


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

    def load_lawgs_geometry(self, lawgsFileName, dirname, name='main', plot=True):
        #key = self.case_keys[self.icase]
        #case = self.result_cases[key]

        skip_reading = self._remove_old_geometry(lawgsFileName)
        if skip_reading:
            return

        model = LaWGS(lawgsFileName)
        self.model_type = model.model_type
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
        #self. turn_text_on()
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
        region_res = GuiResult(ID, header='Region', title='Region',
                               location='centroid', scalar=regions)
        eid_res = GuiResult(ID, header='ElementID', title='ElementID',
                            location='centroid', scalar=eids)
        nid_res = GuiResult(ID, header='NodeID', title='NodeID',
                            location='node', scalar=nids)
        cases[icase] = (region_res, (ID, 'Region'))
        cases[icase + 1] = (eid_res, (ID, 'ElementID'))
        cases[icase + 2] = (nid_res, (ID, 'NodeID'))

        #cases[(ID, icase, 'Region', 1, 'centroid', '%i', '')] = regions
        #cases[(ID, icase + 1, 'ElementID', 1, 'centroid', '%i', '')] = eids
        #cases[(ID, icase + 2, 'NodeID', 1, 'node', '%i', '')] = nids

        #nnids = len(nids)
        neids = len(elements)

        a = nodes[elements[:, 2], :] - nodes[elements[:, 0], :]
        b = nodes[elements[:, 3], :] - nodes[elements[:, 1], :]
        normals = cross(a, b, axis=1)

        assert normals.shape[0] == neids, normals.shape
        assert normals.shape[1] == 3, normals.shape

        x_res = GuiResult(ID, header='X', title='X',
                          location='node', scalar=nodes[:, 0])
        y_res = GuiResult(ID, header='X', title='X',
                          location='node', scalar=nodes[:, 1])
        z_res = GuiResult(ID, header='X', title='X',
                          location='node', scalar=nodes[:, 2])

        nx_res = GuiResult(ID, header='NormalX', title='NormalX',
                           location='node', scalar=normals[:, 0])
        ny_res = GuiResult(ID, header='NormalY', title='NormalY',
                           location='node', scalar=normals[:, 1])
        nz_res = GuiResult(ID, header='NormalZ', title='NormalZ',
                           location='node', scalar=normals[:, 2])

        cases[icase + 3] = (x_res, (ID, 'X'))
        cases[icase + 4] = (y_res, (ID, 'Y'))
        cases[icase + 5] = (z_res, (ID, 'Z'))

        cases[icase + 6] = (nx_res, (ID, 'NormalX'))
        cases[icase + 7] = (ny_res, (ID, 'NormalY'))
        cases[icase + 8] = (nz_res, (ID, 'NormalZ'))

        #cases[(ID, icase + 3, 'X', 1, 'node', '%.3f', '')] = nodes[:, 0]
        #cases[(ID, icase + 4, 'Y', 1, 'node', '%.3f', '')] = nodes[:, 1]
        #cases[(ID, icase + 5, 'Z', 1, 'node', '%.3f', '')] = nodes[:, 2]

        #cases[(ID, icase + 6, 'NormalX', 1, 'centroid', '%.3f', '')] = normals[:, 0]
        #cases[(ID, icase + 7, 'NormalY', 1, 'centroid', '%.3f', '')] = normals[:, 1]
        #cases[(ID, icase + 8, 'NormalZ', 1, 'centroid', '%.3f', '')] = normals[:, 2]
        return geometry_form, cases

