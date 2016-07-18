from six import iteritems
import numpy as np
from numpy import zeros, array, cross, amax, amin
from numpy.linalg import norm

import vtk
from vtk import vtkQuad

from pyNastran.converters.shabp.shabp import SHABP
from pyNastran.gui.gui_objects.gui_result import GuiResult

is_shabp = True


class ShabpIO(object):
    def __init__(self):
        pass

    def get_shabp_wildcard_geometry_results_functions(self):
        data = ('S/HABP',
                'Shabp (*.geo; *.mk5; *.inp)', self.load_shabp_geometry,
                'Shabp (*.out)', self.load_shabp_results)
        return data

    def load_shabp_geometry(self, shabpFilename, dirname, name='main', plot=True):
        self.nid_map = {}

        #key = self.case_keys[self.icase]
        #case = self.result_cases[key]

        skip_reading = self._remove_old_geometry(shabpFilename)
        if skip_reading:
            return

        self.model = SHABP(log=self.log, debug=self.debug)
        self.model_type = 'shabp' # model.model_type
        self.model.read_shabp(shabpFilename)

        nodes, elements, patches, components, impact, shadow = self.model.getPointsElementsRegions()
        #for nid,node in enumerate(nodes):
            #print "node[%s] = %s" %(nid,str(node))

        self.nNodes = len(nodes)
        self.nElements = len(elements)

        #print("nNodes = ",self.nNodes)
        #print("nElements = ", self.nElements)

        self.grid.Allocate(self.nElements, 1000)
        #self.gridResult.SetNumberOfComponents(self.nElements)

        points = vtk.vtkPoints()
        points.SetNumberOfPoints(self.nNodes)
        #self.gridResult.Allocate(self.nNodes, 1000)
        #vectorReselt.SetNumberOfComponents(3)
        #elem.SetNumberOfPoints(nNodes)
        if 0:
            fraction = 1. / nNodes  # so you can color the nodes by ID
            for nid, node in sorted(iteritems(nodes)):
                points.InsertPoint(nid - 1, *node)
                self.gridResult.InsertNextValue(nid * fraction)
                #print str(element)

                #elem = vtk.vtkVertex()
                #elem.GetPointIds().SetId(0, i)
                #self.aQuadGrid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                #vectorResult.InsertTuple3(0, 0.0, 0.0, 1.0)

        assert len(nodes) > 0
        mmax = amax(nodes, axis=0)
        mmin = amin(nodes, axis=0)
        dim_max = (mmax - mmin).max()
        self.create_global_axes(dim_max)
        for nid, node in enumerate(nodes):
            points.InsertPoint(nid, *node)

        assert len(elements) > 0
        for eid, element in enumerate(elements):
            (p1, p2, p3, p4) = element
            elem = vtkQuad()
            elem.GetPointIds().SetId(0, p1)
            elem.GetPointIds().SetId(1, p2)
            elem.GetPointIds().SetId(2, p3)
            elem.GetPointIds().SetId(3, p4)
            self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())

        self.grid.SetPoints(points)
        #self.grid.GetPointData().SetScalars(self.gridResult)
        #print dir(self.grid) #.SetNumberOfComponents(0)
        #self.grid.GetCellData().SetNumberOfTuples(1);
        #self.grid.GetCellData().SetScalars(self.gridResult)
        self.grid.Modified()
        if hasattr(self.grid, 'Update'):
            self.grid.Update()

        # loadCart3dResults - regions/loads
        self. turn_text_on()
        self.scalarBar.VisibilityOn()
        self.scalarBar.Modified()

        self.iSubcaseNameMap = {1: ['S/HABP', '']}
        cases = {}
        ID = 1

        self.log.debug("nNodes=%i nElements=%i" % (self.nNodes, self.nElements))
        form, cases = self._fill_shabp_geometry_case(cases, ID, nodes, elements, patches, components, impact, shadow)
        self._finish_results_io2(form, cases)

    def clear_shabp(self):
        del self.elements
        del self.model

    def _fill_shabp_geometry_case(self, cases, ID, nodes, elements, patches, components, impact, shadow):
        self.elements = elements

        icase = 0
        location_form = [
            ('CentroidX', icase + 5, []),
            ('CentroidY', icase + 6, []),
            ('CentroidZ', icase + 7, []),

            ('NodeX', icase + 8, []),
            ('NodeY', icase + 9, []),
            ('NodeZ', icase + 10, []),
        ]
        normal_form = [
            ('NormalX', icase + 11, []),
            ('NormalY', icase + 12, []),
            ('NormalZ', icase + 13, []),
        ]

        geometry_form = [
            ('Component', icase, []),
            ('PatchID', icase + 1, []),
            ('Impact', icase + 2, []),
            ('Shadow', icase + 3, []),
            ('Area', icase + 4, []),
            ('Location', None, location_form),
            ('Normal', None, normal_form),
        ]
        form = [
            ('Geometry', None, geometry_form),
        ]
        itime = 0
        ID = 0
        components_res = GuiResult(ID, header='Component', title='Component',
                                   location='centroid', scalar=components)
        patch_res = GuiResult(ID, header='PatchID', title='PatchID',
                              location='centroid', scalar=patches)
        impact_res = GuiResult(ID, header='Impact', title='Impact',
                               location='centroid', scalar=impact)
        shadow_res = GuiResult(ID, header='Shadow', title='Shadow',
                               location='centroid', scalar=shadow)

        cases[icase] = (components_res, (itime, 'Component'))
        cases[icase+1] = (patch_res, (itime, 'PatchID'))
        cases[icase+2] = (impact_res, (itime, 'Impact'))
        cases[icase+3] = (shadow_res, (itime, 'Shadow'))

        XYZc = zeros((len(elements), 3), dtype='float32')
        Normal = zeros((len(elements),3), dtype='float32')
        area = zeros(len(elements), dtype='float32')

        if 0:
            elements = np.array(elements, dtype='int32')
            p1 = nodes[elements[:, 0]]
            p2 = nodes[elements[:, 1]]
            p3 = nodes[elements[:, 2]]
            p4 = nodes[elements[:, 3]]
            a = p3 - p1
            b = p4 - p2
            n = np.cross(a, b)
            ni = np.linalg.norm(n, axis=1)
            assert len(ni) == n.shape[0]
            i = np.where(ni != 0.0)[0]
            #n[i] /= ni[:, i]
            area[i] = 0.5 * ni[i]
            assert p1.shape == n.shape, n.shape

        for i, element in enumerate(elements):
            n1, n2, n3, n4 = element
            p1 = nodes[n1, :]
            p2 = nodes[n2, :]
            p3 = nodes[n3, :]
            p4 = nodes[n4, :]
            a = p3 - p1
            b = p4 - p2
            n = cross(a, b)
            nnorm = norm(n)

            XYZc[i, :] = (p1 + p2 + p3 + p4) / 4.0
            if nnorm == 0.:
                print('p1=%s p2=%s p3=%s p4=%s; area=0' % (p1, p2, p3, p4))
                continue
            normal = n / nnorm
            A = 0.5 * nnorm

            Normal[i, :] = normal
            area[i] = A

        area_res = GuiResult(ID, header='Area', title='Area',
                             location='centroid', scalar=area) # data_format='%.2f

        cenx_res = GuiResult(ID, header='CentroidX', title='CentroidX',
                           location='centroid', scalar=XYZc[:,0]) # data_format='%.2f
        ceny_res = GuiResult(ID, header='CentroidY', title='CentroidY',
                             location='centroid', scalar=XYZc[:,1]) # data_format='%.2f
        cenz_res = GuiResult(ID, header='CentroidZ', title='CentroidZ',
                             location='centroid', scalar=XYZc[:,2]) # data_format='%.2f

        nodex_res = GuiResult(ID, header='NodeX', title='NodeX',
                           location='node', scalar=nodes[:, 0]) # data_format='%.2f
        nodey_res = GuiResult(ID, header='NodeY', title='NodeY',
                           location='node', scalar=nodes[:, 1]) # data_format='%.2f
        nodez_res = GuiResult(ID, header='NodeZ', title='NodeZ',
                           location='node', scalar=nodes[:, 2]) # data_format='%.2f

        nx_res = GuiResult(ID, header='NormalX', title='NormalX',
                           location='centroid', scalar=Normal[:,0]) # data_format='%.2f
        ny_res = GuiResult(ID, header='NormalY', title='NormalY',
                           location='centroid', scalar=Normal[:,1]) # data_format='%.2f
        nz_res = GuiResult(ID, header='NormalZ', title='NormalZ',
                           location='centroid', scalar=Normal[:,2]) # data_format='%.2f

        cases[icase+4] = (area_res, (itime, 'Area'))
        cases[icase+5] = (cenx_res, (itime, 'CentroidX'))
        cases[icase+6] = (ceny_res, (itime, 'CentroidY'))
        cases[icase+7] = (cenz_res, (itime, 'CentroidZ'))
        cases[icase+8] = (nodex_res, (itime, 'NodeX'))
        cases[icase+9] = (nodey_res, (itime, 'NodeY'))
        cases[icase+10] = (nodez_res, (itime, 'NodeZ'))

        cases[icase+11] = (nx_res, (itime, 'NormalX'))
        cases[icase+12] = (ny_res, (itime, 'NormalY'))
        cases[icase+13] = (nz_res, (itime, 'NormalZ'))
        return form, cases

    def load_shabp_results(self, shabp_filename, dirname):
        Cpd, deltad = self.model.read_shabp_out(shabp_filename)

        cases = self.result_cases
        icase = len(cases)
        mach_results = []
        form = self.form
        form.append(('Results', None, mach_results))
        #self.result_cases = {}
        mach_forms = {}
        for case_id, Cp in sorted(iteritems(Cpd)):
            Cp = Cpd[case_id]
            #delta = deltad[case_id]

            mach, alpha, beta = self.model.shabp_cases[case_id]
            #name = 'Mach=%g Alpha=%g' % (mach, alpha)
            name = 'Mach=%g Alpha=%g' % (mach, alpha)
            cases[(name, icase, 'Cp', 1, 'centroid', '%.3f', '')] = Cp
            cp_form = [
                ('Cp', icase, [])
            ]
            mach_forms[mach].append(('Cp', None, cp_form))
            #self.result_cases[(name, 'delta', 1, 'centroid', '%.3f')] = delta

        for mach, mach_form in sorted(iteritems(mach_forms)):
            mach_results.append(mach_form)
        self._finish_results_io2(form, cases)

