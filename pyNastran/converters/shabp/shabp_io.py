from six import iteritems
from numpy import zeros, array, cross, amax, amin
from numpy.linalg import norm

import vtk
from vtk import vtkQuad

from pyNastran.converters.shabp.shabp import SHABP

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
            #print "element = ",element
            elem = vtkQuad()
            elem.GetPointIds().SetId(0, p1)
            elem.GetPointIds().SetId(1, p2)
            elem.GetPointIds().SetId(2, p3)
            elem.GetPointIds().SetId(3, p4)
            self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())

        #print("eid = ", eid)
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
            ('centroidX', icase + 5, []),
            ('centroidY', icase + 6, []),
            ('centroidZ', icase + 7, []),

            ('nodeX', icase + 8, []),
            ('nodeY', icase + 9, []),
            ('nodeZ', icase + 10, []),
        ]

        geometry_form = [
            ('Component', icase, []),
            ('PatchID', icase + 1, []),
            ('Impact', icase + 2, []),
            ('Shadow', icase + 3, []),
            ('Area', icase + 4, []),
            ('Location', None, location_form),
        ]
        form = [
            ('Geometry', None, geometry_form),
        ]
        cases[(ID, icase, 'Component', 1, 'centroid', '%i', '')] = components
        cases[(ID, icase + 1, 'PatchID', 1, 'centroid', '%i', '')] = patches
        cases[(ID, icase + 2, 'Impact', 1, 'centroid', '%i', '')] = impact
        cases[(ID, icase + 3, 'Shadow', 1, 'centroid', '%i', '')] = shadow

        XYZc = zeros((len(elements), 3), dtype='float32')
        #Normal = zeros((len(elements),3), dtype='float32')
        area = zeros(len(elements), dtype='float32')

        for i, element in enumerate(elements):
            p1, p2, p3, p4 = element
            P1 = array(nodes[p1])
            P2 = array(nodes[p2])
            P3 = array(nodes[p3])
            P4 = array(nodes[p4])
            a = P3 - P1
            b = P4 - P2
            n = cross(a, b)
            nnorm = norm(n)
            #normal = n / nnorm
            A = 0.5 * nnorm

            XYZc[i, :] = (P1 + P2 + P3 + P4) / 4.0
            #Normal[i, :] = normal
            area[i] = A
        cases[(ID, icase + 4, 'Area', 1, 'centroid', '%.2f')] = area
        cases[(ID, icase + 5, 'centroidX', 1, 'centroid', '%.2f', '')] = XYZc[:, 0]
        cases[(ID, icase + 6, 'centroidY', 1, 'centroid', '%.2f', '')] = XYZc[:, 1]
        cases[(ID, icase + 7, 'centroidZ', 1, 'centroid', '%.2f', '')] = XYZc[:, 2]

        #cases[(ID, 'normalX', 1, 'centroid', '%.2f')] = Normal[:,0]
        #cases[(ID, 'normalY', 1, 'centroid', '%.2f')] = Normal[:,1]
        #cases[(ID, 'normalZ', 1, 'centroid', '%.2f')] = Normal[:,2]

        Xn = zeros(len(nodes), dtype='float32')
        Yn = zeros(len(nodes), dtype='float32')
        Zn = zeros(len(nodes), dtype='float32')
        for i, node in enumerate(nodes):
            Xn[i] = node[0]
            Yn[i] = node[1]
            Zn[i] = node[2]
        cases[(ID, icase + 8, 'nodeX', 1, 'node', '%.2f', '')] = Xn
        cases[(ID, icase + 9, 'nodeY', 1, 'node', '%.2f', '')] = Yn
        cases[(ID, icase + 10, 'nodeZ', 1, 'node', '%.2f', '')] = Zn
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

