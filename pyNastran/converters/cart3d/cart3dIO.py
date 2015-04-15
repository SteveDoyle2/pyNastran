from __future__ import print_function
from six import iteritems
from six.moves import range
from numpy import arange, mean, amax, amin

import vtk
from vtk import vtkTriangle

from pyNastran.converters.cart3d.cart3d_reader import Cart3DReader


class Cart3dIO(object):
    def __init__(self):
        pass

    def get_cart3d_wildcard_geometry_results_functions(self):
        data = ('Cart3d',
                'Cart3d (*.tri; *.triq)', self.load_cart3d_geometry,
                'Cart3d (*.triq)', self.load_cart3d_results)
        return data

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
            #self.gridResult.Reset()
            #self.gridResult.Modified()

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

    def load_cart3d_geometry(self, cart3d_filename, dirname, plot=True):
        #key = self.caseKeys[self.iCase]
        #case = self.resultCases[key]

        skipReading = self.removeOldGeometry(cart3d_filename)
        if skipReading:
            return

        model = Cart3DReader(log=self.log, debug=False)
        self.modelType = 'cart3d'
        #self.modelType = model.modelType
        (nodes, elements, regions, loads) = model.read_cart3d(cart3d_filename)
        self.nNodes = model.nPoints
        self.nElements = model.nElementsRead

        #print("nNodes = ",self.nNodes)
        #print("nElements = ", self.nElements)

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
                #print(str(element))

                #elem = vtk.vtkVertex()
                #elem.GetPointIds().SetId(0, i)
                #self.aQuadGrid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                #vectorResult.InsertTuple3(0, 0.0, 0.0, 1.0)

        assert nodes is not None
        nnodes = nodes.shape[0]

        nid = 0
        #print("nnodes=%s" % nnodes)
        mmax = amax(nodes, axis=0)
        mmin = amin(nodes, axis=0)
        dim_max = (mmax - mmin).max()
        self.update_axes_length(dim_max)

        for i in range(nnodes):
            points.InsertPoint(nid, nodes[i, :])
            nid += 1

        nelements = elements.shape[0]
        elements -= 1
        for eid in range(nelements):
            elem = vtkTriangle()
            node_ids = elements[eid, :]
            elem.GetPointIds().SetId(0, node_ids[0])
            elem.GetPointIds().SetId(1, node_ids[1])
            elem.GetPointIds().SetId(2, node_ids[2])
            self.grid.InsertNextCell(5, elem.GetPointIds())  #elem.GetCellType() = 5  # vtkTriangle

        self.grid.SetPoints(points)
        #self.grid2.SetPoints(points2)
        #self.grid.GetPointData().SetScalars(self.gridResult)
        #print(dir(self.grid) #.SetNumberOfComponents(0))
        #self.grid.GetCellData().SetNumberOfTuples(1);
        #self.grid.GetCellData().SetScalars(self.gridResult)
        self.grid.Modified()
        #self.grid2.Modified()
        if hasattr(self.grid, 'Update'):
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

        #print("nElements = ",nElements)
        form, cases = self._fill_cart3d_case(cases, ID, nodes, elements, regions, loads, model)
        #self._finish_results_io(cases)
        self._finish_results_io2(form, cases)

    def clear_cart3d(self):
        pass

    def load_cart3d_results(self, cart3d_filename, dirname):
        model = Cart3DReader(log=self.log, debug=False)

        self.load_cart3d_geometry(cart3d_filename, dirname)


    def _fill_cart3d_case(self, cases, ID, nodes, elements, regions, loads, model):
        print("is_centroidal=%s isNodal=%s" % (self.is_centroidal, self.is_nodal))
        assert self.is_centroidal != self.is_nodal

        result_names = ['Cp', 'Mach', 'U', 'V', 'W', 'E', 'rho',
                        'rhoU', 'rhoV', 'rhoW', 'rhoE']
        nelements = elements.shape[0]
        nnodes = nodes.shape[0]

        cases_new = []
        new = False
        is_normals = True

        results_form = []
        if self.is_centroidal:
            geometry_form = [
                ('Region', 0, []),
                ('ElementID', 1, []),
            ]

            eids = arange(1, nelements + 1)

            if new:
                cases_new[0] = (ID, regions, 'Region', 'centroid', '%i')
                cases_new[1] = (ID, eids, 'ElementID', 'centroid', '%i')
            else:
                cases[(ID, 0, 'Region', 1, 'centroid', '%i')] = regions
                cases[(ID, 1, 'ElementID', 1, 'centroid', '%i')] = eids

            i = 2
            #if 0:
                #from pyNastran.converters.cart3d.cart3d_to_quad import get_normal_groups
                #normal, normal_groups = get_normal_groups(nodes, elements)
                #groups = zeros(nelements, dtype='int32')
                #for igroup, normal_group in enumerate(normal_groups):
                    #if igroup % 2 == 0:
                        #continue
                    #for ni in normal_group:
                        #groups[ni] = igroup + 1
                #cases[(ID, 2, 'Quad Group', 1, 'centroid', '%i')] = groups
                #geometry_form.append('Quad Group', 2, [])


            if is_normals:
                geometry_form.append(('Normal X', i, []))
                geometry_form.append(('Normal Y', i + 1, []))
                geometry_form.append(('Normal Z', i + 2, []))

                cnormals = model.get_normals(nodes, elements, shift_nodes=False)
                cnnodes = cnormals.shape[0]
                assert cnnodes == nelements, len(cnnodes)

                if new:
                    cases_new[i] = (ID, cnormals[:, 0], 'Normal X', 'centroid', '%.3f')
                    cases_new[i + 1] = (ID, cnormals[:, 1], 'Normal Y', 'centroid', '%.3f')
                    cases_new[i + 2] = (ID, cnormals[:, 2], 'Normal Z', 'centroid', '%.3f')
                else:
                    cases[(ID, i, 'Normal X', 1, 'centroid', '%.3f')] = cnormals[:, 0]
                    cases[(ID, i + 1, 'Normal Y', 1, 'centroid', '%.3f')] = cnormals[:, 1]
                    cases[(ID, i + 2, 'Normal Z', 1, 'centroid', '%.3f')] = cnormals[:, 2]
                i += 3

            # these are actually nodal results, so we convert to the centroid
            # by averaging the data (e.g. the Cp data)
            for result_name in result_names:
                if result_name in loads:
                    nodal_data = loads[result_name]
                    n1 = elements[:, 0]
                    n2 = elements[:, 1]
                    n3 = elements[:, 2]
                    result = (nodal_data[n1] + nodal_data[n2] + nodal_data[n3]) / 3.0
                    if new:
                        cases_new[i] = (result, result_name, 1, 'centroid', '%.3f')
                    else:
                        cases[(ID, i, result_name, 1, 'centroid', '%.3f')] = result
                    results_form.append((result_name, i, []))
                    i += 1

        elif self.is_nodal:
            geometry_form = [
                #('Region', 0, []),
                ('NodeID', 0, []),
            ]
            i = 0
            nids = arange(1, nnodes+1)
            if new:
                #cases_new[0] = (regions, 'Region', 'centroid', '%i')
                cases_new[0] = (ID, nids, 'NodeID', 'node', '%i')
            else:
                cases[(ID, 0, 'NodeID', 1, 'node', '%i')] = nids
            i += 1

            if is_normals:
                geometry_form.append(('Normal X', 1, []))
                geometry_form.append(('Normal Y', 2, []))
                geometry_form.append(('Normal Z', 3, []))

                cnormals = model.get_normals(nodes, elements, shift_nodes=False)
                nnormals = model.get_normals_at_nodes(nodes, elements, cnormals, shift_nodes=False)

                if new:
                    cases_new[i] = (ID, nnormals[:, 0], 'Normal X', 'node', '%.3f')
                    cases_new[i + 1] = (ID, nnormals[:, 1], 'Normal Y', 'node', '%.3f')
                    cases_new[i + 2] = (ID, nnormals[:, 2], 'Normal Z', 'node', '%.3f')
                else:
                    cases[(ID, i, 'Normal X', 1, 'node', '%.3f')] = nnormals[:, 0]
                    cases[(ID, i + 1, 'Normal Y', 1, 'node', '%.3f')] = nnormals[:, 1]
                    cases[(ID, i + 2, 'Normal Z', 1, 'node', '%.3f')] = nnormals[:, 2]
                i += 3

            for result_name in result_names:
                if result_name in loads:
                    nodal_data = loads[result_name]
                    if new:
                        cases_new[i] = (result, result_name, 1, 'node', '%.3f')
                    else:
                        cases[(ID, i, result_name, 1, 'node', '%.3f')] = nodal_data
                    results_form.append((result_name, i, []))
                    i += 1

        form = [
            ('Geometry', None, geometry_form),
        ]
        if len(results_form):
            form.append(('Results', None, results_form))
        return form, cases
        #if new:
            #obj = Cart3dResultsObj(form, new_cases)
            #obj.get_result_index_by_form_key()
            #obj.get_result_by_index(1)
