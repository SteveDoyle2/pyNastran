from __future__ import print_function
from six import iteritems
from six.moves import range
from numpy import arange, mean, amax, amin, zeros, hstack

import vtk
from vtk import vtkTriangle, vtkParametricSpline

from pyNastran.converters.openvsp.adb_reader import ADB_Reader


class ADB_IO(object):
    def __init__(self):
        pass

    def get_openvsp_wildcard_geometry_results_functions(self):
        data = ('VSPAero',
                'VSPAero (*.adb)', self.load_vsp_aero_geometry,
                #'Cart3d (*.triq)', self.load_cart3d_results,
                None, None
               )
        return data

    def _remove_old_adb_geometry(self, fileName):
        self.eid_map = {}
        self.nid_map = {}
        if fileName is None:
            #self.emptyResult = vtk.vtkFloatArray()
            #self.vectorResult = vtk.vtkFloatArray()
            self.scalarBar.VisibilityOff()
            skip_reading = True
        else:
            self.turn_text_off()
            self.grid.Reset()
            #self.gridResult.Reset()
            #self.gridResult.Modified()

            self.result_cases = {}
            self.ncases = 0
            try:
                del self.case_keys
                del self.icase
                del self.iSubcaseNameMap
            except:
                # print("cant delete geo")
                pass

            #print(dir(self))
            skip_reading = False
        #self.scalarBar.VisibilityOff()
        self.scalarBar.Modified()
        return skip_reading

    def load_vsp_aero_geometry(self, adb_filename, dirname, name='main', plot=True):
        #key = self.case_keys[self.icase]
        #case = self.result_cases[key]

        skip_reading = self._remove_old_adb_geometry(adb_filename)
        if skip_reading:
            return

        if 0:
            plot_wakes = False
            #plot_wakes = False  # does this work right?
        else:
            plot_wakes = True # doesn't work perfectly

        model = ADB_Reader(log=self.log, debug=False)
        self.model_type = 'vspaero'
        #self.model_type = model.model_type
        (nodes, elements) = model.read_adb(adb_filename)
        nxyz_nodes = nodes.shape[0]
        nwake_nodes = model.wake_xyz.shape[0]

        nxyz_elements = model.tris.shape[0]
        nwake_elements = max(0, model.wake_xyz.shape[0] - 1)
        if plot_wakes:
            nnodes = nxyz_nodes + nwake_nodes
            nelements = nxyz_elements + nwake_elements
        else:
            nnodes = nxyz_nodes
            nelements = nxyz_elements
        self.nNodes = nnodes
        self.nElements = nelements

        self.grid.Allocate(self.nElements, 1000)
        #self.gridResult.SetNumberOfComponents(self.nElements)

        points = vtk.vtkPoints()
        points.SetNumberOfPoints(self.nNodes)
        #self.gridResult.Allocate(self.nNodes, 1000)
        #vectorReselt.SetNumberOfComponents(3)
        self.nid_map = {}

        assert nodes is not None

        nid = 0
        #print("nxyz_nodes=%s" % nxyz_nodes)
        mmax = amax(nodes, axis=0)
        mmin = amin(nodes, axis=0)
        dim_max = (mmax - mmin).max()
        self.create_global_axes(dim_max)

        for i in range(nxyz_nodes):
            points.InsertPoint(nid, nodes[i, :])
            nid += 1
        self.log.info('nxyz_nodes=%s nwake_nodes=%s total=%s' % (
            nxyz_nodes, nwake_nodes, nxyz_nodes + nwake_nodes))
        self.log.info('nxyz_elements=%s nwake_elements=%s total=%s' % (
            nxyz_elements, nwake_elements, nxyz_elements + nwake_elements))

        elements -= 1
        for eid in range(nxyz_elements):
            elem = vtkTriangle()
            node_ids = elements[eid, :]
            elem.GetPointIds().SetId(0, node_ids[0])
            elem.GetPointIds().SetId(1, node_ids[1])
            elem.GetPointIds().SetId(2, node_ids[2])
            self.grid.InsertNextCell(5, elem.GetPointIds())  #elem.GetCellType() = 5  # vtkTriangle

        if plot_wakes:
            for j in range(nwake_nodes):
                node = model.wake_xyz[j, :]
                points.InsertPoint(nid, node)
                nid += 1

            for wake in model.wake_elements:
                i0, i1 = wake
                #node_ids = range(i0, i1 + 1)
                #elem.SetPoints(points[i0:i1 + 1])
                print(i0, i1)
                for ii in range(i0+1, i1):
                    elem = vtk.vtkLine()
                    elem.GetPointIds().SetId(0, ii - 1)
                    elem.GetPointIds().SetId(1, ii)
                    self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                eid += 1
                #break
            #assert i1 == nelements, 'ii=%s nelements=%s' % (ii, nelements)

        self.grid.SetPoints(points)
        self.grid.Modified()
        if hasattr(self.grid, 'Update'):
            self.grid.Update()
        print("updated grid")

        # load results - regions/loads
        self. turn_text_on()
        self.scalarBar.VisibilityOn()
        self.scalarBar.Modified()

        mach = model.machs[0]
        alpha = model.alphas[0]
        beta = model.betas[0]
        note = ':  Mach=%.2f, alpha=%.1f, beta=%.1f' % (mach, alpha, beta)
        self.iSubcaseNameMap = {1: ['OpenVSP%s' % note, '']}
        cases = {}
        ID = 1

        form, cases = self._fill_adb_case(cases, ID, model, plot_wakes)
        self._finish_results_io2(form, cases)

    #def clear_adb(self):
        #pass

    #def load_adb_results(self, cart3d_filename, dirname):
        #raise NotImplementedError()


    def _fill_adb_case(self, cases, ID, model, plot_wakes=False):
        nxyz_nodes = model.nodes.shape[0]
        nxyz_elements = model.tris.shape[0]
        nwake_nodes = model.wake_xyz.shape[0]
        nwake_elements = max(0, model.wake_xyz.shape[0] - 1)

        if plot_wakes:
            nnodes = nxyz_nodes + nwake_nodes
            nelements = nxyz_elements + nwake_elements
            #assert nnodes == 408 + 1584, nnodes
            #assert nelements == 2223, nelements
        else:
            nnodes = nxyz_nodes
            nelements = nxyz_elements

        Cp = model.Cp
        surf_id = model.surf_id
        area = model.area

        cases_new = []
        new = False
        is_normals = False

        results_form = []
        if 1:
            geometry_form = [
                ('Region', 0, []),
                ('ElementID', 1, []),
                ('Area', 2, []),
            ]

            eids = arange(1, nelements + 1)

            if nwake_elements and plot_wakes:
                fzero_pad = zeros(nwake_elements, dtype='float32')
                izero_pad = zeros(nwake_elements, dtype='int32')
                surf_id = hstack([surf_id, izero_pad])
                area = hstack([area, fzero_pad])

            assert len(eids) == nelements, len(eids)
            assert len(surf_id) == nelements, len(surf_id)
            assert len(area) == nelements, len(area)
            if new:
                cases_new[0] = (ID, surf_id, 'Region', 'centroid', '%i', '')
                cases_new[1] = (ID, eids, 'ElementID', 'centroid', '%i', '')
                cases_new[2] = (ID, area, 'Area', 'centroid', '%.3f', '')
            else:
                # this one...
                cases[(ID, 0, 'Region', 1, 'centroid', '%i', '')] = surf_id
                cases[(ID, 1, 'ElementID', 1, 'centroid', '%i', '')] = eids
                cases[(ID, 2, 'Area', 1, 'centroid', '%.3f', '')] = area

            i = 3
            if is_normals:
                geometry_form.append(('Normal X', i, []))
                geometry_form.append(('Normal Y', i + 1, []))
                geometry_form.append(('Normal Z', i + 2, []))

                cnormals = model.get_normals(nodes, elements, shift_nodes=False)
                cnnodes = cnormals.shape[0]
                assert cnnodes == nelements, len(cnnodes)

                if new:
                    cases_new[i] = (ID, cnormals[:, 0], 'Normal X', 'centroid', '%.3f', '')
                    cases_new[i + 1] = (ID, cnormals[:, 1], 'Normal Y', 'centroid', '%.3f', '')
                    cases_new[i + 2] = (ID, cnormals[:, 2], 'Normal Z', 'centroid', '%.3f', '')
                else:
                    cases[(ID, i, 'Normal X', 1, 'centroid', '%.3f', '')] = cnormals[:, 0]
                    cases[(ID, i + 1, 'Normal Y', 1, 'centroid', '%.3f', '')] = cnormals[:, 1]
                    cases[(ID, i + 2, 'Normal Z', 1, 'centroid', '%.3f', '')] = cnormals[:, 2]
                i += 3

            if new:
                cases_new[i] = (Cp, 'Cp', 1, 'centroid', '%.3f', '')
            else:
                cases[(ID, 3, 'Cp', 1, 'centroid', '%.3f', '')] = Cp
            results_form.append(('Cp', i, []))
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
                cases_new[0] = (ID, nids, 'NodeID', 'node', '%i', '')
            else:
                cases[(ID, 0, 'NodeID', 1, 'node', '%i', '')] = nids
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
                    cases[(ID, i, 'Normal X', 1, 'node', '%.3f', '')] = nnormals[:, 0]
                    cases[(ID, i + 1, 'Normal Y', 1, 'node', '%.3f', '')] = nnormals[:, 1]
                    cases[(ID, i + 2, 'Normal Z', 1, 'node', '%.3f', '')] = nnormals[:, 2]
                i += 3

            # convert element data to nodal data
            if nwake_nodes and plot_wakes:
                fzero_pad = zeros(nwake_nodes, dtype='float32')
                izero_pad = zeros(nwake_nodes, dtype='int32')
                surf_id = hstack([surf_id, izero_pad])
                area = hstack([area, izero_pad])

            Cp_nodes = zeros(nnodes, dtype='float32')
            area_nodes = zeros(nnodes, dtype='float32')
            region_nodes = zeros(nnodes, dtype='int32')
            ncount = zeros(nnodes, dtype='int32')
            tris = model.tris
            for ei in range(nxyz_elements):
                n1, n2, n3 = tris[ei, :]
                cpi = Cp[ei]
                areai = area[ei]
                surf_idi = surf_id[ei]
                ncount[n1] += 1
                ncount[n2] += 1
                ncount[n3] += 1
                Cp_nodes[n1] += cpi
                Cp_nodes[n2] += cpi
                Cp_nodes[n3] += cpi
                area_nodes[n1] += areai
                area_nodes[n2] += areai
                area_nodes[n3] += areai
                region_nodes[n1] += surf_idi
                region_nodes[n2] += surf_idi
                region_nodes[n3] += surf_idi
            for ni in range(nxyz_nodes):
                ncounti = ncount[ni]
                if ncounti > 0:
                    Cp_nodes[ni] /= ncounti
                    area_nodes[ni] /= ncounti
                    region_nodes[ni] /= ncounti
            #print('min', tris.min())

            assert len(Cp_nodes) == nnodes, 'len(Cp)=%s nnodes=%s' % (len(Cp_nodes), nnodes)
            assert len(area_nodes) == nnodes, 'len(area)=%s nnodes=%s' % (len(area_nodes), nnodes)
            assert len(region_nodes) == nnodes, 'len(surf_id)=%s nnodes=%s' % (len(region_nodes), nnodes)
            if new:
                cases_new[i] = (Cp_nodes, 'Cp', 1, 'node', '%.3f', '')
                cases_new[i + 1] = (area_nodes, 'Area', 1, 'node', '%.3f', '')
                cases_new[i + 2] = (region_nodes, 'SurfaceID', 1, 'node', '%i', '')
            else:
                cases[(ID, i, 'Cp', 1, 'node', '%.3f', '')] = Cp_nodes
                cases[(ID, i + 1, 'Area', 1, 'node', '%.3f', '')] = area_nodes
                cases[(ID, i + 2, 'SurfaceID', 1, 'node', '%i', '')] = region_nodes
            results_form.append(('Cp', i, []))
            results_form.append(('Area', i  + 1, []))
            geometry_form.append(('SurfaceID', i + 2, []))
            i += 3

        form = [
            ('Geometry', None, geometry_form),
        ]
        print(form)
        if len(results_form):
            form.append(('Results', None, results_form))
        return form, cases
        #if new:
            #obj = Cart3dResultsObj(form, new_cases)
            #obj.get_result_index_by_form_key()
            #obj.get_result_by_index(1)
