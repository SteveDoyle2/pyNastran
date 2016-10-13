from __future__ import print_function
from six import iteritems
from six.moves import range

import os
from numpy import arange, mean, amax, amin, vstack, zeros, unique, where, sqrt

import vtk
from vtk import vtkLine, vtkTriangle, vtkQuad
from vtk.util.numpy_support import numpy_to_vtk

from pyNastran.gui.gui_objects.gui_result import GuiResult
#from pyNastran.gui.qt_files.result import Result
from pyNastran.converters.dev.abaqus.abaqus import Abaqus
#from pyNastran.converters.cart3d.cart3d_result import Cart3dGeometry, Cart3dResult

#from pyNastran.converters.cart3d.input_c3d_reader import read_input_c3d
#from pyNastran.converters.cart3d.input_cntl_reader import read_input_cntl


class AbaqusIO(object):
    def __init__(self):
        pass

    def get_abaqus_wildcard_geometry_results_functions(self):
        data = ('Abaqus',
                'Abaqus (*.inp)', self.load_abaqus_geometry,
                #'Abaqus (*.triq)', self.load_cart3d_results
                )
        return data

    #def _remove_old_geometry(self, geom_filename):
        #skip_reading = False
        #params_to_delete = (
            #'case_keys', 'icase', 'iSubcaseNameMap',
            #'result_cases', 'eid_map', 'nid_map'
        #)
        #if geom_filename is None or geom_filename is '':
            #skip_reading = True
            #return skip_reading
        #else:
            #self.turn_text_off()
            #self.grid.Reset()

            #self.result_cases = {}
            #self.ncases = 0
            #for param in params_to_delete:
                #if hasattr(self, param):  # TODO: is this correct???
                    #try:
                        #delattr(self, param)
                    #except AttributeError:
                        #print('param =', param, hasattr(self, param))

            #skip_reading = False
        ##self.scalarBar.VisibilityOff()
        #self.scalarBar.Modified()
        #return skip_reading

    #def _remove_old_cart3d_geometry(self, filename):
        ##return self._remove_old_geometry(filename)

        #self.eid_map = {}
        #self.nid_map = {}
        #if filename is None:
            ##self.emptyResult = vtk.vtkFloatArray()
            ##self.vectorResult = vtk.vtkFloatArray()
            #self.scalarBar.VisibilityOff()
            #skip_reading = True
        #else:
            #self.turn_text_off()
            #self.grid.Reset()
            ##self.gridResult.Reset()
            ##self.gridResult.Modified()

            #self.result_cases = {}
            #self.ncases = 0
            #try:
                #del self.case_keys
                #del self.icase
                #del self.iSubcaseNameMap
            #except:
                ## print("cant delete geo")
                #pass

            ##print(dir(self))
            #skip_reading = False
        ##self.scalarBar.VisibilityOff()
        #self.scalarBar.Modified()
        #return skip_reading

    def load_abaqus_geometry(self, abaqus_filename, dirname, name='main', plot=True):
        skip_reading = self._remove_old_cart3d_geometry(abaqus_filename)
        if skip_reading:
            return

        self.eid_map = {}
        self.nid_map = {}
        model = Abaqus(log=self.log, debug=False)
        self.model_type = 'abaqus'
        #self.model_type = model.model_type
        model.read_abaqus_inp(abaqus_filename)
        for part_name, part in model.parts:
            nids = part.nids - 1
            nodes = part.nodes
            break
        nodes = model.nodes
        #elements = model.elements
        #regions = model.regions
        #loads = model.loads

        nnodes = nodes.shape[0]
        n_r2d2 = 0
        n_cpe3 = 0
        n_cpe4 = 0
        n_cpe4r = 0
        n_coh2d4 = 0
        if part.r2d2 is not None:
            n_r2d2 = part.r2d2.shape[0]
        if part.cpe3 is not None:
            n_cpe3 = part.cpe3.shape[0]
        if part.cpe4 is not None:
            n_cpe4 = part.cpe4.shape[0]
        if part.cpe4r is not None:
            n_cpe4r = part.cpe4r.shape[0]
        if part.coh2d4 is not None:
            n_coh2d4 = part.coh2d4.shape[0]
        nelements = n_r2d2 + n_cpe3 + n_cpe4 + n_cpe4r + n_coh2d4

        self.nNodes = nnodes
        self.nElements = nelements

        self.grid.Allocate(self.nElements, 1000)

        points = vtk.vtkPoints()
        points.SetNumberOfPoints(self.nNodes)
        self.nid_map = {}

        assert nodes is not None
        nnodes = nodes.shape[0]

        mmax = amax(nodes, axis=0)
        mmin = amin(nodes, axis=0)
        dim_max = (mmax - mmin).max()
        self.create_global_axes(dim_max)

        data_type = vtk.VTK_FLOAT
        points_array = numpy_to_vtk(
            num_array=nodes,
            deep=True,
            array_type=data_type
        )
        points.SetData(points_array)

        if part.r2d2 is not None:
            n_r2d2 = part.r2d2.shape[0]
        if part.cpe3 is not None:
            n_cpe3 = part.cpe3.shape[0]
        if part.cpe4 is not None:
            n_cpe4 = part.cpe4.shape[0]
        if part.cpe4r is not None:
            n_cpe4r = part.cpe4r.shape[0]
        if part.coh2d4 is not None:
            n_coh2d4 = part.coh2d4.shape[0]

        if part.r2d2:
            part.r2d2[:, 1:] -= 1
            for eid, node_ids in part.r2d2:
                elem = vtkLine()
                elem.GetPointIds().SetId(0, node_ids[0])
                elem.GetPointIds().SetId(1, node_ids[1])
                self.grid.InsertNextCell(5, elem.GetPointIds())

        if part.cpe3:
            part.cpe3[:, 1:] -= 1
            for eid, node_ids in part.cpe3:
                elem = vtkTriangle()
                elem.GetPointIds().SetId(0, node_ids[0])
                elem.GetPointIds().SetId(1, node_ids[1])
                elem.GetPointIds().SetId(2, node_ids[3])
                self.grid.InsertNextCell(5, elem.GetPointIds())

        if part.cpe4:
            part.cpe4[:, 1:] -= 1
            for eid, node_ids in part.cpe4:
                elem = vtkQuad()
                elem.GetPointIds().SetId(0, node_ids[0])
                elem.GetPointIds().SetId(1, node_ids[1])
                elem.GetPointIds().SetId(2, node_ids[2])
                elem.GetPointIds().SetId(3, node_ids[3])
                self.grid.InsertNextCell(5, elem.GetPointIds())

        if part.cpe4r:
            part.cpe4r[:, 1:] -= 1
            for eid, node_ids in part.cpe4r:
                elem = vtkLine()
                elem.GetPointIds().SetId(0, node_ids[0])
                elem.GetPointIds().SetId(1, node_ids[1])
                elem.GetPointIds().SetId(2, node_ids[2])
                elem.GetPointIds().SetId(3, node_ids[3])
                self.grid.InsertNextCell(5, elem.GetPointIds())

        if part.coh2d4:
            part.coh2d4[:, 1:] -= 1
            for eid, node_ids in part.coh2d4:
                elem = vtkLine()
                elem.GetPointIds().SetId(0, node_ids[0])
                elem.GetPointIds().SetId(1, node_ids[1])
                elem.GetPointIds().SetId(2, node_ids[2])
                elem.GetPointIds().SetId(3, node_ids[3])
                self.grid.InsertNextCell(5, elem.GetPointIds())

        nelements = elements.shape[0]
        elements -= 1
        for eid in range(nelements):
            elem = vtkTriangle()
            node_ids = elements[eid, :]
            elem.GetPointIds().SetId(0, node_ids[0])
            elem.GetPointIds().SetId(1, node_ids[1])
            elem.GetPointIds().SetId(2, node_ids[2])
            self.grid.InsertNextCell(5, elem.GetPointIds())

        self.grid.SetPoints(points)
        self.grid.Modified()
        if hasattr(self.grid, 'Update'):
            self.grid.Update()

        # loadCart3dResults - regions/loads
        self. turn_text_on()
        self.scalarBar.VisibilityOn()
        self.scalarBar.Modified()

        note = ''
        self.iSubcaseNameMap = {1: ['Cart3d%s' % note, '']}
        cases = {}
        ID = 1
        form, cases, icase = self._fill_abaqus_case(cases, ID, nodes, elements, regions, model)
        self._fill_cart3d_results(cases, form, icase, ID, loads, model, mach)
        self._finish_results_io2(form, cases)

    def clear_abaqus(self):
        pass

    def load_abaqus_results(self, cart3d_filename, dirname):
        raise NotImplementedError()

    def _fill_abaqus_case(self, cases, ID, nodes, elements, regions, model):
        return [], {}, 0
        nelements = elements.shape[0]
        nnodes = nodes.shape[0]

        eids = arange(1, nelements + 1)
        nids = arange(1, nnodes + 1)
        cnormals = model.get_normals(shift_nodes=False)
        cnnodes = cnormals.shape[0]
        assert cnnodes == nelements, len(cnnodes)

        #print('nnodes =', nnodes)
        #print('nelements =', nelements)
        #print('regions.shape =', regions.shape)
        subcase_id = 0
        labels = ['NodeID', 'ElementID']
        #cart3d_geo = Cart3dGeometry(subcase_id, labels,
                                    #nids, eids, regions, cnormals,
                                    #uname='Cart3dGeometry')

        rho_res = GuiResult(ID, header='NodeID', title='NodeID',
                            location='node', scalar=node_ids)
        rho_res = GuiResult(ID, header='ElementID', title='ElementID',
                            location='centroid', scalar=element_ids)

        cases = {
            0 : (cart3d_geo, (0, 'NodeID')),
            1 : (cart3d_geo, (0, 'ElementID')),
            #2 : (cart3d_geo, (0, 'Region')),
            #3 : (cart3d_geo, (0, 'NormalX')),
            #4 : (cart3d_geo, (0, 'NormalY')),
            #5 : (cart3d_geo, (0, 'NormalZ')),
        }
        geometry_form = [
            ('NodeID', 0, []),
            ('ElementID', 1, []),
            #('Region', 2, []),
            #('Normal X', 3, []),
            #('Normal Y', 4, []),
            #('Normal Z', 5, []),
        ]
        form = [
            ('Geometry', None, geometry_form),
        ]
        icase = 6
        return form, cases, icase
