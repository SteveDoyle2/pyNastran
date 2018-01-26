from __future__ import print_function
from six.moves import range
import numpy as np
from numpy import amax, amin

import vtk
from vtk import vtkQuad

from pyNastran.converters.dev.openvsp.degen_geom import DegenGeom
from pyNastran.gui.gui_objects.gui_result import GuiResult
from pyNastran.gui.utils.vtk.vtk_utils import (
    create_vtk_cells_of_constant_element_type, numpy_to_vtk_points)


class DegenGeomIO(object):
    def __init__(self):
        pass

    def get_degen_geom_wildcard_geometry_results_functions(self):# pragma: no cover
        data = ('DegenGeom',
                'DegenGeom (*.csv)', self.load_degen_geom_geometry,
                #'Cart3d (*.triq)', self.load_cart3d_results,
                None, None
               )
        return data

    def load_degen_geom_geometry(self, csv_filename,
                                 name='main', plot=True):# pragma: no cover
        #key = self.case_keys[self.icase]
        #case = self.result_cases[key]

        skip_reading = self._remove_old_adb_geometry(csv_filename)
        if skip_reading:
            return

        model = DegenGeom(log=self.log, debug=False)
        self.model_type = 'vspaero'
        #self.model_type = model.model_type
        model.read_degen_geom(csv_filename)
        nodes = []
        elements = []
        inid = 0
        for name, comps in sorted(model.components.items()):
            self.log.debug('name = %r' % name)
            #print(comps)
            #print('------------')
            for comp in comps:
                self.log.info(comp)
                nnodes = comp.xyz.shape[0]
                nodes.append(comp.xyz)
                is_elem = np.linalg.norm(comp.elements, axis=1) > 0
                elements.append(comp.elements[is_elem] + inid)
                inid += nnodes

        if len(nodes) == 1:
            nodes = nodes[0]
            elements = elements[0]
        else:
            nodes = np.vstack(nodes)
            elements = np.vstack(elements)


        nnodes = nodes.shape[0]
        nelements = elements.shape[0]
        self.nnodes = nnodes
        self.nelements = nelements

        self.grid.Allocate(self.nelements, 1000)
        #self.gridResult.SetNumberOfComponents(self.nelements)

        #self.gridResult.Allocate(self.nnodes, 1000)
        #vectorReselt.SetNumberOfComponents(3)
        self.nid_map = {}

        assert nodes is not None

        #print("nxyz_nodes=%s" % nxyz_nodes)
        mmax = amax(nodes, axis=0)

        mmin = amin(nodes, axis=0)
        dim_max = (mmax - mmin).max()
        self.create_global_axes(dim_max)

        points = numpy_to_vtk_points(nodes)
        #self.log.info('nxyz_nodes=%s nwake_nodes=%s total=%s' % (
            #nnodes, nwake_nodes, nxyz_nodes + nwake_nodes))
        #self.log.info('nxyz_elements=%s nwake_elements=%s total=%s' % (
            #nxyz_elements, nwake_elements, nxyz_elements + nwake_elements))

        elements -= 1
        etype = 9 # vtkQuad().GetCellType()
        create_vtk_cells_of_constant_element_type(self.grid, elements, etype)

        self.grid.SetPoints(points)
        self.grid.Modified()
        if hasattr(self.grid, 'Update'):
            self.grid.Update()
        #self.log_info("updated grid")

        # load results - regions/loads
        self.scalarBar.VisibilityOn()
        self.scalarBar.Modified()

        #mach = model.machs[0]
        #alpha = model.alphas[0]
        #beta = model.betas[0]
        #note = ':  Mach=%.2f, alpha=%.1f, beta=%.1f' % (mach, alpha, beta)
        note = 'name=%s' % name
        self.isubcase_name_map = {1: ['OpenVSP%s' % note, '']}
        cases = {}
        ID = 1

        form, cases = self._fill_degen_geom_case(cases, ID, model, nnodes, nelements)
        self._finish_results_io2(form, cases)

    #def clear_adb(self):
        #pass

    #def load_adb_results(self, cart3d_filename):
        #raise NotImplementedError()

    def _fill_degen_geom_case(self, cases, ID, model, nnodes, nelements):  # pragma: no cover
        icase = 0
        itime = 0
        form = [
            ('ElementID', icase, []),
            ('NodeID', icase + 1, []),
        ]

        #form = ['Geometry', None, []]
        #form0 = form[2]
        formi = []
        form0 = form

        nodes = np.arange(nnodes, dtype='int32')
        elements = np.arange(nelements, dtype='int32')

        eid_res = GuiResult(0, header='ElementID', title='ElementID',
                            location='centroid', scalar=elements)
        nid_res = GuiResult(0, header='NodeID', title='NodeID',
                            location='node', scalar=nodes)
        cases[icase] = (eid_res, (itime, 'ElementID'))
        cases[icase + 1] = (nid_res, (itime, 'NodeID'))
        return form, cases
