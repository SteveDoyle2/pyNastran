from __future__ import print_function
from six import iteritems
from six.moves import range
import numpy as np
from numpy import arange, mean, amax, amin, zeros, hstack

import vtk
from vtk import vtkQuad, vtkParametricSpline

from pyNastran.converters.dev.openvsp.degen_geom import DegenGeom
from pyNastran.gui.gui_objects.gui_result import GuiResult


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

    def load_degen_geom_geometry(self, csv_filename, dirname,
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
        for name, comps in sorted(model.components.items()):
            print('name = %r' % name)
            #print(comp)
            print('------------')
            for comp in comps:
                nodes = comp.xyz
                elements = comp.elements
                nnodes = nodes.shape[0]
                nelements = elements.shape[0]

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

        for i in range(nnodes):
            points.InsertPoint(nid, nodes[i, :])
            nid += 1
        #self.log.info('nxyz_nodes=%s nwake_nodes=%s total=%s' % (
            #nnodes, nwake_nodes, nxyz_nodes + nwake_nodes))
        #self.log.info('nxyz_elements=%s nwake_elements=%s total=%s' % (
            #nxyz_elements, nwake_elements, nxyz_elements + nwake_elements))

        elements -= 1
        for eid in range(nelements):
            elem = vtkQuad()
            #assert elem.GetCellType() == 9, elem.GetCellType()
            node_ids = elements[eid, :]
            elem.GetPointIds().SetId(0, node_ids[0])
            elem.GetPointIds().SetId(1, node_ids[1])
            elem.GetPointIds().SetId(2, node_ids[2])
            elem.GetPointIds().SetId(3, node_ids[3])
            #elem.GetCellType() = 5  # vtkTriangle
            self.grid.InsertNextCell(9, elem.GetPointIds())

        self.grid.SetPoints(points)
        self.grid.Modified()
        if hasattr(self.grid, 'Update'):
            self.grid.Update()
        self.log_info("updated grid")

        # load results - regions/loads
        self.scalarBar.VisibilityOn()
        self.scalarBar.Modified()

        #mach = model.machs[0]
        #alpha = model.alphas[0]
        #beta = model.betas[0]
        #note = ':  Mach=%.2f, alpha=%.1f, beta=%.1f' % (mach, alpha, beta)
        note = 'name=%s' % name
        self.iSubcaseNameMap = {1: ['OpenVSP%s' % note, '']}
        cases = {}
        ID = 1

        form, cases = self._fill_degen_geom_case(cases, ID, model, nnodes, nelements)
        self._finish_results_io2(form, cases)

    #def clear_adb(self):
        #pass

    #def load_adb_results(self, cart3d_filename, dirname):
        #raise NotImplementedError()


    def _fill_degen_geom_case(self, cases, ID, model, nnodes, nelements):  # pragma: no cover
        icase = 0
        itime = 0
        form = [
            #('ElementID', icase, []),
            ('NodeID', icase, []),
            #('NodeID', icase + 1, []),
        ]

        #form = ['Geometry', None, []]
        #form0 = form[2]
        formi = []
        form0 = form

        nodes = np.arange(nnodes + 1, dtype='int32')
        elements = np.arange(nelements + 1, dtype='int32')

        eid_res = GuiResult(0, header='ElementID', title='ElementID',
                            location='centroid', scalar=elements)
        nid_res = GuiResult(0, header='NodeID', title='NodeID',
                            location='node', scalar=nodes)
        #cases[icase] = (eid_res, (itime, 'ElementID'))
        cases[icase] = (nid_res, (itime, 'NodeID'))
        #cases[icase + 1] = (nid_res, (itime, 'NodeID'))
        return form, cases
