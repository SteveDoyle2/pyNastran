"""
Defines the GUI IO file for SU2.
"""
from __future__ import print_function
from six import iteritems
from six.moves import range
import numpy as np

import vtk
from vtk import vtkTriangle, vtkQuad

from pyNastran.converters.su2.su2_reader import SU2Reader as SU2
from pyNastran.gui.gui_objects.gui_result import GuiResult
from pyNastran.gui.utils.vtk.vtk_utils import numpy_to_vtk_points


class SU2_IO(object):
    def __init__(self):
        pass

    def get_su2_wildcard_geometry_results_functions(self):
        data = ('SU2',
                'Stanford Unstructured (*.SU2)', self.load_su2_geometry,
                None, None)
        return data

    def load_su2_geometry(self, su2_filename, name='main', plot=True):
        #print("load_su2_geometry...")
        skip_reading = self._remove_old_geometry(su2_filename)
        if skip_reading:
            return

        model = SU2(log=self.log, debug=False)
        #self.model_type = model.model_type
        ndim, nodes, elements, regions = model.read_su2(su2_filename)

        nnodes = nodes.shape[0]
        nelements = 0
        for etype, elems in iteritems(elements):
            nsub_elements = elems.shape[0]
            if nsub_elements:
                nelements += nsub_elements
                #print('min of type = %s' % elems.min())
        assert nnodes > 0, nnodes
        assert nelements > 0, nelements

        self.nnodes = nnodes
        self.nelements = nelements

        self.log.info('nnodes=%s nelements=%s' % (self.nnodes, self.nelements))
        self.grid.Allocate(self.nelements, 1000)
        #self.gridResult.SetNumberOfComponents(self.nelements)

        #vectorReselt.SetNumberOfComponents(3)
        self.nid_map = {}

        assert nodes is not None
        nnodes = nodes.shape[0]
        if ndim == 3:
            xmax, ymax, zmax = nodes.max(axis=0)
            xmin, ymin, zmin = nodes.min(axis=0)
            self.log.info('xmax=%s xmin=%s' % (xmax, xmin))
            self.log.info('ymax=%s ymin=%s' % (ymax, ymin))
            self.log.info('zmax=%s zmin=%s' % (zmax, zmin))
            dim_max = max(xmax-xmin, ymax-ymin, zmax-zmin)
        elif ndim == 2:
            xmax, ymax = nodes.max(axis=0)
            xmin, ymin = nodes.min(axis=0)
            self.log.info('xmax=%s xmin=%s' % (xmax, xmin))
            self.log.info('ymax=%s ymin=%s' % (ymax, ymin))
            dim_max = max(xmax-xmin, ymax-ymin)

        self.create_global_axes(dim_max)

        if ndim == 2:
            nodes = np.hstack([nodes, np.zeros((nnodes, 1), dtype=nodes.dtype)])
        #else:
            # ndim=3

        points = numpy_to_vtk_points(nodes)


        #nelements = 0
        #elements = {
            #5 : tris,
            #9 : quads,
        #}
        #print('dict =', elements)
        for etype, elems in iteritems(elements):
            #print(etype, elems)
            if isinstance(elems, list):
                #print('continue')
                continue
            #print(type(elems))
            nsub_elements = elems.shape[0]
            if nsub_elements == 0:
                #print('continue')
                continue
            #print('eid_min =', elems.min())
            assert nsub_elements > 0, nsub_elements
            if etype == 5:
                for eid in range(nsub_elements):
                    elem = vtkTriangle()
                    node_ids = elems[eid, :]
                    elem.GetPointIds().SetId(0, node_ids[0])
                    elem.GetPointIds().SetId(1, node_ids[1])
                    elem.GetPointIds().SetId(2, node_ids[2])
                    self.grid.InsertNextCell(5, elem.GetPointIds())  #elem.GetCellType() = 5  # vtkTriangle
            elif etype == 9:
                for eid in range(nsub_elements):
                    elem = vtk.vtkQuad()
                    node_ids = elems[eid, :]
                    elem.GetPointIds().SetId(0, node_ids[0])
                    elem.GetPointIds().SetId(1, node_ids[1])
                    elem.GetPointIds().SetId(2, node_ids[2])
                    elem.GetPointIds().SetId(3, node_ids[3])
                    self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            else:
                raise NotImplementedError(etype)

        self.grid.SetPoints(points)
        self.grid.Modified()
        if hasattr(self.grid, 'Update'):
            self.grid.Update()
            self.log_info("updated grid")

        # loadSTLResults - regions/loads
        self.scalarBar.VisibilityOff()
        self.scalarBar.Modified()

        cases = {}
        self.isubcase_name_map = {}
        ID = 1

        form, cases = self._fill_su2_case(cases, ID, nelements, nnodes)
        self._finish_results_io2(form, cases)

    def _fill_su2_case(self, cases, ID, nelements, nnodes):
        """adds the sidebar results"""
        self.isubcase_name_map[ID] = ('SU2', '')

        #nelements = elements.shape[0]
        #nnodes = nodes.shape[0]
        icase = 0
        itime = 0
        eids = np.arange(1, nelements + 1, dtype='int32')
        nids = np.arange(1, nnodes + 1, dtype='int32')
        eid_res = GuiResult(ID, header='ElementID', title='ElementID',
                            location='centroid', scalar=eids)
        nid_res = GuiResult(ID, header='NodeID', title='NodeID',
                            location='node', scalar=nids)

        icase = 0
        itime = 0
        cases[icase] = (eid_res, (itime, 'ElementID'))
        cases[icase + 1] = (nid_res, (itime, 'NodeID'))
        form = [
            ('ElementID', icase, []),
            ('NodeID', icase + 1, []),
        ]
        return form, cases
