from collections import OrderedDict
from numpy import arange, amax, amin

import vtk
from vtk import vtkQuad, vtkLine

from pyNastran.gui.gui_objects.gui_result import GuiResult
from pyNastran.converters.avl.avl import read_avl


class AVL_IO:
    def __init__(self, gui):
        self.gui = gui

    def get_avl_wildcard_geometry_results_functions(self):  # pragma: no cover
        data = ('AVL',
                'AVL (*.avl)', self.load_avl_geometry,
                #'Cart3d (*.triq)', self.load_cart3d_results,
                None, None
               )
        return data

    def load_avl_geometry(self, avl_filename,
                          name='main', plot=True):
        model_name = name
        #key = self.case_keys[self.icase]
        #case = self.result_cases[key]

        #skip_reading = self._remove_old_adb_geometry(adb_filename)
        #if skip_reading:
            #return

        log = self.gui.log
        model = read_avl(avl_filename, log=log, debug=False)
        self.gui.model_type = 'avl'
        nodes, elements, line_elements, surfaces = model.get_nodes_elements()
        #self.model_type = model.model_type

        nxyz_nodes = nodes.shape[0]
        nquad_elements = elements.shape[0]
        nline_elements = 0
        if line_elements:
            nline_elements = line_elements.shape[0]
            assert nline_elements > 0, nline_elements

        nnodes = nxyz_nodes
        nelements = nquad_elements + nline_elements
        self.gui.nnodes = nnodes
        self.gui.nelements = nelements

        grid = self.gui.grid
        grid.Allocate(self.gui.nelements, 1000)

        points = vtk.vtkPoints()
        points.SetNumberOfPoints(self.gui.nnodes)
        #vectorReselt.SetNumberOfComponents(3)
        self.gui.nid_map = {}

        assert nodes is not None

        nid = 0
        #print('nxyz_nodes=%s' % nxyz_nodes)
        mmax = amax(nodes, axis=0)
        mmin = amin(nodes, axis=0)
        dim_max = (mmax - mmin).max()
        self.gui.create_global_axes(dim_max)
        for i in range(nxyz_nodes):
            points.InsertPoint(nid, nodes[i, :])
            nid += 1
        log.info('nnodes=%s nquad_elements=%s nline_elements=%s' % (
            nxyz_nodes, nquad_elements, nline_elements))

        #elements -= 1
        for eid in range(nquad_elements):
            elem = vtkQuad()
            node_ids = elements[eid, :]
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, node_ids[0])
            point_ids.SetId(1, node_ids[1])
            point_ids.SetId(2, node_ids[2])
            point_ids.SetId(3, node_ids[3])
            grid.InsertNextCell(elem.GetCellType(), point_ids)

        for eid in range(nline_elements):
            elem = vtkLine()
            node_ids = line_elements[eid, :]
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, node_ids[0])
            point_ids.SetId(1, node_ids[1])
            grid.InsertNextCell(elem.GetCellType(), point_ids)

        grid.SetPoints(points)
        grid.Modified()

        # load results - regions/loads
        self.gui.scalar_bar_actor.VisibilityOn()
        self.gui.scalar_bar_actor.Modified()

        note = ''
        self.gui.isubcase_name_map = {1: ['AVL%s' % note, '']}
        cases = OrderedDict()
        ID = 1

        form, cases, node_ids, element_ids = self._fill_avl_case(
            cases, ID, nnodes, nelements, surfaces)
        self.gui.node_ids = node_ids
        self.gui.element_ids = element_ids
        self.gui._finish_results_io2(model_name, form, cases)

    #def clear_adb(self):
        #pass

    #def load_adb_results(self, cart3d_filename):
        #raise NotImplementedError()


    def _fill_avl_case(self, cases, ID, nnodes, nelements, surfaces):
        #results_form = []
        geometry_form = [
            ('NodeID', 0, []),
            ('ElementID', 1, []),
            ('SurfaceID', 2, []),
        ]

        nids = arange(1, nnodes + 1)
        eids = arange(1, nelements + 1)

        assert len(eids) == nelements, len(eids)

        nid_res = GuiResult(ID, 'NodeID', 'NodeID', 'node', nids)
        eid_res = GuiResult(ID, 'ElementID', 'ElementID', 'centroid', eids)
        surface_res = GuiResult(ID, 'SurfaceID', 'SurfaceID', 'centroid', surfaces)

        i = 0
        cases[i] = (nid_res, (0, 'NodeID'))
        cases[i + 1] = (eid_res, (0, 'ElementID'))
        cases[i + 2] = (surface_res, (0, 'SurfaceID'))

        form = [
            ('Geometry', None, geometry_form),
        ]
        return form, cases, nids, eids
