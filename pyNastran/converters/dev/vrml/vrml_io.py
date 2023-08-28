import numpy as np
import vtk

from pyNastran.gui.gui_objects.gui_result import GuiResult, NormalResult
from pyNastran.gui.utils.vtk.vtk_utils import (
    create_vtk_cells_of_constant_element_types, numpy_to_vtk_points)
from pyNastran.converters.dev.vrml.vrml import read_vrml


class Vrml_io:
    def __init__(self, gui):
        self.gui = gui

    def get_vrml_wildcard_geometry_results_functions(self):
        """gets the loader methods"""
        data = ('VRML',
                'VRML(*.wrl)', self.load_vrml_geometry,
                None, None)
        return data

    def load_vrml_geometry(self, vrml_filename, name='main', plot=True):
        """loads a VRML file into the GUI"""
        model_name = name
        #key = self.case_keys[self.icase]
        #case = self.result_cases[key]
        self.gui.eid_maps[name] = {}
        self.gui.nid_maps[name] = {}

        skip_reading = self.gui._remove_old_geometry(vrml_filename)
        if skip_reading:
            return

        nodes, quads, tris = read_vrml(vrml_filename)
        points = numpy_to_vtk_points(nodes)

        nnodes = nodes.shape[0]
        assert nnodes > 0
        ntris = len(tris)
        nquads = len(quads)
        nelements = ntris + nquads

        self.gui.nelements = nelements
        self.gui.nnodes = nnodes

        self.gui.log.info(f"nnodes={nnodes} nelements={nelements} (nquads={nquads} ntris={ntris})")
        assert nelements > 0, nelements

        grid = self.gui.grid
        grid.Allocate(self.gui.nelements, 1000)

        #self.gui.nid_map = {}

        etypes = []
        elements = []
        if ntris:
            elements.append(tris)
            etypes.append(5) # vtkTriangle().GetCellType()
        if nquads:
            elements.append(quads)
            etypes.append(9) # vtkQuad().GetCellType()

        assert len(etypes) > 0, etypes
        #self.gui.model.elements = elements
        #self.gui.model.etypes = etypes
        create_vtk_cells_of_constant_element_types(grid, elements, etypes)

        self.gui.nelements = nelements
        grid.SetPoints(points)
        grid.Modified()
        #------------------------------------------------
        self.gui.scalar_bar_actor.VisibilityOn()
        self.gui.scalar_bar_actor.Modified()

        self.gui.isubcase_name_map = {1: ['AFLR UGRID Surface', '']}
        cases = {}
        ID = 1

        node_ids = np.arange(1, nnodes + 1, dtype='int32')
        element_ids = np.arange(1, nelements + 1, dtype='int32')
        self.gui.node_ids = node_ids
        self.gui.element_ids = element_ids

        #eids = arange(1, len(elements) + 1, dtype='int32')
        #nids = arange(1, len(nodes) + 1, dtype='int32')
        #regions = np.array(regions, dtype='int32')

        icase = 0
        colormap = self.gui.settings.colormap
        eid_res = GuiResult(ID, header='ElementID', title='ElementID',
                            location='centroid', scalar=element_ids)
        nid_res = GuiResult(ID, header='NodeID', title='NodeID',
                            location='node', scalar=node_ids)
        nxyz_res = NormalResult(0, 'Normals', 'Normals',
                                nlabels=2, labelsize=5, ncolors=2,
                                colormap=colormap, data_format='%.1f',
                                uname='NormalResult')

        geometry_form = [
            ('ElementID', icase, []),
            ('NodeID', icase + 1, []),
            ('Normals', icase + 2, []),
        ]
        cases[icase] = (eid_res, (ID, 'ElementID'))
        cases[icase + 1] = (nid_res, (ID, 'NodeID'))
        cases[icase + 2] = (nxyz_res, (0, 'Normals'))

        form = geometry_form
        if plot:
            self.gui._finish_results_io2(model_name, form, cases)
