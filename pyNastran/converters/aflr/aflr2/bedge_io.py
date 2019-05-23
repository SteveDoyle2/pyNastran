"""
defines:
 - BEdge_IO

"""
import os
from collections import OrderedDict

import numpy as np
import vtk

from pyNastran.converters.aflr.aflr2.aflr2 import read_bedge
from pyNastran.gui.gui_objects.gui_result import GuiResult
from pyNastran.gui.utils.vtk.vtk_utils import (
    create_vtk_cells_of_constant_element_type, numpy_to_vtk_points)
from pyNastran.gui.qt_files.colors import BLACK_FLOAT

class BEdge_IO:
    """creates BEdge_IO"""
    def __init__(self, gui):
        self.gui = gui

    def get_bedge_wildcard_geometry_results_functions(self):
        data = (
            'AFLR2 BEdge',
            'AFLR2 BEdge (*.bedge)', self.load_bedge_geometry,
            None, None)
        return data

    def load_bedge_geometry(self, bedge_filename, name='main', plot=True):
        model_name = name
        #skip_reading = self.remove_old_openfoam_geometry(openfoam_filename)
        #if skip_reading:
        #    return

        self.gui.modelType = 'bedge'

        model = read_bedge(bedge_filename)
        self.gui.log.info('bedge_filename = %s' % bedge_filename)
        nnodes = model.nodes.shape[0]
        nbars = model.bars.shape[0]
        nelements = nbars

        nodes = model.nodes
        self.gui.nelements = nelements
        self.gui.nnodes = nnodes

        self.gui.log.debug("nNodes = %s" % self.gui.nnodes)
        self.gui.log.debug("nElements = %s" % self.gui.nelements)
        assert nelements > 0, nelements

        self.gui.create_alternate_vtk_grid(
            'nodes', color=BLACK_FLOAT, line_width=5, opacity=1., point_size=3,
            representation='point')
        alt_grid = self.gui.alt_grids['nodes']
        alt_grid.Allocate(nnodes, 1000)

        grid = self.gui.grid
        grid.Allocate(self.gui.nelements, 1000)

        points = numpy_to_vtk_points(nodes)

        mmax = np.amax(nodes, axis=0)
        mmin = np.amin(nodes, axis=0)
        unused_dim_max = (mmax - mmin).max()
        self.gui.log.info('max = %s' % mmax)
        self.gui.log.info('min = %s' % mmin)
        #print('dim_max =', dim_max)
        #self.update_axes_length(dim_max)

        etype = 1  # vtk.vtkVertex().GetCellType()
        #elements = np.arange(0, len(nodes), dtype='int32')
        #assert len(elements) == len(nodes)
        #create_vtk_cells_of_constant_element_type(alt_grid, elements, etype)
        for inode, unused_node in enumerate(nodes):
            elem = vtk.vtkVertex()
            elem.GetPointIds().SetId(0, inode)
            alt_grid.InsertNextCell(etype, elem.GetPointIds())

        bars = model.bars

        etype = 3  # vtkLine().GetCellType()
        create_vtk_cells_of_constant_element_type(grid, bars, etype)

        self.gui.nelements = nelements
        alt_grid.SetPoints(points)
        grid.SetPoints(points)
        grid.Modified()

        # loadBedgeResults - regions/loads
        #self.TurnTextOn()
        self.gui.scalar_bar_actor.VisibilityOn()
        self.gui.scalar_bar_actor.Modified()

        self.gui.isubcase_name_map = {1: ['AFLR BEDGE', '']}
        cases = OrderedDict()
        ID = 1

        self.gui._add_alt_actors(self.gui.alt_grids)
        form, cases, node_ids, element_ids = self._fill_bedge_case(
            bedge_filename, cases, ID, nnodes, nelements, model)
        self.gui.node_ids = node_ids
        self.gui.element_ids = element_ids
        if plot:
            self.gui._finish_results_io2(model_name, form, cases)

    def clear_bedge(self):
        pass

    def _fill_bedge_case(self, bedge_filename, cases, unused_ID, nnodes, nelements, model):
        """
        creates the data for the sidebar 'form' and the result 'cases'
        """
        base, ext = os.path.splitext(bedge_filename)
        assert ext == '.bedge', bedge_filename

        unused_tag_filename = base + '.tags'

        #cases_new = []
        has_tag_data = False
        results_form = []

        geometry_form = [
            ('ElementID', 0, []),
            ('NodeID', 1, []),
            #('Region', 1, []),
            ('CurveID', 2, []),
            ('SubcurveID', 3, []),
            ('GridBC', 4, []),

            #('TurnAngle', 5, [])
            #('normSpacing', 6, []),
            #('BL_thick', 7, []),
        ]

        form = [
            ('Geometry', None, geometry_form),
        ]

        eids = np.arange(1, nelements + 1)
        nids = np.arange(1, nnodes + 1)


        eid_res = GuiResult(0, header='ElementID', title='ElementID',
                            location='centroid', scalar=eids)
        nid_res = GuiResult(0, header='NodeID', title='NodeID',
                            location='node', scalar=nids)
        curve_res = GuiResult(0, header='CurveID', title='CurveID',
                              location='centroid', scalar=model.curves)
        subcurve_res = GuiResult(0, header='SubcurveID', title='SubcurveID',
                                 location='centroid', scalar=model.subcurves)
        gridbc_res = GuiResult(0, header='GridBC', title='GridBC',
                               location='centroid', scalar=model.grid_bcs)

        icase = 0
        cases[icase] = (eid_res, (0, 'ElementID'))
        cases[icase + 1] = (nid_res, (0, 'NodeID'))
        cases[icase + 2] = (curve_res, (0, 'CurveID'))
        cases[icase + 3] = (subcurve_res, (0, 'SubcurveID'))
        cases[icase + 4] = (gridbc_res, (0, 'GridBC'))
        icase += 5


        #if 0:
            #surf_ids = element_props[:, 0]
            #recon_flags = element_props[:, 1]
            #grid_bcs = element_props[:, 2]

        if hasattr(model, 'turn_angle'):
            gf = ('TurnAngle', 5, [])
            geometry_form.append(gf)
            turnangle_res = GuiResult(0, header='TurnAngle', title='TurnAngle',
                                      location='centroid',
                                      scalar=np.degrees(np.abs(model.turn_angle)))
            cases[icase] = (turnangle_res, (0, 'TurnAngle'))
            icase += 1

        #norm_spacing = model.node_props[:, 0]
        #bl_thickness = model.node_props[:, 1]
        #cases[(ID, 1, 'normSpacing', 1, 'node', '%.3e')] = norm_spacing
        #cases[(ID, 2, 'BL_thick',    1, 'node', '%.3e')] = bl_thickness

        form = [
            ('Geometry', None, geometry_form),
        ]
        if has_tag_data:
            tag_form = []
            form.append(('Tag Data', None, tag_form),)

        results_form = []
        if len(results_form):
            form.append(('Results', None, results_form))
        return form, cases, nids, eids
