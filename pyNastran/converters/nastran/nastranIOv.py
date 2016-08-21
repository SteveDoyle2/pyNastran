# pylint: disable=E1101
"""
Defines the GUI IO file for Nastran.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import os
import sys
from copy import deepcopy
from collections import defaultdict, OrderedDict
import traceback
from six import iteritems, itervalues, StringIO
from six.moves import zip, range


#VTK_TRIANGLE = 5
#VTK_QUADRATIC_TRIANGLE = 22

#VTK_QUAD = 9
#VTK_QUADRATIC_QUAD = 23

#VTK_TETRA = 10
#VTK_QUADRATIC_TETRA = 24

#VTK_WEDGE = 13
#VTK_QUADRATIC_WEDGE = 26

#VTK_HEXAHEDRON = 12
#VTK_QUADRATIC_HEXAHEDRON = 25

import numpy as np
from numpy.linalg import norm

#: makes vtk work on certain builds of vtk
#: we have to call this before vtk; you can't just try-except it
#: unused_import
import PyQt4

import vtk
from vtk import (vtkTriangle, vtkQuad, vtkTetra, vtkWedge, vtkHexahedron,
                 vtkQuadraticTriangle, vtkQuadraticQuad, vtkQuadraticTetra,
                 vtkQuadraticWedge, vtkQuadraticHexahedron,
                 vtkPyramid) #vtkQuadraticPyramid

#from pyNastran import is_release
from pyNastran.utils import integer_types
from pyNastran.bdf.bdf import (BDF,
                               CAERO1, CAERO2, CAERO3, CAERO4, CAERO5,
                               CQUAD4, CQUAD8, CQUADR, CSHEAR,
                               CTRIA3, CTRIA6, CTRIAR, CTRIAX6,
                               CONM2,
                               LOAD)

from pyNastran.bdf.cards.elements.shell import ShellElement
from pyNastran.bdf.cards.elements.bars import LineElement
from pyNastran.bdf.cards.elements.springs import SpringElement
from pyNastran.bdf.cards.elements.solid import (
    CTETRA4, CTETRA10, CPENTA6, CPENTA15,
    CHEXA8, CHEXA20, CIHEX1,
    CPYRAM5, CPYRAM13,
)

from pyNastran.converters.nastran.displacements import NastranDisplacementResults
from pyNastran.gui.gui_objects.gui_result import GuiResult

from pyNastran.op2.op2 import OP2
#from pyNastran.f06.f06_formatting import get_key0
try:
    from pyNastran.op2.op2_geom import OP2Geom
    is_geom = True
except ImportError:
    is_geom = False

piover2 = np.pi / 2.
piover3 = np.pi / 3.

class NastranIO(object):
    """
    Defines the GUI class for Nastran.
    """
    def __init__(self):
        # new options, no way to access them through the gui
        # they control results generation
        self.make_xyz = False
        self.make_offset_normals_dim = True
        self.make_nnodes_result = False  # make_offset_normals_dim must be True for this to work
        self.make_released_dofs1 = False
        self.make_released_dofs2 = False
        self.plot_applied_loads = True
        self.plot_pressures = True  # make_offset_normals_dim must be True for this to work

        # these are local variables
        self.model = None
        self.model_results = None
        self.bar_lines = None
        self.bar_eids = None
        self.xyz_cid0 = None
        self.normals = None
        self.save_data = True # was False
        #--------------------------------------

        #: flips the nastran CAERO subpaneling
        #:   False -> borders of CAEROs can be seen
        #:   True  -> individual subpanels can be seen
        self.show_caero_sub_panels = False

        #: coordinate systems can be messy, so this is the
        #: list of coords to show
        self.show_cids = []

        self.show_caero_actor = True  # show the caero mesh
        self.show_control_surfaces = True
        self.show_conm = True

        self.element_ids = None
        self.node_ids = None
        self.nid_map = None
        self.eid_map = None
        self.nNodes = None
        self.nElements = None
        self.model_type = None
        self.iSubcaseNameMap = None
        self.has_caero = False
        self.dependents_nodes = set([])
        self.i_transform = {}
        self.transforms = {}

    def get_nastran_wildcard_geometry_results_functions(self):
        """
        gets the Nastran wildcard loader used in the file load menu
        """
        if is_geom:
            geom_methods = 'Nastran BDF (*.bdf; *.dat; *.nas; *.op2; *.pch)'
        else:
            geom_methods = 'Nastran BDF (*.bdf; *.dat; *.nas; *.pch)'

        data = (
            'Nastran',
            geom_methods, self.load_nastran_geometry,
            'Nastran OP2 (*.op2)', self.load_nastran_results)
        return data

    def _cleanup_nastran_tools_and_menu_items(self):
        """
        hides the Nastran toolbar when loading another format
        """
        #self.menu_help.menuAction().setVisible(True)
        #self.menu_help2.menuAction().setVisible(False)
        self.nastran_toolbar.setVisible(False)
        self.actions['nastran'].setVisible(False)

    def _create_nastran_tools_and_menu_items(self):
        """
        creates the Nastran toolbar when loading a Nastran file
        """
        tools = [
            #('about_nastran', 'About Nastran GUI', 'tabout.png', 'CTRL+H', 'About Nastran GUI and help on shortcuts', self.about_dialog),
            #('about', 'About Orig GUI', 'tabout.png', 'CTRL+H', 'About Nastran GUI and help on shortcuts', self.about_dialog),
        ]
        #self.menu_help2 = self.menubar.addMenu('&HelpMenuNew')
        #self.menu_help.menuAction().setVisible(False)
        if hasattr(self, 'nastran_toolbar'):
            self.nastran_toolbar.setVisible(True)
            self.actions['nastran'].setVisible(True)
        else:
            self.nastran_toolbar = self.addToolBar('Nastran Toolbar')
            self.nastran_toolbar.setObjectName('nastran_toolbar')
            #self.nastran_toolbar.setStatusTip("Show/Hide nastran toolbar")
            self.actions['nastran'] = self.nastran_toolbar.toggleViewAction()
            self.actions['nastran'].setStatusTip("Show/Hide application toolbar")
        #self.file.menuAction().setVisible(False)
        #self.menu_help.

        #self.actions['about'].Disable()

        menu_items = [
            #(self.menu_help2, ('about_nastran',)),
            (self.nastran_toolbar, ('caero', 'caero_subpanels', 'conm2'))
            #(self.menu_window, tuple(menu_window)),
            #(self.menu_help, ('load_geometry', 'load_results', 'script', '', 'exit')),
            #(self.menu_help2, ('load_geometry', 'load_results', 'script', '', 'exit')),
        ]
        return tools, menu_items

    def toggle_caero_panels(self):
        """
        Toggle the visibility of the CAERO panels. The visibility of the sub panels
        or panels will be set according to the current show_caero_sub_panels state.
        """
        if not self.has_caero:
            return
        self.show_caero_actor = not self.show_caero_actor
        if self.show_caero_actor:
            if self.show_caero_sub_panels:
                self.geometry_actors['caero_subpanels'].VisibilityOn()
                self.geometry_properties['caero_subpanels'].is_visble = True
            else:
                self.geometry_actors['caero'].VisibilityOn()
                self.geometry_properties['caero'].is_visble = True
        else:
            self.geometry_actors['caero'].VisibilityOff()
            self.geometry_properties['caero'].is_visble = False
            self.geometry_actors['caero_subpanels'].VisibilityOff()
            self.geometry_properties['caero_subpanels'].is_visble = False
            self.on_update_geometry_properties_window(self.geometry_properties)
        self.vtk_interactor.Render()

    def on_update_geometry_properties_window(self, geometry_properties):
        if self._edit_geometry_properties_window_shown:
            self._edit_geometry_properties.on_update_geometry_properties_window(
                geometry_properties)

    def toggle_caero_sub_panels(self):
        """
        Toggle the visibility of the CAERO sub panels
        """
        if not self.has_caero:
            return
        self.show_caero_sub_panels = not self.show_caero_sub_panels
        if self.show_caero_actor:
            if self.show_caero_sub_panels:
                self.geometry_actors['caero'].VisibilityOff()
                self.geometry_properties['caero'].is_visble = False

                self.geometry_actors['caero_subpanels'].VisibilityOn()
                self.geometry_properties['caero_subpanels'].is_visble = True
            else:
                self.geometry_actors['caero'].VisibilityOn()
                self.geometry_properties['caero'].is_visble = True

                self.geometry_actors['caero_subpanels'].VisibilityOff()
                self.geometry_properties['caero_subpanels'].is_visble = False
        self.vtk_interactor.Render()

    def toggle_conms(self):
        """
        Toggle the visibility of the CONMS
        """
        self.show_conm = not self.show_conm
        if 'conm2' in self.geometry_actors:
            if self.show_conm:
                self.geometry_actors['conm2'].VisibilityOn()
                self.geometry_properties['conm2'].is_visble = True
            else:
                self.geometry_actors['conm2'].VisibilityOff()
                self.geometry_properties['conm2'].is_visble = False
        self.vtk_interactor.Render()

    def _create_coord(self, dim_max, cid, coord, cid_type):
        origin = coord.origin
        beta = coord.beta().T
        self.create_coordinate_system(dim_max, label='%s' % cid, origin=origin, matrix_3x3=beta,
                                      Type=cid_type)

    def _create_nastran_coords(self, model, dim_max):
        cid_types = {
            'R' : 'xyz',
            'C' : 'Rtz',
            'S' : 'Rtp',
        }
        self.create_global_axes(dim_max)
        self.show_cids = True
        for cid, coord in sorted(iteritems(model.coords)):
            if cid == 0:
                continue
            cid_type = cid_types[coord.Type]
            if self.show_cids is True:
                self._create_coord(dim_max, cid, coord, cid_type)
            elif isinstance(self.show_cids, integer_types):
                if cid == self.show_cids:
                    self._create_coord(dim_max, cid, coord, cid_type)
            elif isinstance(self.show_cids, (list, tuple, np.ndarray)):
                if cid in self.show_cids:
                    # .. todo:: has issues in VTK 6 I think due to lack of self.grid.Update()
                    self._create_coord(dim_max, cid, coord, cid_type)
            else:
                print('skipping cid=%s; use a script and set self.show_cids=[%s] to view' % (
                    cid, cid))

    def _remove_old_nastran_geometry(self, bdf_filename):
        #return self._remove_old_geometry(bdf_filename)

        # skip_reading = self.removeOldGeometry(bdf_filename)
        skip_reading = False
        if bdf_filename is None or bdf_filename is '':
            #self.grid = vtk.vtkUnstructuredGrid()
            #self.gridResult = vtk.vtkFloatArray()
            #self.emptyResult = vtk.vtkFloatArray()
            #self.vectorResult = vtk.vtkFloatArray()
            #self.scalarBar.VisibilityOff()
            skip_reading = True
            return skip_reading
        else:
            self.turn_text_off()
            self.grid.Reset()

            #self.gridResult = vtk.vtkFloatArray()
            #self.gridResult.Reset()
            #self.gridResult.Modified()
            #self.eid_map = {}
            #self.nid_map = {}

            self.result_cases = {}
            self.ncases = 0
        for i in ('case_keys', 'icase', 'iSubcaseNameMap'):
            if hasattr(self, i):  # TODO: is this correct???
                del i
        return skip_reading

    def get_xyz_in_coord(self, model, points, cid=0, dtype='float32'):
        nid_map = self.nid_map
        assert cid == 0, cid
        nnodes = len(model.nodes)
        nspoints = 0
        spoints = None
        if model.spoints:
            spoints = model.spoints.points
            nspoints = len(spoints)

        xyz_cid0 = np.zeros((nnodes + nspoints, 3), dtype=dtype)
        if nspoints:
            nids = model.nodes.keys()
            newpoints = nids + list(spoints)
            newpoints.sort()
            for i, nid in enumerate(newpoints):
                if nid in spoints:
                    nid_map[nid] = i
                else:
                    node = model.nodes[nid]
                    xyz_cid0[i, :] = node.get_position()
                    nid_map[nid] = i
                points.InsertPoint(i, *xyz_cid0[i, :])
        else:
            for i, (nid, node) in enumerate(sorted(iteritems(model.nodes))):
                xyz = node.get_position()
                xyz_cid0[i, :] = xyz
                points.InsertPoint(i, *xyz)
                nid_map[nid] = i
        return xyz_cid0

    def load_nastran_geometry(self, bdf_filename, dirname, name='main', plot=True):
        self.eid_maps[name] = {}
        self.nid_maps[name] = {}
        self.i_transform = {}
        self.transforms = {}
        #print('bdf_filename=%r' % bdf_filename)
        #key = self.case_keys[self.icase]
        #case = self.result_cases[key]

        skip_reading = self._remove_old_nastran_geometry(bdf_filename)
        pink = (0.98, 0.4, 0.93)
        orange = (219/255., 168/255., 13/255.)
        #blue = (0., 0., 1.)
        red = (1., 0., 0.)
        # if 0:
            # yellow = (1., 1., 0.)
            # line_width = 3
            # opacity = 1
            # alt_grids = [
                # ['caero', yellow, line_width, opacity],
                # ['caero_subpanels', yellow, line_width, opacity],
            # ]
            # skip_reading = self._remove_old_geometry2(bdf_filename, alt_grids=alt_grids)
        if skip_reading:
            return

        if plot:
            self.scalarBar.VisibilityOff()
            self.scalarBar.Modified()

        ext = os.path.splitext(bdf_filename)[1].lower()
        punch = False
        if ext == '.pch':
            punch = True

        xref_loads = True
        #print('ext=%r is_geom=%s' % (ext, is_geom))
        if ext == '.op2' and is_geom:
            model = OP2Geom(make_geom=True, debug=False, log=self.log,
                            debug_file=None)
            model.clear_results()
            model.read_op2(op2_filename=bdf_filename)
            model.cross_reference(xref=True, xref_loads=xref_loads,
                                  xref_constraints=False,
                                  xref_nodes_with_elements=False)
            #model.safe_cross_reference(xref=True, xref_loads=xref_loads,
                                       #xref_constraints=False)
        else:  # read the bdf/punch
            model = BDF(log=self.log, debug=True)
            self.model_type = 'nastran'
            model.read_bdf(bdf_filename,
                           punch=punch, xref=False)
            # model.cross_reference(xref=True, xref_loads=xref_loads,
                                  # xref_constraints=False)
            model.safe_cross_reference(xref=True, xref_loads=xref_loads,
                                       xref_constraints=False,
                                       xref_nodes_with_elements=False)

        # get indicies and transformations for displacements
        self.i_transform, self.transforms = model.get_displacement_index_transforms()

        nnodes = len(model.nodes)
        nspoints = 0
        spoints = None
        if model.spoints:
            spoints = model.spoints.points
            nspoints = len(spoints)

        assert nnodes + nspoints > 0, model.card_count
        nelements = model.nelements
        assert nelements > 0

        self.nNodes = nnodes + nspoints
        self.nElements = nelements  # approximate...


        # count caeros
        ncaeros_sub = 0
        ncaero_sub_points = 0
        for caero in itervalues(model.caeros):
            if hasattr(caero, 'panel_points_elements'):
                npoints, ncelements = caero.get_npanel_points_elements()
                ncaeros_sub += npoints
                ncaero_sub_points += ncelements
            elif isinstance(caero, CAERO2):
                pass
            else:
                print('%r doesnt support panel_points_elements' % caero.type)
        #ncaeros = model.ncaeros

        ncaeros_points = 0
        ncaeros = 0
        for caero in model.caeros:
            if isinstance(caero, (CAERO1, CAERO3, CAERO4, CAERO5)):
                ncaeros_points += 4
                ncaeros += 1
            elif isinstance(caero, CAERO2):
                points, elems = caero.get_points_elements_3d()
                ncaeros_points += points.shape[0]
                ncaeros += elems.shape[0]

        box_id_to_caero_element_map = {}
        num_prev = 0
        ncaeros_sub = 0
        if model.caeros:
            caero_points = []
            for eid, caero in sorted(iteritems(model.caeros)):
                if caero.type == 'CAERO1':
                    ncaeros_sub += 1
                    pointsi, elementsi = caero.panel_points_elements()
                    caero_points.append(pointsi)
                    for i, box_id in enumerate(caero.box_ids.flat):
                        box_id_to_caero_element_map[box_id] = elementsi[i, :] + num_prev
                    num_prev += pointsi.shape[0]
                else:
                    print('caero\n%s' % caero)
            if ncaeros_sub:
                caero_points = np.vstack(caero_points)
            self.has_caero = True
        if ncaeros_sub == 0:
            caero_points = np.empty((0, 3))

        # check for any control surfcaes
        has_control_surface = False
        if model.aesurfs:
            cs_box_ids = []
            has_control_surface = True
            ncaeros_cs = 0
            #ncaero_cs_points = 0
            if 'caero_control_surfaces' not in self.alt_grids:
                self.create_alternate_vtk_grid(
                    'caero_control_surfaces', color=pink, line_width=5, opacity=1.0,
                    representation='surface')
            for aid, aesurf in iteritems(model.aesurfs):
                aelist = aesurf.alid1
                ncaeros_cs += len(aelist.elements)
                cs_box_ids.extend(aelist.elements)

                if aesurf.alid2 is not None:
                    aelist = aesurf.alid2
                    ncaeros_cs += len(aelist.elements)
                    cs_box_ids.extend(aelist.elements)

        self.log_info("nNodes=%i nElements=%i" % (self.nNodes, self.nElements))
        msg = model.get_bdf_stats(return_type='list')
        #self.log_info(msg)
        for msgi in msg:
            model.log.debug(msgi)

        if 'CONM2' in model.card_count:
            nconm2 = model.card_count['CONM2']
        else:
            nconm2 = 0
        #self.gridResult.SetNumberOfComponents(self.nElements)
        if nconm2 > 0:
            self.create_alternate_vtk_grid(
                'conm2', color=orange, line_width=5, opacity=1., point_size=4,
                representation='point')

        # Allocate grids
        self.grid.Allocate(self.nElements, 1000)
        if self.has_caero:
            yellow = (1., 1., 0.)
            if 'caero' not in self.alt_grids:
                self.create_alternate_vtk_grid(
                    'caero', color=yellow, line_width=3, opacity=1.0,
                    representation='toggle', is_visible=True)
            if 'caero_subpanels' not in self.alt_grids:
                self.create_alternate_vtk_grid(
                    'caero_subpanels', color=yellow, line_width=3, opacity=1.0,
                    representation='toggle', is_visible=False)

            self.alt_grids['caero'].Allocate(ncaeros, 1000)
            self.alt_grids['caero_subpanels'].Allocate(ncaeros_sub, 1000)
            if has_control_surface:
                self.alt_grids['caero_control_surfaces'].Allocate(ncaeros_cs, 1000)

        if nconm2 > 0:
            self.alt_grids['conm2'].Allocate(nconm2, 1000)

        points = vtk.vtkPoints()
        points.SetNumberOfPoints(self.nNodes)
        #self.gridResult.Allocate(self.nNodes, 1000)
        #vectorReselt.SetNumberOfComponents(3)
        #elem.SetNumberOfPoints(nNodes)

        if self.save_data:
            self.model = model

        xyz_cid0 = self.get_xyz_in_coord(model, points, cid=0, dtype='float32')
        self.xyz_cid0 = xyz_cid0

        maxi = xyz_cid0.max(axis=0)
        mini = xyz_cid0.min(axis=0)
        assert len(maxi) == 3, len(maxi)
        xmax, ymax, zmax = maxi
        xmin, ymin, zmin = mini
        dim_max = max(xmax-xmin, ymax-ymin, zmax-zmin)

        self._create_nastran_coords(model, dim_max)

        self.log_info("xmin=%s xmax=%s dx=%s" % (xmin, xmax, xmax-xmin))
        self.log_info("ymin=%s ymax=%s dy=%s" % (ymin, ymax, ymax-ymin))
        self.log_info("zmin=%s zmax=%s dz=%s" % (zmin, zmax, zmax-zmin))
        self.create_splines(model, box_id_to_caero_element_map, caero_points)

        if model.suport:
            ids = []
            for suport in model.suport:
                #print(suport)
                idsi = suport.node_ids
                #print('idsi =', idsi)
                ids += idsi
            grid_name = 'SUPORT'
            self.create_alternate_vtk_grid(
                grid_name, color=red, opacity=1.0, point_size=42,
                representation='point', is_visible=True)
        j = 0
        nid_to_pid_map, icase, cases, form = self.map_elements(
            points, self.nid_map, model, j, dim_max, plot=plot, xref_loads=xref_loads)

        #if 0:
            #nsprings = 0
            #if 0:
                #for eid, element in sorted(iteritems(model.elements)):
                    #if(isinstance(element, LineElement) or
                       #isinstance(element, SpringElement) or
                       #element.type in ['CBUSH', 'CBUSH1D', 'CFAST', 'CROD', 'CONROD',
                                        #'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
                                        #'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4',
                                        #'CDAMP5', 'CVISC', ]):
                        #node_ids = element.node_ids
                        #if None in node_ids:
                            #nsprings += 1

        # fill grids
        if 'caero' in self.alt_grids:
            self.set_caero_grid(ncaeros_points, model)
            self.set_caero_subpanel_grid(ncaero_sub_points, model)
            if has_control_surface:
                self.set_caero_control_surface_grid(
                    'caero_control_surfaces', cs_box_ids,
                    box_id_to_caero_element_map, caero_points)
        if nconm2 > 0:
            self.set_conm_grid(nconm2, dim_max, model)

        self.set_spc_grid(dim_max, model, nid_to_pid_map)
        icase = self._fill_bar_yz(dim_max, model, icase, cases, form)
        assert icase is not None

        #------------------------------------------------------------
        # loads
        try:
            subcase_ids = model.case_control_deck.get_subcase_list()
        except AttributeError:
            return nid_to_pid_map, icase, cases, form

        #print('dependent_nodes =', self.dependents_nodes)
        form0 = form[2]
        assert icase is not None
        for subcase_id in subcase_ids:
            if subcase_id == 0:
                continue
            self.log_debug('NastranIOv subcase_id = %s' % subcase_id)
            subcase = model.case_control_deck.subcases[subcase_id]

            subtitle = ''
            if 'SUBTITLE' in subcase:
                subtitle, options = subcase.get_parameter('SUBTITLE')
                del options

            load_str = 'Load Case=%i' % subcase_id if subtitle == '' else 'Load Case=%i; %s' % (
                subcase_id, subtitle)
            formi = (load_str, None, [])
            formii = formi[2]
            if self.plot_pressures:
                try:
                    icase = self._plot_pressures(
                        model, cases, formii, icase, subcase_id, subcase)
                except:
                    s = StringIO()
                    traceback.print_exc(file=s)
                    sout = s.getvalue()
                    self.log_error(sout)
                    print(sout)

            assert icase is not None
            if self.plot_applied_loads:
                try:
                    icase = self._plot_applied_loads(
                        model, cases, formii, icase, subcase_id, subcase)
                except:
                    s = StringIO()
                    traceback.print_exc(file=s)
                    sout = s.getvalue()
                    self.log_error(sout)
                    print(sout)

            if len(formii):
                form0.append(formi)

        #------------------------------------------------------------
        # add alternate actors
        self._add_alt_actors(self.alt_grids)

        # set default representation
        self._set_caero_representation(has_control_surface)

        for grid_name in ['suport', 'spc', 'mpc', 'mpc_dependent', 'mpc_independent']:
            if grid_name in self.geometry_actors:
                self.geometry_actors[grid_name].Modified()

        if plot:
            #self.log.info(cases.keys())
            self._finish_results_io2([form], cases)
        else:
            self._set_results([form], cases)

    def _set_caero_representation(self, has_control_surface):
        if 'caero_control_surfaces' in self.geometry_actors:
            self.geometry_properties['caero_control_surfaces'].opacity = 0.5

            if 'caero' not in self.geometry_actors:
                return
            self.geometry_actors['caero'].Modified()
            self.geometry_actors['caero_subpanels'].Modified()
            if has_control_surface:
                self.geometry_actors['caero_control_surfaces'].Modified()
            if hasattr(self.geometry_actors['caero'], 'Update'):
                self.geometry_actors['caero'].Update()
            if hasattr(self.geometry_actors['caero_subpanels'], 'Update'):
                self.geometry_actors['caero_subpanels'].Update()
            if has_control_surface and hasattr(self.geometry_actors['caero_subpanels'], 'Update'):
                self.geometry_actors['caero_control_surfaces'].Update()

    def create_splines(self, model, box_id_to_caero_element_map, caero_points):
        if model.splines:
            blue = (0., 0., 1.)
            # 0 - caero / caero_subpanel
            # 1 - control surface
            iaero = 2
            for spline_id, spline in sorted(model.splines.items()):
                # the control surfaces all lie perfectly on top of each other
                # such that we have z fighting, so based on the aero index,
                # we calculate a z offset.
                setg = spline.setg_ref
                structure_points = setg.get_ids()

                try:
                    aero_box_ids = spline.aero_element_ids
                except:
                    print(spline.object_attributes())
                    print(spline.object_methods())
                    raise

                zfighting_offset = 0.0001 * iaero
                grid_name = 'spline_%s_structure_points' % spline_id
                self.create_alternate_vtk_grid(
                    grid_name, color=blue, opacity=1.0, point_size=5,
                    representation='point', is_visible=False)
                self._add_nastran_nodes_to_grid(grid_name, structure_points, model)


                zfighting_offset = 0.0001 * (iaero + 1)
                grid_name = 'spline_%s_boxes' % spline_id
                self.create_alternate_vtk_grid(
                    grid_name, color=blue, opacity=0.3,
                    line_width=4,
                    representation='toggle', is_visible=False)
                self.set_caero_control_surface_grid(
                    grid_name, aero_box_ids, box_id_to_caero_element_map, caero_points,
                    zfighting_offset=zfighting_offset)
                iaero += 2

    def set_caero_grid(self, ncaeros_points, model, j=0):
        """
        Sets the CAERO panel geometry.

        Returns the current id counter.
        """
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(ncaeros_points)

        max_cpoints = []
        min_cpoints = []
        caero_grid = self.alt_grids['caero']
        for eid, element in sorted(iteritems(model.caeros)):
            if isinstance(element, (CAERO1, CAERO3, CAERO4, CAERO5)):
                # wing panel
                cpoints = element.get_points()
                max_cpoints.append(np.array(cpoints).max(axis=0))
                min_cpoints.append(np.array(cpoints).min(axis=0))

                elem = vtkQuad()
                elem.GetPointIds().SetId(0, j)
                elem.GetPointIds().SetId(1, j + 1)
                elem.GetPointIds().SetId(2, j + 2)
                elem.GetPointIds().SetId(3, j + 3)
                points.InsertPoint(j, *cpoints[0])
                points.InsertPoint(j + 1, *cpoints[1])
                points.InsertPoint(j + 2, *cpoints[2])
                points.InsertPoint(j + 3, *cpoints[3])
                caero_grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                j += 4
            elif isinstance(element, CAERO2):
                # slender body
                if 0:
                    # 1D version
                    cpoints = element.get_points()
                    max_cpoints.append(np.array(cpoints).max(axis=0))
                    min_cpoints.append(np.array(cpoints).min(axis=0))

                    elem = vtk.vtkLine()
                    elem.GetPointIds().SetId(0, j)
                    elem.GetPointIds().SetId(1, j + 1)

                    print(', '.join(dir(elem)))
                    #prop = elem.GetProperty()
                    #print(dir(prop))

                    points.InsertPoint(j, *cpoints[0])
                    points.InsertPoint(j + 1, *cpoints[1])
                    j += 2
                    caero_grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                else:
                    # 3D version
                    xyz, elems = element.get_points_elements_3d()
                    for elemi in elems:
                        elem = vtkQuad()
                        elem.GetPointIds().SetId(0, j)
                        elem.GetPointIds().SetId(1, j + 1)
                        elem.GetPointIds().SetId(2, j + 2)
                        elem.GetPointIds().SetId(3, j + 3)
                        n1, n2, n3, n4 = elemi

                        points.InsertPoint(j, *xyz[n1])
                        points.InsertPoint(j + 1, *xyz[n2])
                        points.InsertPoint(j + 2, *xyz[n3])
                        points.InsertPoint(j + 3, *xyz[n4])
                        caero_grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                        j += 4
                        #print('caero:\n', element)
                        #print('nids:\n', n1, n2, n3, n4)
                        #print('xyz:\n', xyz[n1], xyz[n2], xyz[n3], xyz[n4])
            else:
                self.log_info("skipping %s" % element.type)

        if ncaeros_points:
            self.log_info('CAERO.max = %s' % np.vstack(max_cpoints).max(axis=0))
            self.log_info('CAERO.min = %s' % np.vstack(min_cpoints).min(axis=0))
        self.alt_grids['caero'].SetPoints(points)
        return j

    def set_caero_subpanel_grid(self, ncaero_sub_points, model, j=0):
        """
        Sets the CAERO panel geometry.

        Returns the current id counter.
        """
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(ncaero_sub_points)

        vtk_type = vtkQuad().GetCellType()
        for eid, element in sorted(iteritems(model.caeros)):
            if isinstance(element, (CAERO1, CAERO3, CAERO4, CAERO5)):
                pointsi, elementsi = element.panel_points_elements()
                for ipoint, pointii in enumerate(pointsi):
                    points.InsertPoint(j + ipoint, *pointii)
                if eid == 10100:
                    continue
                elem = vtkQuad()
                for elementi in elementsi:
                    elem = vtkQuad()
                    elem.GetPointIds().SetId(0, j + elementi[0])
                    elem.GetPointIds().SetId(1, j + elementi[1])
                    elem.GetPointIds().SetId(2, j + elementi[2])
                    elem.GetPointIds().SetId(3, j + elementi[3])
                    self.alt_grids['caero_subpanels'].InsertNextCell(vtk_type, elem.GetPointIds())
                j += ipoint + 1
            else:
                self.log_info("skipping %s" % element.type)
        self.alt_grids['caero_subpanels'].SetPoints(points)
        return j

    def set_caero_control_surface_grid(self, name, cs_box_ids,
                                       box_id_to_caero_element_map,
                                       caero_points,
                                       zfighting_offset=0.001, j=0):
        points_list = []
        missing_boxes = []
        for box_id in cs_box_ids:
            try:
                ipoints = box_id_to_caero_element_map[box_id]
            except KeyError:
                missing_boxes.append(box_id)
                continue
            points_list.append(caero_points[ipoints, :])
        if missing_boxes:
            msg = 'Missing CAERO AELIST boxes: ' + str(missing_boxes)
            self.log_error(msg)

        points_list = np.array(points_list)
        ncaero_sub_points = len(np.unique(points_list.ravel()))

        points = vtk.vtkPoints()
        points.SetNumberOfPoints(ncaero_sub_points)

        vtk_type = vtkQuad().GetCellType()
        #print('cs_box_ids =', cs_box_ids)
        for ibox, box_id in enumerate(cs_box_ids):
            try:
                elementi = box_id_to_caero_element_map[box_id]
            except KeyError:
                print('cant find box_id=%i' % box_id)
                continue
            pointsi = caero_points[elementi]
            area = np.linalg.norm(np.cross(pointsi[2] - pointsi[0], pointsi[3] - pointsi[1])) / 2.
            if area == 0.0:
                print('box_id=%i has 0 area' % box_id)
                continue
            for ipoint, point in enumerate(pointsi):
                #print('point[%i, %i] = %s; A=%s' % (ibox, ipoint, point, area))
                # shift z to remove z-fighting with caero in surface representation
                point[1] += zfighting_offset
                point[2] += zfighting_offset
                points.InsertPoint(j + ipoint, *point)
            #print('')
            elem = vtkQuad()
            elem.GetPointIds().SetId(0, j)
            elem.GetPointIds().SetId(1, j + 1)
            elem.GetPointIds().SetId(2, j + 2)
            elem.GetPointIds().SetId(3, j + 3)
            self.alt_grids[name].InsertNextCell(vtk_type, elem.GetPointIds())
            j += ipoint + 1

        #if missing_boxes:
            #msg = 'Missing CAERO AELIST boxes: ' + str(missing_boxes)
            #self.log_error(msg)
        self.alt_grids[name].SetPoints(points)
        return j

    def set_caero_wireframe_points(self, name, aero_box_ids,
                                   box_id_to_caero_element_map,
                                   caero_points,
                                   structure_points, xyz_cid0,
                                   zfighting_offset=0.0, j=0):
        points_list = []
        missing_boxes = []
        for box_id in aero_box_ids:
            try:
                ipoints = box_id_to_caero_element_map[box_id]
            except KeyError:
                missing_boxes.append(box_id)
                continue
            points_list.append(caero_points[ipoints, :])
        if missing_boxes:
            msg = 'Missing CAERO AELIST boxes: ' + str(missing_boxes)
            self.log_error(msg)

        nid_map = self.nid_map
        structure_points2 = []
        for nid in structure_points:
            if nid in nid_map:
                structure_points2.append(nid)

        points_list = np.array(points_list)
        ncaero_sub_points = len(np.unique(points_list.ravel()))

        points = vtk.vtkPoints()
        points.SetNumberOfPoints(ncaero_sub_points)

        vtk_type = vtkQuad().GetCellType()
        for box_id in aero_box_ids:
            try:
                elementi = box_id_to_caero_element_map[box_id]
            except KeyError:
                continue
            pointsi = caero_points[elementi]
            for ipoint, point in enumerate(pointsi):
                # shift z to remove z-fighting with caero in surface representation
                point[1] += zfighting_offset
                point[2] += zfighting_offset
                points.InsertPoint(j + ipoint, *point)
            elem = vtkQuad()
            elem.GetPointIds().SetId(0, j)
            elem.GetPointIds().SetId(1, j + 1)
            elem.GetPointIds().SetId(2, j + 2)
            elem.GetPointIds().SetId(3, j + 3)
            self.alt_grids[name].InsertNextCell(vtk_type, elem.GetPointIds())
            assert ipoint == 3, ipoint
            j += ipoint + 1

        #xyz_cid0[i, :] = node.get_position()
        #nid_map[nid] = i

        vtk_type = vtk.vtkVertex().GetCellType()
        for i, nid in enumerate(structure_points2):
            ipoint = nid_map[nid]
            point = xyz_cid0[i, :]
            points.InsertPoint(j, *point)

            j += 1
        self.alt_grids[name].SetPoints(points)
        return j

    def set_conm_grid(self, nconm2, dim_max, model, j=0):
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(nconm2)

        sphere_size = self._get_sphere_size(dim_max)
        for eid, element in sorted(iteritems(model.masses)):
            if isinstance(element, CONM2):
                #xyz = element.nid.get_position()
                centroid = element.Centroid()
                #d = norm(xyz - c)
                points.InsertPoint(j, *centroid)

                #if 1:
                elem = vtk.vtkVertex()
                elem.GetPointIds().SetId(0, j)
                #else:
                    #elem = vtk.vtkSphere()
                    #elem.SetRadius(sphere_size)
                    #elem.SetCenter(points.GetPoint(j))

                self.alt_grids['conm2'].InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                j += 1
            else:
                self.log_info("skipping %s" % element.type)
        self.alt_grids['conm2'].SetPoints(points)
        #self.alt_grids['conm2'].Set

    def set_spc_grid(self, dim_max, model, nid_to_pid_map):
        #case_control = model.case_control_deck
        for subcase_id, subcase in sorted(iteritems(model.subcases)):
            #print('subcase_id=%s' % subcase_id)
            #print(subcase.params.keys())

            if 'SPC' in subcase:
                spc_id, options = subcase.get_parameter('SPC')
                if spc_id is not None:
                    nspcs = model.card_count['SPC'] if 'SPC' in model.card_count else 0
                    nspc1s = model.card_count['SPC1'] if 'SPC1' in model.card_count else 0
                    nspcds = model.card_count['SPCD'] if 'SPCD' in model.card_count else 0
                    if nspcs + nspc1s:
                        self._fill_spc(spc_id, nspcs, nspc1s, nspcds, dim_max,
                                       model, nid_to_pid_map)

            if 'MPC' in subcase:
                mpc_id, options = subcase.get_parameter('MPC')
                if mpc_id is not None:
                    nmpcs = model.card_count['MPC'] if 'MPC' in model.card_count else 0
                    if nmpcs:
                        self._fill_mpc(mpc_id, dim_max, model, nid_to_pid_map)
            else:
                self._fill_rigid(dim_max, model, nid_to_pid_map)

            if 'SUPORT1' in subcase.params:  ## TODO: should this be SUPORT?
                suport_id, options = subcase.get_parameter('SUPORT1')
                if 'SUPORT' in model.card_count or 'SUPORT1' in model.card_count:
                    if suport_id:
                        self._fill_suport(suport_id, dim_max, model)

    def get_SPCx_node_ids(self, model, spc_id, exclude_spcadd=False):
        """
        Get the SPC/SPCADD/SPC1/SPCAX IDs.

        Parameters
        -----------
        exclude_spcadd : bool
            you can exclude SPCADD if you just want a list of all the
            SPCs in the model.  For example, apply all the SPCs when
            there is no SPC=N in the case control deck, but you don't
            need to apply SPCADD=N twice.
        """
        try:
            spcs = model.spcs[spc_id]
        except KeyError:
            model.log.warning('spc_id=%s not found' % spc_id)
            return []

        node_ids = []
        for card in sorted(spcs):
            if card.type == 'SPC':
                nids = card.node_ids
            elif card.type == 'SPC1':
                nids = card.node_ids
            elif card.type == 'SPCADD':
                nids = []
                for new_spc_id in card.sets:
                    nidsi = self.get_SPCx_node_ids(model, new_spc_id, exclude_spcadd=False)
                    nids += nidsi
            else:
                self.log.warning('get_SPCx_node_ids doesnt supprt %r' % card.type)
                continue
            node_ids += nids
        return node_ids

    def get_SPCx_node_ids_c1(self, model, spc_id, exclude_spcadd=False):
        """
        Get the SPC/SPCADD/SPC1/SPCAX IDs.

        Parameters
        -----------
        exclude_spcadd : bool
            you can exclude SPCADD if you just want a list of all the
            SPCs in the model.  For example, apply all the SPCs when
            there is no SPC=N in the case control deck, but you don't
            need to apply SPCADD=N twice.
        """
        try:
            spcs = model.spcs[spc_id]
        except KeyError:
            model.log.warning('spc_id=%s not found' % spc_id)
            return {}

        node_ids_c1 = defaultdict(str)
        #print('spcs = ', spcs)
        for card in spcs:  # used to be sorted(spcs)
            if card.type == 'SPC':
                for nid, c1 in zip(card.gids, card.constraints):
                    assert nid is not None, card.gids
                    node_ids_c1[nid] += c1
            elif card.type == 'SPC1':
                nids = card.node_ids
                c1 = card.constraints
                for nid in nids:
                    node_ids_c1[nid] += c1
            elif card.type == 'SPCADD':
                nids = []
                for new_spc_id in card.sets:
                    nids_c1i = self.get_SPCx_node_ids_c1(model, new_spc_id, exclude_spcadd=False)
                    for nid, c1 in iteritems(nids_c1i):
                        node_ids_c1[nid] += c1
            else:
                self.log.warning('get_SPCx_node_ids_c1 doesnt supprt %r' % card.type)
                continue
        return node_ids_c1

    def _fill_spc(self, spc_id, nspcs, nspc1s, nspcds, dim_max, model, nid_to_pid_map):
        purple = (1., 0., 1.)
        self.create_alternate_vtk_grid('spc', color=purple, line_width=5, opacity=1.,
                                       point_size=5, representation='point', is_visible=False)

        # node_ids = self.get_SPCx_node_ids(model, spc_id, exclude_spcadd=False)
        node_ids_c1 = self.get_SPCx_node_ids_c1(model, spc_id, exclude_spcadd=False)

        node_ids = []
        for nid, c1 in iteritems(node_ids_c1):
            if nid_to_pid_map is not None:
                plot_node = False
                pids = nid_to_pid_map[nid]
                for pid in pids:
                    if pid == 0:
                        continue
                    if pid is None:
                        print(pid)
                        continue
                    prop = model.properties[pid]
                    if prop.type not in ['PSOLID', 'PLSOLID']:
                        plot_node = True
                if not plot_node:
                    # don't include 456 constraints if they're ONLY on solid elemetns
                    # if we had any bar/plate/etc. elements that use this node, we'll plot the node
                    if not('1' in c1 or '2' in c1 or '3' in c1):
                        continue
            node_ids.append(nid)

        node_ids = np.unique(node_ids)
        self._add_nastran_nodes_to_grid('spc', node_ids, model, nid_to_pid_map)

    def get_MPCx_node_ids_c1(self, model, mpc_id, exclude_mpcadd=False):
        r"""
        Get the MPC/MPCADD IDs.

        Parameters
        -----------
        exclude_spcadd : bool
            you can exclude MPCADD if you just want a list of all the
            MPCs in the model.  For example, apply all the MPCs when
            there is no MPC=N in the case control deck, but you don't
            need to apply MPCADD=N twice.

        I      I
          \   /
        I---D---I
        """
        lines = []
        try:
            mpcs = model.mpcs[mpc_id]
        except:
            model.log.warning('mpc_id=%s not found' % mpc_id)
            return []

        # dependent, independent
        for card in sorted(mpcs):
            if card.type == 'MPC':
                nids = card.node_ids
                nid0 = nids[0]
                #constraint0 = card.constraints[0]
                #enforced0 = card.enforced[0]
                #card.constraints[1:]
                for nid, enforced in zip(nids[1:], card.enforced[1:]):
                    if enforced != 0.0:
                        lines.append([nid0, nid])
            elif card.type == 'MPCADD':
                nids = []
                for new_mpc_id in card.sets:
                    linesi = self.get_MPCx_node_ids_c1(model, new_mpc_id, exclude_mpcadd=False)
                    lines += linesi
            else:
                self.log.warning('get_MPCx_node_ids_c1 doesnt supprt %r' % card.type)
                continue
        return lines

    def _fill_bar_yz(self, dim_max, model, icase, cases, form, debug=False):
        """
        plots the y, z vectors for CBAR & CBEAM elements
        """
        green = (0., 1., 0.)
        blue = (0., 0., 1.)
        card_types = ['CBAR', 'CBEAM']
        out = model.get_card_ids_by_card_types(card_types=card_types)
        bar_beam_eids = out['CBAR'] + out['CBEAM']
        self.bar_eids = {}
        self.bar_lines = {}
        if len(bar_beam_eids) == 0:
            return icase

        scale = 0.15
        lines_bar_y = []
        lines_bar_z = []

        bar_types = {
            # PBAR
            'bar' : [],

            # PBEAML/PBARL
            "ROD": [],
            "TUBE": [],
            "TUBE2" : [],
            "I": [],
            "CHAN": [],
            "T": [],
            "BOX": [],
            "BAR": [],
            "CROSS": [],
            "H": [],
            "T1": [],
            "I1": [],
            "CHAN1": [],
            "Z": [],
            "CHAN2": [],
            "T2": [],
            "BOX1": [],
            "HEXA": [],
            "HAT": [],
            "HAT1": [],
            "DBOX": [],  # was 12

            # PBEAM
            'beam' : [],

            # PBEAML specfic
            "L" : [],
        }  # for GROUP="MSCBML0"
        allowed_types = [
            'BAR', 'BOX', 'BOX1', 'CHAN', 'CHAN1', 'CHAN2', 'CROSS', 'DBOX',
            'H', 'HAT', 'HAT1', 'HEXA', 'I', 'I1', 'L', 'ROD',
            'T', 'T1', 'T2', 'TUBE', 'TUBE2', 'Z', 'bar', 'beam',
        ]

        # bar_types['bar'] = [ [...], [...], [...] ]
        #bar_types = defaultdict(lambda : defaultdict(list))

        found_bar_types = set([])
        #neids = len(self.element_ids)
        for bar_type, data in iteritems(bar_types):
            eids = []
            lines_bar_y = []
            lines_bar_z = []
            bar_types[bar_type] = (eids, lines_bar_y, lines_bar_z)

        no_axial = np.zeros(self.element_ids.shape, dtype='int32')
        no_torsion = np.zeros(self.element_ids.shape, dtype='int32')

        if self.make_released_dofs1:
            no_shear_y = np.zeros(self.element_ids.shape, dtype='int32')
            no_shear_z = np.zeros(self.element_ids.shape, dtype='int32')
            no_bending_y = np.zeros(self.element_ids.shape, dtype='int32')
            no_bending_z = np.zeros(self.element_ids.shape, dtype='int32')

        if self.make_released_dofs2:
            no_bending = np.zeros(self.element_ids.shape, dtype='int32')
            no_bending_bad = np.zeros(self.element_ids.shape, dtype='int32')
            no_6_16 = np.zeros(self.element_ids.shape, dtype='int32')
            no_0_56 = np.zeros(self.element_ids.shape, dtype='int32')
            no_0_456 = np.zeros(self.element_ids.shape, dtype='int32')
            no_56_456 = np.zeros(self.element_ids.shape, dtype='int32')
            no_0_6 = np.zeros(self.element_ids.shape, dtype='int32')
            no_0_16 = np.zeros(self.element_ids.shape, dtype='int32')

        #debug = True
        bar_nids = set([])
        for eid in bar_beam_eids:
            if eid not in self.eid_map:
                self.log.error('eid=%s is not a valid element...' % eid)
                if debug:
                    print('eid=%s is not a valid element...' % eid)
                continue
            ieid = self.eid_map[eid]
            elem = model.elements[eid]
            pid = elem.pid
            assert not isinstance(pid, integer_types), elem
            if pid.type in ['PBAR', 'PBEAM']:
                bar_type = 'bar'
            elif pid.type in ['PBEAM']:
                bar_type = 'beam'
            elif pid.type in ['PBARL', 'PBEAML']:
                bar_type = pid.Type
            else:
                if debug:
                    print('NotImplementedError(pid)')
                raise NotImplementedError(pid)
            #print('bar_type =', bar_type)

            if debug:
                print('%s' % elem)
                print('  bar_type =', bar_type)
            found_bar_types.add(bar_type)

            (nid1, nid2) = elem.node_ids
            bar_nids.update([nid1, nid2])
            node1 = model.nodes[nid1]
            node2 = model.nodes[nid2]
            n1 = node1.get_position()
            n2 = node2.get_position()
            centroid = (n1 + n2) / 2.
            i = n2 - n1
            Li = norm(i)
            ihat = i / Li

            if not self.make_released_dofs2:
                if elem.pa == 0 and elem.pb == 0:
                    pass

                if elem.pa == 1 or elem.pb == 1:
                    no_axial[ieid] = 1
                if elem.pa == 2 or elem.pb == 2:
                    no_axial[ieid] = 1
                if elem.pa == 3 or elem.pb == 3:
                    no_axial[ieid] = 1
                if elem.pa == 4 or elem.pb == 4:
                    no_torsion[ieid] = 1
                if elem.pa == 5 or elem.pb == 5:
                    no_axial[ieid] = 1
                if elem.pa == 6 or elem.pb == 6:
                    no_axial[ieid] = 1

            else:
                if elem.pa == 0 and elem.pb == 0:
                    pass
                elif (elem.pa == 6 and elem.pb == 16) or (elem.pa == 16 and elem.pb == 6):
                    no_axial[ieid] = 1
                    no_6_16[ieid] = 1
                elif (elem.pa == 56 and elem.pb == 0) or (elem.pa == 0 and elem.pb == 56):
                    no_bending[ieid] = 1
                    no_0_56[ieid] = 1
                    #print(elem)
                elif (elem.pa == 0 and elem.pb == 456) or (elem.pa == 456 and elem.pb == 0):
                    no_bending[ieid] = 1
                    no_torsion[ieid] = 1
                    no_0_456[ieid] = 1
                    #print(elem)
                elif (elem.pa == 456 and elem.pb == 56) or (elem.pa == 56 and elem.pb == 456):
                    no_torsion[ieid] = 1
                    no_56_456[ieid] = 1
                elif elem.pa == 6 and elem.pb == 0:
                    no_bending_bad[ieid] = 1
                    no_0_6[ieid] = 1
                    #print(elem)
                elif elem.pa == 0 and elem.pb == 16 or elem.pb == 0 and elem.pa == 16:
                    no_axial[ieid] = 1
                    no_bending_bad[ieid] = 1
                    # print(elem)
                    no_0_16[ieid] = 1
                elif elem.pa == 56 and elem.pb == 45 or elem.pb == 56 and elem.pa == 45:
                    no_torsion[ieid] = 1
                    no_bending[ieid] = 1
                else:
                    msg = 'pa=%r pb=%r; elem=\n%s' % (elem.pa, elem.pb, elem)
                    print(msg)
                    continue
                    #raise NotImplementedError(msg)


            # OFFT flag
            # ---------
            # ABC or A-B-C (an example is G-G-G or B-G-G)
            # while the slots are:
            #  - A -> orientation; values=[G, B]
            #  - B -> End A; values=[G, O]
            #  - C -> End B; values=[G, O]
            #
            # and the values for A,B,C mean:
            #  - B -> basic
            #  - G -> global
            #  - O -> orientation
            #
            # so for example G-G-G, that's global for all terms.
            # BOG means basic orientation, orientation end A, global end B
            #
            # so now we're left with what does basic/global/orientation mean?
            # - basic -> the glboal coordinate system defined by cid=0
            # - global -> the local coordinate system defined by the
            #             CD field on the GRID card, but referenced by
            #             the CBAR/CBEAM
            # - orientation -> ???
            #
            if elem.g0:
                #debug = False
                msg = 'which is required by %s eid=%s\n%s' % (elem.type, elem.g0, str(elem))
                g0_ref = model.Node(elem.g0, msg=msg)
                if debug:
                    print('  g0 = %s' % elem.g0)
                    print('  g0_ref = %s' % g0_ref)
                n0 = g0_ref.get_position()
                v = n0 - n1
            else:
                #debug = False
                ga = model.nodes[elem.Ga()]
                v = ga.cd_ref.transform_node_to_global(elem.x)
                if debug:
                    print('  ga = %s' % elem.ga)
                    if ga.Cd() != 0:
                        print('  cd = %s' % ga.cd_ref)
                    else:
                        print('  cd = 0')

                    print('  x = %s' % elem.x)
                    print('  v = %s' % v)
                #v = elem.x

            offt_vector, offt_end_a, offt_end_b = elem.offt
            if debug:
                print('  offt vector,A,B=%r' % (elem.offt))
            # if offt_end_a == 'G' or (offt_end_a == 'O' and offt_vector == 'G'):

            if offt_vector == 'G':
                # end A
                # global - cid != 0
                if node1.Cd() != 0:
                    v = node1.cd_ref.transform_node_to_global(v)
                    if node1.cd_ref.type not in ['CORD2R', 'CORD1R']:
                        msg = 'invalid Cd type (%r) on Node %i; expected CORDxR' % (
                            node1.cd_ref.type, node1.nid)
                        self.log.error(msg)
                        continue
                        #raise NotImplementedError(node1.cd)
            elif offt_vector == 'B':
                # basic - cid = 0
                pass
            else:
                msg = 'offt_vector=%r is not supported; offt=%s' % (offt_vector, elem.offt)
                self.log.error(msg)
                continue
                #raise NotImplementedError(msg)
            #print('v =', v)

            # rotate wa
            wa = elem.wa
            if offt_end_a == 'G':
                if node1.Cd() != 0:
                    if node1.cd.type not in ['CORD2R', 'CORD1R']:
                        continue # TODO: support CD transform
                        #raise NotImplementedError(node1.cd)
                    wa = node1.cd.transform_node_to_global(wa)  # TODO: fixme
            elif offt_end_a == 'B':
                pass
            elif offt_end_a == 'O':
                wa = node1.cd.transform_node_to_global(n1 - wa)  # TODO: fixme
            else:
                msg = 'offt_end_a=%r is not supported; offt=%s' % (offt_end_a, elem.offt)
                self.log.error(msg)
                continue
                #raise NotImplementedError(msg)

            #print('wa =', wa)
            # rotate wb
            wb = elem.wb
            if offt_end_b == 'G':
                if node2.Cd() != 0:
                    if node2.cd.type not in ['CORD2R', 'CORD1R']:
                        continue # TODO: MasterModelTaxi
                        #raise NotImplementedError(node2.cd)
                    #wb = node2.cd.transform_node_to_local_xyz(wb)
                    wb = node2.cd.transform_node_to_global(wb)  # TODO: fixme

            elif offt_end_b == 'B':
                pass
            elif offt_end_b == 'O':
                wb = node1.cd.transform_node_to_global(n2 - wb)  # TODO: fixme
            else:
                msg = 'offt_end_b=%r is not supported; offt=%s' % (offt_end_b, elem.offt)
                model.log.error(msg)
                continue
                #raise NotImplementedError(msg)

            #print('wb =', wb)
            ## concept has a GOO
            #if not elem.offt in ['GGG', 'BGG']:
                #msg = 'offt=%r for CBAR/CBEAM eid=%s is not supported...skipping' % (elem.offt, eid)
                #self.log.error(msg)
                #continue

            vhat = v / norm(v) # j
            try:
                z = np.cross(ihat, vhat) # k
            except ValueError:
                msg = 'Invalid vector length\n'
                msg += 'n1  =%s\n' % str(n1)
                msg += 'n2  =%s\n' % str(n2)
                msg += 'nid1=%s\n' % str(nid1)
                msg += 'nid2=%s\n' % str(nid2)
                msg += 'i   =%s\n' % str(i)
                msg += 'Li  =%s\n' % str(Li)
                msg += 'ihat=%s\n' % str(ihat)
                msg += 'v   =%s\n' % str(v)
                msg += 'vhat=%s\n' % str(vhat)
                msg += 'z=cross(ihat, vhat)'
                print(msg)
                raise ValueError(msg)

            zhat = z / norm(z)
            yhat = np.cross(zhat, ihat) # j
            if debug:
                print('  centroid = %s' % centroid)
                print('  ihat = %s' % ihat)
                print('  yhat = %s' % yhat)
                print('  zhat = %s' % zhat)
                print('  scale = %s' % scale)
            #if eid == 5570:
                #print('  check - eid=%s yhat=%s zhat=%s v=%s i=%s n%s=%s n%s=%s' % (
                      #eid, yhat, zhat, v, i, nid1, n1, nid2, n2))

            if norm(ihat) == 0.0 or norm(yhat) == 0.0 or norm(z) == 0.0:
                print('  invalid_orientation - eid=%s yhat=%s zhat=%s v=%s i=%s n%s=%s n%s=%s' % (
                    eid, yhat, zhat, v, i, nid1, n1, nid2, n2))
            elif not np.allclose(norm(yhat), 1.0) or not np.allclose(norm(zhat), 1.0) or Li == 0.0:
                print('  length_error        - eid=%s Li=%s Lyhat=%s Lzhat=%s'
                      ' v=%s i=%s n%s=%s n%s=%s' % (
                          eid, Li, norm(yhat), norm(zhat), v, i, nid1, n1, nid2, n2))

            #print('adding bar %s' % bar_type)
            #print('   centroid=%s' % centroid)
            #print('   yhat=%s len=%s' % (yhat, np.linalg.norm(yhat)))
            #print('   zhat=%s len=%s' % (zhat, np.linalg.norm(zhat)))
            #print('   Li=%s scale=%s' % (Li, scale))
            if bar_type not in allowed_types:
                msg = 'bar_type=%r allowed=[%s]' % (bar_type, ', '.join(allowed_types))
                raise RuntimeError(msg)

            bar_types[bar_type][0].append(eid)
            bar_types[bar_type][1].append((centroid, centroid + yhat * Li * scale))
            bar_types[bar_type][2].append((centroid, centroid + zhat * Li * scale))
            #if len(bar_types[bar_type][0]) > 5:
                #break

        #print('bar_types =', bar_types)
        for bar_type in bar_types.keys():
            if len(bar_types[bar_type][0]) == 0:
                del bar_types[bar_type]

        debug = False
        if debug:
            #np.set_printoptions(formatter={'float': '{: 0.3f}'.format})
            for bar_type, data in sorted(iteritems(bar_types)):
                eids, lines_bar_y, lines_bar_z = data
                if len(eids):
                    #print('barsi =', barsi)
                    #print('bar_type = %r' % bar_type)
                    for eid, line_y, line_z  in zip(eids, lines_bar_y, lines_bar_z):
                        print('eid=%s centroid=%s cy=%s cz=%s' % (
                            eid, line_y[0], line_y[1], line_z[1]))

        #print('found_bar_types =', found_bar_types)

        bar_nids = list(bar_nids)
        red = (1., 0., 0.)
        self.create_alternate_vtk_grid(
            'Bar Nodes', color=red, line_width=1, opacity=1.,
            point_size=5, representation='point', bar_scale=0., is_visible=True)
        self._add_nastran_nodes_to_grid('Bar Nodes', bar_nids, model)


        geo_form = form[2]
        bar_form = ('CBAR / CBEAM', None, [])
        #print('geo_form =', geo_form)
        bar_types2 = {}
        for bar_type, data in sorted(iteritems(bar_types)):
            eids, lines_bar_y, lines_bar_z = data
            #print(data)
            if len(eids):
                if debug:
                    print('bar_type = %r' % bar_type)
                # if bar_type not in ['ROD', 'TUBE']:
                bar_y = bar_type + '_y'
                bar_z = bar_type + '_z'

                self.create_alternate_vtk_grid(
                    bar_y, color=green, line_width=5, opacity=1.,
                    point_size=5, representation='bar', bar_scale=scale, is_visible=False)
                self.create_alternate_vtk_grid(
                    bar_z, color=blue, line_width=5, opacity=1.,
                    point_size=5, representation='bar', bar_scale=scale, is_visible=False)

                self._add_nastran_lines_xyz_to_grid(bar_y, lines_bar_y, eids, model)
                self._add_nastran_lines_xyz_to_grid(bar_z, lines_bar_z, eids, model)

                # form = ['Geometry', None, []]
                i = np.searchsorted(self.element_ids, eids)
                is_type = np.zeros(self.element_ids.shape, dtype='int32')
                is_type[i] = 1.
                # print('is-type =', is_type.max())
                bar_form[2].append(['is_%s' % bar_type, icase, []])

                msg = 'is_%s' % bar_type
                type_res = GuiResult(0, header=msg, title=msg,
                                     location='centroid', scalar=is_type)
                cases[icase] = (type_res, (0, msg))
                #cases[(0, icase, msg, 1, 'centroid', '%i', '')] = is_type
                icase += 1

        if self.make_released_dofs2:
            if no_axial.max() == 1:
                bar_form[2].append(['No Axial', icase, []])
                axial_res = GuiResult(0, header='No Axial', title='No Axial',
                                      location='centroid', scalar=no_axial)
                cases[icase] = (axial_res, (0, 'No Axial'))
                icase += 1

            if no_torsion.max() == 1:
                bar_form[2].append(['No Torsion', icase, []])
                torsion_res = GuiResult(0, header='No Torsion', title='No Torsion',
                                        location='centroid', scalar=no_torsion)
                cases[icase] = (torsion_res, (0, 'No Torsion'))
                icase += 1

        if self.make_released_dofs1:
            if no_shear_y.max() == 1:
                bar_form[2].append(['No Shear Y', icase, []])
                shear_y_res = GuiResult(0, header='No Shear Y', title='No Shear Y',
                                        location='centroid', scalar=no_shear_y)
                cases[icase] = (shear_y_res, (0, 'No Shear Y'))
                #cases[(0, icase, 'No Shear Y', 1, 'centroid', '%i', '')] = no_shear_y
                icase += 1
            if no_shear_z.max() == 1:
                bar_form[2].append(['No Shear Z', icase, []])
                shear_z_res = GuiResult(0, header='No Shear Z', title='No Shear Z',
                                        location='centroid', scalar=no_shear_y)
                cases[icase] = (shear_z_res, (0, 'No Shear Z'))
                #cases[(0, icase, 'No Shear Z', 1, 'centroid', '%i', '')] = no_shear_z
                icase += 1
            if no_bending_y.max() == 1:
                bar_form[2].append(['No Bending Y', icase, []])
                bending_y_res = GuiResult(0, header='No Bending Z', title='No Bending Z',
                                          location='centroid', scalar=no_bending_y)
                cases[icase] = (bending_y_res, (0, 'No Bending Z'))
                #cases[(0, icase, 'No Bending Y', 1, 'centroid', '%i', '')] = no_bending_y
                icase += 1
            if no_bending_z.max() == 1:
                bar_form[2].append(['No Bending Z', icase, []])
                bending_z_res = GuiResult(0, header='No Bending Z', title='No Bending Z',
                                          location='centroid', scalar=no_bending_z)
                cases[icase] = (bending_z_res, (0, 'No Bending Z'))
                #cases[(0, icase, 'No Bending Z', 1, 'centroid', '%i', '')] = no_bending_z
                icase += 1

        if self.make_released_dofs2:
            if no_bending.max() == 1:
                bar_form[2].append(['No Bending', icase, []])
                bending_res = GuiResult(0, header='No Bending', title='No Bending',
                                        location='centroid', scalar=no_bending)
                cases[icase] = (bending_res, (0, 'No Bending'))
                #cases[(0, icase, 'No Bending', 1, 'centroid', '%i', '')] = no_bending
                icase += 1

            if no_bending_bad.max() == 1:
                bar_form[2].append(['No Bending (Bad)', icase, []])
                cases[(0, icase, 'No Bending (Bad)', 1, 'centroid', '%i', '')] = no_bending_bad
                icase += 1

            if no_6_16.max() == 1:
                bar_form[2].append(['no_6_16', icase, []])
                cases[(0, icase, 'no_6_16', 1, 'centroid', '%i', '')] = no_6_16
                icase += 1
            if no_0_56.max() == 1:
                bar_form[2].append(['no_0_56', icase, []])
                cases[(0, icase, 'no_0_56', 1, 'centroid', '%i', '')] = no_0_56
                icase += 1
            if no_0_456.max() == 1:
                bar_form[2].append(['no_0_456', icase, []])
                cases[(0, icase, 'no_0_456', 1, 'centroid', '%i', '')] = no_0_456
                icase += 1
            if no_56_456.max() == 1:
                bar_form[2].append(['no_56_456', icase, []])
                cases[(0, icase, 'no_56_456', 1, 'centroid', '%i', '')] = no_56_456
                icase += 1
            if no_0_6.max() == 1:
                bar_form[2].append(['no_0_6', icase, []])
                cases[(0, icase, 'no_0_6', 1, 'centroid', '%i', '')] = no_0_6
                icase += 1
            if no_0_16.max() == 1:
                bar_form[2].append(['no_0_16)', icase, []])
                cases[(0, icase, 'no_0_16', 1, 'centroid', '%i', '')] = no_0_16
                icase += 1

        # print(geo_form)
        if len(bar_form[2]):
            geo_form.append(bar_form)
        return icase

    def _add_nastran_lines_xyz_to_grid(self, name, lines, eids, model):
        """creates the bar orientation vector lines"""
        nlines = len(lines)
        nnodes = nlines * 2
        if nlines == 0:
            return
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(nnodes)

        assert name != u'Bar Nodes', name

        bar_eids = np.array(eids, dtype='int32')
        bar_lines = np.zeros((nlines, 6))
        self.bar_lines[name] = bar_lines
        self.bar_eids[name] = np.asarray(eids, dtype='int32')

        j = 0
        for i, (node1, node2) in enumerate(lines):
            bar_lines[i, :3] = node1
            bar_lines[i, 3:] = node2

            points.InsertPoint(j, *node1)
            points.InsertPoint(j + 1, *node2)
            # print('adding %s %s' % (str(node1), str(node2)))

            elem = vtk.vtkLine()
            elem.GetPointIds().SetId(0, j)
            elem.GetPointIds().SetId(1, j + 1)
            self.alt_grids[name].InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            j += 2

        #n1 = bar_lines[:, :3]
        #n2 = bar_lines[:, 3:]
        #dy = n2 - n1
        #Ly = norm(dy, axis=1)
        # v = dy / Ly *  bar_scale
        # n2 = n1 + v
        # print(Ly)
        self.alt_grids[name].SetPoints(points)


    def _get_rigid(self, dim_max, model):
        """
        dependent = (lines[:, 0])
        independent = np.unique(lines[:, 1])
        """
        lines_rigid = []
        for eid, elem in iteritems(model.rigid_elements):
            if elem.type == 'RBE3':
                if elem.Gmi != []:
                    msg = 'UM is not supported; RBE3 eid=%s Gmi=%s' % (elem.eid, elem.Gmi)
                    raise RuntimeError(msg)
                #list_fields = ['RBE3', elem.eid, None, elem.ref_grid_id, elem.refc]
                n1 = elem.ref_grid_id
                assert isinstance(n1, int), 'RBE3 eid=%s ref_grid_id=%s' % (elem.eid, n1)
                for (_weight, ci, Gij) in elem.WtCG_groups:
                    Giji = elem._nodeIDs(nodes=Gij, allow_empty_nodes=True)
                    # list_fields += [wt, ci] + Giji
                    for n2 in Giji:
                        assert isinstance(n2, int), 'RBE3 eid=%s Giji=%s' % (elem.eid, Giji)
                        lines_rigid.append([n1, n2])
            elif elem.type == 'RBE2':
                #list_fields = ['RBE2', elem.eid, elem.Gn(), elem.cm
                               #] + elem.Gmi_node_ids + [elem.alpha]
                n2 = elem.Gn() # independent
                nids1 = elem.Gmi_node_ids # dependent
                for n1 in nids1:
                    lines_rigid.append([n1, n2])
            elif elem.type in ['RBAR', 'RBAR1', 'RROD']:
                dependent = elem.Ga()
                independent = elem.Gb()
                lines_rigid.append([dependent, independent])
            else:
                print(str(elem))
        return lines_rigid

    def _fill_mpc(self, mpc_id, dim_max, model, nid_to_pid_map):
        """helper for making MPCs"""
        lines = self.get_MPCx_node_ids_c1(model, mpc_id, exclude_mpcadd=False)
        lines += self._get_rigid(dim_max, model)
        self._fill_dependent_independent(dim_max, model, lines, nid_to_pid_map)

    def _fill_rigid(self, dim_max, model, nid_to_pid_map):
        """helper for making rigid elements"""
        lines = self._get_rigid(dim_max, model)
        self._fill_dependent_independent(dim_max, model, lines, nid_to_pid_map)

    def _fill_dependent_independent(self, dim_max, model, lines, nid_to_pid_map):
        if not lines:
            return
        green = (0., 1., 0.)
        dunno = (0.5, 1., 0.5)
        self.create_alternate_vtk_grid(
            'mpc_dependent', color=green, line_width=5, opacity=1.,
            point_size=5, representation='point', is_visible=False)
        self.create_alternate_vtk_grid(
            'mpc_independent', color=dunno, line_width=5, opacity=1.,
            point_size=5, representation='point', is_visible=False)
        self.create_alternate_vtk_grid(
            'mpc_lines', color=dunno, line_width=5, opacity=1.,
            point_size=5, representation='wire', is_visible=False)

        lines2 = []
        for line in lines:
            if line not in lines2:
                lines2.append(line)
        lines = np.array(lines2, dtype='int32')
        dependent = (lines[:, 0])
        independent = np.unique(lines[:, 1])
        self.dependents_nodes.update(dependent)
        node_ids = np.unique(lines.ravel())
        self._add_nastran_nodes_to_grid('mpc_dependent', dependent, model, nid_to_pid_map)
        self._add_nastran_nodes_to_grid('mpc_independent', independent, model, nid_to_pid_map)
        self._add_nastran_lines_to_grid('mpc_lines', lines, model, nid_to_pid_map)

    def _add_nastran_nodes_to_grid(self, name, node_ids, model, nid_to_pid_map=None):
        """used to create MPC independent/dependent nodes"""
        nnodes = len(node_ids)
        if nnodes == 0:
            model.log.warning('0 nodes added for %r' % name)
            return
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(nnodes)

        j = 0
        nid_map = self.nid_map
        for nid in sorted(node_ids):
            try:
                i = nid_map[nid]
            except KeyError:
                model.log.warning('nid=%s does not exist' % nid)
                continue

            if nid not in model.nodes:
                model.log.warning('nid=%s doesnt exist' % nid)
                continue
            # point = self.grid.GetPoint(i)
            # points.InsertPoint(j, *point)

            node = model.nodes[nid]
            point = node.get_position()
            #self.log_info('adding SUPORT1; p=%s' % str(point))
            points.InsertPoint(j, *point)

            #if 1:
            elem = vtk.vtkVertex()
            elem.GetPointIds().SetId(0, j)
            #else:
                #elem = vtk.vtkSphere()
                #dim_max = 1.0
                #sphere_size = self._get_sphere_size(dim_max)
                #elem.SetRadius(sphere_size)
                #elem.SetCenter(points.GetPoint(j))

            self.alt_grids[name].InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            j += 1
        self.alt_grids[name].SetPoints(points)

    def _add_nastran_lines_to_grid(self, name, lines, model, nid_to_pid_map=None):
        """used to create MPC lines"""
        nlines = lines.shape[0]
        #nids = np.unique(lines)
        #nnodes = len(nids)
        nnodes = nlines * 2
        if nnodes == 0:
            return
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(nnodes)

        j = 0
        nid_map = self.nid_map
        for nid1, nid2 in lines:
            try:
                i1 = nid_map[nid1]
            except KeyError:
                model.log.warning('nid=%s does not exist' % nid1)
                continue
            try:
                i2 = nid_map[nid2]
            except KeyError:
                model.log.warning('nid=%s does not exist' % nid2)
                continue

            if nid1 not in model.nodes:
                continue
            if nid2 not in model.nodes:
                continue
            node = model.nodes[nid1]
            point = node.get_position()
            points.InsertPoint(j, *point)

            node = model.nodes[nid2]
            point = node.get_position()
            points.InsertPoint(j + 1, *point)

            elem = vtk.vtkLine()
            elem.GetPointIds().SetId(0, j)
            elem.GetPointIds().SetId(1, j + 1)
            self.alt_grids[name].InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            j += 2
        self.alt_grids[name].SetPoints(points)

    def _fill_suport(self, suport_id, dim_max, model):
        """creates SUPORT and SUPORT1 nodes"""
        #pink = (0.98, 0.4, 0.93)
        red = (1.0, 0., 0.)
        self.create_alternate_vtk_grid(
            'suport', color=red, line_width=5, opacity=1., point_size=4,
            representation='point', is_visible=False)

        node_ids = []

        # list
        #for suport in model.suport:
            #node_ids += suport.IDs

        # dict
        if suport_id in model.suport1:
            suport1 = model.suport1[suport_id]
            node_ids += suport1.IDs
        else:
            for suport in model.suport:
                if suport_id in suport.IDs:
                    node_ids.append(suport_id)

        node_ids = np.unique(node_ids)
        self._add_nastran_nodes_to_grid('suport', node_ids, model)

    def _get_sphere_size(self, dim_max):
        return 0.01 * dim_max

    def map_elements(self, points, nid_map, model, j, dim_max,
                     plot=True, xref_loads=True):
        """
        element checks
        http://www.altairuniversity.com/wp-content/uploads/2012/04/Student_Guide_211-233.pdf

        Skew:
          Skew in trias is calculated by finding the minimum angle
          between the vector from each node to the opposing mid-side
          and the vector between the two adjacent mid-sides at each
          node of the element.  Ninety degrees minus the minimum angle
          found is reported.

          Skew in quads is calculated by finding the minimum angle
          between two lines joining opposite midsides of the element.
          Ninety degrees minus the minimum angle found is reported.

        Aspect Ratio:
          Aspect ratio in two-dimensional elements is calculated by
          dividing the maximum length side of an element by the minimum
          length side of the element. The aspect ratio check is
          performed in the same fashion on all faces of 3D elements.

        Warpage:
          Warpage in two-dimensional elements is calculated by splitting
          a quad into two trias and finding the angle between the two
          planes which the trias form. The quad is then split again,
          this time using the opposite corners and forming the second
          set of trias. The angle between the two planes which the trias
          form is then found. The maximum angle found between the planes
          is the warpage of the element.

          Warpage in three-dimensional elements is performed in the same
          fashion on all faces of the element.

        Jacobian:
          determinant of Jacobian matrix (-1.0 to 1.0; 1.0 is ideal)

        2D Checks:

        Warp angle:
          Warp angle is the out of plane angle
          Ideal value = 0 degrees (Acceptable < 100).
          Warp angle is not applicable for triangular elements.
          It is defined as the angle between the normals to two planes
          formed by splitting the quad element along the diagonals.
          The maximum angle of the two possible angles is reported as
          the warp angle.

        Aspect:
          Aspect = maximum element edge length / minimum element edge length
          Ideal value = 1 (Acceptable < 5).

        Skew:
          Ideal value = 0 degrees (Acceptable < 45)
          Skew for quadrilateral element = 90
          minus the minimum angle between the two lines joining the
          opposite mid-sides of the element (alpha).

          Skew for triangular element = 90
          minus the minimum angle between the lines from each node to
          the opposing mid-side and between the two adjacent mid-sides
          at each node of the element

        Jacobian:
          Ideal value = 1.0 (Acceptable > 0.6)
          In simple terms, the jacobian is a scale factor arising
          because of the transformation of the coordinate system.
          Elements are tansformed from the global coordinates to
          local coordinates (defined at the centroid of every
          element), for faster analysis times.

        Distortion:
          Ideal value = 1.0 (Acceptable > 0.6)
          Distortion is defined as:
             d = |Jacobian| * AreaLCS / AreaGCS
          LCS - Local Coordinate system
          GCS - Global Coordinate system

        Stretch:
          Ideal value: 1.0 (Acceptable > 0.2)
          For quadrilateral elements stretch = Lmin * sqrt(2) / dmax
          Stretch for triangular element = R * sqrt(12) / Lmax

        Included angles:
          Skew is based on the overall shape of the element and it does
          not take into account the individual angles of a quadrilateral
          or triangular element.  Included or interior angle check is
          applied for individual angles.
          Quad: Ideal value = 90 (Acceptable = 45 < theta <135)
          Tria: Ideal value = 60 (Acceptable = 20 < theta < 120)

        Taper:
          Ideal value = 0 (Acceptable < 0.5)
          Taper = sum( | (Ai - Aavg) / Aavg |)
          Aavg = (A1 + A2 + A3 + A4) / 4
          A1,A2 are one split form of the CQUAD4 and A3,A4 are the quad
          split in the other direction.
        """
        xyz_cid0 = self.xyz_cid0
        sphere_size = self._get_sphere_size(dim_max)

        # :param i: the element id in grid
        # :param j: the element id in grid2
        i = 0

        #nids = self.eid_to_nid_map[eid]
        self.eid_to_nid_map = {}

        # the list of all pids
        #pids = []

        # pid = pids_dict[eid]
        pids_dict = {}
        nelements = len(model.elements)
        pids = np.zeros(nelements, 'int32')
        mids = np.zeros(nelements, 'int32')
        material_coord = np.zeros(nelements, 'int32')
        min_interior_angle = np.zeros(nelements, 'float32')
        max_interior_angle = np.zeros(nelements, 'float32')
        dideal_theta = np.zeros(nelements, 'float32')
        max_skew_angle = np.zeros(nelements, 'float32')
        max_warp_angle = np.zeros(nelements, 'float32')
        max_aspect_ratio = np.zeros(nelements, 'float32')
        area = np.zeros(nelements, 'float32')
        area_ratio = np.zeros(nelements, 'float32')
        taper_ratio = np.zeros(nelements, 'float32')

        # pids_good = []
        # pids_to_keep = []
        # pids_btm = []
        # pids_to_drop = []

        # 3
        # | \
        # |   \
        # |     \
        # 1------2


        # these normals point inwards
        #      4
        #    / | \
        #   /  |  \
        #  3-------2
        #   \  |   /
        #    \ | /
        #      1
        _ctetra_faces = (
            (0, 1, 2), # (1, 2, 3),
            (0, 3, 1), # (1, 4, 2),
            (0, 3, 2), # (1, 4, 3),
            (1, 3, 2), # (2, 4, 3),
        )

        # these normals point inwards
        #
        #
        #
        #
        #        /4-----3
        #       /       /
        #      /  5    /
        #    /    \   /
        #   /      \ /
        # 1---------2
        _cpyram_faces = (
            (0, 1, 2, 3), # (1, 2, 3, 4),
            (1, 4, 2), # (2, 5, 3),
            (2, 4, 3), # (3, 5, 4),
            (0, 3, 4), # (1, 4, 5),
            (0, 4, 1), # (1, 5, 2),
        )

        # these normals point inwards
        #       /6
        #     /  | \
        #   /    |   \
        # 3\     |     \
        # |  \   /4-----5
        # |    \/       /
        # |   /  \     /
        # |  /    \   /
        # | /      \ /
        # 1---------2
        _cpenta_faces = (
            (0, 2, 1), # (1, 3, 2),
            (3, 4, 5), # (4, 5, 6),

            (0, 1, 4, 3), # (1, 2, 5, 4), # bottom
            (1, 2, 5, 4), # (2, 3, 6, 5), # right
            (0, 3, 5, 2), # (1, 4, 6, 3), # left
        )

        # these normals point inwards
        #      8----7
        #     /|   /|
        #    / |  / |
        #   /  5-/--6
        # 4-----3   /
        # |  /  |  /
        # | /   | /
        # 1-----2
        _chexa_faces = (
            (4, 5, 6, 7), # (5, 6, 7, 8),
            (0, 3, 2, 1), # (1, 4, 3, 2),
            (1, 2, 6, 5), # (2, 3, 7, 6),
            (2, 3, 7, 6), # (3, 4, 8, 7),
            (0, 4, 7, 3), # (1, 5, 8, 4),
            (0, 6, 5, 4), # (1, 7, 6, 5),
        )
        nid_to_pid_map = defaultdict(list)
        for (eid, element) in sorted(iteritems(model.elements)):
            # if element.Pid() >= 82:
                # continue
            # if element.Pid() in pids_to_drop:
                # continue
            # if element.Pid() not in pids_to_keep:
                # continue
            # if element.pid.type == 'PSOLID':
                # continue
            self.eid_map[eid] = i
            pid = 0
            dideal_thetai = 0.0
            min_thetai = 0.0
            max_thetai = 0.0
            max_skew = 0.0
            max_warp = 0.0
            aspect_ratio = 1.0
            areai = 0.
            area_ratioi = 1.
            taper_ratioi = 0.
            if isinstance(element, (CTRIA3, CTRIAR)):
                material_coord[i] = 0 if isinstance(element.thetaMcid, float) else element.thetaMcid
                elem = vtkTriangle()
                node_ids = element.node_ids
                pid = element.Pid()
                self.eid_to_nid_map[eid] = node_ids
                for nid in node_ids:
                    if nid is not None:
                        nid_to_pid_map[nid].append(pid)

                n1, n2, n3 = [nid_map[nid] for nid in node_ids]
                p1 = xyz_cid0[n1, :]
                p2 = xyz_cid0[n2, :]
                p3 = xyz_cid0[n3, :]
                e1 = (p1 + p2) / 2.
                e2 = (p2 + p3) / 2.
                e3 = (p3 + p1) / 2.

                #    3
                #    / \
                # e3/   \ e2
                #  /    /\
                # /    /  \
                # 1---/----2
                #    e1
                e21 = e2 - e1
                e31 = e3 - e1
                e32 = e3 - e2

                e3_p2 = e3 - p2
                e2_p1 = e2 - p1
                e1_p3 = e1 - p3

                v21 = p2 - p1
                v32 = p3 - p2
                v13 = p1 - p3
                areai = 0.5 * np.linalg.norm(np.cross(v21, v13))

                cos_skew1 = np.dot(e2_p1, e31) / (np.linalg.norm(e2_p1) * np.linalg.norm(e31))
                cos_skew2 = np.dot(e2_p1, -e31) / (np.linalg.norm(e2_p1) * np.linalg.norm(e31))
                cos_skew3 = np.dot(e3_p2, e21) / (np.linalg.norm(e3_p2) * np.linalg.norm(e21))
                cos_skew4 = np.dot(e3_p2, -e21) / (np.linalg.norm(e3_p2) * np.linalg.norm(e21))
                cos_skew5 = np.dot(e1_p3, e32) / (np.linalg.norm(e1_p3) * np.linalg.norm(e32))
                cos_skew6 = np.dot(e1_p3, -e32) / (np.linalg.norm(e1_p3) * np.linalg.norm(e32))
                max_skew = np.pi / 2. - np.abs(np.arccos(np.clip([
                    cos_skew1, cos_skew2, cos_skew3,
                    cos_skew4, cos_skew5, cos_skew6], -1., 1.))).min()
                lengths = np.linalg.norm([v21, v32, v13], axis=1)
                #assert len(lengths) == 3, lengths
                aspect_ratio = lengths.max() / lengths.min()

                cos_theta1 = np.dot(v21, -v13) / (np.linalg.norm(v21) * np.linalg.norm(v13))
                cos_theta2 = np.dot(v32, -v21) / (np.linalg.norm(v32) * np.linalg.norm(v21))
                cos_theta3 = np.dot(v13, -v32) / (np.linalg.norm(v13) * np.linalg.norm(v32))
                thetas = np.arccos(np.clip([cos_theta1, cos_theta2, cos_theta3], -1., 1.))
                min_thetai = thetas.min()
                max_thetai = thetas.max()
                dideal_thetai = max(max_thetai - piover3, piover3 - min_thetai)

                #theta_deg = np.degrees(np.arccos(max_cos_theta))
                #if theta_deg < 60.:
                    #print('p1=%s' % xyz_cid0[p1, :])
                    #print('p2=%s' % xyz_cid0[p2, :])
                    #print('p3=%s' % xyz_cid0[p3, :])
                    #print('theta1=%s' % np.degrees(np.arccos(cos_theta1)))
                    #print('theta2=%s' % np.degrees(np.arccos(cos_theta2)))
                    #print('theta3=%s' % np.degrees(np.arccos(cos_theta3)))
                    #print('max_theta=%s' % theta_deg)
                    #asdf

                elem.GetPointIds().SetId(0, n1)
                elem.GetPointIds().SetId(1, n2)
                elem.GetPointIds().SetId(2, n3)
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element, CTRIA6):
                material_coord[i] = element.thetaMcid
                node_ids = element.node_ids
                pid = element.Pid()
                self.eid_to_nid_map[eid] = node_ids[:3]
                for nid in node_ids:
                    if nid is not None:
                        nid_to_pid_map[nid].append(pid)
                if None not in node_ids:
                    elem = vtkQuadraticTriangle()
                    elem.GetPointIds().SetId(3, nid_map[node_ids[3]])
                    elem.GetPointIds().SetId(4, nid_map[node_ids[4]])
                    elem.GetPointIds().SetId(5, nid_map[node_ids[5]])
                else:
                    elem = vtkTriangle()
                elem.GetPointIds().SetId(0, nid_map[node_ids[0]])
                elem.GetPointIds().SetId(1, nid_map[node_ids[1]])
                elem.GetPointIds().SetId(2, nid_map[node_ids[2]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element, CTRIAX6):
                material_coord[i] = element.thetaMcid
                # midside nodes are required, nodes out of order
                node_ids = element.node_ids
                pid = element.Pid()
                for nid in node_ids:
                    if nid is not None:
                        nid_to_pid_map[nid].append(pid)

                if None not in node_ids:
                    elem = vtkQuadraticTriangle()
                    elem.GetPointIds().SetId(3, nid_map[node_ids[1]])
                    elem.GetPointIds().SetId(4, nid_map[node_ids[3]])
                    elem.GetPointIds().SetId(5, nid_map[node_ids[5]])
                else:
                    elem = vtkTriangle()
                elem.GetPointIds().SetId(0, nid_map[node_ids[0]])
                elem.GetPointIds().SetId(1, nid_map[node_ids[2]])
                elem.GetPointIds().SetId(2, nid_map[node_ids[4]])
                self.eid_to_nid_map[eid] = [node_ids[0], node_ids[2], node_ids[4]]
                #a = [0, 2, 4]
                #msg = "CTRIAX6 %i %i %i" %(nid_map[node_ids[a[0]]],
                #                           nid_map[node_ids[a[1]]],
                #                           nid_map[node_ids[a[2]]])
                #raise RuntimeError(msg)
                #sys.stdout.flush()

                #elem.GetPointIds().SetId(0, nid_map[node_ids[0]])
                #elem.GetPointIds().SetId(1, nid_map[node_ids[1]])
                #elem.GetPointIds().SetId(2, nid_map[node_ids[2]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())

            elif isinstance(element, (CQUAD4, CSHEAR, CQUADR)):
                material_coord[i] = 0 if isinstance(element.thetaMcid, float) else element.thetaMcid
                #print('eid=%s theta=%s' % (eid, material_coord[i]))
                node_ids = element.node_ids
                pid = element.Pid()
                for nid in node_ids:
                    if nid is not None:
                        nid_to_pid_map[nid].append(pid)

                self.eid_to_nid_map[eid] = node_ids

                n1, n2, n3, n4 = [nid_map[nid] for nid in node_ids]
                p1 = xyz_cid0[n1, :]
                p2 = xyz_cid0[n2, :]
                p3 = xyz_cid0[n3, :]
                p4 = xyz_cid0[n4, :]

                v21 = p2 - p1
                v32 = p3 - p2
                v43 = p4 - p3
                v14 = p1 - p4

                v42 = p4 - p2
                v31 = p3 - p1
                p12 = (p1 + p2) / 2.
                p23 = (p2 + p3) / 2.
                p34 = (p3 + p4) / 2.
                p14 = (p4 + p1) / 2.
                v31 = p3 - p1
                v42 = p4 - p2
                normal = np.cross(v31, v42)
                areai = 0.5 * np.linalg.norm(normal)

                # still kind of in development
                #
                # the ratio of the ideal area to the actual area
                # this is an hourglass check
                areas = [
                    np.linalg.norm(np.cross(-v14, v21)), # v41 x v21
                    np.linalg.norm(np.cross(v32, -v21)), # v32 x v12
                    np.linalg.norm(np.cross(v43, -v32)), # v43 x v23
                    np.linalg.norm(np.cross(v14, v43)),  # v14 x v43
                ]
                #
                # for:
                #   area=1; area1=0.5 -> area_ratioi1=2.0; area_ratio=2.0
                #   area=1; area1=2.0 -> area_ratioi2=2.0; area_ratio=2.0
                area_ratioi1 = areai / min(areas)
                area_ratioi2 = max(areas) / areai
                area_ratioi = max(area_ratioi1, area_ratioi2)

                area1 = 0.5 * np.linalg.norm(np.cross(-v14, v21)) # v41 x v21
                area2 = 0.5 * np.linalg.norm(np.cross(-v21, v32)) # v12 x v32
                area3 = 0.5 * np.linalg.norm(np.cross(v43, v32)) # v43 x v32
                area4 = 0.5 * np.linalg.norm(np.cross(v14, -v43)) # v14 x v34
                aavg = (area1 + area2 + area3 + area4) / 4.
                taper_ratioi = (abs(area1 - aavg) + abs(area2 - aavg) +
                                abs(area3 - aavg) + abs(area4 - aavg)) / aavg

                #    e3
                # 4-------3
                # |       |
                # |e4     |  e2
                # 1-------2
                #     e1
                e13 = p34 - p12
                e42 = p23 - p14
                cos_skew1 = np.dot(e13, e42) / (np.linalg.norm(e13) * np.linalg.norm(e42))
                cos_skew2 = np.dot(e13, -e42) / (np.linalg.norm(e13) * np.linalg.norm(e42))
                max_skew = np.pi / 2. - np.abs(np.arccos(
                    np.clip([cos_skew1, cos_skew2], -1., 1.))).min()
                #aspect_ratio = max(p12, p23, p34, p14) / max(p12, p23, p34, p14)
                lengths = np.linalg.norm([v21, v32, v43, v14], axis=1)
                #assert len(lengths) == 3, lengths
                aspect_ratio = lengths.max() / lengths.min()

                cos_theta1 = np.dot(v21, -v14) / (np.linalg.norm(v21) * np.linalg.norm(v14))
                cos_theta2 = np.dot(v32, -v21) / (np.linalg.norm(v32) * np.linalg.norm(v21))
                cos_theta3 = np.dot(v43, -v32) / (np.linalg.norm(v43) * np.linalg.norm(v32))
                cos_theta4 = np.dot(v14, -v43) / (np.linalg.norm(v14) * np.linalg.norm(v43))
                #max_thetai = np.arccos([cos_theta1, cos_theta2, cos_theta3, cos_theta4]).max()

                # dot the local normal with the normal vector
                # then take the norm of that to determine the angle relative to the normal
                # then take the sign of that to see if we're pointing roughly towards the normal

                # np.sign(np.linalg.norm(np.dot(
                # a x b = ab sin(theta)
                # a x b / ab = sin(theta)
                # sin(theta) < 0. -> normal is flipped
                normal2 = np.sign(np.dot(np.cross(v21, v32), normal))
                normal3 = np.sign(np.dot(np.cross(v32, v43), normal))
                normal4 = np.sign(np.dot(np.cross(v43, v14), normal))
                normal1 = np.sign(np.dot(np.cross(v14, v21), normal))
                n = np.array([normal1, normal2, normal3, normal4])
                theta_additional = np.where(n < 0, 2*np.pi, 0.)

                theta = n * np.arccos(np.clip(
                    [cos_theta1, cos_theta2, cos_theta3, cos_theta4], -1., 1.)) + theta_additional
                min_thetai = theta.min()
                max_thetai = theta.max()
                dideal_thetai = max(max_thetai - piover2, piover2 - min_thetai)
                #print('theta_max = ', theta_max)

                elem = vtkQuad()
                elem.GetPointIds().SetId(0, n1)
                elem.GetPointIds().SetId(1, n2)
                elem.GetPointIds().SetId(2, n3)
                elem.GetPointIds().SetId(3, n4)
                self.grid.InsertNextCell(9, elem.GetPointIds())

                if 0:
                    # warp
                    v31 = xyz_cid0[p3, :] - xyz_cid0[p1, :]
                    n1a = np.cross(v21, v31) # v21 x v31
                    n1b = np.cross(v31, -v14) # v31 x v41
                    warp1 = np.dot(n1a, n1b) / (np.linalg.norm(n1a) * np.linalg.norm(n1b))

                    v42 = xyz_cid0[p4, :] - xyz_cid0[p2, :]
                    n2a = np.cross(v32, v42) # v32 x v42
                    n2b = np.cross(v42, -v21) # v42 x v12
                    warp2 = np.dot(n2a, n2b) / (np.linalg.norm(n2a) * np.linalg.norm(n2b))
                    max_warp = max(np.arccos(warp1), np.arccos(warp2))

            elif isinstance(element, CQUAD8):
                material_coord[i] = element.thetaMcid
                node_ids = element.node_ids
                pid = element.Pid()
                for nid in node_ids:
                    if nid is not None:
                        nid_to_pid_map[nid].append(pid)
                self.eid_to_nid_map[eid] = node_ids[:4]
                if None not in node_ids:
                    elem = vtkQuadraticQuad()
                    elem.GetPointIds().SetId(4, nid_map[node_ids[4]])
                    elem.GetPointIds().SetId(5, nid_map[node_ids[5]])
                    elem.GetPointIds().SetId(6, nid_map[node_ids[6]])
                    elem.GetPointIds().SetId(7, nid_map[node_ids[7]])
                else:
                    elem = vtkQuad()
                elem.GetPointIds().SetId(0, nid_map[node_ids[0]])
                elem.GetPointIds().SetId(1, nid_map[node_ids[1]])
                elem.GetPointIds().SetId(2, nid_map[node_ids[2]])
                elem.GetPointIds().SetId(3, nid_map[node_ids[3]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element, CTETRA4):
                elem = vtkTetra()
                node_ids = element.node_ids
                pid = element.Pid()
                for nid in node_ids:
                    nid_to_pid_map[nid].append(pid)
                self.eid_to_nid_map[eid] = node_ids[:4]
                elem.GetPointIds().SetId(0, nid_map[node_ids[0]])
                elem.GetPointIds().SetId(1, nid_map[node_ids[1]])
                elem.GetPointIds().SetId(2, nid_map[node_ids[2]])
                elem.GetPointIds().SetId(3, nid_map[node_ids[3]])
                self.grid.InsertNextCell(10, elem.GetPointIds())
                #elem_nid_map = {nid:nid_map[nid] for nid in node_ids[:4]}
                min_thetai, max_thetai, dideal_thetai = get_min_max_theta(
                    _ctetra_faces, node_ids[:4], nid_map, xyz_cid0)
            elif isinstance(element, CTETRA10):
                node_ids = element.node_ids
                pid = element.Pid()
                for nid in node_ids:
                    if nid is not None:
                        nid_to_pid_map[nid].append(pid)
                self.eid_to_nid_map[eid] = node_ids[:4]
                if None not in node_ids:
                    elem = vtkQuadraticTetra()
                    elem.GetPointIds().SetId(4, nid_map[node_ids[4]])
                    elem.GetPointIds().SetId(5, nid_map[node_ids[5]])
                    elem.GetPointIds().SetId(6, nid_map[node_ids[6]])
                    elem.GetPointIds().SetId(7, nid_map[node_ids[7]])
                    elem.GetPointIds().SetId(8, nid_map[node_ids[8]])
                    elem.GetPointIds().SetId(9, nid_map[node_ids[9]])
                else:
                    elem = vtkTetra()
                elem.GetPointIds().SetId(0, nid_map[node_ids[0]])
                elem.GetPointIds().SetId(1, nid_map[node_ids[1]])
                elem.GetPointIds().SetId(2, nid_map[node_ids[2]])
                elem.GetPointIds().SetId(3, nid_map[node_ids[3]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                min_thetai, max_thetai, dideal_thetai = get_min_max_theta(
                    _ctetra_faces, node_ids[:4], nid_map, xyz_cid0)
            elif isinstance(element, CPENTA6):
                elem = vtkWedge()
                node_ids = element.node_ids
                pid = element.Pid()
                for nid in node_ids:
                    nid_to_pid_map[nid].append(pid)
                self.eid_to_nid_map[eid] = node_ids[:6]
                elem.GetPointIds().SetId(0, nid_map[node_ids[0]])
                elem.GetPointIds().SetId(1, nid_map[node_ids[1]])
                elem.GetPointIds().SetId(2, nid_map[node_ids[2]])
                elem.GetPointIds().SetId(3, nid_map[node_ids[3]])
                elem.GetPointIds().SetId(4, nid_map[node_ids[4]])
                elem.GetPointIds().SetId(5, nid_map[node_ids[5]])
                self.grid.InsertNextCell(13, elem.GetPointIds())
                min_thetai, max_thetai, dideal_thetai = get_min_max_theta(
                    _cpenta_faces, node_ids[:6], nid_map, xyz_cid0)

            elif isinstance(element, CPENTA15):
                node_ids = element.node_ids
                pid = element.Pid()
                for nid in node_ids:
                    if nid is not None:
                        nid_to_pid_map[nid].append(pid)
                self.eid_to_nid_map[eid] = node_ids[:6]
                if None not in node_ids:
                    elem = vtkQuadraticWedge()
                    elem.GetPointIds().SetId(6, nid_map[node_ids[6]])
                    elem.GetPointIds().SetId(7, nid_map[node_ids[7]])
                    elem.GetPointIds().SetId(8, nid_map[node_ids[8]])
                    elem.GetPointIds().SetId(9, nid_map[node_ids[9]])
                    elem.GetPointIds().SetId(10, nid_map[node_ids[10]])
                    elem.GetPointIds().SetId(11, nid_map[node_ids[11]])
                    elem.GetPointIds().SetId(12, nid_map[node_ids[12]])
                    elem.GetPointIds().SetId(13, nid_map[node_ids[13]])
                    elem.GetPointIds().SetId(14, nid_map[node_ids[14]])
                else:
                    elem = vtkWedge()
                elem.GetPointIds().SetId(0, nid_map[node_ids[0]])
                elem.GetPointIds().SetId(1, nid_map[node_ids[1]])
                elem.GetPointIds().SetId(2, nid_map[node_ids[2]])
                elem.GetPointIds().SetId(3, nid_map[node_ids[3]])
                elem.GetPointIds().SetId(4, nid_map[node_ids[4]])
                elem.GetPointIds().SetId(5, nid_map[node_ids[5]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                min_thetai, max_thetai, dideal_thetai = get_min_max_theta(
                    _cpenta_faces, node_ids[:6], nid_map, xyz_cid0)
            elif isinstance(element, (CHEXA8, CIHEX1)):
                node_ids = element.node_ids
                pid = element.Pid()
                for nid in node_ids:
                    nid_to_pid_map[nid].append(pid)
                self.eid_to_nid_map[eid] = node_ids[:8]
                elem = vtkHexahedron()
                elem.GetPointIds().SetId(0, nid_map[node_ids[0]])
                elem.GetPointIds().SetId(1, nid_map[node_ids[1]])
                elem.GetPointIds().SetId(2, nid_map[node_ids[2]])
                elem.GetPointIds().SetId(3, nid_map[node_ids[3]])
                elem.GetPointIds().SetId(4, nid_map[node_ids[4]])
                elem.GetPointIds().SetId(5, nid_map[node_ids[5]])
                elem.GetPointIds().SetId(6, nid_map[node_ids[6]])
                elem.GetPointIds().SetId(7, nid_map[node_ids[7]])
                self.grid.InsertNextCell(12, elem.GetPointIds())
                min_thetai, max_thetai, dideal_thetai = get_min_max_theta(
                    _chexa_faces, node_ids[:8], nid_map, xyz_cid0)
            elif isinstance(element, CHEXA20):
                node_ids = element.node_ids
                pid = element.Pid()
                for nid in node_ids:
                    if nid is not None:
                        nid_to_pid_map[nid].append(pid)
                if None not in node_ids:
                    elem = vtkQuadraticHexahedron()
                    elem.GetPointIds().SetId(8, nid_map[node_ids[8]])
                    elem.GetPointIds().SetId(9, nid_map[node_ids[9]])
                    elem.GetPointIds().SetId(10, nid_map[node_ids[10]])
                    elem.GetPointIds().SetId(11, nid_map[node_ids[11]])
                    elem.GetPointIds().SetId(12, nid_map[node_ids[12]])
                    elem.GetPointIds().SetId(13, nid_map[node_ids[13]])
                    elem.GetPointIds().SetId(14, nid_map[node_ids[14]])
                    elem.GetPointIds().SetId(15, nid_map[node_ids[15]])
                    elem.GetPointIds().SetId(16, nid_map[node_ids[16]])
                    elem.GetPointIds().SetId(17, nid_map[node_ids[17]])
                    elem.GetPointIds().SetId(18, nid_map[node_ids[18]])
                    elem.GetPointIds().SetId(19, nid_map[node_ids[19]])
                else:
                    elem = vtkHexahedron()

                self.eid_to_nid_map[eid] = node_ids[:8]
                elem.GetPointIds().SetId(0, nid_map[node_ids[0]])
                elem.GetPointIds().SetId(1, nid_map[node_ids[1]])
                elem.GetPointIds().SetId(2, nid_map[node_ids[2]])
                elem.GetPointIds().SetId(3, nid_map[node_ids[3]])
                elem.GetPointIds().SetId(4, nid_map[node_ids[4]])
                elem.GetPointIds().SetId(5, nid_map[node_ids[5]])
                elem.GetPointIds().SetId(6, nid_map[node_ids[6]])
                elem.GetPointIds().SetId(7, nid_map[node_ids[7]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                min_thetai, max_thetai, dideal_thetai = get_min_max_theta(
                    _chexa_faces, node_ids[:8], nid_map, xyz_cid0)

            elif isinstance(element, CPYRAM5):
                node_ids = element.node_ids
                pid = element.Pid()
                for nid in node_ids:
                    nid_to_pid_map[nid].append(pid)
                self.eid_to_nid_map[eid] = node_ids[:5]
                elem = vtkPyramid()
                elem.GetPointIds().SetId(0, nid_map[node_ids[0]])
                elem.GetPointIds().SetId(1, nid_map[node_ids[1]])
                elem.GetPointIds().SetId(2, nid_map[node_ids[2]])
                elem.GetPointIds().SetId(3, nid_map[node_ids[3]])
                elem.GetPointIds().SetId(4, nid_map[node_ids[4]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                min_thetai, max_thetai, dideal_thetai = get_min_max_theta(
                    _cpyram_faces, node_ids[:5], nid_map, xyz_cid0)
            elif isinstance(element, CPYRAM13):
                node_ids = element.node_ids
                pid = element.Pid()
                #if None not in node_ids:
                    #print(' node_ids =', node_ids)
                    #elem = vtkQuadraticPyramid()
                    #elem.GetPointIds().SetId(5, nid_map[node_ids[5]])
                    #elem.GetPointIds().SetId(6, nid_map[node_ids[6]])
                    #elem.GetPointIds().SetId(7, nid_map[node_ids[7]])
                    #elem.GetPointIds().SetId(8, nid_map[node_ids[8]])
                    #elem.GetPointIds().SetId(9, nid_map[node_ids[9]])
                    #elem.GetPointIds().SetId(10, nid_map[node_ids[10]])
                    #elem.GetPointIds().SetId(11, nid_map[node_ids[11]])
                    #elem.GetPointIds().SetId(12, nid_map[node_ids[12]])
                #else:
                elem = vtkPyramid()
                #print('*node_ids =', node_ids[:5])

                self.eid_to_nid_map[eid] = node_ids[:5]

                elem.GetPointIds().SetId(0, nid_map[node_ids[0]])
                elem.GetPointIds().SetId(1, nid_map[node_ids[1]])
                elem.GetPointIds().SetId(2, nid_map[node_ids[2]])
                elem.GetPointIds().SetId(3, nid_map[node_ids[3]])
                elem.GetPointIds().SetId(4, nid_map[node_ids[4]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                min_thetai, max_thetai, dideal_thetai = get_min_max_theta(
                    _cpyram_faces, node_ids[:5], nid_map, xyz_cid0)

            elif (isinstance(element, (LineElement, SpringElement)) or
                  element.type in ['CBUSH', 'CBUSH1D', 'CFAST', 'CROD', 'CONROD',
                                   'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
                                   'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
                                   'CVISC', 'CGAP']):

                # TODO: verify
                # CBUSH, CBUSH1D, CFAST, CROD, CELAS1, CELAS3
                # CDAMP1, CDAMP2, CDAMP3, CDAMP4, CDAMP5, CVISC
                if hasattr(element, 'pid'):
                    pid = element.Pid()
                else:
                    # CONROD
                    # CELAS2, CELAS4?
                    pid = 0
                node_ids = element.node_ids
                for nid in node_ids:
                    if nid is not None:
                        nid_to_pid_map[nid].append(pid)

                if node_ids[0] is None and  node_ids[0] is None: # CELAS2
                    print('removing CELASx eid=%i -> no node %s' % (eid, node_ids[0]))
                    del self.eid_map[eid]
                    continue
                if None in node_ids:  # used to be 0...
                    if node_ids[0] is None:
                        slot = 1
                    elif node_ids[1] is None:
                        slot = 0
                    #print('node_ids=%s slot=%s' % (str(node_ids), slot))
                    self.eid_to_nid_map[eid] = node_ids[slot]
                    nid = node_ids[slot]
                    if nid not in nid_map:
                        # SPOINT
                        print('removing CELASx eid=%i -> SPOINT %i' % (eid, nid))
                        continue

                    #c = nid_map[nid]
                    elem = vtk.vtkVertex()
                    elem.GetPointIds().SetId(0, j)

                    elem = vtk.vtkSphere()
                    #if d == 0.:
                        #d = sphere_size
                    elem.SetRadius(sphere_size)
                else:
                    # 2 points
                    #d = norm(element.nodes[0].get_position() - element.nodes[1].get_position())
                    self.eid_to_nid_map[eid] = node_ids
                    elem = vtk.vtkLine()
                    try:
                        elem.GetPointIds().SetId(0, nid_map[node_ids[0]])
                        elem.GetPointIds().SetId(1, nid_map[node_ids[1]])
                    except KeyError:
                        print("node_ids =", node_ids)
                        print(str(element))
                        continue

                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            else:
                print('removing\n%s' % (elem))
                print('removing eid=%s; %s' % (eid, element.type))
                del self.eid_map[eid]
                self.log_info("skipping %s" % element.type)
                continue
            # what about MPCs, RBE2s (rigid elements)?
            #   are they plotted as elements?
            #   and thus do they need a property?

            if pid is None:
                # CONROD
                #print(element)
                #pids[i] = 0
                #pids_dict[eid] = 0
                pass
            else:
                pids[i] = pid
                pids_dict[eid] = pid
            if np.isnan(max_thetai):
                print('eid=%s theta=%s...setting to 360. deg' % (eid, max_thetai))
                print(str(element).rstrip())
                for node in element.nodes:
                    print(str(node).rstrip())
                max_thetai = 2 * np.pi
            min_interior_angle[i] = min_thetai
            max_interior_angle[i] = max_thetai
            dideal_theta[i] = dideal_thetai
            max_skew_angle[i] = max_skew
            max_warp_angle[i] = max_warp
            max_aspect_ratio[i] = aspect_ratio
            area[i] = areai
            area_ratio[i] = area_ratioi
            taper_ratio[i] = taper_ratioi
            i += 1
        assert len(self.eid_map) > 0, self.eid_map

        nelements = i
        self.nElements = nelements
        #print('nelements=%s pids=%s' % (nelements, list(pids)))
        pids = pids[:nelements]
        self.grid.SetPoints(points)

        self.grid.Modified()
        if hasattr(self.grid, 'Update'):
            self.grid.Update()
        #self.log_info("updated grid")

        cases = OrderedDict()
        #pids = array(pids, 'int32')
        #print('eid_map')
        #for key, value in sorted(iteritems(self.eid_map)):
            #print('  %s %s' % (key, value))

        #if 0:
            #if not len(pids) == len(self.eid_map):
                #msg = 'ERROR:  len(pids)=%s len(eid_map)=%s\n' % (len(pids), len(self.eid_map))
                #for eid, pid in sorted(iteritems(pids_dict)):
                    ##self.eid_map[eid] = i
                    ##pids_dict[eid] = pid
                    #if eid not in self.eid_map:
                        #msg += 'eid=%s %s' % (eid, str(model.elements[eid]))
                #raise RuntimeError(msg)
        del pids_dict


        self.iSubcaseNameMap = {1: ['Nastran', '']}
        #nelements = len(self.eid_map)
        icase = 0
        form = ['Geometry', None, []]
        form0 = form[2]

        #new_cases = True
        # set to True to enable node_ids as an result
        nids_set = True
        if nids_set:
            nids = np.zeros(self.nNodes, dtype='int32')
            cds = np.zeros(self.nNodes, dtype='int32')
            for (nid, nid2) in iteritems(self.nid_map):
                nids[nid2] = nid
                node = model.Node(nid)
                try:
                    cds[nid2] = node.Cd()
                except AttributeError:
                    # SPOINTs
                    msg = 'nid=%s does not have a Node Cd\n%s' % (nid, node)
                    #raise AttributeError(msg)
                    continue

            nid_res = GuiResult(0, header='NodeID', title='NodeID',
                                location='node', scalar=nids)
            cases[icase] = (nid_res, (0, 'NodeID'))
            form0.append(('NodeID', icase, []))
            icase += 1

            if len(np.unique(cds)) > 1:
                cd_res = GuiResult(0, header='NodeCd', title='NodeCd',
                                   location='node', scalar=cds)
                cases[icase] = (cd_res, (0, 'NodeCd'))
                form0.append(('NodeCd', icase, []))
                icase += 1
            self.node_ids = nids

        # set to True to enable elementIDs as a result
        eids_set = True
        if eids_set:
            eids = np.zeros(nelements, dtype='int32')
            for (eid, eid2) in iteritems(self.eid_map):
                eids[eid2] = eid

            #if new_cases:
            eid_res = GuiResult(0, header='ElementID', title='ElementID',
                                location='centroid', scalar=eids)
            cases[icase] = (eid_res, (0, 'ElementID'))
            #else:
                #cases[(0, icase, 'ElementID', 1, 'centroid', '%i', '')] = eids
            form0.append(('ElementID', icase, []))
            icase += 1
            self.element_ids = eids

        prop_types_with_mid = [
            'PSOLID', 'PSHEAR',
            'PROD', 'CROD', 'PTUBE', 'PBAR', 'PBARL', 'PBEAM', 'PBEAML',
        ]
        # subcase_id, resultType, vector_size, location, dataFormat
        if len(model.properties):

            pid_res = GuiResult(0, header='PropertyID', title='PropertyID',
                                location='centroid', scalar=pids)
            cases[icase] = (pid_res, (0, 'PropertyID'))
            form0.append(('PropertyID', icase, []))
            icase += 1

            upids = np.unique(pids)
            mids = np.zeros(nelements, dtype='int32')
            thickness = np.zeros(nelements, dtype='float32')
            mid_eids_skip = []
            for pid in upids:
                if pid == 0:
                    print('skipping pid=0')
                    continue
                prop = model.properties[pid]
                #try:
                if prop.type in prop_types_with_mid:
                    i = np.where(pids == pid)[0]
                    #print('pid=%s i=%s' % (pid, i))
                    #if isinstance(prop.mid, (int, int32)):
                        #mid = prop.mid
                    #else:
                        #try:
                    mid = prop.mid_ref.mid
                        #except AttributeError:
                            #print('pid=%s prop.type=%s' % (pid, prop.type))
                            #raise
                    mids[i] = mid
                elif prop.type == 'PSHELL':
                    # TODO: only considers mid1
                    i = np.where(pids == pid)[0]
                    mid = prop.Mid1()
                    t = prop.Thickness()
                    mids[i] = mid
                    thickness[i] = t
                elif prop.type == 'PCOMP':
                    # TODO: only considers iply=0
                    i = np.where(pids == pid)[0]
                    mid = prop.Mid(0)
                    t = prop.Thickness()
                    mids[i] = mid
                    thickness[i] = t
                elif prop.type in ['PELAS', 'PBUSH']:
                    i = np.where(pids == pid)[0]
                    mid_eids_skip.append(i)
                else:
                    print('material for pid=%s type=%s not considered' % (pid, prop.type))

            #print('mids =', mids)
            if len(mid_eids_skip):
                mid_eids_skip = np.hstack(mid_eids_skip)
                if mids.min() == 0:
                    i = np.where(mids == 0)[0]
                    diff_ids = np.setdiff1d(i, mid_eids_skip)
                    #eids_missing_material_id = eids[i]
                    not_skipped_eids_missing_material_id = eids[diff_ids]
                    if len(not_skipped_eids_missing_material_id):
                        print('eids=%s dont have materials' %
                              not_skipped_eids_missing_material_id)

            e11 = np.zeros(mids.shape, dtype='float32')
            e22 = np.zeros(mids.shape, dtype='float32')
            for umid in np.unique(mids):
                if umid == 0:
                    continue
                try:
                    mat = model.materials[umid]
                except KeyError:
                    print('mids = %s' % mids)
                    print('mids = %s' % model.materials.keys())
                    continue
                    #raise
                if mat.type == 'MAT1':
                    e11i = e22i = mat.e
                elif mat.type == 'MAT8':
                    e11i = mat.e11
                    e22i = mat.e22
                else:
                    print('skipping %s' % mat)
                    continue
                    #raise NotImplementedError(mat)
                #print('mid=%s e11=%e e22=%e' % (umid, e11i, e22i))
                i = np.where(umid == mids)[0]
                e11[i] = e11i
                e22[i] = e22i

            if thickness.max() > 0.0:
                t_res = GuiResult(0, header='Thickness', title='Thickness',
                                  location='centroid', scalar=thickness)
                cases[icase] = (t_res, (0, 'Thickness'))
                form0.append(('Thickness', icase, []))
                icase += 1

            mid_res = GuiResult(0, header='MaterialID', title='MaterialID',
                                location='centroid', scalar=mids)
            #e11_res = GuiResult(0, header='E_11', title='E_11',
                                #location='centroid', scalar=e11)
            cases[icase] = (mid_res, (0, 'MaterialID'))
            form0.append(('MaterialID', icase, []))

            is_orthotropic = not np.array_equal(e11, e22)
            if is_orthotropic:
                e11_res = GuiResult(0, header='E_11', title='E_11',
                                    location='centroid', scalar=e11, data_format='%.3e')
                e22_res = GuiResult(0, header='E_22', title='E_22',
                                    location='centroid', scalar=e22, data_format='%.3e')
                cases[icase + 1] = (e11_res, (0, 'E_11'))
                cases[icase + 2] = (e22_res, (0, 'E_22'))
                form0.append(('E_11', icase + 1, []))
                form0.append(('E_22', icase + 2, []))
                icase += 3
            else:
                # isotropic
                e11_res = GuiResult(0, header='E', title='E',
                                    location='centroid', scalar=e11, data_format='%.3e')
                cases[icase + 1] = (e11_res, (0, 'E'))
                form0.append(('E', icase + 1, []))
                icase += 2

        try:
            icase = self._build_optimization(model, pids, upids, nelements, cases, form0, icase)
        except:
            s = StringIO()
            traceback.print_exc(file=s)
            sout = s.getvalue()
            self.log_error(sout)
            print(sout)
            #traceback.print_exc(file=sys.stdout)
            #etype, value, tb = sys.exc_info
            #print(etype, value, tb)
            #raise RuntimeError('Optimization Parsing Error') from e
            #traceback.print_tb(e)
            #print(e)

            #mid_eids_skip = []
            #for pid in upids:

        if self.make_offset_normals_dim:
            #ielement = 0
            nelements = self.element_ids.shape[0]
            normals = np.zeros((nelements, 3), dtype='float32')
            offset = np.zeros(nelements, dtype='float32')
            xoffset = np.zeros(nelements, dtype='float32')
            yoffset = np.zeros(nelements, dtype='float32')
            zoffset = np.zeros(nelements, dtype='float32')
            element_dim = np.zeros(nelements, dtype='int32')
            nnodes_array = np.zeros(nelements, dtype='int32')
            for eid, element in sorted(iteritems(model.elements)):
                if isinstance(element, ShellElement):
                    element_dimi = 2
                    try:
                        normali = element.Normal()
                    except RuntimeError:
                        normali = np.ones(3) * 2.
                    #pid = element.pid
                    pid = element.pid
                    pid_type = pid.type
                    if pid_type == 'PSHELL':
                        z0 = element.pid.z1
                    elif pid_type == 'PCOMP':
                        z0 = element.pid.z0
                    else:
                        raise NotImplementedError(pid_type) # PSHEAR, PCOMPG

                    if z0 is None:
                        if element.type in ['CTRIA3', 'CTRIAR', 'CTRIAX']:
                            z0 = (element.T1 + element.T2 + element.T3) / 3.
                            nnodesi = 3
                        elif element.type in ['CTRIA6', 'CTRIAX6']:
                            z0 = (element.T1 + element.T2 + element.T3) / 3.
                            nnodesi = 6

                        elif element.type in ['CQUAD4', 'CQUADR']:
                            z0 = (element.T1 + element.T2 + element.T3 + element.T4) / 4.
                            nnodesi = 4
                        elif element.type == 'CQUAD8':
                            z0 = (element.T1 + element.T2 + element.T3 + element.T4) / 4.
                            nnodesi = 8
                        elif element.type == 'CQUAD':
                            z0 = (element.T1 + element.T2 + element.T3 + element.T4) / 4.
                            nnodesi = 9
                        else:
                            raise NotImplementedError(element.type)
                    else:
                        if element.type in ['CTRIA3', 'CTRIAR', 'CTRIAX']:
                            nnodesi = 3
                        elif element.type in ['CTRIA6', 'CTRIAX6']:
                            nnodesi = 6

                        elif element.type in ['CQUAD4', 'CQUADR']:
                            nnodesi = 4
                        elif element.type == 'CQUAD8':
                            nnodesi = 8
                        elif element.type == 'CQUAD':
                            nnodesi = 9
                        else:
                            raise NotImplementedError(element.type)
                    zi = element.zOffset + z0

                    ie = self.eid_map[eid]
                    normals[ie, :] = normali
                    offset[ie] = zi
                    xoffset[ie] = zi * normali[0]
                    yoffset[ie] = zi * normali[1]
                    zoffset[ie] = zi * normali[2]

                elif element.type == 'CTETRA':
                    ie = self.eid_map[eid]
                    element_dimi = 3
                    nnodesi = 4
                elif element.type == 'CPENTA':
                    ie = self.eid_map[eid]
                    element_dimi = 3
                    nnodesi = 6
                elif element.type == 'CPYRAM':
                    ie = self.eid_map[eid]
                    element_dimi = 3
                    nnodesi = 5
                elif element.type in ['CHEXA', 'CIHEX1']:
                    ie = self.eid_map[eid]
                    element_dimi = 3
                    nnodesi = 8

                elif element.type in ['CROD', 'CONROD', 'CBEND', 'CBAR', 'CBEAM', 'CGAP']:
                    ie = self.eid_map[eid]
                    element_dimi = 1
                    nnodesi = 2
                elif element.type in ['CBUSH', 'CBUSH1D', 'CFAST', 'CVISC',
                                      'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
                                      'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5']:
                    ie = self.eid_map[eid]
                    element_dimi = 0
                    nnodesi = 2
                else:
                    ie = self.eid_map[eid]
                    element_dimi = -1
                    nnodesi = -1
                    print('element.type=%s doesnt have a dimension' % element.type)

                element_dim[ie] = element_dimi
                nnodes_array[ie] = nnodesi
                #ielement += 1

            # if not a flat plate
            #if min(nxs) == max(nxs) and min(nxs) != 0.0:
            eid_dim_res = GuiResult(0, header='ElementDim', title='ElementDim',
                                    location='centroid', scalar=element_dim)
            cases[icase] = (eid_dim_res, (0, 'ElementDim'))

            is_shell = np.abs(normals).max() > 0.
            is_solid = np.abs(max_interior_angle).max() > 0.
            #print('is_shell=%s is_solid=%s' % (is_shell, is_solid))
            if is_shell:
                nx_res = GuiResult(0, header='NormalX', title='NormalX',
                                   location='centroid', scalar=normals[:, 0], data_format='%.2f')
                ny_res = GuiResult(0, header='NormalY', title='NormalY',
                                   location='centroid', scalar=normals[:, 1], data_format='%.2f')
                nz_res = GuiResult(0, header='NormalZ', title='NormalZ',
                                   location='centroid', scalar=normals[:, 2], data_format='%.2f')

                # this is just for testing nan colors that doesn't work
                #max_interior_angle[:1000] = np.nan
                area_res = GuiResult(0, header='Area', title='Area',
                                     location='centroid', scalar=area)
                min_theta_res = GuiResult(0, header='Min Interior Angle', title='Min Interior Angle',
                                          location='centroid', scalar=np.degrees(min_interior_angle))
                max_theta_res = GuiResult(0, header='Max Interior Angle', title='Max Interior Angle',
                                          location='centroid', scalar=np.degrees(max_interior_angle))
                dideal_theta_res = GuiResult(0, header='Delta Ideal Angle', title='Delta Ideal Angle',
                                             location='centroid', scalar=np.degrees(dideal_theta))
                eid_dim_res = GuiResult(0, header='ElementDim', title='ElementDim',
                                        location='centroid', scalar=element_dim)

                skew = np.degrees(max_skew_angle)
                skew_res = GuiResult(0, header='Max Skew Angle', title='MaxSkewAngle',
                                     location='centroid', scalar=skew)
                aspect_res = GuiResult(0, header='Aspect Ratio', title='AspectRatio',
                                       location='centroid', scalar=max_aspect_ratio)

                form_checks = []
                form0.append(('Element Checks', None, form_checks))
                form_checks.append(('ElementDim', icase, []))

                if self.make_nnodes_result:
                    nnodes_res = GuiResult(0, header='NNodes/Elem', title='NNodes/Elem',
                                           location='centroid', scalar=nnodes_array)
                    form_checks.append(('NNodes', icase + 1, []))
                    cases[icase + 1] = (nnodes_res, (0, 'NNodes'))
                    icase += 1

                cases[icase + 1] = (nx_res, (0, 'NormalX'))
                cases[icase + 2] = (ny_res, (0, 'NormalY'))
                cases[icase + 3] = (nz_res, (0, 'NormalZ'))
                cases[icase + 4] = (area_res, (0, 'Area'))
                cases[icase + 5] = (min_theta_res, (0, 'Min Interior Angle'))
                cases[icase + 6] = (max_theta_res, (0, 'Max Interior Angle'))
                cases[icase + 7] = (dideal_theta_res, (0, 'Delta Ideal Angle'))
                cases[icase + 8] = (skew_res, (0, 'Max Skew Angle'))
                cases[icase + 9] = (aspect_res, (0, 'Aspect Ratio'))

                form_checks.append(('NormalX', icase + 1, []))
                form_checks.append(('NormalY', icase + 2, []))
                form_checks.append(('NormalZ', icase + 3, []))
                form_checks.append(('Area', icase + 4, []))
                form_checks.append(('Min Interior Angle', icase + 5, []))
                form_checks.append(('Max Interior Angle', icase + 6, []))
                form_checks.append(('Delta Ideal Angle', icase + 7, []))
                form_checks.append(('Max Skew Angle', icase + 8, []))
                form_checks.append(('Aspect Ratio', icase + 9, []))
                icase += 10

                if area_ratio.max() > 1.:
                    arearatio_res = GuiResult(0, header='Area Ratio', title='Area Ratio',
                                              location='centroid', scalar=area_ratio)
                    cases[icase] = (arearatio_res, (0, 'Area Ratio'))
                    form_checks.append(('Area Ratio', icase, []))
                    icase += 1

                if taper_ratio.max() > 1.:
                    taperratio_res = GuiResult(0, header='Taper Ratio', title='Taper Ratio',
                                               location='centroid', scalar=taper_ratio)
                    cases[icase] = (taperratio_res, (0, 'Taper Ratio'))
                    form_checks.append(('Taper Ratio', icase, []))
                    icase += 1

                if max_warp_angle.max() > 0.0:
                    warp_res = GuiResult(0, header='Max Warp Angle', title='MaxWarpAngle',
                                         location='centroid', scalar=np.degrees(max_warp_angle))
                    cases[icase + 4] = (warp_res, (0, 'Max Warp Angle'))
                    form_checks.append(('Max Warp Angle', icase, []))
                    icase += 1

                #if (np.abs(xoffset).max() > 0.0 or np.abs(yoffset).max() > 0.0 or
                    #np.abs(zoffset).max() > 0.0):
                # offsets
                offset_res = GuiResult(0, header='Offset', title='Offset',
                                       location='centroid', scalar=offset, data_format='%g')
                offset_x_res = GuiResult(0, header='OffsetX', title='OffsetX',
                                         location='centroid', scalar=xoffset, data_format='%g')
                offset_y_res = GuiResult(0, header='OffsetY', title='OffsetY',
                                         location='centroid', scalar=yoffset, data_format='%g')
                offset_z_res = GuiResult(0, header='OffsetZ', title='OffsetZ',
                                         location='centroid', scalar=zoffset, data_format='%g')

                cases[icase] = (offset_res, (0, 'Offset'))
                cases[icase + 1] = (offset_x_res, (0, 'OffsetX'))
                cases[icase + 2] = (offset_y_res, (0, 'OffsetY'))
                cases[icase + 3] = (offset_z_res, (0, 'OffsetZ'))

                form_checks.append(('Offset', icase, []))
                form_checks.append(('OffsetX', icase + 1, []))
                form_checks.append(('OffsetY', icase + 2, []))
                form_checks.append(('OffsetZ', icase + 3, []))
                icase += 4

                if self.make_xyz:
                    x_res = GuiResult(0, header='X', title='X',
                                      location='node', scalar=xyz_cid0[:, 0], data_format='%g')
                    y_res = GuiResult(0, header='Y', title='Y',
                                      location='node', scalar=xyz_cid0[:, 1], data_format='%g')
                    z_res = GuiResult(0, header='Z', title='Z',
                                      location='node', scalar=xyz_cid0[:, 2], data_format='%g')
                    cases[icase] = (x_res, (0, 'X'))
                    cases[icase + 1] = (y_res, (0, 'Y'))
                    cases[icase + 2] = (z_res, (0, 'Z'))
                    form_checks.append(('X', icase + 0, []))
                    form_checks.append(('Y', icase + 1, []))
                    form_checks.append(('Z', icase + 2, []))
                    icase += 3

            elif is_solid:
                # only solid elements
                form_checks = []
                form0.append(('Element Checks', None, form_checks))
                min_theta_res = GuiResult(0, header='Min Interior Angle', title='Min Interior Angle',
                                          location='centroid', scalar=np.degrees(min_interior_angle))
                max_theta_res = GuiResult(0, header='Max Interior Angle', title='Max Interior Angle',
                                          location='centroid', scalar=np.degrees(max_interior_angle))
                #skew = 90. - np.degrees(max_skew_angle)
                #skew_res = GuiResult(0, header='Max Skew Angle', title='MaxSkewAngle',
                                     #location='centroid', scalar=skew)
                form_checks.append(('ElementDim', icase, []))
                form_checks.append(('Min Interior Angle', icase + 1, []))
                form_checks.append(('Max Interior Angle', icase + 2, []))
                #form_checks.append(('Max Skew Angle', icase + 2, []))
                cases[icase + 1] = (min_theta_res, (0, 'Min Interior Angle'))
                cases[icase + 2] = (max_theta_res, (0, 'Max Interior Angle'))
                #cases[icase + 3] = (skew_res, (0, 'Max Interior Angle'))
                icase += 3

            else:
                form0.append(('ElementDim', icase, []))
                icase += 1

            if np.abs(material_coord).max() > 0:
                material_coord_res = GuiResult(0, header='MaterialCoord', title='MaterialCoord',
                                               location='centroid',
                                               scalar=material_coord, data_format='%i')
                cases[icase] = (material_coord_res, (0, 'MaterialCoord'))
                form0.append(('MaterialCoord', icase, []))
                icase += 1

            self.normals = normals
        return nid_to_pid_map, icase, cases, form

    def _build_optimization(self, model, pids, upids, nelements, cases, form0, icase):
        if len(model.properties) and len(model.dvprels):
            # len(model.dvprels) + len(model.dvcrels) + len(model.dvmrels) + len(model.desvars)
            #dvmrel_init = np.zeros(nelements, dtype='int32')
            #dvgrel_init = np.zeros(nelements, dtype='int32')
            dvprel_t_init = np.zeros(nelements, dtype='float32')
            dvprel_t_min = np.zeros(nelements, dtype='float32')
            dvprel_t_max = np.zeros(nelements, dtype='float32')
            design_region = np.zeros(nelements, dtype='int32')

            for key, dvprel in iteritems(model.dvprels):
                if dvprel.type == 'DVPREL1':
                    prop_type = dvprel.Type
                    desvars = dvprel.dvids
                    coeffs = dvprel.coeffs
                    pid = dvprel.pid.pid
                    var_to_change = dvprel.pname_fid
                    assert len(desvars) == 1, len(desvars)

                    if prop_type == 'PSHELL':
                        i = np.where(pids == pid)
                        design_region[i] = dvprel.oid
                        assert len(i) > 0, i
                        if var_to_change == 'T':
                            value = 0.
                            lower_bound = 0.
                            upper_bound = 0.
                            for desvar, coeff in zip(desvars, coeffs):
                                xiniti = desvar.xinit
                                if desvar.xlb != -1e20:
                                    xiniti = max(xiniti, desvar.xlb)
                                    lower_bound = desvar.xlb
                                if desvar.xub != 1e20:
                                    xiniti = min(xiniti, desvar.xub)
                                    upper_bound = desvar.xub
                                if desvar.delx is not None and desvar.delx != 1e20:
                                    desvar.delx
                                if desvar.ddval is not None:
                                    msg = 'DESVAR id=%s DDVAL is not None\n%s' % str(desvar)
                                assert desvar.ddval is None, desvar
                                xinit = coeff * xiniti
                            dvprel_t_init[i] = xinit
                            dvprel_t_min[i] = lower_bound
                            dvprel_t_max[i] = upper_bound

                        else:
                            raise NotImplementedError(dvprel)
                    else:
                        raise NotImplementedError(dvprel)
                else:
                    raise NotImplementedError(dvprel)
                if dvprel.p_max != 1e20:
                    dvprel.p_max
                if dvprel.p_min is not None:
                    dvprel.p_min

            region_res = GuiResult(0, header='DV Region', title='DV Region',
                                   location='centroid', scalar=design_region)
            t_init_res = GuiResult(0, header='DVPREL Init - t', title='DVPREL Init - t',
                                   location='centroid', scalar=dvprel_t_init)
            t_min_res = GuiResult(0, header='DVPREL Min - t', title='DVPREL Min - t',
                                  location='centroid', scalar=dvprel_t_min)
            t_max_res = GuiResult(0, header='DVPREL Max - t', title='DVPREL Max - t',
                                  location='centroid', scalar=dvprel_t_max)
            cases[icase] = (region_res, (0, 'DV Region'))
            cases[icase + 1] = (t_init_res, (0, 'DVPREL Init - t'))
            cases[icase + 2] = (t_min_res, (0, 'DVPREL Min - t'))
            cases[icase + 3] = (t_max_res, (0, 'DVPREL Max - t'))
            form0.append(('DV Region', icase, []))
            form0.append(('DVPREL Init - t', icase + 1, []))
            form0.append(('DVPREL Min - t', icase + 2, []))
            form0.append(('DVPREL Max - t', icase + 3, []))
            icase += 4
        return icase

    def _plot_pressures(self, model, cases, form0, icase, subcase_id, subcase):
        """
        pressure act normal to the face (as opposed to anti-normal)
        """
        eids = self.element_ids

        try:
            load_case_id = subcase.get_parameter('LOAD')[0]
        except KeyError:
            print('no LOAD for isubcase=%s' % subcase_id)
            return icase

        try:
            load_case = model.loads[load_case_id]
        except KeyError:
            self.log.warning('LOAD=%s not found' % load_case_id)
            return icase

        # account for scale factors
        loads2 = []
        scale_factors2 = []
        for load in load_case:
            if isinstance(load, LOAD):
                scale_factors, loads = load.get_reduced_loads()
                scale_factors2 += scale_factors
                loads2 += loads
            else:
                scale_factors2.append(1.)
                loads2.append(load)

        pressures = np.zeros(len(model.elements), dtype='float32')

        iload = 0
        nloads = len(loads2)
        show_nloads = nloads > 5000
        # loop thru scaled loads and plot the pressure
        for load, scale in zip(loads2, scale_factors2):
            if show_nloads and iload % 5000 == 0:
                self.log_debug('  NastranIOv iload=%s/%s' % (iload, nloads))
            if load.type == 'PLOAD4':
                #print(load.object_attributes())
                for elem in load.eids:
                    #elem = self.model.elements[eid]
                    if elem.type in ['CTRIA3', 'CTRIA6', 'CTRIA', 'CTRIAR',
                                     'CQUAD4', 'CQUAD8', 'CQUAD', 'CQUADR', 'CSHEAR']:
                        pressure = load.pressures[0] * scale

                        # single element per PLOAD
                        #eid = elem.eid
                        #pressures[eids.index(eid)] = pressure

                        # multiple elements
                        #for elem in load.eids:
                        ie = np.searchsorted(eids, elem.eid)
                        #pressures[ie] += p  # correct; we can't assume model orientation
                        nz = self.normals[ie, 2]  # considers normal of shell
                        pressures[ie] += pressure * nz

                    #elif elem.type in ['CTETRA', 'CHEXA', 'CPENTA']:
                        #A, centroid, normal = elem.getFaceAreaCentroidNormal(
                            #load.g34.nid, load.g1.nid)
                        #r = centroid - p
            iload += 1
        # if there is no applied pressure, don't make a plot
        if np.abs(pressures).max():
            case_name = 'Pressure'
            # print('iload=%s' % iload)
            # print(case_name)
            # subcase_id, resultType, vector_size, location, dataFormat
            cases[(0, icase, case_name, 1, 'centroid', '%.1f', '')] = pressures
            form0.append((case_name, icase, []))
            icase += 1
        return icase

    def _plot_applied_loads(self, model, cases, form0, icase, subcase_id, subcase):
        found_load = False
        found_temperature = False

        form = []
        load_keys = (
            'LOAD', 'TEMPERATURE(MATERIAL)', 'TEMPERATURE(INITIAL)',
            'TEMPERATURE(LOAD)', 'TEMPERATURE(BOTH)')
        temperature_keys = (
            'TEMPERATURE(MATERIAL)', 'TEMPERATURE(INITIAL)',
            'TEMPERATURE(LOAD)', 'TEMPERATURE(BOTH)')

        for key in load_keys:
            try:
                load_case_id = subcase.get_parameter(key)[0]
            except KeyError:
                # print('no %s for isubcase=%s' % (key, subcase_id))
                continue
            try:
                load_case = model.loads[load_case_id]
            except KeyError:
                self.log.warning('LOAD=%s not found' % load_case_id)
                continue

            if key == 'LOAD':
                p0 = np.array([0., 0., 0.], dtype='float32')
                centroidal_pressures, forces, spcd = self._get_forces_moments_array(
                    model, p0, load_case_id, include_grav=False)
                found_load = True
            elif key in temperature_keys:
                temperatures = self._get_temperatures_array(model, load_case_id)
                found_temperature = True
                temperature_key = key
            else:
                raise NotImplementedError(key)

        if found_load:
            if np.abs(centroidal_pressures).max():
                # print('iload=%s' % iload)
                # print(case_name)
                # subcase_id, resultType, vector_size, location, dataFormat
                cases[(0, icase, 'Pressure', 1, 'centroid', '%.1f', '')] = centroidal_pressures
                form0.append(('Pressure', icase, []))
                icase += 1

            if np.abs(forces.max() - forces.min()) > 0.0:
                # if forces[:, 0].min() != forces[:, 0].max():
                cases[(subcase_id, icase, 'LoadX', 1, 'node', '%.1f', '')] = forces[:, 0]
                # if forces[:, 1].min() != forces[:, 1].max():
                cases[(subcase_id, icase + 1, 'LoadY', 1, 'node', '%.1f', '')] = forces[:, 1]
                # if forces[:, 2].min() != forces[:, 2].max():
                cases[(subcase_id, icase + 2, 'LoadZ', 1, 'node', '%.1f', '')] = forces[:, 2]

                form0.append(('Total Load FX', icase, []))
                form0.append(('Total Load FY', icase + 1, []))
                form0.append(('Total Load FZ', icase + 2, []))
                icase += 3

            if np.abs(spcd.max() - spcd.min()) > 0.0:
                cases[(subcase_id, icase, 'SPCDx', 1, 'node', '%.3g')] = spcd[:, 0]
                form0.append(('SPCDx', icase, []))
                icase += 1

                cases[(subcase_id, icase, 'SPCDy', 1, 'node', '%.3g')] = spcd[:, 1]
                form0.append(('SPCDy', icase, []))
                icase += 1

                #cases[(subcase_id, icase, name + 'Z', 1, 'node', '%g', header)] = t3
                cases[(subcase_id, icase, 'SPCDz', 1, 'node', '%.3g')] = spcd[:, 2]
                form0.append(('SPCDz', icase, []))
                icase += 1

                t123 = spcd[:, :3]
                tnorm = norm(t123, axis=1)
                assert len(tnorm) == len(spcd[:, 2]), len(spcd[:, 2])
                cases[(subcase_id, icase, 'SPCD XYZ', 1, 'node', '%.3g')] = tnorm
                form0.append(('SPCD XYZ', icase, []))
                icase += 1
        if found_temperature:
            cases[(subcase_id, icase, temperature_key, 1, 'node', '%.3g')] = temperatures
            form.append((temperature_key, icase, []))
            icase += 1
        return icase

    def _get_loads_and_scale_factors(self, load_case):
        # account for scale factors
        loads2 = []
        scale_factors2 = []
        for load in load_case:
            if isinstance(load, LOAD):
                scale_factors, loads = load.get_reduced_loads()
                scale_factors2 += scale_factors
                loads2 += loads
            else:
                scale_factors2.append(1.)
                loads2.append(load)
        return loads2, scale_factors2

    def _get_temperatures_array(self, model, load_case_id):
        """builds the temperature array based on thermal cards"""
        nids = sorted(model.nodes.keys())

        load_case = model.loads[load_case_id]
        loads2, scale_factors2 = self._get_loads_and_scale_factors(load_case)
        tempd = model.tempds[load_case_id].temperature if load_case_id in model.tempds else 0.
        temperatures = np.ones(len(model.nodes), dtype='float32') * tempd
        for load, scale in zip(loads2, scale_factors2):
            if load.type == 'TEMP':
                temps_dict = load.temperatures
                for nid, val in iteritems(temps_dict):
                    nidi = nids.index(nid)
                    temperatures[nidi] = val
            else:
                print(load.type)
        return temperatures

    def _get_forces_moments_array(self, model, p0, load_case_id, include_grav=False):
        nids = sorted(model.nodes.keys())
        nnodes = len(nids)
        nid_map = self.nid_map
        eid_map = self.eid_map

        load_case = model.loads[load_case_id]
        loads2, scale_factors2 = self._get_loads_and_scale_factors(load_case)

        eids = sorted(model.elements.keys())
        centroidal_pressures = np.zeros(len(model.elements), dtype='float32')
        nodal_pressures = np.zeros(len(self.node_ids), dtype='float32')

        forces = np.zeros((nnodes, 3), dtype='float32')
        spcd = np.zeros((nnodes, 3), dtype='float32')
        # loop thru scaled loads and plot the pressure
        cards_ignored = {}

        assert self.normals is not None
        for load, scale in zip(loads2, scale_factors2):
            if load.type == 'FORCE':
                scale2 = load.mag * scale  # does this need a magnitude?
                nid = load.node
                if nid in self.dependents_nodes:
                    print('    nid=%s is a dependent node and has an FORCE applied\n%s' % (
                        nid, str(load)))
                forces[nid_map[nid]] += load.xyz * scale2

            elif load.type == 'PLOAD2':
                pressure = load.pressures[0] * scale  # there are 4 pressures, but we assume p0
                for eid in load.eids:
                    elem = self.elements[eid]
                    if elem.type in ['CTRIA3',
                                     'CQUAD4', 'CSHEAR']:
                        node_ids = elem.node_ids
                        nnodes = len(node_ids)
                        normal = elem.Normal()
                        area = elem.Area()
                        forcei = pressure * normal * area / nnodes
                        # r = elem.Centroid() - p0
                        # m = cross(r, f)
                        for nid in node_ids:
                            if nid in self.dependents_nodes:
                                print('    nid=%s is a dependent node and has an PLOAD2 applied\n'
                                      '%s' % (nid, str(load)))
                            forces[nid_map[nid]] += forcei
                        forces += forcei
                        # F += f
                        # M += m
                    else:
                        self.log.debug('    case=%s etype=%r loadtype=%r not supported' % (
                            load_case_id, elem.type, load.type))

            elif load.type == 'PLOAD4':
                # continue  ## TODO: should be removed
                # elem = load.eid
                # area = elem.get_area()
                if 0:
                    if elem.type in ['CTRIA3', 'CTRIA6', 'CTRIA', 'CTRIAR',]:
                        eid = elem.eid
                        node_ids = elem.node_ids
                        k = load.pressures[0] * scale / 3.
                        # TODO: doesn't consider load.eids for distributed pressures???
                        for nid in node_ids[3:]:
                            centroidal_pressures[nid_map[nid]] += k
                    elif elem.type in ['CQUAD4', 'CQUAD8', 'CQUAD', 'CQUADR', 'CSHEAR']:
                        eid = elem.eid
                        node_ids = elem.node_ids
                        k = load.pressures[0] * scale / 4.
                        # TODO: doesn't consider load.eids for distributed pressures???
                        for nid in node_ids[4:]:
                            if nid in self.dependents_nodes:
                                print('    nid=%s is a dependent node and has an PLOAD4 applied\n'
                                      '%s' % (nid, str(load)))
                            centroidal_pressures[nid_map[nid]] += k
                    else:
                        print('    PLOAD4 is unhandled\n%s' % str(load))

                else:
                    # single element per PLOAD
                    #eid = elem.eid
                    #pressures[eids.index(eid)] = p

                    pressure = load.pressures[0] * scale

                    # multiple elements
                    for elem in load.eids:
                        ie = eid_map[elem.eid]
                        nz = self.normals[ie, :]
                        # pressures[eids.index(elem.eid)] += p
                        if elem.type in ['CTRIA3', 'CTRIA6', 'CTRIA', 'CTRIAR',
                                         # TODO: this was split in bdf_methods...
                                         'CQUAD4', 'CQUAD8', 'CQUAD', 'CQUADR', 'CSHEAR']:
                            area = elem.get_area()
                            elem_node_ids = elem.node_ids
                            elem_nnodes = len(elem_node_ids)
                            forcei = pressure * area / elem_nnodes
                            for nid in elem_node_ids:
                                if nid in self.dependents_nodes:
                                    print('    nid=%s is a dependent node and has a'
                                          ' PLOAD4 applied\n%s' % (nid, str(load)))
                                #forces[nids.index(nid)] += F
                                i = nid_map[nid]
                                try:
                                    forces[i, :] += forcei * nz
                                except IndexError:
                                    print('i = %s' % i)
                                    print('normals.shape = %s' %  str(self.normals.shape))
                                    print('forces.shape = %s' % str(forces.shape))
                                    print('nz = ', self.normals[i, :])
                                    print('forces[i, :] = ', forces[i, :])
                                    raise
                        else:
                            elem_node_ids = elem.node_ids
                            if elem.type == 'CTETRA':
                                #face1 = elem.get_face(load.g1.nid, load.g34.nid)
                                face, area, centroid, normal = elem.getFaceAreaCentroidNormal(
                                    load.g1.nid, load.g34.nid)
                                #assert face == face1
                                nface = 3
                            elif elem.type == 'CHEXA':
                                #face1 = elem.get_face(load.g34.nid, load.g1.nid)
                                face, area, centroid, normal = elem.getFaceAreaCentroidNormal(
                                    load.g34.nid, load.g1.nid)
                                #assert face == face1
                                nface = 4
                            elif elem.type == 'CPENTA':
                                g1 = load.g1.nid
                                if load.g34 is None:
                                    #face1 = elem.get_face(g1)
                                    face, area, centroid, normal = elem.getFaceAreaCentroidNormal(g1)
                                    nface = 3
                                else:
                                    #face1 = elem.get_face(g1, load.g34.nid)
                                    face, area, centroid, normal = elem.getFaceAreaCentroidNormal(
                                        g1, load.g34.nid)
                                    nface = 4
                                #assert face == face1
                            else:
                                msg = ('case=%s eid=%s etype=%r loadtype=%r not supported'
                                       % (load_case_id, eid, elem.type, load.type))
                                self.log.debug(msg)
                                continue
                            pressures = load.pressures[:nface]
                            assert len(pressures) == nface
                            if min(pressures) != max(pressures):
                                pressure = np.mean(pressures)
                                #msg = '%s%s\npressure.min=%s != pressure.max=%s using average of %%s; load=%s eid=%%s'  % (
                                    #str(load), str(elem), min(pressures), max(pressures), load.sid)
                                #print(msg % (pressure, eid))
                            else:
                                pressure = pressures[0]
                            #centroidal_pressures
                            f = pressure * area * normal * scale
                            for inid in face:
                                inidi = nid_map[elem_node_ids[inid]]
                                nodal_pressures[inid] += pressure * scale / nface
                                forces[inidi, :] += f / nface
                            centroidal_pressures[ie] += pressure

                            #r = centroid - p
                            #load.cid.transformToGlobal()
                            #m = cross(r, f)
                            #M += m


                #elif elem.type in ['CTETRA', 'CHEXA', 'CPENTA']:
            elif load.type == 'SPCD':
                #self.gids = [integer(card, 2, 'G1'),]
                #self.constraints = [components_or_blank(card, 3, 'C1', 0)]
                #self.enforced = [double_or_blank(card, 4, 'D1', 0.0)]
                for nid, c1, d1 in zip(load.node_ids, load.constraints, load.enforced):
                    if nid in self.dependents_nodes:
                        self.log_warning('    nid=%s is a dependent node and has an'
                                         ' SPCD applied\n%s' % (nid, str(load)))
                    c1 = int(c1)
                    assert c1 in [1, 2, 3, 4, 5, 6], c1
                    if c1 < 4:
                        spcd[nid_map[nid], c1 - 1] = d1
            else:
                if load.type not in cards_ignored:
                    cards_ignored[load.type] = True
                    self.log_info('  NastranIOv _get_forces_moments_array - unsupported '
                                  'load.type = %s' % load.type)

        return centroidal_pressures, forces, spcd

    def load_nastran_results(self, op2_filename, dirname):
        """
        Loads the Nastran results into the GUI
        """
        #gridResult.SetNumberOfComponents(self.nElements)
        self. turn_text_on()
        self.scalarBar.VisibilityOn()
        self.scalarBar.Modified()
        #self.show_caero_mesh()

        print("trying to read...%s" % op2_filename)
        ext = os.path.splitext(op2_filename)[1].lower()

        if ext == '.op2':
            model = OP2(log=self.log, debug=True)

            if 0:
                model._results.saved = set([])
                all_results = model.get_all_results()
                desired_results = [
                    # nodal
                    # ---------
                    'displacements', 'velocities', 'accelerations', 'temperatures',
                    'constraint_forces', 'spc_forces', 'mpc_forces', 'eigenvectors',

                    #'gridPointForces',
                    #'stress',

                    # untested
                    'load_vectors',
                    'applied_loads',
                    'force_vectors',

                    # ---------
                    # centroidal
                    'stress',
                    'chexa_stress', 'cpenta_stress', 'ctetra_stress',

                    'ctria3_stress', 'ctria3_stress',
                    'cquad8_stress''cquad4_stress',

                    'ctria3_composite_stress', 'ctria3_composite_stress',
                    'cquad8_composite_stress''cquad4_composite_stress',

                    'cbar_stress', 'cbeam_stress',
                    'crod_stress', 'conrod_stress', 'ctube_stress',
                    'celas1_stress', 'celas2_stress', 'celas3_stress', 'celas4_stress',
                    #=================================================
                    'strain',
                    'chexa_strain', 'cpenta_strain', 'ctetra_strein',

                    'ctria3_strain', 'ctria3_strain',
                    'cquad8_strain', 'cquad4_strain',

                    'ctria3_composite_strain', 'ctria3_composite_strain',
                    'cquad8_composite_strain', 'cquad4_composite_strain',

                    'cbar_strain', 'cbeam_strain',
                    'crod_strain', 'conrod_strain', 'ctube_strain',
                    'celas1_strain', 'celas2_strain', 'celas3_strain', 'celas4_strain',
                ]
                for result in desired_results:
                    if result in all_results:
                        model._results.saved.add(result)
            model.read_op2(op2_filename, combine=False)

            if not self.is_testing:
                self.log.info(model.get_op2_stats())
            # print(model.get_op2_stats())

        elif ext == '.pch':
            raise NotImplementedError('*.pch is not implemented; filename=%r' % op2_filename)
        #elif ext == '.f06':
            #model = F06(log=self.log, debug=True)
            #model.set_vectorization(True)
            #model.read_f06(op2_filename)
        else:
            #print("error...")
            msg = 'extension=%r is not supported; filename=%r' % (ext, op2_filename)
            raise NotImplementedError(msg)

        if self.save_data:
            self.model_results = model

        #print(model.print_results())
        #self.iSubcaseNameMap[self.isubcase] = [Subtitle, Label]

        # tansform displacements into global coordinates
        try:
            i_transform = self.i_transform
            transforms = self.transforms
        except AttributeError:
            self.log.error('Skipping displacment transformation')
        else:
            model.transform_displacements_to_global(i_transform, transforms)

        #if 0:
            #cases = OrderedDict()
            #self.iSubcaseNameMap = {}
            #form = []
            #icase = 0
        #else:
        cases = self.result_cases
        form = self.get_form()
        icase = len(cases)
        # form = self.res_widget.get_form()

        subcase_ids = model.iSubcaseNameMap.keys()
        #self.iSubcaseNameMap = model.iSubcaseNameMap
        # self.iSubcaseNameMap = model.subcase_key
        #print(self.iSubcaseNameMap)
        for isubcase, values in iteritems(model.iSubcaseNameMap):
            if not isinstance(isubcase, integer_types):
                print('isubcase type =', type(isubcase))
                continue
            if isinstance(values, str):
                # eigenvalue???
                label = values
                print('label??? =', label)
            elif isinstance(values, list):
                print(values)
                subtitle, analysis_code, label = values
            else:
                print(values)
                print(type(values))
                raise RuntimeError(values)

            self.iSubcaseNameMap[isubcase] = [subtitle, label]
            del subtitle, label
        # self.iSubcaseNameMap = {subcase_id : label for
                                # in iteritems(model.iSubcaseNameMap)}

        form = self._fill_op2_output(op2_filename, cases, model, form, icase)
        self._finish_results_io2(form, cases)

        #name = 'spike'
        #eids = np.arange(10, 40)
        #self.create_group_with_name(name, eids)
        #self.post_group_by_name(name)


    def _fill_op2_output(self, op2_filename, cases, model, form, icase):
        """
        SOL 101 (Static)
        ----------------
        Subcase 1
         - DisplacementXYZ
         - SPCForceX
         - ...
         - Stress
           - oxx
         - Strain

        SOL 103 (modal)
        ---------------
        Subcase 1
         - mode 1; eigr=123.4
          - EigenvectorXYZ
          - Stress
        - mode 2: eigr=156.3
          - EigenvectorXYZ
          - Stress

        SOL 109 (Freq)
        --------------
        Subcase 1
         - freq=123.4
          - DisplacementXYZ
          - Stress

        SOL 105 (Buckling)
        ------------------
        Subcase 1
         - Preload
          - DisplacementXYZ
         - mode 1; eigr=123.4
          - EigenvectorXYZ
          - Stress
        """
        keys = self._get_nastran_key_order(model)
        assert keys is not None, keys
        #print('keys_order =', keys)

        disp_dict = defaultdict(list)
        stress_dict = defaultdict(list)
        strain_dict = defaultdict(list)
        force_dict = defaultdict(list)
        strain_energy_dict = defaultdict(list)

        header_dict = {}
        key_itime = []


        for key in keys:
            #print('key = %r' % str(key))
            header_dict[(key, 0)] = '; Static'

            formi = []
            form_time = []
            is_data, is_static, is_real, times = self._get_times(model, key)

            ncases_old = icase
            icase = self._fill_op2_oug_oqg(cases, model, key, icase,
                                           disp_dict, header_dict)
            ncases = icase - ncases_old
            #print('ncases=%s icase=%s' % (ncases, icase))
            #assert ncases > 0, ncases

            if ncases:
                # can potentially make a double listing, but we need it
                # eigenvector only cases
                for itime, dt in enumerate(times):
                    new_key = (key, itime)
                    key_itime.append(new_key)

            for itime, dt in enumerate(times):
                ncases_old = icase
                icase = self._fill_op2_stress(
                    cases, model, key, icase, itime,
                    stress_dict, header_dict, is_static)
                icase = self._fill_op2_strain(
                    cases, model, key, icase, itime,
                    strain_dict, header_dict, is_static)
                icase = self._fill_op2_force(
                    cases, model, key, icase, itime,
                    force_dict, header_dict, is_static)
                icase = self._fill_op2_time_centroidal_strain_energy(
                    cases, model, key, icase, itime,
                    strain_energy_dict, header_dict, is_static)
                ncases = icase - ncases_old
                new_key = (key, itime)
                if ncases and new_key not in key_itime:
                    key_itime.append(new_key)
            #icase = self._fill_op2_oug_oqg(cases, model, key, icase,
                                           #disp_dict, header_dict)
            #print('******form_time =', form_time)
            #print(header)

        form_out = []
        is_results = False

        form_resultsi = []
        form_resultsi_subcase = []
        basename = os.path.basename(op2_filename).rstrip()
        form_results = (basename + '-Results', None, form_resultsi)

        if len(key_itime) == 0:
            print('header_dict =', header_dict)
            print('key_itime =', key_itime)
            return form
        key_itime0 = key_itime[0]
        key0 = key_itime0[0]
        # isubcase, analysis_code, sort_method, count, subtitle
        subcase_id_old = key0[0]
        count_old = key0[3]
        subtitle_old = key0[4]

        for key, itime in key_itime:
            #print('key =', key)
            subcase_id = key[0]
            count = key[3]
            subtitle = key[4]
            if subcase_id != subcase_id_old or subtitle != subtitle_old:
                count_str = '' if count == 0 else ' ; opt_count=%s' % count_old
                subcase_str = 'Subcase %s; %s%s' % (subcase_id_old, subtitle_old, count_str)
                res = (
                    subcase_str.rstrip('; '),
                    None,
                    form_resultsi_subcase
                )
                form_resultsi.append(res)
                form_resultsi_subcase = []
                subcase_id_old = subcase_id
                subtitle_old = subtitle
                count_old = count


            try:
                header = header_dict[(key, itime)]
            except KeyError:
                msg = 'keys =\n'
                for keyi in header_dict:
                    msg += '  %s\n' % str(keyi)
                print(msg.rstrip())
                print('expected = (%s, %r)\n' % (str(key), itime))
                continue
                #raise KeyError(msg)
            try:
                header = header.strip()
            except:
                print('header =', header)
                raise


            form_outi = []
            form_out = (header, None, form_outi)
            disp_formi = disp_dict[(key, itime)]
            stress_formi = stress_dict[(key, itime)]
            strain_formi = strain_dict[(key, itime)]
            force_formi = force_dict[(key, itime)]
            strain_energy_formi = strain_energy_dict[(key, itime)]
            if disp_formi:
                form_outi += disp_formi
                #form_outi.append(('Disp', None, disp_formi))
            if stress_formi:
                form_outi.append(('Stress', None, stress_formi))
                is_results = True
            if strain_formi:
                form_out[2].append(('Strain', None, strain_formi))
                is_results = True
            if force_formi:
                form_out[2].append(('Force', None, force_formi))
                is_results = True
            if strain_energy_formi:
                form_out[2].append(('Strain Energy', None, strain_energy_formi))
                is_results = True

            if form_outi:
                is_results = True
                form_resultsi_subcase.append(form_out)
                #break

        if subcase_id:
            count_str = '' if count == 0 else ' ; opt_count=%s' % count_old
            subcase_str = 'Subcase %s; %s%s' % (subcase_id, subtitle, count_str)
            res = (
                subcase_str.strip('; '),
                None,
                form_resultsi_subcase
            )
            form_resultsi.append(res)
            assert len(form_out) > 0, form_out
            form_resultsi_subcase = []

        if is_results:
            form.append(form_results)
            assert len(form_out) > 0, form_out
            #print('formi =', formi)
            #print('form_out =', form_out)
        #print('form_resultsi =', form_resultsi)
        #print('form_results =', form_results)
            #print(form)
        #if len(formi):
            #form.append(form0)
        #print(form)
        #aa
        return form


    def _get_nastran_key_order(self, model):
        displacement_like = [
            model.displacements,
            model.velocities,
            model.accelerations,
            model.eigenvectors,
            model.spc_forces,
            model.mpc_forces,

            # untested
            model.load_vectors,
            model.applied_loads,
            model.force_vectors,
            #[model.grid_point_forces, 'GridPointForces'],  # TODO: this is buggy...
        ]
        temperature_like = [
            model.temperatures,
        ]
        stress = [
            model.crod_stress, model.conrod_stress, model.ctube_stress,
            model.cbar_stress,
            model.cbeam_stress,

            model.ctria3_stress, model.cquad4_stress,
            model.ctria6_stress, model.cquad8_stress,
            model.ctriar_stress, model.cquadr_stress,

            model.ctria3_composite_stress, model.cquad4_composite_stress,
            model.ctria6_composite_stress, model.cquad8_composite_stress,

            model.ctetra_stress, model.cpenta_stress, model.chexa_stress,
        ]
        strain = [
            model.crod_strain, model.conrod_strain, model.ctube_strain,
            model.cbar_strain,
            model.cbeam_strain,

            model.ctria3_strain, model.cquad4_strain,
            model.ctria6_strain, model.cquad8_strain,
            model.ctriar_strain, model.cquadr_strain,

            model.ctria3_composite_strain, model.cquad4_composite_strain,
            model.ctria6_composite_strain, model.cquad8_composite_strain,

            model.ctetra_strain, model.cpenta_strain, model.chexa_strain,
        ]
        strain_energy = [
            #model.strain_energy,
            model.cquad4_strain_energy, model.cquad8_strain_energy,
            model.cquadr_strain_energy, model.cquadx_strain_energy,
            model.ctria3_strain_energy, model.ctria6_strain_energy,
            model.ctriar_strain_energy, model.ctriax_strain_energy,
            model.ctriax6_strain_energy,
            model.ctetra_strain_energy, model.cpenta_strain_energy,
            model.chexa_strain_energy,
            model.crod_strain_energy, model.ctube_strain_energy,
            model.conrod_strain_energy,
            model.cbar_strain_energy, model.cbeam_strain_energy,
            model.cgap_strain_energy, model.cbush_strain_energy,
            model.celas1_strain_energy, model.celas2_strain_energy,
            model.celas3_strain_energy, model.celas4_strain_energy,
            model.cdum8_strain_energy, model.dmig_strain_energy,
            model.cbend_strain_energy,
            model.genel_strain_energy, model.cshear_strain_energy,
        ]

        result_groups = [
            displacement_like, temperature_like, stress, strain, strain_energy,
        ]

        nids = self.node_ids
        eids = self.element_ids
        keys_order = []
        # model = OP2()

        # subcase_ids = model.subcase_key.keys()

        #self.iSubcaseNameMap[self.isubcase] = [self.subtitle, self.analysis_code, self.label]
        subcase_ids = list(model.iSubcaseNameMap.keys())
        subcase_ids.sort()
        #print('subcase_ids =', subcase_ids)


        # isubcase, analysis_code, sort_method, count, subtitle
        #(1, 2, 1, 0, 'SUPERELEMENT 0') : result1

        #subcase_ids = subcase_ids
        #print('subcase_idsB =' % subcase_ids)
        for isubcase in sorted(subcase_ids):
            if isubcase == 0:
                # beam_modes
                self.log.error('*****isubcase=0')
                continue
            # value = (analysis_codei, sort_methodi, counti, isubtitle)
            #print('subcase_key =', model.subcase_key)
            keys = model.subcase_key[isubcase]
            #print('keys[%s] =%s' % (isubcase, keys))
            key0 = keys[0]

            # this while loop lets us make sure we pull the analysis codes in the expected order
            # TODO: doesn't pull count in the right order
            # TODO: doesn't pull subtitle in right order
            keys2 = deepcopy(keys)
            while keys2:
                key = keys2[-1]
                #print('while keys ->', key)
                (analysis_code, sort_method, count, subtitle) = key
                #assert isubcase == isubcasei, 'isubcase=%s isubcasei=%s' % (isubcase, isubcasei)
                assert analysis_code < 12, analysis_code
                for ianalysis_code in range(12):
                    keyi = (ianalysis_code, sort_method, count, subtitle)
                    if keyi in keys2:
                        #print(keyi)
                        keyi2 = (isubcase, ianalysis_code, sort_method, count, subtitle)
                        #print(keyi2)
                        keys_order.append(keyi2)
                        keys2.remove(keyi)
                #keys_order += keys
        return keys_order

    def _fill_gpforces(self, model):
        pass
        #[model.grid_point_forces, 'GridPointForces'],  # TODO: this is not really an OUG table

    def _fill_op2_oug_oqg(self, cases, model, key, icase,
                          form_dict, header_dict):
        """
        loads the nodal dispalcements/velocity/acceleration/eigenvector/spc/mpc forces
        """
        new_cases = True
        nnodes = self.nNodes
        displacement_like = [
            # slot, name, deflects
            (model.displacements, 'Displacement', True),
            (model.velocities, 'Velocity', False),
            (model.accelerations, 'Acceleration', False),
            (model.eigenvectors, 'Eigenvectors', True),
            (model.spc_forces, 'SPC Forces', False),
            (model.mpc_forces, 'MPC Forces', False),

            (model.load_vectors, 'LoadVectors', False),
            (model.applied_loads, 'AppliedLoads', False),
            (model.force_vectors, 'ForceVectors', False),
        ]
        temperature_like = [
            (model.temperatures, 'Temperature'),
        ]
        nids = self.node_ids
        for (result, name, deflects) in displacement_like:
            if key not in result:
                continue

            title = name + 'XYZ'
            case = result[key]
            subcase_idi = case.isubcase
            if not hasattr(case, 'data'):
                print('str(%s) has no data...' % case.__class.__name__)
                continue
            if not case.is_real():
                print('complex results not supported...')
                continue
            # transient
            if case.nonlinear_factor is not None:
                #code_name = case.data_code['name']
                has_cycle = hasattr(case, 'mode_cycle')
            else:
                has_cycle = False
                code_name = None
            assert case.is_sort1(), case.is_sort1()

            itime0 = 0
            t1 = case.data[itime0, :, 0]
            ndata = t1.shape[0]
            if nnodes != ndata:
                nidsi = case.node_gridtype[:, 0]
                assert len(nidsi) == nnodes
                j = np.searchsorted(nids, nidsi)  # searching for nidsi

                try:
                    if not np.allclose(nids[j], nidsi):
                        msg = 'nids[j]=%s nidsi=%s' % (nids[j], nidsi)
                        raise RuntimeError(msg)
                except IndexError:
                    msg = 'node_ids = %s\n' % list(nids)
                    msg += 'nidsi in disp = %s\n' % list(nidsi)
                    raise IndexError(msg)

            t123 = case.data[:, :, :3]
            if nnodes != ndata:
                t123i = np.zeros((nnodes, 3), dtype='float32')
                t123i[j, :] = t123
                t123 = t123i
            tnorm = norm(t123, axis=1)
            assert len(tnorm) == t123.shape[0]
            ntimes = case.ntimes
            titles = []
            scales = []
            headers = []
            if deflects:
                nastran_res = NastranDisplacementResults(subcase_idi, titles, headers,
                                                         self.xyz_cid0, t123, tnorm,
                                                         scales, deflects=deflects,
                                                         uname='NastranResult')

                for itime in range(ntimes):
                    dt = case._times[itime]
                    header = self._get_nastran_header(case, dt, itime)
                    header_dict[(key, itime)] = header

                    tnorm_abs_max = tnorm.max()
                    if tnorm_abs_max == 0.0:
                        scale = self.displacement_scale_factor
                    else:
                        scale = self.displacement_scale_factor / tnorm_abs_max
                    scale = self.dim_max / tnorm_abs_max * 0.25
                    scales.append(scale)
                    titles.append(title)
                    headers.append(header)
                    cases[icase] = (nastran_res, (itime, title))  # do I keep this???
                    formii = (title, icase, [])
                    form_dict[(key, itime)].append(formii)
                    icase += 1
                nastran_res.save_defaults()
            else:
                for itime in range(ntimes):
                    dt = case._times[itime]
                    header = self._get_nastran_header(case, dt, itime)
                    header_dict[(key, itime)] = header

                    loads = case.data[itime, :, :]
                    nxyz = norm(loads[:, :3], axis=1)
                    assert len(nxyz) == nnodes, 'len(nxyz)=%s nnodes=%s' % (
                        len(nxyz), nnodes)

                    #cases[(subcase_id, icase, word + 'XX', 1, 'node', '%.3f')] = oxx
                    tx_res = GuiResult(subcase_idi, header=name + 'Tx', title=name + 'Tx',
                                       location='node', scalar=loads[:, 0])
                    ty_res = GuiResult(subcase_idi, header=name + 'Ty', title=name + 'Ty',
                                       location='node', scalar=loads[:, 1])
                    tz_res = GuiResult(subcase_idi, header=name + 'Tz', title=name + 'Tz',
                                       location='node', scalar=loads[:, 2])
                    rx_res = GuiResult(subcase_idi, header=name + 'Rx', title=name + 'Rx',
                                       location='node', scalar=loads[:, 3])
                    ry_res = GuiResult(subcase_idi, header=name + 'Ry', title=name + 'Ry',
                                       location='node', scalar=loads[:, 4])
                    rz_res = GuiResult(subcase_idi, header=name + 'Rz', title=name + 'Rz',
                                       location='node', scalar=loads[:, 5])
                    txyz_res = GuiResult(subcase_idi, header=name + 'Txyz',
                                         title=name + name + 'Txyz', location='node', scalar=nxyz)

                    cases[icase] = (tx_res, (0, name + 'Tx'))
                    cases[icase + 1] = (ty_res, (0, name + 'Ty'))
                    cases[icase + 2] = (tz_res, (0, name + 'Tz'))
                    cases[icase + 3] = (rx_res, (0, name + 'Rx'))
                    cases[icase + 4] = (ry_res, (0, name + 'Ry'))
                    cases[icase + 5] = (rz_res, (0, name + 'Rz'))
                    cases[icase + 6] = (txyz_res, (0, name  + 'Txyz'))

                    form_dict[(key, itime)].append((name + 'Tx', icase, []))
                    form_dict[(key, itime)].append((name + 'Ty', icase + 1, []))
                    form_dict[(key, itime)].append((name + 'Tz', icase + 2, []))
                    form_dict[(key, itime)].append((name + 'Rx', icase + 3, []))
                    form_dict[(key, itime)].append((name + 'Ry', icase + 4, []))
                    form_dict[(key, itime)].append((name + 'Rz', icase + 5, []))
                    form_dict[(key, itime)].append((name + 'Txyz', icase + 6, []))
                    icase += 7

        for (result, name) in temperature_like:
            if key not in result:
                continue
            case = result[key]
            subcase_idi = case.isubcase
            if not hasattr(case, 'data'):
                continue

            ntimes = case.ntimes
            for itime in range(ntimes):
                dt = case._times[itime]
                header = self._get_nastran_header(case, dt, itime)
                header_dict[(key, itime)] = header

                loads = case.data[itime, :, :]
                nxyz = norm(loads[:, :3], axis=1)
                assert len(nxyz) == nnodes, 'len(nxyz)=%s nnodes=%s' % (
                    len(nxyz), nnodes)

                temp_res = GuiResult(subcase_idi, header=name, title=name + name,
                                     location='node', scalar=loads[:, 0])
                cases[icase] = (temp_res, (0, name))
                form_dict[(key, itime)].append((name, icase, []))
                icase += 1
        return icase

    def clear_nastran(self):
        """cleans up variables specific to Nastran"""
        self.eid_map = {}
        self.nid_map = {}
        self.eid_to_nid_map = {}
        self.element_ids = None
        self.node_ids = None

    def _fill_op2_force(self, cases, model, key, icase, itime,
                        form_dict, header_dict, is_static):
        """creates the force plots"""
        #assert isinstance(key, int), key
        assert isinstance(icase, int), icase
        assert isinstance(form_dict, dict), form_dict
        icase = self._fill_op2_time_centroidal_force(
            cases, model, key, icase, itime,
            form_dict, header_dict, is_static)
        return icase

    def _fill_op2_stress(self, cases, model, key, icase, itime,
                         form_dict, header_dict, is_static, is_stress=True):
        """creates the stress plots"""
        assert isinstance(icase, int), icase
        assert isinstance(form_dict, dict), form_dict
        icase = self._fill_op2_time_centroidal_stress(
            cases, model, key, icase, itime, form_dict, header_dict,
            is_static, is_stress=is_stress)
        return icase

    def _fill_op2_strain(self, cases, model, key, icase, itime,
                         form_dict, header_dict, is_static):
        """creates the strain plots"""
        return self._fill_op2_stress(cases, model, key, icase, itime,
                                     form_dict, header_dict,
                                     is_static, is_stress=False)

    def _get_times(self, model, isubcase):
        """
        Get the times/frequencies/eigenvalues/loadsteps used on a given
        subcase
        """
        table_types = model.get_table_types()
        is_real = True
        is_data = False
        is_static = False
        times = None
        #print('isubcase =', isubcase)
        #print('table_types =', table_types)
        #print('model.eigenvectors.keys() =', model.eigenvectors.keys())
        for table_type in table_types:
            if not hasattr(model, table_type):
                print('no table_type=%s' % table_type)
                continue
            table = getattr(model, table_type)
            if len(table) == 0:
                continue
            #print(table)
            if isubcase in table:
                is_data = True
                case = table[isubcase]
                #print(case)
                is_real = case.is_real()
                if case.nonlinear_factor is not None:
                    times = case._times
                    is_static = False
                else:
                    is_static = True
                    times = np.zeros(1, dtype='int32')
                #print('times = ', times)
                break
                #return is_data, is_static, is_real, times
        #print('isubcase =', isubcase)
        return is_data, is_static, is_real, times

    def _get_stress_times(self, model, isubcase):
        table_types = self._get_stress_table_types()
        is_real = True
        is_data = False
        is_static = False
        times = None
        for table_type in table_types:
            if not hasattr(model, table_type):
                # print('no table_type=%s' % table_type)
                continue
            table = getattr(model, table_type)
            if isubcase in table:
                is_data = True
                case = table[isubcase]
                is_real = case.is_real()
                if case.nonlinear_factor is not None:
                    times = case._times
                    is_static = False
                else:
                    is_static = True
                    times = np.zeros(1, dtype='int32')
                break
                #return is_data, is_static, is_real, times
        return is_data, is_static, is_real, times

    def _get_stress_table_types(self):
        """
        Gets the list of Nastran stress objects that the GUI supports
        """
        table_types = [
            # OES - tCode=5 thermal=0 s_code=0,1 (stress/strain)
            # OES - CELAS1/CELAS2/CELAS3/CELAS4 stress
            'celas1_stress',
            'celas2_stress',
            'celas3_stress',
            'celas4_stress',

            # OES - CELAS1/CELAS2/CELAS3/CELAS4 strain
            'celas1_strain',
            'celas2_strain',
            'celas3_strain',
            'celas4_strain',

            # OES - isotropic CROD/CONROD/CTUBE stress
            'crod_stress',
            'conrod_stress',
            'ctube_stress',

            # OES - isotropic CROD/CONROD/CTUBE strain
            'crod_strain',
            'conrod_strain',
            'ctube_strain',

            # OES - isotropic CBAR stress
            'cbar_stress',
            # OES - isotropic CBAR strain
            'cbar_strain',
            # OES - isotropic CBEAM stress
            'cbeam_stress',
            # OES - isotropic CBEAM strain
            'cbeam_strain',

            # OES - isotropic CTRIA3/CQUAD4 stress
            'ctria3_stress',
            'cquad4_stress',

            # OES - isotropic CTRIA3/CQUAD4 strain
            'ctria3_strain',
            'cquad4_strain',

            # OES - isotropic CTETRA/CHEXA/CPENTA stress
            'ctetra_stress',
            'chexa_stress',
            'cpenta_stress',

            # OES - isotropic CTETRA/CHEXA/CPENTA strain
            'ctetra_strain',
            'chexa_strain',
            'cpenta_strain',

            # OES - CSHEAR stress
            'cshear_stress',
            # OES - CSHEAR strain
            'cshear_strain',
            # OES - CEALS1 224, CELAS3 225
            'nonlinear_spring_stress',
            # OES - GAPNL 86
            'nonlinear_cgap_stress',
            # OES - CBUSH 226
            'nolinear_cbush_stress',
        ]

        table_types += [
            # OES - CTRIAX6
            'ctriax_stress',
            'ctriax_strain',

            'cbush_stress',
            'cbush_strain',
            'cbush1d_stress_strain',

            # OES - nonlinear CROD/CONROD/CTUBE stress
            'nonlinear_rod_stress',
            'nonlinear_rod_strain',

            # OESNLXR - CTRIA3/CQUAD4 stress
            'nonlinear_plate_stress',
            'nonlinear_plate_strain',
            #'hyperelastic_plate_stress',
            'hyperelastic_cquad4_strain',

            # OES - composite CTRIA3/CQUAD4 stress
            'cquad4_composite_stress',
            'cquad8_composite_stress',
            'ctria3_composite_stress',
            'ctria6_composite_stress',

            'cquad4_composite_strain',
            'cquad8_composite_strain',
            'ctria3_composite_strain',
            'ctria6_composite_strain',

            # OGS1 - grid point stresses
            'grid_point_stresses',        # tCode=26
            'grid_point_volume_stresses',  # tCode=27
        ]
        return table_types

    def _get_nastran_header(self, case, dt, itime):
        #if case is None:
            #return None
        try:
            code_name = case.data_code['name']
        except KeyError:
            return 'Static'

        if isinstance(dt, float):
            header = ' %s = %.4E' % (code_name, dt)
        else:
            header = ' %s = %i' % (code_name, dt)

        if hasattr(case, 'mode_cycle'):
            freq = case.eigrs[itime]
            #msg.append('%16s = %13E\n' % ('EIGENVALUE', freq))
            cycle = np.sqrt(np.abs(freq)) / (2. * np.pi)
            header += '; freq=%g' % cycle
        elif hasattr(case, 'eigrs'):
            freq = case.eigrs[itime]
            #msg.append('%16s = %13E\n' % ('EIGENVALUE', freq))
            cycle = np.sqrt(np.abs(freq)) / (2. * np.pi)
            header += '; freq=%g' % cycle
        return header.strip('; ')

    def _fill_op2_time_centroidal_strain_energy(self, cases, model,
                                                key, icase, itime,
                                                form_dict, header_dict, is_static):
        """
        Creates the time accurate strain energy objects for the pyNastranGUI
        """
        oxx = np.zeros(self.nElements, dtype='float32')
        oyy = np.zeros(self.nElements, dtype='float32')
        ozz = np.zeros(self.nElements, dtype='float32')
        fmt = '%g'
        header = ''
        form0 = ('Element Strain Energy', None, [])

        #op2.strain_energy[1]
          #type=StrainEnergyObject ntimes=3 nelements=16
          #energy, percent, density
          #modes = [1, 2, 3]
        case = None
        subcase_id = key[2]
        if not hasattr(model, 'strain_energy'):
            return icase
        if key in model.strain_energy:
            #print('key =', key)
            ese = model.strain_energy[key]
            #print(ese)
            times = sorted(ese.energy.keys())  # TODO: not vectorized
            #assert times[itime] == dt, 'actual=%s expected=%s' % (times[itime], dt)
            dt = times[itime]

            try:
                header = self._get_nastran_header(case, dt, itime)
                header_dict[(key, itime)] = header
            except AttributeError:
                pass

            if is_static:
                percent = ese.percent
                energy = ese.energy
                density = ese.density
            else:
                percent = ese.percent[dt]
                energy = ese.energy[dt]
                density = ese.density[dt]
            for eid, percenti in sorted(iteritems(percent)):
                if eid not in self.eid_map:
                    continue
                i = self.eid_map[eid]
                oxx[i] = energy[eid]
                oyy[i] = percenti
                ozz[i] = density[eid]

            case = ese
            fmt = '%.4f'
            cases[(subcase_id, icase, 'StrainEnergy', 1, 'centroid', fmt, header)] = oxx
            form_dict[(key, itime)].append(('StrainEnergy', icase, []))
            icase += 1

            cases[(subcase_id, icase, 'PercentOfTotal', 1, 'centroid', fmt, header)] = oyy
            form_dict[(key, itime)].append(('PercentOfTotal', icase, []))
            icase += 1

            cases[(subcase_id, icase, 'Density', 1, 'centroid', fmt, header)] = ozz
            form_dict[(key, itime)].append(('Density', icase, []))
            icase += 1
        return icase

    #icase = self._fill_op2_time_centroidal_force(
        #cases, model, subcase_id, icase, itime, form_dict,
        #is_static)

    def _fill_op2_time_centroidal_force(self, cases, model,
                                        key, icase, itime,
                                        form_dict, header_dict, is_static):
        """
        Creates the time accurate strain energy objects for the pyNastranGUI
        """
        #new_cases = True
        nelements = self.nElements
        fx = np.zeros(nelements, dtype='float32') # axial
        fy = np.zeros(nelements, dtype='float32') # shear_y
        fz = np.zeros(nelements, dtype='float32') # shear_z

        rx = np.zeros(nelements, dtype='float32') # torque
        ry = np.zeros(nelements, dtype='float32') # bending_y
        rz = np.zeros(nelements, dtype='float32') # bending_z

        is_element_on = np.zeros(nelements, dtype='float32') # torque
        fmt = '%g'
        header = ''
        form0 = ('Force', None, [])

        case = None
        found_force = False
        for res_type in (model.conrod_force, model.crod_force, model.ctube_force):
            if key in res_type:
                found_force = True
                case = res_type[key]
                if case.is_complex():
                    continue
                data = case.data
                if case.nonlinear_factor is None:
                    ntimes = data.shape[:1]
                    eids = case.element
                    dt = case._times[itime]
                    header = self._get_nastran_header(case, dt, itime)
                    header_dict[(key, itime)] = header
                    #eids_to_find = intersect1d(self.element_ids, eids)
                    i = np.searchsorted(self.element_ids, eids)
                    assert np.array_equal(self.element_ids[i], eids)
                    fxi = data[itime, :, 0]
                    rxi = data[itime, :, 1]
                    if fxi.size != i.size:
                        msg = 'fx.size=%s i.size=%s fx=%s eids_to_find=%s' % (
                            fxi.size, i.size, fxi, eids)
                        raise RuntimeError(msg)
                    fx[i] = fxi
                    rx[i] = rxi
                    is_element_on[i] = 1.
                else:
                    continue

        if key in model.cbar_force:
            found_force = True
            ## CBAR-34
            case = model.cbar_force[key]
            if case.is_real():
                eids = case.element
                i = np.searchsorted(self.element_ids, eids)
                is_element_on[i] = 1.

                dt = case._times[itime]
                header = self._get_nastran_header(case, dt, itime)
                header_dict[(key, itime)] = header

                #[bending_moment_a1, bending_moment_a2, bending_moment_b1, bending_moment_b2,
                # shear1, shear2, axial, torque]
                #fx[i] = case.data[:, :, 6]
                #fy[i] = case.data[:, :, 4]
                #fz[i] = case.data[:, :, 5]

                if i.size == 1:
                    rxi = case.data[itime, :, 7].max()
                    ryi = np.vstack([case.data[itime, :, 0], case.data[itime, :, 2]]).max()
                    rzi = np.vstack([case.data[itime, :, 1], case.data[itime, :, 3]]).max()
                else:
                    rxi = case.data[itime, :, 7]#.max(axis=0)
                    ryi = np.vstack([case.data[itime, :, 0], case.data[itime, :, 2]]).max(axis=0)
                    rzi = np.vstack([case.data[itime, :, 1], case.data[itime, :, 3]]).max(axis=0)
                    rzv = rzi

                    # rza = array([case.data[itime, :, 1], case.data[itime, :, 3]])#.max(axis=0)
                    # rzh = hstack([case.data[itime, :, 1], case.data[itime, :, 3]])#.max(axis=0)
                    # print(rzv.shape, rzv.shape, rzv.shape)
                assert rxi.size == i.size, 'rx.size=%s i.size=%s rx=%s' % (rxi.size, i.size, rxi)
                assert ryi.size == i.size, 'ry.size=%s i.size=%s ry=%s' % (ryi.size, i.size, ryi)
                assert rzi.size == i.size, 'rz.size=%s i.size=%s rz=%s' % (rzi.size, i.size, rzi)

                rx[i] = rxi
                ry[i] = ryi
                rz[i] = rzi

        if key in model.cbar_force_10nodes:
            found_force = True
            ## CBAR-100
            case = model.cbar_force_10nodes[key]
            eids = case.element
            ueids = np.unique(eids)

            dt = case._times[itime]
            header = self._get_nastran_header(case, dt, itime)
            header_dict[(key, itime)] = header

            j = np.searchsorted(self.element_ids, ueids)
            is_element_on[j] = 1.
            di = j[1:-1] - j[0:-2]
            if di.max() != 2:
                print('di =', np.unique(di))
                # [station, bending_moment1, bending_moment2, shear1, shear2, axial, torque]
                ii = 0
                eid_old = eids[0]
                fxi = defaultdict(list)
                fyi = defaultdict(list)
                fzi = defaultdict(list)
                rxi = defaultdict(list)
                ryi = defaultdict(list)
                rzi = defaultdict(list)
                for ii, eid in enumerate(eids):
                    fxi[eid].append(case.data[:, ii, 5])
                    fyi[eid].append(case.data[:, ii, 3])
                    fzi[eid].append(case.data[:, ii, 4])

                    rxi[eid].append(case.data[:, ii, 6])
                    ryi[eid].append(case.data[:, ii, 1])
                    rzi[eid].append(case.data[:, ii, 2])
                    #if eidi == eid_old:
                    #    fx[ii] = array([case.data[:, j, 5], case.data[:, j, 5]]).max(axis=0)
                    #else:
                for ii, eidi in zip(j, eids[j]):
                    fx[ii] = max(fxi[eidi])
                    fy[ii] = max(fyi[eidi])
                    fz[ii] = max(fyi[eidi])
                    rx[ii] = max(rxi[eidi])
                    ry[ii] = max(ryi[eidi])
                    rz[ii] = max(rzi[eidi])
            else:
                # [station, bending_moment1, bending_moment2, shear1, shear2, axial, torque]
                neids = len(np.unique(eids)) * 2
                if len(eids) != len(np.unique(eids)) * 2:
                    msg = 'CBAR-100 Error: len(eids)=%s neids=%s' % (len(eids), neids)
                    raise RuntimeError(msg)
                fx[i] = np.array(
                    [case.data[itime, ::-1, 5],
                     case.data[itime, 1::-1, 5]]).max(axis=0)
                fy[i] = np.array(
                    [case.data[itime, ::-1, 3],
                     case.data[itime, 1::-1, 3]]).max(axis=0)
                fz[i] = np.array(
                    [case.data[itime, ::-1, 4],
                     case.data[itime, 1::-1, 4]]).max(axis=0)
                rx[i] = np.array(
                    [case.data[itime, ::-1, 6],
                     case.data[itime, 1::-1, 6]]).max(axis=0)
                ry[i] = np.array(
                    [case.data[itime, ::-1, 1],
                     case.data[itime, 1::-1, 1]]).max(axis=0)
                rz[i] = np.array(
                    [case.data[itime, ::-1, 2],
                     case.data[itime, 1::-1, 2]]).max(axis=0)

        subcase_id = key[2]
        if found_force:
            fmt = '%.4f'
            # header = self._get_nastran_header(case, dt, itime)

            #num_on = nelements
            num_off = 0
            if itime == 0 and is_element_on.min() == 0.0:
                ioff = np.where(is_element_on == 0)[0]
                num_off = len(ioff)
                print('force_eids_off = %s; n=%s' % (self.element_ids[ioff], num_off))
                self.log_error('force_eids_off = %s; n=%s' % (self.element_ids[ioff], num_off))
                cases[(subcase_id, icase, 'Force\nIsElementOn', 1, 'centroid', '%i', header)] = is_element_on
                form_dict[(key, itime)].append(('Force - IsElementOn', icase, []))
                #num_on -= num_off
                icase += 1


            if fx.min() != fx.max() or rx.min() != rx.max() and not num_off == nelements:
                fx_res = GuiResult(subcase_id, header='Axial', title='Axial',
                                   location='centroid', scalar=fx)
                fy_res = GuiResult(subcase_id, header='ShearY', title='ShearY',
                                   location='centroid', scalar=fy)
                fz_res = GuiResult(subcase_id, header='ShearZ', title='ShearZ',
                                   location='centroid', scalar=fz)
                mx_res = GuiResult(subcase_id, header='Torsion', title='Torsion',
                                   location='centroid', scalar=rx)
                my_res = GuiResult(subcase_id, header='BendingY', title='BendingY',
                                   location='centroid', scalar=ry)
                mz_res = GuiResult(subcase_id, header='BendingZ', title='BendingZ',
                                   location='centroid', scalar=rz)
                cases[icase] = (fx_res, (subcase_id, 'Axial'))
                cases[icase + 1] = (fy_res, (subcase_id, 'ShearY'))
                cases[icase + 2] = (fz_res, (subcase_id, 'ShearZ'))
                cases[icase + 3] = (mx_res, (subcase_id, 'Torsion'))
                cases[icase + 4] = (my_res, (subcase_id, 'BendingY'))
                cases[icase + 5] = (mz_res, (subcase_id, 'BendingZ'))

                form_dict[(key, itime)].append(('Axial', icase, []))
                form_dict[(key, itime)].append(('ShearY', icase + 1, []))
                form_dict[(key, itime)].append(('ShearZ', icase + 2, []))
                form_dict[(key, itime)].append(('Torque', icase + 3, []))
                form_dict[(key, itime)].append(('BendingY', icase + 4, []))
                form_dict[(key, itime)].append(('BendingZ', icase + 5, []))
                icase += 6

                is_axial = np.zeros(self.nElements, dtype='int8')
                is_shear_y = np.zeros(self.nElements, dtype='int8')
                is_shear_z = np.zeros(self.nElements, dtype='int8')
                is_torsion = np.zeros(self.nElements, dtype='int8')
                is_bending_y = np.zeros(self.nElements, dtype='int8')
                is_bending_z = np.zeros(self.nElements, dtype='int8')
                is_axial[np.where(np.abs(fx) > 0.0)[0]] = 1
                is_shear_y[np.where(np.abs(fy) > 0.0)[0]] = 1
                is_shear_z[np.where(np.abs(fz) > 0.0)[0]] = 1
                is_torsion[np.where(np.abs(rx) > 0.0)[0]] = 1
                is_bending_y[np.where(np.abs(ry) > 0.0)[0]] = 1
                is_bending_z[np.where(np.abs(rz) > 0.0)[0]] = 1
                #is_bending[where(abs(rx) > 0.0)[0]] = 1

                is_fx_res = GuiResult(subcase_id, header='IsAxial', title='IsAxial',
                                      location='centroid', scalar=is_axial, data_format=fmt)
                is_fy_res = GuiResult(subcase_id, header='IsShearY', title='IsShearY',
                                      location='centroid', scalar=is_shear_y, data_format=fmt)
                is_fz_res = GuiResult(subcase_id, header='IsShearZ', title='IsShearZ',
                                      location='centroid', scalar=is_shear_z, data_format=fmt)
                is_mx_res = GuiResult(subcase_id, header='IsTorsion', title='IsTorsion',
                                      location='centroid', scalar=is_torsion, data_format=fmt)
                is_my_res = GuiResult(subcase_id, header='IsBendingY', title='IsBendingY',
                                      location='centroid', scalar=is_bending_y, data_format=fmt)
                is_mz_res = GuiResult(subcase_id, header='IsBendingZ', title='IsBendingZ',
                                      location='centroid', scalar=is_bending_z, data_format=fmt)

                cases[icase] = (is_fx_res, (subcase_id, 'IsAxial'))
                cases[icase + 1] = (is_fy_res, (subcase_id, 'IsShearY'))
                cases[icase + 2] = (is_fz_res, (subcase_id, 'IsShearZ'))
                cases[icase + 3] = (is_mx_res, (subcase_id, 'IsTorsion'))
                cases[icase + 4] = (is_my_res, (subcase_id, 'IsBendingY'))
                cases[icase + 5] = (is_mz_res, (subcase_id, 'IsBendingZ'))

                form_dict[(key, itime)].append(('IsAxial', icase, []))
                form_dict[(key, itime)].append(('IsShearY', icase + 1, []))
                form_dict[(key, itime)].append(('IsShearZ', icase + 2, []))
                form_dict[(key, itime)].append(('IsTorsion', icase + 3, []))
                form_dict[(key, itime)].append(('IsBendingY', icase + 4, []))
                form_dict[(key, itime)].append(('IsBendingZ', icase + 5, []))
                icase += 6
        return icase

    def _fill_op2_time_centroidal_stress(self, cases, model, key, icase, itime,
                                         form_dict, header_dict, is_static, is_stress=True):
        """
        Creates the time accurate stress objects for the pyNastranGUI
        """
        #new_cases = True
        case = None
        #assert isinstance(subcase_id, int), type(subcase_id)
        assert isinstance(icase, int), icase
        #assert isinstance(itime, int), type(itime)
        assert is_stress in [True, False], is_stress
        eids = self.element_ids
        assert len(eids) > 0, eids
        nelements = self.nElements
        dt = None

        is_element_on = np.zeros(nelements, dtype='int8')  # is the element supported
        oxx = np.zeros(nelements, dtype='float32')
        oyy = np.zeros(nelements, dtype='float32')
        ozz = np.zeros(nelements, dtype='float32')

        txy = np.zeros(nelements, dtype='float32')
        tyz = np.zeros(nelements, dtype='float32')
        txz = np.zeros(nelements, dtype='float32')

        max_principal = np.zeros(nelements, dtype='float32')  # max
        mid_principal = np.zeros(nelements, dtype='float32')  # mid
        min_principal = np.zeros(nelements, dtype='float32')  # min
        ovm = np.zeros(nelements, dtype='float32')

        vm_word = None
        if is_stress:
            rods = [model.crod_stress, model.conrod_stress, model.ctube_stress,]
        else:
            rods = [model.crod_strain, model.conrod_strain, model.ctube_strain,]

        for result in rods:
            if key not in result:
                continue

            case = result[key]
            if case.is_complex():
                continue
            eidsi = case.element
            i = np.searchsorted(eids, eidsi)
            if len(i) != len(np.unique(i)):
                msg = 'irod=%s is not unique\n' % str(i)
                print('eids = %s\n' % str(list(eids)))
                print('eidsi = %s\n' % str(list(eidsi)))
                raise RuntimeError(msg)

            is_element_on[i] = 1
            dt = case._times[itime]
            header = self._get_nastran_header(case, dt, itime)
            header_dict[(key, itime)] = header

            # data=[1, nnodes, 4] where 4=[axial, SMa, torsion, SMt]
            oxx[i] = case.data[itime, :, 0]
            txy[i] = case.data[itime, :, 2]
            ovm[i] = np.sqrt(oxx[i]**2 + 3*txy[i]**2) # plane stress
            # max_principal[i] = sqrt(oxx[i]**2 + txy[i]**2)
            # min_principal[i] = max_principal[i] - 2 * txy[i]
            # simplification of:
            #   eig(A) = [oxx, txy]
            #            [txy, 0.0]
            # per Equation 7: http://www.soest.hawaii.edu/martel/Courses/GG303/Eigenvectors.pdf
            max_principal[i] = (oxx[i] + np.sqrt(oxx[i]**2 + 4 * txy[i]**2)) / 2.
            min_principal[i] = (oxx[i] - np.sqrt(oxx[i]**2 + 4 * txy[i]**2)) / 2.
        del rods


        if is_stress:
            bars = model.cbar_stress
        else:
            bars = model.cbar_strain

        if key in bars:
            case = bars[key]
            if case.is_complex():
                pass
            else:
                dt = case._times[itime]
                header = self._get_nastran_header(case, dt, itime)
                header_dict[(key, itime)] = header
                #s1a = case.data[itime, :, 0]
                #s2a = case.data[itime, :, 1]
                #s3a = case.data[itime, :, 2]
                #s4a = case.data[itime, :, 3]

                axial = case.data[itime, :, 4]
                smaxa = case.data[itime, :, 5]
                smina = case.data[itime, :, 6]
                #MSt = case.data[itime, :, 7]

                #s1b = case.data[itime, :, 8]
                #s2b = case.data[itime, :, 9]
                #s3b = case.data[itime, :, 10]
                #s4b = case.data[itime, :, 11]

                smaxb = case.data[itime, :, 12]
                sminb = case.data[itime, :, 13]
                #MSc   = case.data[itime, :, 14]

                eidsi = case.element # [:, 0]

                i = np.searchsorted(eids, eidsi)
                if len(i) != len(np.unique(i)):
                    print('ibar = %s' % i)
                    print('eids = %s' % eids)
                    msg = 'ibar=%s is not unique' % str(i)
                    raise RuntimeError(msg)

                is_element_on[i] = 1.
                oxx[i] = axial

                ## TODO :not sure if this block is general for multiple CBAR elements
                samax = np.amax([smaxa, smaxb], axis=0)
                samin = np.amin([smaxa, smaxb], axis=0)
                assert len(samax) == len(i), len(samax)
                assert len(samin) == len(i)
                savm = np.amax(np.abs(
                    [smina, sminb,
                     smaxa, smaxb, axial]), axis=0)

                max_principal[i] = samax
                min_principal[i] = samin
                ovm[i] = savm
                del axial, smaxa, smina, smaxb, sminb, eidsi, i, samax, samin, savm
        del bars


        if is_stress:
            bars2 = model.cbar_stress_10nodes
        else:
            bars2 = model.cbar_strain_10nodes

        if key in bars2:
            case = bars[key]
            if case.is_complex():
                pass
            else:
                dt = case._times[itime]
                header = self._get_nastran_header(case, dt, itime)
                header_dict[(key, itime)] = header
                #s1a = case.data[itime, :, 0]
                #s2a = case.data[itime, :, 1]
                #s3a = case.data[itime, :, 2]
                #s4a = case.data[itime, :, 3]

                axial = case.data[itime, :, 4]
                smaxa = case.data[itime, :, 5]
                smina = case.data[itime, :, 6]
                #MSt = case.data[itime, :, 7]

                #s1b = case.data[itime, :, 8]
                #s2b = case.data[itime, :, 9]
                #s3b = case.data[itime, :, 10]
                #s4b = case.data[itime, :, 11]

                smaxb = case.data[itime, :, 12]
                sminb = case.data[itime, :, 13]
                #MSc   = case.data[itime, :, 14]

                eidsi = case.element # [:, 0]

                i = np.searchsorted(eids, eidsi)
                if len(i) != len(np.unique(i)):
                    print('ibar = %s' % i)
                    print('eids = %s' % eids)
                    msg = 'ibar=%s is not unique' % str(i)
                    raise RuntimeError(msg)

                is_element_on[i] = 1.
                oxx[i] = axial

                ## TODO :not sure if this block is general for multiple CBAR elements
                samax = np.amax([smaxa, smaxb], axis=0)
                samin = np.amin([smaxa, smaxb], axis=0)
                assert len(samax) == len(i), len(samax)
                assert len(samin) == len(i)
                savm = np.amax(np.abs(
                    [smina, sminb,
                     smaxa, smaxb, axial]), axis=0)

                max_principal[i] = samax
                min_principal[i] = samin
                ovm[i] = savm
                #del axial, smaxa, smina, smaxb, sminb, eidsi, i, samax, samin, savm
        del bars2


        if is_stress:
            beams = model.cbeam_stress
        else:
            beams = model.cbeam_strain

        if key in beams:
            case = beams[key]
            if case.is_complex():
                pass
            else:
                eidsi = case.element_node[:, 0]
                ueids = np.unique(eidsi)
                #neids = len(ueids)

                # sxc, sxd, sxe, sxf
                # smax, smin, MSt, MSc
                dt = case._times[itime]
                header = self._get_nastran_header(case, dt, itime)
                header_dict[(key, itime)] = header
                sxc = case.data[itime, :, 0]
                sxd = case.data[itime, :, 1]
                sxe = case.data[itime, :, 2]
                sxf = case.data[itime, :, 3]
                smax = case.data[itime, :, 4]
                smin = case.data[itime, :, 5]

                imin = np.searchsorted(eidsi, ueids)
                imax = np.searchsorted(eidsi, ueids, side='right')
                #sxxi = smax[imin:imax]
                for eid, imini, imaxi in zip(ueids, imin, imax):
                    oxxi = 0.
                    smaxi = 0.
                    smini = 0.
                    eid2 = self.eid_map[eid]
                    is_element_on[eid2] = 1.
                    oxxi = max(
                        sxc[imini:imaxi].max(),
                        sxd[imini:imaxi].max(),
                        sxe[imini:imaxi].max(),
                        sxf[imini:imaxi].max(),
                    )
                    smaxi = smax[imini:imaxi].max()
                    smini = smin[imini:imaxi].min()
                    ovmi = max(np.abs(smaxi), np.abs(smini))
                    oxxi = oxx[eid2]
                    max_principal[eid2] = smaxi
                    min_principal[eid2] = smini
                    ovm[eid2] = ovmi
                del eidsi, ueids, sxc, sxd, sxe, sxf, smax, smin, oxxi, smaxi, smini, ovmi
        del beams


        if is_stress:
            plates = [
                model.ctria3_stress, model.cquad4_stress,
                model.ctria6_stress, model.cquad8_stress,
                model.ctriar_stress, model.cquadr_stress,
            ]
        else:
            plates = [
                model.ctria3_strain, model.cquad4_strain,
                model.ctria6_strain, model.cquad8_strain,
                model.ctriar_strain, model.cquadr_strain,
            ]

        for result in plates:
            if key not in result:
                continue
            case = result[key]
            if case.is_complex():
                continue

            if case.is_von_mises():
                vm_word = 'vonMises'
            else:
                vm_word = 'maxShear'

            nnodes_per_element = case.nnodes
            nlayers_per_element = nnodes_per_element * 2  # *2 for every other layer
            eidsi = case.element_node[::nlayers_per_element, 0]  # ::2 is for layer skipping

            i = np.searchsorted(eids, eidsi)
            if len(i) != len(np.unique(i)):
                print('iplate = %s' % i)
                print('eids = %s' % eids)
                print('eidsiA = %s' % case.element_node[:, 0])
                print('eidsiB = %s' % eidsi)
                msg = 'iplate=%s is not unique' % str(i)
                raise RuntimeError(msg)
            #self.data[self.itime, self.itotal, :] = [fd, oxx, oyy,
            #                                         txy, angle,
            #                                         majorP, minorP, ovm]
            is_element_on[i] = 1.
            ntotal = case.data.shape[1]  # (ndt, ntotal, nresults)
            if nlayers_per_element == 1:
                j = None
            else:
                j = np.arange(ntotal)[::nlayers_per_element]

            #self.data[self.itime, self.itotal, :] = [fd, oxx, oyy,
            #                                         txy, angle,
            #                                         majorP, minorP, ovm]
            dt = case._times[itime]
            header = self._get_nastran_header(case, dt, itime)
            header_dict[(key, itime)] = header
            oxxi = case.data[itime, j, 1]
            oyyi = case.data[itime, j, 2]
            txyi = case.data[itime, j, 3]
            o1i = case.data[itime, j, 5]
            o3i = case.data[itime, j, 6]
            ovmi = case.data[itime, j, 7]

            for inode in range(1, nlayers_per_element):
                #print('%s - ilayer = %s' % (case.element_name, inode))
                oxxi = np.amax(np.vstack([oxxi, case.data[itime, j + inode, 1]]), axis=0)
                oyyi = np.amax(np.vstack([oyyi, case.data[itime, j + inode, 2]]), axis=0)
                txyi = np.amax(np.vstack([txyi, case.data[itime, j + inode, 3]]), axis=0)
                o1i = np.amax(np.vstack([o1i, case.data[itime, j + inode, 5]]), axis=0)
                o3i = np.amin(np.vstack([o3i, case.data[itime, j + inode, 6]]), axis=0)
                ovmi = np.amax(np.vstack([ovmi, case.data[itime, j + inode, 7]]), axis=0)
                assert len(oxxi) == len(j)

            oxx[i] = oxxi
            oyy[i] = oyyi
            txy[i] = txyi
            max_principal[i] = o1i
            min_principal[i] = o3i
            ovm[i] = ovmi


        if is_stress:
            cplates = [
                ('CTRIA3', model.ctria3_composite_stress),
                ('CQUAD4', model.cquad4_composite_stress),
                ('CTRIA6', model.ctria6_composite_stress),
                ('CQUAD8', model.cquad8_composite_stress),
                #model.ctriar_composite_stress,
                #model.cquadr_composite_stress,
            ]
        else:
            cplates = [
                ('CTRIA3', model.ctria3_composite_strain),
                ('CQUAD4', model.cquad4_composite_strain),
                ('CTRIA6', model.ctria6_composite_strain),
                ('CQUAD8', model.cquad8_composite_strain),
                #model.ctriar_composite_strain,
                #model.cquadr_composite_strain,
            ]

        for cell_type, result in cplates:
            if key not in result:
                continue
            case = result[key]
            if case.is_complex():
                continue

            if case.is_von_mises():
                vm_word = 'vonMises'
            else:
                vm_word = 'maxShear'

            dt = case._times[itime]
            header = self._get_nastran_header(case, dt, itime)
            header_dict[(key, itime)] = header
            eidsi = case.element_layer[:, 0]
            layers = case.element_layer[:, 1]
            ntotal = case.data.shape[1]

            #[o11, o22, t12, t1z, t2z, angle, major, minor, max_shear]
            oxxs = case.data[itime, :, 0]
            oyys = case.data[itime, :, 1]
            txys = case.data[itime, :, 2]
            txzs = case.data[itime, :, 3]
            tyzs = case.data[itime, :, 4]
            # angle
            omaxs = case.data[itime, :, 6]
            omins = case.data[itime, :, 7]
            ovms = case.data[itime, :, 8]

            j = 0
            for eid in np.unique(eidsi):
                ieid = np.where(eidsi == eid)[0]
                ieid.sort()
                layersi = layers[ieid]
                eid2 = self.eid_map[eid]
                is_element_on[eid2] = 1.

                oxxi = 0.
                oyyi = 0.
                txyi = 0.
                tyzi = 0.
                txzi = 0.
                omaxi = 0.
                omini = 0.
                ovmi = 0.
                nlayers = len(layersi)
                for ilayer in range(nlayers):
                    oxxi = max(oxxs[j], oxxi)
                    oyyi = max(oyys[j], oyyi)
                    txyi = max(txys[j], txyi)
                    tyzi = max(tyzs[j], tyzi)
                    txzi = max(txzs[j], txzi)

                    omaxi = max(omaxs[j], omaxi)
                    omini = min(omins[j], omini)
                    ovmi = max(ovms[j], ovmi)
                    j += 1

                oxx[eid2] = oxxi
                oyy[eid2] = oyyi
                txy[eid2] = txyi
                tyz[eid2] = tyzi
                txz[eid2] = txzi
                max_principal[eid2] = omaxi
                min_principal[eid2] = omini
                ovm[eid2] = ovmi
            del oxxi, oyyi, txyi, tyzi, txzi, omaxi, omini, ovmi, eid2, j, layers, eidsi
        del cplates


        if is_stress:
            solids = [(model.ctetra_stress),
                      (model.cpenta_stress),
                      (model.chexa_stress),]
        else:
            solids = [(model.ctetra_strain),
                      (model.cpenta_strain),
                      (model.chexa_strain),]

        for result in solids:
            if key not in result:
                continue
            case = result[key]
            if case.is_complex():
                continue

            if case.is_von_mises():
                vm_word = 'vonMises'
            else:
                vm_word = 'maxShear'

            nnodes_per_element = case.nnodes
            eidsi = case.element_cid[:, 0]
            ntotal = len(eidsi)  * nnodes_per_element

            i = np.searchsorted(eids, eidsi)
            if len(i) != len(np.unique(i)):
                print('isolid = %s' % str(i))
                print('eids = %s' % eids)
                print('eidsi = %s' % eidsi)
                assert len(i) == len(np.unique(i)), 'isolid=%s is not unique' % str(i)

            is_element_on[i] = 1
            #self.data[self.itime, self.itotal, :] = [oxx, oyy, ozz,
            #                                         txy, tyz, txz,
            #                                         o1, o2, o3, ovm]

            if nnodes_per_element == 1:
                j = None
            else:
                j = np.arange(ntotal)[::nnodes_per_element]
                ueidsi = np.unique(eidsi)
                assert len(j) == len(ueidsi), 'j=%s ueidsi=%s' % (j, ueidsi)

            dt = case._times[itime]
            header = self._get_nastran_header(case, dt, itime)
            header_dict[(key, itime)] = header
            oxxi = case.data[itime, j, 0]
            oyyi = case.data[itime, j, 1]
            ozzi = case.data[itime, j, 2]
            txyi = case.data[itime, j, 3]
            tyzi = case.data[itime, j, 4]
            txzi = case.data[itime, j, 5]
            o1i = case.data[itime, j, 6]
            o2i = case.data[itime, j, 7]
            o3i = case.data[itime, j, 8]
            ovmi = case.data[itime, j, 9]

            for inode in range(1, nnodes_per_element):
                oxxi = np.amax(np.vstack([oxxi, case.data[itime, j + inode, 0]]), axis=0)
                oyyi = np.amax(np.vstack([oyyi, case.data[itime, j + inode, 1]]), axis=0)
                ozzi = np.amax(np.vstack([ozzi, case.data[itime, j + inode, 2]]), axis=0)
                txyi = np.amax(np.vstack([txyi, case.data[itime, j + inode, 3]]), axis=0)
                tyzi = np.amax(np.vstack([tyzi, case.data[itime, j + inode, 4]]), axis=0)
                txzi = np.amax(np.vstack([txzi, case.data[itime, j + inode, 2]]), axis=0)

                o1i = np.amax(np.vstack([o1i, case.data[itime, j + inode, 6]]), axis=0)
                o2i = np.amax(np.vstack([o2i, case.data[itime, j + inode, 7]]), axis=0)
                o3i = np.amin(np.vstack([o3i, case.data[itime, j + inode, 8]]), axis=0)
                ovmi = np.amax(np.vstack([ovmi, case.data[itime, j + inode, 9]]), axis=0)
                assert len(oxxi) == len(j)

            oxx[i] = oxxi
            oyy[i] = oyyi
            ozz[i] = ozzi
            txy[i] = txyi
            tyz[i] = tyzi
            txz[i] = txzi
            max_principal[i] = o1i
            mid_principal[i] = o2i
            min_principal[i] = o3i
            ovm[i] = ovmi
        del solids


        if is_stress:
            word = 'Stress'
            fmt = '%.3f'
        else:
            word = 'Strain'
            fmt = '%.4e'

        # a form is the table of output...
        # Subcase 1         <--- formi  - form_isubcase
        #    Time 1
        #        Stress     <--- form0  - the root level
        #            oxx    <--- formis - form_itime_stress
        #            oyy
        #            ozz

        if dt is None:
            return icase

        header = ''
        if not is_static:
            #print('is_static = %s' % is_static)
            if case is None:
                formis = None
                return icase
            header = self._get_nastran_header(case, dt, itime)
            #form_time[0] = header

        form0 = (word, None, [])
        formis = form0[2]
        # subcase_id, icase, resultType, vector_size, location, dataFormat
        if is_stress and itime == 0:
            if is_element_on.min() == 0:  # if all elements aren't on
                ioff = np.where(is_element_on == 0)[0]
                print('stress_eids_off = %s' % np.array(self.element_ids[ioff]))
                self.log_error('stress_eids_off = %s' % self.element_ids[ioff])
                #stress_res = GuiResult(subcase_id, header='Stress - isElementOn', title='Stress\nisElementOn',
                                    #location='centroid', scalar=oxx, data_format=fmt)
                #cases[icase] = (stress_res, (subcase_id, 'Stress - isElementOn'))

                cases[(1, icase, 'Stress\nisElementOn', 1, 'centroid', '%i', header)] = is_element_on
                form_dict[(key, itime)].append(('Stress - IsElementOn', icase, []))
                icase += 1

        subcase_id = key[2]
        if oxx.min() != oxx.max():
            oxx_res = GuiResult(subcase_id, header=word + 'XX', title=word + 'XX',
                                location='centroid', scalar=oxx, data_format=fmt)
            cases[icase] = (oxx_res, (subcase_id, word + 'XX'))
            form_dict[(key, itime)].append((word + 'XX', icase, []))
            icase += 1
        if oyy.min() != oyy.max():
            oyy_res = GuiResult(subcase_id, header=word + 'YY', title=word + 'YY',
                                location='centroid', scalar=oyy, data_format=fmt)
            cases[icase] = (oyy_res, (subcase_id, word + 'YY'))
            form_dict[(key, itime)].append((word + 'YY', icase, []))
            icase += 1
        if ozz.min() != ozz.max():
            ozz_res = GuiResult(subcase_id, header=word + 'ZZ', title=word + 'ZZ',
                                location='centroid', scalar=ozz, data_format=fmt)
            cases[icase] = (ozz_res, (subcase_id, word + 'ZZ'))
            form_dict[(key, itime)].append((word + 'ZZ', icase, []))
            icase += 1
        if txy.min() != txy.max():
            oxy_res = GuiResult(subcase_id, header=word + 'XY', title=word + 'XY',
                                location='centroid', scalar=txy, data_format=fmt)
            cases[icase] = (oxy_res, (subcase_id, word + 'XY'))
            form_dict[(key, itime)].append((word + 'XY', icase, []))
            icase += 1
        if tyz.min() != tyz.max():
            oyz_res = GuiResult(subcase_id, header=word + 'YZ', title=word + 'YZ',
                                location='centroid', scalar=tyz, data_format=fmt)
            cases[icase] = (oyz_res, (subcase_id, word + 'YZ'))
            form_dict[(key, itime)].append((word + 'YZ', icase, []))
            icase += 1
        if txz.min() != txz.max():
            oxz_res = GuiResult(subcase_id, header=word + 'XZ', title=word + 'XZ',
                                location='centroid', scalar=txz, data_format=fmt)
            cases[icase] = (oxz_res, (subcase_id, word + 'XZ'))
            form_dict[(key, itime)].append((word + 'XZ', icase, []))
            icase += 1
        if max_principal.min() != max_principal.max():
            maxp_res = GuiResult(subcase_id, header='MaxPrincipal', title='MaxPrincipal',
                                 location='centroid', scalar=max_principal, data_format=fmt)
            cases[icase] = (maxp_res, (subcase_id, 'MaxPrincipal'))
            form_dict[(key, itime)].append(('Max Principal', icase, []))
            icase += 1
        if mid_principal.min() != mid_principal.max():
            midp_res = GuiResult(subcase_id, header='MidPrincipal', title='MidPrincipal',
                                 location='centroid', scalar=mid_principal, data_format=fmt)
            cases[icase] = (midp_res, (subcase_id, 'MidPrincipal'))
            form_dict[(key, itime)].append(('Mid Principal', icase, []))
            icase += 1
        if min_principal.min() != min_principal.max():
            minp_res = GuiResult(subcase_id, header='MinPrincipal', title='MinPrincipal',
                                 location='centroid', scalar=min_principal, data_format=fmt)
            cases[icase] = (minp_res, (subcase_id, 'MinPrincipal'))
            form_dict[(key, itime)].append(('Min Principal', icase, []))
            icase += 1
        if vm_word is not None:
            if not is_stress:
                max_min = max(ovm.max(), np.abs(ovm.min()))
                if max_min > 100:
                    raise RuntimeError('vm strain = %s' % ovm)

            ovm_res = GuiResult(subcase_id, header=vm_word, title=vm_word,
                                location='centroid', scalar=ovm, data_format=fmt)
            cases[icase] = (ovm_res, (subcase_id, 'MinPrincipal'))
            form_dict[(key, itime)].append((vm_word, icase, []))
            icase += 1

        #, case, header, form0
        return icase

def get_min_max_theta(faces, all_node_ids, nid_map, xyz_cid0):
    cos_thetas = []
    ideal_theta = []
    #print('faces =', faces)
    #assert len(faces) > 0, 'faces=%s nids=%s' % (faces, all_node_ids)
    for face in faces:
        if len(face) == 3:
            node_ids = all_node_ids[face[0]], all_node_ids[face[1]], all_node_ids[face[2]]
            p1, p2, p3 = nid_map[node_ids[0]], nid_map[node_ids[1]], nid_map[node_ids[2]]
            v21 = xyz_cid0[p2, :] - xyz_cid0[p1, :]
            v32 = xyz_cid0[p3, :] - xyz_cid0[p2, :]
            v13 = xyz_cid0[p1, :] - xyz_cid0[p3, :]

            cos_theta1 = np.dot(v21, -v13) / (np.linalg.norm(v21) * np.linalg.norm(v13))
            cos_theta2 = np.dot(v32, -v21) / (np.linalg.norm(v32) * np.linalg.norm(v21))
            cos_theta3 = np.dot(v13, -v32) / (np.linalg.norm(v13) * np.linalg.norm(v32))
            cos_thetas.extend([cos_theta1, cos_theta2, cos_theta3])
            ideal_theta.extend([piover3, piover3, piover3])
        elif len(face) == 4:
            try:
                node_ids = all_node_ids[face[0]], all_node_ids[face[1]], all_node_ids[face[2]], all_node_ids[face[3]]
            except:
                print(face)
                print(node_ids)
                raise
            p1, p2, p3, p4 = nid_map[node_ids[0]], nid_map[node_ids[1]], nid_map[node_ids[2]], nid_map[node_ids[3]]
            v21 = xyz_cid0[p2, :] - xyz_cid0[p1, :]
            v32 = xyz_cid0[p3, :] - xyz_cid0[p2, :]
            v43 = xyz_cid0[p4, :] - xyz_cid0[p3, :]
            v14 = xyz_cid0[p1, :] - xyz_cid0[p4, :]
            cos_theta1 = np.dot(v21, -v14) / (np.linalg.norm(v21) * np.linalg.norm(v14))
            cos_theta2 = np.dot(v32, -v21) / (np.linalg.norm(v32) * np.linalg.norm(v21))
            cos_theta3 = np.dot(v43, -v32) / (np.linalg.norm(v43) * np.linalg.norm(v32))
            cos_theta4 = np.dot(v14, -v43) / (np.linalg.norm(v14) * np.linalg.norm(v43))
            cos_thetas.extend([cos_theta1, cos_theta2, cos_theta3, cos_theta4])
            ideal_theta.extend([piover2, piover2, piover2, piover2])
        else:
            raise NotImplementedError(face)
    thetas = np.arccos(cos_thetas)
    ideal_theta = np.array(ideal_theta)
    ideal_thetai = max((thetas - ideal_theta).max(), (ideal_theta - thetas).min())

    min_thetai = thetas.min()
    max_thetai = thetas.max()
    return min_thetai, max_thetai, ideal_thetai
