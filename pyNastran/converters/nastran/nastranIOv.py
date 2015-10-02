from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
"""
Defines the GUI IO file for Nastran.
"""
# pylint: disable=C0103,E1101
from __future__ import print_function
from six import iteritems, itervalues
from six.moves import zip, range
from copy import deepcopy

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

import os
from collections import defaultdict, OrderedDict

from numpy import zeros, abs, mean, where, nan_to_num, amax, amin, vstack, array, empty, ones
from numpy import searchsorted, sqrt, pi, arange, unique, allclose, ndarray, int32, cross
from numpy.linalg import norm

import vtk
from vtk import (vtkTriangle, vtkQuad, vtkTetra, vtkWedge, vtkHexahedron,
                 vtkQuadraticTriangle, vtkQuadraticQuad, vtkQuadraticTetra,
                 vtkQuadraticWedge, vtkQuadraticHexahedron,
                 vtkPyramid) #vtkQuadraticPyramid

#from pyNastran import is_release
from pyNastran.bdf.bdf import (BDF,
                               CAERO1, CAERO3, CAERO4, CAERO5, # CAERO2,
                               CQUAD4, CQUAD8, CQUADR, CSHEAR,
                               CTRIA3, CTRIA6, CTRIAR, CTRIAX6,
                               CTETRA4, CTETRA10, CPENTA6, CPENTA15,
                               CHEXA8, CHEXA20,
                               CPYRAM5, CPYRAM13,
                               CONM2,
                               LOAD)
from pyNastran.bdf.cards.elements.shell import ShellElement
from pyNastran.bdf.cards.elements.bars import LineElement
from pyNastran.bdf.cards.elements.springs import SpringElement

from pyNastran.op2.op2 import OP2
from pyNastran.f06.f06_formatting import get_key0
try:
    from pyNastran.op2.op2_geom import OP2Geom
    is_geom = True
except ImportError:
    is_geom = False
#from pyNastran.f06.f06 import F06
is_geom = False

class NastranComplexDisplacementResults(object):
    def __init__(self, subcase_id, titles, xyz, dxyz, scalar,
                 default_scale=40., uname='NastranGeometry'):
        self.subcase_id = subcase_id
        self.data_formats = ['%g'] * len(titles)
        self.xyz = xyz
        self.dxyz = dxyz
        self.dxyz_norm = norm(dxyz, axis=1)
        self.titles = titles
        self.scale = default_scale

        self.default_scale = default_scale
        self.titles_default = deepcopy(titles)
        self.data_formats_default = deepcopy(self.data_formats)
        self.default_scale = default_scales

        #theta = (2*np.pi * i/frame) % (2 * pi)
        theta = 0.0

        # calculate deflections
        eigvs = model.eigenvectors[1000].data[6, :, :]
        defl = scale * (np.real(eigvs[:,:3]) * np.cos(theta) + np.imag(eigvs[:,:3]) * sin(theta))

class NastranDisplacementResults(object):
    def __init__(self, subcase_id, titles, xyz, dxyz, scalar,
                 scales, uname='NastranGeometry'):
        self.subcase_id = subcase_id
        assert self.subcase_id > 0, self.subcase_id
        self.xyz = xyz
        self.dxyz = dxyz
        self.dxyz_norm = norm(dxyz, axis=1)
        self.titles = titles
        self.scales = scales
        self.subcase_id = subcase_id

    def save_defaults(self):
        self.data_formats = ['%g'] * len(self.titles)
        self.titles_default = deepcopy(self.titles)
        self.data_formats_default = deepcopy(self.data_formats)
        self.scales_default = deepcopy(self.scales)

    def get_location(self, i, name):
        return 'node'

    def get_title(self, i, name):
        #j = self.titles_default.index(name)
        #return self.titles[j]
        return self.titles[i]

    def get_default_title(self, i, name):
        #j = self.titles_default.index(name)
        return self.titles_default[i]

    def get_data_format(self, i, name):
        #print(self.titles_default, i)
        #j = self.titles_default.index(name)
        #return self.data_formats[j]
        return self.data_formats[i]

    def get_vector_size(self, i, name):
        #print(i)
        #j = self.titles_default.index(name)
        return 3

    def get_methods(self, i):
        return ['node']

    def get_plot_value(self, i, name):
        return self.dxyz[i, :]
    def get_result(self, i, name):
        return self.dxyz[i, :]

    def get_scalar(self, i, name):
        return self.dxyz_norm

    def get_vector_result(self, i, name):
        xyz = self.xyz + self.scales[i] * self.dxyz[i, :]
        return self.xyz, xyz

    def get_scale(self, i, name):
        #j = self.titles_default.index(name)
        return self.scales[i]

    def get_default_scale(self, i, name):
        #j = self.titles_default.index(name)
        return self.scales_default[i]

    def set_scale(self, i, name, scale):
        j = self.titles_default.index(name)
        self.scales[i] = scale

    def __repr__(self):
        msg = 'NastranDisplacementResults\n'
        msg += '    uname=%r\n' % self.uname
        return msg

class NastranIO(object):
    """
    Defines the GUI class for Nastran.
    """
    def __init__(self):
        #: flips the nastran CAERO subpaneling
        #:   False -> borders of CAEROs can be seen
        #:   True  -> individual subpanels can be seen
        self.show_caero_sub_panels = False

        #: coordinate systems can be messy, so this is the
        #: list of coords to show
        self.show_cids = []
        self.save_data = False
        self.show_caero_actor = True  # show the caero mesh
        self.show_control_surfaces = True
        self.show_conm = True

        self.element_ids = None
        self.node_ids = None
        self.nidMap = None
        self.eidMap = None
        self.nNodes = None
        self.nElements = None
        self.modelType = None
        self.iSubcaseNameMap = None
        self.has_caero = False

    def get_nastran_wildcard_geometry_results_functions(self):
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
        #self.menu_help.menuAction().setVisible(True)
        #self.menu_help2.menuAction().setVisible(False)
        self.nastran_toolbar.setVisible(False)
        self.actions['nastran'].setVisible(False)

    def _create_nastran_tools_and_menu_items(self):
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
            (self.nastran_toolbar, ('caero', 'caero_sub', 'conm'))
            #(self.menu_window, tuple(menu_window)),
            #(self.menu_help, ('load_geometry', 'load_results', 'script', '', 'exit')),
            #(self.menu_help2, ('load_geometry', 'load_results', 'script', '', 'exit')),
        ]
        return tools, menu_items

    def show_caero_mesh(self, is_shown=None):
        """
        :param is_shown: should the mesh be shown/hidden
                         (default=None -> flip between shown/not shown)
        """
        msg = 'self.show_alt_actor=True/False and self.is_sub_panels=True/False may be used'
        self.log.info(msg)
        if is_shown is None:
            is_shown = not self.show_caero_actor

        self.show_caero_actor = is_shown
        if is_shown:
            if not self.show_caero_actor:
                return
            self.geometry_actors['caero'].VisibilityOn()
            if self.show_caero_sub_panels:
                self.geometry_actors['caero_sub'].VisibilityOn()
        else:
            if self.show_caero_actor:
                return
            self.geometry_actors['caero'].VisibilityOff()
            self.geometry_actors['caero_sub'].VisibilityOff()

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
                self.geometry_actors['caero_sub'].VisibilityOn()
            else:
                self.geometry_actors['caero'].VisibilityOn()
        else:
            self.geometry_actors['caero'].VisibilityOff()
            self.geometry_actors['caero_sub'].VisibilityOff()
        self.vtk_interactor.Render()

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
                self.geometry_actors['caero_sub'].VisibilityOn()
            else:
                self.geometry_actors['caero'].VisibilityOn()
                self.geometry_actors['caero_sub'].VisibilityOff()
        self.vtk_interactor.Render()

    def toggle_conms(self):
        """
        Toggle the visibility of the CONMS
        """
        self.show_conm = not self.show_conm
        if 'conm' in self.geometry_actors:
            if self.show_conm:
                self.geometry_actors['conm'].VisibilityOn()
            else:
                self.geometry_actors['conm'].VisibilityOff()
        self.vtk_interactor.Render()

    def _create_coord(self, cid, coord, Type):
        origin = coord.origin
        beta = coord.beta()
        self.create_coordinate_system(label=cid, origin=origin, matrix_3x3=beta, Type=Type)

    def _create_nastran_coords(self, model):
        cid_types = {
            'R' : 'xyz',
            'C' : 'Rtz',
            'S' : 'Rtp',
        }
        for cid, coord in sorted(iteritems(model.coords)):
            if cid == 0:
                continue
            Type = cid_types[coord.Type]
            if self.show_cids is True:
                self._create_coord(cid, coord, Type)
            elif isinstance(self.show_cids, (int, int32)):
                if cid == self.show_cids:
                    self._create_coord(cid, coord, Type)
            elif isinstance(self.show_cids, (list, tuple, ndarray)):
                if cid in self.show_cids:
                    # .. todo:: has issues in VTK 6 I think due to lack of self.grid.Update()
                    self._create_coord(cid, coord, Type)
            else:
                print('skipping cid=%s; use a script and set self.show_cids=[%s] to view' % (cid, cid))

    def _remove_old_nastran_geometry(self, bdf_filename):
        # skip_reading = self.removeOldGeometry(bdf_filename)
        skip_reading = False
        if bdf_filename is None or bdf_filename is '':
            #self.grid = vtk.vtkUnstructuredGrid()
            #self.gridResult = vtk.vtkFloatArray()
            #self.emptyResult = vtk.vtkFloatArray()
            #self.vectorResult = vtk.vtkFloatArray()
            #self.grid2 = vtk.vtkUnstructuredGrid()
            #self.scalarBar.VisibilityOff()
            skip_reading = True
            return skip_reading
        else:
            self.TurnTextOff()
            self.grid.Reset()
            #self.grid2.Reset()

            #self.gridResult = vtk.vtkFloatArray()
            #self.gridResult.Reset()
            #self.gridResult.Modified()
            #self.eidMap = {}
            #self.nidMap = {}

            self.resultCases = {}
            self.nCases = 0
        for i in ('caseKeys', 'iCase', 'iSubcaseNameMap'):
            if hasattr(self, i):  # TODO: is this correct???
                del i
        return skip_reading

    def _remove_old_geometry_old(self, filename, alt_grids):
        skip_reading = False
        if filename is None or filename is '':
            #self.grid = vtk.vtkUnstructuredGrid()
            #self.gridResult = vtk.vtkFloatArray()
            #self.emptyResult = vtk.vtkFloatArray()
            #self.vectorResult = vtk.vtkFloatArray()
            #self.grid2 = vtk.vtkUnstructuredGrid()
            #self.scalarBar.VisibilityOff()
            skip_reading = True
            return skip_reading
        else:
            self.TurnTextOff()
            self.grid.Reset()

            # create alt grids
            yellow = (1., 1., 0.)
            pink = (0.98, 0.4, 0.93)
            if 'caero' not in self.alt_grids:
                self.create_alternate_vtk_grid('caero', color=yellow, line_width=3, opacity=1.0, representation='surface')
            if 'caero_sub' not in self.alt_grids:
                self.create_alternate_vtk_grid('caero_sub', color=yellow, line_width=3, opacity=1.0, representation='surface')
            if 'conm' not in self.alt_grids:
                self.create_alternate_vtk_grid('conm', color=orange, line_width=3, opacity=1.0, point_size=4, representation='point')

            #print('alt_grids', self.alt_grids.keys())

            #self.gridResult = vtk.vtkFloatArray()
            #self.gridResult.Reset()
            #self.gridResult.Modified()
            #self.eidMap = {}
            #self.nidMap = {}

            self.resultCases = {}
            self.nCases = 0
        for i in ('caseKeys', 'iCase', 'iSubcaseNameMap'):
            if hasattr(self, i):  # TODO: is this correct???
                del i
        return skip_reading

    def load_nastran_geometry(self, bdf_filename, dirname, plot=True):
        self.eidMap = {}
        self.nidMap = {}
        #print('bdf_filename=%r' % bdf_filename)
        #key = self.caseKeys[self.iCase]
        #case = self.resultCases[key]

        skip_reading = self._remove_old_nastran_geometry(bdf_filename)
        pink = (0.98, 0.4, 0.93)
        orange = (219/255., 168/255., 13/255.)
        # if 0:
            # yellow = (1., 1., 0.)
            # line_width = 3
            # opacity = 1
            # alt_grids = [
                # ['caero', yellow, line_width, opacity],
                # ['caero_sub', yellow, line_width, opacity],
            # ]
            # skip_reading = self._remove_old_geometry2(bdf_filename, alt_grids=alt_grids)
        if skip_reading:
            return

        if plot:
            self.scalarBar.VisibilityOff()
            self.scalarBar.Modified()

        ext = os.path.splitext(bdf_filename)[0].lower()
        punch = False
        if ext == '.pch':
            punch = True

        xref_loads = True
        if ext == '.op2' and 0 and is_geom:
            model = OP2Geom(make_geom=True, debug=False, log=self.log,
                            debug_file=None)
            model._clear_results()
            model.read_op2(op2_filename=bdf_filename)
            model.cross_reference(xref=True, xref_loads=xref_loads,
                                  xref_constraints=False)
        else:  # read the bdf/punch
            model = BDF(log=self.log, debug=True)
            self.modelType = model.modelType
            model.read_bdf(bdf_filename,
                           punch=punch, xref=False)
            model.cross_reference(xref=True, xref_loads=xref_loads,
                                  xref_constraints=False)

        nnodes = model.nnodes
        assert nnodes > 0
        nelements = model.nelements
        assert nelements > 0

        # count caeros
        ncaeros_sub = 0
        ncaero_sub_points = 0
        for caero in itervalues(model.caeros):
            if hasattr(caero, 'panel_points_elements'):
                npoints, ncelements = caero.get_npanel_points_elements()
                ncaeros_sub += npoints
                ncaero_sub_points += ncelements
            else:
                print('%r doesnt support panel_points_elements' % caero.type)
        ncaeros = model.ncaeros
        ncaeros_points = ncaeros * 4

        box_id_to_caero_element_map = {}
        num_prev = 0
        if model.caeros:
            caero_points = []
            for eid, caero in sorted(iteritems(model.caeros)):
                if caero.type == 'CAERO1':
                    pointsi, elementsi = caero.panel_points_elements()
                    caero_points.append(pointsi)
                    for i, box_id in enumerate(caero.box_ids.flat):
                        box_id_to_caero_element_map[box_id] = elementsi[i, :] + num_prev
                    num_prev += pointsi.shape[0]
            caero_points = vstack(caero_points)
            self.has_caero = True
        else:
            caero_points = empty((0, 3))

        # check for any control surfcaes
        has_control_surface = False
        if model.aesurfs:
            cs_box_ids = []
            has_control_surface = True
            ncaeros_cs = 0
            ncaero_cs_points = 0
            if 'caero_cs' not in self.alt_grids:
                self.create_alternate_vtk_grid('caero_cs', color=pink, line_width=5, opacity=1.0, representation='surface')
            for aid, aesurf in iteritems(model.aesurfs):
                aelist = aesurf.alid1
                ncaeros_cs += len(aelist.elements)
                cs_box_ids.extend(aelist.elements)

                if aesurf.alid2 is not None:
                    aelist = aesurf.alid2
                    ncaeros_cs += len(aelist.elements)
                    cs_box_ids.extend(aelist.elements)

        self.nNodes = nnodes
        self.nElements = nelements  # approximate...

        self.log_info("nNodes=%i nElements=%i" % (self.nNodes, self.nElements))
        msg = model.get_bdf_stats(return_type='list')
        #self.log_info(msg)
        for msgi in msg:
            model.log.debug(msgi)

        if 'CONM2' in model.card_count:
            nCONM2 = model.card_count['CONM2']
        else:
            nCONM2 = 0
        #self.gridResult.SetNumberOfComponents(self.nElements)
        if nCONM2 > 0:
            self.create_alternate_vtk_grid('conm', color=orange, line_width=5, opacity=1., point_size=4, representation='point')

        # Allocate grids
        self.grid.Allocate(self.nElements, 1000)
        if self.has_caero:
            yellow = (1., 1., 0.)
            if 'caero' not in self.alt_grids:
                self.create_alternate_vtk_grid('caero', color=yellow, line_width=3, opacity=1.0, representation='surface')
            if 'caero_sub' not in self.alt_grids:
                self.create_alternate_vtk_grid('caero_sub', color=yellow, line_width=3, opacity=1.0, representation='surface')

            self.alt_grids['caero'].Allocate(ncaeros, 1000)
            self.alt_grids['caero_sub'].Allocate(ncaeros_sub, 1000)
            if has_control_surface:
                self.alt_grids['caero_cs'].Allocate(ncaeros_cs, 1000)

        if nCONM2 > 0:
            self.alt_grids['conm'].Allocate(nCONM2, 1000)

        points = vtk.vtkPoints()
        points.SetNumberOfPoints(self.nNodes)
        #self.gridResult.Allocate(self.nNodes, 1000)
        #vectorReselt.SetNumberOfComponents(3)
        #elem.SetNumberOfPoints(nNodes)


        # add the nodes
        node0 = get_key0(model.nodes)
        position0 = model.nodes[node0].get_position()
        xmin = position0[0]
        xmax = position0[0]

        ymin = position0[1]
        ymax = position0[1]

        zmin = position0[2]
        zmax = position0[2]
        if self.save_data:
            self.model = model

        n = len(model.nodes)
        xyz_cid0 = zeros((n, 3), dtype='float32')
        for i, (nid, node) in enumerate(sorted(iteritems(model.nodes))):
            xyz = node.get_position()
            xyz_cid0[i, :] = xyz
        self.xyz_cid0 = xyz_cid0

        self._create_nastran_coords(model)

        for i, (nid, node) in enumerate(sorted(iteritems(model.nodes))):
            point = node.get_position()
            xmin = min(xmin, point[0])
            xmax = max(xmax, point[0])

            ymin = min(ymin, point[1])
            ymax = max(ymax, point[1])

            zmin = min(zmin, point[2])
            zmax = max(zmax, point[2])
            points.InsertPoint(i, *point)
            self.nidMap[nid] = i

        dim_max = max(xmax-xmin, ymax-ymin, zmax-zmin)
        self.update_axes_length(dim_max)

        self.log_info("xmin=%s xmax=%s dx=%s" % (xmin, xmax, xmax-xmin))
        self.log_info("ymin=%s ymax=%s dy=%s" % (ymin, ymax, ymax-ymin))
        self.log_info("zmin=%s zmax=%s dz=%s" % (zmin, zmax, zmax-zmin))

        j = 0
        nid_to_pid_map, cases, form = self.mapElements(points, self.nidMap, model, j, dim_max, plot=plot, xref_loads=xref_loads)

        if 0:
            nsprings = 0
            if 0:
                for eid, element in sorted(iteritems(model.elements)):
                    if(isinstance(element, LineElement) or
                       isinstance(element, SpringElement) or
                       element.type in ['CBUSH', 'CBUSH1D', 'CFAST', 'CROD', 'CONROD',
                                        'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
                                        'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5', 'CVISC', ]):
                        node_ids = element.node_ids
                        if None in node_ids:
                            nsprings += 1

        # fill grids
        if 'caero' in self.alt_grids:
            self.set_caero_grid(ncaeros_points, model)
            self.set_caero_subpanel_grid(ncaero_sub_points, model)
            if has_control_surface:
                self.set_caero_control_surface_grid(cs_box_ids, box_id_to_caero_element_map, caero_points)

        if nCONM2 > 0:
            self.set_conm_grid(nCONM2, dim_max, model)
        self.set_spc_grid(dim_max, model, nid_to_pid_map)

        # add alternate actors
        self._add_alt_actors(self.alt_grids)

        # set default representation
        if 'caero_cs' in self.geometry_actors:
            self.geometry_properties['caero_cs'].opacity = 0.5

        #for (name, size) in [('conm', 4), ('spc', 4), ('suport', 4)]:
            #if size == 0:
                #continue
            #if name in self.geometry_actors:
                #actor = self.geometry_actors[name]
                #prop = actor.GetProperty()
                #prop.SetRepresentationToPoints()
                #prop.SetPointSize(size)

        # set initial caero visibility
        if 'caero' in self.alt_grids:
            if self.show_caero_actor:
                if self.show_caero_sub_panels:
                    self.geometry_actors['caero'].VisibilityOff()
                    self.geometry_actors['caero_sub'].VisibilityOn()
                else:
                    self.geometry_actors['caero'].VisibilityOn()
                    self.geometry_actors['caero_sub'].VisibilityOff()
            else:
                self.geometry_actors['caero'].VisibilityOff()
                self.geometry_actors['caero_sub'].VisibilityOff()

            if has_control_surface:
                if self.show_control_surfaces:
                    self.geometry_actors['caero_cs'].VisibilityOn()
                else:
                    self.geometry_actors['caero_cs'].VisibilityOn()

            self.geometry_actors['caero'].Modified()
            self.geometry_actors['caero_sub'].Modified()
            if has_control_surface:
                self.geometry_actors['caero_cs'].Modified()
            if hasattr(self.geometry_actors['caero'], 'Update'):
                self.geometry_actors['caero'].Update()
            if hasattr(self.geometry_actors['caero_sub'], 'Update'):
                self.geometry_actors['caero_sub'].Update()
            if has_control_surface and hasattr(self.geometry_actors['caero_sub'], 'Update'):
                    self.geometry_actors['caero_cs'].Update()

        if 'conm' in self.geometry_actors:
            if nCONM2 > 0:
                self.geometry_actors['conm'].VisibilityOn()
            else:
                self.geometry_actors['conm'].VisibilityOff()
            self.geometry_actors['conm'].Modified()

        for name in ['suport', 'spc', 'mpc', 'mpc_dependent', 'mpc_independent']:
            if name in self.geometry_actors:
                self.geometry_actors[name].Modified()

        if plot:
            self.log.info(cases.keys())
            self._finish_results_io2([form], cases)

    def set_caero_grid(self, ncaeros_points, model, j=0):
        """
        Sets the CAERO panel geometry.

        Returns the current id counter.
        """
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(ncaeros_points)

        for eid, element in sorted(iteritems(model.caeros)):
            if(isinstance(element, CAERO1) or isinstance(element, CAERO3) or
               isinstance(element, CAERO4) or isinstance(element, CAERO5)):
                cpoints = element.get_points()
                elem = vtkQuad()
                elem.GetPointIds().SetId(0, j)
                elem.GetPointIds().SetId(1, j + 1)
                elem.GetPointIds().SetId(2, j + 2)
                elem.GetPointIds().SetId(3, j + 3)
                points.InsertPoint(j, *cpoints[0])
                points.InsertPoint(j + 1, *cpoints[1])
                points.InsertPoint(j + 2, *cpoints[2])
                points.InsertPoint(j + 3, *cpoints[3])
                self.alt_grids['caero'].InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                j += 4
            else:
                self.log_info("skipping %s" % element.type)
        self.alt_grids['caero'].SetPoints(points)
        return j

    def set_caero_subpanel_grid(self, ncaero_sub_points, model, j=0):
        """
        Sets the CAERO panel geometry.

        Returns the current id counter.
        """
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(ncaero_sub_points)

        eType = vtkQuad().GetCellType()
        for eid, element in sorted(iteritems(model.caeros)):
            if(isinstance(element, CAERO1) or isinstance(element, CAERO3) or
               isinstance(element, CAERO4) or isinstance(element, CAERO5)):
                pointsi, elementsi = element.panel_points_elements()
                for ipoint, pointii in enumerate(pointsi):
                    points.InsertPoint(j + ipoint, *pointii)

                elem = vtkQuad()
                for elementi in elementsi:
                    elem = vtkQuad()
                    elem.GetPointIds().SetId(0, j + elementi[0])
                    elem.GetPointIds().SetId(1, j + elementi[1])
                    elem.GetPointIds().SetId(2, j + elementi[2])
                    elem.GetPointIds().SetId(3, j + elementi[3])
                    self.alt_grids['caero_sub'].InsertNextCell(eType, elem.GetPointIds())
                j += ipoint + 1
            else:
                self.log_info("skipping %s" % element.type)
        self.alt_grids['caero_sub'].SetPoints(points)
        return j

    def set_caero_control_surface_grid(self, cs_box_ids, box_id_to_caero_element_map, caero_points, j=0):
        points_list = []
        missing_boxes = []
        for ibox, box_id in enumerate(cs_box_ids):
            try:
                ipoints = box_id_to_caero_element_map[box_id]
            except KeyError:
                missing_boxes.append(box_id)
                continue
            points_list.append(caero_points[ipoints, :])
        if missing_boxes:
            msg = 'Missing CAERO AELIST boxes: ' + str(missing_boxes)
            self.log_error(msg)

        points_list = array(points_list)
        ncaero_sub_points = len(unique(points_list.ravel()))

        points = vtk.vtkPoints()
        points.SetNumberOfPoints(ncaero_sub_points)

        eType = vtkQuad().GetCellType()
        for ibox, box_id in enumerate(cs_box_ids):
            try:
                elementi = box_id_to_caero_element_map[box_id]
            except KeyError:
                continue
            pointsi = caero_points[elementi]
            for ipoint, point in enumerate(pointsi):
                # shift z to remove z-fighting with caero in surface representation
                point[2] += 0.001
                points.InsertPoint(j + ipoint, *point)
            elem = vtkQuad()
            elem.GetPointIds().SetId(0, j)
            elem.GetPointIds().SetId(1, j + 1)
            elem.GetPointIds().SetId(2, j + 2)
            elem.GetPointIds().SetId(3, j + 3)
            self.alt_grids['caero_cs'].InsertNextCell(eType, elem.GetPointIds())
            j += ipoint + 1
        self.alt_grids['caero_cs'].SetPoints(points)
        return j

    def set_conm_grid(self, nCONM2, dim_max, model, j=0):
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(nCONM2)

        sphere_size = self._get_sphere_size(dim_max)
        for eid, element in sorted(iteritems(model.masses)):
            if isinstance(element, CONM2):
                xyz = element.nid.get_position()
                c = element.Centroid()
                #d = norm(xyz - c)
                points.InsertPoint(j, *c)

                if 1:
                    elem = vtk.vtkVertex()
                    elem.GetPointIds().SetId(0, j)
                else:
                    elem = vtk.vtkSphere()
                    elem.SetRadius(sphere_size)
                    elem.SetCenter(points.GetPoint(j))

                self.alt_grids['conm'].InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                j += 1
            else:
                self.log_info("skipping %s" % element.type)
        self.alt_grids['conm'].SetPoints(points)
        #self.alt_grids['conm'].Set

    def set_spc_grid(self, dim_max, model, nid_to_pid_map):
        case_control = model.case_control_deck
        for subcase_id, subcase in sorted(iteritems(case_control.subcases)):
            #print('subcase_id=%s' % subcase_id)
            #print(subcase.params.keys())

            if 'SPC' in subcase:
                spc_id, options = subcase.get_parameter('SPC')
                if spc_id is not None:
                    nspcs = model.card_count['SPC'] if 'SPC' in model.card_count else 0
                    nspc1s = model.card_count['SPC1'] if 'SPC1' in model.card_count else 0
                    nspcds = model.card_count['SPCD'] if 'SPCD' in model.card_count else 0
                    if nspcs + nspc1s:
                        self._fill_spc(spc_id, nspcs, nspc1s, nspcds, dim_max, model, nid_to_pid_map)

            if 'MPC' in subcase:
                mpc_id, options = subcase.get_parameter('MPC')
                if spc_id is not None:
                    nmpcs = model.card_count['MPC'] if 'MPC' in model.card_count else 0
                    if nmpcs:
                        self._fill_mpc(mpc_id, dim_max, model, nid_to_pid_map)

            if 'SUPORT1' in subcase.params:  ## TODO: should this be SUPORT?
                # print('suport in subcase %s' % subcase_id)
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
        except:
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
        except:
            model.log.warning('spc_id=%s not found' % spc_id)
            return {}

        node_ids_c1 = defaultdict(str)
        for card in sorted(spcs):
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
        self.create_alternate_vtk_grid('spc', color=purple, line_width=5, opacity=1., point_size=5, representation='point')

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
                    prop = model.properties[pid]
                    if prop.type not in ['PSOLID', 'PLSOLID']:
                        plot_node = True
                if not plot_node:
                    # don't include 456 constraints if they're ONLY on solid elemetns
                    # if we had any bar/plate/etc. elements that use this node, we'll plot the node
                    if not('1' in c1 or '2' in c1 or '3' in c1):
                        continue
            node_ids.append(nid)

        node_ids = unique(node_ids)
        self._add_nastran_nodes_to_grid('spc', node_ids, model, nid_to_pid_map)

    def get_MPCx_node_ids_c1(self, model, mpc_id, exclude_mpcadd=False):
        """
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
        try:
            mpcs = model.mpcs[mpc_id]
        except:
            model.log.warning('mpc_id=%s not found' % mpc_id)
            return []

        # dependent, independent
        lines = []
        for card in sorted(mpcs):
            if card.type == 'MPC':
                nids = card.node_ids
                nid0 = nids[0]
                constraint0 = card.constraints[0]
                enforced0 = card.enforced[0]
                for nid, constraint, enforced in zip(nids[1:], card.constraints[1:], card.enforced[1:]):
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

    def _fill_mpc(self, mpc_id, dim_max, model, nid_to_pid_map):
        green = (0., 1., 0.)
        dunno = (0.5, 1., 0.5)
        self.create_alternate_vtk_grid('mpc_dependent', color=green, line_width=5, opacity=1., point_size=5, representation='point')
        self.create_alternate_vtk_grid('mpc_independent', color=dunno, line_width=5, opacity=1., point_size=5, representation='point')
        self.create_alternate_vtk_grid('mpc_lines', color=dunno, line_width=5, opacity=1., point_size=5, representation='wire')

        lines = self.get_MPCx_node_ids_c1(model, mpc_id, exclude_mpcadd=False)
        lines2 = []
        for line in lines:
            if line not in lines2:
                lines2.append(line)
        lines = array(lines2, dtype='int32')
        dependent = (lines[:, 0])
        independent = unique(lines[:, 1])

        node_ids = unique(lines.ravel())

        self._add_nastran_nodes_to_grid('mpc_dependent', dependent, model, nid_to_pid_map)
        self._add_nastran_nodes_to_grid('mpc_independent', independent, model, nid_to_pid_map)
        self._add_nastran_lines_to_grid('mpc_lines', lines, model, nid_to_pid_map)

    def _add_nastran_nodes_to_grid(self, name, node_ids, model, nid_to_pid_map=None):
        nnodes = len(node_ids)
        if nnodes == 0:
            return
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(nnodes)

        j = 0
        for nid in sorted(node_ids):
            try:
                i = self.nidMap[nid]
            except KeyError:
                model.log.warning('nid=%s does not exist' % nid)
                continue

            # point = self.grid.GetPoint(i)
            # points.InsertPoint(j, *point)

            node = model.nodes[nid]
            point = node.get_position()
            #self.log_info('adding SUPORT1; p=%s' % str(point))
            points.InsertPoint(j, *point)

            if 1:
                elem = vtk.vtkVertex()
                elem.GetPointIds().SetId(0, j)
            else:
                elem = vtk.vtkSphere()
                sphere_size = self._get_sphere_size(dim_max)
                elem.SetRadius(sphere_size)
                elem.SetCenter(points.GetPoint(j))

            self.alt_grids[name].InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            j += 1
        self.alt_grids[name].SetPoints(points)

    def _add_nastran_lines_to_grid(self, name, lines, model, nid_to_pid_map=None):
        nlines = lines.shape[0]
        #nids = unique(lines)
        #nnodes = len(nids)
        nnodes = nlines * 2
        if nnodes == 0:
            return
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(nnodes)

        j = 0
        for nid1, nid2 in lines:
            try:
                i1 = self.nidMap[nid1]
            except KeyError:
                model.log.warning('nid=%s does not exist' % nid1)
                continue
            try:
                i2 = self.nidMap[nid2]
            except KeyError:
                model.log.warning('nid=%s does not exist' % nid2)
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
            # print(nid1, nid2)
            self.alt_grids[name].InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            j += 2
        self.alt_grids[name].SetPoints(points)

    def set_quad_grid(self, name, nodes, elements, color, line_width=1, opacity=1.):
        """
        Makes a CQUAD4 grid
        """
        self.create_alternate_vtk_grid(name, color=color, line_width=5, opacity=1., representation='wire')

        nnodes = nodes.shape[0]
        nquads = elements.shape[0]
        #print(nodes)
        if nnodes == 0:
            return
        if nquads == 0:
            return

        print('adding quad_grid %s; nnodes=%s nquads=%s' % (name, nnodes, nquads))
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(nnodes)
        for nid, node in enumerate(nodes):
            #print(nid, node)
            points.InsertPoint(nid, *list(node))

        #elem = vtkQuad()
        #assert elem.GetCellType() == 9, elem.GetCellType()
        self.alt_grids[name].Allocate(nquads, 1000)
        for element in elements:
            elem = vtkQuad()
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, element[0])
            point_ids.SetId(1, element[1])
            point_ids.SetId(2, element[2])
            point_ids.SetId(3, element[3])
            self.alt_grids[name].InsertNextCell(9, elem.GetPointIds())
        self.alt_grids[name].SetPoints(points)

        self._add_alt_actors({name : self.alt_grids[name]})

        #if name in self.geometry_actors:
        self.geometry_actors[name].Modified()

    def _fill_suport(self, suport_id, dim_max, model):
        #pink = (0.98, 0.4, 0.93)
        red = (1.0, 0., 0.)
        self.create_alternate_vtk_grid('suport', color=red, line_width=5, opacity=1., point_size=4, representation='point')

        node_ids = []

        # list
        #for suport in model.suport:
            #node_ids += suport.IDs

        # dict
        suport1 = model.suport1[suport_id]
        node_ids += suport1.IDs

        node_ids = unique(node_ids)
        self._add_nastran_nodes_to_grid('suport', node_ids, model)

    def _get_sphere_size(self, dim_max):
        return 0.01 * dim_max

    def mapElements(self, points, nidMap, model, j, dim_max,
                    plot=True, xref_loads=True):
        sphere_size = self._get_sphere_size(dim_max)
        #self.eidMap = {}

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
        pids = zeros(nelements, 'int32')

        # pids_good = []
        # pids_to_keep = []
        # pids_btm = []
        # pids_to_drop = []
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
            self.eidMap[eid] = i
            pid = 0
            if isinstance(element, CTRIA3) or isinstance(element, CTRIAR):
                elem = vtkTriangle()
                nodeIDs = element.node_ids
                pid = element.Pid()
                self.eid_to_nid_map[eid] = nodeIDs
                for nid in nodeIDs:
                    if nid is not None:
                        nid_to_pid_map[nid].append(pid)

                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element, CTRIA6):
                nodeIDs = element.node_ids
                pid = element.Pid()
                self.eid_to_nid_map[eid] = nodeIDs[:3]
                for nid in nodeIDs:
                    if nid is not None:
                        nid_to_pid_map[nid].append(pid)
                if None not in nodeIDs:
                    elem = vtkQuadraticTriangle()
                    elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                    elem.GetPointIds().SetId(4, nidMap[nodeIDs[4]])
                    elem.GetPointIds().SetId(5, nidMap[nodeIDs[5]])
                else:
                    elem = vtkTriangle()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element, CTRIAX6):
                # midside nodes are required, nodes out of order
                nodeIDs = element.node_ids
                pid = element.Pid()
                for nid in nodeIDs:
                    if nid is not None:
                        nid_to_pid_map[nid].append(pid)

                if None not in nodeIDs:
                    elem = vtkQuadraticTriangle()
                    elem.GetPointIds().SetId(3, nidMap[nodeIDs[1]])
                    elem.GetPointIds().SetId(4, nidMap[nodeIDs[3]])
                    elem.GetPointIds().SetId(5, nidMap[nodeIDs[5]])
                else:
                    elem = vtkTriangle()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[4]])
                self.eid_to_nid_map[eid] = [nodeIDs[0], nodeIDs[2], nodeIDs[4]]
                #a = [0, 2, 4]
                #msg = "CTRIAX6 %i %i %i" %(nidMap[nodeIDs[a[0]]],
                #                           nidMap[nodeIDs[a[1]]],
                #                           nidMap[nodeIDs[a[2]]])
                #raise RuntimeError(msg)
                #sys.stdout.flush()

                #elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                #elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                #elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())

            elif (isinstance(element, CQUAD4) or isinstance(element, CSHEAR) or
                  isinstance(element, CQUADR)):
                nodeIDs = element.node_ids
                pid = element.Pid()
                for nid in nodeIDs:
                    if nid is not None:
                        nid_to_pid_map[nid].append(pid)
                self.eid_to_nid_map[eid] = nodeIDs[:4]
                elem = vtkQuad()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element, CQUAD8):
                nodeIDs = element.node_ids
                pid = element.Pid()
                for nid in nodeIDs:
                    if nid is not None:
                        nid_to_pid_map[nid].append(pid)
                self.eid_to_nid_map[eid] = nodeIDs[:4]
                if None not in nodeIDs:
                    elem = vtkQuadraticQuad()
                    elem.GetPointIds().SetId(4, nidMap[nodeIDs[4]])
                    elem.GetPointIds().SetId(5, nidMap[nodeIDs[5]])
                    elem.GetPointIds().SetId(6, nidMap[nodeIDs[6]])
                    elem.GetPointIds().SetId(7, nidMap[nodeIDs[7]])
                else:
                    elem = vtkQuad()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element, CTETRA4):
                elem = vtkTetra()
                nodeIDs = element.node_ids
                pid = element.Pid()
                for nid in nodeIDs:
                    nid_to_pid_map[nid].append(pid)
                self.eid_to_nid_map[eid] = nodeIDs[:4]
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element, CTETRA10):
                nodeIDs = element.node_ids
                pid = element.Pid()
                for nid in nodeIDs:
                    if nid is not None:
                        nid_to_pid_map[nid].append(pid)
                self.eid_to_nid_map[eid] = nodeIDs[:4]
                if None not in nodeIDs:
                    elem = vtkQuadraticTetra()
                    elem.GetPointIds().SetId(4, nidMap[nodeIDs[4]])
                    elem.GetPointIds().SetId(5, nidMap[nodeIDs[5]])
                    elem.GetPointIds().SetId(6, nidMap[nodeIDs[6]])
                    elem.GetPointIds().SetId(7, nidMap[nodeIDs[7]])
                    elem.GetPointIds().SetId(8, nidMap[nodeIDs[8]])
                    elem.GetPointIds().SetId(9, nidMap[nodeIDs[9]])
                else:
                    elem = vtkTetra()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element, CPENTA6):
                elem = vtkWedge()
                nodeIDs = element.node_ids
                pid = element.Pid()
                for nid in nodeIDs:
                    nid_to_pid_map[nid].append(pid)
                self.eid_to_nid_map[eid] = nodeIDs[:6]
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                elem.GetPointIds().SetId(4, nidMap[nodeIDs[4]])
                elem.GetPointIds().SetId(5, nidMap[nodeIDs[5]])
                self.grid.InsertNextCell(elem.GetCellType(),
                                         elem.GetPointIds())

            elif isinstance(element, CPENTA15):
                nodeIDs = element.node_ids
                pid = element.Pid()
                for nid in nodeIDs:
                    if nid is not None:
                        nid_to_pid_map[nid].append(pid)
                self.eid_to_nid_map[eid] = nodeIDs[:6]
                if None not in nodeIDs:
                    elem = vtkQuadraticWedge()
                    elem.GetPointIds().SetId(6, nidMap[nodeIDs[6]])
                    elem.GetPointIds().SetId(7, nidMap[nodeIDs[7]])
                    elem.GetPointIds().SetId(8, nidMap[nodeIDs[8]])
                    elem.GetPointIds().SetId(9, nidMap[nodeIDs[9]])
                    elem.GetPointIds().SetId(10, nidMap[nodeIDs[10]])
                    elem.GetPointIds().SetId(11, nidMap[nodeIDs[11]])
                    elem.GetPointIds().SetId(12, nidMap[nodeIDs[12]])
                    elem.GetPointIds().SetId(13, nidMap[nodeIDs[13]])
                    elem.GetPointIds().SetId(14, nidMap[nodeIDs[14]])
                else:
                    elem = vtkWedge()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                elem.GetPointIds().SetId(4, nidMap[nodeIDs[4]])
                elem.GetPointIds().SetId(5, nidMap[nodeIDs[5]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element, CHEXA8):
                nodeIDs = element.node_ids
                pid = element.Pid()
                for nid in nodeIDs:
                    nid_to_pid_map[nid].append(pid)
                self.eid_to_nid_map[eid] = nodeIDs[:8]
                elem = vtkHexahedron()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                elem.GetPointIds().SetId(4, nidMap[nodeIDs[4]])
                elem.GetPointIds().SetId(5, nidMap[nodeIDs[5]])
                elem.GetPointIds().SetId(6, nidMap[nodeIDs[6]])
                elem.GetPointIds().SetId(7, nidMap[nodeIDs[7]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element, CHEXA20):
                nodeIDs = element.node_ids
                pid = element.Pid()
                for nid in nodeIDs:
                    if nid is not None:
                        nid_to_pid_map[nid].append(pid)
                if None not in nodeIDs:
                    elem = vtkQuadraticHexahedron()
                    elem.GetPointIds().SetId(8, nidMap[nodeIDs[8]])
                    elem.GetPointIds().SetId(9, nidMap[nodeIDs[9]])
                    elem.GetPointIds().SetId(10, nidMap[nodeIDs[10]])
                    elem.GetPointIds().SetId(11, nidMap[nodeIDs[11]])
                    elem.GetPointIds().SetId(12, nidMap[nodeIDs[12]])
                    elem.GetPointIds().SetId(13, nidMap[nodeIDs[13]])
                    elem.GetPointIds().SetId(14, nidMap[nodeIDs[14]])
                    elem.GetPointIds().SetId(15, nidMap[nodeIDs[15]])
                    elem.GetPointIds().SetId(16, nidMap[nodeIDs[16]])
                    elem.GetPointIds().SetId(17, nidMap[nodeIDs[17]])
                    elem.GetPointIds().SetId(18, nidMap[nodeIDs[18]])
                    elem.GetPointIds().SetId(19, nidMap[nodeIDs[19]])
                else:
                    elem = vtkHexahedron()

                self.eid_to_nid_map[eid] = nodeIDs[:8]
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                elem.GetPointIds().SetId(4, nidMap[nodeIDs[4]])
                elem.GetPointIds().SetId(5, nidMap[nodeIDs[5]])
                elem.GetPointIds().SetId(6, nidMap[nodeIDs[6]])
                elem.GetPointIds().SetId(7, nidMap[nodeIDs[7]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())

            elif isinstance(element, CPYRAM5):
                nodeIDs = element.node_ids
                pid = element.Pid()
                for nid in nodeIDs:
                    nid_to_pid_map[nid].append(pid)
                self.eid_to_nid_map[eid] = nodeIDs[:5]
                elem = vtkPyramid()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                elem.GetPointIds().SetId(4, nidMap[nodeIDs[4]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element, CPYRAM13):
                nodeIDs = element.node_ids
                pid = element.Pid()
                #if None not in nodeIDs:
                    #print(' node_ids =', nodeIDs)
                    #elem = vtkQuadraticPyramid()
                    #elem.GetPointIds().SetId(5, nidMap[nodeIDs[5]])
                    #elem.GetPointIds().SetId(6, nidMap[nodeIDs[6]])
                    #elem.GetPointIds().SetId(7, nidMap[nodeIDs[7]])
                    #elem.GetPointIds().SetId(8, nidMap[nodeIDs[8]])
                    #elem.GetPointIds().SetId(9, nidMap[nodeIDs[9]])
                    #elem.GetPointIds().SetId(10, nidMap[nodeIDs[10]])
                    #elem.GetPointIds().SetId(11, nidMap[nodeIDs[11]])
                    #elem.GetPointIds().SetId(12, nidMap[nodeIDs[12]])
                #else:
                elem = vtkPyramid()
                #print('*node_ids =', nodeIDs[:5])

                self.eid_to_nid_map[eid] = nodeIDs[:5]

                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                elem.GetPointIds().SetId(4, nidMap[nodeIDs[4]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())

            elif (isinstance(element, LineElement) or
                  isinstance(element, SpringElement) or
                  element.type in ['CBUSH', 'CBUSH1D', 'CFAST', 'CROD', 'CONROD',
                                   'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
                                   'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5', 'CVISC', ]):

                # TODO: verify
                # CBUSH, CBUSH1D, CFAST, CROD, CELAS1, CELAS3
                # CDAMP1, CDAMP2, CDAMP3, CDAMP4, CDAMP5, CVISC
                if hasattr(element, 'pid'):
                    pid = element.Pid()
                else:
                    # CONROD
                    # CELAS2, CELAS4?
                    pid = 0
                nodeIDs = element.node_ids
                for nid in nodeIDs:
                    if nid is not None:
                        nid_to_pid_map[nid].append(pid)

                if nodeIDs[0] is None and  nodeIDs[0] is None: # CELAS2
                    print('removing CELASx eid=%i -> no node %i' % (eid, nodeIDs[0]))
                    del self.eidMap[eid]
                    continue
                if None in nodeIDs:  # used to be 0...
                    if nodeIDs[0] is None:
                        slot = 1
                    elif nodeIDs[1] is None:
                        slot = 0
                    #print('nodeIDs=%s slot=%s' % (str(nodeIDs), slot))
                    self.eid_to_nid_map[eid] = nodeIDs[slot]
                    nid = nodeIDs[slot]
                    if nid not in nidMap:
                        # SPOINT
                        print('removing CELASx eid=%i -> SPOINT %i' % (eid, nid))
                        continue

                    #c = nidMap[nid]
                    elem = vtk.vtkVertex()
                    elem.GetPointIds().SetId(0, j)

                    elem = vtk.vtkSphere()
                    #if d == 0.:
                        #d = sphere_size
                    elem.SetRadius(sphere_size)
                else:
                    # 2 points
                    #d = norm(element.nodes[0].get_position() - element.nodes[1].get_position())
                    self.eid_to_nid_map[eid] = nodeIDs
                    elem = vtk.vtkLine()
                    try:
                        elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                        elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                    except KeyError:
                        print("nodeIDs =", nodeIDs)
                        print(str(element))
                        continue

                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            else:
                print('removing eid=%s' % eid)
                del self.eidMap[eid]
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
            i += 1
        assert len(self.eidMap) > 0, self.eidMap

        nelements = i
        self.nElements = nelements
        #print('nelements=%s pids=%s' % (nelements, list(pids)))
        pids = pids[:nelements]
        #print('len(pids) = ', len(pids))
        self.grid.SetPoints(points)

        self.grid.Modified()
        if hasattr(self.grid, 'Update'):
            self.grid.Update()
        #self.log_info("updated grid")

        cases = {}
        #pids = array(pids, 'int32')
        #print('eid_map')
        #for key, value in sorted(iteritems(self.eidMap)):
            #print('  %s %s' % (key, value))

        if 0:
            if not len(pids) == len(self.eidMap):
                msg = 'ERROR:  len(pids)=%s len(eidMap)=%s\n' % (len(pids), len(self.eidMap))
                for eid, pid in sorted(iteritems(pids_dict)):
                    #self.eidMap[eid] = i
                    #pids_dict[eid] = pid
                    if eid not in self.eidMap:
                        msg += 'eid=%s %s' % (eid, str(model.elements[eid]))
                raise RuntimeError(msg)
        del pids_dict


        self.iSubcaseNameMap = {1: ['Nastran', '']}
        #nelements = len(self.eidMap)
        icase = 0
        form = ['Geometry', None, []]
        form0 = form[2]

        # set to True to enable nodeIDs as an result
        nidsSet = True
        if nidsSet:
            nids = zeros(self.nNodes, dtype='int32')
            for (nid, nid2) in iteritems(self.nidMap):
                nids[nid2] = nid
            cases[(0, icase, 'NodeID', 1, 'node', '%i')] = nids
            form0.append(('NodeID', icase, []))
            icase += 1
            self.node_ids = nids
            nidsSet = True

        # set to True to enable elementIDs as a result
        eidsSet = True
        if eidsSet:
            eids = zeros(nelements, dtype='int32')
            for (eid, eid2) in iteritems(self.eidMap):
                eids[eid2] = eid
            cases[(0, icase, 'ElementID', 1, 'centroid', '%i')] = eids
            form0.append(('ElementID', icase, []))
            icase += 1
            self.element_ids = eids
            eidsSet = True

        # subcase_id, resultType, vector_size, location, dataFormat
        if len(model.properties):
            cases[(0, icase, 'PropertyID', 1, 'centroid', '%i')] = pids
            form0.append(('PropertyID', icase, []))
            icase += 1

        icase = self._plot_pressures(model, cases, form0, icase, xref_loads)
        icase = self._plot_applied_loads(model, cases, form0, icase)

        if 0:
            nxs = []
            nys = []
            nzs = []
            i = 0

            for eid, element in sorted(iteritems(model.elements)):
                if isinstance(element, ShellElement):
                    (nx, ny, nz) = element.Normal()
                else:
                    nx = ny = nz = 0.0
                nxs.append(nx)
                nys.append(ny)
                nzs.append(nz)

            # if not a flat plate
            #if min(nxs) == max(nxs) and min(nxs) != 0.0:
            # subcase_id, resultType, vector_size, location, dataFormat
            cases[(0, icase, 'Normalx', 1, 'centroid', '%.1f')] = nxs
            form0.append(('Normalx', icase, []))
            icase += 1

            cases[(0, icase, 'Normaly', 1, 'centroid', '%.1f')] = nys
            form0.append(('Normaly', icase, []))
            icase += 1

            cases[(0, icase, 'Normalz', 1, 'centroid', '%.1f')] = nzs
            form0.append(('Normalz', icase, []))
            icase += 1

        return nid_to_pid_map, cases, form

    def _plot_pressures(self, model, cases, form0, icase, xref_loads):
        """
        pressure act normal to the face (as opposed to anti-normal)
        """
        assert xref_loads is True, 'xref_loads must be set to True; change it above near the read_bdf'
        try:
            sucaseIDs = model.case_control_deck.get_subcase_list()
        except AttributeError:
            return icase

        print('_plot_pressures')
        for subcase_id in sucaseIDs:
            if subcase_id == 0:
                continue
            try:
                load_case_id, options = model.case_control_deck.get_subcase_parameter(subcase_id, 'LOAD')
            except KeyError:
                continue
            try:
                loadCase = model.loads[load_case_id]
            except KeyError:
                self.log.warning('LOAD=%s not found' % load_case_id)
                continue

            # account for scale factors
            loads2 = []
            scale_factors2 = []
            for load in loadCase:
                if isinstance(load, LOAD):
                    scale_factors, loads = load.get_reduced_loads()
                    scale_factors2 += scale_factors
                    loads2 += loads
                else:
                    scale_factors2.append(1.)
                    loads2.append(load)

            eids = sorted(model.elements.keys())
            pressures = zeros(len(model.elements), dtype='float32')

            iload = 0
            # loop thru scaled loads and plot the pressure
            for load, scale in zip(loads2, scale_factors2):
                if iload % 5000 == 0:
                    print('NastranIOv iload=%s' % iload)
                if load.type == 'PLOAD4':
                    elem = load.eid
                    if elem.type in ['CTRIA3', 'CTRIA6', 'CTRIA', 'CTRIAR',
                                     'CQUAD4', 'CQUAD8', 'CQUAD', 'CQUADR', 'CSHEAR']:
                        p = load.pressures[0] * scale

                        # single element per PLOAD
                        #eid = elem.eid
                        #pressures[eids.index(eid)] = p

                        # multiple elements
                        for el in load.eids:
                            pressures[eids.index(el.eid)] += p
                    #elif elem.type in ['CTETRA', 'CHEXA', 'CPENTA']:
                        #A, centroid, normal = elem.getFaceAreaCentroidNormal(load.g34.nid, load.g1.nid)
                        #r = centroid - p
                iload += 1
            # if there is no applied pressure, don't make a plot
            if abs(pressures).max():
                case_name = 'Pressure Case=%i' % subcase_id
                # print('iload=%s' % iload)
                # print(case_name)
                # subcase_id, resultType, vector_size, location, dataFormat
                cases[(0, case_name, 1, 'centroid', '%.1f')] = pressures
                form0.append((case_name, icase, []))
                icase += 1
        return icase

    def _plot_applied_loads(self, model, cases, form0, icase):
        print('_plot_applied_loads')
        try:
            sucase_ids = model.case_control_deck.get_subcase_list()
        except AttributeError:
            print('no subcases....')
            return icase

        for subcase_id in sucase_ids:
            if subcase_id == 0:
                continue

            found_load = False
            found_temperature = False

            form = []
            for key in ('LOAD', 'TEMPERATURE(INITIAL)'):
                try:
                    load_case_id, options = model.case_control_deck.get_subcase_parameter(subcase_id, key)
                except KeyError:
                    print('no load for isubcase=%s' % subcase_id)
                    continue
                try:
                    load_case = model.loads[load_case_id]
                except KeyError:
                    self.log.warning('LOAD=%s not found' % load_case_id)
                    continue

                if key == 'LOAD':
                    p0 = array([0., 0., 0.], dtype='float32')
                    pressures, forces, spcd = self._get_forces_moments_array(model, p0, load_case_id, include_grav=False)
                    found_load = True
                elif key == 'TEMPERATURE(INITIAL)':
                    temperatures = self._get_temperatures_array(model, load_case_id)
                    found_temperature = True
                else:
                    raise NotImplementedError(key)

            if found_load:
                if abs(pressures).max():
                    case_name = 'Pressure Case=%i' % subcase_id
                    # print('iload=%s' % iload)
                    # print(case_name)
                    # subcase_id, resultType, vector_size, location, dataFormat
                    cases[(0, case_name, 1, 'centroid', '%.1f')] = pressures
                    form.append((case_name, icase, []))
                    icase += 1

                if abs(forces.max() - forces.min()) > 0.0:
                    print('plot applied loads loadcase =', load_case_id)
                    # if forces[:, 0].min() != forces[:, 0].max():
                    cases[(subcase_id, icase, 'LoadX Case=%i' % subcase_id, 1, 'node', '%.1f')] = forces[:, 0]
                    # if forces[:, 1].min() != forces[:, 1].max():
                    cases[(subcase_id, icase + 1, 'LoadY Case=%i' % subcase_id, 1, 'node', '%.1f')] = forces[:, 1]
                    # if forces[:, 2].min() != forces[:, 2].max():
                    cases[(subcase_id, icase + 2, 'LoadZ Case=%i' % subcase_id, 1, 'node', '%.1f')] = forces[:, 2]

                    form.append(('Total Load FX', icase, []))
                    form.append(('Total Load FY', icase + 1, []))
                    form.append(('Total Load FZ', icase + 2, []))
                    icase += 2

                if abs(spcd).max():
                    cases[(subcase_id, icase, 'SPCDx', 1, 'node', '%.3g')] = spcd[:, 0]
                    form.append(('SPCDx', icase, []))
                    icase += 1

                    cases[(subcase_id, icase, 'SPCDy', 1, 'node', '%.3g')] = spcd[:, 1]
                    form.append(('SPCDy', icase, []))
                    icase += 1

                    #cases[(subcase_id, icase, name + 'Z', 1, 'node', '%g', header)] = t3
                    cases[(subcase_id, icase, 'SPCDz', 1, 'node', '%.3g')] = spcd[:, 2]
                    form.append(('SPCDz', icase, []))
                    icase += 1

                    t123 = spcd[:, :3]
                    tnorm = norm(t123, axis=1)
                    assert len(tnorm) == len(spcd[:, 2]), len(spcd[:, 2])
                    cases[(subcase_id, icase, 'SPCD XYZ', 1, 'node', '%.3g')] = tnorm
                    form.append(('SPCD XYZ', icase, []))
                    icase += 1
            if found_temperature:
                cases[(subcase_id, icase, 'Temperature(Initial)', 1, 'node', '%.3g')] = temperatures
                form.append(('Temperature(Initial)', icase, []))
                icase += 1
            if form:
                form0.append(('Load Case=%i' % subcase_id, None, form))
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
        print('_get_forces_moments_array')
        nids = sorted(model.nodes.keys())
        nnodes = len(nids)

        load_case = model.loads[load_case_id]
        loads2, scale_factors2 = self._get_loads_and_scale_factors(load_case)
        tempd = model.tempds[load_case_id].temperature if load_case_id in model.tempds else 0.
        temperatures = ones(len(model.nodes), dtype='float32') * tempd
        for load, scale in zip(loads2, scale_factors2):
            if load.type == 'TEMP':
                #print(dir(load))
                temps_dict = load.temperatures
                for nid, val in iteritems(temps_dict):
                    nidi = nids.index(nid)
                    temperatures[nidi] = val
            else:
                print(load.type)
        return temperatures

    def _get_forces_moments_array(self, model, p0, load_case_id, include_grav=False):
        print('_get_forces_moments_array')
        nids = sorted(model.nodes.keys())
        nnodes = len(nids)

        load_case = model.loads[load_case_id]
        loads2, scale_factors2 = self._get_loads_and_scale_factors(load_case)

        eids = sorted(model.elements.keys())
        pressures = zeros(len(model.elements), dtype='float32')

        forces = zeros((nnodes, 3), dtype='float32')
        spcd = zeros((nnodes, 3), dtype='float32')
        # loop thru scaled loads and plot the pressure
        cards_ignored = {}
        for load, scale in zip(loads2, scale_factors2):
            if load.type == 'FORCE':
                scale2 = load.mag * scale  # does this need a magnitude?
                nid = load.node
                forces[nids.index(nid)] += load.xyz * scale2

            elif load.type == 'PLOAD2':
                pressure = load.pressures[0] * scale  # there are 4 pressures, but we assume p0
                for eid in load.eids:
                    elem = self.elements[eid]
                    if elem.type in ['CTRIA3',
                                     'CQUAD4', 'CSHEAR']:
                        node_ids = elem.node_ids
                        nnodes = len(node_ids)
                        n = elem.Normal()
                        area = elem.Area()
                        f = pressure * n * area / nnodes
                        # r = elem.Centroid() - p0
                        # m = cross(r, f)
                        for nid in node_ids:
                            forces[nids.index(nid)] += f
                        forces += f
                        # F += f
                        # M += m
                    else:
                        self.log.debug('case=%s etype=%r loadtype=%r not supported' % (load_case_id, elem.type, load.type))

            elif load.type == 'PLOAD4':
                continue  ## TODO: should be removed
                # elem = load.eid
                # A = elem.get_area()
                if 0:
                    if elem.type in ['CTRIA3', 'CTRIA6', 'CTRIA', 'CTRIAR',]:
                        eid = elem.eid
                        node_ids = elem.node_ids
                        k = load.pressures[0] * scale / 3.
                        # TODO: doesn't consider load.eids for distributed pressures???
                        for nid in node_ids[3:]:
                            pressures[eids.index(nid)] += k
                    elif elem.type in ['CQUAD4', 'CQUAD8', 'CQUAD', 'CQUADR', 'CSHEAR']:
                        eid = elem.eid
                        node_ids = elem.node_ids
                        k = load.pressures[0] * scale / 4.
                        # TODO: doesn't consider load.eids for distributed pressures???
                        for nid in node_ids[4:]:
                            pressures[eids.index(nid)] += k

                else:
                    # single element per PLOAD
                    #eid = elem.eid
                    #pressures[eids.index(eid)] = p

                    p = load.pressures[0] * scale

                    # multiple elements
                    for el in load.eids:
                        # pressures[eids.index(el.eid)] += p
                        A = el.get_area()
                        F = p * A * scale
                        for nid in el.node_ids:
                            forces[nids.index(nid)] += F
                #elif elem.type in ['CTETRA', 'CHEXA', 'CPENTA']:
            elif load.type == 'SPCD':
                #self.gids = [integer(card, 2, 'G1'),]
                #self.constraints = [components_or_blank(card, 3, 'C1', 0)]
                #self.enforced = [double_or_blank(card, 4, 'D1', 0.0)]
                for nid, c1, d1 in zip(load.node_ids, load.constraints, load.enforced):
                    c1 = int(c1)
                    assert c1 in [1, 2, 3, 4, 5, 6], c1
                    if c1 < 4:
                        spcd[nids.index(nid), c1 - 1] = d1
            else:
                if load.type not in cards_ignored:
                    cards_ignored[load.type] = True
                    print('  _get_forces_moment_array - load.type = %s' % load.type)

        return pressures, forces, spcd

    def load_nastran_results(self, op2_filename, dirname):
        """
        Loads the Nastran results into the GUI
        """
        #gridResult.SetNumberOfComponents(self.nElements)
        self.TurnTextOn()
        self.scalarBar.VisibilityOn()
        self.scalarBar.Modified()
        #self.show_caero_mesh()

        print("tring to read...%s" % op2_filename)
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
            model.read_op2(op2_filename, combine=True)

            self.log.info(model.get_op2_stats())
            print(model.get_op2_stats())

        elif ext == '.pch':
            raise NotImplementedError('*.pch is not implemented; filename=%r' % op2_filename)
        #elif ext == '.f06':
            #model = F06(log=self.log, debug=True)
            #model.set_vectorization(True)
            #model.read_f06(op2_filename)
        else:
            print("error...")
            msg = 'extension=%r is not supported; filename=%r' % (ext, op2_filename)
            raise NotImplementedError(msg)

        #print(model.print_results())
        #self.iSubcaseNameMap[self.isubcase] = [Subtitle, Label]

        cases = OrderedDict()
        subcase_ids = model.iSubcaseNameMap.keys()
        #self.iSubcaseNameMap = model.iSubcaseNameMap
        # self.iSubcaseNameMap = model.subcase_key
        #print(self.iSubcaseNameMap)
        self.iSubcaseNameMap = {}
        for isubcase, values in iteritems(model.iSubcaseNameMap):
            if not isinstance(isubcase, (int, int32)):
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

        form = []
        icase = 0
        for subcase_id in subcase_ids:
            if subcase_id == 0:
                # subcase id can be 0...what...see ISAT....
                continue
            subcase_name = 'Subcase %i' % subcase_id
            form0 = (subcase_name, None, [])
            formi = form0[2]
            icase = self.fill_oug_oqg(cases, model, subcase_id, formi, icase)
            icase = self.fill_stress(cases, model, subcase_id, formi, icase)
            if len(formi):
                form.append(form0)

        self._finish_results_io2(form, cases)

    def __fill_cart3d_case2(self, cases, ID, nodes, elements, regions, model):
        print('_fill_cart3d_case2')
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
        labels = ['NodeID', 'ElementID', 'Region',
                  'Normal X', 'Normal Y', 'Normal Z']
        cart3d_geo = Cart3dGeometry(subcase_id, labels,
                                    nids, eids, regions, cnormals,
                                    uname='Cart3dGeometry')

        cases = {
            0 : (cart3d_geo, (0, 'NodeID')),
            1 : (cart3d_geo, (0, 'ElementID')),
            2 : (cart3d_geo, (0, 'Region')),
            3 : (cart3d_geo, (0, 'NormalX')),
            4 : (cart3d_geo, (0, 'NormalY')),
            5 : (cart3d_geo, (0, 'NormalZ')),
        }
        geometry_form = [
            ('NodeID', 0, []),
            ('ElementID', 1, []),
            ('Region', 2, []),
            ('Normal X', 3, []),
            ('Normal Y', 4, []),
            ('Normal Z', 5, []),
        ]
        form = [
            ('Geometry', None, geometry_form),
        ]
        icase = 6
        return form, cases, icase

    def fill_oug_oqg(self, cases, model, subcase_id, formi, icase):
        """
        loads the nodal dispalcements/velocity/acceleration/eigenvector/spc/mpc forces
        """
        nnodes = self.nNodes
        displacement_like = [
            (model.displacements, 'Displacement'),
            (model.velocities, 'Velocity'),
            (model.accelerations, 'Acceleration'),
            (model.eigenvectors, 'Eigenvectors'),
            (model.spc_forces, 'SPC Forces'),
            (model.mpc_forces, 'MPC Forces'),

            # untested
            (model.load_vectors, 'LoadVectors'),
            (model.applied_loads, 'AppliedLoads'),
            (model.force_vectors, 'ForceVectors'),
            #[model.grid_point_forces, 'GridPointForces'],  # TODO: this is buggy...
        ]
        temperature_like = [
            (model.temperatures, 'Temperature'),
        ]
        nids = self.node_ids

        keys2 = []
        for (result, name) in displacement_like:
            keys = result.keys()
            if len(keys) == 0:
                continue

            for key in keys:
                if isinstance(key, (int, int32)):
                    if key != subcase_id:
                        continue
                else:
                    subcase_idi = key[0]
                    if subcase_id != subcase_idi:
                        continue

                if key not in result:
                    continue
                if key not in keys2:
                    keys2.append(key)

        for key in keys2:
            for (result, name) in displacement_like:
                if key not in result:
                    continue
                if isinstance(key, (int, int32)):
                    if key != subcase_id:
                        continue
                else:
                    subcase_idi = key[0]
                    if subcase_id != subcase_idi:
                        continue

                #print('dkey =', key)
                if key not in result:
                    continue
                case = result[key]
                subcase_idi = case.isubcase
                if not hasattr(case, 'data'):
                    continue
                if not case.is_real():
                    continue
                if case.nonlinear_factor is not None: # transient
                    code_name = case.data_code['name']
                    has_cycle = hasattr(case, 'mode_cycle')
                    assert case.is_sort1(), case.is_sort1()

                    itime0 = 0
                    t1 = case.data[itime0, :, 0]
                    ndata = t1.shape[0]
                    if nnodes != ndata:
                        nidsi = case.node_gridtype[:, 0]
                        assert len(nidsi) == nnodes
                        j = searchsorted(nids, nidsi)  # searching for nidsi

                        try:
                            if not allclose(nids[j], nidsi):
                                msg = 'nids[j]=%s nidsi=%s' % (nids[j], nidsi)
                                raise RuntimeError(msg)
                        except IndexError:
                            msg = 'node_ids = %s\n' % list(nids)
                            msg += 'nidsi in disp = %s\n' % list(nidsi)
                            raise IndexError(msg)


                    if name in ['Displacement', 'Eigenvectors']:
                        assert case.is_sort1(), case.is_sort1()
                        # we'll pass all the times in
                        t123 = case.data[:, :, :3]
                        if nnodes != ndata:
                            t123i = zeros((nnodes, 3), dtype='float32')
                            t123i[j, :] = t123
                            t123 = t123i
                        tnorm = norm(t123, axis=1)
                        assert len(tnorm) == t123.shape[0]
                        ntimes = case.ntimes
                        titles = []
                        scales = []
                        nastran_res = NastranDisplacementResults(subcase_idi, titles,
                                                                 self.xyz_cid0, t123, tnorm,
                                                                 scales,
                                                                 uname='NastranResult')

                        for itime in range(ntimes):
                            dt = case._times[itime]

                            if isinstance(dt, float):
                                header = ' %s = %.4E' % (code_name, dt)
                            else:
                                header = ' %s = %i' % (code_name, dt)

                            if has_cycle:
                                freq = case.eigrs[itime]
                                #msg.append('%16s = %13E\n' % ('EIGENVALUE', freq))
                                cycle = sqrt(abs(freq)) / (2. * pi)
                                header += '; freq=%g' % cycle

                            form0 = (header, None, [])
                            formi2 = form0[2]

                            print('*name = %r' % name)
                            tnorm_abs_max = tnorm.max()
                            scale = self.displacement_scale_factor / tnorm_abs_max

                            scale = self.dim_max / tnorm_abs_max * 0.25
                            scales.append(scale)

                            title = name + 'XYZ'
                            titles.append(title)

                            cases[icase] = (nastran_res, (itime, title))
                            formi2.append((title, icase, []))
                            icase += 1
                            formi.append(form0)
                        nastran_res.save_defaults()
                    else:
                        # print('name=%r' % name)
                        assert case.is_sort1(), case.is_sort1()
                        for itime in range(case.ntimes):
                            dt = case._times[itime]
                            t1 = case.data[itime, :, 0]
                            t2 = case.data[itime, :, 1]
                            t3 = case.data[itime, :, 2]

                            t123 = case.data[itime, :, :3]
                            #tnorm = norm(t123, axis=1)

                            if(t1.min() == t1.max() and t2.min() == t2.max() and
                               t3.min() == t3.max()):
                                continue
                            if nnodes != ndata:
                                t1i = zeros(nnodes, dtype='float32')
                                t2i = zeros(nnodes, dtype='float32')
                                t3i = zeros(nnodes, dtype='float32')
                                t1i[j] = t1
                                t2i[j] = t2
                                t3i[j] = t3
                                t1 = t1i
                                t2 = t2i
                                t3 = t3i

                            if isinstance(dt, float):
                                header = ' %s = %.4E' % (code_name, dt)
                            else:
                                header = ' %s = %i' % (code_name, dt)

                            if has_cycle:
                                freq = case.eigrs[itime]
                                #msg.append('%16s = %13E\n' % ('EIGENVALUE', freq))
                                cycle = sqrt(abs(freq)) / (2. * pi)
                                header += '; freq=%g' % cycle

                            form0 = (header, None, [])
                            formi2 = form0[2]

                            if 1:
                                cases[(subcase_idi, icase, name + 'X', 1, 'node', '%g', header)] = t1
                                formi2.append((name + 'X', icase, []))
                                icase += 1

                                cases[(subcase_idi, icase, name + 'Y', 1, 'node', '%g', header)] = t2
                                formi2.append((name + 'Y', icase, []))
                                icase += 1

                                cases[(subcase_idi, icase, name + 'Z', 1, 'node', '%g', header)] = t3
                                formi2.append((name + 'Z', icase, []))
                                icase += 1

                            cases[(subcase_idi, icase, name + 'XYZ', 3, 'node', '%g', header)] = t123
                            formi2.append((name + 'XYZ', icase, []))
                            icase += 1

                            formi.append(form0)
                else:
                    assert case.is_sort1(), case.is_sort1()

                    if name == 'Displacement' and 1:
                        t123 = case.data[:, :, :3]
                        tnorm = norm(t123, axis=1)
                        #print('*name = %r' % name)
                        titles = [name + 'XYZ']

                        tnorm_abs_max = tnorm.max()
                        scale = self.displacement_scale_factor / tnorm_abs_max
                        scales = [scale]
                        nastran_res = NastranDisplacementResults(subcase_idi, titles,
                                                                 self.xyz_cid0, t123, tnorm,
                                                                 scales=scales,
                                                                 uname='NastranResult')
                        nastran_res.save_defaults()
                        cases[icase] = (nastran_res, (0, name + 'XYZ'))
                        formi.append((name + 'XYZ', icase, []))
                        icase += 1
                    else:
                        t123 = case.data[0, :, :3]
                        tnorm = norm(t123, axis=1)
                        print('name = %r' % name)
                        t1 = case.data[0, :, 0]
                        t2 = case.data[0, :, 1]
                        t3 = case.data[0, :, 2]

                        if(t1.min() == t1.max() and t2.min() == t2.max() and
                           t3.min() == t3.max() and t123.min() == t123.max()):
                            continue
                        #if t1.min() != t1.max():
                        cases[(subcase_idi, icase, name + 'X', 1, 'node', '%g')] = t1
                        formi.append((name + 'X', icase, []))
                        icase += 1

                        #if t2.min() != t2.max():
                        cases[(subcase_idi, icase, name + 'Y', 1, 'node', '%g')] = t2
                        formi.append((name + 'Y', icase, []))
                        icase += 1

                        #if t3.min() != t3.max():
                        cases[(subcase_idi, icase, name + 'Z', 1, 'node', '%g')] = t3
                        formi.append((name + 'Z', icase, []))
                        icase += 1

                        #if t123.min() != t123.max():
                        cases[(subcase_idi, icase, name + 'XYZ', 1, 'node', '%g')] = case.data[0, :, :3]
                        #cases[(subcase_id, icase, name + 'XYZ', 1, 'node', '%g')] = tnorm
                        formi.append((name + 'XYZ', icase, []))
                        icase += 1

        for (result, name) in temperature_like:
            if subcase_id in result:
                case = result[subcase_id]
                subcase_idi = case.subcase_id
                if not hasattr(case, 'data'):
                    continue
                temperatures = case.data[0, :, 0]
                cases[(subcase_idi, name, 1, 'node', '%g')] = temperatures
                formi.append((name, icase, []))
                icase += 1
        return icase

    def clear_nastran(self):
        self.eidMap = {}
        self.nidMap = {}
        self.eid_to_nid_map = {}
        self.element_ids = None
        self.node_ids = None

    def fill_stress(self, cases, model, subcase_id, formi, icase):
        icase = self._fill_stress_centroidal(cases, model, subcase_id, formi, icase)
        #elif self.is_nodal:
            #icase = self._fill_stress_nodal(cases, model, subcase_id, formi, icase)
        #else:
            #raise RuntimeError('this shouldnt happen...')
        return icase

    def _fill_stress_nodal(self, cases, model, subcase_id, formi, icase):
        """
        disabled...
        """
        return icase

        is_stress = True
        if is_stress:
            word = 'Stress'
        else:
            word = 'Strain'

        oxx_dict = {}
        oyy_dict = {}
        ozz_dict = {}
        o1_dict = {}
        o2_dict = {}
        o3_dict = {}
        ovm_dict = {}

        for nid in self.nidMap:
            oxx_dict[nid] = []
            oyy_dict[nid] = []
            ozz_dict[nid] = []
            o1_dict[nid] = []
            o2_dict[nid] = []
            o3_dict[nid] = []
            ovm_dict[nid] = []

        vm_word = None
        if subcase_id in model.rodStress:
            case = model.rodStress[subcase_id]
            if case.nonlinear_factor is not None: # transient
                return
            for eid in case.axial:
                axial = case.axial[eid]
                torsion = case.torsion[eid]
                node_ids = self.eid_to_nid_map[eid]
                o1i = max(axial, torsion)  # not really
                o3i = min(axial, torsion)
                ovmi = max(abs(axial), abs(torsion))
                for nid in node_ids:
                    oxx_dict[nid].append(axial)
                    oyy_dict[nid].append(torsion)
                    o1_dict[nid].append(o1i)
                    o3_dict[nid].append(o3i)
                    ovm_dict[nid].append(ovmi)

        if is_stress:
            bars = model.barStress
        else:
            bars = model.barStrain

        if subcase_id in bars:
            case = bars[subcase_id]
            if case.nonlinear_factor is not None: # transient
                return
            for eid in case.axial:
                node_ids = self.eid_to_nid_map[eid]
                oxxi = case.axial[eid]
                o1i = max(case.smax[eid])
                o3i = min(case.smin[eid])
                ovmi = max(abs(max(case.smax[eid])),
                           abs(min(case.smin[eid])))
                for nid in node_ids:
                    oxx_dict[nid].append(oxxi)
                    o1_dict[nid].append(o1i)
                    o3_dict[nid].append(o3i)
                    ovm_dict[nid].append(ovmi)

        if subcase_id in model.beamStress:
            case = model.beamStress[subcase_id]
            if case.nonlinear_factor is not None: # transient
                return
            for eid in case.smax:
                node_ids = self.eid_to_nid_map[eid]
                oxxi = max(max(case.sxc[eid]),
                           max(case.sxd[eid]),
                           max(case.sxe[eid]),
                           max(case.sxf[eid]))
                o1i = max(case.smax[eid])
                o3i = min(case.smin[eid])
                ovmi = max(abs(max(case.smax[eid])),
                           abs(min(case.smin[eid])))
                for nid in node_ids:
                    oxx_dict[nid].append(oxxi)
                    o1_dict[nid].append(o1i)
                    o3_dict[nid].append(o3i)
                    ovm_dict[nid].append(ovmi)

        if subcase_id in model.plateStress:
            case = model.plateStress[subcase_id]
            if case.nonlinear_factor is not None: # transient
                return
            if case.is_von_mises():
                vm_word = 'vonMises'
            else:
                vm_word = 'maxShear'
            for eid in case.ovmShear:
                node_ids = self.eid_to_nid_map[eid]

                eType = case.eType[eid]
                if eType in ['CQUAD4', 'CQUAD8']:
                    #cen = 'CEN/%s' % eType[-1]
                    for nid in node_ids:
                        oxxi = max(case.oxx[eid][nid])
                        oyyi = max(case.oyy[eid][nid])
                        ozzi = min(case.oxx[eid][nid], min(case.oyy[eid][nid]))
                        o1i = max(case.majorP[eid][nid])
                        o2i = max(case.minorP[eid][nid])
                        o3i = min(case.majorP[eid][nid], min(case.minorP[eid][nid]))
                        ovmi = max(case.ovmShear[eid][nid])

                        oxx_dict[nid].append(oxxi)
                        oyy_dict[nid].append(oyyi)
                        o1_dict[nid].append(o1i)
                        o3_dict[nid].append(o3i)
                        ovm_dict[nid].append(ovmi)

                elif eType in ['CTRIA3', 'CTRIA6']:
                    cen = 'CEN/%s' % eType[-1]
                    oxxi = case.oxx[eid][cen]

                    oxxi = max(case.oxx[eid][cen])
                    oyyi = max(case.oyy[eid][cen])
                    ozzi = min(case.oxx[eid][cen], min(case.oyy[eid][cen]))

                    o1i = max(case.majorP[eid][cen])
                    o2i = max(case.minorP[eid][cen])
                    o3i = min(case.majorP[eid][cen], min(case.minorP[eid][cen]))
                    ovmi = max(case.ovmShear[eid][cen])

                    for nid in node_ids:
                        oxx_dict[nid].append(oxxi)
                        oyy_dict[nid].append(oyyi)
                        o1_dict[nid].append(o1i)
                        o3_dict[nid].append(o3i)
                        ovm_dict[nid].append(ovmi)

        if subcase_id in model.compositePlateStress:
            case = model.compositePlateStress[subcase_id]
            if case.nonlinear_factor is not None: # transient
                return
            if case.is_von_mises():
                vm_word = 'vonMises'
            else:
                vm_word = 'maxShear'

            for eid in case.ovmShear:
                node_ids = self.eid_to_nid_map[eid]

                oxxi = max(case.o11[eid])
                oyyi = max(case.o22[eid])
                o1i = max(case.majorP[eid])
                o3i = min(case.minorP[eid])
                ovmi = max(case.ovmShear[eid])

                for nid in node_ids:
                    oxx_dict[nid].append(oxxi)
                    oyy_dict[nid].append(oyyi)
                    o1_dict[nid].append(o1i)
                    o3_dict[nid].append(o3i)
                    ovm_dict[nid].append(ovmi)

        if subcase_id in model.solidStress:
            case = model.solidStress[subcase_id]
            if case.nonlinear_factor is not None: # transient
                return
            if case.is_von_mises():
                vm_word = 'vonMises'
            else:
                vm_word = 'maxShear'
            for eid in case.ovmShear:
                node_ids = self.eid_to_nid_map[eid]
                for nid in node_ids:
                    oxxi = case.oxx[eid][nid]
                    oyyi = case.oyy[eid][nid]
                    ozzi = case.ozz[eid][nid]
                    o1i = case.o1[eid][nid]
                    o2i = case.o2[eid][nid]
                    o3i = case.o3[eid][nid]
                    ovmi = case.ovmShear[eid][nid]

                    oxx_dict[nid].append(oxxi)
                    oyy_dict[nid].append(oyyi)
                    ozz_dict[nid].append(ozzi)
                    o1_dict[nid].append(o1i)
                    o2_dict[nid].append(o2i)
                    o3_dict[nid].append(o3i)
                    ovm_dict[nid].append(ovmi)

        nnodes = self.nNodes
        oxx = zeros(nnodes, dtype='float32')
        oyy = zeros(nnodes, dtype='float32')
        ozz = zeros(nnodes, dtype='float32')
        o1 = zeros(nnodes, dtype='float32')
        o2 = zeros(nnodes, dtype='float32')
        o3 = zeros(nnodes, dtype='float32')
        ovm = zeros(nnodes, dtype='float32')
        for i, nid in enumerate(sorted(self.nidMap)):
            oxx[i] = mean(oxx_dict[nid])
            oyy[i] = mean(oyy_dict[nid])
            ozz[i] = mean(ozz_dict[nid])
            o1[i] = mean(o1_dict[nid])
            o2[i] = mean(o2_dict[nid])
            o3[i] = mean(o3_dict[nid])
            ovm[i] = mean(ovm_dict[nid])

        # do this to prevent screwy stresses at points that have no stress
        oxx = nan_to_num(oxx)
        oyy = nan_to_num(oyy)
        ozz = nan_to_num(ozz)
        o1 = nan_to_num(o1)
        o2 = nan_to_num(o2)
        o3 = nan_to_num(o3)
        ovm = nan_to_num(ovm)
        if oxx.min() != oxx.max():
            cases[(subcase_id, icase, word + 'XX', 1, 'node', '%.3f')] = oxx
            icase += 1
        if oyy.min() != oyy.max():
            cases[(subcase_id, icase, word + 'YY', 1, 'node', '%.3f')] = oyy
            icase += 1
        if ozz.min() != ozz.max():
            cases[(subcase_id, icase, word + 'ZZ', 1, 'node', '%.3f')] = ozz
            icase += 1

        if o1.min() != o1.max():
            cases[(subcase_id, icase, word + '1', 1, 'node', '%.3f')] = o1
            icase += 1
        if o2.min() != o2.max():
            cases[(subcase_id, icase, word + '2', 1, 'node', '%.3f')] = o2
            icase += 1
        if o3.min() != o3.max():
            cases[(subcase_id, icase, word + '3', 1, 'node', '%.3f')] = o3
            icase += 1
        if vm_word is not None:
            cases[(subcase_id, icase, vm_word, 1, 'node', '%.3f')] = ovm
            icase += 1
        return icase

    def _is_nonlinear(self, model, isubcase):
        table_types = model.get_table_types()
        for table_type in table_types:
            table = getattr(model, table_type)
            if isubcase in table:
                case = table[isubcase]
                if case.nonlinear_factor:
                    return True
                else:
                    return False
        raise RuntimeError('self._is_nonlinear(...) failed')

    def _get_stress_times(self, model, isubcase):
        table_types = self._get_stress_table_types()
        is_real = True
        is_data = False
        is_static = False
        times = None
        for table_type in table_types:
            if not hasattr(model, table_type):
                print('no table_type=%s' % table_type)
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
                    times = zeros(1, dtype='int32')
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

    def _fill_stress_centroidal(self, cases, model, subcase_id, form, icase):
        assert isinstance(subcase_id, int), type(subcase_id)
        assert isinstance(icase, int), type(icase)
        is_data, is_static, is_real, times = self._get_stress_times(model, subcase_id)

        if not is_data:
            times = []

        #print('times = %s' % times)
        for itime, dt in enumerate(times):
            if is_static:
                formi = form
            else:
                header = 'dummy'
                form_time = [header, None, []]
                formi = form_time[2]
            is_form_time = False

            # stress
            icase, ncase, case, header, form0 = self._get_nastran_time_centroidal_stress(
                cases, model, subcase_id, form, icase, itime, dt,
                is_stress=True, is_real=is_real, is_static=is_static)
            if ncase:
                assert ncase > 0, ncase
                if is_static:
                    formi.append(form0)
                else:
                    form_time[0] = header
                    form_time[2].append(form0)
                    is_form_time = True

            # strain
            icase, ncase, case, header, form0 = self._get_nastran_time_centroidal_stress(
                cases, model, subcase_id, form, icase, itime, dt,
                is_stress=False, is_real=is_real, is_static=is_static)
            if ncase:
                assert ncase > 0, ncase
                if is_static:
                    formi.append(form0)
                else:
                    form_time[0] = header
                    form_time[2].append(form0)
                    is_form_time = True

            # ese
            icase, ncase, case, header, form0 = self._get_nastran_time_centroidal_strain_energy(
                cases, model, subcase_id, form, icase, itime, dt,
                is_real=is_real, is_static=is_static)
            if ncase:
                assert ncase > 0, ncase
                if is_static:
                    formi.append(form0)
                else:
                    #form_time[0] = header
                    form_time[2].append(form0)
                    is_form_time = True

            #--------------------------
            if is_form_time:
                form.append(form_time)

        return icase

    def _get_nastran_header(self, case, dt, itime):
        if case is None:
            return None
        code_name = case.data_code['name']

        if isinstance(dt, float):
            header = ' %s = %.4E' % (code_name, dt)
        else:
            header = ' %s = %i' % (code_name, dt)

        if hasattr(case, 'mode_cycle'):
            freq = case.eigrs[itime]
            #msg.append('%16s = %13E\n' % ('EIGENVALUE', freq))
            cycle = sqrt(abs(freq)) / (2. * pi)
            header += '; freq=%g' % cycle
        elif hasattr(case, 'eigrs'):
            freq = case.eigrs[itime]
            #msg.append('%16s = %13E\n' % ('EIGENVALUE', freq))
            cycle = sqrt(abs(freq)) / (2. * pi)
            header += '; freq=%g' % cycle
        return header

    def _get_nastran_time_centroidal_strain_energy(self, cases, model,
                                                   subcase_id, form, icase, itime, dt,
                                                   is_real=True, is_static=False):
        """
        Creates the time accurate strain energy objects for the pyNastranGUI
        """
        oxx = zeros(self.nElements, dtype='float32')
        oyy = zeros(self.nElements, dtype='float32')
        ozz = zeros(self.nElements, dtype='float32')
        fmt = '%g'
        header = ''
        ncase = 0
        form0 = ('Element Strain Energy', None, [])

        #op2.strain_energy[1]
          #type=StrainEnergyObject ntimes=3 nelements=16
          #energy, percent, density
          #modes = [1, 2, 3]
        case = None
        if subcase_id in model.strain_energy:
            ese = model.strain_energy[subcase_id]
            times = sorted(ese.energy.keys())  # TODO: not vectorized
            assert times[itime] == dt, 'actual=%s expected=%s' % (times[itime], dt)

            if is_static:
                percent = ese.percent
                energy = ese.energy
                density = ese.density
            else:
                percent = ese.percent[dt]
                energy = ese.energy[dt]
                density = ese.density[dt]
            for eid, p in sorted(iteritems(percent)):
                if eid not in self.eidMap:
                    continue
                i = self.eidMap[eid]
                oxx[i] = energy[eid]
                oyy[i] = p
                ozz[i] = density[eid]

            case = ese
            fmt = '%.4f'
            header = self._get_nastran_header(case, dt, itime)
            cases[(subcase_id, icase, 'StrainEnergy', 1, 'centroid', fmt, header)] = oxx
            form0[2].append(('StrainEnergy', icase, []))
            icase += 1
            ncase += 1

            cases[(subcase_id, icase, 'PercentOfTotal', 1, 'centroid', fmt, header)] = oyy
            form0[2].append(('PercentOfTotal', icase, []))
            icase += 1
            ncase += 1

            cases[(subcase_id, icase, 'Density', 1, 'centroid', fmt, header)] = ozz
            form0[2].append(('Density', icase, []))
            icase += 1
            ncase += 1
        return icase, ncase, case, header, form0

    def _get_nastran_time_centroidal_stress(self, cases, model, subcase_id, form, icase, itime, dt,
                                            is_stress=True, is_real=True, is_static=False):
        """
        Creates the time accurate stress objects for the pyNastranGUI
        """
        ncase = 0
        case = None
        assert isinstance(subcase_id, int), type(subcase_id)
        assert isinstance(icase, int), icase
        assert isinstance(itime, int), type(itime)
        assert is_real in [True, False], is_real
        assert is_stress in [True, False], is_stress
        assert is_static in [True, False], is_static
        eids = self.element_ids
        assert len(eids) > 0, eids
        nelements = self.nElements

        isElementOn = zeros(nelements, dtype='int8')  # is the element supported
        oxx = zeros(nelements, dtype='float32')
        oyy = zeros(nelements, dtype='float32')
        ozz = zeros(nelements, dtype='float32')

        txy = zeros(nelements, dtype='float32')
        tyz = zeros(nelements, dtype='float32')
        txz = zeros(nelements, dtype='float32')

        max_principal = zeros(nelements, dtype='float32')  # max
        mid_principal = zeros(nelements, dtype='float32')  # mid
        min_principal = zeros(nelements, dtype='float32')  # min
        ovm = zeros(nelements, dtype='float32')

        vm_word = None
        if is_stress:
            rods = [model.crod_stress, model.conrod_stress, model.ctube_stress,]
        else:
            rods = [model.crod_strain, model.conrod_strain, model.ctube_strain,]

        for result in rods:
            if subcase_id not in result:
                continue

            case = result[subcase_id]
            eidsi = case.element
            i = searchsorted(eids, eidsi)
            if len(i) != len(unique(i)):
                msg = 'irod=%s is not unique\n' % str(i)
                print('eids = %s\n' % str(list(eids)))
                print('eidsi = %s\n' % str(list(eidsi)))
                raise RuntimeError(msg)

            isElementOn[i] = 1

            # data=[1, nnodes, 4] where 4=[axial, SMa, torsion, SMt]
            oxx[i] = case.data[itime, :, 0]
            txy[i] = case.data[itime, :, 2]
            ovm[i] = sqrt(oxx[i]**2 + 3*txy[i]**2)
            max_principal[i] = sqrt(oxx[i]**2 + txy[i]**2)
            min_principal[i] = max_principal[i] - 2 * txy[i]
        del rods

        if is_stress:
            bars = model.cbar_stress
        else:
            bars = model.cbar_strain

        if subcase_id in bars:  # vectorized....
            case = bars[subcase_id]
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

            eidsi = case.element_node # [:, 0]

            i = searchsorted(eids, eidsi)
            if len(i) != len(unique(i)):
                print('ibar = %s' % i)
                print('eids = %s' % eids)
                msg = 'ibar=%s is not unique' % str(i)
                raise RuntimeError(msg)

            isElementOn[i] = 1.
            oxx[i] = axial

            ## TODO :not sure if this block is general for multiple CBAR elements
            samax = amax([smaxa, smaxb], axis=0)
            samin = amin([smaxa, smaxb], axis=0)
            assert len(samax) == len(i), len(samax)
            assert len(samin) == len(i)
            savm = amax(abs([smina, sminb,
                            smaxa, smaxb, axial]), axis=0)

            max_principal[i] = samax
            min_principal[i] = samin
            ovm[i] = savm
            del axial, smaxa, smina, smaxb, sminb, eidsi, i, samax, samin, savm
        del bars

        if is_stress:
            beams = model.cbeam_stress
        else:
            beams = model.cbeam_strain

        if subcase_id in beams:  # vectorized
            case = beams[subcase_id]
            eidsi = case.element_node[:, 0]
            ueids = unique(eidsi)
            #neids = len(ueids)

            j = 0
            # sxc, sxd, sxe, sxf
            # smax, smin, MSt, MSc
            sxc = case.data[itime, :, 0]
            sxd = case.data[itime, :, 1]
            sxe = case.data[itime, :, 2]
            sxf = case.data[itime, :, 3]
            smax = case.data[itime, :, 4]
            smin = case.data[itime, :, 5]
            for ieid, eid in enumerate(ueids):
                oxxi = 0.
                smaxi = 0.
                smini = 0.
                eid2 = self.eidMap[eid]
                isElementOn[eid2] = 1.
                for i in range(11):
                    oxxi = max(sxc[j], sxd[j], sxe[j], sxf[j], oxxi)
                    smaxi = max(smax[j], smaxi)
                    smini = min(smin[j], smini)
                    j += 1
                ovmi = max(abs(smaxi), abs(smini))
                oxxi = oxx[eid2]
                max_principal[eid2] = smaxi
                min_principal[eid2] = smini
                ovm[eid2] = ovmi
            del j, eidsi, ueids, sxc, sxd, sxe, sxf, smax, smin, oxxi, smaxi, smini, ovmi
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
            ## TODO: is tria6, quad8, bilinear quad handled?
            if subcase_id not in result:
                continue

            case = result[subcase_id]
            if case.is_von_mises():
                vm_word = 'vonMises'
            else:
                vm_word = 'maxShear'

            nnodes_per_element = case.nnodes
            nlayers_per_element = nnodes_per_element * 2  # *2 for every other layer
            eidsi = case.element_node[::nlayers_per_element, 0]  # ::2 is for layer skipping

            i = searchsorted(eids, eidsi)
            if len(i) != len(unique(i)):
                print('iplate = %s' % i)
                print('eids = %s' % eids)
                print('eidsiA = %s' % case.element_node[:, 0])
                print('eidsiB = %s' % eidsi)
                msg = 'iplate=%s is not unique' % str(i)
                raise RuntimeError(msg)
            #self.data[self.itime, self.itotal, :] = [fd, oxx, oyy,
            #                                         txy, angle,
            #                                         majorP, minorP, ovm]
            isElementOn[i] = 1.
            ntotal = case.data.shape[1]  # (ndt, ntotal, nresults)
            if nlayers_per_element == 1:
                j = None
            else:
                j = arange(ntotal)[::nlayers_per_element]

            #self.data[self.itime, self.itotal, :] = [fd, oxx, oyy,
            #                                         txy, angle,
            #                                         majorP, minorP, ovm]
            oxxi = case.data[itime, j, 1]
            oyyi = case.data[itime, j, 2]
            txyi = case.data[itime, j, 3]
            o1i = case.data[itime, j, 5]
            o3i = case.data[itime, j, 6]
            ovmi = case.data[itime, j, 7]

            for inode in range(1, nlayers_per_element):
                #print('%s - ilayer = %s' % (case.element_name, inode))
                oxxi = amax(vstack([oxxi, case.data[itime, j + inode, 1]]), axis=0)
                oyyi = amax(vstack([oyyi, case.data[itime, j + inode, 2]]), axis=0)
                txyi = amax(vstack([txyi, case.data[itime, j + inode, 3]]), axis=0)
                o1i = amax(vstack([o1i, case.data[itime, j + inode, 5]]), axis=0)
                o3i = amin(vstack([o3i, case.data[itime, j + inode, 6]]), axis=0)
                ovmi = amax(vstack([ovmi, case.data[itime, j + inode, 7]]), axis=0)
                assert len(oxxi) == len(j)

            oxx[i] = oxxi
            oyy[i] = oyyi
            txy[i] = txyi
            max_principal[i] = o1i
            min_principal[i] = o3i
            ovm[i] = ovmi

        if is_stress:
            cplates = [
                ('CTRIA3', model.ctria3_composite_stress), ('CQUAD4', model.cquad4_composite_stress),
                ('CTRIA6', model.ctria6_composite_stress), ('CQUAD8', model.cquad8_composite_stress),
                #model.ctriar_composite_stress, model.cquadr_composite_stress,
            ]
        else:
            cplates = [
                ('CTRIA3', model.ctria3_composite_strain), ('CQUAD4', model.cquad4_composite_strain),
                ('CTRIA6', model.ctria6_composite_strain), ('CQUAD8', model.cquad8_composite_strain),
                #model.ctriar_composite_strain, model.cquadr_composite_strain,
            ]

        for cell_type, result in cplates:
            if subcase_id not in result:
                continue

            case = result[subcase_id]
            if case.is_von_mises():
                vm_word = 'vonMises'
            else:
                vm_word = 'maxShear'

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
            for eid in unique(eidsi):
                ieid = where(eidsi == eid)[0]
                ieid.sort()
                layersi = layers[ieid]
                eid2 = self.eidMap[eid]
                isElementOn[eid2] = 1.

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
            if subcase_id not in result:
                continue

            case = result[subcase_id]
            if case.is_von_mises():
                vm_word = 'vonMises'
            else:
                vm_word = 'maxShear'

            nnodes_per_element = case.nnodes
            eidsi = case.element_cid[:, 0]
            ntotal = len(eidsi)  * nnodes_per_element

            i = searchsorted(eids, eidsi)
            if len(i) != len(unique(i)):
                print('isolid = %s' % str(i))
                print('eids = %s' % eids)
                print('eidsi = %s' % eidsi)
                assert len(i) == len(unique(i)), 'isolid=%s is not unique' % str(i)

            isElementOn[i] = 1
            #self.data[self.itime, self.itotal, :] = [oxx, oyy, ozz,
            #                                         txy, tyz, txz,
            #                                         o1, o2, o3, ovm]

            if nnodes_per_element == 1:
                j = None
            else:
                j = arange(ntotal)[::nnodes_per_element]
                ueidsi = unique(eidsi)
                assert len(j) == len(ueidsi), 'j=%s ueidsi=%s' % (j, ueidsi)

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
                #print('%s - inode = %s' % (case.element_name, inode))
                oxxi = amax(vstack([oxxi, case.data[itime, j + inode, 0]]), axis=0)
                oyyi = amax(vstack([oyyi, case.data[itime, j + inode, 1]]), axis=0)
                ozzi = amax(vstack([ozzi, case.data[itime, j + inode, 2]]), axis=0)
                txyi = amax(vstack([txyi, case.data[itime, j + inode, 3]]), axis=0)
                tyzi = amax(vstack([tyzi, case.data[itime, j + inode, 4]]), axis=0)
                txzi = amax(vstack([txzi, case.data[itime, j + inode, 2]]), axis=0)

                o1i = amax(vstack([o1i, case.data[itime, j + inode, 6]]), axis=0)
                o2i = amax(vstack([o2i, case.data[itime, j + inode, 7]]), axis=0)
                o3i = amin(vstack([o3i, case.data[itime, j + inode, 8]]), axis=0)
                ovmi = amax(vstack([ovmi, case.data[itime, j + inode, 9]]), axis=0)
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

        header = ''
        if not is_static:
            #print('is_static = %s' % is_static)
            if case is None:
                formis = None
                return icase, ncase, case, header, formis
            header = self._get_nastran_header(case, dt, itime)
            #form_time[0] = header

        form0 = (word, None, [])
        formis = form0[2]
        # subcase_id, icase, resultType, vector_size, location, dataFormat
        if is_stress and itime == 0:
            if isElementOn.min() == 0:  # if all elements aren't on
                ioff = where(isElementOn == 0)[0]
                print('eids_off = %s' % self.element_ids[ioff])
                self.log_error('eids_off = %s' % self.element_ids[ioff])
                cases[(1, icase, 'isElementOn', 1, 'centroid', '%i')] = isElementOn
                form.append(('IsElementOn', icase, []))
                icase += 1
                ncase += 1

        if oxx.min() != oxx.max():
            cases[(subcase_id, icase, word + 'XX', 1, 'centroid', fmt, header)] = oxx
            formis.append((word + 'XX', icase, []))
            icase += 1
            ncase += 1
        if oyy.min() != oyy.max():
            cases[(subcase_id, icase, word + 'YY', 1, 'centroid', fmt, header)] = oyy
            formis.append((word + 'YY', icase, []))
            icase += 1
            ncase += 1
        if ozz.min() != ozz.max():
            cases[(subcase_id, icase, word + 'ZZ', 1, 'centroid', fmt, header)] = ozz
            formis.append((word + 'ZZ', icase, []))
            ncase += 1
            icase += 1

        if txy.min() != txy.max():
            cases[(subcase_id, icase, word + 'XY', 1, 'centroid', fmt, header)] = txy
            formis.append((word + 'XY', icase, []))
            icase += 1
            ncase += 1
        if tyz.min() != tyz.max():
            cases[(subcase_id, icase, word + 'YZ', 1, 'centroid', fmt, header)] = tyz
            formis.append((word + 'YZ', icase, []))
            icase += 1
            ncase += 1
        if txz.min() != txz.max():
            cases[(subcase_id, icase, word + 'XZ', 1, 'centroid', fmt, header)] = txz
            formis.append((word + 'XZ', icase, []))
            icase += 1
            ncase += 1

        if max_principal.min() != max_principal.max():
            cases[(subcase_id, icase, 'MaxPrincipal', 1, 'centroid', fmt, header)] = max_principal
            formis.append(('Max Principal', icase, []))
            icase += 1
            ncase += 1
        if mid_principal.min() != mid_principal.max():
            cases[(subcase_id, icase, 'MidPrincipal', 1, 'centroid', fmt, header)] = mid_principal
            formis.append(('Mid Principal', icase, []))
            icase += 1
            ncase += 1
        if min_principal.min() != min_principal.max():
            cases[(subcase_id, icase, 'MinPrincipal', 1, 'centroid', fmt, header)] = min_principal
            formis.append(('Min Principal', icase, []))
            icase += 1
            ncase += 1
        if vm_word is not None:
            if not is_stress:
                max_min = max(ovm.max(), abs(ovm.min()))
                if max_min > 100:
                    raise RuntimeError('vm strain = %s' % ovm)
            cases[(subcase_id, icase, vm_word, 1, 'centroid', fmt, header)] = ovm
            formis.append((vm_word, icase, []))
            icase += 1
            ncase += 1
        return icase, ncase, case, header, form0

