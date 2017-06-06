# pylint: disable=E1101
"""
Defines the GUI IO file for Nastran.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import os
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
from pyNastran.gui.qt_version import qt_version
if qt_version == 4:
    import PyQt4
elif qt_version == 5:
    import PyQt5
else:
    raise NotImplementedError(qt_version)

import vtk
from vtk import (vtkTriangle, vtkQuad, vtkTetra, vtkWedge, vtkHexahedron,
                 vtkQuadraticTriangle, vtkQuadraticQuad, vtkQuadraticTetra,
                 vtkQuadraticWedge, vtkQuadraticHexahedron,
                 vtkPyramid) #vtkQuadraticPyramid

#from pyNastran import is_release
from pyNastran.utils import integer_types
from pyNastran.bdf.bdf import (BDF,
                               CAERO1, CAERO2, CAERO3, CAERO4, CAERO5,
                               CQUAD4, CQUAD8, CQUAD, CQUADR, CSHEAR,
                               CTRIA3, CTRIA6, CTRIAR,
                               CPLSTN3, CPLSTN4, CPLSTN6, CPLSTN8,
                               CTRAX3, CTRIAX6, CTRIAX, #CTRAX6,
                               CQUADX4, CQUADX8, CQUADX,
                               CONM2)

from pyNastran.bdf.cards.elements.shell import ShellElement
from pyNastran.bdf.cards.elements.solid import (
    CTETRA4, CTETRA10, CPENTA6, CPENTA15,
    CHEXA8, CHEXA20, CIHEX1, CIHEX2,
    CPYRAM5, CPYRAM13,
)

from pyNastran.gui.errors import NoGeometry
from pyNastran.gui.gui_objects.gui_result import GuiResult, NormalResult
from pyNastran.converters.nastran.geometry_helper import (
    NastranGeometryHelper, tri_quality, quad_quality, get_min_max_theta)
from pyNastran.converters.nastran.results_helper import NastranGuiResults
from pyNastran.converters.nastran.displacements import (
    ForceTableResults)

from pyNastran.op2.op2 import OP2
#from pyNastran.f06.f06_formatting import get_key0
try:
    from pyNastran.op2.op2_geom import OP2Geom
    is_geom = True
except ImportError:
    is_geom = False

GREEN = (0., 1., 0.)
BLUE = (0., 0., 1.)
LIGHT_GREEN = (0.5, 1., 0.5)
PINK = (0.98, 0.4, 0.93)
ORANGE = (219/255., 168/255., 13/255.)
RED = (1., 0., 0.)
YELLOW = (1., 1., 0.)
PURPLE = (1., 0., 1.)


class NastranIO(NastranGuiResults, NastranGeometryHelper):
    """
    Defines the GUI class for Nastran.
    """
    def __init__(self):
        super(NastranIO, self).__init__()

    def get_nastran_wildcard_geometry_results_functions(self):
        """
        gets the Nastran wildcard loader used in the file load menu
        """
        if is_geom:
            geom_methods1 = 'Nastran BDF ''(*.bdf; *.dat; *.nas; *.ecd; *.op2; *.pch)'
            geom_methods2 = 'Nastran Punch (*.bdf; *.dat; *.nas; *.ecd; *.op2; *.pch)'
        else:
            geom_methods1 = 'Nastran BDF ''(*.bdf; *.dat; *.nas; *.ecd; *.pch)'
            geom_methods2 = 'Nastran Punch (*.bdf; *.dat; *.nas; *.ecd; *.pch)'

        data_bdf = (
            'nastran',
            geom_methods1, self.load_nastran_geometry,
            'Nastran OP2 (*.op2)', self.load_nastran_results)
        data_pch = (
            'nastran',
            geom_methods2, self.load_nastran_geometry,
            'Nastran OP2 (*.op2)', self.load_nastran_results)

        return [data_bdf, data_pch]

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
        """
        Create a coordinate system

        Parameters
        ----------
        dim_max : float
            the max model dimension; 10% of the max will be used for the coord length
        cid : int
           the coordinate system id
        coord : Coord()
           the Nastran coord object
        cid_type : str
            a string of 'xyz', 'Rtz', 'Rtp' (xyz, cylindrical, spherical)
            that changes the axis names
        """
        origin = coord.origin
        beta = coord.beta().T
        self.create_coordinate_system(dim_max, label='%s' % cid, origin=origin,
                                      matrix_3x3=beta, Type=cid_type)

    def _create_nastran_coords(self, model, dim_max):
        """
        Creates the Nastran coordinate systems.

        Parameters
        ----------
        model : BDF()
            the BDF object
        dim_max : float
            the max model dimension; 10% of the max will be used for the coord length
        """
        cid_types = {
            'R' : 'xyz',
            'C' : 'Rtz',
            'S' : 'Rtp',
        }
        self.create_global_axes(dim_max)
        for cid, coord in sorted(iteritems(model.coords)):
            if cid == 0:
                continue
            cid_type = cid_types[coord.Type]
            self._create_coord(dim_max, cid, coord, cid_type)

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

    def get_xyz_in_coord(self, model, cid=0, fdtype='float32'):
        """Creates the grid points efficiently"""
        #import time
        #t0 = time.time()

        if 1:
            # t=.578
            #print("get_displacement_index_xyz_cp_cd")
            out = model.get_displacement_index_xyz_cp_cd(
                fdtype=fdtype, idtype='int32', sort_ids=True)
            icd_transform, icp_transform, xyz_cp, nid_cp_cd = out
            self.i_transform = icd_transform

            #print("transform_xyzcp_to_xyz_cid")
            xyz_cid0 = model.transform_xyzcp_to_xyz_cid(xyz_cp, icp_transform, cid=0,
                                                        in_place=False)

            nid_map = self.nid_map
            for i, nid in enumerate(nid_cp_cd[:, 0]):
                nid_map[nid] = i
        elif 0:  # pragma: no cover
            # t=.573
            out = model.get_displacement_index_xyz_cp_cd(
                fdtype='float32', idtype='int32', sort_ids=True)
            icd_transform, icp_transform, xyz_cp, nid_cp_cd = out
            self.i_transform = icd_transform
            xyz_cid0 = model.transform_xyzcp_to_xyz_cid(xyz_cp, icp_transform, cid=0)

            nid_map = self.nid_map
            for i, nid in enumerate(nid_cp_cd[:, 0]):
                nid_map[nid] = i
        else:  # pragma: no cover
            # t=.75
            nid_map = self.nid_map
            assert cid == 0, cid
            nnodes = len(model.nodes)
            nspoints = 0
            spoints = None
            if model.spoints:
                spoints = model.spoints.points
                nspoints = len(spoints)

            xyz_cid0 = np.zeros((nnodes + nspoints, 3), dtype=fdtype)
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
            else:
                for i, (nid, node) in enumerate(sorted(iteritems(model.nodes))):
                    xyz = node.get_position()
                    xyz_cid0[i, :] = xyz
                    nid_map[nid] = i

            # get indicies and transformations for displacements
            #self.i_transform, self.transforms = model.get_displacement_index_transforms()
            self.i_transform = model.get_displacement_index()

        self._add_nastran_spoints_to_grid(model)
        #print('dt_nastran_xyz =', time.time() - t0)
        return xyz_cid0, nid_cp_cd

    def _get_model(self, bdf_filename, xref_loads=True):
        """Loads the BDF/OP2 geometry"""
        #print('get_model')
        ext = os.path.splitext(bdf_filename)[1].lower()
        punch = False
        if ext == '.pch':
            punch = True

        self.model_type = 'nastran'
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
            model.read_bdf(bdf_filename,
                           punch=punch, xref=False,
                           validate=True)
            #print('done with read_bdf')

            #xref_loads = False
            xref_aero = len(model.caeros) > 0
            #model.cross_reference(
            model.safe_cross_reference(
                xref=True,
                xref_nodes=True,
                xref_elements=True,
                xref_nodes_with_elements=False,
                xref_properties=True,
                xref_masses=True,
                xref_materials=False,
                xref_loads=xref_loads,
                xref_constraints=False,
                xref_optimization=False,
                xref_aero=True,
                xref_sets=False,
            )
            #print('done with cross_reference')
        return model

    def load_nastran_geometry(self, bdf_filename, dirname, name='main', plot=True):
        """
        The entry point for Nastran geometry loading.

        Parameters
        ----------
        bdf_filename : str
            the Nastran filename to load
        dirname : str
            ???
        name : str
            the name of the "main" actor for the GUI
        plot : bool; default=True
            should the model be generated or should we wait until
            after the results are loaded
        """
        self.eid_maps[name] = {}
        self.nid_maps[name] = {}
        self.i_transform = {}
        #self.transforms = {}
        #print('bdf_filename=%r' % bdf_filename)
        #key = self.case_keys[self.icase]
        #case = self.result_cases[key]

        skip_reading = self._remove_old_nastran_geometry(bdf_filename)
        # if 0:
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

        xref_loads = True # should be True
        model = self._get_model(bdf_filename, xref_loads=xref_loads)

        nnodes = len(model.nodes)
        nspoints = 0
        spoints = None
        if model.spoints:
            spoints = model.spoints.points
            nspoints = len(spoints)

        if nnodes + nspoints == 0:
            msg = 'nnodes + nspoints = 0\n'
            msg += 'card_count = %r' % str(model.card_count)
            raise NoGeometry(msg)

        nelements = len(model.elements)
        nplotels = len(model.plotels)
        ncaero_cards = len(model.caeros)
        if nelements + ncaero_cards + nplotels == 0:
            msg = 'nelements + ncaero_cards + nplotels = 0\n'
            msg += 'card_count = %r' % str(model.card_count)
            raise NoGeometry(msg)

        self.nNodes = nnodes + nspoints
        self.nElements = nelements  # approximate...

        out = self.make_caeros(model)
        (caero_points, ncaeros, ncaeros_sub, ncaeros_cs,
         ncaeros_points, ncaero_sub_points,
         has_control_surface, box_id_to_caero_element_map, cs_box_ids) = out

        self.log_info("nNodes=%i nElements=%i" % (self.nNodes, self.nElements))
        msg = model.get_bdf_stats(return_type='string')
        self.log_debug(msg)
        msg = model.get_bdf_stats(return_type='list')

        # this call will break the GUI if there are a lot of lines and
        # by a lot I mean 37641.  It's fine for a single call.
        #for msgi in msg:
            #model.log.debug(msgi)

        nconm2 = 0
        if 'CONM2' in model.card_count:
            nconm2 += model.card_count['CONM2']
        if 'CMASS1' in model.card_count:
            nconm2 += model.card_count['CMASS1']
        if 'CMASS2' in model.card_count:
            nconm2 += model.card_count['CMASS2']

        if nconm2 > 0:
            self.create_alternate_vtk_grid(
                'conm2', color=ORANGE, line_width=5, opacity=1., point_size=4,
                representation='point')

        # Allocate grids
        self.grid.Allocate(self.nElements, 1000)
        self._create_caero_actors(ncaeros, ncaeros_sub, ncaeros_cs, has_control_surface)
        if nconm2 > 0:
            self.alt_grids['conm2'].Allocate(nconm2, 1000)

        #vectorReselt.SetNumberOfComponents(3)
        #elem.SetNumberOfPoints(nNodes)

        if self.save_data:
            self.model = model

        #print('get_xyz_in_coord')
        xyz_cid0, nid_cp_cd = self.get_xyz_in_coord(model, cid=0, fdtype='float32')
        points = self.numpy_to_vtk_points(xyz_cid0)
        self.xyz_cid0 = xyz_cid0

        maxi = xyz_cid0.max(axis=0)
        mini = xyz_cid0.min(axis=0)
        assert len(maxi) == 3, len(maxi)
        xmax, ymax, zmax = maxi
        xmin, ymin, zmin = mini
        dim_max = max(xmax-xmin, ymax-ymin, zmax-zmin)

        #print('_create_nastran_coords')
        self._create_nastran_coords(model, dim_max)
        #print('done _create_nastran_coords')

        self.log_info("xmin=%s xmax=%s dx=%s" % (xmin, xmax, xmax-xmin))
        self.log_info("ymin=%s ymax=%s dy=%s" % (ymin, ymax, ymax-ymin))
        self.log_info("zmin=%s zmax=%s dz=%s" % (zmin, zmax, zmax-zmin))

        j = 0
        nid_to_pid_map, icase, cases, form = self.map_elements(
            points, self.nid_map, model, j, dim_max, nid_cp_cd,
            plot=plot, xref_loads=xref_loads)

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
        self.create_splines(model, box_id_to_caero_element_map, caero_points)
        if 'caero' in self.alt_grids:
            self.set_caero_grid(ncaeros_points, model, j=0)
            self.set_caero_subpanel_grid(ncaero_sub_points, model, j=0)
            if has_control_surface:
                self.set_caero_control_surface_grid(
                    'caero_control_surfaces', cs_box_ids,
                    box_id_to_caero_element_map, caero_points)
        if nconm2 > 0:
            self.set_conm_grid(nconm2, dim_max, model)

        make_spc_mpc_supports = True
        if make_spc_mpc_supports:
            geometry_names = self.set_spc_mpc_suport_grid(dim_max, model, nid_to_pid_map)
        else:
            geometry_names = []
        icase = self._fill_bar_yz(dim_max, model, icase, cases, form)
        assert icase is not None

        #------------------------------------------------------------

        #print('dependent_nodes =', self.dependents_nodes)
        form0 = form[2]
        assert icase is not None
        for subcase_id, subcase in sorted(iteritems(model.subcases)):
            if subcase_id == 0:
                continue
            self.log_debug('NastranIOv subcase_id = %s' % subcase_id)

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
                        model, cases, formii, icase, subcase_id)
                except KeyError:
                    s = StringIO()
                    traceback.print_exc(file=s)
                    sout = s.getvalue()
                    self.log_error(sout)
                    print(sout)

            assert icase is not None
            if self.plot_applied_loads:
                try:
                    icase = self._plot_applied_loads(
                        model, cases, formii, icase, subcase_id, xref_loads=xref_loads)
                except KeyError:
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

        for grid_name in geometry_names:
            if grid_name in self.geometry_actors:
                self.geometry_actors[grid_name].Modified()

        #self.grid_mapper.SetResolveCoincidentTopologyToPolygonOffset()
        if plot:
            #self.log.info(cases.keys())
            self._finish_results_io2([form], cases)
        else:
            self._set_results([form], cases)

    def _create_caero_actors(self, ncaeros, ncaeros_sub, ncaeros_cs, has_control_surface):
        if self.has_caero:
            if 'caero' not in self.alt_grids:
                self.create_alternate_vtk_grid(
                    'caero', color=YELLOW, line_width=3, opacity=1.0,
                    representation='toggle', is_visible=True, is_pickable=False)
            if 'caero_subpanels' not in self.alt_grids:
                self.create_alternate_vtk_grid(
                    'caero_subpanels', color=YELLOW, line_width=3, opacity=1.0,
                    representation='toggle', is_visible=False, is_pickable=False)

            self.alt_grids['caero'].Allocate(ncaeros, 1000)
            self.alt_grids['caero_subpanels'].Allocate(ncaeros_sub, 1000)
            if has_control_surface:
                self.alt_grids['caero_control_surfaces'].Allocate(ncaeros_cs, 1000)

    def make_caeros(self, model):
        ncaeros = 0
        ncaeros_sub = 0
        ncaeros_cs = 0
        ncaeros_points = 0
        ncaero_sub_points = 0
        has_control_surface = False
        box_id_to_caero_element_map = {}
        cs_box_ids = []

        # when caeros is empty, SPLINEx/AESURF cannot be defined
        if len(model.caeros) == 0:
            #print('returning early')
            caero_points = np.empty((0, 3))
            out = (
                caero_points, ncaeros, ncaeros_sub, ncaeros_cs,
                ncaeros_points, ncaero_sub_points,
                has_control_surface, box_id_to_caero_element_map, cs_box_ids,
            )
            return out

        # count caeros
        # sorting doesn't matter here because we're just trying to size the array
        for caero in itervalues(model.caeros):
            if hasattr(caero, 'panel_points_elements'):
                npoints, ncelements = caero.get_npanel_points_elements()
                ncaeros_sub += npoints
                ncaero_sub_points += ncelements
            elif isinstance(caero, CAERO2):
                pass
            else:
                print('%r doesnt support panel_points_elements' % caero.type)

        for eid, caero in sorted(iteritems(model.caeros)):
            if isinstance(caero, (CAERO1, CAERO3, CAERO4, CAERO5)):
                ncaeros_points += 4
                ncaeros += 1
            elif isinstance(caero, CAERO2):
                points, elems = caero.get_points_elements_3d()
                ncaeros_points += points.shape[0]
                ncaeros += elems.shape[0]

        num_prev = 0
        ncaeros_sub = 0
        if model.caeros:
            caero_points = []
            for eid, caero in sorted(iteritems(model.caeros)):
                if caero.type in ['CAERO1', 'CAERO4']:
                    ncaeros_sub += 1
                    pointsi, elementsi = caero.panel_points_elements()
                    caero_points.append(pointsi)
                    for i, box_id in enumerate(caero.box_ids.flat):
                        box_id_to_caero_element_map[box_id] = elementsi[i, :] + num_prev
                    num_prev += pointsi.shape[0]
                elif caero.type == 'CAERO2':
                    pass
                else:
                    print('caero\n%s' % caero)
            if ncaeros_sub:
                caero_points = np.vstack(caero_points)
            self.has_caero = True
        if ncaeros_sub == 0:
            caero_points = np.empty((0, 3))

        # check for any control surfcaes
        if model.aesurfs:
            has_control_surface = True
            #ncaero_cs_points = 0
            if 'caero_control_surfaces' not in self.alt_grids:
                self.create_alternate_vtk_grid(
                    'caero_control_surfaces', color=PINK, line_width=5, opacity=1.0,
                    representation='surface')
            for aid, aesurf in iteritems(model.aesurfs):
                aelist = aesurf.alid1
                ncaeros_cs += len(aelist.elements)
                cs_box_ids.extend(aelist.elements)

                if aesurf.alid2 is not None:
                    aelist = aesurf.alid2
                    ncaeros_cs += len(aelist.elements)
                    cs_box_ids.extend(aelist.elements)
        out = (
            caero_points, ncaeros, ncaeros_sub, ncaeros_cs,
            ncaeros_points, ncaero_sub_points,
            has_control_surface, box_id_to_caero_element_map, cs_box_ids,
        )
        return out

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
        """
        Sets the following actors:
          - spline_%s_structure_points % spline_id
          - spline_%s_boxes % spline_id

        Parameters
        ----------
        model : BDF()
            the bdf model
        box_id_to_caero_element_map : dict[key] : value
            ???
        caero_points : ???
            ???
        """
        if model.splines:
            # 0 - caero / caero_subpanel
            # 1 - control surface
            iaero = 2
            for spline_id, spline in sorted(iteritems(model.splines)):
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

                zfighting_offset = 0.0001 * (iaero + 1)
                grid_name = 'spline_%s_structure_points' % spline_id
                self.create_alternate_vtk_grid(
                    grid_name, color=BLUE, opacity=1.0, point_size=5,
                    representation='point', is_visible=False)
                msg = ', which is required by %r' % grid_name
                self._add_nastran_nodes_to_grid(grid_name, structure_points, model, msg)


                zfighting_offset = 0.0001 * (iaero + 2)
                grid_name = 'spline_%s_boxes' % spline_id
                self.create_alternate_vtk_grid(
                    grid_name, color=BLUE, opacity=0.3,
                    line_width=4,
                    representation='toggle', is_visible=False)
                self.set_caero_control_surface_grid(
                    grid_name, aero_box_ids, box_id_to_caero_element_map, caero_points,
                    zfighting_offset=zfighting_offset)
                iaero += 2

    def set_caero_grid(self, ncaeros_points, model, j=0):
        """
        Sets the CAERO panel geometry.

        Parameters
        ----------
        ncaeros_points : int
            number of points used by the 'caero' actor
        model : BDF()
            the bdf model
        j : int; default=0
            the current id counter (ID of what???)

        Returns
        -------
        j : int
            the current id counter (ID of what???)
        """
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(ncaeros_points)

        max_cpoints = []
        min_cpoints = []

        zfighting_offset = 0.0001
        caero_grid = self.alt_grids['caero']
        for eid, element in sorted(iteritems(model.caeros)):
            if isinstance(element, (CAERO1, CAERO3, CAERO4, CAERO5)):
                # wing panel
                cpoints = element.get_points()
                cpoints[0][2] += zfighting_offset
                cpoints[1][2] += zfighting_offset
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
                    cpoints[:, 2] +=  zfighting_offset
                    max_cpoints.append(np.array(cpoints).max(axis=0))
                    min_cpoints.append(np.array(cpoints).min(axis=0))

                    elem = vtk.vtkLine()
                    elem.GetPointIds().SetId(0, j)
                    elem.GetPointIds().SetId(1, j + 1)

                    #print(', '.join(dir(elem)))
                    #prop = elem.GetProperty()

                    points.InsertPoint(j, *cpoints[0])
                    points.InsertPoint(j + 1, *cpoints[1])
                    j += 2
                    caero_grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                else:
                    # 3D version
                    xyz, elems = element.get_points_elements_3d()
                    xyz[:, 2] +=  zfighting_offset
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
        #self.alt_grids['caero']
        #edge_mapper.SetResolveCoincidentTopologyToPolygonOffset()
        return j

    def set_caero_subpanel_grid(self, ncaero_sub_points, model, j=0):
        """
        Sets the CAERO sub-panel geometry.

        Parameters
        ----------
        ncaero_sub_points : int
            number of points used by the 'caero_subpanels' actor
        model : BDF()
            the bdf model
        j : int; default=0
            the current id counter (ID of what???)

        Returns
        -------
        j : int
            the current id counter (ID of what???)
        """
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(ncaero_sub_points)

        vtk_type = vtkQuad().GetCellType()
        for eid, element in sorted(iteritems(model.caeros)):
            if isinstance(element, (CAERO1, CAERO3, CAERO4, CAERO5)):
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
        """
        Creates a single CAERO control surface?

        Parameters
        ----------
        name : str
            ???
        aero_box_ids : List[int]
            the ids of the box as seen on the AESURF? SET card?
        box_id_to_caero_element_map : Dict[key]=value
            key : ???
                ???
            value : ???
                ???
        caero_points : (ncaero_points, 3)
            the xyz coordinates used by the CAEROx actor
        zfighting_offset : float
            z-fighting is when two elements "fight" for who is in front
            leading.  The standard way to fix this is to bump the
            element.
        j : int; default=0
            ???

        Returns
        -------
        j : int
            ???
        """
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
        for box_id in cs_box_ids:
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

    def set_conm_grid(self, nconm2, dim_max, model, j=0):
        """
        creates the mass secondary actor called:
         - conm2

        which includes:
         - CONM2
         - CMASS1
         - CMASS2

        because it's really a "mass" actor
        """
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(nconm2)

        #sphere_size = self._get_sphere_size(dim_max)
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
            elif element.type in ['CMASS1', 'CMASS2']:
                n1 = element.G1()
                n2 = element.G1()
                centroid = element.Centroid()
                #print('n1=%s n2=%s centroid=%s' % (n1, n2, centroid))
                points.InsertPoint(j, *centroid)

                elem = vtk.vtkVertex()
                elem.GetPointIds().SetId(0, j)
                self.alt_grids['conm2'].InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                j += 1

            else:
                self.log_info("skipping %s" % element.type)
        self.alt_grids['conm2'].SetPoints(points)

    def set_spc_mpc_suport_grid(self, dim_max, model, nid_to_pid_map):
        """
        for each subcase, make secondary actors including:
         - spc_id=spc_id
         - mpc_id=mpc_id              (includes rigid elements)
         - mpc_dependent_id=mpc_id    (includes rigid elements)
         - mpc_independent_id=mpc_id  (includes rigid elements)
         - suport_id=suport1_id       (includes SUPORT/SUPORT1)

        TODO: consider changing the varying ids to huh???
        """
        spc_names = []
        mpc_names = []
        suport_names = []

        #print('getting rigid')
        rigid_lines = model._get_rigid()

        spc_ids_used = set()
        mpc_ids_used = set()
        suport1_ids_used = set()
        for subcase_id, subcase in sorted(iteritems(model.subcases)):
            #print(subcase.params.keys())
            if 'SPC' in subcase:
                #print('getting spcs')
                spc_id = subcase.get_parameter('SPC')[0]
                if spc_id is not None and spc_id not in spc_ids_used:
                    spc_ids_used.add(spc_id)
                    nspcs = model.card_count['SPC'] if 'SPC' in model.card_count else 0
                    nspc1s = model.card_count['SPC1'] if 'SPC1' in model.card_count else 0
                    nspcds = model.card_count['SPCD'] if 'SPCD' in model.card_count else 0

                    ## TODO: this line seems too loose...
                    ## TODO: why isn't SPCDs included?
                    if nspcs + nspc1s + nspcds:
                        spc_names += self._fill_spc(spc_id, nspcs, nspc1s, nspcds, dim_max,
                                                    model, nid_to_pid_map)

            # rigid body elements and MPCs
            if 'MPC' in subcase:
                #print('getting mpc')
                mpc_id = subcase.get_parameter('MPC')[0]
                if mpc_id is not None and mpc_id not in mpc_ids_used:
                    mpc_ids_used.add(mpc_id)

                    ## TODO: this line seems too loose
                    nmpcs = model.card_count['MPC'] if 'MPC' in model.card_count else 0
                    if nmpcs:
                        lines = model.get_MPCx_node_ids_c1(mpc_id, exclude_mpcadd=False)
                        lines2 = list(lines) + rigid_lines
                        mpc_names += self._fill_dependent_independent(
                            mpc_id, dim_max, model, lines2, nid_to_pid_map)

            # SUPORTs are node/dofs that deconstrained to allow rigid body motion
            # SUPORT1s are subcase-specific SUPORT cards
            if 'SUPORT1' in subcase.params:  ## TODO: should this be SUPORT?
                suport_id = subcase.get_parameter('SUPORT1')[0]

                # TODO: is this line correct???
                if 'SUPORT' in model.card_count or 'SUPORT1' in model.card_count:

                    # TODO: this "if block" seems unnecessary
                    if suport_id is not None and suport_id not in suport1_ids_used:
                        # SUPORT1 / SUPORT
                        suport1_ids_used.add(suport_id)
                        suport_name = self._fill_suport(suport_id, subcase_id, dim_max, model)
                        suport_names.append(suport_name)

        # create a SUPORT actor if there are no SUPORT1s
        # otherwise, we already included it in suport_id=suport_id
        if len(suport_names) == 0 and model.suport:
            # handle SUPORT without SUPORT1
            ids = []
            for suport in model.suport:
                idsi = suport.node_ids
                ids += idsi
            grid_name = 'SUPORT'
            self.create_alternate_vtk_grid(
                grid_name, color=RED, opacity=1.0, point_size=4,
                representation='point', is_visible=True)

        if len(mpc_names) == 0 and len(rigid_lines):
            # handle RBEs without MPCs
            mpc_id = 0
            mpc_names += self._fill_dependent_independent(
                mpc_id, dim_max, model, rigid_lines, nid_to_pid_map)

        geometry_names = spc_names + mpc_names + suport_names
        return geometry_names

    def _fill_spc(self, spc_id, nspcs, nspc1s, nspcds, dim_max, model, nid_to_pid_map):
        """creates the spc secondary actors"""
        spc_name = 'spc_id=%i' % spc_id
        spc_names = [spc_name]
        self.create_alternate_vtk_grid(spc_name, color=PURPLE, line_width=5, opacity=1.,
                                       point_size=5, representation='point', is_visible=False)

        # node_ids = model.get_SPCx_node_ids(spc_id, exclude_spcadd=False)
        node_ids_c1 = model.get_SPCx_node_ids_c1(
            spc_id, exclude_spcadd=False, stop_on_failure=False)

        node_ids = []
        for nid, c1 in iteritems(node_ids_c1):
            if nid_to_pid_map is not None:
                plot_node = False
                pids = nid_to_pid_map[nid]
                for pid in pids:
                    if pid == 0:
                        # CONROD
                        continue
                    if pid is None:
                        print('pid is None in _fill_spc...')
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
        msg = ', which is required by %r' % spc_name
        self._add_nastran_nodes_to_grid(spc_name, node_ids, model, msg, nid_to_pid_map)
        return spc_names

    def create_bar_pin_flag_text(self, pin_flag=None):
        """TODO: needs a better interface in the gui"""
        nids = []
        text = []
        #result_name = self.icase
        result_name = str('ElementID')
        for nid, data in sorted(iteritems(self.nid_release_map)):
            sub_release_map = defaultdict(str)
            for (eid, pin_flagi) in data:
                sub_release_map[pin_flagi] += (str(eid) + ', ')
            texti = '\n'.join(['%s-%s'  % (pin_flagi, msg.rstrip(', '))
                               for (pin_flagi, msg) in sorted(iteritems(sub_release_map))])

            # super messy
            #texti = ', '.join(['%s-%s'  % (pin_flagi, eid) for (eid, pin_flagi) in data])
            nids.append(nid)
            text.append(texti)
        self.mark_nodes(nids, result_name, text)

    def _fill_bar_yz(self, dim_max, model, icase, cases, form, debug=False):
        """
        plots the y, z vectors for CBAR & CBEAM elements
        """
        card_types = ['CBAR', 'CBEAM']
        out = model.get_card_ids_by_card_types(card_types=card_types)
        bar_beam_eids = out['CBAR'] + out['CBEAM']

        self.bar_eids = {}
        self.bar_lines = {}
        if len(bar_beam_eids) == 0:
            return icase
        scale = 0.15

        # TODO: this should be reworked
        bar_nids, bar_types, nid_release_map = self._get_bar_yz_arrays(
            model, bar_beam_eids, scale, debug)
        self.nid_release_map = nid_release_map

        bar_nids = list(bar_nids)
        self.create_alternate_vtk_grid(
            'Bar Nodes', color=RED, line_width=1, opacity=1.,
            point_size=5, representation='point', bar_scale=0., is_visible=False)
        msg = ", which is required by 'Bar Nodes'"
        self._add_nastran_nodes_to_grid('Bar Nodes', bar_nids, model, msg)


        geo_form = form[2]
        bar_form = ('CBAR / CBEAM', None, [])
        #print('geo_form =', geo_form)
        bar_types2 = {}
        for bar_type, data in sorted(iteritems(bar_types)):
            eids, lines_bar_y, lines_bar_z = data
            if len(eids):
                if debug: # pragma: no cover
                    print('bar_type = %r' % bar_type)
                    print('eids     = %r' % eids)
                    print('all_eids = %r' % self.element_ids.tolist())
                # if bar_type not in ['ROD', 'TUBE']:
                bar_y = bar_type + '_y'
                bar_z = bar_type + '_z'

                self.create_alternate_vtk_grid(
                    bar_y, color=GREEN, line_width=5, opacity=1.,
                    point_size=5, representation='bar', bar_scale=scale, is_visible=False)
                self.create_alternate_vtk_grid(
                    bar_z, color=BLUE, line_width=5, opacity=1.,
                    point_size=5, representation='bar', bar_scale=scale, is_visible=False)

                self._add_nastran_lines_xyz_to_grid(bar_y, lines_bar_y, eids, model)
                self._add_nastran_lines_xyz_to_grid(bar_z, lines_bar_z, eids, model)

                # form = ['Geometry', None, []]
                i = np.searchsorted(self.element_ids, eids)
                is_type = np.zeros(self.element_ids.shape, dtype='int32')
                try:
                    is_type[i] = 1.
                except:
                    #print('self.element_ids =', self.element_ids)
                    #print('eids =', eids)
                    ii = np.where(i == len(self.element_ids))[0]
                    print('ii = %s' % ii)
                    print('failed eids =', eids[ii])
                    #assert self.element_ids[i] == eids
                    raise
                bar_form[2].append(['is_%s' % bar_type, icase, []])

                msg = 'is_%s' % bar_type
                type_res = GuiResult(0, header=msg, title=msg,
                                     location='centroid', scalar=is_type)
                cases[icase] = (type_res, (0, msg))
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

        assert name != u'Bar Nodes', name
        grid = self.alt_grids[name]

        bar_eids = np.array(eids, dtype='int32')
        bar_lines = np.asarray(lines, dtype='float32').reshape(nlines, 6)
        self.bar_eids[name] = np.asarray(eids, dtype='int32')
        self.bar_lines[name] = bar_lines

        nodes = bar_lines.reshape(nlines * 2, 3)
        points = self.numpy_to_vtk_points(nodes)
        elements = np.arange(0, nnodes, dtype='int32').reshape(nlines, 2)

        etype = 3 # vtk.vtkLine().GetCellType()
        self.create_vtk_cells_of_constant_element_type(grid, elements, etype)
        grid.SetPoints(points)

    def _fill_dependent_independent(self, mpc_id, dim_max, model, lines, nid_to_pid_map):
        """creates the mpc actors"""
        if not lines:
            return

        #print('_fill_dependent_independent')
        if mpc_id == 0:
            depname = 'rigid_dependent'
            indname = 'rigid_independent'
            linename = 'rigid_lines'
        else:
            depname = 'mpc_id=%i_dependent' % mpc_id
            indname = 'mpc_id=%i_independent' % mpc_id
            linename = 'mpc_id=%i_lines' % mpc_id
        self.create_alternate_vtk_grid(
            depname, color=GREEN, line_width=5, opacity=1.,
            point_size=5, representation='point', is_visible=False)
        self.create_alternate_vtk_grid(
            indname, color=LIGHT_GREEN, line_width=5, opacity=1.,
            point_size=5, representation='point', is_visible=False)
        self.create_alternate_vtk_grid(
            linename, color=LIGHT_GREEN, line_width=5, opacity=1.,
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

        msg = ', which is required by %r' % depname
        self._add_nastran_nodes_to_grid(depname, dependent, model, msg, nid_to_pid_map)

        msg = ', which is required by %r' % indname
        self._add_nastran_nodes_to_grid(indname, independent, model, msg, nid_to_pid_map)

        msg = ', which is required by %r' % linename
        self._add_nastran_lines_to_grid(linename, lines, model, nid_to_pid_map)

        mpc_names = [depname, indname, linename]
        return mpc_names

    def _add_nastran_nodes_to_grid(self, name, node_ids, model, msg, nid_to_pid_map=None):
        """used to create MPC independent/dependent nodes"""
        nnodes = len(node_ids)
        if nnodes == 0:
            model.log.warning('0 nodes added for %r' % name)
            return
        self.follower_nodes[name] = node_ids

        #self.numpy_to_vtk_points(nodes)
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(nnodes)

        j = 0
        nid_map = self.nid_map
        for nid in sorted(node_ids):
            try:
                i = nid_map[nid]
            except KeyError:
                model.log.warning('nid=%s does not exist%s' % (nid, msg))
                continue

            if nid not in model.nodes:
                # I think this hits for SPOINTs
                model.log.warning('nid=%s doesnt exist%s' % (nid, msg))
                continue
            # point = self.grid.GetPoint(i)
            # points.InsertPoint(j, *point)

            node = model.nodes[nid]
            point = node.get_position()
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

    def _add_nastran_spoints_to_grid(self, model, nid_to_pid_map=None):
        """used to create SPOINTs"""
        if model.spoints is None:
            return
        spoint_ids = list(model.spoints.points) # set -> list
        assert isinstance(spoint_ids, list), type(spoint_ids)

        nspoints = len(spoint_ids)
        name = 'SPoints'
        if nspoints == 0:
            model.log.warning('0 spoints added for %r' % name)
            return
        self.create_alternate_vtk_grid(
            name, color=BLUE, line_width=1, opacity=1.,
            point_size=5, representation='point', bar_scale=0., is_visible=True)

        self.follower_nodes[name] = spoint_ids
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(nspoints)

        j = 0
        nid_map = self.nid_map
        for spointi in sorted(spoint_ids):
            try:
                i = nid_map[spointi]
            except KeyError:
                model.log.warning('spointi=%s does not exist' % spointi)
                continue

            if spointi not in model.spoints.points:
                model.log.warning('spointi=%s doesnt exist' % spointi)
                continue
            # point = self.grid.GetPoint(i)
            # points.InsertPoint(j, *point)

            points.InsertPoint(j, 0., 0., 0.)

            elem = vtk.vtkVertex()
            elem.GetPointIds().SetId(0, j)
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
        self.follower_nodes[name] = lines.ravel()
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

    def _fill_suport(self, suport_id, subcase_id, dim_max, model):
        """creates SUPORT and SUPORT1 nodes"""
        suport_name = 'suport1_id=%i' % suport_id
        self.create_alternate_vtk_grid(
            suport_name, color=RED, line_width=5, opacity=1., point_size=4,
            representation='point', is_visible=False)
        suport_nids = self._get_suport_node_ids(model, suport_id)
        msg = ', which is required by %r' % suport_name
        self._add_nastran_nodes_to_grid(suport_name, suport_nids, model, msg)
        return suport_name

    def _get_sphere_size(self, dim_max):
        return 0.01 * dim_max

    def map_elements2(self, points, nid_map, model, j, dim_max,
                      nid_cp_cd, plot=True, xref_loads=True):  # pragma: no cover
        #model = BDF()
        cards_to_consider = [
            'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
            'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4',
            'CBAR', 'CBEAM', 'CROD', 'CONROD', 'CTUBE',
            'CTRIA3', 'CTRIAR', 'CQUAD4',
            'CTRIA6', 'CQUAD8', 'CQUAD',
            'CTRIAX', 'CTRIAX6',
            'CQUADX', 'CQUADX4', 'CQUADX8',
            'CTETRA', 'CPENTA', 'CHEXA', 'CPYRAM',
        ]
        non_element_cards = [
            'PBUSH', 'PELAS', 'PROD', 'PSHELL', 'PCOMP', 'PBARL', 'PSOLID', 'PBAR',
            'PBEAM', 'PBEAML',

            'MAT1', 'MAT2', 'MAT3', 'MAT4', 'MAT5', 'MAT8', 'MAT9', 'MAT10', 'MAT11',
            'MAT3D', 'MATHE', 'MATHP',
            'SUPORT', 'SUPORT1', 'EIGR', 'EIGRL', 'EIGB', 'EIGC',
            'GRID', 'CORD1R', 'CORD1C', 'CORD1S', 'CORD2R', 'CORD2C', 'CORD2S',
        ]
        grid = self.grid
        all_eids = []
        all_pids = []
        all_nids = nid_cp_cd[:, 0]
        #print(model._type_to_id_map)
        for card_type in model._type_to_id_map:
            #print('card_type = %r' % card_type)
            if card_type in cards_to_consider:
                eids = model._type_to_id_map[card_type]
                if card_type in ['CBAR', 'CBEAM', 'CROD', 'CTUBE']:
                    nid = np.array([model.elements[eid].node_ids for eid in eids],
                                   dtype='int32')
                    pids = np.array([model.elements[eid].Pid() for eid in eids],
                                    dtype='int32')
                    inids = np.searchsorted(all_nids, nid)
                    nnodes = 2

                    for elem_nid in inids:
                        elem = vtk.vtkLine()
                        pts = elem.GetPointIds()
                        pts.SetId(0, elem_nid[0])
                        pts.SetId(1, elem_nid[1])
                        grid.InsertNextCell(elem.GetCellType(), pts)

                elif card_type in ['CTRIA3', 'CTRIAR']:
                    nnodes = 3
                    nid = np.array([model.elements[eid].node_ids for eid in eids],
                                   dtype='int32')
                    pids = np.array([model.elements[eid].Pid() for eid in eids],
                                    dtype='int32')
                    theta_mcid = [model.elements[eid].theta_mcid for eid in eids]
                    inids = np.searchsorted(all_nids, nid)

                    for elem_nid in inids:
                        elem = vtkTriangle()
                        pts = elem.GetPointIds()
                        pts.SetId(0, elem_nid[0])
                        pts.SetId(1, elem_nid[1])
                        pts.SetId(2, elem_nid[2])
                        grid.InsertNextCell(elem.GetCellType(), pts)

                elif card_type in ['CQUAD4']:
                    nnodes = 4
                    nid = np.array([model.elements[eid].node_ids for eid in eids],
                                    dtype='int32')
                    pids = np.array([model.elements[eid].Pid() for eid in eids],
                                    dtype='int32')
                    inids = np.searchsorted(all_nids, nid)

                    for elem_nid in inids:
                        elem = vtkQuad()
                        pts = elem.GetPointIds()
                        pts.SetId(0, elem_nid[0])
                        pts.SetId(1, elem_nid[1])
                        pts.SetId(2, elem_nid[2])
                        pts.SetId(3, elem_nid[3])
                        grid.InsertNextCell(elem.GetCellType(), pts)
                elif card_type == 'CTETRA':
                    pids_list = []
                    nids1 = []
                    nids2 = []
                    for eid in eids:
                        elem = model.elements[eid]
                        pid = elem.Pid()
                        node_ids = elem.node_ids
                        if len(node_ids) == 4:
                            nids1.append(node_ids)
                        else:
                            nids2.append(node_ids)
                        pids_list.append(pid)
                    pids = np.array(pids_list, dtype='int32')

                    if nids1:
                        nnodes = 4
                        nids1 = np.array(nids1, dtype='int32')
                        inids1 = np.searchsorted(all_nids, nids1)
                        for elem_nid in inids1:
                            elem = vtkTetra()
                            pts = elem.GetPointIds()
                            pts.SetId(0, elem_nid[0])
                            pts.SetId(1, elem_nid[1])
                            pts.SetId(2, elem_nid[2])
                            pts.SetId(3, elem_nid[3])
                            grid.InsertNextCell(elem.GetCellType(), pts)

                    if nids2:
                        nnodes = 10
                        nids2 = np.array(nids2, dtype='int32')
                        inids2 = np.searchsorted(all_nids, nids2)
                        for elem_nid in inids2:
                            elem = vtkQuadraticTetra()
                            pts = elem.GetPointIds()
                            pts.SetId(0, elem_nid[0])
                            pts.SetId(1, elem_nid[1])
                            pts.SetId(2, elem_nid[2])
                            pts.SetId(3, elem_nid[3])
                            pts.SetId(4, elem_nid[4])
                            pts.SetId(5, elem_nid[5])
                            pts.SetId(6, elem_nid[6])
                            pts.SetId(7, elem_nid[7])
                            pts.SetId(8, elem_nid[8])
                            pts.SetId(9, elem_nid[9])
                            grid.InsertNextCell(elem.GetCellType(), pts)


                elif card_type in ['CHEXA']:
                    #continue
                    pids_list = []
                    nids1 = []
                    nids2 = []
                    eids1 = []
                    eids2 = []
                    for eid in eids:
                        elem = model.elements[eid]
                        pid = elem.Pid()
                        node_ids = elem.node_ids
                        if len(node_ids) == 8:
                            eids1.append(eid)
                            nids1.append(node_ids)
                        else:
                            eids2.append(eid)
                            nids2.append(node_ids)
                        pids_list.append(pid)

                    if eids1 and eids2:
                        eids = eids1 + eids2

                    pids = np.array(pids_list, dtype='int32')

                    if nids1:
                        nids1 = np.array(nids1, dtype='int32')
                        inids1 = np.searchsorted(all_nids, nids1)
                        for elem_nid in inids1:
                            nnodes = 8
                            elem = vtkHexahedron()
                            pts = elem.GetPointIds()
                            pts.SetId(0, elem_nid[0])
                            pts.SetId(1, elem_nid[1])
                            pts.SetId(2, elem_nid[2])
                            pts.SetId(3, elem_nid[3])
                            pts.SetId(4, elem_nid[4])
                            pts.SetId(5, elem_nid[5])
                            pts.SetId(6, elem_nid[6])
                            pts.SetId(7, elem_nid[7])
                            grid.InsertNextCell(elem.GetCellType(), pts)

                    if nids2:
                        nids2 = np.array(nids2, dtype='int32')
                        print(nids2)
                        inids2 = np.searchsorted(all_nids, nids2)
                        for elem_nid in inids2:
                            nnodes = 20
                            elem = vtkQuadraticHexahedron()
                            pts = elem.GetPointIds()
                            pts.SetId(0, elem_nid[0])
                            pts.SetId(1, elem_nid[1])
                            pts.SetId(2, elem_nid[2])
                            pts.SetId(3, elem_nid[3])
                            pts.SetId(4, elem_nid[4])
                            pts.SetId(5, elem_nid[5])
                            pts.SetId(6, elem_nid[6])
                            pts.SetId(7, elem_nid[7])
                            pts.SetId(8, elem_nid[8])
                            pts.SetId(9, elem_nid[9])
                            pts.SetId(10, elem_nid[10])
                            pts.SetId(11, elem_nid[11])
                            pts.SetId(12, elem_nid[12])
                            pts.SetId(13, elem_nid[13])
                            pts.SetId(14, elem_nid[14])
                            pts.SetId(15, elem_nid[15])
                            pts.SetId(16, elem_nid[16])
                            pts.SetId(17, elem_nid[17])
                            pts.SetId(18, elem_nid[18])
                            pts.SetId(19, elem_nid[19])
                            grid.InsertNextCell(elem.GetCellType(), pts)

                #elif card_type in ['CPENTA']:
                #elif card_type in ['CPYRAM']:

                elif card_type in non_element_cards:
                    continue
                else:
                    print('card_type = %r' % card_type)
                    continue
                all_eids.append(eids)
                all_pids.append(pids)

        if len(all_eids) == 0:
            raise RuntimeError('all_eids=0 ... huh?')
        if len(all_pids) == 0:
            raise RuntimeError('all_pids=0 ... huh?')

        if len(all_eids) == 1:
            all_eids = all_eids[0]
            all_pids = all_pids[0]
        else:
            all_eids = np.hstack(all_eids)
            all_pids = np.hstack(all_pids)

        eids = all_eids
        pids = all_pids
        ieids = np.argsort(eids)
        #print(ieids)
        self.ieids = ieids
        # we need some
        #eids.sort()
        #-------------------------------------------------------------
        for i, eid in enumerate(eids):
            self.eid_map[eid] = i

        # sorted = unsorted[isort]
        eids = eids[ieids]
        pids = pids[ieids]

        nelements = len(all_eids)
        nid_to_pid_map = None
        #celas1s = model._type_to_id_map['CELAS1']
        #celas2s = model._type_to_id_map['CELAS2']
        #celas3s = model._type_to_id_map['CELAS3']
        #celas4s = model._type_to_id_map['CELAS4']

        #cdamp1s = model._type_to_id_map['CDAMP1']
        #cdamp2s = model._type_to_id_map['CDAMP2']
        #cdamp3s = model._type_to_id_map['CDAMP3']
        #cdamp4s = model._type_to_id_map['CDAMP4']

        #cbars = model._type_to_id_map['CBAR']
        #cbeams = model._type_to_id_map['CBEAM']
        #crods = model._type_to_id_map['CROD']
        #conrods = model._type_to_id_map['CONROD']
        #ctubes = model._type_to_id_map['CTUBE']

        ## simple
        #ctria3s = model._type_to_id_map['CTRIA3']
        #ctriars = model._type_to_id_map['CTRIAR']
        #cquad4s = model._type_to_id_map['CQUAD4']

        ## curved
        #ctria6s = model._type_to_id_map['CTRIA6']
        #cquad8s = model._type_to_id_map['CQUAD8']
        #cquads = model._type_to_id_map['CQUAD']

        ## axisymmetric
        #ctriaxs = model._type_to_id_map['CTRIAX']   # the sane CTRIA6-X
        #ctriax6s = model._type_to_id_map['CTRIAX6'] # the dumb CTRIA6-X

        #cquadxs = model._type_to_id_map['CQUADX']
        #cquadx4s = model._type_to_id_map['CQUADX4']
        #cquadx8s = model._type_to_id_map['CQUADX8']

        #ctetras = model._type_to_id_map['CTETRA']
        #cpentas = model._type_to_id_map['CPENTA']
        #chexas = model._type_to_id_map['CHEXA']
        #cpyrams = model._type_to_id_map['CPYRAM']

        #nelements = i
        self.nElements = nelements
        #print('nelements=%s pids=%s' % (nelements, list(pids)))
        #pids = pids[:nelements]

        grid.Modified()
        if hasattr(grid, 'Update'):
            grid.Update()
        #self.log_info("updated grid")

        cases = OrderedDict()
        #del pids_dict

        self.iSubcaseNameMap = {1: ['Nastran', '']}
        icase = 0
        form = ['Geometry', None, []]
        form0 = form[2]

        #new_cases = True
        # set to True to enable node_ids as an result
        nids_set = True
        if nids_set:
            # this intentionally makes a deepcopy
            nids = np.array(nid_cp_cd[:, 0])
            cds = np.array(nid_cp_cd[:, 2])

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
        if eids_set and nelements:
            #eids = np.zeros(nelements, dtype='int32')
            #for (eid, eid2) in iteritems(self.eid_map):
                #eids[eid2] = eid
            assert isinstance(eids, np.ndarray), type(eids)

            eid_res = GuiResult(0, header='ElementID', title='ElementID',
                                location='centroid', scalar=eids)
            cases[icase] = (eid_res, (0, 'ElementID'))
            form0.append(('ElementID', icase, []))
            icase += 1
            self.element_ids = eids

        # subcase_id, resultType, vector_size, location, dataFormat
        if len(model.properties) and 0:
            icase, upids, mids, thickness, nplies, is_pshell_pcomp = self._build_properties(
                model, nelements, eids, pids, cases, form0, icase)
            icase = self._build_materials(model, mids, thickness, nplies, is_pshell_pcomp,
                                          cases, form0, icase)

            try:
                icase = self._build_optimization(model, pids, upids, nelements, cases, form0, icase)
            except:
                #raise
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

        #self.glyphs.SetScaleFactor(min_edge_length.mean())
        #if self.make_offset_normals_dim and nelements:
            #icase, normals = self._build_normals_quality(
                #model, nelements, cases, form0, icase,
                #xyz_cid0, material_coord,
                #min_interior_angle, max_interior_angle, dideal_theta,
                #area, max_skew_angle, taper_ratio,
                #max_warp_angle, area_ratio, min_edge_length, max_aspect_ratio)
            #self.normals = normals
        return nid_to_pid_map, icase, cases, form


    def map_elements(self, points, nid_map, model, j, dim_max,
                     nid_cp_cd, plot=True, xref_loads=True):
        """
        Creates the elements
        """
        grid = self.grid
        grid.SetPoints(points)
        #return self.map_elements2(points, nid_map, model, j, dim_max,
                                  #nid_cp_cd, plot=plot, xref_loads=xref_loads)

        out = self._map_elements(model, dim_max, nid_map, j)
        (nid_to_pid_map, xyz_cid0, pids, nelements, material_coord,
         area, min_interior_angle, max_interior_angle, max_aspect_ratio,
         max_skew_angle, taper_ratio, dideal_theta,
         area_ratio, min_edge_length, max_warp_angle) = out

        #self.grid_mapper.SetResolveCoincidentTopologyToPolygonOffset()
        grid.Modified()
        if hasattr(grid, 'Update'):
            grid.Update()
        #self.log_info("updated grid")

        cases = OrderedDict()


        self.iSubcaseNameMap = {1: ['Nastran', '']}
        icase = 0
        form = ['Geometry', None, []]
        form0 = form[2]

        #new_cases = True
        # set to True to enable node_ids as an result
        nids_set = True
        if nids_set:
            # this intentionally makes a deepcopy
            nids = np.array(nid_cp_cd[:, 0])
            cds = np.array(nid_cp_cd[:, 2])

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
        if eids_set and nelements:
            eids = np.zeros(nelements, dtype='int32')
            for (eid, eid2) in iteritems(self.eid_map):
                eids[eid2] = eid

            eid_res = GuiResult(0, header='ElementID', title='ElementID',
                                location='centroid', scalar=eids)
            cases[icase] = (eid_res, (0, 'ElementID'))
            form0.append(('ElementID', icase, []))
            icase += 1
            self.element_ids = eids

        # subcase_id, resultType, vector_size, location, dataFormat
        if len(model.properties):
            icase, upids, mids, thickness, nplies, is_pshell_pcomp = self._build_properties(
                model, nelements, eids, pids, cases, form0, icase)
            icase = self._build_materials(model, mids, thickness, nplies, is_pshell_pcomp,
                                          cases, form0, icase)

            try:
                icase = self._build_optimization(model, pids, upids, nelements, cases, form0, icase)
            except:
                #raise
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

        #print('nelements=%s eid_map=%s' % (nelements, self.eid_map))
        self.set_glyph_scale_factor(np.nanmean(min_edge_length) * 2.5)  # was 1.5
        if self.make_offset_normals_dim and nelements:
            icase, normals = self._build_normals_quality(
                model, nelements, cases, form0, icase,
                xyz_cid0, material_coord,
                min_interior_angle, max_interior_angle, dideal_theta,
                area, max_skew_angle, taper_ratio,
                max_warp_angle, area_ratio, min_edge_length, max_aspect_ratio)
            self.normals = normals
        return nid_to_pid_map, icase, cases, form

    def _map_elements(self, model, dim_max, nid_map, j):
        """
        Helper for map_elements


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

        Aspect Ratio:
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
        min_edge_length = np.zeros(nelements, 'float32')

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
            (0, 3, 2), # (1, 3, 4),
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
        pid = 0

        grid = self.grid
        nplotels = len(model.plotels)
        if nplotels:
            lines = []
            for (eid, element) in sorted(iteritems(model.plotels)):
                node_ids = element.node_ids
                lines.append(node_ids)
            lines = np.array(lines, dtype='int32')

            self.create_alternate_vtk_grid(
                'plotel', color=RED, line_width=2, opacity=0.8,
                point_size=5, representation='wire', is_visible=True)
            self._add_nastran_lines_to_grid('plotel', lines, model)

        #print("map_elements...")
        for (eid, element) in sorted(iteritems(model.elements)):
            self.eid_map[eid] = i
            if i % 5000 == 0 and i > 0:
                print('  map_elements = %i' % i)
            etype = element.type
            # if element.Pid() >= 82:
                # continue
            # if element.Pid() in pids_to_drop:
                # continue
            # if element.Pid() not in pids_to_keep:
                # continue
            # if element.pid.type == 'PSOLID':
                # continue
            if 1:
                pid = np.nan
                dideal_thetai = np.nan
                min_thetai = np.nan
                max_thetai = 0.0
                #max_thetai = np.nan
                max_skew = np.nan
                #max_warp = np.nan
                max_warp = 0.0
                aspect_ratio = np.nan
                areai = np.nan
                area_ratioi = np.nan
                taper_ratioi = np.nan
                min_edge_lengthi = np.nan
            else:
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
                min_edge_lengthi = 0.
            if isinstance(element, (CTRIA3, CTRIAR, CTRAX3, CPLSTN3)):
                if isinstance(element, (CTRIA3, CTRIAR)):
                    material_coord[i] = 0 if isinstance(element.theta_mcid, float) else element.theta_mcid
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
                out = tri_quality(p1, p2, p3)
                (areai, max_skew, aspect_ratio,
                 min_thetai, max_thetai, dideal_thetai, min_edge_lengthi) = out

                elem.GetPointIds().SetId(0, n1)
                elem.GetPointIds().SetId(1, n2)
                elem.GetPointIds().SetId(2, n3)
                grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element, (CTRIA6, CPLSTN6, CTRIAX)):
                # the CTRIAX is a standard 6-noded element
                if isinstance(element, CTRIA6):
                    material_coord[i] = 0 if isinstance(element.theta_mcid, float) else element.theta_mcid
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

                n1, n2, n3 = [nid_map[nid] for nid in node_ids[:3]]
                p1 = xyz_cid0[n1, :]
                p2 = xyz_cid0[n2, :]
                p3 = xyz_cid0[n3, :]
                out = tri_quality(p1, p2, p3)
                (areai, max_skew, aspect_ratio,
                 min_thetai, max_thetai, dideal_thetai, min_edge_lengthi) = out
                elem.GetPointIds().SetId(0, n1)
                elem.GetPointIds().SetId(1, n2)
                elem.GetPointIds().SetId(2, n3)
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element, CTRIAX6):
                # the CTRIAX6 is not a standard second-order triangle
                #
                # 5
                # |\
                # |  \
                # 6    4
                # |     \
                # |       \
                # 1----2----3
                #
                material_coord[i] = element.theta_mcid
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

                n1 = nid_map[node_ids[0]]
                n2 = nid_map[node_ids[2]]
                n3 = nid_map[node_ids[4]]
                p1 = xyz_cid0[n1, :]
                p2 = xyz_cid0[n2, :]
                p3 = xyz_cid0[n3, :]
                out = tri_quality(p1, p2, p3)
                (areai, max_skew, aspect_ratio,
                 min_thetai, max_thetai, dideal_thetai, min_edge_lengthi) = out
                elem.GetPointIds().SetId(0, n1)
                elem.GetPointIds().SetId(1, n2)
                elem.GetPointIds().SetId(2, n3)
                self.eid_to_nid_map[eid] = [node_ids[0], node_ids[2], node_ids[4]]
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())

            elif isinstance(element, (CQUAD4, CSHEAR, CQUADR, CPLSTN4, CQUADX4)):
                if isinstance(element, (CQUAD4, CQUADR)):
                    material_coord[i] = 0 if isinstance(element.theta_mcid, float) else element.theta_mcid
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
                out = quad_quality(p1, p2, p3, p4)
                (areai, taper_ratioi, area_ratioi, max_skew, aspect_ratio,
                 min_thetai, max_thetai, dideal_thetai, min_edge_lengthi) = out

                elem = vtkQuad()
                elem.GetPointIds().SetId(0, n1)
                elem.GetPointIds().SetId(1, n2)
                elem.GetPointIds().SetId(2, n3)
                elem.GetPointIds().SetId(3, n4)
                self.grid.InsertNextCell(9, elem.GetPointIds())

            elif isinstance(element, (CQUAD8, CPLSTN8, CQUADX8)):
                if isinstance(element, CQUAD8):
                    material_coord[i] = 0 if isinstance(element.theta_mcid, float) else element.theta_mcid
                node_ids = element.node_ids
                pid = element.Pid()
                for nid in node_ids:
                    if nid is not None:
                        nid_to_pid_map[nid].append(pid)
                self.eid_to_nid_map[eid] = node_ids[:4]

                n1, n2, n3, n4 = [nid_map[nid] for nid in node_ids[:4]]
                p1 = xyz_cid0[n1, :]
                p2 = xyz_cid0[n2, :]
                p3 = xyz_cid0[n3, :]
                p4 = xyz_cid0[n4, :]
                out = quad_quality(p1, p2, p3, p4)
                (areai, taper_ratioi, area_ratioi, max_skew, aspect_ratio,
                 min_thetai, max_thetai, dideal_thetai, min_edge_lengthi) = out
                if None not in node_ids:
                    elem = vtkQuadraticQuad()
                    elem.GetPointIds().SetId(4, nid_map[node_ids[4]])
                    elem.GetPointIds().SetId(5, nid_map[node_ids[5]])
                    elem.GetPointIds().SetId(6, nid_map[node_ids[6]])
                    elem.GetPointIds().SetId(7, nid_map[node_ids[7]])
                else:
                    elem = vtkQuad()
                    elem.GetPointIds().SetId(0, n1)
                    elem.GetPointIds().SetId(1, n2)
                    elem.GetPointIds().SetId(2, n3)
                    elem.GetPointIds().SetId(3, n4)
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element, (CQUAD, CQUADX)):
                # CQUAD, CQUADX are 9 noded quads
                material_coord[i] = 0 if isinstance(element.theta_mcid, float) else element.theta_mcid

                node_ids = element.node_ids
                pid = element.Pid()
                for nid in node_ids:
                    if nid is not None:
                        nid_to_pid_map[nid].append(pid)
                self.eid_to_nid_map[eid] = node_ids[:4]

                n1, n2, n3, n4 = [nid_map[nid] for nid in node_ids[:4]]
                p1 = xyz_cid0[n1, :]
                p2 = xyz_cid0[n2, :]
                p3 = xyz_cid0[n3, :]
                p4 = xyz_cid0[n4, :]
                out = quad_quality(p1, p2, p3, p4)
                (areai, taper_ratioi, area_ratioi, max_skew, aspect_ratio,
                 min_thetai, max_thetai, dideal_thetai, min_edge_lengthi) = out
                if None in node_ids:
                    elem = vtkQuad()
                    elem.GetPointIds().SetId(0, n1)
                    elem.GetPointIds().SetId(1, n2)
                    elem.GetPointIds().SetId(2, n3)
                    elem.GetPointIds().SetId(3, n4)
                else:
                    elem = vtk.vtkBiQuadraticQuad()
                    elem.GetPointIds().SetId(4, nid_map[node_ids[4]])
                    elem.GetPointIds().SetId(5, nid_map[node_ids[5]])
                    elem.GetPointIds().SetId(6, nid_map[node_ids[6]])
                    elem.GetPointIds().SetId(7, nid_map[node_ids[7]])
                    elem.GetPointIds().SetId(8, nid_map[node_ids[8]])
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
                min_thetai, max_thetai, dideal_thetai, min_edge_lengthi = get_min_max_theta(
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
                min_thetai, max_thetai, dideal_thetai, min_edge_lengthi = get_min_max_theta(
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
                min_thetai, max_thetai, dideal_thetai, min_edge_lengthi = get_min_max_theta(
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
                min_thetai, max_thetai, dideal_thetai, min_edge_lengthi = get_min_max_theta(
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
                min_thetai, max_thetai, dideal_thetai, min_edge_lengthi = get_min_max_theta(
                    _chexa_faces, node_ids[:8], nid_map, xyz_cid0)
            elif isinstance(element, (CHEXA20, CIHEX2)):
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

                    # these two blocks are flipped
                    elem.GetPointIds().SetId(12, nid_map[node_ids[16]])
                    elem.GetPointIds().SetId(13, nid_map[node_ids[17]])
                    elem.GetPointIds().SetId(14, nid_map[node_ids[18]])
                    elem.GetPointIds().SetId(15, nid_map[node_ids[19]])

                    elem.GetPointIds().SetId(16, nid_map[node_ids[12]])
                    elem.GetPointIds().SetId(17, nid_map[node_ids[13]])
                    elem.GetPointIds().SetId(18, nid_map[node_ids[14]])
                    elem.GetPointIds().SetId(19, nid_map[node_ids[15]])
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
                min_thetai, max_thetai, dideal_thetai, min_edge_lengthi = get_min_max_theta(
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
                # etype = 14
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                min_thetai, max_thetai, dideal_thetai, min_edge_lengthi = get_min_max_theta(
                    _cpyram_faces, node_ids[:5], nid_map, xyz_cid0)
            elif isinstance(element, CPYRAM13):
                node_ids = element.node_ids
                pid = element.Pid()
                #if None not in node_ids:
                    #print(' node_ids =', node_ids)
                    #elem = vtkQuadraticPyramid()
                    # etype = 27
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
                min_thetai, max_thetai, dideal_thetai, min_edge_lengthi = get_min_max_theta(
                    _cpyram_faces, node_ids[:5], nid_map, xyz_cid0)

            #elif (isinstance(element, (LineElement, SpringElement)) or
                  #etype in ['CBUSH', 'CBUSH1D', 'CFAST', 'CROD', 'CONROD',
                                   #'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
                                   #'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
                                   #'CVISC', 'CGAP']):
            elif etype in ['CBUSH', 'CBUSH1D', 'CFAST',
                           'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
                           'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
                           'CVISC', 'CGAP']:

                # TODO: verify
                # CBUSH, CBUSH1D, CFAST, CELAS1, CELAS3
                # CDAMP1, CDAMP2, CDAMP3, CDAMP4, CDAMP5, CVISC
                if hasattr(element, 'pid'):
                    pid = element.Pid()
                else:
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

                    if 1:
                        elem = vtk.vtkVertex()
                        elem.GetPointIds().SetId(0, j)
                    else:
                        elem = vtk.vtkSphere()
                        #elem = vtk.vtkSphereSource()
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

            elif etype in ('CBAR', 'CBEAM', 'CROD', 'CONROD', 'CTUBE'):
                if etype == 'CONROD':
                    pid = 0
                    areai = element.Area()
                else:
                    pid = element.Pid()
                    try:
                        areai = element.pid_ref.Area()
                    except:
                        print(element)
                        raise

                node_ids = element.node_ids
                for nid in node_ids:
                    nid_to_pid_map[nid].append(pid)

                # 2 points
                min_edge_lengthi = norm(element.nodes[0].get_position() -
                                        element.nodes[1].get_position())
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
            elif etype == 'CHBDYG':
                node_ids = element.node_ids
                pid = element.Pid()
                for nid in node_ids:
                    if nid is not None:
                        nid_to_pid_map[nid].append(pid)

                if element.Type in ['AREA4', 'AREA8']:
                    self.eid_to_nid_map[eid] = node_ids[:4]

                    n1, n2, n3, n4 = [nid_map[nid] for nid in node_ids[:4]]
                    p1 = xyz_cid0[n1, :]
                    p2 = xyz_cid0[n2, :]
                    p3 = xyz_cid0[n3, :]
                    p4 = xyz_cid0[n4, :]
                    out = quad_quality(p1, p2, p3, p4)
                    (areai, taper_ratioi, area_ratioi, max_skew, aspect_ratio,
                     min_thetai, max_thetai, dideal_thetai, min_edge_lengthi) = out
                    if element.Type == 'AREA4' or None in node_ids:
                        elem = vtkQuad()
                    else:
                        elem = vtkQuadraticQuad()
                        elem.GetPointIds().SetId(4, nid_map[node_ids[4]])
                        elem.GetPointIds().SetId(5, nid_map[node_ids[5]])
                        elem.GetPointIds().SetId(6, nid_map[node_ids[6]])
                        elem.GetPointIds().SetId(7, nid_map[node_ids[7]])

                    elem.GetPointIds().SetId(0, n1)
                    elem.GetPointIds().SetId(1, n2)
                    elem.GetPointIds().SetId(2, n3)
                    elem.GetPointIds().SetId(3, n4)
                    self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                elif element.Type in ['AREA3', 'AREA6']:
                    self.eid_to_nid_map[eid] = node_ids[:3]
                    if element.Type == 'AREA3' or None in node_ids:
                        elem = vtkTriangle()
                    else:
                        elem = vtkQuadraticTriangle()
                        elem.GetPointIds().SetId(3, nid_map[node_ids[3]])
                        elem.GetPointIds().SetId(4, nid_map[node_ids[4]])
                        elem.GetPointIds().SetId(5, nid_map[node_ids[5]])

                    n1, n2, n3 = [nid_map[nid] for nid in node_ids[:3]]
                    p1 = xyz_cid0[n1, :]
                    p2 = xyz_cid0[n2, :]
                    p3 = xyz_cid0[n3, :]
                    out = tri_quality(p1, p2, p3)
                    (areai, max_skew, aspect_ratio,
                     min_thetai, max_thetai, dideal_thetai, min_edge_lengthi) = out
                    elem.GetPointIds().SetId(0, n1)
                    elem.GetPointIds().SetId(1, n2)
                    elem.GetPointIds().SetId(2, n3)
                    self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                else:
                    #print('removing\n%s' % (element))
                    print('removing eid=%s; %s' % (eid, element.type))
                    del self.eid_map[eid]
                    self.log_info("skipping %s" % element.type)
                    continue
            else:
                print('removing\n%s' % (element))
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
            #print(eid, min_thetai, max_thetai, '\n', element)
            #asdf
            min_interior_angle[i] = min_thetai
            max_interior_angle[i] = max_thetai
            dideal_theta[i] = dideal_thetai
            max_skew_angle[i] = max_skew
            max_warp_angle[i] = max_warp
            max_aspect_ratio[i] = aspect_ratio
            area[i] = areai
            area_ratio[i] = area_ratioi
            taper_ratio[i] = taper_ratioi
            min_edge_length[i] = min_edge_lengthi
            i += 1
        #assert len(self.eid_map) > 0, self.eid_map
        #print('mapped elements')

        nelements = i
        self.nElements = nelements
        #print('nelements=%s pids=%s' % (nelements, list(pids)))
        pids = pids[:nelements]

        out = (
            nid_to_pid_map, xyz_cid0, pids, nelements, material_coord,
            area, min_interior_angle, max_interior_angle, max_aspect_ratio,
            max_skew_angle, taper_ratio, dideal_theta,
            area_ratio, min_edge_length, max_warp_angle,
        )
        return out

    def _build_normals_quality(self, model, nelements, cases, form0, icase,
                               xyz_cid0, material_coord,
                               min_interior_angle, max_interior_angle, dideal_theta,
                               area, max_skew_angle, taper_ratio,
                               max_warp_angle, area_ratio, min_edge_length, max_aspect_ratio):
        """
        creates:
         - ElementDim
         - Normal X/Y/Z
         - NNodes/Elem
         - Area
         - Min/Max Interior Angle
         - Skew Angle
         - Taper Ratio
         - Area Ratio
         - MaterialCoord
        """
        #ielement = 0
        nelements = self.element_ids.shape[0]
        normals = np.zeros((nelements, 3), dtype='float32')
        offset = np.full(nelements, np.nan, dtype='float32')
        xoffset = np.full(nelements, np.nan, dtype='float32')
        yoffset = np.full(nelements, np.nan, dtype='float32')
        zoffset = np.full(nelements, np.nan, dtype='float32')
        element_dim = np.full(nelements, -1, dtype='int32')
        nnodes_array = np.full(nelements, np.nan, dtype='int32')
        for eid, element in sorted(iteritems(model.elements)):
            etype = element.type
            if isinstance(element, ShellElement):
                ie = None
                element_dimi = 2
                try:
                    normali = element.Normal()
                #except AttributeError:
                    #raise
                except RuntimeError:
                    # this happens when you have a degenerate tri
                    msg = 'eid=%i normal=nan...setting to [2, 2, 2]\n'
                    msg += '%s' % (element)
                    for node in element.nodes:
                        msg += str(node)
                    self.log.error(msg)
                    normali = np.ones(3) * 2.
                    #raise

                prop = element.pid
                ptype = prop.type
                if ptype == 'PSHELL':
                    z0 = prop.z1
                elif ptype in ['PCOMP', 'PCOMPG']:
                    z0 = prop.z0
                elif ptype == 'PLPLANE':
                    z0 = 0.
                elif ptype == 'PSHEAR':
                    z0 = 0.
                elif ptype in ['PSOLID', 'PLSOLID']:
                    z0 = 0.
                else:
                    raise NotImplementedError(ptype) # PSHEAR, PCOMPG

                if z0 is None:
                    if etype in ['CTRIA3', 'CTRIAR']:
                        #node_ids = self.nodes[3:]
                        z0 = (element.T1 + element.T2 + element.T3) / 3.
                        nnodesi = 3
                    elif etype == 'CTRIA6':
                        #node_ids = self.nodes[3:]
                        z0 = (element.T1 + element.T2 + element.T3) / 3.
                        nnodesi = 6
                    elif etype in ['CQUAD4', 'CQUADR']:
                        #node_ids = self.nodes[4:]
                        z0 = (element.T1 + element.T2 + element.T3 + element.T4) / 4.
                        nnodesi = 4
                    elif etype == 'CQUAD8':
                        #node_ids = self.nodes[4:]
                        z0 = (element.T1 + element.T2 + element.T3 + element.T4) / 4.
                        nnodesi = 8
                    elif etype == 'CQUAD':
                        #node_ids = self.nodes[4:]
                        z0 = (element.T1 + element.T2 + element.T3 + element.T4) / 4.
                        nnodesi = 9

                    # axisymmetric
                    elif etype == 'CTRAX3':
                        #node_ids = self.nodes[3:]
                        nnodesi = 3
                        z0 = 0.
                    elif etype == 'CTRAX6':
                        #node_ids = self.nodes[3:]
                        nnodesi = 6
                        z0 = 0.
                    elif etype in ['CTRIAX', 'CTRIAX6']:
                        # the CTRIAX6 uses a non-standard node orientation
                        #node_ids = self.nodes[3:]
                        z0 = 0.
                        nnodesi = 6
                    elif etype == 'CQUADX':
                        #node_ids = self.nodes[4:]
                        nnodesi = 9
                        z0 = 0.
                    elif etype == 'CQUADX4':
                        #node_ids = self.nodes[4:]
                        nnodesi = 4
                        z0 = 0.
                    elif etype == 'CQUADX8':
                        #node_ids = self.nodes[4:]
                        nnodesi = 8
                        z0 = 0.
                    else:
                        raise NotImplementedError(element)
                else:
                    if etype in ['CTRIA3', 'CTRIAR', 'CTRAX3', 'CPLSTN3']:
                        nnodesi = 3
                    elif etype in ['CTRIA6', 'CTRIAX', 'CTRIAX6', 'CPLSTN6', 'CTRAX6']:
                        # no a CTRIAX really has 6 nodes because reasons...
                        nnodesi = 6

                    elif etype in ['CQUAD4', 'CQUADR', 'CPLSTN4', 'CSHEAR', 'CQUADX4']:
                        nnodesi = 4
                    elif etype in ['CQUAD8', 'CPLSTN8', 'CQUADX8']:
                        nnodesi = 8
                    elif etype in ['CQUAD', 'CQUADX']:
                        nnodesi = 9
                    else:
                        raise NotImplementedError(element)

                ie = self.eid_map[eid]
                normals[ie, :] = normali
                if element.type in ['CPLSTN3', 'CPLSTN4', 'CPLSTN6', 'CPLSTN8']:
                    element_dim[ie] = element_dimi
                    nnodes_array[ie] = nnodesi
                    self.log.debug('continue...element.type=%r' % element.type)
                    continue

                offset[ie] = z0
                xoffset[ie] = z0 * normali[0]
                yoffset[ie] = z0 * normali[1]
                zoffset[ie] = z0 * normali[2]

            elif etype == 'CTETRA':
                ie = self.eid_map[eid]
                element_dimi = 3
                nnodesi = 4
            elif etype == 'CPENTA':
                ie = self.eid_map[eid]
                element_dimi = 3
                nnodesi = 6
            elif etype == 'CPYRAM':
                ie = self.eid_map[eid]
                element_dimi = 3
                nnodesi = 5
            elif etype in ['CHEXA', 'CIHEX1', 'CIHEX2']:
                ie = self.eid_map[eid]
                element_dimi = 3
                nnodesi = 8

            elif etype in ['CROD', 'CONROD', 'CBEND', 'CBAR', 'CBEAM', 'CGAP', 'CTUBE']:
                ie = self.eid_map[eid]
                element_dimi = 1
                nnodesi = 2
            elif etype in ['CBUSH', 'CBUSH1D', 'CBUSH2D',
                           'CFAST', 'CVISC',
                           'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
                           'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5']:
                ie = self.eid_map[eid]
                element_dimi = 0
                nnodesi = 2
            elif etype == 'CHBDYG':
                ie = self.eid_map[eid]
                if element.Type == 'AREA3':
                    nnodesi = 3
                    element_dimi = 2
                elif element.Type == 'AREA4':
                    nnodesi = 4
                    element_dimi = 2
                elif element.Type == 'AREA6':
                    nnodesi = 6
                    element_dimi = 2
                elif element.Type == 'AREA8':
                    nnodesi = 8
                    element_dimi = 2
                #elif element.Type == 'REV':
                    #nnodesi = 2 # ???
                    #element_dimi = 1 # ???
                else:
                    element_dimi = -1
                    nnodesi = -1
                    print('element.type=%s doesnt have a dimension' % element.type)

            else:
                ie = self.eid_map[eid]
                element_dimi = -1
                nnodesi = -1
                print('element.type=%s doesnt have a dimension' % element.type)
            assert ie is not None
            element_dim[ie] = element_dimi
            nnodes_array[ie] = nnodesi
            #ielement += 1

        # if not a flat plate
        #if min(nxs) == max(nxs) and min(nxs) != 0.0:
        is_element_dim = np.max(element_dim) != np.min(element_dim)
        is_element_dim = True
        if is_element_dim:
            eid_dim_res = GuiResult(0, header='ElementDim', title='ElementDim',
                                    location='centroid', scalar=element_dim, mask_value=-1)
            cases[icase] = (eid_dim_res, (0, 'ElementDim'))

        is_shell = np.abs(normals).max() > 0.
        is_solid = np.any(np.isfinite(max_interior_angle)) and np.nanmax(np.abs(max_interior_angle)) > 0.
        #print('is_shell=%s is_solid=%s' % (is_shell, is_solid))
        if is_shell:
            nx_res = GuiResult(
                0, header='NormalX', title='NormalX',
                location='centroid', scalar=normals[:, 0], data_format='%.2f')
            ny_res = GuiResult(
                0, header='NormalY', title='NormalY',
                location='centroid', scalar=normals[:, 1], data_format='%.2f')
            nz_res = GuiResult(
                0, header='NormalZ', title='NormalZ',
                location='centroid', scalar=normals[:, 2], data_format='%.2f')
            nxyz_res = NormalResult(0, 'Normals', 'Normals',
                                    nlabels=2, labelsize=5, ncolors=2,
                                    colormap='jet', data_format='%.1f',
                                    uname='NormalResult')
            # this is just for testing nan colors that doesn't work
            #max_interior_angle[:1000] = np.nan
            area_res = GuiResult(0, header='Area', title='Area',
                                 location='centroid', scalar=area)
            min_edge_length_res = GuiResult(
                0, header='Min Edge Length', title='Min Edge Length',
                location='centroid', scalar=min_edge_length)

            min_theta_res = GuiResult(
                0, header='Min Interior Angle', title='Min Interior Angle',
                location='centroid', scalar=np.degrees(min_interior_angle))
            max_theta_res = GuiResult(
                0, header='Max Interior Angle', title='Max Interior Angle',
                location='centroid', scalar=np.degrees(max_interior_angle))
            dideal_theta_res = GuiResult(
                0, header='Delta Ideal Angle', title='Delta Ideal Angle',
                location='centroid', scalar=np.degrees(dideal_theta))

            skew = np.degrees(max_skew_angle)
            skew_res = GuiResult(
                0, header='Max Skew Angle', title='MaxSkewAngle',
                location='centroid', scalar=skew)
            aspect_res = GuiResult(
                0, header='Aspect Ratio', title='AspectRatio',
                location='centroid', scalar=max_aspect_ratio)

            form_checks = []
            form0.append(('Element Checks', None, form_checks))
            if is_element_dim:
                form_checks.append(('ElementDim', icase, []))

            if self.make_nnodes_result:
                nnodes_res = GuiResult(
                    0, header='NNodes/Elem', title='NNodes/Elem',
                    location='centroid', scalar=nnodes_array)
                form_checks.append(('NNodes', icase + 1, []))
                cases[icase + 1] = (nnodes_res, (0, 'NNodes'))
                icase += 1

            cases[icase + 1] = (nx_res, (0, 'NormalX'))
            cases[icase + 2] = (ny_res, (0, 'NormalY'))
            cases[icase + 3] = (nz_res, (0, 'NormalZ'))
            cases[icase + 4] = (nxyz_res, (0, 'Normal'))
            cases[icase + 5] = (area_res, (0, 'Area'))
            cases[icase + 6] = (min_edge_length_res, (0, 'Min Edge Length'))
            cases[icase + 7] = (min_theta_res, (0, 'Min Interior Angle'))
            cases[icase + 8] = (max_theta_res, (0, 'Max Interior Angle'))
            cases[icase + 9] = (dideal_theta_res, (0, 'Delta Ideal Angle'))
            cases[icase + 10] = (skew_res, (0, 'Max Skew Angle'))
            cases[icase + 11] = (aspect_res, (0, 'Aspect Ratio'))

            form_checks.append(('NormalX', icase + 1, []))
            form_checks.append(('NormalY', icase + 2, []))
            form_checks.append(('NormalZ', icase + 3, []))
            form_checks.append(('Normal', icase + 4, []))
            form_checks.append(('Area', icase + 5, []))
            form_checks.append(('Min Edge Length', icase + 6, []))
            form_checks.append(('Min Interior Angle', icase + 7, []))
            form_checks.append(('Max Interior Angle', icase + 8, []))
            form_checks.append(('Delta Ideal Angle', icase + 9, []))
            form_checks.append(('Max Skew Angle', icase + 10, []))
            form_checks.append(('Aspect Ratio', icase + 11, []))
            icase += 12

            if np.nanmax(area_ratio) > 1.:
                arearatio_res = GuiResult(
                    0, header='Area Ratio', title='Area Ratio',
                    location='centroid', scalar=area_ratio)
                cases[icase] = (arearatio_res, (0, 'Area Ratio'))
                form_checks.append(('Area Ratio', icase, []))
                icase += 1

            if np.nanmax(taper_ratio) > 1.:
                taperratio_res = GuiResult(
                    0, header='Taper Ratio', title='Taper Ratio',
                    location='centroid', scalar=taper_ratio)
                cases[icase] = (taperratio_res, (0, 'Taper Ratio'))
                form_checks.append(('Taper Ratio', icase, []))
                icase += 1

            if np.nanmax(max_warp_angle) > 0.0:
                warp_res = GuiResult(
                    0, header='Max Warp Angle', title='MaxWarpAngle',
                    location='centroid', scalar=np.degrees(max_warp_angle))
                cases[icase + 4] = (warp_res, (0, 'Max Warp Angle'))
                form_checks.append(('Max Warp Angle', icase, []))
                icase += 1

            #if (np.abs(xoffset).max() > 0.0 or np.abs(yoffset).max() > 0.0 or
                #np.abs(zoffset).max() > 0.0):
            # offsets
            offset_res = GuiResult(
                0, header='Offset', title='Offset',
                location='centroid', scalar=offset, data_format='%g')
            offset_x_res = GuiResult(
                0, header='OffsetX', title='OffsetX',
                location='centroid', scalar=xoffset, data_format='%g')
            offset_y_res = GuiResult(
                0, header='OffsetY', title='OffsetY',
                location='centroid', scalar=yoffset, data_format='%g')
            offset_z_res = GuiResult(
                0, header='OffsetZ', title='OffsetZ',
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
                x_res = GuiResult(
                    0, header='X', title='X',
                    location='node', scalar=xyz_cid0[:, 0], data_format='%g')
                y_res = GuiResult(
                    0, header='Y', title='Y',
                    location='node', scalar=xyz_cid0[:, 1], data_format='%g')
                z_res = GuiResult(
                    0, header='Z', title='Z',
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

            min_edge_length_res = GuiResult(
                0, header='Min Edge Length', title='Min Edge Length',
                location='centroid', scalar=min_edge_length)
            min_theta_res = GuiResult(
                0, header='Min Interior Angle', title='Min Interior Angle',
                location='centroid', scalar=np.degrees(min_interior_angle))
            max_theta_res = GuiResult(
                0, header='Max Interior Angle', title='Max Interior Angle',
                location='centroid', scalar=np.degrees(max_interior_angle))
            #skew = 90. - np.degrees(max_skew_angle)
            #skew_res = GuiResult(0, header='Max Skew Angle', title='MaxSkewAngle',
                                    #location='centroid', scalar=skew)
            if is_element_dim:
                form_checks.append(('ElementDim', icase, []))
            form_checks.append(('Min Edge Length', icase + 1, []))
            form_checks.append(('Min Interior Angle', icase + 2, []))
            form_checks.append(('Max Interior Angle', icase + 3, []))
            #form_checks.append(('Max Skew Angle', icase + 4, []))
            cases[icase + 1] = (min_edge_length_res, (0, 'Min Edge Length'))
            cases[icase + 2] = (min_theta_res, (0, 'Min Interior Angle'))
            cases[icase + 3] = (max_theta_res, (0, 'Max Interior Angle'))
            #cases[icase + 4] = (skew_res, (0, 'Max Skew Angle'))
            icase += 4

        else:
            form0.append(('ElementDim', icase, []))
            icase += 1

        if np.abs(material_coord).max() > 0:
            material_coord_res = GuiResult(
                0, header='MaterialCoord', title='MaterialCoord',
                location='centroid',
                scalar=material_coord, data_format='%i')
            cases[icase] = (material_coord_res, (0, 'MaterialCoord'))
            form0.append(('MaterialCoord', icase, []))
            icase += 1
        return icase, normals

    def _build_properties(self, model, nelements, eids, pids,
                          cases, form0, icase):
        """
        creates:
          - PropertyID

        TODO: CONROD
        """
        prop_types_with_mid = [
            'PSOLID', 'PSHEAR',
            'PROD', 'PTUBE', 'PBAR', 'PBARL', 'PBEAM', 'PBEAML',
        ]

        upids = None
        pid_res = GuiResult(0, header='PropertyID', title='PropertyID',
                            location='centroid', scalar=pids)
        cases[icase] = (pid_res, (0, 'PropertyID'))
        form0.append(('PropertyID', icase, []))
        icase += 1

        upids = np.unique(pids)
        mid_eids_skip = []

        pcomp_nplies = 0
        nplies = 1
        is_pshell = False
        is_pcomp = False
        if 'PSHELL' in model.card_count:
            nplies = 4
            is_pshell = True
        for pid in model.get_card_ids_by_card_types(['PCOMP', 'PCOMPG'], combine=True):
            prop = model.properties[pid]
            pcomp_nplies = max(pcomp_nplies, prop.nplies)
            is_pcomp = True
        nplies = max(nplies, pcomp_nplies + 1)

        mids = np.zeros((nelements, nplies), dtype='int32')
        thickness = np.full((nelements, nplies), np.nan, dtype='float32')
        #rho = np.full((nelements, nplies), np.nan, dtype='float32')
        nplies = np.zeros(nelements, dtype='int32')
        for pid in upids:
            if pid == 0:
                print('skipping pid=0')
                continue
            prop = model.properties[pid]
            #try:
            if prop.type in prop_types_with_mid:
                # simple types
                i = np.where(pids == pid)[0]
                mid = prop.mid_ref.mid
                mids[i, 0] = mid
            elif prop.type == 'PSHELL':
                # TODO: only considers mid1
                i = np.where(pids == pid)[0]
                mid1 = prop.Mid1()
                mid2 = prop.Mid2()
                mid3 = prop.Mid3()
                mid4 = prop.Mid4()
                mids[i, 0] = mid1 if mid1 is not None else 0
                mids[i, 1] = mid2 if mid2 is not None else 0
                mids[i, 2] = mid3 if mid3 is not None else 0
                mids[i, 3] = mid4 if mid4 is not None else 0
                thickness[i, 0] = prop.Thickness()
            elif prop.type in ['PCOMP', 'PCOMPG']:
                # TODO: only considers iply=0
                i = np.where(pids == pid)[0]
                npliesi = prop.nplies
                nplies[i] = npliesi
                for iply in range(npliesi):
                    mids[i, iply+1] = prop.Mid(iply)
                    thickness[i, iply+1] = prop.Thickness(iply)
                thickness[i, 0] = thickness[i[0], 1:].sum()

                #mids[i, 0] = mids[i, 1]
            elif prop.type in ['PELAS', 'PBUSH', 'PDAMP', 'PDAMPT']:
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
        return icase, upids, mids, thickness, nplies, (is_pshell, is_pcomp)

    def _build_materials(self, model, mids, thickness, nplies, is_pshell_pcomp,
                         cases, form0, icase):
        """
        creates:
          - Thickness
          - nPlies (composite only)
          - Material ID
          - E_11
          - E_22
          - E_33
          - Is Isotropic?
        """
        nlayers = mids.shape[1]

        if nplies.max() > 0:
            nplies_res = GuiResult(0, header='Number of Plies', title='nPlies',
                              location='centroid', scalar=nplies)
            cases[icase] = (nplies_res, (0, 'Number of Plies'))
            form0.append(('Number of Plies', icase, []))
            icase += 1

        for ilayer in range(nlayers):
            midsi = mids[:, ilayer]
            if midsi.max() == 0:
                print('cant find anything in ilayer=%s' % ilayer)
                continue
            thicknessi = thickness[:, ilayer]

            form_layer = []
            has_mat8, has_mat9, e11, e22, e33 = self._get_material_arrays(model, midsi)
            if np.any(np.isfinite(thicknessi)) and np.nanmax(thicknessi) > 0.0:
                t_res = GuiResult(0, header='Thickness', title='Thickness',
                                  location='centroid', scalar=thicknessi)
                cases[icase] = (t_res, (0, 'Thickness'))
                form_layer.append(('Thickness', icase, []))
                icase += 1

            mid_res = GuiResult(0, header='MaterialID', title='MaterialID',
                                location='centroid', scalar=midsi, mask_value=0)
            cases[icase] = (mid_res, (0, 'MaterialID'))
            form_layer.append(('MaterialID', icase, []))
            icase += 1

            if has_mat9: # also implicitly has_mat8
                is_orthotropic = not (np.array_equal(e11, e22) and np.array_equal(e11, e33))
            elif has_mat8:
                is_orthotropic = not np.array_equal(e11, e22)
            else:
                is_orthotropic = False

            if is_orthotropic:
                e11_res = GuiResult(0, header='E_11', title='E_11',
                                    location='centroid', scalar=e11, data_format='%.3e')
                e22_res = GuiResult(0, header='E_22', title='E_22',
                                    location='centroid', scalar=e22, data_format='%.3e')
                cases[icase] = (e11_res, (0, 'E_11'))
                cases[icase + 1] = (e22_res, (0, 'E_22'))
                form_layer.append(('E_11', icase, []))
                form_layer.append(('E_22', icase + 1, []))
                icase += 2

                is_isotropic = np.zeros(len(e11), dtype='int8')
                if has_mat9:
                    is_isotropic[(e11 == e22) | (e11 == e33)] = 1
                    e33_res = GuiResult(0, header='E_33', title='E_33',
                                        location='centroid', scalar=e33, data_format='%.3e')
                    cases[icase] = (e33_res, (0, 'E_33'))
                    form_layer.append(('E_33', icase, []))
                    icase += 1
                else:
                    #is_isotropic_map = e11 == e22
                    is_isotropic[e11 == e22] = 1

                iso_res = GuiResult(0, header='IsIsotropic?', title='IsIsotropic?',
                                    location='centroid', scalar=is_isotropic, data_format='%i')
                cases[icase] = (iso_res, (0, 'Is Isotropic?'))
                form_layer.append(('Is Isotropic?', icase, []))
                icase += 1
            elif np.nanmax(e11) > 0.:
                # isotropic
                e11_res = GuiResult(0, header='E', title='E',
                                    location='centroid', scalar=e11, data_format='%.3e')
                cases[icase] = (e11_res, (0, 'E'))
                form_layer.append(('E', icase, []))
                icase += 1

            if nlayers == 1:
                form0 += form_layer
            else:
                word = self._get_nastran_gui_layer_word(ilayer, is_pshell_pcomp)
                form0.append((word, None, form_layer))
        return icase

    def _get_nastran_gui_layer_word(self, ilayer, is_pshell_pcomp):
        """gets the PSHELL/PCOMP layer word"""
        is_pshell, is_pcomp = is_pshell_pcomp
        word = ''
        if ilayer == 0:
            if is_pshell:
                word += 'PSHELL: mid%i; ' % (ilayer + 1)
            if is_pcomp:
                word += 'PCOMP: Total; '
            word = word.rstrip('; ') + ' & others'

        elif ilayer in [1, 2, 3]:
            if is_pshell:
                word += 'PSHELL: mid%i; ' % (ilayer + 1)
            if is_pcomp:
                word += 'PCOMP: ilayer=%i' % (ilayer)
            word = word.rstrip('; ')
        else:
            assert is_pcomp, ilayer
            word += 'PCOMP: ilayer=%i' % (ilayer)
        return word

    def _build_optimization(self, model, pids, upids, nelements, cases, form0, icase):
        """
        Creates the optimization visualization.  Supports:
          - DVPREL1/2 shell thickness:
             - DV Region
             - DVPREL Init - t
             - DVPREL Min - t
             - DVPREL Max - ts
        """
        if upids is None:
            return icase
        if len(model.properties) and len(model.dvprels):
            # len(model.dvprels) + len(model.dvcrels) + len(model.dvmrels) + len(model.desvars)
            #dvmrel_init = np.zeros(nelements, dtype='int32')
            #dvgrel_init = np.zeros(nelements, dtype='int32')
            out = model._get_dvprel_ndarrays(nelements, pids)
            dvprel_t_init, dvprel_t_min, dvprel_t_max, design_region = out

            region_res = GuiResult(
                0, header='DV Region', title='DV Region',
                location='centroid', scalar=design_region)
            t_init_res = GuiResult(
                0, header='DVPREL Init - t', title='DVPREL Init - t',
                location='centroid', scalar=dvprel_t_init)
            t_min_res = GuiResult(
                0, header='DVPREL Min - t', title='DVPREL Min - t',
                location='centroid', scalar=dvprel_t_min)
            t_max_res = GuiResult(
                0, header='DVPREL Max - t', title='DVPREL Max - t',
                location='centroid', scalar=dvprel_t_max)
            cases[icase] = (region_res, (0, 'DV Region'))
            cases[icase + 1] = (t_init_res, (0, 'DVPREL Init - t'))
            cases[icase + 2] = (t_min_res, (0, 'DVPREL Min - t'))
            cases[icase + 3] = (t_max_res, (0, 'DVPREL Max - t'))
            opt = []
            opt.append(('DV Region', icase, []))
            opt.append(('DVPREL Init - t', icase + 1, []))
            opt.append(('DVPREL Min - t', icase + 2, []))
            opt.append(('DVPREL Max - t', icase + 3, []))
            form0.append(('Optimization', '', opt))
            icase += 4
        return icase

    def _plot_pressures(self, model, cases, form0, icase, subcase_id):
        """
        pressure act normal to a shell (as opposed to anti-normal to a solid face)
        """
        # quit out if we're going to make pressure plots anyways
        if self.plot_applied_loads:
            return icase

        # quit out if we don't have pressures
        if 'PLOAD4' not in model.card_count:
            return icase

        subcase = model.subcases[subcase_id]

        try:
            load_case_id = subcase.get_parameter('LOAD')[0]
        except KeyError:
            return icase

        try:
            load_case = model.loads[load_case_id]
        except KeyError:
            self.log.warning('LOAD=%s not found' % load_case_id)
            return icase

        is_pressure, pressures = model.get_pressure_array(
            load_case, eids=self.element_ids, normals=self.normals)
        if not is_pressure:
            return icase

        # if there is no applied pressure, don't make a plot
        if np.abs(pressures).max():
            case_name = 'Pressure'
            # print('iload=%s' % iload)
            # print(case_name)
            pressure_res = GuiResult(
                subcase_id, header='Pressure', title='Pressure',
                location='centroid', scalar=pressures)
            cases[icase] = (pressure_res, (0, 'Pressure'))
            form0.append((case_name, icase, []))
            icase += 1
        return icase

    def _plot_applied_loads(self, model, cases, form0, icase, subcase_id,
                            xref_loads=True):
        """
        Applied loads include:
        ----------------------
         - Centroidal Pressure
         - Fx, Fy, Fz
         - SPCDx, SPCDy, SPCDz, SPCDxyz
         - Temperature(MATERIAL)
         - Temperature(INITIAL)
         - Temperature(LOAD)
         - Temperature(BOTH)
        """
        if not xref_loads:
            print('returning from plot_applied_loads_early')
            return icase

        form = []
        out = model.get_load_arrays(
            subcase_id, nid_map=self.nid_map,
            eid_map=self.eid_map, node_ids=self.node_ids,
            normals=self.normals)
        is_loads, is_temperatures, temperature_data, load_data = out

        if is_loads:
            centroidal_pressures, forces, spcd = load_data
            if np.abs(centroidal_pressures).max():
                pressure_res = GuiResult(subcase_id, header='Pressure', title='Pressure',
                                         location='centroid', scalar=centroidal_pressures)
                cases[icase] = (pressure_res, (0, 'Pressure'))
                form0.append(('Pressure', icase, []))
                icase += 1

            if np.abs(forces.max() - forces.min()) > 0.0:
                fxyz = forces[:, :3]
                mxyz = forces[:, 3:]
                fscalar = np.linalg.norm(fxyz, axis=1)
                mscalar = np.linalg.norm(mxyz, axis=1)
                if fscalar.max() > 0:
                    titles = ['Force XYZ']
                    headers = titles
                    assert fxyz.shape[1] == 3, fxyz.shape
                    assert fxyz.shape[0] == len(fscalar)
                    scales = [1.0]

                    force_xyz_res = ForceTableResults(
                        subcase_id, titles, headers, fxyz, fscalar,
                        scales, data_formats=None,
                        nlabels=None, labelsize=None, ncolors=None, colormap='jet',
                        set_max_min=False, uname='NastranGeometry')
                    force_xyz_res.save_defaults()

                    cases[icase] = (force_xyz_res, (0, 'Force XYZ'))
                    form0.append(('Force XYZ', icase, []))
                    icase += 1

                if mscalar.max() > 0:
                    titles = ['Moment XYZ']
                    headers = titles
                    assert mxyz.shape[1] == 3, mxyz.shape
                    assert mxyz.shape[0] == len(mscalar)
                    scales = [1.0]

                    moment_xyz_res = ForceTableResults(
                        subcase_id, titles, headers, mxyz, mscalar,
                        scales, data_formats=None,
                        nlabels=None, labelsize=None, ncolors=None, colormap='jet',
                        set_max_min=False, uname='NastranGeometry')
                    moment_xyz_res.save_defaults()

                    cases[icase] = (moment_xyz_res, (0, 'Moment XYZ'))
                    form0.append(('Moment XYZ', icase, []))
                    icase += 1

            if np.abs(spcd.max() - spcd.min()) > 0.0:
                t123 = spcd[:, :3]
                tnorm = norm(t123, axis=1)
                assert len(tnorm) == len(spcd[:, 2]), len(spcd[:, 2])

                spcd_x_res = GuiResult(subcase_id, header='SPCDx', title='SPCDx',
                                       location='node', scalar=forces[:, 0])
                spcd_y_res = GuiResult(subcase_id, header='SPCDy', title='SPCDy',
                                       location='node', scalar=forces[:, 1])
                spcd_z_res = GuiResult(subcase_id, header='SPCDz', title='SPCDz',
                                       location='node', scalar=forces[:, 2])
                spcd_xyz_res = GuiResult(subcase_id, header='SPCD XYZ', title='SPCD XYZ',
                                         location='node', scalar=tnorm)

                cases[icase] = (spcd_x_res, (0, 'SPCDx'))
                form0.append(('SPCDx', icase, []))
                icase += 1

                cases[icase] = (spcd_y_res, (0, 'SPCDy'))
                form0.append(('SPCDy', icase, []))
                icase += 1

                cases[icase] = (spcd_z_res, (0, 'SPCDz'))
                form0.append(('SPCDz', icase, []))
                icase += 1

                cases[icase] = (spcd_xyz_res, (0, 'SPCD XYZ'))
                form0.append(('SPCD XYZ', icase, []))
                icase += 1

        if is_temperatures:
            temperature_key, temperatures = temperature_data
            temperature_res = GuiResult(subcase_id, header=temperature_key, title=temperature_key,
                                        location='node', scalar=temperatures)
            cases[icase] = (temperature_res, (0, temperature_key))
            form.append((temperature_key, icase, []))
            icase += 1
        return icase

    def load_nastran_results(self, op2_filename, dirname):
        """
        Loads the Nastran results into the GUI
        """
        #gridResult.SetNumberOfComponents(self.nElements)
        self.scalarBar.VisibilityOn()
        self.scalarBar.Modified()

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
            #transforms = self.transforms
        except AttributeError:
            self.log.error('Skipping displacment transformation')
        else:
            model.transform_displacements_to_global(
                i_transform, self.model.coords, xyz_cid0=self.xyz_cid0)

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
                # stress
                icase = self._fill_op2_stress(
                    cases, model, key, icase, itime,
                    stress_dict, header_dict, is_static)

                # strain
                icase = self._fill_op2_strain(
                    cases, model, key, icase, itime,
                    strain_dict, header_dict, is_static)

                # force
                icase = self._fill_op2_force(
                    cases, model, key, icase, itime,
                    force_dict, header_dict, is_static)

                # strain energy
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
                print('header = %r' % header)
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

    def clear_nastran(self):
        """cleans up variables specific to Nastran"""
        self.eid_map = {}
        self.nid_map = {}
        self.eid_to_nid_map = {}
        self.element_ids = None
        self.node_ids = None
