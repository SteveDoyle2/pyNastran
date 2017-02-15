# pylint: disable=E1101
"""
Defines the GUI IO file for Nastran.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import os
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
from vtk.util.numpy_support import numpy_to_vtk

#from pyNastran import is_release
from pyNastran.utils import integer_types
from pyNastran.bdf.bdf import (BDF,
                               CAERO1, CAERO2, CAERO3, CAERO4, CAERO5,
                               CQUAD4, CQUAD8, CQUADR, CSHEAR,
                               CTRIA3, CTRIA6, CTRIAR,
                               CPLSTN3, CPLSTN4, CPLSTN6, CPLSTN8,
                               CTRAX3, CTRIAX6, # CTRIAX, CTRAX6,
                               CQUADX4, CQUADX8, # CQUADX,
                               CONM2,
                               LOAD)

from pyNastran.bdf.cards.elements.shell import ShellElement
#from pyNastran.bdf.cards.elements.bars import LineElement
#from pyNastran.bdf.cards.elements.springs import SpringElement
from pyNastran.bdf.cards.elements.solid import (
    CTETRA4, CTETRA10, CPENTA6, CPENTA15,
    CHEXA8, CHEXA20, CIHEX1, CIHEX2,
    CPYRAM5, CPYRAM13,
)

from pyNastran.gui.gui_objects.gui_result import GuiResult
from pyNastran.converters.nastran.geometry_helper import (
    NastranGeometryHelper, tri_quality, quad_quality, get_min_max_theta)
from pyNastran.converters.nastran.results_helper import NastranGuiResults

from pyNastran.op2.op2 import OP2
#from pyNastran.f06.f06_formatting import get_key0
try:
    from pyNastran.op2.op2_geom import OP2Geom
    is_geom = True
except ImportError:
    is_geom = False

green = (0., 1., 0.)
blue = (0., 0., 1.)
dunno = (0.5, 1., 0.5)
pink = (0.98, 0.4, 0.93)
orange = (219/255., 168/255., 13/255.)
red = (1., 0., 0.)
yellow = (1., 1., 0.)
purple = (1., 0., 1.)


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

    def get_xyz_in_coord(self, model, points, cid=0, dtype='float32'):
        #import time
        #t0 = time.time()

        if 1:
            # t=.578
            out = model.get_displacement_index_xyz_cp_cd(
                dtype='float32')
            icd_transform, icp_transform, xyz_cp, nid_cp_cd = out
            self.i_transform = icd_transform
            xyz_cid0 = model.transform_xyzcp_to_xyz_cid(xyz_cp, icp_transform, cid=0)

            data_type = vtk.VTK_FLOAT
            points_array = numpy_to_vtk(
                num_array=xyz_cid0,
                deep=True,
                array_type=data_type
            )
            points.SetData(points_array)
            nid_map = self.nid_map
            for i, nid in enumerate(nid_cp_cd[:, 0]):
                nid_map[nid] = i
        elif 0:
            # t=.573
            out = model.get_displacement_index_xyz_cp_cd(
                dtype='float32')
            icd_transform, icp_transform, xyz_cp, nid_cp_cd = out
            self.i_transform = icd_transform
            xyz_cid0 = model.transform_xyzcp_to_xyz_cid(xyz_cp, icp_transform, cid=0)

            data_type = vtk.VTK_FLOAT
            points_array = numpy_to_vtk(
                num_array=xyz_cid0,
                deep=False,
                array_type=data_type
            )
            points.SetData(points_array)
            nid_map = self.nid_map
            for i, nid in enumerate(nid_cp_cd[:, 0]):
                nid_map[nid] = i
        else:
            # t=.75
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

            # get indicies and transformations for displacements
            #self.i_transform, self.transforms = model.get_displacement_index_transforms()
            self.i_transform = model.get_displacement_index()

        self._add_nastran_spoints_to_grid(model)
        #print('dt_nastran_xyz =', time.time() - t0)
        return xyz_cid0

    def load_nastran_geometry(self, bdf_filename, dirname, name='main', plot=True):
        self.eid_maps[name] = {}
        self.nid_maps[name] = {}
        self.i_transform = {}
        self.spc_names = []
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

        nnodes = len(model.nodes)
        nspoints = 0
        spoints = None
        if model.spoints:
            spoints = model.spoints.points
            nspoints = len(spoints)

        assert nnodes + nspoints > 0, model.card_count
        nelements = model.nelements
        nplotels = len(model.plotels)
        ncaero_cards = len(model.caeros)
        assert nelements + ncaero_cards + nplotels > 0

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

        nconm2 = 0
        if 'CONM2' in model.card_count:
            nconm2 += model.card_count['CONM2']
        if 'CMASS1' in model.card_count:
            nconm2 += model.card_count['CMASS1']
        if 'CMASS2' in model.card_count:
            nconm2 += model.card_count['CMASS2']

        if nconm2 > 0:
            self.create_alternate_vtk_grid(
                'conm2', color=orange, line_width=5, opacity=1., point_size=4,
                representation='point')

        # Allocate grids
        self.grid.Allocate(self.nElements, 1000)
        if self.has_caero:
            if 'caero' not in self.alt_grids:
                self.create_alternate_vtk_grid(
                    'caero', color=yellow, line_width=3, opacity=1.0,
                    representation='toggle', is_visible=True, is_pickable=False)
            if 'caero_subpanels' not in self.alt_grids:
                self.create_alternate_vtk_grid(
                    'caero_subpanels', color=yellow, line_width=3, opacity=1.0,
                    representation='toggle', is_visible=False, is_pickable=False)

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
                idsi = suport.node_ids
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

        for grid_name in ['suport', 'mpc', 'mpc_dependent', 'mpc_independent'] + self.spc_names:
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

            # rigid body elements and MPCS
            lines = self._get_rigid(model)
            if 'MPC' in subcase:
                mpc_id, options = subcase.get_parameter('MPC')
                if mpc_id is not None:
                    nmpcs = model.card_count['MPC'] if 'MPC' in model.card_count else 0
                    if nmpcs:
                        lines += model.get_MPCx_node_ids_c1(mpc_id, exclude_mpcadd=False)
            self._fill_dependent_independent(dim_max, model, lines, nid_to_pid_map)

            if 'SUPORT1' in subcase.params:  ## TODO: should this be SUPORT?
                suport_id, options = subcase.get_parameter('SUPORT1')
                if 'SUPORT' in model.card_count or 'SUPORT1' in model.card_count:
                    if suport_id:
                        self._fill_suport(suport_id, dim_max, model)

    def _fill_spc(self, spc_id, nspcs, nspc1s, nspcds, dim_max, model, nid_to_pid_map):
        spc_name = 'spc_subcase=%i' % spc_id
        self.spc_names.append(spc_name)
        self.create_alternate_vtk_grid(spc_name, color=purple, line_width=5, opacity=1.,
                                       point_size=5, representation='point', is_visible=False)

        # node_ids = model.get_SPCx_node_ids(spc_id, exclude_spcadd=False)
        node_ids_c1 = model.get_SPCx_node_ids_c1(spc_id, exclude_spcadd=False,
                                                 stop_on_failure=False)

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
        self._add_nastran_nodes_to_grid(spc_name, node_ids, model, nid_to_pid_map)

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
        bar_nids, bar_types, out = self._get_bar_yz_arrays(model, bar_beam_eids, scale, debug)

        bar_nids = list(bar_nids)
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
                icase += 1
            if no_shear_z.max() == 1:
                bar_form[2].append(['No Shear Z', icase, []])
                shear_z_res = GuiResult(0, header='No Shear Z', title='No Shear Z',
                                        location='centroid', scalar=no_shear_z)
                cases[icase] = (shear_z_res, (0, 'No Shear Z'))
                icase += 1
            if no_bending_y.max() == 1:
                bar_form[2].append(['No Bending Y', icase, []])
                bending_y_res = GuiResult(0, header='No Bending Z', title='No Bending Z',
                                          location='centroid', scalar=no_bending_y)
                cases[icase] = (bending_y_res, (0, 'No Bending Z'))
                icase += 1
            if no_bending_z.max() == 1:
                bar_form[2].append(['No Bending Z', icase, []])
                bending_z_res = GuiResult(0, header='No Bending Z', title='No Bending Z',
                                          location='centroid', scalar=no_bending_z)
                cases[icase] = (bending_z_res, (0, 'No Bending Z'))
                icase += 1

        if self.make_released_dofs2 and 0:
            if no_bending.max() == 1:
                bar_form[2].append(['No Bending', icase, []])
                bending_res = GuiResult(0, header='No Bending', title='No Bending',
                                        location='centroid', scalar=no_bending)
                cases[icase] = (bending_res, (0, 'No Bending'))
                icase += 1

            if no_bending_bad.max() == 1:
                bar_form[2].append(['No Bending (Bad)', icase, []])
                type_res = GuiResult(0, header='No Bending (Bad)', title='No Bending (Bad)',
                                     location='centroid', scalar=no_bending_bad)
                cases[icase] = (type_res, (0, 'No Bending (Bad)'))
                icase += 1

            if no_6_16.max() == 1:
                bar_form[2].append(['no_6_16', icase, []])
                type_res = GuiResult(0, header='no_6_16', title='no_6_16',
                                     location='centroid', scalar=no_6_16)
                cases[icase] = (type_res, (0, 'no_6_16'))
                icase += 1
            if no_0_56.max() == 1:
                bar_form[2].append(['no_0_56', icase, []])
                type_res = GuiResult(0, header='no_0_56', title='no_0_56',
                                     location='centroid', scalar=no_0_56)
                cases[icase] = (type_res, (0, 'no_0_56'))
                icase += 1
            if no_0_456.max() == 1:
                bar_form[2].append(['no_0_456', icase, []])
                type_res = GuiResult(0, header='no_0_456', title='no_0_456',
                                     location='centroid', scalar=no_0_456)
                cases[icase] = (type_res, (0, 'no_0_456'))
                icase += 1
            if no_56_456.max() == 1:
                bar_form[2].append(['no_56_456', icase, []])
                type_res = GuiResult(0, header='no_56_456', title='no_56_456',
                                     location='centroid', scalar=no_56_456)
                cases[icase] = (type_res, (0, 'no_56_456'))
                icase += 1
            if no_0_6.max() == 1:
                bar_form[2].append(['no_0_6', icase, []])
                type_res = GuiResult(0, header='no_0_6', title='no_0_6',
                                     location='centroid', scalar=no_0_6)
                cases[icase] = (type_res, (0, 'no_0_6'))
                icase += 1
            if no_0_16.max() == 1:
                bar_form[2].append(['no_0_16)', icase, []])
                type_res = GuiResult(0, header='no_0_16', title='no_0_16',
                                     location='centroid', scalar=no_0_16)
                cases[icase] = (type_res, (0, 'no_0_16'))
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

    def _fill_dependent_independent(self, dim_max, model, lines, nid_to_pid_map):
        if not lines:
            return
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
        self.follower_nodes[name] = node_ids
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

    def _add_nastran_spoints_to_grid(self, model, nid_to_pid_map=None):
        """used to create SPOINTs"""
        if model.spoints is None:
            return
        spoint_ids = list(model.spoints.points) # se t-> list
        assert isinstance(spoint_ids, list), type(spoint_ids)

        nspoints = len(spoint_ids)
        name = 'SPoints'
        if nspoints == 0:
            model.log.warning('0 spoints added for %r' % name)
            return
        self.create_alternate_vtk_grid(
            name, color=blue, line_width=1, opacity=1.,
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

    def _fill_suport(self, suport_id, dim_max, model):
        """creates SUPORT and SUPORT1 nodes"""
        self.create_alternate_vtk_grid(
            'suport', color=red, line_width=5, opacity=1., point_size=4,
            representation='point', is_visible=False)
        suport_nids = self._get_suport_node_ids(model, suport_id)
        self._add_nastran_nodes_to_grid('suport', suport_nids, model)

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

        nplotels = len(model.plotels)
        if nplotels:
            lines = []
            for (eid, element) in sorted(iteritems(model.plotels)):
                node_ids = element.node_ids
                lines.append(node_ids)
            lines = np.array(lines, dtype='int32')

            self.create_alternate_vtk_grid(
                'plotel', color=red, line_width=2, opacity=0.8,
                point_size=5, representation='wire', is_visible=True)
            self._add_nastran_lines_to_grid('plotel', lines, model)

        for (eid, element) in sorted(iteritems(model.elements)):
            self.eid_map[eid] = i
            etype = element.type
            # if element.Pid() >= 82:
                # continue
            # if element.Pid() in pids_to_drop:
                # continue
            # if element.Pid() not in pids_to_keep:
                # continue
            # if element.pid.type == 'PSOLID':
                # continue
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
                areai, max_skew, aspect_ratio, min_thetai, max_thetai, dideal_thetai = tri_quality(p1, p2, p3)

                elem.GetPointIds().SetId(0, n1)
                elem.GetPointIds().SetId(1, n2)
                elem.GetPointIds().SetId(2, n3)
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element, (CTRIA6, CTRIAX6, CPLSTN6)):
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
                areai, max_skew, aspect_ratio, min_thetai, max_thetai, dideal_thetai = tri_quality(p1, p2, p3)
                elem.GetPointIds().SetId(0, n1)
                elem.GetPointIds().SetId(1, n2)
                elem.GetPointIds().SetId(2, n3)
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element, CTRIAX6):
                # the CTRIAX6 is not a standard second-order element
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
                areai, max_skew, aspect_ratio, min_thetai, max_thetai, dideal_thetai = tri_quality(p1, p2, p3)
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
                areai, taper_ratioi, area_ratioi, max_skew, aspect_ratio, min_thetai, max_thetai, dideal_thetai = quad_quality(p1, p2, p3, p4)

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
                areai, taper_ratioi, area_ratioi, max_skew, aspect_ratio, min_thetai, max_thetai, dideal_thetai = quad_quality(p1, p2, p3, p4)
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
                # etype = 14
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                min_thetai, max_thetai, dideal_thetai = get_min_max_theta(
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
                min_thetai, max_thetai, dideal_thetai = get_min_max_theta(
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

            elif etype in ('CBAR', 'CBEAM', 'CROD', 'CONROD'):
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
                lengthi = norm(element.nodes[0].get_position() -
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
            i += 1
        #assert len(self.eid_map) > 0, self.eid_map

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
        del pids_dict


        self.iSubcaseNameMap = {1: ['Nastran', '']}
        icase = 0
        form = ['Geometry', None, []]
        form0 = form[2]

        #new_cases = True
        # set to True to enable node_ids as an result
        nids_set = True
        if nids_set:
            nids = np.zeros(self.nNodes, dtype='int32')
            cds = np.zeros(self.nNodes, dtype='int32')
            for (nid, nid2) in iteritems(self.nid_map):  # map node ids to index
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
            icase, upids, mids, thickness = self._build_properties(
                model, nelements, eids, pids, cases, form0, icase)
            icase = self._build_materials(model, mids, thickness, cases, form0, icase)

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

        #print('nelements=%s eid_map=%s' % (nelements, self.eid_map))
        if self.make_offset_normals_dim and nelements:
            icase, normals = self._build_normals_quality(
                model, nelements, cases, form0, icase,
                xyz_cid0, material_coord,
                min_interior_angle, max_interior_angle, dideal_theta,
                area, max_skew_angle, taper_ratio,
                max_warp_angle, area_ratio, max_aspect_ratio)
            self.normals = normals
        return nid_to_pid_map, icase, cases, form

    def _build_normals_quality(self, model, nelements, cases, form0, icase,
                               xyz_cid0, material_coord,
                               min_interior_angle, max_interior_angle, dideal_theta,
                               area, max_skew_angle, taper_ratio,
                               max_warp_angle, area_ratio, max_aspect_ratio):
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
                elif pid_type in ['PCOMP', 'PCOMPG']:
                    z0 = element.pid.z0
                elif pid_type == 'PLPLANE':
                    z0 = 0.
                elif pid_type == 'PSHEAR':
                    z0 = 0.
                elif pid_type in ['PSOLID', 'PLSOLID']:
                    z0 = 0.
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
                    elif element.type == 'CTRAX3':
                        nnodesi = 3
                        z0 = 0.
                    elif element.type == 'CTRAX6':
                        nnodesi = 6
                        z0 = 0.
                    elif element.type == 'CQUADX4':
                        nnodesi = 4
                        z0 = 0.
                    elif element.type == 'CQUADX8':
                        nnodesi = 8
                        z0 = 0.
                    elif element.type == 'CQUADX':
                        nnodesi = 9
                        z0 = 0.
                    else:
                        raise NotImplementedError(element.type)
                else:
                    if element.type in ['CTRIA3', 'CTRIAR', 'CTRAX3', 'CTRIAX', 'CPLSTN3']:
                        nnodesi = 3
                    elif element.type in ['CTRIA6', 'CTRIAX6', 'CPLSTN6', 'CTRAX6']:
                        nnodesi = 6

                    elif element.type in ['CQUAD4', 'CQUADR', 'CPLSTN4', 'CSHEAR', 'CQUADX4']:
                        nnodesi = 4
                    elif element.type in ['CQUAD8', 'CPLSTN8', 'CQUADX8']:
                        nnodesi = 8
                    elif element.type == ['CQUAD', 'CQUADX']:
                        nnodesi = 9
                    else:
                        raise NotImplementedError(element.type)

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
            elif element.type in ['CBUSH', 'CBUSH1D', 'CBUSH2D',
                                  'CFAST', 'CVISC',
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
        is_element_dim = element_dim.max() != element_dim.min()
        is_element_dim = True
        if is_element_dim:
            eid_dim_res = GuiResult(0, header='ElementDim', title='ElementDim',
                                    location='centroid', scalar=element_dim)
            cases[icase] = (eid_dim_res, (0, 'ElementDim'))

        is_shell = np.abs(normals).max() > 0.
        is_solid = np.abs(max_interior_angle).max() > 0.
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

            # this is just for testing nan colors that doesn't work
            #max_interior_angle[:1000] = np.nan
            area_res = GuiResult(0, header='Area', title='Area',
                                 location='centroid', scalar=area)
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
                arearatio_res = GuiResult(
                    0, header='Area Ratio', title='Area Ratio',
                    location='centroid', scalar=area_ratio)
                cases[icase] = (arearatio_res, (0, 'Area Ratio'))
                form_checks.append(('Area Ratio', icase, []))
                icase += 1

            if taper_ratio.max() > 1.:
                taperratio_res = GuiResult(
                    0, header='Taper Ratio', title='Taper Ratio',
                    location='centroid', scalar=taper_ratio)
                cases[icase] = (taperratio_res, (0, 'Taper Ratio'))
                form_checks.append(('Taper Ratio', icase, []))
                icase += 1

            if max_warp_angle.max() > 0.0:
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
            material_coord_res = GuiResult(
                0, header='MaterialCoord', title='MaterialCoord',
                location='centroid',
                scalar=material_coord, data_format='%i')
            cases[icase] = (material_coord_res, (0, 'MaterialCoord'))
            form0.append(('MaterialCoord', icase, []))
            icase += 1
        return icase, normals

    def _build_properties(self, model, nelements, eids, pids, cases, form0, icase):
        """
        creates:
          - PropertyID
        """
        prop_types_with_mid = [
            'PSOLID', 'PSHEAR',
            'PROD', 'CROD', 'PTUBE', 'PBAR', 'PBARL', 'PBEAM', 'PBEAML',
        ]

        upids = None
        pid_res = GuiResult(0, header='PropertyID', title='PropertyID',
                            location='centroid', scalar=pids)
        cases[icase] = (pid_res, (0, 'PropertyID'))
        form0.append(('PropertyID', icase, []))
        icase += 1

        upids = np.unique(pids)
        mid_eids_skip = []

        nplies = 1
        if 'PSHELL' in model.card_count:
            nplies = 4
        for pid in model.get_card_ids_by_card_types(['PCOMP', 'PCOMPG'], combine=True):
            prop = model.properties[pid]
            nplies = max(nplies, prop.nplies)
        if nplies > 1:
            nplies += 1

        mids = np.zeros((nelements, nplies), dtype='int32')
        thickness = np.zeros((nelements, nplies), dtype='float32')
        rho = np.zeros((nelements, nplies), dtype='float32')
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
                mids[i] = mid
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
                for iply in range(prop.nplies):
                    mids[i, iply+1] = prop.Mid(iply)
                    thickness[i, iply+1] = prop.Thickness(iply)
                thickness[i, 0] = thickness[i, :].sum()
                mids[i, 0] = mids[i, 1]
            elif prop.type in ['PELAS', 'PBUSH', 'PDAMP', 'PDAMPT']:
                i = np.where(pids == pid)[0]
                mid_eids_skip.append(i)
            else:
                print('material for pid=%s type=%s not considered' % (pid, prop.type))

        #print('mids =', mids)
        mids = mids[:, 0]
        thickness = thickness[:, 0]
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
        return icase, upids, mids, thickness

    def _build_materials(self, model, mids, thickness, cases, form0, icase):
        """
        creates:
          - Material ID
          - E_11
          - E_22
          - E_33
          - Is Isotropic?
        """
        has_mat8, has_mat9, e11, e22, e33 = self._get_material_arrays(model, mids)

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
            cases[icase + 1] = (e11_res, (0, 'E_11'))
            cases[icase + 2] = (e22_res, (0, 'E_22'))
            form0.append(('E_11', icase + 1, []))
            form0.append(('E_22', icase + 2, []))
            icase += 3

            is_isotropic = np.zeros(len(e11), dtype='int8')
            if has_mat9:
                is_isotropic[(e11 == e22) | (e11 == e33)] = 1
                e33_res = GuiResult(0, header='E_33', title='E_33',
                                    location='centroid', scalar=e33, data_format='%.3e')
                cases[icase] = (e33_res, (0, 'E_33'))
                form0.append(('E_33', icase, []))
                icase += 1
            else:
                #is_isotropic_map = e11 == e22
                is_isotropic[e11 == e22] = 1

            iso_res = GuiResult(0, header='IsIsotropic?', title='IsIsotropic?',
                                location='centroid', scalar=is_isotropic, data_format='%i')
            cases[icase] = (iso_res, (0, 'Is Isotropic?'))
            form0.append(('Is Isotropic?', icase, []))
            icase += 1
        else:
            # isotropic
            e11_res = GuiResult(0, header='E', title='E',
                                location='centroid', scalar=e11, data_format='%.3e')
            cases[icase + 1] = (e11_res, (0, 'E'))
            form0.append(('E', icase + 1, []))
            icase += 2
        return icase

    def _build_optimization(self, model, pids, upids, nelements, cases, form0, icase):
        if upids is None:
            return icase
        if len(model.properties) and len(model.dvprels):
            # len(model.dvprels) + len(model.dvcrels) + len(model.dvmrels) + len(model.desvars)
            #dvmrel_init = np.zeros(nelements, dtype='int32')
            #dvgrel_init = np.zeros(nelements, dtype='int32')
            out = self._get_dvprel_ndarrays(model, nelements, pids)
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

    def _plot_pressures(self, model, cases, form0, icase, subcase_id, subcase):
        """
        pressure act normal to the face (as opposed to anti-normal)
        """

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

        pressures = self.get_pressure_array(model, load_case)

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

    def _plot_applied_loads(self, model, cases, form0, icase, subcase_id, subcase):
        form = []
        found_load, found_temperature, temperature_data, load_data = self.get_load_array(model, subcase)
        if found_load:
            centroidal_pressures, forces, spcd = load_data
            if np.abs(centroidal_pressures).max():
                # print('iload=%s' % iload)
                # print(case_name)
                pressure_res = GuiResult(subcase_id, header='Pressure', title='Pressure',
                                         location='centroid', scalar=centroidal_pressures)
                cases[icase] = (pressure_res, (0, 'Pressure'))
                form0.append(('Pressure', icase, []))
                icase += 1

            if np.abs(forces.max() - forces.min()) > 0.0:
                load_x_res = GuiResult(subcase_id, header='LoadX', title='LoadX',
                                       location='node', scalar=forces[:, 0])
                load_y_res = GuiResult(subcase_id, header='LoadY', title='LoadY',
                                       location='node', scalar=forces[:, 1])
                load_z_res = GuiResult(subcase_id, header='LoadZ', title='LoadZ',
                                       location='node', scalar=forces[:, 2])
                cases[icase] = (load_x_res, (0, 'LoadX'))
                cases[icase + 1] = (load_y_res, (0, 'LoadY'))
                cases[icase + 2] = (load_z_res, (0, 'LoadZ'))

                # if forces[:, 0].min() != forces[:, 0].max():
                # if forces[:, 1].min() != forces[:, 1].max():
                # if forces[:, 2].min() != forces[:, 2].max():

                form0.append(('Total Load FX', icase, []))
                form0.append(('Total Load FY', icase + 1, []))
                form0.append(('Total Load FZ', icase + 2, []))
                icase += 3

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
        if found_temperature:
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

    def clear_nastran(self):
        """cleans up variables specific to Nastran"""
        self.eid_map = {}
        self.nid_map = {}
        self.eid_to_nid_map = {}
        self.element_ids = None
        self.node_ids = None
