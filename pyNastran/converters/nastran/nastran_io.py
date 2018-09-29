# pylint: disable=E1101,C1801,C0103
"""
Defines the GUI IO file for Nastran.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import os
import sys
from collections import defaultdict, OrderedDict
import traceback
from six import iteritems, StringIO, string_types
from pyNastran.op2.result_objects.stress_object import StressObject

SIDE_MAP = {}
SIDE_MAP['CHEXA'] = {
    1 : [4, 3, 2, 1],
    2 : [1, 2, 6, 5],
    3 : [2, 3, 7, 6],
    4 : [3, 4, 8, 7],
    5 : [4, 1, 5, 8],
    6 : [5, 6, 7, 8],
}

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
from numpy.linalg import norm  # type: ignore

#: makes vtk work on certain builds of vtk
#: we have to call this before vtk; you can't just try-except it
#: unused_import
from pyNastran.gui.qt_version import qt_version
if qt_version == 'pyqt4':
    import PyQt4
elif qt_version == 'pyqt5':
    import PyQt5
elif qt_version == 'pyside':
    import PySide
elif qt_version == 'pyside2':
    import PySide2
else:
    raise NotImplementedError(qt_version)

import vtk
from vtk import (vtkTriangle, vtkQuad, vtkTetra, vtkWedge, vtkHexahedron,
                 vtkQuadraticTriangle, vtkQuadraticQuad, vtkQuadraticTetra,
                 vtkQuadraticWedge, vtkQuadraticHexahedron,
                 vtkPyramid) #vtkQuadraticPyramid

#from pyNastran import is_release
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.utils.log import properties as log_properties
from pyNastran.utils.numpy_utils import isfinite_and_nonzero, isfinite_and_greater_than, isfinite
from pyNastran.bdf.bdf import (BDF,
                               CAERO1, CAERO2, CAERO3, CAERO4, CAERO5,
                               CQUAD4, CQUAD8, CQUAD, CQUADR, CSHEAR,
                               CTRIA3, CTRIA6, CTRIAR,
                               CPLSTN3, CPLSTN4, CPLSTN6, CPLSTN8,
                               CTRAX3, CTRIAX6, CTRIAX, #CTRAX6,
                               CQUADX4, CQUADX8, CQUADX,
                               CONM2)
from pyNastran.bdf.cards.aero.zona import CAERO7, BODY7
from pyNastran.bdf.cards.elements.shell import ShellElement
from pyNastran.bdf.cards.elements.solid import (
    CTETRA4, CTETRA10, CPENTA6, CPENTA15,
    CHEXA8, CHEXA20, CIHEX1, CIHEX2,
    CPYRAM5, CPYRAM13,
)

from pyNastran.gui.utils.vtk.vtk_utils import (
    get_numpy_idtype_for_vtk, numpy_to_vtk_points,
    create_vtk_cells_of_constant_element_type,
    numpy_to_vtk, numpy_to_vtkIdTypeArray)
from pyNastran.gui.errors import NoGeometry
from pyNastran.gui.gui_objects.gui_result import GuiResult, NormalResult
from pyNastran.converters.nastran.geometry_helper import (
    NastranGeometryHelper, tri_quality, quad_quality, get_min_max_theta,
    get_material_arrays, get_suport_node_ids)
from pyNastran.converters.nastran.results_helper import NastranGuiResults, _get_times
from pyNastran.converters.nastran.displacements import (
    ForceTableResults, ElementalTableResults)

from pyNastran.op2.op2 import OP2
try:
    import h5py
    IS_H5PY = True
except ImportError:
    IS_H5PY = False
#from pyNastran.f06.f06_formatting import get_key0
from pyNastran.op2.op2_geom import OP2Geom

from pyNastran.gui.qt_files.colors import (
    RED_FLOAT, BLUE_FLOAT, GREEN_FLOAT, LIGHT_GREEN_FLOAT, PINK_FLOAT, PURPLE_FLOAT,
    YELLOW_FLOAT, ORANGE_FLOAT)
#GREEN = (0., 1., 0.)
#BLUE = (0., 0., 1.)
#LIGHT_GREEN = (0.5, 1., 0.5)
#PINK = (0.98, 0.4, 0.93)
#ORANGE = (219/255., 168/255., 13/255.)
#RED = (1., 0., 0.)
#YELLOW = (1., 1., 0.)
#PURPLE = (1., 0., 1.)

NO_THETA = [
    'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
    'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
    'CBAR', 'CBEAM', 'CBEAM3',
    'CBUSH', 'CBUSH1D', 'CBUSH2D', 'CVISC',
    'CONROD', 'CROD', 'PLOTEL',
]

IS_TESTING = 'test' in sys.argv[0]
class NastranIO(NastranGuiResults, NastranGeometryHelper):
    """
    Defines the GUI class for Nastran.
    """
    def __init__(self):
        super(NastranIO, self).__init__()
        self.nid_release_map = {}
        self.stress = {}
        self.strain = {}

    @property
    def gui(self):
        return self

    #def __init__(self, gui):
        #super(NastranIO, self).__init__()
        #self.gui = gui  # make sure to comment out the property on line 124
        #self.nid_release_map = {}
        #self.stress = {}
        #self.strain = {}


    def get_nastran_wildcard_geometry_results_functions(self):
        """
        gets the Nastran wildcard loader used in the file load menu
        """
        geom_methods_bdf = 'Nastran Geometry - BDF (*.bdf; *.dat; *.nas; *.ecd; *.op2; *.pch)'
        geom_methods_pch = 'Nastran Geometry - Punch (*.bdf; *.dat; *.nas; *.ecd; *.pch)'
        combined_methods_op2 = 'Nastran Geometry + Results - OP2 (*.op2)'

        results_fmts = ['Nastran OP2 (*.op2)',]
        if IS_H5PY:
            results_fmts.append('pyNastran H5 (*.h5)')
        results_fmts.append('Patran nod (*.nod)')
        results_fmt = ';;'.join(results_fmts)
        #results_fmt = 'Nastran OP2 (*.op2)'

        data_geom = (
            'nastran',
            geom_methods_bdf, self.load_nastran_geometry,
            results_fmt, self.load_nastran_results)

        data_geom_pch = (
            'nastran',
            geom_methods_pch, self.load_nastran_geometry,
            results_fmt, self.load_nastran_results)

        unused_data_geom_results = (
            'nastran',
            combined_methods_op2, self.load_nastran_geometry_and_results,
            results_fmt, self.load_nastran_results)

        return [data_geom, data_geom_pch]
        #return [data_geom, data_geom_pch, data_geom_results]

    def load_nastran_geometry_and_results(self, op2_filename, name='main', plot=True):
        """loads geometry and results, so you don't have to double define the same BDF/OP2"""
        self.load_nastran_geometry(op2_filename, name='main', plot=False)
        self.load_nastran_results(self.model) # name='main', plot=True

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
            #('about_nastran', 'About Nastran GUI', 'tabout.png', 'CTRL+H',
            #'About Nastran GUI and help on shortcuts', self.about_dialog),
            #('about', 'About Orig GUI', 'tabout.png', 'CTRL+H',
            #'About Nastran GUI and help on shortcuts', self.about_dialog),
        ]
        #self.gui.menu_help2 = self.gui.menubar.addMenu('&HelpMenuNew')
        #self.gui.menu_help.menuAction().setVisible(False)
        if hasattr(self, 'nastran_toolbar'):
            self.gui.nastran_toolbar.setVisible(True)
            self.gui.actions['nastran'].setVisible(True)
        else:
            self.gui.nastran_toolbar = self.addToolBar('Nastran Toolbar')
            self.gui.nastran_toolbar.setObjectName('nastran_toolbar')
            #self.gui.nastran_toolbar.setStatusTip("Show/Hide nastran toolbar")
            self.gui.actions['nastran'] = self.nastran_toolbar.toggleViewAction()
            self.gui.actions['nastran'].setStatusTip("Show/Hide application toolbar")
        #self.gui.file.menuAction().setVisible(False)
        #self.gui.menu_help.

        #self.gui.actions['about'].Disable()
        menu_items = {}
        menu_items['nastran_toolbar'] = (self.gui.nastran_toolbar, ('caero', 'caero_subpanels', 'conm2'))
        #menu_items = [
            #(self.menu_help2, ('about_nastran',)),
            #(self.gui.nastran_toolbar, ('caero', 'caero_subpanels', 'conm2'))
            #(self.menu_window, tuple(menu_window)),
            #(self.menu_help, ('load_geometry', 'load_results', 'script', '', 'exit')),
            #(self.menu_help2, ('load_geometry', 'load_results', 'script', '', 'exit')),

        return tools, menu_items

    def toggle_caero_panels(self):
        """
        Toggle the visibility of the CAERO panels. The visibility of the
        sub panels or panels will be set according to the current
        show_caero_sub_panels state.
        """
        if not self.has_caero:
            return
        self.show_caero_actor = not self.show_caero_actor

        names = ['caero', 'caero_subpanels', 'caero_control_surfaces']
        geometry_properties = self.gui._get_geometry_properties_by_name(names)

        if self.show_caero_actor:
            try:
                geometry_properties['caero_control_surfaces'].is_visible = True
            except KeyError:
                pass
            if self.show_caero_sub_panels:
                geometry_properties['caero_subpanels'].is_visible = True
            else:
                geometry_properties['caero'].is_visible = True
        else:
            try:
                geometry_properties['caero_control_surfaces'].is_visible = False
            except KeyError:
                pass
            geometry_properties['caero'].is_visible = False
            geometry_properties['caero_subpanels'].is_visible = False
        self.gui.on_update_geometry_properties_override_dialog(geometry_properties)

    def _get_geometry_properties_by_name(self, names):
        """
        Get a subset of the self.geometry_properties dict specified by
        names.  Any names not in the dict will be ignored.

        Parameters
        -----------
        names : list [str, ...]
            List of names.

        Returns
        --------
        geometry_properties : dict {str : AltGeometry or CoordProperties}
            Dictonairy from name to property object.

        """
        geometry_properties = {}
        for name in names:
            try:
                prop = self.gui.geometry_properties[name]
            except KeyError:
                continue
            geometry_properties[name] = prop
        return geometry_properties

    def on_update_geometry_properties_window(self, geometry_properties):
        """updates the 'Edit Geometry Properties' window"""
        self.gui.on_update_geometry_properties_window(geometry_properties)

    def toggle_caero_sub_panels(self):
        """
        Toggle the visibility of the CAERO sub panels
        """
        if not self.has_caero:
            return

        names = ['caero', 'caero_subpanels']
        geometry_properties = self.gui._get_geometry_properties_by_name(names)

        self.show_caero_sub_panels = not self.show_caero_sub_panels
        if self.show_caero_actor:
            if self.show_caero_sub_panels:
                geometry_properties['caero'].is_visible = False
                geometry_properties['caero_subpanels'].is_visible = True
            else:
                geometry_properties['caero'].is_visible = True
                geometry_properties['caero_subpanels'].is_visible = False
        self.gui.on_update_geometry_properties_override_dialog(geometry_properties)

    def toggle_conms(self):
        """
        Toggle the visibility of the CONMS
        """
        name = 'conm2'
        if name in self.gui.geometry_actors:
            geometry_properties_change = {name : self.gui.geometry_properties[name]}
            visibility_prev = geometry_properties_change[name].is_visible
            geometry_properties_change[name].is_visible = not visibility_prev

            self.gui.on_update_geometry_properties_override_dialog(geometry_properties_change)

    def _create_coord(self, dim_max, cid, coord, coord_type):
        """
        Create a coordinate system

        Parameters
        ----------
        dim_max : float
            the max model dimension; 10% of the max will be used for the
            coord length
        cid : int
           the coordinate system id
        coord : Coord()
           the Nastran coord object
        coord_type : str
            a string of 'xyz', 'Rtz', 'Rtp' (xyz, cylindrical, spherical)
            that changes the axis names
        """
        origin = coord.origin
        beta = coord.beta().T
        ## TODO: support FEMAP syntax
        self.gui.create_coordinate_system(
            cid, dim_max, label='%s' % cid, origin=origin,
            matrix_3x3=beta, coord_type=coord_type)

    def _create_nastran_coords(self, model, dim_max):
        """
        Creates the Nastran coordinate systems.

        Parameters
        ----------
        model : BDF()
            the BDF object
        dim_max : float
            the max model dimension; 10% of the max will be used for the
            coord length
        """
        cid_types = {
            'R' : 'xyz',
            'C' : 'Rtz',
            'S' : 'Rtp',
        }
        self.gui.create_global_axes(dim_max)
        for cid, coord in sorted(model.coords.items()):
            if cid in [0, -1]:
                continue
            cid_type = cid_types[coord.Type]
            self.gui._create_coord(dim_max, cid, coord, cid_type)

    def _remove_old_nastran_geometry(self, bdf_filename):
        """cleans up the nastran model"""
        #return self._remove_old_geometry(bdf_filename)

        # skip_reading = self.removeOldGeometry(bdf_filename)
        skip_reading = False
        if bdf_filename is None or bdf_filename is '':
            #self.grid = vtk.vtkUnstructuredGrid()
            #self.scalarBar.VisibilityOff()
            skip_reading = True
            return skip_reading
        else:
            self.gui.turn_text_off()
            self.gui.grid.Reset()

            #self.gui.eid_map = {}
            #self.gui.nid_map = {}

            self.gui.result_cases = {}
            self.gui.ncases = 0

        # TODO: is this doing anything?
        for name in ('case_keys', 'icase', 'isubcase_name_map'):
            if hasattr(self, name):
                del name
        return skip_reading

    def get_xyz_in_coord(self, model, cid=0, fdtype='float32'):
        """
        Creates the grid points efficiently

        Used by ``load_nastran_geometry_unvectorized``
        """
        #import time
        #time0 = time.time()

        # t=.578
        #print("get_displacement_index_xyz_cp_cd")
        out = model.get_displacement_index_xyz_cp_cd(
            fdtype=fdtype, idtype='int32', sort_ids=True)
        icd_transform, icp_transform, xyz_cp, nid_cp_cd = out
        self.i_transform = icd_transform

        #print("transform_xyzcp_to_xyz_cid")
        xyz_cid0 = model.transform_xyzcp_to_xyz_cid(
            xyz_cp, nid_cp_cd[:, 0], icp_transform, cid=cid,
            in_place=False)

        #print('model.spoints =', model.spoints)
        #import json
        #for spoint_id, spoint in model.spoints.items():
            #if spoint.comment: # or spoint._comment?
                #print('SPOINT comment=%r _comment=%r' % (spoint.comment, spoint._comment))
                #comment_lower = spoint.comment.lower()
                #print('comment_lower = %r' % comment_lower)
                ## pyNastran: SPOINT={'id':10, 'xyz':[10.,10.,10.]}
                #if 'pynastran' in comment_lower and 'spoint' in comment_lower:
                    #dict_str = comment_lower.split('=')[1].rstrip().replace("'", '"').replace('}', ',}').replace(',,}', ',}')
                    #print('dict_str = %r' % dict_str)
                    #dicti = json.loads(dict_str)
                    #print(dicti)
        #for epoint_id, epoint in model.epoints.items():
            #if epoints.comment:
                #print('EPOINT comment=%r _comment=%r' % (spoint.comment, spoint._comment))
        #sys.stdout.flush()

        nid_map = self.gui.nid_map
        for i, nid in enumerate(nid_cp_cd[:, 0]):
            nid_map[nid] = i

        self._add_nastran_spoints_to_grid(model)
        #print('dt_nastran_xyz =', time.time() - time0)
        return xyz_cid0, nid_cp_cd

    def get_xyz_in_coord_vectorized(self, model, cid=0, fdtype='float32'):
        """
        Creates the grid points efficiently

        Used by ``load_nastran_geometry_vectorized``
        """
        #import time
        #time0 = time.time()

        # t=.578
        #print("get_displacement_index_xyz_cp_cd")
        xyz_cid0 = None
        nid_cp_cd = None
        if self.gui.nnodes > 0:
            out = model.get_displacement_index_xyz_cp_cd(
                fdtype=fdtype, idtype='int32')
            icd_transform, icp_transform, xyz_cp, nid_cp_cd = out
            self.i_transform = icd_transform

            #print("transform_xyzcp_to_xyz_cid")
            #model.nodes.cp = nid_cp_cd[:, 1]
            xyz_cid0 = model.transform_xyzcp_to_xyz_cid(
                xyz_cp, nid_cp_cd[:, 0], icp_transform, cid=cid,
                in_place=False)
            model.nodes.xyz_cid0 = xyz_cid0
            model.nodes.nids = nid_cp_cd[:, 0]

            nid_map = self.gui.nid_map
            for i, nid in enumerate(nid_cp_cd[:, 0]):
                nid_map[nid] = i

            self._add_nastran_spoints_to_grid(model)
        #print('dt_nastran_xyz =', time.time() - time0)
        return xyz_cid0, nid_cp_cd

    def _get_model_unvectorized(self, bdf_filename, xref_loads=True):
        """Loads the BDF/OP2 geometry"""
        ext = '.bdf'
        if isinstance(bdf_filename, str):
            ext = os.path.splitext(bdf_filename)[1].lower()

        punch = False
        if ext == '.pch':
            punch = True

        self.model_type = 'nastran'
        if ext == '.op2':
            model = OP2Geom(make_geom=True, debug=False, log=self.gui.log,
                            debug_file=None)
            model.clear_results()
            model.read_op2(op2_filename=bdf_filename)
        else:  # read the bdf/punch
            model = BDF(log=self.gui.log, debug=True)
            model.read_bdf(bdf_filename,
                           punch=punch, xref=False,
                           validate=True)
            #print('done with read_bdf')
            #xref_loads = False
        #xref_aero = len(model.caeros) > 0

        xref_nodes = True
        #model.cross_reference()
        model.safe_cross_reference(
            xref=True,
            xref_nodes=xref_nodes,
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
        return model, xref_nodes

    def load_nastran_geometry(self, bdf_filename, name='main', plot=True, **kwargs):
        """
        The entry point for Nastran geometry loading.

        Parameters
        ----------
        bdf_filename : str
            the Nastran filename to load
        name : str
            the name of the "main" actor for the GUI
        plot : bool; default=True
            should the model be generated or should we wait until
            after the results are loaded

        kwargs:
        -------
        is_geometry_results : bool; default=True
            code is being called from load_nastran_geometry_and_results
            not used...
        """
        self.gui.eid_maps[name] = {}
        self.gui.nid_maps[name] = {}
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

        #load_geom = True
        if isinstance(bdf_filename, str) and bdf_filename.lower().endswith(('.bdf', '.dat', '.pch',)): # '.op2'
            if IS_TESTING or self.gui.is_testing_flag:
                self.load_nastran_geometry_vectorized(bdf_filename, plot=plot)
                self.load_nastran_geometry_unvectorized(bdf_filename, plot=plot)
            else:
                self.load_nastran_geometry_unvectorized(bdf_filename, plot=plot)
                #self.load_nastran_geometry_vectorized(bdf_filename, plot=plot)
        else:
            self.load_nastran_geometry_unvectorized(bdf_filename, plot=plot)

    def load_nastran_geometry_vectorized(self, bdf_filename, plot=True):
        """
        The entry point for Nastran geometry loading.

        Parameters
        ----------
        bdf_filename : str
            the Nastran filename to load
        plot : bool; default=True
            should the model be generated or should we wait until
            after the results are loaded
        """
        #self.isubcase_name_map[None] = ['a', 'b']
        reset_labels = True
        if plot:
            self.gui.scalarBar.VisibilityOff()
            self.gui.scalarBar.Modified()

        model = self._get_model_vectorized(bdf_filename)

        nnodes = len(model.grid)
        nspoints = len(model.spoints)
        nepoints = len(model.epoints)
        ncaero_cards = len(model.caeros)
        ngridb = len(model.gridb)
        #if model.spoints:
            #spoints = sorted([spoint.nid for spoint in itervalues(model.spoints)])
        #if model.epoints:
            #epoints = sorted([epoint.nid for epoint in itervalues(model.epoints)])

        ngui_nodes = nnodes + nspoints + nepoints + ngridb
        if ngui_nodes + ncaero_cards == 0:
            msg = 'nnodes + nspoints + nepoints = 0\n'
            msg += 'card_count = %r' % str(model.card_count)
            raise NoGeometry(msg)

        nelements2 = len(model.elements2)
        #nelements = len(model.elements) + nelements2
        nelements = nelements2

        nmasses = len(model.masses)
        nplotels = len(model.plotels)
        nrigid = len(model.rigid_elements)

        #if 'CTRIAX6' in model.card_count and nelements == 0: # skipping...
            #return

        #nmpc = len(model.mpcs)  # really should only be allowed if we have it in a subcase
        if nelements + nmasses + ncaero_cards + nplotels + nrigid == 0:
            msg = 'nelements + nmasses + ncaero_cards + nplotels + nrigid = 0\n'
            msg += 'card_count = %r' % str(model.card_count)
            raise NoGeometry(msg)

        self.gui.nnodes = ngui_nodes
        self.gui.nelements = nelements  # approximate...

        self.gui.log_info("nnodes=%i nelements=%i" % (self.nnodes, self.nelements))
        msg = model.get_bdf_stats(return_type='string')
        self.gui.log_debug(msg)
        msg = model.get_bdf_stats(return_type='list')

        # this call will break the GUI if there are a lot of lines and
        # by a lot I mean 37641.  It's fine for a single call.
        #for msgi in msg:
            #model.log.debug(msgi)

        nconm2 = 0
        #if 'CONM2' in model.card_count:
            #nconm2 += model.card_count['CONM2']
        #if 'CMASS1' in model.card_count:
            #nconm2 += model.card_count['CMASS1']
        #if 'CMASS2' in model.card_count:
            #nconm2 += model.card_count['CMASS2']

        if nconm2 > 0:
            self.gui.create_alternate_vtk_grid(
                'conm2', color=ORANGE_FLOAT, line_width=5, opacity=1., point_size=4,
                representation='point', follower_function=None)

        # Allocate grids
        self.gui.grid.Allocate(self.nelements, 1000)
        #self._create_caero_actors(ncaeros, ncaeros_sub, ncaeros_cs, has_control_surface)
        #if nconm2 > 0:
            #self.gui.alt_grids['conm2'].Allocate(nconm2, 1000)

        if self.save_data:
            self.model = model

        #-----------------------------------------------------------------------
        # nodes/coords
        #print('get_xyz_in_coord')
        dim_max = 1.0
        xyz_cid0, nid_cp_cd = self.get_xyz_in_coord_vectorized(model, cid=0, fdtype='float32')
        if xyz_cid0 is not None:
            dim_max = self._points_to_vtkpoints_coords(model, xyz_cid0)
        #-----------------------------------------------------------------------

        #------------------------------------------------------------
        # TEMP
        j = 0
        results = self._map_elements_vectorized(self.nid_map, model, j, dim_max,
                                                nid_cp_cd, plot=True, xref_loads=True)
        has_control_surface = False
        geometry_names = []
        #------------------------------------------------------------
        cases = OrderedDict()
        form = ['Geometry', None, []]
        form0 = form[2]

        subcase_id = 0

        colormap = self.gui.settings.colormap
        if self.gui.nnodes > 0:
            icase = 0
            all_nids = nid_cp_cd[:, 0]
            self.gui.node_ids = all_nids

            nid_res = GuiResult(subcase_id, 'NodeID', 'NodeID', 'node', all_nids,
                                mask_value=0,
                                nlabels=None,
                                labelsize=None,
                                ncolors=None,
                                colormap=colormap,
                                data_format=None,
                                uname='GuiResult')
            cases[icase] = (nid_res, (0, 'Node ID'))
            form0.append(('Node ID', icase, []))
            icase += 1

            nid_res = GuiResult(subcase_id, 'iNode', 'iNode', 'node',
                                np.arange(len(all_nids), dtype='int32'),
                                mask_value=0,
                                nlabels=None,
                                labelsize=None,
                                ncolors=None,
                                colormap=colormap,
                                data_format=None,
                                uname='GuiResult')
            cases[icase] = (nid_res, (0, 'Node ID'))
            form0.append(('iNode', icase, []))
            icase += 1

            # this intentionally makes a deepcopy
            cds = np.array(nid_cp_cd[:, 2])
            if cds.max() > 0:
                cd_res = GuiResult(0, header='NodeCd', title='NodeCd',
                                   location='node', scalar=cds, colormap=colormap)
                cases[icase] = (cd_res, (0, 'NodeCd'))
                form0.append(('NodeCd', icase, []))
                icase += 1

        if self.gui.nelements > 0:
            eids_array = results['eid']
            eid_res = GuiResult(subcase_id, 'ElementID', 'ElementID', 'centroid', eids_array,
                                mask_value=0,
                                nlabels=None,
                                labelsize=None,
                                ncolors=None,
                                colormap=colormap,
                                data_format=None,
                                uname='GuiResult')
            cases[icase] = (eid_res, (0, 'ElementID'))
            form0.append(('ElementID', icase, []))
            icase += 1

            eids_array = results['eid']
            eid_res = GuiResult(subcase_id, 'iElement', 'iElement', 'centroid',
                                np.arange(len(eids_array), dtype='int32'),
                                mask_value=-1,
                                nlabels=None,
                                labelsize=None,
                                ncolors=None,
                                colormap=colormap,
                                data_format=None,
                                uname='GuiResult')
            cases[icase] = (eid_res, (0, 'iElement'))
            form0.append(('iElement', icase, []))
            icase += 1

            #is_element_dim = True
            dim_array = results['dim']
            if len(np.unique(dim_array)) > 1:
                dim_res = GuiResult(subcase_id, 'ElementDim', 'ElementDim', 'centroid', dim_array,
                                    mask_value=-1,
                                    nlabels=None,
                                    labelsize=None,
                                    ncolors=None,
                                    colormap=colormap,
                                    data_format=None,
                                    uname='GuiResult')
                cases[icase] = (dim_res, (0, 'ElementDim'))
                form0.append(('ElementDim', icase, []))
                icase += 1

            nnodes_array = results['nnodes']
            if nnodes_array.max() > -1:
                nnodes_res = GuiResult(subcase_id, 'NNodes/Elem', 'NNodes/Elem',
                                       'centroid', nnodes_array,
                                       mask_value=0,
                                       nlabels=None,
                                       labelsize=None,
                                       ncolors=None,
                                       colormap=colormap,
                                       data_format=None,
                                       uname='GuiResult')
                cases[icase] = (nnodes_res, (0, 'NNodes/Elem'))
                form0.append(('NNodes/Elem', icase, []))
                icase += 1

            pids_array = results['pid']
            pid_res = GuiResult(0, header='PropertyID', title='PropertyID',
                                location='centroid', scalar=pids_array, mask_value=0)
            cases[icase] = (pid_res, (0, 'PropertyID'))
            form0.append(('PropertyID', icase, []))
            icase += 1

            #upids = np.unique(pids_array)
            unused_mid_eids_skip = []


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
            is_pshell_pcomp = (is_pshell, is_pcomp)
            nplies = max(nplies, pcomp_nplies + 1)
            mids = np.zeros((nelements, nplies), dtype='int32')
            thickness = np.full((nelements, nplies), np.nan, dtype='float32')
            #rho = np.full((nelements, nplies), np.nan, dtype='float32')
            nplies = np.zeros(nelements, dtype='int32')

            # materials
            upids = np.unique(pids_array)
            ipids = np.zeros(len(pids_array), dtype='int32')
            iupid = 0
            for upid in upids: # upid_old
                if upid == 0:
                    # elements w/o properties
                    continue

                ipid = np.where(pids_array == upid)[0]
                ipids[ipid] = iupid
                if len(ipid):
                    try:
                        prop = model.properties[upid]
                    except KeyError:
                        raise KeyError('pid=%r properties=%s' % (upid, str(model.properties)))
                    if prop.type == 'PSHELL':
                        nplies[ipid] = 4
                        thickness[ipid, 0] = prop.Thickness()
                    elif prop.type in ['PCOMP', 'PCOMPG']:
                        nplies[ipid] = prop.nplies
                        for iply in range(prop.nplies):
                            mids[ipid, iply+1] = prop.Mid(iply)
                            thickness[ipid, iply+1] = prop.Thickness(iply)
                    #else:
                        #self.log.error('nplies\n%s' % str(prop))
                    iupid += 1

            if len(model.conrod):
                #mids[ieid, 0] = 42
                pass

            pid_res = GuiResult(0, header='iProperty', title='iProperty',
                                location='centroid', scalar=ipids, colormap=colormap)
            cases[icase] = (pid_res, (0, 'iProperty'))
            form0.append(('iProperty', icase, []))
            icase += 1

            #if nplies.max() > 0:
                #nplies_res = GuiResult(0, header='Number of Plies', title='nPlies',
                                       #location='centroid', scalar=nplies, mask_value=0)
                #cases[icase] = (nplies_res, (0, 'Number of Plies'))
                #form0.append(('Number of Plies', icase, []))
                #icase += 1

            pshell = {
                'mids' : mids,
                'thickness' : nplies,
            }
            pcomp = {
                'mids' : mids,
                'thickness' : nplies,
                'nplies' : nplies,
            }
            icase = self._build_materials(model, pshell, pcomp, is_pshell_pcomp,
                                          cases, form0, icase)

        #------------------------------------------------------------
        # add alternate actors
        self.gui._add_alt_actors(self.gui.alt_grids)

        # set default representation
        self._set_caero_representation(has_control_surface)

        for grid_name in geometry_names:
            if grid_name in self.gui.geometry_actors:
                self.gui.geometry_actors[grid_name].Modified()

        #self.gui.grid_mapper.SetResolveCoincidentTopologyToPolygonOffset()
        if 0:
            if plot:
                #self.log.info(cases.keys())
                self.gui._finish_results_io2([form], cases, reset_labels=reset_labels)
            else:
                self.gui._set_results([form], cases)


    def _map_elements_vectorized(self, unused_nid_map, model, unused_j, unused_dim_max,
                                 unused_nid_cp_cd, plot=True, xref_loads=True):
        """
        Much, much faster way to add elements that directly builds the
        VTK objects rather than using for loops.

        Parameters
        ----------
        nid_map : ???
           ???
        model : BDF()
            the BDF model object
        j : int
            ???
        dim_max : float
            ???
        nid_cp_cd : ???
            ???
        plot : bool; default=True
            ???
        xref_loads : bool; default=True
            ???

        Returns
        -------
        nid_to_pid_map  : dict
            node to property id map
            used to show SPC constraints (we don't want to show constraints
            on 456 DOFs)
        icase : int
            the result number
        cases : dict
            the GuiResult objects
        form : List[???, ???, ???]
            the Results sidebar data

        TDOO: Not quite done on:
               - ???
        """
        self.gui.isubcase_name_map = {1: ['Nastran', '']}
        grid = self.gui.grid

        nelements = self.nelements
        dim_array = np.zeros(nelements, 'int32')
        pids_array = np.zeros(nelements, 'int32')
        nnodes_array = np.zeros(nelements, 'int32')
        #mids = np.zeros(nelements, 'int32')
        #material_coord = np.zeros(nelements, 'int32')
        #min_interior_angle = np.zeros(nelements, 'float32')
        #max_interior_angle = np.zeros(nelements, 'float32')
        #dideal_theta = np.zeros(nelements, 'float32')
        #max_skew_angle = np.zeros(nelements, 'float32')
        #max_warp_angle = np.zeros(nelements, 'float32')
        #max_aspect_ratio = np.zeros(nelements, 'float32')
        #area = np.zeros(nelements, 'float32')
        #area_ratio = np.zeros(nelements, 'float32')
        #taper_ratio = np.zeros(nelements, 'float32')
        #min_edge_length = np.zeros(nelements, 'float32')
        if nelements == 0:
            return None

        if len(model.ctria3):
            model.ctria3.quality()
        #if len(model.tria6):
            #model.tria6.quality()
        #if len(model.quad4):
            #model.quad4.quality()
        #if len(model.cquad8):
            #model.cquad8.quality()
        #if len(model.cquad):
            #model.cquad.quality()

        nids_list = []
        unused_ieid = 0
        unused_cell_offset = 0
        dtype = get_numpy_idtype_for_vtk()

        nelements = self.nelements
        eids_array = np.zeros(nelements, dtype=dtype)
        cell_types_array = np.zeros(nelements, dtype=dtype)
        cell_offsets_array = np.zeros(nelements, dtype=dtype)

        results = {
            'pid' : pids_array,
            'eid' : eids_array,
            'nnodes' : nnodes_array,
            'dim' : dim_array,
        }

        cell_type_point = vtk.vtkVertex().GetCellType()
        cell_type_line = vtk.vtkLine().GetCellType()
        cell_type_tri3 = vtkTriangle().GetCellType()
        cell_type_tri6 = vtkQuadraticTriangle().GetCellType()
        cell_type_quad4 = vtkQuad().GetCellType()
        cell_type_quad8 = vtkQuadraticQuad().GetCellType()
        cell_type_tetra4 = vtkTetra().GetCellType()
        cell_type_tetra10 = vtkQuadraticTetra().GetCellType()
        cell_type_pyram5 = vtkPyramid().GetCellType()
        cell_type_pyram13 = vtk.vtkQuadraticPyramid().GetCellType()
        cell_type_penta6 = vtkWedge().GetCellType()
        cell_type_penta15 = vtkQuadraticWedge().GetCellType()
        cell_type_hexa8 = vtkHexahedron().GetCellType()
        cell_type_hexa20 = vtkQuadraticHexahedron().GetCellType()

        unused_all_eids = model.elements2.eids
        #print('type(eids) =', type(all_eids)) # list
        #print('all_eids =', all_eids)

        #ncelas1 = len(model.celas1)
        #ncelas2 = len(model.celas2)
        #ncelas3 = len(model.celas3)
        #ncelas4 = len(model.celas4)

        #ncdamp1 = len(model.cdamp1)
        #ncdamp2 = len(model.cdamp2)
        #ncdamp3 = len(model.cdamp3)
        #ncdamp4 = len(model.cdamp4)
        #ncdamp5 = len(model.cdamp5)

        #nconrod = len(model.conrod)
        #ncrod = len(model.crod)
        #nctube = len(model.ctube)

        #ncbar = len(model.cbar)
        #ncbeam = len(model.cbeam)

        #ncshear = len(model.cshear)

        #nctria3 = len(model.ctria3)
        #ncquad4 = len(model.cquad4)
        #nctria6 = len(model.ctria6)
        #ncquad8 = len(model.cquad8)
        #ncquad = len(model.cquad)

        #nctetra4 = len(model.ctetra4)
        #ncpenta6 = len(model.cpenta6)
        #nchexa8 = len(model.chexa8)
        #ncpyram5 = len(model.cpyram5)
        #nsolids = nctetra4 + ncpenta6 + nchexa8

        ieid0 = 0
        cell_offset0 = 0
        nids_list = []

        nodes = model.nodes

        ieid0, cell_offset0 = self._map_elements_vectorized_fill_spring(
            ieid0, cell_offset0,
            nodes, nids_list,
            eids_array, pids_array, nnodes_array, dim_array,
            cell_types_array, cell_offsets_array,
            model.celas1, cell_type_line, cell_type_point)
        ieid0, cell_offset0 = self._map_elements_vectorized_fill_spring(
            ieid0, cell_offset0,
            nodes, nids_list,
            eids_array, pids_array, nnodes_array, dim_array,
            cell_types_array, cell_offsets_array,
            model.celas2, cell_type_line, cell_type_point)
        ieid0, cell_offset0 = self._map_elements_vectorized_fill_spring(
            ieid0, cell_offset0,
            nodes, nids_list,
            eids_array, pids_array, nnodes_array, dim_array,
            cell_types_array, cell_offsets_array,
            model.celas3, cell_type_line, cell_type_point)
        ieid0, cell_offset0 = self._map_elements_vectorized_fill_spring(
            ieid0, cell_offset0,
            nodes, nids_list,
            eids_array, pids_array, nnodes_array, dim_array,
            cell_types_array, cell_offsets_array,
            model.celas4, cell_type_line, cell_type_point)

        ieid0, cell_offset0 = self._map_elements_vectorized_fill_spring(
            ieid0, cell_offset0,
            nodes, nids_list,
            eids_array, pids_array, nnodes_array, dim_array,
            cell_types_array, cell_offsets_array,
            model.cdamp1, cell_type_line, cell_type_point)
        ieid0, cell_offset0 = self._map_elements_vectorized_fill_spring(
            ieid0, cell_offset0,
            nodes, nids_list,
            eids_array, pids_array, nnodes_array, dim_array,
            cell_types_array, cell_offsets_array,
            model.cdamp2, cell_type_line, cell_type_point)
        ieid0, cell_offset0 = self._map_elements_vectorized_fill_spring(
            ieid0, cell_offset0,
            nodes, nids_list,
            eids_array, pids_array, nnodes_array, dim_array,
            cell_types_array, cell_offsets_array,
            model.cdamp3, cell_type_line, cell_type_point)
        ieid0, cell_offset0 = self._map_elements_vectorized_fill_spring(
            ieid0, cell_offset0,
            nodes, nids_list,
            eids_array, pids_array, nnodes_array, dim_array,
            cell_types_array, cell_offsets_array,
            model.cdamp4, cell_type_line, cell_type_point)

        ieid0, cell_offset0 = self._map_elements_vectorized_fill_spring(
            ieid0, cell_offset0,
            nodes, nids_list,
            eids_array, pids_array, nnodes_array, dim_array,
            cell_types_array, cell_offsets_array,
            model.cvisc, cell_type_line, cell_type_point)

        ieid0, cell_offset0 = self._map_elements_vectorized_fill(
            ieid0, cell_offset0,
            nodes, nids_list,
            eids_array, pids_array, nnodes_array, dim_array,
            cell_types_array, cell_offsets_array,
            model.plotel, cell_type_line, nnodesi=2, dimi=2)

        ieid0, cell_offset0 = self._map_elements_vectorized_fill_spring(
            ieid0, cell_offset0,
            nodes, nids_list,
            eids_array, pids_array, nnodes_array, dim_array,
            cell_types_array, cell_offsets_array,
            model.cbush, cell_type_line, cell_type_point)

        ieid0, cell_offset0 = self._map_elements_vectorized_fill(
            ieid0, cell_offset0,
            nodes, nids_list,
            eids_array, pids_array, nnodes_array, dim_array,
            cell_types_array, cell_offsets_array,
            model.conrod, cell_type_line, nnodesi=2, dimi=2)
        ieid0, cell_offset0 = self._map_elements_vectorized_fill(
            ieid0, cell_offset0,
            nodes, nids_list,
            eids_array, pids_array, nnodes_array, dim_array,
            cell_types_array, cell_offsets_array,
            model.crod, cell_type_line, nnodesi=2, dimi=2)
        ieid0, cell_offset0 = self._map_elements_vectorized_fill(
            ieid0, cell_offset0,
            nodes, nids_list,
            eids_array, pids_array, nnodes_array, dim_array,
            cell_types_array, cell_offsets_array,
            model.ctube, cell_type_line, nnodesi=2, dimi=2)

        ieid0, cell_offset0 = self._map_elements_vectorized_fill(
            ieid0, cell_offset0,
            nodes, nids_list,
            eids_array, pids_array, nnodes_array, dim_array,
            cell_types_array, cell_offsets_array,
            model.cbar, cell_type_line, nnodesi=2, dimi=1)

        ieid0, cell_offset0 = self._map_elements_vectorized_fill(
            ieid0, cell_offset0,
            nodes, nids_list,
            eids_array, pids_array, nnodes_array, dim_array,
            cell_types_array, cell_offsets_array,
            model.cbeam, cell_type_line, nnodesi=2, dimi=1)

        #model.cbend
        ieid0, cell_offset0 = self._map_elements_vectorized_fill(
            ieid0, cell_offset0,
            nodes, nids_list,
            eids_array, pids_array, nnodes_array, dim_array,
            cell_types_array, cell_offsets_array,
            model.cshear, cell_type_quad4, nnodesi=4, dimi=2)

        ieid0, cell_offset0 = self._map_elements_vectorized_fill(
            ieid0, cell_offset0,
            nodes, nids_list,
            eids_array, pids_array, nnodes_array, dim_array,
            cell_types_array, cell_offsets_array,
            model.ctria3, cell_type_tri3, nnodesi=3, dimi=2)
        ieid0, cell_offset0 = self._map_elements_vectorized_fill(
            ieid0, cell_offset0,
            nodes, nids_list,
            eids_array, pids_array, nnodes_array, dim_array,
            cell_types_array, cell_offsets_array,
            model.cquad4, cell_type_quad4, nnodesi=4, dimi=2)

        ieid0, cell_offset0 = self._map_elements_vectorized_fill(
            ieid0, cell_offset0,
            nodes, nids_list,
            eids_array, pids_array, nnodes_array, dim_array,
            cell_types_array, cell_offsets_array,
            model.ctria6, cell_type_tri6, nnodesi=6, dimi=2,
            allow0=True, cell_type_allow=cell_type_tri3)
        ieid0, cell_offset0 = self._map_elements_vectorized_fill(
            ieid0, cell_offset0,
            nodes, nids_list,
            eids_array, pids_array, nnodes_array, dim_array,
            cell_types_array, cell_offsets_array,
            model.cquad8, cell_type_quad8, nnodesi=8, dimi=2,
            allow0=True, cell_type_allow=cell_type_quad8)

        ieid0, cell_offset0 = self._map_elements_vectorized_fill(
            ieid0, cell_offset0,
            nodes, nids_list,
            eids_array, pids_array, nnodes_array, dim_array,
            cell_types_array, cell_offsets_array,
            model.ctriar, cell_type_tri6, nnodesi=6, dimi=2,
            allow0=True, cell_type_allow=cell_type_tri3)
        ieid0, cell_offset0 = self._map_elements_vectorized_fill(
            ieid0, cell_offset0,
            nodes, nids_list,
            eids_array, pids_array, nnodes_array, dim_array,
            cell_types_array, cell_offsets_array,
            model.cquadr, cell_type_tri6, nnodesi=6, dimi=2,
            allow0=True, cell_type_allow=cell_type_quad4)

        ieid0, cell_offset0 = self._map_elements_vectorized_fill(
            ieid0, cell_offset0,
            nodes, nids_list,
            eids_array, pids_array, nnodes_array, dim_array,
            cell_types_array, cell_offsets_array,
            model.ctetra4, cell_type_tetra4, nnodesi=4, dimi=3)
        ieid0, cell_offset0 = self._map_elements_vectorized_fill(
            ieid0, cell_offset0,
            nodes, nids_list,
            eids_array, pids_array, nnodes_array, dim_array,
            cell_types_array, cell_offsets_array,
            model.cpenta6, cell_type_penta6, nnodesi=8, dimi=3)
        ieid0, cell_offset0 = self._map_elements_vectorized_fill(
            ieid0, cell_offset0,
            nodes, nids_list,
            eids_array, pids_array, nnodes_array, dim_array,
            cell_types_array, cell_offsets_array,
            model.chexa8, cell_type_hexa8, nnodesi=8, dimi=3)
        ieid0, cell_offset0 = self._map_elements_vectorized_fill(
            ieid0, cell_offset0,
            nodes, nids_list,
            eids_array, pids_array, nnodes_array, dim_array,
            cell_types_array, cell_offsets_array,
            model.cpyram5, cell_type_pyram5, nnodesi=8, dimi=3)

        ieid0, cell_offset0 = self._map_elements_vectorized_fill(
            ieid0, cell_offset0,
            nodes, nids_list,
            eids_array, pids_array, nnodes_array, dim_array,
            cell_types_array, cell_offsets_array,
            model.ctetra10, cell_type_tetra10, nnodesi=4, dimi=3,
            allow0=True, cell_type_allow=cell_type_tetra4)
        ieid0, cell_offset0 = self._map_elements_vectorized_fill(
            ieid0, cell_offset0,
            nodes, nids_list,
            eids_array, pids_array, nnodes_array, dim_array,
            cell_types_array, cell_offsets_array,
            model.cpenta15, cell_type_penta15, nnodesi=8, dimi=3,
            allow0=True, cell_type_allow=cell_type_penta6)
        ieid0, cell_offset0 = self._map_elements_vectorized_fill(
            ieid0, cell_offset0,
            nodes, nids_list,
            eids_array, pids_array, nnodes_array, dim_array,
            cell_types_array, cell_offsets_array,
            model.chexa20, cell_type_hexa20, nnodesi=8, dimi=3,
            allow0=True, cell_type_allow=cell_type_hexa8)
        ieid0, cell_offset0 = self._map_elements_vectorized_fill(
            ieid0, cell_offset0,
            nodes, nids_list,
            eids_array, pids_array, nnodes_array, dim_array,
            cell_types_array, cell_offsets_array,
            model.cpyram13, cell_type_pyram13, nnodesi=8, dimi=3,
            allow0=True, cell_type_allow=cell_type_pyram5)

        # model.chbdyg
        # model.chbdye
        # model.chbdyp

        if cell_types_array.min() == 0:

            # all the non-elemental cards should be listed
            # it's not hugely important, but it cleans up dev error messages
            skip_cards = [
                'CONM2',
                #'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4', 'PLOTEL',
                'PARAM',
                #'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CVISC',
                'TABLEM1', 'TABLEM2', 'TABLEM3', 'TABLEM4',
                'TABLED1', 'TABLED2', 'TABLED3', 'TABLED4', 'TABLEST',
                'MAT1', 'MAT2', 'MAT4', 'MAT5', 'MAT8', 'MAT9', 'MAT10',
                'MATT1', 'MATT2', 'MATT8',
                'MATS1',

                'PLOAD', 'PLOAD1', 'PLOAD2', 'FORCE', 'PLOAD4', 'LOAD',
                'MAT1', 'PSHEAR', 'PSHELL', 'PTUBE', 'PDAMP',
                'PELAST', 'PBEND', 'PBEAM', 'PCOMP', 'PCOMPG', 'PBAR', 'PSOLID',
                'PROD', 'PELAS', 'PVISC', 'PBUSH1D', 'PBUSH2D',
                #'EPOINT',
                #'CQUADR', 'CTRIAR', 'SPOINT',
                #'CQUAD8', 'CTRIA6',
                'ENDDATA',
                'CORD2R', 'CORD2C', 'CORD2S', 'CORD1R', 'CORD1C', 'CORD1S',
                'GRID', 'SPOINT', 'EPOINT', 'TF',

                'RFORCE', 'RFORCE1', 'RFORCE2', 'FORCE', 'FORCE1', 'FORCE2',
                'MOMENT', 'MOMENT1', 'MOMENT2', 'PLOAD', 'PLOAD1', 'PLOAD2', 'PLOAD4',
                'LOAD', 'TLOAD1', 'TLOAD2', 'DLOAD', 'LSEQ', 'DAREA',
                'RLOAD1', 'RLOAD2',

                'SUPORT', 'SUPORT1', 'MPC', 'MPCADD', 'RBE1', 'RBE2', 'RBE3', 'RBAR', 'RCROSS',
                'SPCADD', 'SPC', 'SPC1', 'SPCD', 'SPCAX', 'DMIG', 'DMI', 'DMIJ', 'DMIJI', 'DMIK',

                'AELIST', 'AELINK', 'AESURF', 'AESURFS', 'AERO', 'AEROS', 'TRIM',
                'FLUTTER', 'DIVERG',
                'CAERO1', 'CAERO2', 'CAERO3', 'CAERO4', 'CAERO5',
                'PAERO1', 'PAERO2', 'PAERO3', 'PAERO4', 'PAERO5',
                'SPLINE1', 'SPLINE2', 'SPLINE3', 'SPLINE4', 'SPLINE5', 'SPLINE6', 'SPLINE7',
            ]
            potential_elements_found = [key for key in model.card_count if key not in skip_cards]
            etypes = [
                'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
                'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5', 'CVISC',
                'CBUSH', 'CBUSH1D', 'CBUSH2D',
                'CONROD', 'CROD', 'CTUBE', 'PLOTEL',
                'CBAR', 'CBEAM', 'CBEND',
                'CSHEAR',
                'CTRIA3', 'CQUAD4', 'CTRIA6', 'CQUAD8', 'CTRIAR', 'CQUADR',
                'CTETRA', 'CPENTA', 'CHEXA', 'CPYRAM',
                'CHBDYG', 'CHBDYE', 'CHBDYP',
            ]
            for key in potential_elements_found:
                if key not in etypes:
                    self.log.warning('is %s an element?' % key)

            msg = (
                'Cell Type is not defined (cell_type=0).\n'
                '  cell_types_array = %s\n'
                '  potential_elements_found=[%s]\n'
                '  nelements=%s\n\n'
                '%s\n\n' % (
                    cell_types_array,
                    ', '.join(potential_elements_found),
                    len(cell_types_array),
                    '', #str(model.elements2),
                )
            )
            print(str(model.elements2))
            #msg += model.get_bdf_stats()
            raise RuntimeError(msg)

        deep = 1
        if len(nids_list) == 1:
            nids_array = nids_list[0].ravel()
        else:
            #raise NotImplementedError(len(nids_list))
            nids_array = np.hstack([nid_list.flatten() for nid_list in nids_list])
            #nids_array = np.array(nids_list, dtype=dtype)

        #-----------------------------------------------------------------
        # saving some data members
        self.gui.element_ids = eids_array

        #-----------------------------------------------------------------
        # build the grid

        #self.log.info('nids_array = %s' % nids_array)
        #self.log.info('cell_offsets_array = %s' % cell_offsets_array)
        #self.log.info('cell_types_array = %s' % cell_types_array)

        # Create the array of cells
        #print('nids_array =', nids_array)
        cells_id_type = numpy_to_vtkIdTypeArray(nids_array, deep=1)
        vtk_cells = vtk.vtkCellArray()
        vtk_cells.SetCells(nelements, cells_id_type)

        # Cell types
        vtk_cell_types = numpy_to_vtk(
            cell_types_array, deep=deep,
            array_type=vtk.vtkUnsignedCharArray().GetDataType())

        vtk_cell_offsets = numpy_to_vtk(cell_offsets_array, deep=deep,
                                        array_type=vtk.VTK_ID_TYPE)

        grid = self.gui.grid
        #grid = vtk.vtkUnstructuredGrid()
        grid.SetCells(vtk_cell_types, vtk_cell_offsets, vtk_cells)
        return results

    def _map_elements_vectorized_fill_spring(self,
                                             ieid0, cell_offset0,
                                             nodes, nids_list,
                                             eids_array, pids_array, nnodes_array, dim_array,
                                             cell_types_array, cell_offsets_array,
                                             model_obj, cell_type_line, cell_type_point):
        """helper method for ``map_elements_vectorized``"""
        nelems = len(model_obj)
        if nelems:
            ieid0_in = ieid0
            dimi = 1
            nnodesi = 2
            elem = model_obj
            elem.make_current()
            self.log.info('%s; nelem=%s nnodes=%s' % (
                model_obj.card_name, nelems, model_obj.nids.shape[1]))

            ieid = np.arange(ieid0, ieid0 + nelems)
            eid_type = np.full(nelems, cell_type_line, dtype='int32')
            dim = np.full(nelems, dimi, dtype='int32')
            nnodes = np.full((nelems, 1), nnodesi, dtype='int32')
            dim_array[ieid] = dim

            eid = elem.eid
            eids_array[ieid] = eid

            nids = elem.nids
            inonzero = np.where(nids[:, 1] > 0)[0]
            izero = np.where(nids[:, 1] == 0)[0]

            if len(inonzero):
                cell_types_array[ieid[inonzero]] = cell_type_line

            nzero = len(izero)
            if nzero:
                #raise NotImplementedError('izero')
                eid_type[izero] = cell_type_point
                nnodes[izero] = 1
                cell_types_array[ieid[izero]] = cell_type_point
            assert len(ieid) == len(nnodes), 'len(ieid)=%s len(nnodes)=%s' % (len(ieid), nnodes)
            nnodes_array[ieid] = nnodes.ravel()

            if hasattr(elem, 'pid'):
                pid = elem.pid
                pids_array[ieid] = pid
            inids = nodes.get_node_index(nids, allow0=True)

            nnodes_inids = np.hstack([nnodes, inids])#.ravel()
            #print("  nnodes_inids=", nnodes_inids)
            if nzero:
                nnodes_inids_ravel = nnodes_inids.flatten()
                index_minus1 = np.where(nnodes_inids_ravel != -1)[0]
                nnodes_inids_cleaned = nnodes_inids_ravel[index_minus1]
                assert -1 not in nnodes_inids_cleaned, nnodes_inids_cleaned
                nids_list.append(nnodes_inids_cleaned)
            else:
                nids_list.append(nnodes_inids.ravel())

            nnodes[0] = -1
            cumsum = cell_offset0 + np.cumsum(nnodes + 1)
            #print("nnodes =", nnodes)
            #print('offset0=%s cumsum=%s' % (cell_offset0, cumsum))
            self.log.debug('  ieid = %s' % ieid)
            assert len(ieid) == nelems
            assert len(ieid) == len(cumsum)
            cell_offsets_array[ieid] = cumsum
            ieid0 += nelems
            self.log.debug("  ieid0_in=%s ieid0_out=%s" % (ieid0_in, ieid0))
            cell_offset0 += cumsum[-1]
            #print('nids_list[%s] = %s' % (model_obj.card_name, nids_list[-1]))
            #print('cell_offsets_array =', cell_offsets_array)
        return ieid0, cell_offset0

    def _map_elements_vectorized_fill(self,
                                      ieid0, cell_offset0,
                                      nodes, nids_list,
                                      eids_array, pids_array, nnodes_array, dim_array,
                                      cell_types_array, cell_offsets_array,
                                      model_obj, cell_type, nnodesi=None, dimi=None,
                                      allow0=False, cell_type_allow=None):
        """helper method for ``map_elements_vectorized``"""
        assert nnodesi is not None
        assert dimi is not None

        nelems = len(model_obj)
        #print('nelements =', nelems)
        if nelems:
            ieid0_in = ieid0
            elem = model_obj
            elem.make_current()
            try:
                nnodes = model_obj.nids.shape[1]
            except IndexError:
                self.log.info('%s; nelem=%s nnodes=???' % (
                    model_obj.card_name, nelems))
                raise

            self.log.info('%s; nelem=%s nnodes=%s' % (
                model_obj.card_name, nelems, nnodes))
            eid = elem.eid
            ieid = np.arange(ieid0, ieid0 + nelems)
            dim = np.full(nelems, dimi, dtype='int32')
            dim_array[ieid] = dim

            nids = elem.nids
            #print('nids =', nids)
            if hasattr(elem, 'pid'):
                pid = elem.pid
                pids_array[ieid] = pid

            try:
                inids = nodes.get_node_index(nids, allow0=allow0)
            except:
                bdf_file = StringIO()
                model_obj.write_card(bdf_file=bdf_file)
                print(bdf_file.getvalue())
                raise
            if inids.min() == -1:
                nrequired = elem.nrequired
                #print('nrequired =', nrequired)
                #print('nids.shape =', nids.shape)
                #print('inids.shape =', inids.shape)
                inids = inids[:, :nrequired]
                nnodesi = nrequired
                assert cell_type_allow is not None
                unused_eid_type = np.full(nelems, cell_type_allow, dtype='int32')
            else:
                unused_eid_type = np.full(nelems, cell_type, dtype='int32')

            nnodes = np.full((nelems, 1), nnodesi, dtype='int32')
            if len(ieid) != len(nnodes):
                msg = 'len(ieid)=%s len(nnodes)=%s' % (len(ieid), len(nnodes))
                raise RuntimeError(msg)

            # nnodes is a 2D column matrix and nnodes_array is a 1D array,
            # so we ravel it
            nnodes_array[ieid] = nnodes.ravel()
            assert -1 not in inids, inids
            nnodes_inids = np.hstack([nnodes, inids])
            nids_list.append(nnodes_inids)
            self.log.debug('  ieid = %s' % ieid)
            assert len(ieid) == nelems

            nnodes[0] = -1
            cumsum = cell_offset0 + np.cumsum(nnodes + 1)
            assert len(ieid) == len(cumsum)
            try:
                eids_array[ieid] = eid
            except IndexError:
                # this hits when an element doesn't correctly implement make_current()
                self.log.error('nelems = %s' % nelems)
                self.log.error('card_count = %s' % model_obj.model.card_count)
                self.log.error('eids_array = %s' % eids_array)
                self.log.error('eids_array.shape = %s' % str(eids_array.shape))
                self.log.error('ieid.max() = %s' % ieid.max())
                raise
            cell_types_array[ieid] = cell_type
            cell_offsets_array[ieid] = cumsum
            ieid0 += nelems
            cell_offset0 += cumsum[-1]
            #print('nids_list[%s] = %s' % (model_obj.card_name, nids_list[-1]))
            #print('cell_offsets_array =', cell_offsets_array)
            self.log.debug("  ieid0_in=%s ieid0_out=%s" % (ieid0_in, ieid0))
        return ieid0, cell_offset0

    def _get_model_vectorized(self, bdf_filename):
        """Loads the BDF/OP2 geometry"""
        ext = os.path.splitext(bdf_filename)[1].lower()
        punch = False
        if ext == '.pch':
            punch = True

        self.model_type = 'nastran'
        if ext == '.op2':
            from pyNastran.dev.bdf_vectorized2.op2_geom_vectorized import (
                OP2Geom as OP2Geom_)
            model = OP2Geom_(make_geom=True, debug=False, log=self.log,
                             debug_file=None)
            model.clear_results()
            model.read_op2(op2_filename=bdf_filename)
        else:  # read the bdf/punch
            from pyNastran.dev.bdf_vectorized2.bdf_vectorized import BDF as BDF_
            model = BDF_(log=self.log, debug=True)

            # static_elements.bdf
            #skip_cards = [
                #'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4', 'PLOTEL', 'PARAM',
                #'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CVISC',
                #'TABLEM1', 'TABLEM2', 'TABLEM3', 'TABLEM4',
                #'TABLED1', 'TABLED2', 'TABLED3', 'TABLED4',
                #'PLOAD', 'PLOAD1', 'PLOAD2', 'FORCE', 'PLOAD4', 'LOAD',
                #'SPCADD', 'MAT1', 'PSHEAR', 'PSHELL', 'PTUBE', 'PDAMP',
                #'SPC1', 'CONM2', 'PELAST', 'PBEND', 'PBEAM', 'PCOMP', 'PCOMPG', 'PBAR', 'PSOLID',
                #'PBUSH1D',
                #'EPOINT',
                #'CQUADR', 'CTRIAR', 'SPOINT', 'PROD', 'PELAS', 'PVISC',
                #'CQUAD8', 'CTRIA6',
                #]
            #model.disable_cards(skip_cards)
            model.read_bdf(bdf_filename,
                           punch=punch, xref=False,
                           validate=True)
            #print(list(key for key in model.card_count.keys() if key not in skip_cards))

        #xref_loads = False
        #xref_aero = len(model.caeros) > 0
        #model.cross_reference(
            #xref=True,
            #xref_nodes=True,
            #xref_elements=False,
            #xref_nodes_with_elements=False,
            #xref_properties=True,
            #xref_masses=True,
            #xref_materials=False,
            #xref_loads=xref_loads,
            #xref_constraints=False,
            #xref_optimization=False,
            #xref_aero=False,
            #xref_sets=False,
        #)
        return model

    def _points_to_vtkpoints_coords(self, model, xyz_cid0):
        """
        helper method for:
         - load_nastran_geometry_unvectorized
         - load_nastran_geometry_vectorized
        """
        points = numpy_to_vtk_points(xyz_cid0)
        self.gui.grid.SetPoints(points)

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

        self.gui.log_info("xmin=%s xmax=%s dx=%s" % (xmin, xmax, xmax-xmin))
        self.gui.log_info("ymin=%s ymax=%s dy=%s" % (ymin, ymax, ymax-ymin))
        self.gui.log_info("zmin=%s zmax=%s dz=%s" % (zmin, zmax, zmax-zmin))
        return dim_max

    def load_nastran_geometry_unvectorized(self, bdf_filename, plot=True):
        """
        The entry point for Nastran geometry loading.

        Parameters
        ----------
        bdf_filename : str
            the Nastran filename to load
        plot : bool; default=True
            should the model be generated or should we wait until
            after the results are loaded
        """
        reset_labels = True
        if plot:
            self.gui.scalarBar.VisibilityOff()
            self.gui.scalarBar.Modified()

        xref_loads = True # should be True
        model, xref_nodes = self._get_model_unvectorized(bdf_filename, xref_loads=xref_loads)

        nnodes = len(model.nodes)
        nspoints = len(model.spoints)
        nepoints = len(model.epoints)
        ngridb = len(model.gridb)
        ncaero_cards = len(model.caeros)

        ngui_nodes = nnodes + nspoints + nepoints + ngridb
        if ngui_nodes + ncaero_cards == 0:
            msg = 'nnodes + nspoints + nepoints = 0\n'
            msg += 'card_count = %r' % str(model.card_count)
            raise NoGeometry(msg)

        nelements = len(model.elements)
        nmasses = len(model.masses)
        nplotels = len(model.plotels)
        nrigid = len(model.rigid_elements)
        #nmpc = len(model.mpcs)  # really should only be allowed if we have it in a subcase
        if nelements + nmasses + ncaero_cards + nplotels + nrigid == 0:
            msg = 'nelements + nmasses + ncaero_cards + nplotels + nrigid = 0\n'
            msg += 'card_count = %r' % str(model.card_count)
            raise NoGeometry(msg)

        self.nnodes = ngui_nodes
        self.nelements = nelements  # approximate...

        out = self.make_caeros(model)
        (has_caero, caero_points, ncaeros, ncaeros_sub, ncaeros_cs,
         ncaeros_points, ncaero_sub_points,
         has_control_surface, box_id_to_caero_element_map, cs_box_ids) = out
        self.has_caero = has_caero

        self.gui.log_info("nnodes=%i nelements=%i" % (self.nnodes, self.nelements))
        msg = model.get_bdf_stats(return_type='string')
        self.gui.log_debug(msg)
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
        # CMASS3, CMASS4 are applied to SPOINTs

        update_conm2s_function = None
        if nconm2 > 0:
            def update_conm2s_function(unused_nid_map, unused_ugrid, points, nodes):
                j2 = 0
                mass_grid = self.gui.alt_grids['conm2']
                for unused_eid, element in sorted(model.masses.items()):
                    if isinstance(element, CONM2):
                        nid = element.nid
                        inid = np.searchsorted(self.node_ids, nid)
                        xyz_nid = nodes[inid, :]
                        centroid = element.offset(xyz_nid)
                        points.SetPoint(j2, *centroid)

                    elif element.type in ['CMASS1', 'CMASS2']:
                        n1, n2 = element.nodes
                        factor = 0.
                        if element.nodes[0] is not None:
                            inid = np.searchsorted(self.node_ids, n1)
                            p1 = nodes[inid, :]
                            factor += 1.
                        if element.nodes[1] is not None:
                            inid = np.searchsorted(self.node_ids, n2)
                            p2 = nodes[inid, :]
                            factor += 1.
                        centroid = (p1 + p2) / factor
                        points.SetPoint(j2, *centroid)

                        elem = vtk.vtkVertex()
                        elem.GetPointIds().SetId(0, j2)
                        mass_grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                    else:
                        continue
                        #self.gui.log_info("skipping %s" % element.type)
                    j2 += 1
                return

            self.gui.create_alternate_vtk_grid(
                'conm2', color=ORANGE_FLOAT, line_width=5, opacity=1., point_size=4,
                follower_function=update_conm2s_function,
                representation='point')

        # Allocate grids
        self.gui.grid.Allocate(self.nelements, 1000)
        self._create_caero_actors(ncaeros, ncaeros_sub, ncaeros_cs, has_control_surface)
        if nconm2 > 0:
            self.gui.alt_grids['conm2'].Allocate(nconm2, 1000)

        if self.save_data:
            self.model = model

        #-----------------------------------------------------------------------
        # nodes/coords

        #print('get_xyz_in_coord')
        dim_max = 1.0
        xyz_cid0 = None
        nid_cp_cd = None
        if self.gui.nnodes:
            xyz_cid0, nid_cp_cd = self.get_xyz_in_coord(model, cid=0, fdtype='float32')
            dim_max = self._points_to_vtkpoints_coords(model, xyz_cid0)
        #-----------------------------------------------------------------------

        j = 0
        nid_to_pid_map, icase, cases, form = self.map_elements(
            xyz_cid0, nid_cp_cd, self.gui.nid_map, model, j, dim_max,
            plot=plot, xref_loads=xref_loads)

        # fill grids
        zfighting_offset0 = 0.001
        zfighting_offset = zfighting_offset0
        self._create_splines(model, box_id_to_caero_element_map, caero_points)
        if 'caero' in self.gui.alt_grids:
            self.set_caero_grid(ncaeros_points, model, j=0)
            self.set_caero_subpanel_grid(ncaero_sub_points, model, j=0)
            if has_control_surface:
                cs_name = 'caero_control_surfaces'
                self.set_caero_control_surface_grid(
                    cs_name, cs_box_ids[cs_name],
                    box_id_to_caero_element_map, caero_points,
                    zfighting_offset=zfighting_offset)
                zfighting_offset += zfighting_offset0

                # sort the control surfaces
                labels_to_aesurfs = {aesurf.label: aesurf for aesurf in model.aesurf.values()}
                if len(labels_to_aesurfs) != len(model.aesurf):
                    msg = (
                        'Expected same number of label->aesurf as aid->aesurf\n'
                        'labels_to_aesurfs = %r\n'
                        'model.aesurf = %r\n' % (labels_to_aesurfs, model.aesurf))
                    raise RuntimeError(msg)

                for unused_label, aesurf in sorted(labels_to_aesurfs.items()):
                    reset_labels = False
                    cs_name = '%s_control_surface' % aesurf.label
                    self.set_caero_control_surface_grid(
                        cs_name, cs_box_ids[cs_name],
                        box_id_to_caero_element_map, caero_points, note=aesurf.label,
                        zfighting_offset=zfighting_offset)
                    zfighting_offset += zfighting_offset0

        if nconm2 > 0 and xref_nodes:
            self._set_conm_grid(nconm2, model)


        make_spc_mpc_supports = True
        geometry_names = []
        if make_spc_mpc_supports and xref_nodes:
            geometry_names = self.set_spc_mpc_suport_grid(model, nid_to_pid_map)

        if xref_nodes:
            icase = self._fill_bar_yz(dim_max, model, icase, cases, form)
        assert icase is not None

        #------------------------------------------------------------
        #print('dependent_nodes =', self.dependents_nodes)
        form0 = form[2]
        assert icase is not None
        nsubcases = len(model.subcases)
        for subcase_idi, subcase in sorted(model.subcases.items()):
            if not xref_nodes:
                continue

            subcase_id = subcase_idi
            if subcase_id == 0 and nsubcases == 1:
                subcase_id = 1
            elif subcase_id == 0:
                continue
            self.gui.log_debug('NastranIOv subcase_id = %s' % subcase_id)

            subtitle = ''
            if 'SUBTITLE' in subcase:
                subtitle, options = subcase.get_parameter('SUBTITLE')
                del options

            load_str = 'Load Case=%i' % subcase_id if subtitle == '' else 'Load Case=%i; %s' % (
                subcase_id, subtitle)
            formi = (load_str, None, [])
            formii = formi[2]

            assert icase is not None
            if self.normals is not None and self.plot_applied_loads:
                icase = self._plot_applied_loads(
                    model, cases, formii, icase, subcase_idi, xref_loads=xref_loads)
                #plot_pressures = False
                plot_pressures = True
            else:
                plot_pressures = True

            if plot_pressures: # and self._plot_pressures:
                try:
                    icase = self._plot_pressures(
                        model, cases, formii, icase, subcase_idi)
                except KeyError:
                    s = StringIO()
                    traceback.print_exc(file=s)
                    sout = s.getvalue()
                    self.gui.log_error(sout)
                    print(sout)

            if len(formii):
                form0.append(formi)

        name = 'main_copy'
        self.gui.duplicate_alternate_vtk_grid(
            name, 'main', color=(0., 0., 0.), line_width=5,
            opacity=0.1, is_visible=False)

        #------------------------------------------------------------
        # add alternate actors
        self.gui._add_alt_actors(self.gui.alt_grids)

        # set default representation
        self._set_caero_representation(has_control_surface)

        for grid_name in geometry_names:
            if grid_name in self.gui.geometry_actors:
                self.gui.geometry_actors[grid_name].Modified()

        #self.grid_mapper.SetResolveCoincidentTopologyToPolygonOffset()
        if plot:
            self.gui._finish_results_io2([form], cases, reset_labels=reset_labels)
        else:
            self.gui._set_results([form], cases)

    def _create_caero_actors(self, ncaeros, ncaeros_sub, ncaeros_cs, has_control_surface):
        """
        This just creates the following actors.  It does not fill them.
        These include:
         - caero
         - caero_subpanels
         - caero_control_surfaces
        """
        if self.has_caero:
            if 'caero' not in self.gui.alt_grids:
                self.gui.create_alternate_vtk_grid(
                    'caero', color=YELLOW_FLOAT, line_width=3, opacity=1.0,
                    representation='toggle', is_visible=True, is_pickable=False)
            if 'caero_subpanels' not in self.gui.alt_grids:
                self.gui.create_alternate_vtk_grid(
                    'caero_subpanels', color=YELLOW_FLOAT, line_width=3, opacity=1.0,
                    representation='toggle', is_visible=False, is_pickable=False)

            self.gui.alt_grids['caero'].Allocate(ncaeros, 1000)
            self.gui.alt_grids['caero_subpanels'].Allocate(ncaeros_sub, 1000)
            if has_control_surface:
                self.gui.alt_grids['caero_control_surfaces'].Allocate(ncaeros_cs, 1000)


    def _set_caero_representation(self, has_control_surface):
        """
        Parameters
        ----------
        has_control_surface : bool
            is there a control surface
        """
        if 'caero_control_surfaces' in self.gui.geometry_actors:
            self.geometry_properties['caero_control_surfaces'].opacity = 0.5

            if 'caero' not in self.gui.geometry_actors:
                return
            self.gui.geometry_actors['caero'].Modified()
            self.gui.geometry_actors['caero_subpanels'].Modified()
            if has_control_surface:
                self.gui.geometry_actors['caero_control_surfaces'].Modified()
            if hasattr(self.gui.geometry_actors['caero'], 'Update'):
                self.gui.geometry_actors['caero'].Update()
            if hasattr(self.gui.geometry_actors['caero_subpanels'], 'Update'):
                self.gui.geometry_actors['caero_subpanels'].Update()
            if has_control_surface and hasattr(self.gui.geometry_actors['caero_subpanels'], 'Update'):
                self.gui.geometry_actors['caero_control_surfaces'].Update()

    def _create_splines(self, model, box_id_to_caero_element_map, caero_points):
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
        stored_msg = []
        if model.splines:
            # 0 - caero / caero_subpanel
            # 1 - control surface
            # 3/5/7/... - spline points
            # 2/4/6/... - spline panels
            iaero = 2
            for spline_id, spline in sorted(model.splines.items()):
                setg_ref = spline.setg_ref
                if setg_ref is None:
                    msg = 'error cross referencing SPLINE:\n%s' % spline.rstrip()
                    #n, filename = log_properties(1)
                    #print(filename, n)
                    #stored_msg.append(msg)
                    self.log.error(msg)
                    #raise RuntimeError(msg)
                    continue
                else:
                    structure_points = setg_ref.get_ids()

                try:
                    aero_box_ids = spline.aero_element_ids
                except:
                    print(spline.object_attributes())
                    print(spline.object_methods())
                    raise
                if spline.type != 'SPLINE3_ZAERO':
                    assert len(aero_box_ids) > 0, spline
                # the control surfaces all lie perfectly on top of each other
                # such that we have z fighting, so based on the aero index,
                # we calculate a z offset.
                zfighting_offset = 0.0001 * (iaero + 1)
                grid_name = 'spline_%s_structure_points' % spline_id
                self.gui.create_alternate_vtk_grid(
                    grid_name, color=BLUE_FLOAT, opacity=1.0, point_size=5,
                    representation='point', is_visible=False)
                msg = ', which is required by %r' % grid_name
                stored_msgi = self._add_nastran_nodes_to_grid(
                    grid_name, structure_points, model, msg, store_msg=True)

                zfighting_offset = 0.0001 * (iaero + 2)
                grid_name = 'spline_%s_boxes' % spline_id
                self.gui.create_alternate_vtk_grid(
                    grid_name, color=BLUE_FLOAT, opacity=0.3,
                    line_width=4,
                    representation='toggle', is_visible=False)
                stored_msgi2 = self.set_caero_control_surface_grid(
                    grid_name, aero_box_ids,
                    box_id_to_caero_element_map, caero_points,
                    zfighting_offset=zfighting_offset, store_msg=True)
                iaero += 2
                if stored_msgi:
                    stored_msg.append(stored_msgi)
                if stored_msgi2:
                    stored_msg.append(stored_msgi2)
        if stored_msg:
            model.log.warning('\n' + '\n'.join(stored_msg))

    def make_caeros(self, model):
        """
        Creates the CAERO panel inputs including:
         - caero
         - caero_subpanels
         - caero_control_surfaces
         - N control surfaces

        Parameters
        ----------
        model : BDF()
            the bdf model

        Returns
        -------
        caero_points : (N_aero_points, 3) float ndarray
            the xyz points for the aero panels
            N_aero_points can be 0
        ncaeros : int
            the number of aero sub-panels?
        ncaeros_sub : int
            ???
        ncaeros_cs : int
            ???
        ncaeros_points : int
            number of points for the caero coarse grid
        ncaero_sub_points : int
            number of points for the caero fine/subpanel grid
        has_control_surface : bool
            is there a control surface
        box_id_to_caero_element_map : dict[box_id] = box_index
            used to map the CAEROx box id to index in the ???
            (aero panel elements) array, which will be used with
            cs_box_ids
        cs_box_ids : dict[control_surface_name] : List[panel ids]
            list of panels used by each aero panel
        """
        has_caero = False
        ncaeros = 0
        ncaeros_sub = 0
        ncaeros_cs = 0
        ncaeros_points = 0
        ncaero_sub_points = 0
        has_control_surface = False
        box_id_to_caero_element_map = {}
        cs_box_ids = defaultdict(list)

        # when caeros is empty, SPLINEx/AESURF cannot be defined
        if len(model.caeros) == 0:
            #print('returning early')
            caero_points = np.empty((0, 3))
            out = (
                has_caero, caero_points, ncaeros, ncaeros_sub, ncaeros_cs,
                ncaeros_points, ncaero_sub_points,
                has_control_surface, box_id_to_caero_element_map, cs_box_ids,
            )
            return out

        # count caeros
        # sorting doesn't matter here because we're just trying to size the array
        for caero in model.caeros.values():
            if hasattr(caero, 'panel_points_elements'):
                npoints, ncelements = caero.get_npanel_points_elements()
                ncaeros_sub += npoints
                ncaero_sub_points += ncelements
            elif isinstance(caero, (CAERO2, BODY7)):
                pass
            else:  # pragma: no cover
                msg = '%r doesnt support panel_points_elements\n%s' % (caero.type, caero.rstrip())
                raise NotImplementedError(msg)

        for unused_eid, caero in sorted(model.caeros.items()):
            if isinstance(caero, (CAERO1, CAERO3, CAERO4, CAERO5, CAERO7)):
                ncaeros_points += 4
                ncaeros += 1
            elif isinstance(caero, (CAERO2, BODY7)):
                points, elems = caero.get_points_elements_3d()
                if points is None:
                    continue
                ncaeros_points += points.shape[0]
                ncaeros += elems.shape[0]
            else:  # pragma: no cover
                msg = '%r doesnt support panel counter\n%s' % (caero.type, caero.rstrip())
                raise NotImplementedError(msg)

        num_prev = 0
        ncaeros_sub = 0
        if model.caeros:
            caero_points = []
            for unused_eid, caero in sorted(model.caeros.items()):
                if caero.type in ['CAERO1', 'CAERO4', 'CAERO7']:
                    box_ids = caero.box_ids
                    nboxes = len(box_ids.ravel())
                    if nboxes > 1000:
                        print('skipping nboxes=%s for:\n%s' % (nboxes, str(caero)))
                        continue

                    ncaeros_sub += 1
                    pointsi, elementsi = caero.panel_points_elements()
                    caero_points.append(pointsi)

                    for i, box_id in enumerate(caero.box_ids.flat):
                        box_id_to_caero_element_map[box_id] = elementsi[i, :] + num_prev
                    num_prev += pointsi.shape[0]
                elif caero.type in ('CAERO2', 'BODY7'):
                    pass
                else:
                    print('caero\n%s' % caero)
            if ncaeros_sub:
                caero_points = np.vstack(caero_points)
            has_caero = True

        if ncaeros_sub == 0:
            caero_points = np.empty((0, 3))

        # check for any control surfcaes
        if model.aesurf:
            has_control_surface = True
            #ncaero_cs_points = 0
            if 'caero_control_surfaces' not in self.gui.alt_grids:
                self.gui.create_alternate_vtk_grid(
                    'caero_control_surfaces', color=PINK_FLOAT, line_width=5, opacity=1.0,
                    representation='surface', is_visible=False)

            # sort the control surfaces
            labels_to_aesurfs = {aesurf.label: aesurf for aesurf in model.aesurf.values()}
            if len(labels_to_aesurfs) != len(model.aesurf):
                msg = (
                    'Expected same number of label->aesurf as aid->aesurf\n'
                    'labels_to_aesurfs = %r\n'
                    'model.aesurf = %r\n' % (labels_to_aesurfs, model.aesurf))
                raise RuntimeError(msg)

            for unused_label, aesurf in sorted(model.aesurf.items()):
                if aesurf.type == 'AESURFZ':
                    aero_element_ids = aesurf.aero_element_ids
                    ncaeros_cs += len(aero_element_ids)

                    cs_name = '%s_control_surface' % aesurf.label
                    if cs_name not in self.alt_grids:
                        self.gui.create_alternate_vtk_grid(
                            cs_name, color=PINK_FLOAT, line_width=5, opacity=0.5,
                            representation='surface')

                    cs_box_ids['caero_control_surfaces'].extend(aero_element_ids)
                    cs_box_ids[cs_name].extend(aero_element_ids)
                else:
                    aelist_ref = aesurf.alid1_ref
                    if aelist_ref is None:
                        self.log.error('AESURF does not reference an AELIST\n%s' % (
                            aesurf.rstrip()))
                        continue
                    ncaeros_cs += len(aelist_ref.elements)

                    cs_name = '%s_control_surface' % aesurf.label
                    if cs_name not in self.alt_grids:
                        self.gui.create_alternate_vtk_grid(
                            cs_name, color=PINK_FLOAT, line_width=5, opacity=0.5,
                            representation='surface')

                    cs_box_ids['caero_control_surfaces'].extend(aelist_ref.elements)
                    cs_box_ids[cs_name].extend(aelist_ref.elements)
                    if aesurf.alid2 is not None:
                        aelist_ref = aesurf.alid2_ref
                        ncaeros_cs += len(aelist_ref.elements)
                        cs_box_ids[cs_name].extend(aelist_ref.elements)
                        cs_box_ids['caero_control_surfaces'].extend(aelist_ref.elements)
        out = (
            has_caero, caero_points, ncaeros, ncaeros_sub, ncaeros_cs,
            ncaeros_points, ncaero_sub_points,
            has_control_surface, box_id_to_caero_element_map, cs_box_ids,
        )
        return out

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
        caero_grid = self.gui.alt_grids['caero']
        for unused_eid, element in sorted(model.caeros.items()):
            if isinstance(element, (CAERO1, CAERO3, CAERO4, CAERO5, CAERO7)):
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
            elif isinstance(element, (CAERO2, BODY7)):
                # slender body

                #if 0:  # pragma: no cover
                # 1D version
                #cpoints = element.get_points()
                #cpoints[:, 2] += zfighting_offset
                #max_cpoints.append(np.array(cpoints).max(axis=0))
                #min_cpoints.append(np.array(cpoints).min(axis=0))

                #elem = vtk.vtkLine()
                #elem.GetPointIds().SetId(0, j)
                #elem.GetPointIds().SetId(1, j + 1)
                #points.InsertPoint(j, *cpoints[0])
                #points.InsertPoint(j + 1, *cpoints[1])
                #j += 2
                #caero_grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())

                #else:
                # 3D version
                xyz, elems = element.get_points_elements_3d()
                assert xyz is not None, element
                xyz[:, 2] += zfighting_offset
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

                    #cpoints = element.get_points()
                    #cpoints[0][2] += zfighting_offset
                    #cpoints[1][2] += zfighting_offset
                    #max_cpoints.append(np.array(cpoints).max(axis=0))
                    #min_cpoints.append(np.array(cpoints).min(axis=0))

                    caero_grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                    j += 4
            else:
                self.gui.log_info("skipping %s" % element.type)

        if ncaeros_points and len(max_cpoints):
            self.gui.log_info('CAERO.max = %s' % np.vstack(max_cpoints).max(axis=0))
            self.gui.log_info('CAERO.min = %s' % np.vstack(min_cpoints).min(axis=0))
        caero_grid.SetPoints(points)
        #self.gui.alt_grids['caero']
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
        for unused_eid, element in sorted(model.caeros.items()):
            if isinstance(element, (CAERO1, CAERO3, CAERO4, CAERO5, CAERO7)):
                pointsi, elementsi = element.panel_points_elements()

                ipoint = 0
                for ipoint, pointii in enumerate(pointsi):
                    points.InsertPoint(j + ipoint, *pointii)

                elem = vtkQuad()
                for elementi in elementsi:
                    elem = vtkQuad()
                    elem.GetPointIds().SetId(0, j + elementi[0])
                    elem.GetPointIds().SetId(1, j + elementi[1])
                    elem.GetPointIds().SetId(2, j + elementi[2])
                    elem.GetPointIds().SetId(3, j + elementi[3])
                    self.gui.alt_grids['caero_subpanels'].InsertNextCell(
                        vtk_type, elem.GetPointIds())
                j += ipoint + 1
            else:
                self.gui.log_info("skipping %s" % element.type)
        self.gui.alt_grids['caero_subpanels'].SetPoints(points)
        return j

    def set_caero_control_surface_grid(self, name, cs_box_ids,
                                       box_id_to_caero_element_map,
                                       caero_points, note=None,
                                       zfighting_offset=0.001, store_msg=False):
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
        label : str / None
            None : no label will be used
            str : the name of the control surface card will be placed
            at the centroid of the panel
        zfighting_offset : float
            z-fighting is when two elements "fight" for who is in front
            leading.  The standard way to fix this is to bump the
            element.

        Returns
        -------
        j : int
            ???
        """
        j = 0
        boxes_to_show, stored_msg = check_for_missing_control_surface_boxes(
            name, cs_box_ids, box_id_to_caero_element_map, self.gui.log,
            store_msg=store_msg)
        #if not boxes_to_show:
            #print('*%s' % name)
            #print('*%s' % boxes_to_show)
            #return

        areas = []
        centroids = []
        vtk_type = vtkQuad().GetCellType()

        all_points = []
        #if name not in self.gui.alt_grids:
            #print('**%s' % name)
            #return

        grid = self.gui.alt_grids[name]
        for box_id in boxes_to_show:
            elementi = box_id_to_caero_element_map[box_id]
            pointsi = caero_points[elementi]
            centroid = (pointsi[0] + pointsi[1] + pointsi[2] + pointsi[3]) / 4.
            area = np.linalg.norm(np.cross(pointsi[2] - pointsi[0], pointsi[3] - pointsi[1])) / 2.
            if area == 0.0:
                print('box_id=%i has 0 area' % box_id)
                continue

            elem = vtkQuad()
            elem.GetPointIds().SetId(0, j)
            elem.GetPointIds().SetId(1, j + 1)
            elem.GetPointIds().SetId(2, j + 2)
            elem.GetPointIds().SetId(3, j + 3)
            grid.InsertNextCell(vtk_type, elem.GetPointIds())
            all_points.append(pointsi)
            centroids.append(centroid)
            areas.append(area)
            j += 4

        if len(all_points) == 0:
            self.log.error('deleting %r' % name)

            # name = spline_1000_boxes
            sname = name.split('_')
            sname[-1] = 'structure_points'

            # points_name = spline_1000_structure_points
            points_name = '_'.join(sname)
            self.log.error('deleting %r' % points_name)
            if name in self.gui.alt_grids:
                del self.gui.alt_grids[name]

            if points_name in self.gui.alt_grids:
                del self.gui.alt_grids[points_name]

            if name in self.gui.geometry_properties:
                del self.gui.geometry_properties[name]

            if points_name in self.gui.geometry_properties:
                del self.gui.geometry_properties[points_name]
            return stored_msg
        # combine all the points
        all_points_array = np.vstack(all_points)

        # shift z to remove z-fighting with caero in surface representation
        all_points_array[:, [1, 2]] += zfighting_offset

        # get the vtk object
        points = numpy_to_vtk_points(all_points_array, deep=0)

        #if missing_boxes:
            #msg = 'Missing CAERO AELIST boxes: ' + str(missing_boxes)
            #self.gui.log_error(msg)
        if note:
            # points_list (15, 4, 3) = (elements, nodes, 3)
            x, y, z = np.average(centroids, weights=areas, axis=0)
            text = str(note)
            #slot = self.gui.label_actors[-1]
            slot = self.gui.geometry_properties[name].label_actors
            self.gui.create_annotation(text, slot, x, y, z)

        self.gui.alt_grids[name].SetPoints(points)
        return stored_msg

    def _set_conm_grid(self, nconm2, model):
        """
        creates the mass secondary actor called:
         - conm2

        which includes:
         - CONM2
         - CMASS1
         - CMASS2

        because it's really a "mass" actor
        """
        j = 0
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(nconm2)

        #sphere_size = self._get_sphere_size(dim_max)
        alt_grid = self.gui.alt_grids['conm2']
        for unused_eid, element in sorted(model.masses.items()):
            if isinstance(element, CONM2):
                xyz_nid = element.nid_ref.get_position()
                centroid = element.offset(xyz_nid)
                #centroid_old = element.Centroid()
                #assert np.all(np.allclose(centroid_old, centroid)), 'centroid_old=%s new=%s' % (centroid_old, centroid)
                #d = norm(xyz - c)
                points.InsertPoint(j, *centroid)

                #if 1:
                elem = vtk.vtkVertex()
                elem.GetPointIds().SetId(0, j)
                #else:
                    #elem = vtk.vtkSphere()
                    #elem.SetRadius(sphere_size)
                    #elem.SetCenter(points.GetPoint(j))

                alt_grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                j += 1
            elif element.type in ['CMASS1', 'CMASS2']:
                centroid = element.Centroid()
                #n1 = element.G1()
                #n2 = element.G2()
                #print('n1=%s n2=%s centroid=%s' % (n1, n2, centroid))
                points.InsertPoint(j, *centroid)

                elem = vtk.vtkVertex()
                elem.GetPointIds().SetId(0, j)
                alt_grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                j += 1
            else:
                self.gui.log_info("skipping %s" % element.type)
        alt_grid.SetPoints(points)

    def set_spc_mpc_suport_grid(self, model, nid_to_pid_map):
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
        for subcase_id, subcase in sorted(model.subcases.items()):
            if 'SPC' in subcase:
                spc_id = subcase.get_parameter('SPC')[0]
                if spc_id is not None and spc_id not in spc_ids_used:
                    spc_ids_used.add(spc_id)
                    nspcs = model.card_count['SPC'] if 'SPC' in model.card_count else 0
                    nspc1s = model.card_count['SPC1'] if 'SPC1' in model.card_count else 0
                    nspcds = model.card_count['SPCD'] if 'SPCD' in model.card_count else 0

                    ## TODO: this line seems too loose...
                    ## TODO: why aren't SPCDs included?
                    if nspcs + nspc1s + nspcds:
                        spc_names += self._fill_spc(spc_id, model, nid_to_pid_map)

            # rigid body elements and MPCs
            if 'MPC' in subcase:
                mpc_id = subcase.get_parameter('MPC')[0]
                if mpc_id is not None and mpc_id not in mpc_ids_used:
                    mpc_ids_used.add(mpc_id)

                    ## TODO: this line seems too loose
                    nmpcs = model.card_count['MPC'] if 'MPC' in model.card_count else 0
                    if nmpcs:
                        lines = model.get_MPCx_node_ids(mpc_id, stop_on_failure=False)
                        lines2 = list(lines) + rigid_lines
                        mpc_names += self._fill_dependent_independent(
                            mpc_id, model, lines2)

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
                        suport_name = self._fill_suport(suport_id, subcase_id, model)
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
            self.gui.create_alternate_vtk_grid(
                grid_name, color=RED_FLOAT, opacity=1.0, point_size=4,
                representation='point', is_visible=True)

        if len(mpc_names) == 0 and len(rigid_lines):
            # handle RBEs without MPCs
            mpc_id = 0
            mpc_names += self._fill_dependent_independent(
                mpc_id, model, rigid_lines)

        geometry_names = spc_names + mpc_names + suport_names
        return geometry_names

    def _fill_spc(self, spc_id, model, nid_to_pid_map):
        """creates the spc secondary actors"""
        spc_name = 'spc_id=%i' % spc_id
        spc_names = [spc_name]
        self.gui.create_alternate_vtk_grid(
            spc_name, color=PURPLE_FLOAT, line_width=5, opacity=1.,
            point_size=5, representation='point', is_visible=False)

        # node_ids = model.get_SPCx_node_ids(spc_id)
        node_ids_c1 = model.get_SPCx_node_ids_c1(
            spc_id, stop_on_failure=False)

        node_ids = []
        for nid, c1 in node_ids_c1.items():
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
        self._add_nastran_nodes_to_grid(spc_name, node_ids, model, msg)
        return spc_names

    def create_bar_pin_flag_text(self, unused_pin_flag=None):
        """
        Lists the pin flag for each element (that has a pin flag)
        self.nid_release_map is set by ``_fill_bar_yz``

        TODO: needs a better interface in the gui
        """
        nids = []
        text = []
        #result_name = self.icase
        result_name = str('ElementID')
        for nid, data in sorted(self.nid_release_map.items()):
            sub_release_map = defaultdict(str)
            for (eid, pin_flagi) in data:
                sub_release_map[pin_flagi] += (str(eid) + ', ')
            texti = '\n'.join(['%s-%s'  % (pin_flagi, msg.rstrip(', '))
                               for (pin_flagi, msg) in sorted(iteritems(sub_release_map))])

            # super messy
            #texti = ', '.join(['%s-%s'  % (pin_flagi, eid) for (eid, pin_flagi) in data])
            nids.append(nid)
            text.append(texti)
        self.gui.mark_nodes(nids, result_name, text)

    def _fill_bar_yz(self, unused_dim_max, model, icase, cases, form, debug=False):
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
        self.gui.create_alternate_vtk_grid(
            'Bar Nodes', color=RED_FLOAT, line_width=1, opacity=1.,
            point_size=5, representation='point', bar_scale=0., is_visible=False)
        msg = ", which is required by 'Bar Nodes'"
        self._add_nastran_nodes_to_grid('Bar Nodes', bar_nids, model, msg)


        geo_form = form[2]
        bar_form = ('CBAR / CBEAM', None, [])
        #print('geo_form =', geo_form)
        #bar_types2 = {}
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

                self.gui.create_alternate_vtk_grid(
                    bar_y, color=GREEN_FLOAT, line_width=5, opacity=1.,
                    point_size=5, representation='bar', bar_scale=scale, is_visible=False)
                self.gui.create_alternate_vtk_grid(
                    bar_z, color=BLUE_FLOAT, line_width=5, opacity=1.,
                    point_size=5, representation='bar', bar_scale=scale, is_visible=False)

                self._add_nastran_lines_xyz_to_grid(bar_y, lines_bar_y, eids)
                self._add_nastran_lines_xyz_to_grid(bar_z, lines_bar_z, eids)

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

    def _add_nastran_lines_xyz_to_grid(self, name, lines, eids):
        """creates the bar orientation vector lines"""
        nlines = len(lines)
        nnodes = nlines * 2
        if nlines == 0:
            return

        assert name != u'Bar Nodes', name
        grid = self.gui.alt_grids[name]

        bar_eids = np.asarray(eids, dtype='int32')
        bar_lines = np.asarray(lines, dtype='float32').reshape(nlines, 6)
        self.bar_eids[name] = bar_eids
        self.bar_lines[name] = bar_lines

        nodes = bar_lines.reshape(nlines * 2, 3)
        points = numpy_to_vtk_points(nodes)
        elements = np.arange(0, nnodes, dtype='int32').reshape(nlines, 2)

        etype = 3 # vtk.vtkLine().GetCellType()
        create_vtk_cells_of_constant_element_type(grid, elements, etype)
        grid.SetPoints(points)

    def _fill_dependent_independent(self, mpc_id, model, lines):
        """creates the mpc actors"""
        if not lines:
            return []

        if mpc_id == 0:
            depname = 'rigid_dependent'
            indname = 'rigid_independent'
            linename = 'rigid_lines'
        else:
            depname = 'mpc_id=%i_dependent' % mpc_id
            indname = 'mpc_id=%i_independent' % mpc_id
            linename = 'mpc_id=%i_lines' % mpc_id
        self.gui.create_alternate_vtk_grid(
            depname, color=GREEN_FLOAT, line_width=5, opacity=1.,
            point_size=5, representation='point', is_visible=False)
        self.gui.create_alternate_vtk_grid(
            indname, color=LIGHT_GREEN_FLOAT, line_width=5, opacity=1.,
            point_size=5, representation='point', is_visible=False)
        self.gui.create_alternate_vtk_grid(
            linename, color=LIGHT_GREEN_FLOAT, line_width=5, opacity=1.,
            point_size=5, representation='wire', is_visible=False)

        lines2 = []
        for line in lines:
            if line not in lines2:
                lines2.append(line)
        lines = np.array(lines2, dtype='int32')
        dependent = (lines[:, 0])
        independent = np.unique(lines[:, 1])
        self.dependents_nodes.update(dependent)
        unused_node_ids = np.unique(lines.ravel())

        msg = ', which is required by %r' % depname
        self._add_nastran_nodes_to_grid(depname, dependent, model, msg)

        msg = ', which is required by %r' % indname
        self._add_nastran_nodes_to_grid(indname, independent, model, msg)

        msg = ', which is required by %r' % linename
        self._add_nastran_lines_to_grid(linename, lines, model)

        mpc_names = [depname, indname, linename]
        return mpc_names

    def _add_nastran_nodes_to_grid(self, name, node_ids, model, msg, store_msg=False):
        """used to create MPC independent/dependent nodes"""
        nnodes = len(node_ids)
        stored_msg = []
        if nnodes == 0:
            msg = '0 nodes added for %r' % name
            out_msg = store_warning(model.log, store_msg, msg)
            return out_msg
        self.gui.follower_nodes[name] = node_ids

        #numpy_to_vtk_points(nodes)
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(nnodes)

        j = 0
        nid_map = self.gui.nid_map
        alt_grid = self.gui.alt_grids[name]
        missing_nodes = []
        for nid in sorted(node_ids):
            try:
                unused_i = nid_map[nid]
            except KeyError:
                missing_nodes.append(str(nid))
                continue

            if nid not in model.nodes:
                # I think this hits for SPOINTs
                missing_nodes.append(str(nid))
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

            alt_grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            j += 1

        out_msg = ''
        if missing_nodes:
            stored_msg = 'nids=[%s] do not exist%s' % (', '.join(missing_nodes), msg)
        alt_grid.SetPoints(points)

        if stored_msg:
            out_msg = store_warning(model.log, store_msg, stored_msg)
        return out_msg

    def _add_nastran_spoints_to_grid(self, model):
        """used to create SPOINTs"""
        if not model.spoints:
            return
        spoint_ids = list(model.spoints.keys())
        assert isinstance(spoint_ids, list), type(spoint_ids)

        nspoints = len(spoint_ids)
        name = 'SPoints'
        if nspoints == 0:
            model.log.warning('0 spoints added for %r' % name)
            return
        self.gui.create_alternate_vtk_grid(
            name, color=BLUE_FLOAT, line_width=1, opacity=1.,
            point_size=5, representation='point', bar_scale=0., is_visible=True)

        self.gui.follower_nodes[name] = spoint_ids
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(nspoints)

        j = 0
        nid_map = self.gui.nid_map
        alt_grid = self.gui.alt_grids[name]
        for spointi in sorted(spoint_ids):
            try:
                unused_i = nid_map[spointi]
            except KeyError:
                model.log.warning('spointi=%s does not exist' % spointi)
                continue

            if spointi not in model.spoints:
                model.log.warning('spointi=%s doesnt exist' % spointi)
                continue
            # point = self.grid.GetPoint(i)
            # points.InsertPoint(j, *point)

            points.InsertPoint(j, 0., 0., 0.)

            elem = vtk.vtkVertex()
            elem.GetPointIds().SetId(0, j)
            alt_grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            j += 1
        alt_grid.SetPoints(points)

    def _add_nastran_lines_to_grid(self, name, lines, model, nid_to_pid_map=None):
        """used to create MPC lines"""
        nlines = lines.shape[0]
        #nids = np.unique(lines)
        #nnodes = len(nids)
        nnodes = nlines * 2
        if nnodes == 0:
            return
        self.gui.follower_nodes[name] = lines.ravel()
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(nnodes)

        j = 0
        nid_map = self.gui.nid_map
        alt_grid = self.gui.alt_grids[name]
        for nid1, nid2 in lines:
            try:
                unused_i1 = nid_map[nid1]
            except KeyError:
                model.log.warning('nid=%s does not exist' % nid1)
                continue
            try:
                unused_i2 = nid_map[nid2]
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
            alt_grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            j += 2
        alt_grid.SetPoints(points)

    def _fill_suport(self, suport_id, unused_subcase_id, model):
        """creates SUPORT and SUPORT1 nodes"""
        suport_name = 'suport1_id=%i' % suport_id
        self.gui.create_alternate_vtk_grid(
            suport_name, color=RED_FLOAT, line_width=5, opacity=1., point_size=4,
            representation='point', is_visible=False)
        suport_nids = get_suport_node_ids(model, suport_id)
        msg = ', which is required by %r' % suport_name
        self._add_nastran_nodes_to_grid(suport_name, suport_nids, model, msg)
        return suport_name

    def _get_sphere_size(self, dim_max):
        return 0.01 * dim_max

    def _map_elements3(self, nid_map, model, unused_j, unused_dim_max,
                       nid_cp_cd, xref_loads=True):
        """
        Much, much faster way to add elements that directly builds the VTK objects
        rather than using for loops.

        Returns
        -------
        nid_to_pid_map  : dict
            node to property id map
            used to show SPC constraints (we don't want to show constraints on 456 DOFs)
        icase : int
            the result number
        cases : dict
            the GuiResult objects
        form : List[???, ???, ???]
            the Results sidebar data

        TDOO: Not quite done on:
               - ???
        """
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

        nelements = len(model.elements)
        xyz_cid0 = self.xyz_cid0
        pids_array = np.zeros(nelements, dtype='int32')
        eids_array = np.zeros(nelements, dtype='int32')
        mcid_array = np.full(nelements, -1, dtype='int32')
        theta_array = np.full(nelements, np.nan, dtype='float32')
        dim_array = np.full(nelements, -1, dtype='int32')
        nnodes_array = np.full(nelements, -1, dtype='int32')

        # quality
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
        normals = np.full((nelements, 3), np.nan, 'float32')


        nids_list = []
        ieid = 0
        cell_offset = 0

        dtype = get_numpy_idtype_for_vtk()

        cell_types_array = np.zeros(nelements, dtype=dtype)
        cell_offsets_array = np.zeros(nelements, dtype=dtype)

        cell_type_point = vtk.vtkVertex().GetCellType()
        cell_type_line = vtk.vtkLine().GetCellType()
        cell_type_tri3 = vtkTriangle().GetCellType()
        cell_type_tri6 = vtkQuadraticTriangle().GetCellType()
        cell_type_quad4 = vtkQuad().GetCellType()
        #cell_type_quad8 = vtkQuadraticQuad().GetCellType()
        cell_type_tetra4 = vtkTetra().GetCellType()
        cell_type_tetra10 = vtkQuadraticTetra().GetCellType()
        cell_type_pyram5 = vtkPyramid().GetCellType()
        #cell_type_pyram13 = vtk.vtkQuadraticPyramid().GetCellType()
        cell_type_penta6 = vtkWedge().GetCellType()
        cell_type_penta15 = vtkQuadraticWedge().GetCellType()
        cell_type_hexa8 = vtkHexahedron().GetCellType()
        cell_type_hexa20 = vtkQuadraticHexahedron().GetCellType()

        # per gui/testing_methods.py/create_vtk_cells_of_constant_element_type
        #1  = vtk.vtkVertex().GetCellType()
        #3  = vtkLine().GetCellType()
        #5  = vtkTriangle().GetCellType()
        #9  = vtk.vtkQuad().GetCellType()
        #10 = vtkTetra().GetCellType()
        #vtkPenta().GetCellType()
        #vtkHexa().GetCellType()
        #vtkPyram().GetCellType()

        skipped_etypes = set([])
        all_nids = nid_cp_cd[:, 0]
        ieid = 0
        for eid, elem in sorted(iteritems(model.elements)):
            if ieid % 5000 == 0 and ieid > 0:
                print('  map_elements = %i' % ieid)
            etype = elem.type
            nnodes = None
            nids = None
            pid = None
            cell_type = None
            inids = None

            dideal_thetai = np.nan
            min_thetai = np.nan
            max_thetai = np.nan
            #max_thetai = np.nan
            max_skew = np.nan
            #max_warp = np.nan
            max_warp = np.nan
            aspect_ratio = np.nan
            areai = np.nan
            area_ratioi = np.nan
            taper_ratioi = np.nan
            min_edge_lengthi = np.nan
            normali = np.nan
            if etype in ['CTRIA3', 'CTRIAR', 'CTRAX3', 'CPLSTN3']:
                nids = elem.nodes
                pid = elem.pid
                cell_type = cell_type_tri3 # 5
                inids = np.searchsorted(all_nids, nids)
                p1, p2, p3 = xyz_cid0[inids, :]
                out = tri_quality(p1, p2, p3)
                (areai, max_skew, aspect_ratio,
                 min_thetai, max_thetai, dideal_thetai, min_edge_lengthi) = out
                normali = np.cross(p1 - p2, p1 - p3)
                if isinstance(elem.theta_mcid, float):
                    theta_array[ieid] = elem.theta_mcid
                else:
                    mcid_array[ieid] = elem.theta_mcid
                nnodes = 3
                dim = 2

            elif etype in ['CQUAD4', 'CQUADR', 'CPLSTN4', 'CQUADX4']:
                nids = elem.nodes
                pid = elem.pid
                cell_type = cell_type_quad4 #9
                inids = np.searchsorted(all_nids, nids)
                p1, p2, p3, p4 = xyz_cid0[inids, :]
                out = quad_quality(p1, p2, p3, p4)
                (areai, taper_ratioi, area_ratioi, max_skew, aspect_ratio,
                 min_thetai, max_thetai, dideal_thetai, min_edge_lengthi) = out
                normali = np.cross(p1 - p3, p2 - p4)
                if isinstance(elem.theta_mcid, float):
                    theta_array[ieid] = elem.theta_mcid
                else:
                    mcid_array[ieid] = elem.theta_mcid
                nnodes = 4
                dim = 2

            elif etype == 'CTRIA6':
                nids = elem.nodes
                pid = elem.pid
                if None in nids:
                    cell_type = cell_type_tri3
                    inids = np.searchsorted(all_nids, nids[:3])
                    nids = nids[:3]
                    p1, p2, p3 = xyz_cid0[inids, :]
                    nnodes = 3
                else:
                    cell_type = cell_type_tri6
                    inids = np.searchsorted(all_nids, nids)
                    p1, p2, p3, p4, unused_p5, unused_p6 = xyz_cid0[inids, :]
                    nnodes = 6
                out = tri_quality(p1, p2, p3)
                (areai, max_skew, aspect_ratio,
                 min_thetai, max_thetai, dideal_thetai, min_edge_lengthi) = out
                normali = np.cross(p1 - p2, p1 - p3)
                if isinstance(elem.theta_mcid, float):
                    theta_array[ieid] = elem.theta_mcid
                else:
                    mcid_array[ieid] = elem.theta_mcid
                dim = 2
            elif etype == 'CQUAD8':
                nids = elem.nodes
                pid = elem.pid
                if None in nids:
                    cell_type = cell_type_tri3
                    inids = np.searchsorted(all_nids, nids[:4])
                    nids = nids[:4]
                    p1, p2, p3, p4 = xyz_cid0[inids, :]
                    nnodes = 4
                else:
                    cell_type = cell_type_tri6
                    inids = np.searchsorted(all_nids, nids)
                    p1, p2, p3, p4 = xyz_cid0[inids[:4], :]
                    nnodes = 8
                out = quad_quality(p1, p2, p3, p4)
                (areai, taper_ratioi, area_ratioi, max_skew, aspect_ratio,
                 min_thetai, max_thetai, dideal_thetai, min_edge_lengthi) = out
                normali = np.cross(p1 - p3, p2 - p4)
                if isinstance(elem.theta_mcid, float):
                    theta_array[ieid] = elem.theta_mcid
                else:
                    mcid_array[ieid] = elem.theta_mcid
                nnodes = 4
                dim = 2

            elif etype == 'CSHEAR':
                nids = elem.nodes
                pid = elem.pid
                cell_type = cell_type_quad4 #9
                inids = np.searchsorted(all_nids, nids)
                p1, p2, p3, p4 = xyz_cid0[inids, :]
                out = quad_quality(p1, p2, p3, p4)
                (areai, taper_ratioi, area_ratioi, max_skew, aspect_ratio,
                 min_thetai, max_thetai, dideal_thetai, min_edge_lengthi) = out
                normali = np.cross(p1 - p3, p2 - p4)
                nnodes = 4
                dim = 2

            elif etype == 'CTETRA':
                nids = elem.nodes
                pid = elem.pid
                if None in nids:
                    cell_type = cell_type_tetra4
                    nids = nids[:4]
                    nnodes = 4
                else:
                    cell_type = cell_type_tetra10
                    nnodes = 10
                inids = np.searchsorted(all_nids, nids)
                min_thetai, max_thetai, dideal_thetai, min_edge_lengthi = get_min_max_theta(
                    _ctetra_faces, nids, nid_map, xyz_cid0)
                dim = 3

            elif etype == 'CHEXA':
                nids = elem.nodes
                pid = elem.pid
                if None in nids:
                    cell_type = cell_type_hexa8
                    nids = nids[:8]
                    nnodes = 8
                else:
                    cell_type = cell_type_hexa20
                    nnodes = 20
                inids = np.searchsorted(all_nids, nids)
                min_thetai, max_thetai, dideal_thetai, min_edge_lengthi = get_min_max_theta(
                    _chexa_faces, nids, nid_map, xyz_cid0)
                dim = 3

            elif etype == 'CPENTA':
                nids = elem.nodes
                pid = elem.pid

                if None in nids:
                    cell_type = cell_type_penta6
                    nids = nids[:6]
                    nnodes = 6
                else:
                    cell_type = cell_type_penta15
                    nnodes = 15

                inids = np.searchsorted(all_nids, nids)
                min_thetai, max_thetai, dideal_thetai, min_edge_lengthi = get_min_max_theta(
                    _cpenta_faces, nids, nid_map, xyz_cid0)
                dim = 3
            elif etype == 'CPYRAM':
                # TODO: assuming 5
                nids = elem.nodes
                pid = elem.pid
                if None in nids:
                    cell_type = cell_type_pyram5
                    nids = nids[:5]
                    nnodes = 5
                else:
                    cell_type = cell_type_penta15
                    nnodes = 15
                inids = np.searchsorted(all_nids, nids)
                min_thetai, max_thetai, dideal_thetai, min_edge_lengthi = get_min_max_theta(
                    _cpyram_faces, nids, nid_map, xyz_cid0)
                dim = 3
            elif etype in ['CELAS2', 'CELAS4', 'CDAMP4']:
                # these can have empty nodes and have no property
                # CELAS1: 1/2 GRID/SPOINT and pid
                # CELAS2: 1/2 GRID/SPOINT, k, ge, and s
                # CELAS3: 1/2 SPOINT and pid
                # CELAS4: 1/2 SPOINT and k
                nids = elem.nodes
                assert nids[0] != nids[1]
                if None in nids:
                    assert nids[0] is not None, nids
                    assert nids[1] is None, nids
                    nids = [nids[0]]
                    cell_type = cell_type_point
                    nnodes = 1
                else:
                    nids = elem.nodes
                    assert nids[0] != nids[1]
                    cell_type = cell_type_line
                    nnodes = 2
                inids = np.searchsorted(all_nids, nids)
                pid = 0
                dim = 0
            elif etype in ['CBUSH', 'CBUSH1D', 'CBUSH2D',
                           'CELAS1', 'CELAS3',
                           'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP5',
                           'CFAST', 'CGAP', 'CVISC']:
                nids = elem.nodes
                assert nids[0] != nids[1]
                assert None not in nids, 'nids=%s\n%s' % (nids, elem)
                pid = elem.pid
                cell_type = cell_type_line
                inids = np.searchsorted(all_nids, nids)
                nnodes = 2
                dim = 0
            elif etype in ['CBAR', 'CBEAM']:
                nids = elem.nodes
                pid = elem.pid
                pid_ref = model.Property(pid)
                areai = pid_ref.Area()
                cell_type = cell_type_line
                inids = np.searchsorted(all_nids, nids)
                p1, p2 = xyz_cid0[inids, :]
                min_edge_lengthi = norm(p2 - p1)
                nnodes = 2
                dim = 1
            elif etype in ['CROD', 'CTUBE']:
                nids = elem.nodes
                pid = elem.pid
                pid_ref = model.Property(pid)
                areai = pid_ref.Area()
                cell_type = cell_type_line
                inids = np.searchsorted(all_nids, nids)
                p1, p2 = xyz_cid0[inids, :]
                min_edge_lengthi = norm(p2 - p1)
                nnodes = 2
                dim = 1
            elif etype == 'CONROD':
                nids = elem.nodes
                areai = elem.Area()
                pid = 0
                cell_type = cell_type_line
                inids = np.searchsorted(all_nids, nids)
                p1, p2 = xyz_cid0[inids, :]
                min_edge_lengthi = norm(p2 - p1)
                nnodes = 2
                dim = 1
            #------------------------------
            # rare
            #elif etype == 'CIHEX1':
                #nids = elem.nodes
                #pid = elem.pid
                #cell_type = cell_type_hexa8
                #inids = np.searchsorted(all_nids, nids)
                #min_thetai, max_thetai, dideal_thetai, min_edge_lengthi = get_min_max_theta(
                    #_chexa_faces, nids, nid_map, xyz_cid0)
                #nnodes = 8
                #dim = 3
            elif etype == 'CHBDYE':
                #self.eid_map[eid] = ieid
                eid_solid = elem.eid2
                side = elem.side
                element_solid = model.elements[eid_solid]

                mapped_inids = SIDE_MAP[element_solid.type][side]
                side_inids = [nid - 1 for nid in mapped_inids]
                nodes = element_solid.node_ids

                pid = 0
                nnodes = len(side_inids)
                nids = [nodes[inid] for inid in side_inids]
                inids = np.searchsorted(all_nids, nids)

                if len(side_inids) == 4:
                    cell_type = cell_type_quad4
                else:
                    msg = 'element_solid:\n%s' % (str(element_solid))
                    msg += 'mapped_inids = %s\n' % mapped_inids
                    msg += 'side_inids = %s\n' % side_inids
                    msg += 'nodes = %s\n' % nodes
                    #msg += 'side_nodes = %s\n' % side_nodes
                    raise NotImplementedError(msg)
            else:
                #raise NotImplementedError(elem)
                skipped_etypes.add(etype)
                nelements -= 1
                continue
            #for nid in nids:
                #assert isinstance(nid, integer_types), 'not an integer. nids=%s\n%s' % (nids, elem)
                #assert nid != 0, 'not a positive integer. nids=%s\n%s' % (nids, elem)

            assert inids is not None
            if not np.array_equal(all_nids[inids], nids):
                msg = 'all_nids[inids]=%s nids=%s\n%s' % (all_nids[inids], nids, elem)
                raise RuntimeError(msg)

            assert cell_type is not None
            assert cell_offset is not None
            assert eid is not None
            assert pid is not None
            assert dim is not None
            assert nnodes is not None
            nids_list.append(nnodes)
            nids_list.extend(inids)
            normals[ieid] = normali
            eids_array[ieid] = eid
            pids_array[ieid] = pid
            dim_array[ieid] = dim
            cell_types_array[ieid] = cell_type
            cell_offsets_array[ieid] = cell_offset  # I assume the problem is here
            cell_offset += nnodes + 1
            self.eid_map[eid] = ieid

            min_interior_angle[ieid] = min_thetai
            max_interior_angle[ieid] = max_thetai
            dideal_theta[ieid] = dideal_thetai
            max_skew_angle[ieid] = max_skew
            max_warp_angle[ieid] = max_warp
            max_aspect_ratio[ieid] = aspect_ratio
            area[ieid] = areai
            area_ratio[ieid] = area_ratioi
            taper_ratio[ieid] = taper_ratioi
            min_edge_length[ieid] = min_edge_lengthi
            ieid += 1

        #print('self.eid_map =', self.eid_map)

        icells_zero = np.where(cell_types_array == 0)[0]
        # TODO: I'd like to get rid of deep=1, but it'll crash the edges
        deep = 1
        if len(icells_zero):
            icells = np.where(cell_types_array != 0)[0]
            if len(icells) == 0:
                self.log.error('skipped_etypes = %s' % skipped_etypes)
                raise RuntimeError('there are no elements...')
            eids_array = eids_array[icells]
            pids_array = pids_array[icells]
            #dim_array = pids_array[dim_array]
            cell_types_array = cell_types_array[icells]
            cell_offsets_array = cell_offsets_array[icells]
            nnodes_array = nnodes_array[icells]
            normals = normals[icells, :]
            #deep = 1
        #print('deep = %s' % deep)
        if skipped_etypes:
            self.log.error('skipped_etypes = %s' % list(skipped_etypes))
            #print('skipped_etypes = %s' % skipped_etypes)
        if len(pids_array) != nelements:
            msg = 'nelements=%s len(pids_array)=%s' % (nelements, len(pids_array))
            raise RuntimeError(msg)
        if len(cell_offsets_array) != nelements:
            msg = 'nelements=%s len(cell_offsets_array)=%s' % (nelements, len(cell_offsets_array))
            raise RuntimeError(msg)

        nids_array = np.array(nids_list, dtype=dtype)

        #-----------------------------------------------------------------
        # saving some data members
        self.element_ids = eids_array

        #print('cell_types_array* = ', cell_types_array.tolist())
        #print('cell_offsets_array* = ', cell_offsets_array.tolist())

        #-----------------------------------------------------------------
        # build the grid

        #self.log.info('nids_array = %s' % nids_array)
        #self.log.info('cell_offsets_array = %s' % cell_offsets_array)
        #self.log.info('cell_types_array = %s' % cell_types_array)

        # Create the array of cells
        cells_id_type = numpy_to_vtkIdTypeArray(nids_array, deep=1)
        vtk_cells = vtk.vtkCellArray()
        vtk_cells.SetCells(nelements, cells_id_type)

        # Cell types
        vtk_cell_types = numpy_to_vtk(
            cell_types_array, deep=deep,
            array_type=vtk.vtkUnsignedCharArray().GetDataType())

        vtk_cell_offsets = numpy_to_vtk(cell_offsets_array, deep=deep,
                                        array_type=vtk.VTK_ID_TYPE)

        grid = self.grid
        #grid = vtk.vtkUnstructuredGrid()
        grid.SetCells(vtk_cell_types, vtk_cell_offsets, vtk_cells)

        #-----------------------------------------------------------------
        # fill the results
        nid_to_pid_map = None
        self.isubcase_name_map = {1: ['Nastran', '']}
        icase = 0
        cases = OrderedDict()
        form = ['Geometry', None, []]
        form0 = form[2]

        subcase_id = 0

        #nids_set = True
        #if nids_set:
        # this intentionally makes a deepcopy
        #nids = np.array(nid_cp_cd[:, 0])

        # this intentionally makes a deepcopy
        cds = np.array(nid_cp_cd[:, 2])
        colormap = self.gui.settings.colormap
        nid_res = GuiResult(subcase_id, 'NodeID', 'NodeID', 'node', all_nids,
                            mask_value=0,
                            nlabels=None,
                            labelsize=None,
                            ncolors=None,
                            colormap=colormap,
                            data_format=None,
                            uname='GuiResult')
        cases[icase] = (nid_res, (0, 'Node ID'))
        form0.append(('Node ID', icase, []))
        icase += 1

        if cds.max() > 0:
            cd_res = GuiResult(0, header='NodeCd', title='NodeCd',
                               location='node', scalar=cds)
            cases[icase] = (cd_res, (0, 'NodeCd'))
            form0.append(('NodeCd', icase, []))
            icase += 1

        eid_res = GuiResult(subcase_id, 'ElementID', 'ElementID', 'centroid', eids_array,
                            mask_value=0,
                            nlabels=None,
                            labelsize=None,
                            ncolors=None,
                            colormap=colormap,
                            data_format=None,
                            uname='GuiResult')
        cases[icase] = (eid_res, (0, 'ElementID'))
        form0.append(('ElementID', icase, []))
        icase += 1

        is_element_dim = True
        #if len(np.unique(dim_array)) > 1:
            #dim_res = GuiResult(subcase_id, 'ElementDim', 'ElementDim', 'centroid', dim_array,
                                   #mask_value=-1,
                                   #nlabels=None,
                                   #labelsize=None,
                                   #ncolors=None,
                                   #colormap=colormap,
                                   #data_format=None,
                                   #uname='GuiResult')
            #cases[icase] = (dim_res, (0, 'ElementDim'))
            #form0.append(('ElementDim', icase, []))
            #icase += 1

        if nnodes_array.max() > -1:
            nnodes_res = GuiResult(subcase_id, 'NNodes/Elem', 'NNodes/Elem',
                                   'centroid', nnodes_array,
                                   mask_value=0,
                                   nlabels=None,
                                   labelsize=None,
                                   ncolors=None,
                                   colormap=colormap,
                                   data_format=None,
                                   uname='GuiResult')
            cases[icase] = (nnodes_res, (0, 'NNodes/Elem'))
            form0.append(('NNodes/Elem', icase, []))
            icase += 1

        #pid_res = GuiResult(subcase_id, 'PropertyID', 'PropertyID', 'centroid', pids_array,
                            #mask_value=0,
                            #nlabels=None,
                            #labelsize=None,
                            #ncolors=None,
                            #colormap=colormap,
                            #data_format=None,
                            #uname='GuiResult')
        #cases[icase] = (pid_res, (0, 'PropertyID'))
        #form0.append(('PropertyID', icase, []))
        #icase += 1

        if len(model.properties) and nelements:
            icase, upids, pcomp, pshell, is_pshell_pcomp = self._build_properties(
                model, nelements, eids_array, pids_array, cases, form0, icase)
            icase = self._build_materials(model, pcomp, pshell, is_pshell_pcomp,
                                          cases, form0, icase)
            try:
                icase = self._build_optimization(model, pids_array, upids,
                                                 nelements, cases, form0, icase)
            except:
                #raise
                s = StringIO()
                traceback.print_exc(file=s)
                sout = s.getvalue()
                self.gui.log_error(sout)
                print(sout)


        #if mcid_array.max() > -1:
            #mcid_res = GuiResult(subcase_id, 'Material Coordinate System', 'MaterialCoord',
                                 #'centroid', mcid_array,
                                 #mask_value=-1,
                                 #nlabels=None,
                                 #labelsize=None,
                                 #ncolors=None,
                                 #colormap=colormap,
                                 #data_format=None,
                                 #uname='GuiResult')
            #cases[icase] = (mcid_res, (0, 'Material Coordinate System'))
            #form0.append(('Material Coordinate System', icase, []))
            #icase += 1

        #if np.isfinite(theta_array).any():
            #print('np.nanmax(theta_array) =', np.nanmax(theta_array))
            #theta_res = GuiResult(subcase_id, 'Theta', 'Theta', 'centroid', theta_array,
                                  #mask_value=None,
                                  #nlabels=None,
                                  #labelsize=None,
                                  #ncolors=None,
                                  #colormap=colormap,
                                  #data_format=None,
                                  #uname='GuiResult')
            #cases[icase] = (theta_res, (0, 'Theta'))
            #form0.append(('Theta', icase, []))
            #icase += 1

        normal_mag = np.linalg.norm(normals, axis=1)
        assert len(normal_mag) == nelements
        normals /= normal_mag.reshape(nelements, 1)
        i_not_nan = np.isnan(normal_mag)

        #if self.make_offset_normals_dim and nelements:
            #material_coord = None
            #icase, normals = self._build_normals_quality(
                #model, nelements, cases, form0, icase,
                #xyz_cid0, material_coord,
                #min_interior_angle, max_interior_angle, dideal_theta,
                #area, max_skew_angle, taper_ratio,
                #max_warp_angle, area_ratio, min_edge_length, max_aspect_ratio)
            #self.normals = normals

        #----------------------------------------------------------

        is_shell = False
        if False in i_not_nan:
            #max_normal = np.nanmax(normal_mag[i_not_nan])
            #is_shell = np.abs(max_normal) > 0.
            is_shell = True
        is_solid = isfinite_and_nonzero(max_interior_angle)
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
                                    colormap=colormap, data_format='%.1f',
                                    uname='NormalResult')


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

            if self.make_nnodes_result and 0:  # pragma: no cover
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

            if np.any(np.isfinite(area_ratio)) and np.nanmax(area_ratio) > 1.:
                arearatio_res = GuiResult(
                    0, header='Area Ratio', title='Area Ratio',
                    location='centroid', scalar=area_ratio)
                cases[icase] = (arearatio_res, (0, 'Area Ratio'))
                form_checks.append(('Area Ratio', icase, []))
                icase += 1

            if np.any(np.isfinite(taper_ratio)) and np.nanmax(taper_ratio) > 1.:
                taperratio_res = GuiResult(
                    0, header='Taper Ratio', title='Taper Ratio',
                    location='centroid', scalar=taper_ratio)
                cases[icase] = (taperratio_res, (0, 'Taper Ratio'))
                form_checks.append(('Taper Ratio', icase, []))
                icase += 1

            if isfinite_and_nonzero(max_warp_angle):
                warp_res = GuiResult(
                    0, header='Max Warp Angle', title='MaxWarpAngle',
                    location='centroid', scalar=np.degrees(max_warp_angle))
                cases[icase + 4] = (warp_res, (0, 'Max Warp Angle'))
                form_checks.append(('Max Warp Angle', icase, []))
                icase += 1

            #if (np.abs(xoffset).max() > 0.0 or np.abs(yoffset).max() > 0.0 or
                #np.abs(zoffset).max() > 0.0):
            # offsets
            #offset_res = GuiResult(
                #0, header='Offset', title='Offset',
                #location='centroid', scalar=offset, data_format='%g')
            #offset_x_res = GuiResult(
                #0, header='OffsetX', title='OffsetX',
                #location='centroid', scalar=xoffset, data_format='%g')
            #offset_y_res = GuiResult(
                #0, header='OffsetY', title='OffsetY',
                #location='centroid', scalar=yoffset, data_format='%g')
            #offset_z_res = GuiResult(
                #0, header='OffsetZ', title='OffsetZ',
                #location='centroid', scalar=zoffset, data_format='%g')

            #cases[icase] = (offset_res, (0, 'Offset'))
            #cases[icase + 1] = (offset_x_res, (0, 'OffsetX'))
            #cases[icase + 2] = (offset_y_res, (0, 'OffsetY'))
            #cases[icase + 3] = (offset_z_res, (0, 'OffsetZ'))

            #form_checks.append(('Offset', icase, []))
            #form_checks.append(('OffsetX', icase + 1, []))
            #form_checks.append(('OffsetY', icase + 2, []))
            #form_checks.append(('OffsetZ', icase + 3, []))
            #icase += 4

            if self.make_xyz or IS_TESTING:
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
            skew = 90. - np.degrees(max_skew_angle)
            #skew_res = GuiResult(0, header='Max Skew Angle', title='MaxSkewAngle',
                                    #location='centroid', scalar=skew)
            if is_element_dim:
                form_checks.append(('ElementDim', icase, []))
            form_checks.append(('Min Edge Length', icase + 1, []))
            form_checks.append(('Min Interior Angle', icase + 2, []))
            form_checks.append(('Max Interior Angle', icase + 3, []))
            form_checks.append(('Max Skew Angle', icase + 4, []))
            cases[icase + 1] = (min_edge_length_res, (0, 'Min Edge Length'))
            cases[icase + 2] = (min_theta_res, (0, 'Min Interior Angle'))
            cases[icase + 3] = (max_theta_res, (0, 'Max Interior Angle'))
            #cases[icase + 4] = (skew_res, (0, 'Max Skew Angle'))
            icase += 4

        else:
            form0.append(('ElementDim', icase, []))
            icase += 1

        if mcid_array is not None and np.abs(mcid_array).max() > 0:
            material_coord_res = GuiResult(
                0, header='MaterialCoord', title='MaterialCoord',
                location='centroid',
                scalar=mcid_array, data_format='%i')
            cases[icase] = (material_coord_res, (0, 'MaterialCoord'))
            form0.append(('MaterialCoord', icase, []))
            icase += 1
        #print(normals)
        #----------------------------------------------------------
        # finishing up vtk
        self.set_glyph_scale_factor(np.nanmean(min_edge_length) * 2.5)  # was 1.5

        grid.Modified()
        if hasattr(grid, 'Update'):  # pragma: no cover
            grid.Update()
        #----------------------------------------------------------
        # finishing up parameters
        self.node_ids = all_nids
        self.normals = normals

        return nid_to_pid_map, icase, cases, form

    def map_elements(self, xyz_cid0, nid_cp_cd, nid_map, model, j, dim_max,
                     plot=True, xref_loads=True):
        """
        Creates the elements

        nid_cp_cd : (nnodes, 3) int ndarray
            the node_id and coordinate systems corresponding to xyz_cid0
            used for setting the NodeID and CD coordinate results
        xyz_cid0 : (nnodes, 3) float ndarray
            the global xyz locations
        nid_map : dict[nid] : nid_index
            nid : int
                the GRID/SPOINT/EPOINT id
            nid_index : int
                the index for the GRID/SPOINT/EPOINT in xyz_cid0
        model : BDF()
            the model object
        j : int
            ???
        dim_max : float
            the max(dx, dy, dz) dimension
            use for ???
        """
        grid = self.gui.grid

        if IS_TESTING:
            self._map_elements3(nid_map, model, j, dim_max,
                                nid_cp_cd, xref_loads=xref_loads)
        #return self._map_elements3(nid_map, model, j, dim_max,
                                   #nid_cp_cd, plot=plot, xref_loads=xref_loads)

        out = self._map_elements1(model, xyz_cid0, nid_cp_cd, dim_max, nid_map, j)
        (nid_to_pid_map, xyz_cid0, pids, nelements, material_coord,
         area, min_interior_angle, max_interior_angle, max_aspect_ratio,
         max_skew_angle, taper_ratio, dideal_theta,
         area_ratio, min_edge_length, max_warp_angle) = out

        #self.grid_mapper.SetResolveCoincidentTopologyToPolygonOffset()
        grid.Modified()
        if hasattr(grid, 'Update'):  # pragma: no cover
            grid.Update()
        #self.gui.log_info("updated grid")

        cases = OrderedDict()


        self.gui.isubcase_name_map = {1: ['Nastran', '']}
        icase = 0
        form = ['Geometry', None, []]
        form0 = form[2]

        #new_cases = True
        # set to True to enable node_ids as an result
        nids_set = True
        if nids_set and self.gui.nnodes > 0:
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
            eid_map = self.gui.eid_map
            for (eid, eid2) in iteritems(eid_map):
                eids[eid2] = eid

            eid_res = GuiResult(0, header='ElementID', title='ElementID',
                                location='centroid', scalar=eids, mask_value=0)
            cases[icase] = (eid_res, (0, 'ElementID'))
            form0.append(('ElementID', icase, []))
            icase += 1
            self.element_ids = eids

        # subcase_id, resultType, vector_size, location, dataFormat
        if len(model.properties) and nelements:
            icase, upids, pcomp, pshell, is_pshell_pcomp = self._build_properties(
                model, nelements, eids, pids, cases, form0, icase)
            icase = self._build_materials(model, pcomp, pshell, is_pshell_pcomp,
                                          cases, form0, icase)

            try:
                icase = self._build_optimization(model, pids, upids, nelements, cases, form0, icase)
            except:
                if IS_TESTING or self.is_testing_flag:
                    raise
                s = StringIO()
                traceback.print_exc(file=s)
                sout = s.getvalue()
                self.gui.log_error(sout)
                print(sout)
                #traceback.print_exc(file=sys.stdout)
                #etype, value, tb = sys.exc_info
                #print(etype, value, tb)
                #raise RuntimeError('Optimization Parsing Error') from e
                #traceback.print_tb(e)
                #print(e)

        #print('nelements=%s eid_map=%s' % (nelements, self.eid_map))
        if nelements and isfinite(min_edge_length):
            glyph_scale = np.nanmean(min_edge_length) * 2.5
            self.gui.set_glyph_scale_factor(glyph_scale)  # was 1.5

        if self.make_offset_normals_dim and nelements:
            icase, normals = self._build_normals_quality(
                model, nelements, cases, form0, icase,
                xyz_cid0, material_coord,
                min_interior_angle, max_interior_angle, dideal_theta,
                area, max_skew_angle, taper_ratio,
                max_warp_angle, area_ratio, min_edge_length, max_aspect_ratio)
            self.normals = normals
        return nid_to_pid_map, icase, cases, form

    def _build_plotels(self, model):
        """creates the plotel actor"""
        nplotels = len(model.plotels)
        if nplotels:
            lines = []
            for unused_eid, element in sorted(model.plotels.items()):
                node_ids = element.node_ids
                lines.append(node_ids)
            lines = np.array(lines, dtype='int32')

            self.gui.create_alternate_vtk_grid(
                'plotel', color=RED_FLOAT, line_width=2, opacity=0.8,
                point_size=5, representation='wire', is_visible=True)
            self._add_nastran_lines_to_grid('plotel', lines, model)

    def _map_elements1(self, model, xyz_cid0, nid_cp_cd, unused_dim_max, nid_map, j):
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
        assert nid_map is not None
        if xyz_cid0 is None:
            nid_to_pid_map = None
            pids = None
            nelements = None
            material_coord = None
            area = None
            min_interior_angle = None
            max_interior_angle = None
            max_aspect_ratio = None
            max_skew_angle = None
            taper_ratio = None
            dideal_theta = None
            area_ratio = None
            min_edge_length = None
            max_warp_angle = None
            out = (
                nid_to_pid_map, xyz_cid0, pids, nelements, material_coord,
                area, min_interior_angle, max_interior_angle, max_aspect_ratio,
                max_skew_angle, taper_ratio, dideal_theta,
                area_ratio, min_edge_length, max_warp_angle,
            )
            return out

        xyz_cid0 = self.xyz_cid0
        nids = nid_cp_cd[:, 0]
        #sphere_size = self._get_sphere_size(dim_max)

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

        grid = self.gui.grid
        self._build_plotels(model)

        #print("map_elements...")
        eid_to_nid_map = self.eid_to_nid_map
        eid_map = self.gui.eid_map
        for (eid, element) in sorted(iteritems(model.elements)):
            eid_map[eid] = i
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

            pid = np.nan
            dideal_thetai = np.nan
            min_thetai = np.nan
            max_thetai = np.nan
            #max_thetai = np.nan
            max_skew = np.nan
            #max_warp = np.nan
            max_warp = np.nan
            aspect_ratio = np.nan
            areai = np.nan
            area_ratioi = np.nan
            taper_ratioi = np.nan
            min_edge_lengthi = np.nan

            if isinstance(element, (CTRIA3, CTRIAR, CTRAX3, CPLSTN3)):
                if isinstance(element, (CTRIA3, CTRIAR)):
                    material_coord[i] = _get_shell_material_coord_int(element)
                elem = vtkTriangle()
                node_ids = element.node_ids
                pid = element.Pid()
                eid_to_nid_map[eid] = node_ids
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
                    material_coord[i] = _get_shell_material_coord_int(element)
                node_ids = element.node_ids
                pid = element.Pid()
                eid_to_nid_map[eid] = node_ids[:3]
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
                grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
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
                #material_coord[i] = element.theta # TODO: no mcid
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
                eid_to_nid_map[eid] = [node_ids[0], node_ids[2], node_ids[4]]
                grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())

            elif isinstance(element, (CQUAD4, CSHEAR, CQUADR, CPLSTN4, CQUADX4)):
                if isinstance(element, (CQUAD4, CQUADR)):
                    material_coord[i] = _get_shell_material_coord_int(element)
                #print('eid=%s theta=%s' % (eid, material_coord[i]))
                node_ids = element.node_ids
                pid = element.Pid()
                for nid in node_ids:
                    if nid is not None:
                        nid_to_pid_map[nid].append(pid)

                eid_to_nid_map[eid] = node_ids

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
                grid.InsertNextCell(9, elem.GetPointIds())

            elif isinstance(element, (CQUAD8, CPLSTN8, CQUADX8)):
                if isinstance(element, CQUAD8):
                    material_coord[i] = _get_shell_material_coord_int(element)
                node_ids = element.node_ids
                pid = element.Pid()
                for nid in node_ids:
                    if nid is not None:
                        nid_to_pid_map[nid].append(pid)
                self.eid_to_nid_map[eid] = node_ids[:4]

                n1, n2, n3, n4 = [nid_map[nid] for nid in node_ids[:4]]
                print('node_ids = ', node_ids)
                print(xyz_cid0)
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
                grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element, (CQUAD, CQUADX)):
                # CQUAD, CQUADX are 9 noded quads
                material_coord[i] = _get_shell_material_coord_int(element)

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
                    elem = vtk.vtkBiQuadraticQuad()
                    elem.GetPointIds().SetId(4, nid_map[node_ids[4]])
                    elem.GetPointIds().SetId(5, nid_map[node_ids[5]])
                    elem.GetPointIds().SetId(6, nid_map[node_ids[6]])
                    elem.GetPointIds().SetId(7, nid_map[node_ids[7]])
                    elem.GetPointIds().SetId(8, nid_map[node_ids[8]])
                else:
                    elem = vtkQuad()
                elem.GetPointIds().SetId(0, n1)
                elem.GetPointIds().SetId(1, n2)
                elem.GetPointIds().SetId(2, n3)
                elem.GetPointIds().SetId(3, n4)

                grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element, CTETRA4):
                elem = vtkTetra()
                node_ids = element.node_ids
                pid = element.Pid()
                for nid in node_ids:
                    nid_to_pid_map[nid].append(pid)
                eid_to_nid_map[eid] = node_ids[:4]
                elem.GetPointIds().SetId(0, nid_map[node_ids[0]])
                elem.GetPointIds().SetId(1, nid_map[node_ids[1]])
                elem.GetPointIds().SetId(2, nid_map[node_ids[2]])
                elem.GetPointIds().SetId(3, nid_map[node_ids[3]])
                grid.InsertNextCell(10, elem.GetPointIds())
                #elem_nid_map = {nid:nid_map[nid] for nid in node_ids[:4]}
                min_thetai, max_thetai, dideal_thetai, min_edge_lengthi = get_min_max_theta(
                    _ctetra_faces, node_ids[:4], nid_map, xyz_cid0)
            elif isinstance(element, CTETRA10):
                node_ids = element.node_ids
                pid = element.Pid()
                for nid in node_ids:
                    if nid is not None:
                        nid_to_pid_map[nid].append(pid)
                eid_to_nid_map[eid] = node_ids[:4]
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
                grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                min_thetai, max_thetai, dideal_thetai, min_edge_lengthi = get_min_max_theta(
                    _ctetra_faces, node_ids[:4], nid_map, xyz_cid0)
            elif isinstance(element, CPENTA6):
                elem = vtkWedge()
                node_ids = element.node_ids
                pid = element.Pid()
                for nid in node_ids:
                    nid_to_pid_map[nid].append(pid)
                eid_to_nid_map[eid] = node_ids[:6]
                elem.GetPointIds().SetId(0, nid_map[node_ids[0]])
                elem.GetPointIds().SetId(1, nid_map[node_ids[1]])
                elem.GetPointIds().SetId(2, nid_map[node_ids[2]])
                elem.GetPointIds().SetId(3, nid_map[node_ids[3]])
                elem.GetPointIds().SetId(4, nid_map[node_ids[4]])
                elem.GetPointIds().SetId(5, nid_map[node_ids[5]])
                grid.InsertNextCell(13, elem.GetPointIds())
                min_thetai, max_thetai, dideal_thetai, min_edge_lengthi = get_min_max_theta(
                    _cpenta_faces, node_ids[:6], nid_map, xyz_cid0)

            elif isinstance(element, CPENTA15):
                node_ids = element.node_ids
                pid = element.Pid()
                for nid in node_ids:
                    if nid is not None:
                        nid_to_pid_map[nid].append(pid)
                eid_to_nid_map[eid] = node_ids[:6]
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
                grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                min_thetai, max_thetai, dideal_thetai, min_edge_lengthi = get_min_max_theta(
                    _cpenta_faces, node_ids[:6], nid_map, xyz_cid0)
            elif isinstance(element, (CHEXA8, CIHEX1)):
                node_ids = element.node_ids
                pid = element.Pid()
                for nid in node_ids:
                    nid_to_pid_map[nid].append(pid)
                eid_to_nid_map[eid] = node_ids[:8]
                elem = vtkHexahedron()
                elem.GetPointIds().SetId(0, nid_map[node_ids[0]])
                elem.GetPointIds().SetId(1, nid_map[node_ids[1]])
                elem.GetPointIds().SetId(2, nid_map[node_ids[2]])
                elem.GetPointIds().SetId(3, nid_map[node_ids[3]])
                elem.GetPointIds().SetId(4, nid_map[node_ids[4]])
                elem.GetPointIds().SetId(5, nid_map[node_ids[5]])
                elem.GetPointIds().SetId(6, nid_map[node_ids[6]])
                elem.GetPointIds().SetId(7, nid_map[node_ids[7]])
                grid.InsertNextCell(12, elem.GetPointIds())
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

                eid_to_nid_map[eid] = node_ids[:8]
                elem.GetPointIds().SetId(0, nid_map[node_ids[0]])
                elem.GetPointIds().SetId(1, nid_map[node_ids[1]])
                elem.GetPointIds().SetId(2, nid_map[node_ids[2]])
                elem.GetPointIds().SetId(3, nid_map[node_ids[3]])
                elem.GetPointIds().SetId(4, nid_map[node_ids[4]])
                elem.GetPointIds().SetId(5, nid_map[node_ids[5]])
                elem.GetPointIds().SetId(6, nid_map[node_ids[6]])
                elem.GetPointIds().SetId(7, nid_map[node_ids[7]])
                grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                min_thetai, max_thetai, dideal_thetai, min_edge_lengthi = get_min_max_theta(
                    _chexa_faces, node_ids[:8], nid_map, xyz_cid0)

            elif isinstance(element, CPYRAM5):
                node_ids = element.node_ids
                pid = element.Pid()
                for nid in node_ids:
                    nid_to_pid_map[nid].append(pid)
                eid_to_nid_map[eid] = node_ids[:5]
                elem = vtkPyramid()
                elem.GetPointIds().SetId(0, nid_map[node_ids[0]])
                elem.GetPointIds().SetId(1, nid_map[node_ids[1]])
                elem.GetPointIds().SetId(2, nid_map[node_ids[2]])
                elem.GetPointIds().SetId(3, nid_map[node_ids[3]])
                elem.GetPointIds().SetId(4, nid_map[node_ids[4]])
                # etype = 14
                grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
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

                eid_to_nid_map[eid] = node_ids[:5]

                elem.GetPointIds().SetId(0, nid_map[node_ids[0]])
                elem.GetPointIds().SetId(1, nid_map[node_ids[1]])
                elem.GetPointIds().SetId(2, nid_map[node_ids[2]])
                elem.GetPointIds().SetId(3, nid_map[node_ids[3]])
                elem.GetPointIds().SetId(4, nid_map[node_ids[4]])
                grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                min_thetai, max_thetai, dideal_thetai, min_edge_lengthi = get_min_max_theta(
                    _cpyram_faces, node_ids[:5], nid_map, xyz_cid0)

            elif etype in ['CBUSH', 'CBUSH1D', 'CFAST',
                           'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
                           'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
                           'CVISC', 'CGAP']:

                # TODO: verify
                # CBUSH, CBUSH1D, CFAST, CELAS1, CELAS3
                # CDAMP1, CDAMP3, CDAMP4, CDAMP5, CVISC
                if hasattr(element, 'pid'):
                    pid = element.pid
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
                    eid_to_nid_map[eid] = node_ids[slot]
                    nid = node_ids[slot]
                    if nid not in nid_map:
                        # SPOINT
                        print('removing CELASx eid=%i -> SPOINT %i' % (eid, nid))
                        continue

                    #c = nid_map[nid]

                    #if 1:
                    elem = vtk.vtkVertex()
                    elem.GetPointIds().SetId(0, j)
                    #else:
                        #elem = vtk.vtkSphere()
                        #elem = vtk.vtkSphereSource()
                        #if d == 0.:
                        #d = sphere_size
                        #elem.SetRadius(sphere_size)
                else:
                    # 2 points
                    #d = norm(element.nodes[0].get_position() - element.nodes[1].get_position())
                    eid_to_nid_map[eid] = node_ids
                    elem = vtk.vtkLine()
                    try:
                        elem.GetPointIds().SetId(0, nid_map[node_ids[0]])
                        elem.GetPointIds().SetId(1, nid_map[node_ids[1]])
                    except KeyError:
                        print("node_ids =", node_ids)
                        print(str(element))
                        continue

                grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())

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
                #min_edge_lengthi = norm(element.nodes_ref[0].get_position() -
                                        #element.nodes_ref[1].get_position())
                n1, n2 = np.searchsorted(nids, element.nodes)
                xyz1 = xyz_cid0[n1, :]
                xyz2 = xyz_cid0[n2, :]
                min_edge_lengthi = norm(xyz2 - xyz1)
                eid_to_nid_map[eid] = node_ids
                elem = vtk.vtkLine()
                try:
                    elem.GetPointIds().SetId(0, nid_map[node_ids[0]])
                    elem.GetPointIds().SetId(1, nid_map[node_ids[1]])
                except KeyError:
                    print("node_ids =", node_ids)
                    print(str(element))
                    continue
                grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())

            elif etype == 'CBEND':
                pid = element.Pid()
                node_ids = element.node_ids
                for nid in node_ids:
                    nid_to_pid_map[nid].append(pid)

                # 2 points
                n1, n2 = np.searchsorted(nids, element.nodes)
                xyz1 = xyz_cid0[n1, :]
                xyz2 = xyz_cid0[n2, :]
                #min_edge_lengthi = norm(element.nodes_ref[0].get_position() -
                                        #element.nodes_ref[1].get_position())
                eid_to_nid_map[eid] = node_ids

                g0 = element.g0 #_vector
                if not isinstance(g0, integer_types):
                    msg = 'CBEND: g0 must be an integer; g0=%s x=%s\n%s' % (
                        g0, element.x, element)
                    raise NotImplementedError(msg)
                # only supports g0 as an integer
                elem = vtk.vtkQuadraticEdge()
                elem.GetPointIds().SetId(0, nid_map[node_ids[0]])
                elem.GetPointIds().SetId(1, nid_map[node_ids[1]])
                elem.GetPointIds().SetId(2, nid_map[g0])
                grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())

            elif etype == 'CHBDYG':
                node_ids = element.node_ids
                pid = 0
                #pid = element.Pid()
                for nid in node_ids:
                    if nid is not None:
                        nid_to_pid_map[nid].append(pid)

                if element.surface_type in ['AREA4', 'AREA8']:
                    eid_to_nid_map[eid] = node_ids[:4]

                    n1, n2, n3, n4 = [nid_map[nid] for nid in node_ids[:4]]
                    p1 = xyz_cid0[n1, :]
                    p2 = xyz_cid0[n2, :]
                    p3 = xyz_cid0[n3, :]
                    p4 = xyz_cid0[n4, :]
                    out = quad_quality(p1, p2, p3, p4)
                    (areai, taper_ratioi, area_ratioi, max_skew, aspect_ratio,
                     min_thetai, max_thetai, dideal_thetai, min_edge_lengthi) = out
                    if element.surface_type == 'AREA4' or None in node_ids:
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
                    grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                elif element.surface_type in ['AREA3', 'AREA6']:
                    eid_to_nid_map[eid] = node_ids[:3]
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
                    grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                else:
                    #print('removing\n%s' % (element))
                    print('removing eid=%s; %s' % (eid, element.type))
                    del self.eid_map[eid]
                    self.gui.log_info("skipping %s" % element.type)
                    continue
            #elif etype == 'CBYDYP':
            elif etype == 'CHBDYE':
                eid_solid = element.eid2
                side = element.side
                element_solid = model.elements[eid_solid]

                try:
                    mapped_inids = SIDE_MAP[element_solid.type][side]
                except KeyError:  # pragma: no cover
                    print('removing\n%s' % (element))
                    print('removing eid=%s; %s' % (eid, element.type))
                    del self.eid_map[eid]
                    self.gui.log_info("skipping %s" % element.type)
                    continue
                side_inids = [nid - 1 for nid in mapped_inids]
                nodes = element_solid.node_ids

                pid = 0
                unused_nnodes = len(side_inids)
                node_ids = [nodes[inid] for inid in side_inids]
                #inids = np.searchsorted(all_nids, node_ids)

                if len(side_inids) == 3:
                    n1, n2, n3 = [nid_map[nid] for nid in node_ids[:3]]
                    p1 = xyz_cid0[n1, :]
                    p2 = xyz_cid0[n2, :]
                    p3 = xyz_cid0[n3, :]
                    out = tri_quality(p1, p2, p3)
                    (areai, max_skew, aspect_ratio,
                     min_thetai, max_thetai, dideal_thetai, min_edge_lengthi) = out

                    elem = vtkTriangle()
                    elem.GetPointIds().SetId(0, n1)
                    elem.GetPointIds().SetId(1, n2)
                    elem.GetPointIds().SetId(2, n3)
                elif len(side_inids) == 4:
                    n1, n2, n3, n4 = [nid_map[nid] for nid in node_ids[:4]]
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
                else:
                    msg = 'element_solid:\n%s' % (str(element_solid))
                    msg += 'mapped_inids = %s\n' % mapped_inids
                    msg += 'side_inids = %s\n' % side_inids
                    msg += 'nodes = %s\n' % nodes
                    #msg += 'side_nodes = %s\n' % side_nodes
                    raise NotImplementedError(msg)
                grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            else:
                print('removing\n%s' % (element))
                print('removing eid=%s; %s' % (eid, element.type))
                del self.eid_map[eid]
                self.gui.log_info("skipping %s" % element.type)
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

            if np.isnan(max_thetai) and etype not in NO_THETA:
                print('eid=%s theta=%s...setting to 360. deg' % (eid, max_thetai))
                print(element.rstrip())
                if isinstance(element.nodes[0], integer_types):
                    print('  nodes = %s' % element.nodes)
                else:
                    for node in element.nodes:
                        print(str(node).rstrip())
                max_thetai = 2 * np.pi
            #print(eid, min_thetai, max_thetai, '\n', element)

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
        self.gui.nelements = nelements
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
        Creates some nastran specific results

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
        colormap = self.gui.settings.colormap
        #ielement = 0
        nelements = self.element_ids.shape[0]
        normals = np.zeros((nelements, 3), dtype='float32')
        offset = np.full(nelements, np.nan, dtype='float32')
        xoffset = np.full(nelements, np.nan, dtype='float32')
        yoffset = np.full(nelements, np.nan, dtype='float32')
        zoffset = np.full(nelements, np.nan, dtype='float32')
        element_dim = np.full(nelements, -1, dtype='int32')
        nnodes_array = np.full(nelements, np.nan, dtype='int32')

        eid_map = self.gui.eid_map
        assert eid_map is not None
        for eid, element in sorted(iteritems(model.elements)):
            etype = element.type
            if isinstance(element, ShellElement):
                ieid = None
                element_dimi = 2
                #assert element.nodes_ref is not None, element.nodes_ref
                try:
                    normali = element.Normal()
                except AttributeError:
                    msg = str(element)
                    msg += 'nodes_ref = %s\n' % element
                    msg += 'nodes = %s' % element.node_ids
                    raise AttributeError(msg)
                except RuntimeError:
                    # this happens when you have a degenerate tri
                    msg = 'eid=%i normal=nan...setting to [2, 2, 2]\n'
                    msg += '%s' % (element)
                    msg += 'nodes = %s' % str(element.nodes)
                    self.log.error(msg)
                    normali = np.ones(3) * 2.
                    #raise

                prop = element.pid_ref
                if prop is None:
                    # F:\work\pyNastran\examples\Dropbox\move_tpl\ehbus69.op2
                    # CTRIAX
                    z0 = np.nan
                else:
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

                ieid = eid_map[eid]
                normals[ieid, :] = normali
                if element.type in ['CPLSTN3', 'CPLSTN4', 'CPLSTN6', 'CPLSTN8']:
                    element_dim[ieid] = element_dimi
                    nnodes_array[ieid] = nnodesi
                    self.log.debug('continue...element.type=%r' % element.type)
                    continue

                offset[ieid] = z0
                xoffset[ieid] = z0 * normali[0]
                yoffset[ieid] = z0 * normali[1]
                zoffset[ieid] = z0 * normali[2]

            elif etype == 'CTETRA':
                ieid = self.eid_map[eid]
                element_dimi = 3
                nnodesi = 4
            elif etype == 'CPENTA':
                ieid = self.eid_map[eid]
                element_dimi = 3
                nnodesi = 6
            elif etype == 'CPYRAM':
                ieid = self.eid_map[eid]
                element_dimi = 3
                nnodesi = 5
            elif etype in ['CHEXA', 'CIHEX1', 'CIHEX2']:
                ieid = self.eid_map[eid]
                element_dimi = 3
                nnodesi = 8

            elif etype in ['CROD', 'CONROD', 'CBEND', 'CBAR', 'CBEAM', 'CGAP', 'CTUBE']:
                ieid = self.eid_map[eid]
                element_dimi = 1
                nnodesi = 2
            elif etype in ['CBUSH', 'CBUSH1D', 'CBUSH2D',
                           'CFAST', 'CVISC',
                           'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
                           'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5']:
                ieid = self.eid_map[eid]
                element_dimi = 0
                nnodesi = 2
            elif etype == 'CHBDYG':
                ieid = self.eid_map[eid]
                if element.surface_type == 'AREA3':
                    nnodesi = 3
                    element_dimi = 2
                elif element.surface_type == 'AREA4':
                    nnodesi = 4
                    element_dimi = 2
                elif element.surface_type == 'AREA6':
                    nnodesi = 6
                    element_dimi = 2
                elif element.surface_type == 'AREA8':
                    nnodesi = 8
                    element_dimi = 2
                #elif element.surface_type == 'REV':
                    #nnodesi = 2 # ???
                    #element_dimi = 1 # ???
                else:
                    element_dimi = -1
                    nnodesi = -1
                    print('element.type=%s doesnt have a dimension' % element.type)
            elif etype == 'CHBDYP':
                continue
            else:
                try:
                    ieid = self.eid_map[eid]
                except KeyError:
                    print('didnt add element to eid_map')
                    print(model.elements[eid])
                    raise
                element_dimi = -1
                nnodesi = -1
                print('element.type=%s doesnt have a dimension' % element.type)
            assert ieid is not None
            element_dim[ieid] = element_dimi
            nnodes_array[ieid] = nnodesi
            #ielement += 1

        # if not a flat plate
        #if min(nxs) == max(nxs) and min(nxs) != 0.0:
        is_element_dim = np.max(element_dim) != np.min(element_dim)
        is_element_dim = True
        if is_element_dim and isfinite_and_greater_than(element_dim, -1):
            eid_dim_res = GuiResult(0, header='ElementDim', title='ElementDim',
                                    location='centroid', scalar=element_dim, mask_value=-1)
            cases[icase] = (eid_dim_res, (0, 'ElementDim'))

        is_shell = np.abs(normals).max() > 0.

        # we have to add the 2nd/3rd lines to make sure bars are getting into this check

        is_solid = (
            isfinite_and_nonzero(min_interior_angle) and
            isfinite_and_nonzero(max_interior_angle)
        )

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
                                    colormap=colormap, data_format='%.1f',
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

            if np.any(np.isfinite(area_ratio)) and np.nanmax(area_ratio) > 1.:
                arearatio_res = GuiResult(
                    0, header='Area Ratio', title='Area Ratio',
                    location='centroid', scalar=area_ratio)
                cases[icase] = (arearatio_res, (0, 'Area Ratio'))
                form_checks.append(('Area Ratio', icase, []))
                icase += 1

            if np.any(np.isfinite(taper_ratio)) and np.nanmax(taper_ratio) > 1.:
                taperratio_res = GuiResult(
                    0, header='Taper Ratio', title='Taper Ratio',
                    location='centroid', scalar=taper_ratio)
                cases[icase] = (taperratio_res, (0, 'Taper Ratio'))
                form_checks.append(('Taper Ratio', icase, []))
                icase += 1

            if isfinite_and_nonzero(max_warp_angle):
                warp_res = GuiResult(
                    0, header='Max Warp Angle', title='MaxWarpAngle',
                    location='centroid', scalar=np.degrees(max_warp_angle))
                cases[icase + 4] = (warp_res, (0, 'Max Warp Angle'))
                form_checks.append(('Max Warp Angle', icase, []))
                icase += 1

            #if (np.abs(xoffset).max() > 0.0 or np.abs(yoffset).max() > 0.0 or
                #np.abs(zoffset).max() > 0.0):
            # offsets
            if isfinite(max_warp_angle):
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

            if 0:  # pragma: no cover
                xyz_offset = np.vstack([xoffset, yoffset, zoffset]).T
                titles = ['Offset XYZ']
                headers = titles
                assert xyz_offset.shape[1] == 3, xyz_offset.shape
                assert xyz_offset.shape[0] == len(offset)
                scales = [1.0]
                subcase_id = 0
                offset_xyz_res = ElementalTableResults(
                    subcase_id, titles, headers, xyz_offset, offset, scales)
                offset_xyz_res.save_defaults()
                cases[icase] = (offset_z_res, (0, 'OffsetZ'))
                form_checks.append(('OffsetXYZ', icase, []))
                icase += 1

            if self.make_xyz or IS_TESTING:
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

        if material_coord is not None and np.abs(material_coord).max() > 0:
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
        upids = None
        pcomp = None
        pshell = None
        is_pcomp = False
        is_pshell = False

        mids_pcomp = None
        thickness_pcomp = None
        nplies_pcomp = None
        pcomp = {
            'mids' : mids_pcomp,
            'thickness' : thickness_pcomp,
            'nplies' : nplies_pcomp,
        }

        mids = None
        thickness = None
        pshell = {
            'mids' : mids,
            'thickness' : thickness,
        }
        if not isfinite_and_greater_than(pids, 0):
            return icase, upids, pcomp, pshell, (is_pshell, is_pcomp)

        prop_types_with_mid = [
            'PSOLID', 'PSHEAR',
            'PROD', 'PTUBE', 'PBAR', 'PBARL', 'PBEAM', 'PBEAML',
            'PBEND',
        ]
        prop_types_without_mid = ['PVISC', 'PELAS', 'PBUSH', 'PDAMP', 'PDAMPT']

        pid_res = GuiResult(0, header='PropertyID', title='PropertyID',
                            location='centroid', scalar=pids, mask_value=0)
        cases[icase] = (pid_res, (0, 'PropertyID'))
        form0.append(('PropertyID', icase, []))
        icase += 1
        upids = np.unique(pids)

        mid_eids_skip = []

        #mids_pshell = None
        #thickness_pshell = None
        if 'PSHELL' in model.card_count:
            is_pshell = True

        pids_pcomp = model.get_card_ids_by_card_types(['PCOMP', 'PCOMPG'], combine=True)
        if pids_pcomp:
            npliesi = 0
            pcomp_nplies = 0
            for pid in pids_pcomp:
                prop = model.properties[pid]
                pcomp_nplies = max(pcomp_nplies, prop.nplies + 1)
            npliesi = max(npliesi, pcomp_nplies)

            nplies_pcomp = np.zeros(nelements, dtype='int32')
            mids_pcomp = np.zeros((nelements, npliesi), dtype='int32')
            thickness_pcomp = np.full((nelements, npliesi), np.nan, dtype='float32')
            mids_pcomp = np.zeros((nelements, npliesi), dtype='int32')
            is_pcomp = True

        #rho = np.full((nelements, nplies), np.nan, dtype='float32')
        mids = np.zeros((nelements, 4), dtype='int32')
        thickness = np.full((nelements, 4), np.nan, dtype='float32')
        for pid in upids:
            if pid == 0:
                print('skipping pid=0')
                continue
            elif pid < 0:
                continue

            prop = model.properties[pid]
            #try:
            if prop.type in prop_types_with_mid:
                # simple types
                i = np.where(pids == pid)[0]
                mid = prop.mid_ref.mid
                mids[i, 0] = mid
            elif prop.type == 'PSHELL':
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
                thickness[i, 1] = prop.twelveIt3
                thickness[i, 2] = prop.tst

            elif prop.type in ['PCOMP', 'PCOMPG']:
                i = np.where(pids == pid)[0]
                npliesi = prop.nplies
                nplies_pcomp[i] = npliesi
                thickness_pcomp[i, 0] = 0.
                for iply in range(npliesi):
                    mids_pcomp[i, iply+1] = prop.Mid(iply)
                    thickniess_ply = prop.Thickness(iply)
                    thickness_pcomp[i, iply+1] = thickniess_ply
                    thickness_pcomp[i, 0] += thickniess_ply
                #mids[i, 0] = mids[i, 1]

            #elif prop.type == 'PSHEAR':  # element has the thickness
                #i = np.where(pids == pid)[0]
                #mids[i, 0] = prop.Mid()
                #thickness[i, 0] = elem.Thickness()
            elif prop.type in prop_types_without_mid:
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

        pcomp = {
            'mids' : mids_pcomp,
            'thickness' : thickness_pcomp,
            'nplies' : nplies_pcomp,
        }
        pshell = {
            'mids' : mids,
            'thickness' : thickness,
        }
        return icase, upids, pcomp, pshell, (is_pshell, is_pcomp)

    def _build_materials(self, model, pcomp, pshell, is_pshell_pcomp,
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
        for i, pshell_pcompi in enumerate([pshell, pcomp]):
            mids = pshell_pcompi['mids']
            thickness = pshell_pcompi['thickness']

            if 'nplies' in pshell_pcompi:
                nplies = pshell_pcompi['nplies']
                if nplies is not None and nplies.max() > 0:
                    nplies_res = GuiResult(0, header='Number of Plies', title='nPlies',
                                           location='centroid', scalar=nplies, mask_value=0)
                    cases[icase] = (nplies_res, (0, 'Number of Plies'))
                    form0.append(('Number of Plies', icase, []))
                    icase += 1

            if mids is None:
                continue
            nlayers = mids.shape[1]
            for ilayer in range(nlayers):
                if len(thickness.shape) == 2:
                    thicknessi = thickness[:, ilayer]
                else:
                    ## TODO: I think this is used by a non-PSHELL/PCOMP case
                    #print('B-shape...i=%s ilayer=%s' % (i, ilayer))
                    thicknessi = thickness

                form_layer = []
                #if i == 1 and ilayer == 0:
                    #print('thicknessi = ', thicknessi)
                if isfinite_and_nonzero(thicknessi):
                    if i == 1 and ilayer == 0:
                        tword = 'Total Thickness'  # thickness is nan
                    elif i == 0 and ilayer == 1:
                        tword = '12/t^3'
                    elif i == 0 and ilayer == 2:
                        tword = 'ts/t'
                    elif i == 0 and ilayer == 3:
                        tword = 'mid4'
                    else:
                        tword = 'Thickness'
                    if tword != 'mid4':
                        t_res = GuiResult(0, header=tword, title=tword,
                                          location='centroid', scalar=thicknessi)
                        cases[icase] = (t_res, (0, tword))
                        form_layer.append((tword, icase, []))
                        icase += 1

                midsi = mids[:, ilayer]
                if midsi.max() == 0:
                    if not(i == 1 and ilayer == 0):
                        print('cant find anything in ilayer=%s' % ilayer)
                    #continue
                else:
                    has_mat8, has_mat11, e11, e22, e33 = get_material_arrays(model, midsi)
                    mid_res = GuiResult(0, header='MaterialID', title='MaterialID',
                                        location='centroid', scalar=midsi, mask_value=0)
                    cases[icase] = (mid_res, (0, 'MaterialID'))
                    form_layer.append(('MaterialID', icase, []))
                    icase += 1

                    if has_mat11: # also implicitly has_mat8
                        is_orthotropic = not (np.array_equal(e11, e22) and np.array_equal(e11, e33))
                    elif has_mat8:
                        is_orthotropic = not np.array_equal(e11, e22)
                    else:
                        is_orthotropic = False

                    # np.nanmax(e11) > 0. can fail if e11=[nan, nan]
                    e112 = np.fmax.reduce(e11)
                    is_e11 = True
                    if np.isnan(e112):
                        is_e11 = False
                        #
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
                        if has_mat11:
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
                    elif is_e11:
                        # isotropic
                        assert np.nanmax(e11) > 0, np.nanmax(e11)
                        e11_res = GuiResult(0, header='E', title='E',
                                            location='centroid', scalar=e11, data_format='%.3e')
                        cases[icase] = (e11_res, (0, 'E'))
                        form_layer.append(('E', icase, []))
                        icase += 1

                #print('form_layer =', form_layer)
                if form_layer:
                    if nlayers == 1:
                        form0 += form_layer
                    else:
                        word = _get_nastran_gui_layer_word(i, ilayer, is_pshell_pcomp)
                        form0.append((word, None, form_layer))
        return icase

    def _build_optimization(self, model, pids, upids, nelements, cases, form0, icase):
        """
        Creates the optimization visualization.  Supports:
          - DVPREL1/2 shell thickness:
             - DV Region
             - DVPREL Init - t
             - DVPREL Min - t
             - DVPREL Max - t
        """
        if upids is None:
            return icase
        if len(model.properties) and len(model.dvprels):
            # len(model.dvprels) + len(model.dvcrels) + len(model.dvmrels) + len(model.desvars)
            #dvmrel_init = np.zeros(nelements, dtype='int32')
            #dvgrel_init = np.zeros(nelements, dtype='int32')
            out_dict = model._get_dvprel_ndarrays(nelements, pids)

            optimization_cases = []
            for key, dvprel_data in out_dict.items():
                design_region, dvprel_init, dvprel_min, dvprel_max = dvprel_data
                if np.nanmax(design_region) == 0:
                    continue

                region_res = GuiResult(
                    0, header='DV Region', title='DV Region',
                    location='centroid', scalar=design_region)
                t_init_res = GuiResult(
                    0, header='DVPREL Init - %s' % key, title='DVPREL Init - %s' % key,
                    location='centroid', scalar=dvprel_init)
                t_min_res = GuiResult(
                    0, header='DVPREL Min - %s' % key, title='DVPREL Min - %s' % key,
                    location='centroid', scalar=dvprel_min)
                t_max_res = GuiResult(
                    0, header='DVPREL Max - %s' % key, title='DVPREL Max - %s' % key,
                    location='centroid', scalar=dvprel_max)
                cases[icase] = (region_res, (0, 'DV Region'))
                cases[icase + 1] = (t_init_res, (0, 'DVPREL Init - %s' % key))
                cases[icase + 2] = (t_min_res, (0, 'DVPREL Min - %s' % key))
                cases[icase + 3] = (t_max_res, (0, 'DVPREL Max - %s' % key))
                opt_cases = []
                opt_cases.append(('DV Region', icase, []))
                opt_cases.append(('DVPREL Init - %s' % key, icase + 1, []))
                opt_cases.append(('DVPREL Min - %s' % key, icase + 2, []))
                opt_cases.append(('DVPREL Max - %s' % key, icase + 3, []))
                optimization_cases.append((key, None, opt_cases))
                icase += 4
            if optimization_cases:
                form0.append(('Optimization', None, optimization_cases))
        return icase

    def _plot_pressures(self, model, cases, form0, icase, subcase_id):
        """
        pressure act normal to a shell (as opposed to anti-normal to a solid face)
        """
        # quit out if we're going to make pressure plots anyways
        #if self.plot_applied_loads:
            #return icase

        # quit out if we don't have pressures
        if not any(['PLOAD' in model.card_count, 'PLOAD2' in model.card_count,
                    'PLOAD4' in model.card_count]):
            return icase

        subcase = model.subcases[subcase_id]

        try:
            load_case_id = subcase.get_parameter('LOAD')[0]
        except KeyError:
            self.gui.log.warning('LOAD not found in subcase_id=%s' % (subcase_id))
            return icase

        if load_case_id not in model.loads and load_case_id not in model.load_combinations:
            self.gui.log.warning('LOAD=%s not found' % load_case_id)
            return icase

        is_pressure, pressures = model.get_pressure_array(
            load_case_id, eids=self.element_ids, stop_on_failure=False)
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
        #if not self.plot_applied_loads:
            #model.log.debug('self.plot_applied_loads=False')
            #return icase

        if not xref_loads:
            model.log.debug('returning from plot_applied_loads_early')
            return icase

        colormap = self.gui.settings.colormap
        try:
            #form = []
            out = model.get_load_arrays(
                subcase_id,
                eid_map=self.eid_map, node_ids=self.node_ids,
                normals=self.normals, nid_map=self.nid_map,)
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
                            nlabels=None, labelsize=None, ncolors=None, colormap=colormap,
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
                            nlabels=None, labelsize=None, ncolors=None, colormap=colormap,
                            set_max_min=False, uname='NastranGeometry')
                        moment_xyz_res.save_defaults()

                        cases[icase] = (moment_xyz_res, (0, 'Moment XYZ'))
                        form0.append(('Moment XYZ', icase, []))
                        icase += 1

                if np.abs(spcd.max() - spcd.min()) > 0.0:
                    t123 = spcd[:, :3]
                    tnorm = norm(t123, axis=1)
                    assert len(tnorm) == len(spcd[:, 2]), len(spcd[:, 2])
                    assert len(tnorm) == len(self.nid_map)

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
                assert len(temperatures) == len(self.nid_map)
                temperature_res = GuiResult(
                    subcase_id, header=temperature_key, title=temperature_key,
                    location='node', scalar=temperatures)
                cases[icase] = (temperature_res, (0, temperature_key))
                form0.append((temperature_key, icase, []))
                icase += 1
        except KeyError:
            stringio = StringIO()
            traceback.print_exc(file=stringio)
            sout = stringio.getvalue()
            self.gui.log_error(sout)
            print(sout)
        return icase

    def load_nastran_results(self, results_filename):
        """
        Loads the Nastran results into the GUI
        """
        self.scalarBar.VisibilityOn()
        self.scalarBar.Modified()

        if isinstance(results_filename, string_types):
            print("trying to read...%s" % results_filename)
            ext = os.path.splitext(results_filename)[1].lower()

            if ext == '.op2':
                op2_filename = results_filename
                model = OP2(log=self.log, debug=True)

                if 0:  # pragma: no cover
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

                if not IS_TESTING or self.is_testing_flag:
                    self.log.info(model.get_op2_stats())
                # print(model.get_op2_stats())

            elif ext == '.nod':
                self.load_patran_nod(results_filename)
                self.cycle_results_explicit()  # start at nCase=0
                return
            elif ext == '.h5' and IS_H5PY:
                model = OP2(log=self.log, debug=True)
                hdf5_filename = results_filename
                model.load_hdf5(hdf5_filename, combine=False)
            #elif ext == '.pch':
                #raise NotImplementedError('*.pch is not implemented; filename=%r' % op2_filename)
            #elif ext == '.f06':
                #model = F06(log=self.log, debug=True)
                #model.set_vectorization(True)
                #model.read_f06(op2_filename)
            else:
                #print("error...")
                msg = 'extension=%r is not supported; filename=%r' % (ext, op2_filename)
                raise NotImplementedError(msg)
        else:
            model = op2_filename
            op2_filename = op2_filename.filename

        if self.save_data:
            self.model_results = model

        #print(model.print_results())
        #self.isubcase_name_map[self.isubcase] = [Subtitle, Label]

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
            #self.isubcase_name_map = {}
            #form = []
            #icase = 0
        #else:
        cases = self.result_cases
        form = self.get_form()
        icase = len(cases)
        # form = self.res_widget.get_form()

        #subcase_ids = model.isubcase_name_map.keys()
        #self.isubcase_name_map = model.isubcase_name_map
        # self.isubcase_name_map = model.subcase_key
        #print(self.isubcase_name_map)
        for isubcase, values in model.isubcase_name_map.items():
            if not isinstance(isubcase, integer_types):
                print('isubcase type =', type(isubcase))
                continue
            if isinstance(values, str):
                # eigenvalue???
                label = values
                self.gui.log.debug('label_str = %r' % label)
            elif isinstance(values, list):
                self.gui.log.debug(str(values))
                subtitle, superelement_adaptivity, analysis_code, label = values
                del analysis_code
            else:
                self.gui.log.debug(str(values))
                self.gui.log.debug(str(type(values)))
                raise RuntimeError(values)

            if superelement_adaptivity:
                subcase_name = '%s: %s' % (subtitle, superelement_adaptivity)
            else:
                subcase_name = subtitle
            self.isubcase_name_map[isubcase] = [subcase_name, label]
            del subtitle, label
        # self.isubcase_name_map = {subcase_id : label for
                                # in iteritems(model.isubcase_name_map)}

        form = self._fill_op2_output(results_filename, cases, model, form, icase)
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
        keys = model.get_key_order()
        assert keys is not None, keys
        #print('keys_order =', keys)

        disp_dict = defaultdict(list)
        stress_dict = defaultdict(list)
        strain_dict = defaultdict(list)
        force_dict = defaultdict(list)
        strain_energy_dict = defaultdict(list)

        header_dict = {}
        keys_map = {}
        key_itime = []

        for key in keys:
            unused_is_data, unused_is_static, unused_is_real, times = _get_times(model, key)
            if times is None:
                # we dynamically created the keys and created extra ones
                continue
            #assert times is not None  # gen22x_modes

            #print('--------------')
            #print('key = %r' % str(key))
            self.stress[key] = StressObject(model, key, self.element_ids, is_stress=True)
            self.strain[key] = StressObject(model, key, self.element_ids, is_stress=False)

            #header_dict[(key, 0)] = '; Static'

            unused_formi = []
            unused_form_time = []

            ncases_old = icase
            icase = self._fill_op2_oug_oqg(cases, model, key, icase,
                                           disp_dict, header_dict, keys_map)
            ncases = icase - ncases_old
            #print('ncases=%s icase=%s' % (ncases, icase))
            #assert ncases > 0, ncases

            if ncases:
                # can potentially make a double listing, but we need it
                # eigenvector only cases
                for itime, unused_dt in enumerate(times):
                    new_key = (key, itime)
                    key_itime.append(new_key)


            for itime, unused_dt in enumerate(times):
                ncases_old = icase
                # stress
                icase = self._fill_op2_stress(
                    cases, model, key, icase, itime,
                    stress_dict, header_dict, keys_map)

                # strain
                icase = self._fill_op2_strain(
                    cases, model, key, icase, itime,
                    strain_dict, header_dict, keys_map)

                # force
                icase = self._fill_op2_force(
                    cases, model, key, icase, itime,
                    force_dict, header_dict, keys_map)

                # strain energy
                icase = self._fill_op2_time_centroidal_strain_energy(
                    cases, model, key, icase, itime,
                    strain_energy_dict, header_dict, keys_map)
                ncases = icase - ncases_old
                new_key = (key, itime)
                if ncases and new_key not in key_itime:
                    key_itime.append(new_key)
            #print('******form_time =', form_time)
            #print(header)

        #print('Key,itime:')
        #for key_itimei in key_itime:
            #print('  %s' % str(key_itimei))

        unused_form_out = []
        is_results = False

        form_resultsi = []
        form_resultsi_subcase = []
        basename = os.path.basename(op2_filename).rstrip()
        form_results = (basename + '-Results', None, form_resultsi)

        if len(key_itime) == 0:
            print('header_dict =', header_dict)
            print('key_itime =', key_itime)
            return form

        # (isubcase, analysis_code, sort_method,
        #  count, ogs, superelement_adaptivity_index) = key
        key_itime0 = key_itime[0]
        key0 = key_itime0[0]
        # (isubcase, analysis_code, sort_method,
        #  count, ogs, superelement_adaptivity_index, pval_step) = key
        subcase_id_old = key0[0]
        count_old = key0[3]
        unused_ogs_old = key0[4]
        subtitle_old = key0[5]
        subtitle_old, label_old, superelement_adaptivity_index_old, unused_pval_step_old = keys_map[key0]
        del label_old
        del superelement_adaptivity_index_old

        # now that we have the data built, we put it in the form
        # in sorted order
        #
        # TODO: consider pval_step
        for key, itime in key_itime:
            # (isubcase, analysis_code, sort_method,
            #  count, ogs, superelement_adaptivity_index, pval_step) = key
            #print('key =', key)
            subcase_id = key[0]
            count = key[3]
            #subtitle = key[4]
            try:
                subtitle, unused_label, superelement_adaptivity_index, unused_pval_step = keys_map[key]
            except:
                subcase_id = subcase_id_old
                subtitle = subtitle_old + '?'
                superelement_adaptivity_index = '?'
                raise

            if subcase_id != subcase_id_old or subtitle != subtitle_old:
                count_str = '' if count == 0 else ' ; opt_count=%s' % count_old
                subcase_str = 'Subcase %s; %s%s%s' % (
                    subcase_id_old, subtitle_old, superelement_adaptivity_index, count_str)
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
                form_outi.append(('Strain', None, strain_formi))
                is_results = True
            if force_formi:
                form_outi.append(('Force', None, force_formi))
                is_results = True
            if strain_energy_formi:
                form_outi.append(('Strain Energy', None, strain_energy_formi))
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
        #print('form', form)
        #print('form_results =', form_results)
        return form

    def clear_nastran(self):
        """cleans up variables specific to Nastran"""
        self.eid_map = {}
        self.nid_map = {}
        self.eid_to_nid_map = {}
        self.element_ids = None
        self.node_ids = None

def _get_shell_material_coord_int(element):
    """
    used by:
     - CQUAD4, CQUADR, CQUAD8, CQUAD, CQUADX
     - CTRIA3, CTRIAR, CTRIA6
    """
    return 0 if isinstance(element.theta_mcid, float) else element.theta_mcid

def _get_nastran_gui_layer_word(i, ilayer, is_pshell_pcomp):
    """gets the PSHELL/PCOMP layer word"""
    ## TODO: this makes no sense...
    is_pshell, unused_is_pcomp = is_pshell_pcomp
    word = ''
    if i == 0:
        if ilayer == 0:
            if is_pshell:
                word += 'PSHELL: ilayer=1 & others'
            else:
                word += 'Other Properties'
        else:
            word += 'PSHELL: ilayer=%i' % (ilayer + 1)
    else:
        if ilayer == 0:
            word += 'PCOMP: Total'
        else:
            word += 'PCOMP: ilayer=%i' % (ilayer)
    return word

def check_for_missing_control_surface_boxes(name, cs_box_ids,
                                            box_id_to_caero_element_map,
                                            log, store_msg=False):
    """helper method for creating control surface"""
    boxes_to_show = []
    missing_boxes = []
    for box_id in cs_box_ids:
        try:
            unused_ipoints = box_id_to_caero_element_map[box_id]
        except KeyError:
            missing_boxes.append(box_id)
            continue
        boxes_to_show.append(box_id)

    out_msg = ''
    if missing_boxes:
        #print('\nboxes_to_show =', boxes_to_show)
        msg = 'Missing CAERO AELIST/SPLINE control surface %r boxes: %s\n' % (
            name, str(missing_boxes))
        if not boxes_to_show:
            msg += 'boxes_to_show=%s\n' % boxes_to_show

        out_msg = store_error(log, store_msg, msg.rstrip())
    return boxes_to_show, out_msg

def store_error(log, store_msg, msg):
    out_msg = ''
    if store_msg:
        n, filename = log_properties(nframe=2)
        out_msg = 'ERROR: %s:%s %s\n' % (filename, n, msg)
    else:
        log.warning(msg)
    return out_msg

def store_warning(log, store_msg, msg):
    out_msg = ''
    if store_msg:
        n, filename = log_properties(nframe=2)
        out_msg = 'WARNING: %s:%s %s\n' % (filename, n, msg)
    else:
        log.warning(msg)
    return out_msg
