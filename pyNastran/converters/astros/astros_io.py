from typing import Optional, Any, TYPE_CHECKING
import numpy as np
import sys
from itertools import chain
from collections import defaultdict
from io import StringIO
import traceback
from pyNastran.gui.qt_files.colors import RED_FLOAT, BLUE_FLOAT, GREEN_FLOAT, PINK_FLOAT, PURPLE_FLOAT
from pyNastran.gui.vtk_common_core import vtkPoints
from pyNastran.gui.vtk_interface import vtkVertex, vtkLine, vtkUnstructuredGrid, vtkQuad
from pyNastran.gui.utils.vtk.vtk_utils import create_vtk_cells_of_constant_element_type, create_vtk_cells_of_constant_element_types
from pyNastran.bdf.bdf import BDF, CAERO1, CAERO2, CAERO3, CAERO4, CAERO5
from pyNastran.bdf.cards.aero.zona import CAERO7, BODY7
from pyNastran.bdf.cards.aero.aero import get_caero_subpanel_grid
from pyNastran.femutils.nan import isfinite, isfinite_and_greater_than
from pyNastran.femutils.utils import duplicates, is_monotonic
from pyNastran.converters.nastran.gui.nastran_io import (_create_caero_actors, _create_masses, set_acoustic_grid, get_beam_sections_map, get_bar_nids, build_map_centroidal_result,
                                                         numpy_to_vtk_points, _build_materials, _build_optimization, build_normals_quality, _set_conm_grid, build_caero_paneling,
                                                         get_mpc_node_ids, _apply_colon_set, _prepare_superelement_model, plotels_to_groups)
from pyNastran.converters.nastran.gui.nastran_io_utils import (build_superelement_model, map_elements1_no_quality_helper, map_elements1_quality_helper, create_monpnt, get_pcomp_nplies, get_caero_control_surface_grid)
from pyNastran.converters.nastran.gui.utils import make_nid_map, store_warning, check_for_missing_control_surface_boxes 
if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger
    from pyNastran.gui.gui_objects.settings import Settings, NastranSettings
    from pyNastran.gui.main_window import MainWindow
from pyNastran.gui.gui_objects.settings import Settings, NastranSettings
from pyNastran.gui.gui_objects.gui_result import GuiResult, NormalResult
from pyNastran.converters.astros.astros import Astros

IS_TESTING = 'test' in sys.argv[0]

class AstrosIO():
    """Defines the GUI class for Astros."""
    def __init__(self, gui):
        self.gui = gui
        self.create_secondary_actors = True
        self.save_data = True
        self.make_offset_normals_dim = True
        self.nid_release_map = {}
        self.stress = {}
        self.strain = {}

    def get_astros_wildcard_geometry_results_functions(self):
        """dynamic named method for loading astros input files"""
        data = (
            'Astros',
            'Astros (*.in)', self.load_astros_geometry,
            None, self.load_astros_results
        )
        return data
    
    def load_astros_geometry(self, astros_filename: str,
                            name: str='main', plot: bool=True):
        """loads astros input files into the gui"""
        gui: MainWindow = self.gui
        gui.eid_maps[name] = {}
        gui.nid_maps[name] = {}

        self.clear_astros()

        model_name = 'main'
        reset_labels = True
        settings: Settings = gui.settings
        nastran_settings: NastranSettings = settings.nastran_settings
        if plot:
            gui.scalar_bar_actor.VisibilityOff()
            gui.scalar_bar_actor.Modified()

        xref_loads = True # should be True
        model = Astros(log=self.gui.log, debug=False)
        self.gui.model_type = 'astros'
        model.read_astros_inp(astros_filename)
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
        create_superelement_geometry=True,
    )
        self.log = model.log

        # ids = []
        # for id, node in model.nodes.items():
        #     if node.xyz[1] < 0.01:
        #         ids.append(id)


        nnodes = len(model.nodes)
        # nspoints = len(model.spoints)
        # nepoints = len(model.epoints)
        # ngridb = len(model.gridb)
        # ncaero_cards = len(model.caeros)

        # for superelement in model.superelement_models.values():
        #     nnodes += len(superelement.nodes)
        #     nspoints += len(superelement.spoints)
        #     nepoints += len(superelement.epoints)
        #     ngridb += len(superelement.gridb)
        #     ncaero_cards += len(superelement.caeros)

        # ngui_nodes = nnodes + nspoints + nepoints + ngridb
        ngui_nodes = nnodes
        self.gui.nnodes = nnodes
        # if ngui_nodes + ncaero_cards == 0:
        #     msg = 'nnodes + nspoints + nepoints = 0\n'
        #     msg += 'card_count = %r' % str(model.card_count)
        #     raise NoGeometry(msg)

        nelements = len(model.elements)

        # nmasses = len(model.masses)
        # nplotels = len(model.plotels)
        # nrigid = len(model.rigid_elements)
        # for superelement in model.superelement_models.values():
        #     nelements += len(superelement.elements)
        #     nmasses += len(superelement.masses)
        #     nplotels += len(superelement.plotels)
        #     nrigid += len(superelement.rigid_elements)

        #nmpc = len(model.mpcs)  # really should only be allowed if we have it in a subcase
        # if nelements + nmasses + ncaero_cards + nplotels + nrigid == 0:
        #     msg = 'nelements + nmasses + ncaero_cards + nplotels + nrigid = 0\n'
        #     msg += 'card_count = %r' % str(model.card_count)
        #     raise NoGeometry(msg)

        self.nnodes = ngui_nodes
        self.nelements = nelements  # approximate...

        out = self.make_caeros(model)
        (has_caero, caero_points, ncaeros, ncaeros_sub, ncaeros_cs,
         ncaeros_points, ncaero_sub_points,
         has_control_surface, box_id_to_caero_element_map, cs_box_ids) = out
        self.has_caero = has_caero

        #-----------------------------------------------------------------------
        gui.log_info(f'nnodes={self.nnodes:d} nelements={self.nelements:d}')
        # msg = model.get_bdf_stats(return_type='string')
        # gui.log_debug(msg)
        # msg = model.get_bdf_stats(return_type='list')

        # this call will break the GUI if there are a lot of lines and
        # by a lot I mean 37641.  It's fine for a single call.
        #for msgi in msg:
            #model.log.debug(msgi)
        #-----------------------------------------------------------------------
        # nodes/coords

        #print('get_xyz_in_coord')
        dim_max = 1.0
        xyz_cid0 = None
        nid_cp_cd = None
        if gui.nnodes:
            xyz_cid0, nid_cp_cd = self.get_xyz_in_coord(model, cid=0, fdtype='float32')
            xyz_cid0 = gui.scale_length(xyz_cid0)
            dim_max = self._points_to_vtkpoints_coords(model, xyz_cid0)
            self.node_ids = nid_cp_cd[:, 0]
        else:
            self.node_ids = np.zeros(0, dtype='int32')
        gui.node_ids = self.node_ids
        #-----------------------------------------------------------------------
        nconm2 = _create_masses(
            gui, model, self.node_ids,
            create_secondary_actors=self.create_secondary_actors)

        # Allocate grids
        gui.grid.Allocate(self.nelements, 1000)
        _create_caero_actors(gui, ncaeros, ncaeros_sub, ncaeros_cs,
                             has_control_surface, self.has_caero)
        if nconm2 > 0:
            gui.alt_grids['conm2'].Allocate(nconm2, 1000)

        if self.save_data:
            self.model = model

        #-----------------------------------------------------------------------
        j = 0
        nid_map = gui.nid_map
        idtype = self.node_ids.dtype
        self.element_ids = np.zeros(0, dtype='int32')
        nid_to_pid_map, icase, cases, form = self.map_elements(
            xyz_cid0, nid_cp_cd, nid_map, model, j, dim_max,
            plot=plot, xref_loads=xref_loads)

        if xyz_cid0 is not None:
            create_monpnt(self, model, xyz_cid0, nid_cp_cd)
        self._create_aero(model, box_id_to_caero_element_map, cs_box_ids,
                          caero_points, ncaeros_points, ncaero_sub_points,
                          has_control_surface)

        if nconm2 > 0:
            _set_conm_grid(gui, nconm2, model,
                           self.create_secondary_actors)

        geometry_names = []
        make_spc_mpc_supports = nastran_settings.is_constraints
        make_acoustic = True
        if self.create_secondary_actors:
            if make_spc_mpc_supports:
                geometry_names = self.set_spc_mpc_suport_grid(
                    model, nid_to_pid_map, idtype)
            if make_acoustic:
                set_acoustic_grid(gui, model, xyz_cid0, nid_cp_cd, nid_map)

            if nastran_settings.is_bar_axes:
                icase = self._fill_bar_yz(dim_max, model, icase, cases, form)
        assert icase is not None

        #------------------------------------------------------------
        #print('dependent_nodes =', self.dependents_nodes)
        icase = self._set_subcases_unvectorized(model, form, cases, icase,
                                                True, xref_loads)

        name = 'main_copy'
        gui.duplicate_alternate_vtk_grid(
            name, 'main', color=(0., 0., 0.), line_width=5,
            opacity=0.1, is_visible=False)

        #------------------------------------------------------------
        # add alternate actors
        gui._add_alt_actors(gui.alt_grids)

        # set default representation
        self._set_caero_representation(has_control_surface)

        for grid_name in geometry_names:
            if grid_name in gui.geometry_actors:
                gui.geometry_actors[grid_name].Modified()

        #self.grid_mapper.SetResolveCoincidentTopologyToPolygonOffset()
        stop_on_failure = False
        build_map_centroidal_result(
            model, nid_map, stop_on_failure=stop_on_failure)

        #if self.create_secondary_actors and not IS_TESTING and 'dev' in __version__:
            #self.sidebar_nastran = ModelSidebar(self.gui, nastran_io=self)
            #self.sidebar_nastran.set_model(model)

            #self.res_dock_nastran = QDockWidget("Nastran Model", self)
            #self.res_dock_nastran.setObjectName("nastran_model")
            #self.res_dock_nastran.setWidget(self.sidebar_nastran)
            #self.addDockWidget(QtCore.Qt.RightDockWidgetArea, self.res_dock_nastran)

        #self.res_dock.setWidget(self.res_widget)
        if plot:
            gui._finish_results_io2(model_name, [form], cases, reset_labels=reset_labels)
        else:
            gui._set_results([form], cases)
        gui.format = 'astros'

    def make_caeros(self, model: BDF) -> tuple[np.ndarray, int, int, int, int, bool,
                                               dict[int, int], list[int]]:
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
        cs_box_ids : dict[control_surface_name] : list[panel ids]
            list of panels used by each aero panel

        """
        gui: MainWindow = self.gui
        nastran_settings: NastranSettings = gui.settings.nastran_settings
        is_aero = nastran_settings.is_aero

        # when caeros is empty, SPLINEx/AESURF cannot be defined
        if not is_aero or len(model.caeros) == 0:
            caero_points = np.empty((0, 3))
            has_caero = False
            ncaeros = 0
            ncaeros_sub = 0
            ncaeros_cs = 0
            ncaeros_points = 0
            ncaero_sub_points = 0
            has_control_surface = False
            box_id_to_caero_element_map = {}
            cs_box_ids = defaultdict(list)
            all_control_surface_name = ''
            #caero_control_surface_names = []
            out = (
                has_caero, caero_points, ncaeros, ncaeros_sub, ncaeros_cs,
                ncaeros_points, ncaero_sub_points,
                has_control_surface, box_id_to_caero_element_map, cs_box_ids,
            )
            return out

        all_control_surface_name, caero_control_surfaces, out = build_caero_paneling(model)
        if all_control_surface_name:
            gui.create_alternate_vtk_grid(
                'caero_control_surfaces', color=PINK_FLOAT, line_width=5, opacity=1.0,
                representation='surface', is_visible=False)

        for cs_name in caero_control_surfaces:
            gui.create_alternate_vtk_grid(
                cs_name, color=PINK_FLOAT, line_width=5, opacity=0.5,
                representation='surface')
        return out
    
    def map_elements(self, xyz_cid0: np.ndarray, nid_cp_cd: np.ndarray,
                     nid_map: dict[int, int],
                     model: BDF, j: int, dim_max: float,
                     plot: bool=True, xref_loads: bool=True):
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
        gui: MainWindow = self.gui
        grid = gui.grid
        settings: Settings = gui.settings
        nastran_settings: NastranSettings = settings.nastran_settings

        if nastran_settings.is_element_quality:
            out = self._map_elements1_quality(model, xyz_cid0, nid_cp_cd, dim_max, nid_map, j)
        else:
            out = self._map_elements1_no_quality(model, xyz_cid0, nid_cp_cd, dim_max, nid_map, j)
        (nid_to_pid_map, xyz_cid0, superelements, pids, nelements,
         material_coord, material_theta,
         area, min_interior_angle, max_interior_angle, max_aspect_ratio,
         max_skew_angle, taper_ratio, dideal_theta,
         area_ratio, min_edge_length, max_warp_angle) = out

        #self.grid_mapper.SetResolveCoincidentTopologyToPolygonOffset()
        grid.Modified()

        cases = {}

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
            eids = np.zeros(nelements, dtype=nid_cp_cd.dtype)
            eid_map = self.gui.eid_map
            for (eid, eid2) in eid_map.items():
                eids[eid2] = eid

            eid_res = GuiResult(0, header='ElementID', title='ElementID',
                                location='centroid', scalar=eids, mask_value=0)
            cases[icase] = (eid_res, (0, 'ElementID'))
            form0.append(('ElementID', icase, []))
            icase += 1
            self.element_ids = eids
            gui.element_ids = eids

        if superelements is not None:
            nid_res = GuiResult(0, header='SuperelementID', title='SuperelementID',
                                location='centroid', scalar=superelements)
            cases[icase] = (nid_res, (0, 'SuperelementID'))
            form0.append(('SuperelementID', icase, []))
            icase += 1

        # subcase_id, resultType, vector_size, location, dataFormat
        if len(model.properties) and nelements and nastran_settings.is_properties:
            icase, upids, pcomp, pshell, is_pshell_pcomp = self._build_properties(
                model, nelements, eids, pids, cases, form0, icase)
            icase = _build_materials(model, pcomp, pshell, is_pshell_pcomp,
                                     cases, form0, icase)

            try:
                icase = _build_optimization(model, pids, upids, nelements,
                                            cases, form0, icase)
            except Exception:
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
            mean_edge_length = np.nanmean(min_edge_length) * 2.5
            self.gui.set_glyph_scale_factor(mean_edge_length)  # was 1.5

        if (self.make_offset_normals_dim or nastran_settings.is_element_quality) and nelements:
            icase, normals = build_normals_quality(
                settings, model, self.gui.eid_map, nelements, cases, form0, icase,
                xyz_cid0,
                material_coord, material_theta,
                min_interior_angle, max_interior_angle, dideal_theta,
                area, max_skew_angle, taper_ratio,
                max_warp_angle, area_ratio, min_edge_length, max_aspect_ratio,
                make_offset_normals_dim=self.make_offset_normals_dim,
                is_testing=IS_TESTING)
            self.normals = normals
        return nid_to_pid_map, icase, cases, form
    
    def _map_elements1_no_quality(self,
                                  model: BDF,
                                  xyz_cid0: Optional[np.ndarray],
                                  nid_cp_cd: np.ndarray,
                                  unused_dim_max,
                                  nid_map: dict[int, int],
                                  j: int) -> tuple:
        """
        Helper for map_elements

        No element quality
        """
        print('_map_elements1_no_quality')
        assert nid_map is not None

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

        if xyz_cid0 is None:
            superelements = None
            nid_to_pid_map = None
            pids = None
            nelements = None
            material_coord = None
            material_theta = None
        else:
            xyz_cid0 = self.xyz_cid0
            (nid_to_pid_map, xyz_cid0, superelements,
             pids, nelements, material_coord,
             material_theta) = map_elements1_no_quality_helper(
                 self, xyz_cid0, nid_cp_cd, model, nid_map, j)
            self._build_plotels(model)
            self.gui.nelements = nelements
        out = (
            nid_to_pid_map, xyz_cid0, superelements, pids, nelements,
            material_coord, material_theta,
            area, min_interior_angle, max_interior_angle, max_aspect_ratio,
            max_skew_angle, taper_ratio, dideal_theta,
            area_ratio, min_edge_length, max_warp_angle,
        )
        return out
    
    def _map_elements1_quality(self,
                               model: BDF,
                               xyz_cid0: np.ndarray,
                               nid_cp_cd: np.ndarray,
                               unused_dim_max,
                               nid_map: dict[int, int],
                               j: int):
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
            superelements = None
            pids = None
            nelements = None
            material_coord = None
            material_theta = None
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
                nid_to_pid_map, xyz_cid0, superelements, pids, nelements,
                material_coord, material_theta,
                area, min_interior_angle, max_interior_angle, max_aspect_ratio,
                max_skew_angle, taper_ratio, dideal_theta,
                area_ratio, min_edge_length, max_warp_angle,
            )
        else:
            xyz_cid0 = self.xyz_cid0
            out = map_elements1_quality_helper(
                self, model, xyz_cid0, nid_cp_cd, nid_map, j)
        return out
    
    def _create_aero(self, model: BDF,
                     box_id_to_caero_element_map: dict[int, Any],
                     cs_box_ids,
                     caero_points,
                     ncaeros_points: int,
                     ncaero_sub_points: int,
                     has_control_surface: bool):
        # fill grids
        zfighting_offset0 = 0.001
        zfighting_offset = zfighting_offset0
        self._create_splines(model, box_id_to_caero_element_map, caero_points)
        gui: MainWindow = self.gui
        if 'caero' in gui.alt_grids:
            self.set_caero_grid(model)
            self.set_caero_subpanel_grid(model)
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
                    #reset_labels = False
                    cs_name = '%s_control_surface' % aesurf.label
                    self.set_caero_control_surface_grid(
                        cs_name, cs_box_ids[cs_name],
                        box_id_to_caero_element_map, caero_points, note=aesurf.label,
                        zfighting_offset=zfighting_offset)
                    zfighting_offset += zfighting_offset0

    def _create_splines(self, model: BDF,
                        box_id_to_caero_element_map: dict[int, int],
                        caero_points):
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
        if not len(model.splines):
            return

        gui: MainWindow = self.gui
        # 0 - caero / caero_subpanel
        # 1 - control surface
        # 3/5/7/... - spline points
        # 2/4/6/... - spline panels
        iaero = 2
        nsplines = 0
        all_structure_points_list = []
        for spline_id, spline in sorted(model.splines.items()):
            # SPLINE2 -> spline2
            # SPLINE3_ZONEA -> spline3
            spline_type = spline.type.lower().split('_')[0]

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

            all_structure_points_list.extend(structure_points)
            nsplines += 1
            try:
                aero_box_ids = spline.aero_element_ids
            except Exception:
                print(spline.object_attributes())
                print(spline.object_methods())
                raise
            if spline.type != 'SPLINE3_ZAERO':
                assert len(aero_box_ids) > 0, spline
            # the control surfaces all lie perfectly on top of each other
            # such that we have z fighting, so based on the aero index,
            # we calculate a z offset.
            zfighting_offset = 0.0001 * (iaero + 1)
            grid_name = f'{spline_type}_{spline_id:d}_structure_points'
            gui.create_alternate_vtk_grid(
                grid_name, color=BLUE_FLOAT, opacity=1.0, point_size=5,
                representation='point', is_visible=False)
            msg = f', which is required by {grid_name!r}'
            stored_msgi = self._add_astros_nodes_to_grid(
                grid_name, structure_points, model, msg, store_msg=True)

            zfighting_offset = 0.0001 * (iaero + 2)
            grid_name = f'{spline_type}_{spline_id:d}_boxes'
            gui.create_alternate_vtk_grid(
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

        all_structure_points = np.unique(all_structure_points_list)
        if nsplines > 1:
            grid_name = 'all_spline_points'
            gui.create_alternate_vtk_grid(
                grid_name, color=BLUE_FLOAT, opacity=1.0, point_size=5,
                representation='point', is_visible=False)
            msg = f', which is required by {grid_name!r}'
            stored_msgi = self._add_astros_nodes_to_grid(
                grid_name, all_structure_points, model, msg, store_msg=True)
            if stored_msgi:
                stored_msg.append(stored_msgi)

        if stored_msg:
            model.log.warning('\n' + '\n'.join(stored_msg))

    def set_spc_mpc_suport_grid(self, model: BDF,
                                nid_to_pid_map: dict[int, int],
                                idtype: str):
        """
        for each subcase, make secondary actors including:
         - spc_id=spc_id
         - mpc_id=mpc_id              (includes rigid elements)
         - mpc_dependent_id=mpc_id    (includes rigid elements)
         - mpc_independent_id=mpc_id  (includes rigid elements)
         - suport_id=suport1_id       (includes SUPORT/SUPORT1)

        TODO: consider changing the varying ids to huh???
        """
        gui: MainWindow = self.gui
        nastran_settings: NastranSettings = gui.settings.nastran_settings
        spc_names = []
        mpc_names = []
        suport_names = []

        #print('getting rigid')
        rigid_lines = model._get_rigid()

        spc_ids_used = set()
        mpc_ids_used = set()
        suport1_ids_used = set()

        spc_to_subcase = defaultdict(list)
        mpc_to_subcase = defaultdict(list)
        #suport1_to_subcase = defaultdict(list)
        for subcase_id, subcase in sorted(model.subcases.items()):
            if 'SPC' in subcase:
                spc_id = subcase.get_parameter('SPC')[0]
                if spc_id is not None:
                    nspcs = model.card_count['SPC'] if 'SPC' in model.card_count else 0
                    nspc1s = model.card_count['SPC1'] if 'SPC1' in model.card_count else 0
                    nspcds = model.card_count['SPCD'] if 'SPCD' in model.card_count else 0

                    ## TODO: this line seems too loose...
                    ## TODO: why aren't SPCDs included?
                    if nspcs + nspc1s + nspcds:
                        spc_to_subcase[spc_id].append(subcase_id)

            if 'MPC' in subcase:
                mpc_id = subcase.get_parameter('MPC')[0]
                if mpc_id is not None:
                    ## TODO: this line seems too loose
                    nmpcs = model.card_count['MPC'] if 'MPC' in model.card_count else 0
                    if nmpcs:
                        mpc_to_subcase[mpc_id].append(subcase_id)

        for spc_id in chain(model.spcs, model.spcadds):
            spc_name = f'SPC={spc_id:d}'
            if spc_id in mpc_to_subcase:
                subcases = spc_to_subcase[spc_id]
                spc_name += ': Subcases='
                spc_name += ', '.join(str(subcase_id) for subcase_id in subcases)
            spc_names += self._fill_spc(spc_id, spc_name, model, nid_to_pid_map)

        if nastran_settings.is_rbe:
            for mpc_id in chain(model.mpcs, model.mpcadds):
                depname = f'MPC={mpc_id:d}_dependent'
                indname = f'MPC={mpc_id:d}_independent'
                linename = f'MPC={mpc_id:d}_lines'
                if mpc_id in mpc_to_subcase:
                    subcases = mpc_to_subcase[mpc_id]
                    mpc_name = ': Subcases='
                    mpc_name += ', '.join(str(subcase_id) for subcase_id in subcases)
                    depname += mpc_name
                    indname += mpc_name
                    linename += mpc_name

                lines = get_mpc_node_ids(model, mpc_id, stop_on_failure=False)
                lines2 = list(lines)
                mpc_names += self._fill_dependent_independent(
                    mpc_id, model, lines2,
                    depname, indname, linename, idtype)

        if 0:  # pragma: no cover
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
                            spc_name = 'spc_id=%i' % spc_id
                            spc_names += self._fill_spc(spc_id, spc_name, model, nid_to_pid_map)

                # rigid body elements and MPCs
                if 'MPC' in subcase:
                    mpc_id = subcase.get_parameter('MPC')[0]
                    if mpc_id is not None and mpc_id not in mpc_ids_used:
                        mpc_ids_used.add(mpc_id)

                        ## TODO: this line seems too loose
                        nmpcs = model.card_count['MPC'] if 'MPC' in model.card_count else 0
                        if nmpcs:
                            lines = get_mpc_node_ids(model, mpc_id, stop_on_failure=False)
                            lines2 = list(lines)
                            depname = 'mpc_id=%i_dependent' % mpc_id
                            indname = 'mpc_id=%i_independent' % mpc_id
                            linename = 'mpc_id=%i_lines' % mpc_id
                            mpc_names += self._fill_dependent_independent(
                                mpc_id, model, lines2,
                                depname, indname, linename, idtype)

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
            gui.create_alternate_vtk_grid(
                grid_name, color=RED_FLOAT, opacity=1.0, point_size=4,
                representation='point', is_visible=True)

        if nastran_settings.is_rbe and len(rigid_lines):
            # handle RBEs without MPCs
            mpc_id = 0
            depname = 'rigid_dependent'
            indname = 'rigid_independent'
            linename = 'rigid_lines'
            mpc_names += self._fill_dependent_independent(
                mpc_id, model, rigid_lines,
                depname, indname, linename, idtype)

        for grid_name, group in model.model_groups.items():
            if group.nodes is None or not len(group.nodes):
                continue
            nids_colon = []
            for groupi in group.nodes:
                colon_groupi = '%d:%d:%d' % groupi
                expanded_nids = _apply_colon_set(colon_groupi)
                nids_colon.extend(expanded_nids)
            nids = np.unique(np.hstack(nids_colon))
            msg = f', which is required by {grid_name!r}'
            gui.create_alternate_vtk_grid(
                grid_name, color=RED_FLOAT, opacity=1.0, point_size=4,
                representation='point', is_visible=group.is_visible)
            self._add_astros_nodes_to_grid(grid_name, nids, model, msg)
            del nids, nids_colon

        geometry_names = spc_names + mpc_names + suport_names
        return geometry_names

    def _fill_bar_yz(self, unused_dim_max: float,
                     model: BDF, icase: int,
                     cases, form,
                     debug: bool=False) -> int:
        """
        plots the y, z vectors for CBAR & CBEAM elements
        """
        gui: MainWindow = self.gui
        include_rod = False
        card_types = ['CBAR', 'CBEAM']
        if include_rod:
            card_types.append('CROD')
        out = model.get_card_ids_by_card_types(card_types=card_types)
        bar_beam_eids = out['CBAR'] + out['CBEAM']
        if include_rod:
            bar_beam_eids += out['CROD']

        bar_pid_to_eids = get_beam_sections_map(model, bar_beam_eids)
        bar_nids = get_bar_nids(model, bar_beam_eids)
        #ugrid_temp = create_3d_beams(model, bar_pid_to_eids)

        self.bar_eids = {}
        self.bar_lines = {}
        if len(bar_beam_eids) == 0:
            return icase
        scale = 0.15

        # TODO: this should be reworked
        bar_nids, bar_types, nid_release_map = self._get_bar_yz_arrays(
            model, bar_beam_eids, bar_pid_to_eids,
            scale, debug)
        self.nid_release_map = nid_release_map

        bar_nids = list(bar_nids)
        gui.create_alternate_vtk_grid(
            'Bar Nodes', color=RED_FLOAT, line_width=1, opacity=1.,
            point_size=5, representation='point', bar_scale=0., is_visible=False)
        msg = ", which is required by 'Bar Nodes'"
        self._add_astros_nodes_to_grid('Bar Nodes', bar_nids, model, msg)


        geo_form = form[2]
        bar_form = ('CBAR / CBEAM', None, [])
        #print('geo_form =', geo_form)
        #bar_types2 = {}
        bar_eids = []
        for bar_type, data in sorted(bar_types.items()):
            eids, lines_bar_y, lines_bar_z = data
            if len(eids):
                bar_eids.append(eids)
        ibars = 0
        if bar_eids:
            bar_eids = np.hstack(bar_eids)
            ibars = np.searchsorted(self.element_ids, bar_eids)

        for bar_type, data in sorted(bar_types.items()):
            eids, lines_bar_y, lines_bar_z = data
            if len(eids):
                if debug: # pragma: no cover
                    print('bar_type = %r' % bar_type)
                    print('eids     = %r' % eids)
                    print('all_eids = %r' % self.element_ids.tolist())
                # if bar_type not in ['ROD', 'TUBE']:
                bar_y = bar_type + '_y'
                bar_z = bar_type + '_z'

                gui.create_alternate_vtk_grid(
                    bar_y, color=GREEN_FLOAT, line_width=5, opacity=1.,
                    point_size=5, representation='bar', bar_scale=scale, is_visible=False)
                gui.create_alternate_vtk_grid(
                    bar_z, color=BLUE_FLOAT, line_width=5, opacity=1.,
                    point_size=5, representation='bar', bar_scale=scale, is_visible=False)

                self._add_astros_lines_xyz_to_grid(bar_y, lines_bar_y, eids)
                self._add_astros_lines_xyz_to_grid(bar_z, lines_bar_z, eids)

                # form = ['Geometry', None, []]
                i = np.searchsorted(self.element_ids, eids)
                is_type = np.full(self.element_ids.shape, -1, dtype='int32')
                is_type[ibars] = 0
                try:
                    is_type[i] = 1
                except Exception:
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
                                     location='centroid', scalar=is_type, mask_value=-1)
                cases[icase] = (type_res, (0, msg))
                icase += 1

        # print(geo_form)
        if len(bar_form[2]):
            geo_form.append(bar_form)
        return icase

    def _fill_spc(self, spc_id: int, spc_name: str, model: BDF,
                  nid_to_pid_map: dict[int, list[int]]) -> list[str]:
        """creates the spc secondary actors"""
        spc_names = [spc_name]
        gui: MainWindow = self.gui
        gui.create_alternate_vtk_grid(
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
                    if pid < 0:
                        print('pid=%s in _fill_spc...' % pid)
                        continue
                    prop = model.properties[pid]
                    if prop.type not in ['PSOLID', 'PLSOLID']:
                        plot_node = True
                if not plot_node:
                    # don't include 456 constraints if they're ONLY on solid elements
                    # if we had any bar/plate/etc. elements that use this node, we'll plot the node
                    if not('1' in c1 or '2' in c1 or '3' in c1):
                        continue
            node_ids.append(nid)

        node_ids = np.unique(node_ids)
        msg = ', which is required by %r' % spc_name
        self._add_astros_nodes_to_grid(spc_name, node_ids, model, msg)
        return spc_names

    def _set_subcases_unvectorized(self, model: BDF,
                                   form,
                                   cases: dict[int, Any],
                                   icase: int,
                                   xref_nodes: bool,
                                   xref_loads: bool) -> None:
        """helper for ``load_nastran_geometry_unvectorized``"""
        gui: MainWindow = self.gui
        settings: Settings = gui.settings
        colormap = settings.colormap
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
            gui.log_debug('NastranIOv subcase_id = %s' % subcase_id)

            subtitle = ''
            if 'SUBTITLE' in subcase:
                subtitle, options = subcase.get_parameter('SUBTITLE')
                del options

            load_str = f'Load Case={subcase_id:d}' if subtitle == '' else 'Load Case=%i; %s' % (
                subcase_id, subtitle)
            formi = (load_str, None, [])
            formii = formi[2]

            assert icase is not None
            if self.normals is not None and self.plot_applied_loads:
                icase = self._plot_applied_loads(
                    model, cases, formii, icase, subcase_idi, xref_loads=xref_loads,
                    colormap=colormap,
                )
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
                    gui.log_error(sout)
                    print(sout)

            if len(formii):
                form0.append(formi)
        return icase

    def _set_caero_representation(self, has_control_surface: bool) -> None:
        """
        Parameters
        ----------
        has_control_surface : bool
            is there a control surface

        """
        gui: MainWindow = self.gui
        geometry_actors = gui.geometry_actors
        if 'caero_control_surfaces' in geometry_actors:
            gui.geometry_properties['caero_control_surfaces'].opacity = 0.5
            geometry_actors['caero_control_surfaces'].Modified()

        if 'caero' not in geometry_actors:
            return
        geometry_actors['caero'].Modified()
        geometry_actors['caero_subpanels'].Modified()

    def set_caero_control_surface_grid(self, name: str, cs_box_ids: list[int],
                                       box_id_to_caero_element_map: dict[int, Any],
                                       caero_points: np.ndarray,
                                       note: Optional[str]=None,
                                       zfighting_offset: float=0.001,
                                       store_msg: bool=False) -> str:
        """
        Creates a single CAERO control surface?

        Parameters
        ----------
        name : str
            ???
        aero_box_ids : list[int]
            the ids of the box as seen on the AESURF? SET card?
        box_id_to_caero_element_map : dict[key]=value
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
        stored_msg : str
            ???

        """
        gui: MainWindow = self.gui
        log: SimpleLogger = gui.log
        boxes_to_show, stored_msg = check_for_missing_control_surface_boxes(
            name, cs_box_ids, box_id_to_caero_element_map, log,
            store_msg=store_msg)
        #if not boxes_to_show:
            #print('*%s' % name)
            #print('*%s' % boxes_to_show)
            #return

        #if name not in gui.alt_grids:
            #print('**%s' % name)
            #return

        grid = gui.alt_grids[name]
        grid.Reset()

        all_points, elements, centroids, areas = get_caero_control_surface_grid(
            grid,
            box_id_to_caero_element_map,
            caero_points, boxes_to_show, log)

        if len(all_points) == 0:
            log.error('deleting %r' % name)

            # name = spline_1000_boxes
            sname = name.split('_')
            sname[-1] = 'structure_points'

            # points_name = spline_1000_structure_points
            points_name = '_'.join(sname)
            log.error(f'deleting {points_name!r}')

            gui.remove_alt_grid(name, remove_geometry_property=True)
            gui.remove_alt_grid(points_name, remove_geometry_property=True)
            return stored_msg

        # combine all the points
        all_points_array = np.vstack(all_points)

        #vtk_etype = 9 # vtkQuad
        #create_vtk_cells_of_constant_element_type(grid, elements, vtk_etype)

        # shift z to remove z-fighting with caero in surface representation
        all_points_array[:, [1, 2]] += zfighting_offset

        # get the vtk object
        vtk_points = numpy_to_vtk_points(all_points_array, deep=0)
        grid.SetPoints(vtk_points)

        #if missing_boxes:
            #msg = 'Missing CAERO AELIST boxes: ' + str(missing_boxes)
            #gui.log_error(msg)
        if note:
            # points_list (15, 4, 3) = (elements, nodes, 3)
            x, y, z = np.average(centroids, weights=areas, axis=0)
            text = str(note)
            #slot = gui.label_actors[-1]

            slot = gui.reset_label_actors(name)
            annotation = gui.create_annotation(text, x, y, z)
            slot.append(annotation)
        return stored_msg

    def set_caero_grid(self, model: BDF) -> None:
        """
        Sets the CAERO panel geometry.

        Parameters
        ----------
        model : BDF()
            the bdf model

        """
        gui: MainWindow = self.gui
        xyzs = []
        max_cpoints = []
        min_cpoints = []

        zfighting_offset = 0.0001
        caero_grid = gui.alt_grids['caero']
        j = 0
        for unused_eid, element in sorted(model.caeros.items()):
            if isinstance(element, (CAERO1, CAERO3, CAERO4, CAERO5, CAERO7)):
                # wing panel
                cpoints = element.get_points()
                cpoints[0][2] += zfighting_offset
                cpoints[1][2] += zfighting_offset
                max_cpoints.append(np.array(cpoints).max(axis=0))
                min_cpoints.append(np.array(cpoints).min(axis=0))

                elem = vtkQuad()
                point_ids = elem.GetPointIds()
                point_ids.SetId(0, j)
                point_ids.SetId(1, j + 1)
                point_ids.SetId(2, j + 2)
                point_ids.SetId(3, j + 3)
                xyzs.append(cpoints[0])
                xyzs.append(cpoints[1])
                xyzs.append(cpoints[2])
                xyzs.append(cpoints[3])
                caero_grid.InsertNextCell(elem.GetCellType(), point_ids)
                j += 4
            elif isinstance(element, (CAERO2, BODY7)):
                # slender body

                #if 0:  # pragma: no cover
                # 1D version
                #cpoints = element.get_points()
                #cpoints[:, 2] += zfighting_offset
                #max_cpoints.append(np.array(cpoints).max(axis=0))
                #min_cpoints.append(np.array(cpoints).min(axis=0))

                #elem = vtkLine()
                #point_ids = elem.GetPointIds()
                #point_ids.SetId(0, j)
                #point_ids.SetId(1, j + 1)
                #xyzs.append(cpoints[0])
                #xyzs.append(cpoints[1])
                #j += 2
                #caero_grid.InsertNextCell(elem.GetCellType(), point_ids)

                #else:
                # 3D version
                xyz, elems = element.get_points_elements_3d()
                assert xyz is not None, element
                xyz[:, 2] += zfighting_offset
                for elemi in elems:
                    elem = vtkQuad()
                    point_ids = elem.GetPointIds()
                    point_ids.SetId(0, j)
                    point_ids.SetId(1, j + 1)
                    point_ids.SetId(2, j + 2)
                    point_ids.SetId(3, j + 3)
                    n1, n2, n3, n4 = elemi

                    xyzs.append(xyz[n1, :])
                    xyzs.append(xyz[n2, :])
                    xyzs.append(xyz[n3, :])
                    xyzs.append(xyz[n4, :])

                    #cpoints = element.get_points()
                    #cpoints[0][2] += zfighting_offset
                    #cpoints[1][2] += zfighting_offset
                    #max_cpoints.append(np.array(cpoints).max(axis=0))
                    #min_cpoints.append(np.array(cpoints).min(axis=0))
                    caero_grid.InsertNextCell(elem.GetCellType(), point_ids)
                    j += 4
            else:
                gui.log_info(f'skipping {element.type}')

        if len(xyzs) and len(max_cpoints):
            gui.log_info('CAERO.max = %s' % np.vstack(max_cpoints).max(axis=0))
            gui.log_info('CAERO.min = %s' % np.vstack(min_cpoints).min(axis=0))

        xyz = np.array(xyzs)
        xyz = gui.scale_length(xyz)
        vtk_points = numpy_to_vtk_points(xyz, points=None, dtype='<f', deep=1)
        caero_grid.SetPoints(vtk_points)
        #gui.alt_grids['caero']
        #edge_mapper.SetResolveCoincidentTopologyToPolygonOffset()

    def set_caero_subpanel_grid(self, model: BDF) -> None:
        """
        Sets the CAERO sub-panel geometry.

        Parameters
        ----------
        model : BDF()
            the bdf model

        """
        nodes, elements = get_caero_subpanel_grid(model)
        if elements.shape[0] == 0:  # pragma: no cover
            return
        gui: MainWindow = self.gui
        grid = gui.alt_grids['caero_subpanels']
        quad_etype = 9
        create_vtk_cells_of_constant_element_type(grid, elements, quad_etype)

        nodes = gui.scale_length(nodes)
        vtk_points = numpy_to_vtk_points(nodes, points=None, dtype='<f', deep=1)
        grid.SetPoints(vtk_points)
        return

    def get_xyz_in_coord(self, model: BDF,
                         cid: int=0,
                         fdtype: str='float32',
                         check_mirror: bool=True) -> tuple[np.ndarray, np.ndarray]:
        """
        Creates the grid points efficiently

        Used by ``load_nastran_geometry_unvectorized``
        """
        xyz_cid0, nid_cp_cd, icd_transform = build_superelement_model(
            model, cid=cid, fdtype=fdtype)

        if len(xyz_cid0) == 1:
            super_id = 0
            nid_mapi = self.gui.nid_map
            make_nid_map(nid_mapi, nid_cp_cd[super_id][:, 0])
            self._add_astros_spoints_to_grid(model.spoints, nid_mapi)

            self.icd_transform = icd_transform[super_id]
            return xyz_cid0[super_id], nid_cp_cd[super_id]

        # superelements
        self.icd_transform = icd_transform
        xyz_cid0_full = []
        nid_cp_cd_full = []
        for super_id, xyz_cid0i in sorted(xyz_cid0.items()):
            xyz_cid0_full.append(xyz_cid0[super_id])
            nid_cp_cd_full.append(nid_cp_cd[super_id])

        xyz_cid0_out = np.vstack(xyz_cid0_full)
        nid_cp_cd_out = np.vstack(nid_cp_cd_full)

        all_nids = nid_cp_cd_out[:, 0]
        unids = np.unique(all_nids)

        log = self.log
        if not len(all_nids) == len(unids):
            if model.sebulk and check_mirror:
                _prepare_superelement_model(model, log)
                return self.get_xyz_in_coord(model, cid=0, fdtype=fdtype, check_mirror=False)

            msg = ('superelement nodes are not unique; use superelement_renumber\n'
                   'renumbering; duplicate nids=\n%s' % duplicates(all_nids))
            raise NotImplementedError(msg)

        if not is_monotonic(all_nids):
            #msg = ('superelement nodes are not monotonic; use superelement_renumber\n'
                   #'renumbering; nids=\n%s' % all_nids)
            #self.log.warning(msg)
            isort = np.argsort(all_nids)
            xyz_cid0_out = xyz_cid0_out[isort, :]
            nid_cp_cd_out = nid_cp_cd_out[isort, :]

        make_nid_map(self.gui.nid_map, nid_cp_cd_out[:, 0])
        return xyz_cid0_out, nid_cp_cd_out

    def _add_astros_nodes_to_grid(self, name: str,
                                   node_ids: list[int],
                                   model: BDF,
                                   msg: str, store_msg: bool=False):
        """used to create MPC independent/dependent nodes"""
        nnodes = len(node_ids)
        stored_msg = []
        if nnodes == 0:
            msg = f'0 nodes added for {name!r}'
            out_msg = store_warning(model.log, store_msg, msg)
            return out_msg
        gui: MainWindow = self.gui
        gui.follower_nodes[name] = node_ids

        #numpy_to_vtk_points(nodes)
        points = vtkPoints()
        #print(name, 'nnodes', nnodes)
        points.SetNumberOfPoints(nnodes)

        j = 0
        nid_map = gui.nid_map
        alt_grid = gui.alt_grids[name]
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
            elem = vtkVertex()
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, j)
            #else:
                #elem = vtkSphere()
                #dim_max = 1.0
                #sphere_size = self._get_sphere_size(dim_max)
                #elem.SetRadius(sphere_size)
                #elem.SetCenter(points.GetPoint(j))

            alt_grid.InsertNextCell(elem.GetCellType(), point_ids)
            j += 1

        out_msg = ''
        if missing_nodes:
            stored_msg = 'nids=[%s] do not exist%s' % (', '.join(missing_nodes), msg)
        alt_grid.SetPoints(points)

        if stored_msg:
            out_msg = store_warning(model.log, store_msg, stored_msg)
        return out_msg

    def _add_astros_spoints_to_grid(self, spoints: list[int],
                                     nid_map: dict[int, int]) -> None:
        """used to create SPOINTs"""
        if not spoints:  # pragma: no cover
            return
        spoint_ids = list(spoints.keys())
        assert isinstance(spoint_ids, list), type(spoint_ids)

        nspoints = len(spoint_ids)
        name = 'SPoints'
        if nspoints == 0:
            self.log.warning(f'0 spoints added for {name!r}')
            return
        gui: MainWindow = self.gui
        gui.create_alternate_vtk_grid(
            name, color=BLUE_FLOAT, line_width=1, opacity=1.,
            point_size=5, representation='point', bar_scale=0., is_visible=True)

        gui.follower_nodes[name] = spoint_ids
        points = vtkPoints()
        points.SetNumberOfPoints(nspoints)

        j = 0
        alt_grid = gui.alt_grids[name]
        for spointi in sorted(spoint_ids):
            try:
                unused_i = nid_map[spointi]
            except KeyError:
                self.log.warning('spointi=%s does not exist' % spointi)
                continue

            if spointi not in spoints:
                self.log.warning('spointi=%s doesnt exist' % spointi)
                continue
            # point = self.grid.GetPoint(i)
            # points.InsertPoint(j, *point)

            points.InsertPoint(j, 0., 0., 0.)

            elem = vtkVertex()
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, j)
            alt_grid.InsertNextCell(elem.GetCellType(), point_ids)
            j += 1
        alt_grid.SetPoints(points)

    def _add_astros_lines_to_grid(self, name: str,
                                   lines: np.ndarray,
                                   model: BDF,
                                   nid_to_pid_map=None):
        """used to create MPC lines"""
        nlines = lines.shape[0]
        #nids = np.unique(lines)
        #nnodes = len(nids)
        nnodes = nlines * 2
        if nnodes == 0:  # pragma: no cover
            return
        gui: MainWindow = self.gui
        gui.follower_nodes[name] = lines.ravel()
        points = vtkPoints()
        points.SetNumberOfPoints(nnodes)

        j = 0
        etype = 3 # vtkLine
        nid_map = gui.nid_map
        alt_grid = gui.alt_grids[name]
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

            if nid1 not in model.nodes or nid2 not in model.nodes:
                continue
            node = model.nodes[nid1]
            point = node.get_position()
            points.InsertPoint(j, *point)

            node = model.nodes[nid2]
            point = node.get_position()
            points.InsertPoint(j + 1, *point)

            elem = vtkLine()
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, j)
            point_ids.SetId(1, j + 1)
            alt_grid.InsertNextCell(etype, point_ids)
            j += 2
        alt_grid.SetPoints(points)

    def _add_astros_lines_xyz_to_grid(self, name: str, lines, eids) -> None:
        """creates the bar orientation vector lines"""
        nlines = len(lines)
        nnodes = nlines * 2
        if nlines == 0:  # pragma: no cover
            return

        assert name != 'Bar Nodes', name
        gui: MainWindow = self.gui
        grid = gui.alt_grids[name]

        bar_eids = np.asarray(eids, dtype='int32')
        bar_lines = np.asarray(lines, dtype='float32').reshape(nlines, 6)
        self.bar_eids[name] = bar_eids
        self.bar_lines[name] = bar_lines

        nodes = bar_lines.reshape(nlines * 2, 3)
        points = numpy_to_vtk_points(nodes)
        elements = np.arange(0, nnodes, dtype='int32').reshape(nlines, 2)

        etype = 3 # vtkLine().GetCellType()
        create_vtk_cells_of_constant_element_type(grid, elements, etype)
        grid.SetPoints(points)

    def _points_to_vtkpoints_coords(self, model: BDF, xyz_cid0: np.ndarray) -> float:
        """
        helper method for:
         - load_astros_geometry_unvectorized
        """
        gui: MainWindow = self.gui
        points = numpy_to_vtk_points(xyz_cid0)
        gui.grid.SetPoints(points)

        self.xyz_cid0 = xyz_cid0

        maxi = xyz_cid0.max(axis=0)
        mini = xyz_cid0.min(axis=0)
        assert len(maxi) == 3, len(maxi)
        xmax, ymax, zmax = maxi
        xmin, ymin, zmin = mini
        dim_max = max(xmax-xmin, ymax-ymin, zmax-zmin)

        self._create_astros_coords(model, dim_max)

        gui.log_info("xmin=%s xmax=%s dx=%s" % (xmin, xmax, xmax-xmin))
        gui.log_info("ymin=%s ymax=%s dy=%s" % (ymin, ymax, ymax-ymin))
        gui.log_info("zmin=%s zmax=%s dz=%s" % (zmin, zmax, zmax-zmin))
        return dim_max

    def _create_astros_coords(self, model: BDF, dim_max: float) -> None:
        """
        Creates the astros coordinate systems.

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
        gui: MainWindow = self.gui
        gui.create_global_axes(dim_max)
        if not gui.settings.nastran_settings.create_coords:  # pragma: no cover
            return
        for cid, coord in sorted(model.coords.items()):
            if cid in [0, -1]:
                continue
            cid_type = cid_types[coord.Type]
            gui._create_coord(dim_max, cid, coord, cid_type)

    def _build_plotels(self, model: BDF) -> None:
        """creates the plotel actor"""
        nplotels = len(model.plotels)
        if not nplotels:  # pragma: no cover
            return

        gui: MainWindow = self.gui
        nastran_settings: NastranSettings = gui.settings.nastran_settings
        if not nastran_settings.is_plotel:
            return

        lines, tris, quads = plotels_to_groups(model)
        color = nastran_settings.plotel_color

        elements_list = []
        etypes_list = []
        if len(lines):
            etypes_list.append('line')
            elements_list.append(np.array(lines))
        if len(tris):
            etypes_list.append('tri3')
            elements_list.append(np.array(tris))
        if len(quads):
            etypes_list.append('quad4')
            elements_list.append(np.array(quads))

        if len(lines) and (len(tris) or len(quads)):
            representation = 'wire+surf'
        elif len(lines):
            representation = 'wire'
        else:
            representation = 'surface'

        name = 'plotel'
        gui.create_alternate_vtk_grid(
            name, color=color, line_width=2, opacity=0.8,
            point_size=5, representation=representation,
            is_visible=False)
        alt_grid: vtkUnstructuredGrid = gui.alt_grids[name]
        vtk_points: vtkPoints = alt_grid.GetPoints()

        stacked_nids = np.hstack([
            nids.ravel() for nids in elements_list])
        unids = np.unique(stacked_nids)

        nnodes = len(unids)
        xyz_cid0 = np.zeros((nnodes, 3))
        for inid, nid in enumerate(unids):
            node = model.nodes[nid]
            xyz_cid0[inid, :] = node.get_position()

        elements_list2 = [np.searchsorted(unids, nids)
                          for nids in elements_list]
        create_vtk_cells_of_constant_element_types(
            alt_grid, elements_list2, etypes_list)
        #gui.follower_nodes[name] = unids

        xyz_cid0 = gui.scale_length(xyz_cid0)
        vtk_points = numpy_to_vtk_points(xyz_cid0, points=vtk_points)
        alt_grid.SetPoints(vtk_points)

    def _build_properties(self, model: BDF, nelements: int,
                          eids: np.ndarray, pids: np.ndarray,
                          cases, form0,
                          icase: int) -> tuple[int, np.ndarray,
                                               dict[str, Any],
                                               dict[str, Any],
                                               tuple[bool, bool]]:
        """
        creates:
          - PropertyID

        Returns
        -------
        icase : int
            ???
        upids : np.ndarray
            the property ids
        pcomp : dict[str, Any]
        pshell : dict[str, Any]
        flags : tuple[is_pshell, is_pcomp]
            do these exist

        TODO: CONROD
        """
        settings: Settings = self.gui.settings
        nastran_settings: NastranSettings = settings.nastran_settings
        upids = None
        pcomp = None
        pshell = None
        is_pcomp = False
        is_pshell = False

        mids_pcomp = None
        thickness_pcomp = None
        theta_pcomp = None
        nplies_pcomp = None
        pcomp = {
            'mids' : mids_pcomp,
            'thickness' : thickness_pcomp,
            'theta': theta_pcomp,
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

        prop_types_with_mid = (
            'PSOLID',
            'PROD', 'PTUBE', 'PBAR', 'PBARL', 'PBEAM', 'PBEAML',
            'PBEND',
        )
        prop_types_without_mid = ('PVISC', 'PELAS', 'PBUSH', 'PDAMP', 'PDAMPT')

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

        composite_properties = ['PCOMP', 'PCOMPG', 'PCOMPS', 'PCOMPLS']
        pids_pcomp = model.get_card_ids_by_card_types(composite_properties, combine=True)
        properties = model.properties
        for superelement in model.superelement_models.values():
            properties.update(superelement.properties)

        if pids_pcomp:
            npliesi = get_pcomp_nplies(properties, pids_pcomp)

            nplies_pcomp = np.zeros(nelements, dtype='int32')
            thickness_pcomp = np.full((nelements, npliesi), np.nan, dtype='float32')
            theta_pcomp = np.full((nelements, npliesi), np.nan, dtype='float32')
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

            try:
                prop = properties[pid]
            except KeyError:
                print('skipping pid=%i' % pid)
                continue

            if prop.type in prop_types_with_mid:
                # simple types
                i = np.where(pids == pid)[0]
                mid = prop.mid_ref.mid
                mids[i, 0] = mid
            elif prop.type == 'PSHEAR':
                i = np.where(pids == pid)[0]
                mid = prop.mid_ref.mid
                mids[i, 0] = mid
                thickness[i, 0] = prop.Thickness()
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

            elif prop.type in ['PCOMP', 'PCOMPG', 'PCOMPS', 'PCOMPLS']:
                i = np.where(pids == pid)[0]
                npliesi = prop.nplies
                nplies_pcomp[i] = npliesi
                total_pcomp_thickness = 0.
                for iply in range(npliesi):
                    thickniess_ply = prop.Thickness(iply)
                    total_pcomp_thickness += thickniess_ply

                    mids_pcomp[i, iply+1] = prop.Mid(iply)
                    theta_pcomp[i, iply+1] = prop.Theta(iply)
                    thickness_pcomp[i, iply+1] = thickniess_ply
                thickness_pcomp[i, 0] = total_pcomp_thickness
                del total_pcomp_thickness
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
            'theta': theta_pcomp,
            'nplies' : nplies_pcomp,
        }
        pshell = {
            'mids' : mids,
            'thickness' : thickness,
        }
        nplies = None
        if is_pshell:
            nplies = 1
        if is_pcomp:
            nplies = nplies_pcomp.max()

        if nastran_settings.is_shell_mcids and nplies is not None:
            try:
                self._build_mcid_vectors(model, nplies)
            except Exception as exception:
                model.log.error(str(exception))
                if self.stop_on_failure:
                    raise

        return icase, upids, pcomp, pshell, (is_pshell, is_pcomp)

    def _fill_dependent_independent(self, unused_mpc_id: int,
                                    model: BDF, lines,
                                    depname: str, indname: str, linename: str,
                                    idtype: str) -> list[str]:
        """creates the mpc actors"""
        if not lines:
            return []
        gui: MainWindow = self.gui
        settings: NastranSettings = gui.settings.nastran_settings
        gui.create_alternate_vtk_grid(
            depname, color=RED_FLOAT, line_width=5, opacity=1.,
            point_size=5, representation='point', is_visible=False)
        gui.create_alternate_vtk_grid(
            indname, color=BLUE_FLOAT, line_width=5, opacity=1.,
            point_size=5, representation='point', is_visible=False)
        gui.create_alternate_vtk_grid(
            linename, color=settings.rbe_line_color, line_width=5, opacity=1.,
            point_size=5, representation='wire', is_visible=False)

        lines2 = []
        for line in lines:
            if line not in lines2:
                lines2.append(line)
        lines = np.array(lines2, dtype=idtype)
        dependent = lines[:, 0]
        independent = np.unique(lines[:, 1])
        self.dependents_nodes.update(dependent)
        unused_node_ids = np.unique(lines.ravel())

        msg = ', which is required by %r' % depname
        self._add_astros_nodes_to_grid(depname, dependent, model, msg)

        msg = ', which is required by %r' % indname
        self._add_astros_nodes_to_grid(indname, independent, model, msg)

        msg = ', which is required by %r' % linename
        self._add_astros_lines_to_grid(linename, lines, model)

        mpc_names = [depname, indname, linename]
        return mpc_names

    def clear_astros(self):
        """cleans up variables specific to Astros"""
        self.eid_map = {}
        self.nid_map = {}
        self.eid_to_nid_map = {}
        self.element_ids = None
        self.node_ids = None
        self.model = None
        self.nid_release_map = {}
        self.normals = None
        self.icd_transform = {}

        self.dependents_nodes: set[int] = set([])
        #unused_node_ids = np.zeros((0, 2), dtype='int32')
        self.nid_release_map = {}
        self.xyz_cid0 = np.zeros((0, 3), dtype='float64')
        self.nnodes = 0
        self.nelements = 0
        self.bar_eids = {}
        self.bar_lines = {}
        self.gui.isubcase_name_map = {}

    def load_astros_results(self, astros_filename: str):
        """does nothing"""
        raise NotImplementedError()
    
    def _fill_astros_case(self, cases, idi: int, node_ids, nodes, nelements: int,
                          unused_model: Astros) -> tuple[Any, Any, int, np.ndarray, np.ndarray]:
        """creates the result objects for astros"""
        #return [], {}, 0
        #nelements = elements.shape[0]
        #nnodes = nodes.shape[0]

        element_ids = np.arange(1, nelements + 1)
        #print(nodes)
        #node_ids = np.arange(1, nnodes + 1)
        #cnormals = model.get_normals(shift_nodes=False)
        #cnnodes = cnormals.shape[0]
        #assert cnnodes == nelements, len(cnnodes)

        #print('nnodes =', nnodes)
        #print('nelements =', nelements)
        #print('regions.shape =', regions.shape)
        #subcase_id = 0
        #labels = ['NodeID', 'ElementID']
        #cart3d_geo = Cart3dGeometry(subcase_id, labels,
                                    #nids, eids, regions, cnormals,
                                    #uname='Cart3dGeometry')
        colormap = 'jet'
        nid_res = GuiResult(idi, header='NodeID', title='NodeID',
                            location='node', scalar=node_ids)
        eid_res = GuiResult(idi, header='ElementID', title='ElementID',
                            location='centroid', scalar=element_ids)
        nxyz_res = NormalResult(0, 'Normals', 'Normals',
                                nlabels=2, labelsize=5, ncolors=2,
                                colormap=colormap, data_format='%.1f',
                                uname='NormalResult')

        cases[0] = (nid_res, (0, 'NodeID'))
        cases[1] = (eid_res, (0, 'ElementID'))
        cases[2] = (nxyz_res, (0, 'Normal'))

        geometry_form = [
            ('NodeID', 0, []),
            ('ElementID', 1, []),
            ('Normal', 2, []),
        ]
        form = [
            ('Geometry', None, geometry_form),
        ]
        icase = 2
        return form, cases, icase, node_ids, element_ids