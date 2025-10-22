# encoding: utf-8
from __future__ import annotations
#import os
#from itertools import count
import warnings
from pathlib import PurePath
from typing import Optional, Any, TYPE_CHECKING

import numpy as np
#import vtkmodules

from pyNastran.utils import PathLike
from pyNastran.gui.vtk_interface import vtkUnstructuredGrid
#from pyNastran.gui.utils.vtk.base_utils import numpy_to_vtk, numpy_to_vtkIdTypeArray
from pyNastran.gui.gui_objects.gui_result import GuiResult, NormalResult
from pyNastran.gui.gui_objects.types import Formi, Form, Cases #, FormDict, HeaderDict, Case, Cases, KeysMap
from pyNastran.gui.gui_objects.displacements import (
    DisplacementResults, ForceTableResults, ElementalTableResults)

from pyNastran.gui.utils.vtk.vtk_utils import (
    #create_vtk_cells_of_constant_element_type,
    numpy_to_vtk_points)
#from pyNastran.gui.qt_files.colors import (
    #RED_FLOAT, BLUE_FLOAT,
    #LIGHT_GREEN_FLOAT, # GREEN_FLOAT, PINK_FLOAT, PURPLE_FLOAT,
    #YELLOW_FLOAT, ORANGE_FLOAT,
#)
from pyNastran.op2.result_objects.op2_objects import (
    real_modes_to_omega_freq, complex_damping_frequency)

from pyNastran.gui.utils.vtk.vectorized_geometry import (
    create_offset_arrays)
from pyNastran.gui.gui_objects.gui_result import GuiResult# , NormalResult
#from pyNastran.gui.gui_objects.displacements import ForceTableResults, ElementalTableResults
from pyNastran.dev.bdf_vectorized3.bdf import BDF, Subcase
from .alt_actor_builder import (
    build_vtk_geometry,
    create_alt_conm2_grids, create_alt_rbe2_grids, create_alt_rbe3_grids,
    create_alt_spcs, create_alt_axes,
    create_monpnt, create_plotels)
from pyNastran.gui.gui_objects.displacement_results import DisplacementResults2
from pyNastran.gui.gui_objects.force_results import ForceResults2
from pyNastran.op2.result_objects.stress_object import _get_nastran_header
from pyNastran.converters.nastran.gui.types import KeysMap, KeyMap, NastranKey

IS_TABLES = True
if IS_TABLES:
    from pyNastran.dev.op2_vectorized3.op2_hdf5 import OP2, OP2Geom
    from pyNastran.dev.op2_vectorized3.op2_hdf5 import Results
if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger
    from pyNastran.gui.gui_objects.settings import Settings, NastranSettings
    from pyNastran.gui.main_window import MainWindow
    from pyNastran.dev.bdf_vectorized3.bdf_interface.bdf_attributes import (
        CTETRA, CPENTA, CHEXA, CPYRAM,
    )
    SolidElement = CTETRA | CPENTA | CHEXA | CPYRAM


class Nastran3:
    """doesn't have object bleedover like NastranIO"""
    def __init__(self, gui: MainWindow):
        self.gui = gui
        self.data_map = None
        self.save_results_model = False
        self.stop_on_failure = False

        # filters
        self.etypes_filter = {} # TODO: doesn't work
        # TODO: only affects elements
        self.ifile_filter = np.array([], dtype='int32')
        self.properties_filter = np.array([], dtype='int32')

        # TODO: does nothing
        # TODO: need to consider CONROD as well as traced materials
        # TODO: PSHELLs/PCOMPs have multiple materials...how does FEMAP do it?
        # TODO:  - all vs. any filter
        self.materials_filter = np.array([], dtype='int32')

        self.model = BDF(debug=True, log=None, mode='msc')
        #self.model.is_strict_card_parser = False

        # the set of element types that are supported
        self.gui_elements: set[str] = set([])
        self.bar_eids = {}
        self.bar_lines = {}
        self.include_mass_in_geometry = True

        # dictionary of filename to boolean (???)
        self.card_index = {}

    def get_nastran3_wildcard_geometry_results_functions(self) -> tuple:
        """
        gets the Nastran wildcard loader used in the file load menu
        """
        data = (
            'Nastran3',
            'Nastran3 (*.bdf; *.dat; *.ecd, *.h5, *.op2)', self.load_nastran3_geometry,
            #'NastranV (*.h5)', self.load_h5_results,
            'Nastran3 (*.op2)', self.load_op2_results,
        )
        return data

    def load_op2_results(self, op2_filename: PathLike,
                         plot: bool=True) -> None:
        """loads results from an op2 file"""
        model = self.get_op2_results(op2_filename, log=self.gui.log)
        self._load_op2_results(model, plot)

    def get_op2_results(self, op2_filename: PathLike | OP2, log=None) -> OP2:
        from pyNastran.op2.op2 import OP2
        assert isinstance(op2_filename, (str, PurePath)), op2_filename
        model = OP2(debug=True, log=self.gui.log, mode='msc')
        model.read_matpool = False
        model.read_op2(op2_filename, combine=True, build_dataframe=False,
                       skip_undefined_matrices=False, encoding=None)
        return model

    def load_h5_results(self, h5_filename: PathLike, plot: bool=True) -> None:
        """loads results from an h5 file"""
        model = self.get_h5_results(h5_filename, log=self.gui.log)
        self._load_op2_results(model, plot)

    def get_h5_results(self, h5_filename: PathLike | Results,
                       log=None) -> Results:
        assert isinstance(h5_filename, (str, PurePath)), h5_filename
        model = Results(debug=True, log=log, mode='msc')
        model.read_matpool = False
        model.read_h5(h5_filename)
        return model

    def _load_op2_results(self, model: OP2, plot: bool) -> None:
        """loads results from a filled OP2 object (for op2/h5)"""
        name = 'main'
        gui: MainWindow = self.gui
        cases: Cases = gui.result_cases
        form: Form = gui.get_form()
        icase = len(cases)
        icase = self._load_results_from_model(
            model,
            form, cases, icase,
            name, plot)
        gui._finish_results_io2(name, form, cases)

    def load_nastran3_geometry(self, bdf_filename: PathLike,
                               name: str='main', plot: bool=True):
        """loads geometry from a bdf/op2/h5 file"""
        log = self.gui.log
        if isinstance(bdf_filename, PathLike):
            bdf_filename_lower = bdf_filename.lower()
            if bdf_filename_lower.endswith('.op2'):
                geo = self.load_op2_geometry(bdf_filename)
            elif bdf_filename_lower.endswith('.h5'):
                geo = self.load_h5_geometry(bdf_filename)
            else:  # str
                #log.debug(f'start load_bdf_geometry; type={type(bdf_filename)}')
                geo = self.load_bdf_geometry(bdf_filename)
                #log.error(f'end load_bdf_geometry; type={type(bdf_filename)}')
        elif isinstance(bdf_filename, BDF):
            model = bdf_filename
            geo = self.load_bdf_geometry(model)
        else:
            raise RuntimeError(type(bdf_filename))

        if not hasattr(self, 'xyz_cid0'):
            log.error('cant find xyz_cid0')
        self.xyz_cid0
        #log.info('xyz_cid0 exists...')
        return geo

    def load_nastran3_results(self, results_filename: PathLike,
                              name: str='main', plot: bool=True) -> vtkUnstructuredGrid:
        """loads geometry from a op2/h5 file"""
        #self.gui.log.info('load_nastran3_results')
        self.xyz_cid0
        results_filename = str(results_filename)
        results_filename_lower = results_filename.lower()
        if results_filename_lower.endswith('.op2'):
            return self.load_op2_results(results_filename)
        elif results_filename_lower.endswith('.h5'):
            return self.load_h5_results(results_filename)
        raise RuntimeError(f'results_filename={results_filename!r} is not supported')

    def load_op2_geometry(self, op2_filename: PathLike,
                          name: str='main', plot: bool=True) -> vtkUnstructuredGrid:
        """loads geometry from an op2 file"""
        assert isinstance(op2_filename, (str, PurePath)), op2_filename
        model = OP2Geom(debug=True, log=None, mode='msc')
        idtype = model.idtype
        log = model.log

        if not str(op2_filename).lower().endswith('.op2'):
            log.error(f'op2_filename={op2_filename!r} is not a .op2 file')
            assert op2_filename.endswith('.op2'), op2_filename
        model.read_op2(op2_filename)
        model.idtype = idtype
        model.setup(run_geom_check=False)

        ugrid, form, cases, icase = self._load_geometry_from_model(
            model, name, plot)

        self._load_results_from_model(model,
                                      form, cases, icase,
                                      name, plot)
        self.gui._finish_results_io2(name, form, cases)
        return ugrid

    def _load_results_from_model(self, model: OP2,
                                 form: Form, cases: Cases, icase: int,
                                 name: str, plot: bool) -> int:
        #self.xyz_cid0 = model.grid.xyz_cid0() # .astype('float32')
        #self.node_ids = model.grid.node_id
        log = model.log
        if not hasattr(self, 'node_id'):
            log.error(f'missing node_id; load the geometry')
            return

        gui = self.gui
        settings: Settings = gui.settings
        nastran_settings: NastranSettings = settings.nastran_settings

        xyz_cid0 = self.xyz_cid0
        node_id = self.node_id
        element_id = self.element_id
        mmax = xyz_cid0.max(axis=0)
        mmin = xyz_cid0.min(axis=0)
        dim_max = (mmax - mmin).max()
        #nnodes = len(node_id)

        #model.transform_displacements_to_global(icd_transform, coords, xyz_cid0=None, debug=False)
        # stress
        #result_name = 'PropertyID'
        #prop_res = GuiResult(subcase_id, header=result_name, title=result_name,
                             #location='centroid', scalar=property_id, mask_value=0)

        #results = model.res
        name_results = [
            ('Displacement', model.displacements, 'disp'),
            ('Eigenvector', model.eigenvectors, 'disp'),
            ('Temperature', model.temperatures, 'temp'),
            ('Velocity', model.velocities, 'disp'),
            ('Acceleration', model.accelerations, 'disp'),
            ('SPC Forces', model.spc_forces, 'force'),
            ('MPC Forces', model.mpc_forces, 'force'),
        ]

        subcases = set([])
        for (res_name, results, deflects_str) in name_results:
            for key, case in results.items():
                subcases.add(key)

        # SORT1
        #print(f'subcases = {subcases}')
        for subcase in subcases:
            subcase_form: Form = []
            #print(subcase)
            icase = _load_oug(model, name, name_results,
                              key, subcase, subcase_form,
                              cases, icase,
                              node_id, xyz_cid0, dim_max)
            if nastran_settings.stress:
                icase = _load_stress_strain(
                    self.element_cards,
                    self.model, model, name,
                    key, subcase, subcase_form,
                    cases, icase,
                    element_id, is_stress=True)

            if nastran_settings.strain:
                icase = _load_stress_strain(
                    self.element_cards,
                    self.model, model, name,
                    key, subcase, subcase_form,
                    cases, icase,
                    element_id, is_stress=False)

            if len(subcase_form):
                form.append((f'Subcase {subcase}', None, subcase_form))
                print(form)
            #('Geometry', None, geometry_form),

            #e = ('Subcase 1', None,
                 #[('Displacement', None, [('Static', 13, [])]),
                  #[('Stress', None,
                    #[('oxx', 14, []), ('oyy', 15, []), ('ozz', 16, []),
                     #('txy', 17, []), ('txz', 18, []), ('tyz', 19, []),
                     #('omax', 20, []), ('omin', 21, []), ('Von Mises', 22, [])]
                    #)
                   #]])

        if self.save_results_model:
            self.model = model
        return icase

    def load_h5_geometry(self, h5_filename: PathLike,
                         name: str='main', plot: bool=True) -> vtkUnstructuredGrid:
        """loads a geometry only an h5 file"""
        model = OP2Geom(debug=True, log=None, mode='msc')
        model.read_h5(h5_filename)
        ugrid, form, cases, icase = self._load_geometry_from_model(
            model, name, plot)
        self._load_results_from_model(model,
                                      form, cases, icase,
                                      name, plot)
        self.gui._finish_results_io2(name, form, cases)
        return ugrid

    def load_bdf_geometry(self, bdf_filename: PathLike,
                          name: str='main', plot: bool=True) -> vtkUnstructuredGrid:
        """loads a geometry only an h5 file"""
        log = self.gui.log
        model = self.get_bdf_geometry(bdf_filename, log=log)
        ugrid, form, cases, unused_icase = self._load_geometry_from_model(
            model, name, plot)
        self.gui._finish_results_io2(name, form, cases)
        return ugrid

    def get_bdf_geometry(self, bdf_filename: PathLike | BDF,
                         log=None) -> BDF:
        if isinstance(bdf_filename, BDF):
            model = bdf_filename
        else:
            model = BDF(debug=True, log=log, mode='msc')
            model.is_strict_card_parser = False
            model.allow_overwrites_set = {'GRID', 'CONM2'}
            model.idtype = 'int64'
            model.read_bdf(bdf_filename)

        for prop in model.property_cards + model.aero_property_cards:
            # if prop.n == 0:
            #     continue
            try:
                assert len(prop.ifile) == len(prop.property_id)
            except:
                log.error(f'{prop.type}.ifile error; ifile={prop.ifile}')
                raise
        for mat in model.material_cards:
            # if mat.n == 0:
            #     continue
            try:
                assert len(mat.ifile) == len(mat.material_id)
            except:
                log.error(f'{mat.type}.ifile error; ifile={mat.ifile}')
                raise
        for elem in model.aero_element_cards + model.element_cards:
            # if elem.n == 0:
            #     continue
            try:
                assert len(elem.ifile) == len(elem.element_id), elem.type
            except:
                log.error(f'{elem.type}.ifile error; ifile={elem.ifile}')
                raise
        return model

    def _load_geometry_from_model(self, model: BDF, name: str, plot: bool,
                                  ) -> tuple[vtkUnstructuredGrid, Form, Cases, int]:
        self.model = model
        log = model.log
        gui = self.gui
        gui.eid_maps[name] = {}
        gui.nid_maps[name] = {}

        self.xyz_cid0 = model.grid.xyz_cid0() # .astype('float32')
        self.node_id = model.grid.node_id

        xyz_cid0 = self.xyz_cid0
        node_id = self.node_id

        mmax = xyz_cid0.max(axis=0)
        mmin = xyz_cid0.min(axis=0)
        dim_max = (mmax - mmin).max()
        xmax, ymax, zmax = mmax
        xmin, ymin, zmin = mmin
        gui.log_info("xmin=%s xmax=%s dx=%s" % (xmin, xmax, xmax-xmin))
        gui.log_info("ymin=%s ymax=%s dy=%s" % (ymin, ymax, ymax-ymin))
        gui.log_info("zmin=%s zmax=%s dz=%s" % (zmin, zmax, zmax-zmin))
        gui.create_global_axes(dim_max)

        points = numpy_to_vtk_points(xyz_cid0)
        #points_data = numpy_to_vtk(xyz_cid0, deep=1)
        #npoints = len(nids)

        #points = vtkPoints()
        #points.SetNumberOfPoints(npoints)
        #points.SetData(points_data)

        ugrid: vtkUnstructuredGrid = gui.grid
        create_alt_spcs(gui, model, node_id, xyz_cid0)
        create_alt_conm2_grids(gui, model, node_id, xyz_cid0)
        create_alt_rbe2_grids(gui, model, node_id, xyz_cid0)
        create_alt_rbe3_grids(gui, model, node_id, xyz_cid0)
        create_alt_axes(self, gui, model, node_id, xyz_cid0)
        create_monpnt(gui, model, node_id, xyz_cid0)
        create_plotels(gui, model, node_id, xyz_cid0)

        # add alternate actors
        gui._add_alt_actors(gui.alt_grids)

        out = load_elements(
            self.etypes_filter, self.ifile_filter,
            self.properties_filter, self.materials_filter,
            ugrid, model, node_id,
            self.include_mass_in_geometry)
        etypes, element_id, property_id, element_cards, gui_elements, card_index = out
        self.element_types = etypes
        self.element_id = element_id
        self.element_cards = element_cards
        self.gui_elements = gui_elements
        self.card_index = card_index

        #print(f'nelements = {len(element_id)}')
        form, cases, icase = self.save_results(model, node_id, element_id, property_id)

        ugrid.SetPoints(points)
        ugrid.Modified()

        gui.scalar_bar_actor.VisibilityOn()
        gui.scalar_bar_actor.Modified()
        return ugrid, form, cases, icase

    def save_results(self,
                     model: BDF,
                     node_id: np.ndarray,
                     element_id: np.ndarray,
                     property_id: np.ndarray) -> tuple[Form, Cases, int]:
        log = model.log
        node_index = np.arange(len(node_id))
        element_index = np.arange(len(element_id))

        gui = self.gui
        settings: Settings = gui.settings
        nastran_settings: NastranSettings = settings.nastran_settings

        # I think we specifically look up NodeID, ELementID,
        # but I think we want to look up by Index here
        self.node_id = node_id
        self.element_id = element_id

        gui.node_ids = node_id        # TODO: should this be node_id/index
        gui.element_ids = element_id  # TODO: should this be element_id/index
        gui.nnodes = len(node_id)
        gui.nelements = len(element_id)

        self.gui.log_info(f'nnodes={gui.nnodes:d} nelements={gui.nelements:d}')
        msg = model.get_bdf_stats(return_type='string')
        #print(msg)
        self.gui.log_debug(msg)
        msg = model.get_bdf_stats(return_type='list')


        gui.isubcase_name_map = {1: ['Nastran', '']}
        cases: Cases = {}
        quality_form = []

        icase = self.gui.get_new_icase()

        subcase_id = 1
        geometry_form: list[Form] = []
        nelements = len(element_id)
        icase = _add_integer_node_gui_result(icase, cases, geometry_form, subcase_id, 'NodeID', node_id)
        if nelements:
            icase = _add_integer_centroid_gui_result(icase, cases, geometry_form, subcase_id, 'ElementID', element_id)
        icase = _add_integer_node_gui_result(icase, cases, geometry_form, subcase_id, 'NodeIndex', node_index)
        if nelements:
            icase = _add_integer_centroid_gui_result(icase, cases, geometry_form, subcase_id, 'ElementIndex', element_index)

        if len(model.grid):
            node_cp = model.grid.cp
            node_cd = model.grid.cd
            if (node_cp.min(), node_cp.max()) != (0, 0):
                icase = _add_integer_node_gui_result(icase, cases, geometry_form, subcase_id, 'NodeCp', node_cp)
            if (node_cd.min(), node_cd.max()) != (0, 0):
                icase = _add_integer_node_gui_result(icase, cases, geometry_form, subcase_id, 'NodeCd', node_cd)

        if nelements:
            icase = _add_integer_centroid_gui_result(
                icase, cases, geometry_form, subcase_id, 'PropertyID', property_id, mask_value=0)

            uproperty_id = get_property_index(property_id)
            icase = _add_integer_centroid_gui_result(
                icase, cases, geometry_form, subcase_id, 'PropertyIndex', uproperty_id)

            mean_edge_length, icase, quality_form = _set_quality(
                icase, cases, nastran_settings,
                model, subcase_id, self.gui_elements,
                element_id, nelements, self.gui_elements)
            self.mean_edge_length = mean_edge_length

            try:
                # material_id[***ieid, :ilayer***] = material_id_rect
                #
                # ValueError: shape mismatch: value array of shape (636,28) could not be
                # broadcast to indexing result of shape (636,24)
                icase = gui_material_ids(
                    model, icase, cases, geometry_form,
                    element_id, property_id, self.card_index)
            except Exception as e:
                log.error('failed creating material_ids results')
                if self.stop_on_failure:
                    raise

        form = [
            ('Geometry', None, geometry_form),
        ]
        if len(quality_form):
            geometry_form.append(('Quality', None, quality_form))

        return form, cases, icase# , nids, eids, data_map_dict

    def _simple_gui(self, ugrid: vtkUnstructuredGrid) -> None:  # pragma: no cover
        from pyNastran.gui.vtk_rendering_core import (
            vtkActor, vtkDataSetMapper, vtkRenderer,
            vtkRenderWindow, vtkRenderWindowInteractor)
        grid_mapper = vtkDataSetMapper()
        grid_mapper.SetInputData(ugrid)

        geom_actor = vtkActor()
        geom_actor.SetMapper(grid_mapper)

        # Setup renderer
        renderer = vtkRenderer()
        renderer.AddActor(geom_actor)
        #if make_glyphs:
            #renderer.AddActor(arrow_actor)
        renderer.ResetCamera()
        renderer.SetBackground(0.7, 0.8, 1.0)

        # Setup render window
        renderWindow = vtkRenderWindow()
        renderWindow.AddRenderer(renderer)

        # Setup render window
        render_window = vtkRenderWindow()
        render_window.AddRenderer(renderer)

        # Setup render window interactor
        renderWindowInteractor = vtkRenderWindowInteractor()
        #style = vtkInteractorStyleImage()

        # Render and start interaction
        renderWindowInteractor.SetRenderWindow(render_window)
        renderWindowInteractor.Initialize()

        renderWindowInteractor.Start()
        #x = 1


def load_elements(etypes_filter: set[str],
                  ifile_filter: np.ndarray,
                  properties_filter: np.ndarray,
                  materials_filter: np.ndarray,
                  ugrid: vtkUnstructuredGrid,
                  model: BDF,
                  grid_id: np.ndarray,
                  include_mass_in_geometry: bool,
                  ) -> tuple[np.ndarray, np.ndarray, np.ndarray,
                             list[Any], set[str],
                             dict[str, tuple[int, int]]]:
    """
    Fills the vtkUnstructuredGrid.

    Elements are added by the order in model.element_cards, not by element id.
    This makes it easier to fill the geometry/results, but harder to lookup
    a specific element id.
    """
    log = model.log
    property_ids: list[np.ndarray] = []
    element_ids: list[np.ndarray] = []
    etypes: list[np.ndarray] = []
    #nodes_ids: list[np.ndarray] = []
    #nnodes: list[np.ndarray] = []
    n_nodes_ = []
    cell_type_ = []
    cell_offset_ = []
    etypes_ = []

    cell_type_point = 1  # vtkVertex().GetCellType()
    cell_type_line = 3  # vtkLine().GetCellType()
    cell_type_tri3 = 5
    cell_type_tri6 = 22
    cell_type_quad4 = 9
    cell_type_quad8 = 23

    cell_offset0 = 0

    # elements that don't use SPOINTs or
    basic_elements = {
        'CONROD', 'CTUBE', 'CROD',
        'CBEAM', 'CBAR', 'CBUSH',
        'CBAR', 'CBEAM', 'CSHEAR',
        'CTRIA3', 'CQUAD4', 'CTRIAR', 'CQUADR',
    }
    midside_elements = {
        'CTRIA6', 'CQUAD8', # 'CQUAD'
        'CTRIAX6',
    }
    solid_elements = {'CTETRA', 'CPENTA', 'CHEXA', 'CPYRAM'}

    gui_elements = basic_elements | solid_elements | midside_elements
    cards_to_drop_silenced = {
        'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
        'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CVISC',
        'CONM2',
    }

    nelement0 = 0
    card_index = {}
    element_cards = [card for card in model.element_cards if card.n]
    etypes_filter_list = list(etypes_filter)
    etypes_filter_list.sort()
    if len(etypes_filter_list):
        log.debug('filters:\n'
                  f'  etypes_filter = {etypes_filter_list}'
                  f'  ifile_filter = {ifile_filter}\n'
                  f'  properties_filter = {properties_filter}\n'
                  f'  materials_filter = {materials_filter}\n')
    for element in element_cards:
        etype = element.type
        is_removed = len(etypes_filter) and etype not in etypes_filter
        if is_removed:
            log.warning(
            f'\netype={etype!r} etypes_filter={etypes_filter}\n'
            f'  len(etypes_filter)={len(etypes_filter)}\n'
            f'  etype not in etypes_filter={etype not in etypes_filter}')
            log.warning(f'skipping {etype}')
            continue

        is_removed, element = element.filter_by_ifile_filter(ifile_filter)
        if is_removed:
            log.warning(f'skipping {etype} due to ifile_filter')
            continue
        #is_removed, element = element.filter(
            #ifile_filter, properties_filter, materials_filter)
        is_removed, element = element.filter_by_property_filter(
             properties_filter)
        if is_removed:
            #log.warning(f'skipping {etype} due to ifile/prop/mat filter')
            log.warning(f'skipping {etype} due to property_filter')
            continue
        nelement = element.n
        #log.info(f'working on {etype} with nelement={nelement}')

        if etype in basic_elements:
            #log.debug('  basic')
            # basic elements
            # missing nodes are not allowed
            if etype in {'CONROD'}:
                cell_type = cell_type_line
                property_id = np.full(nelement, -1)
            elif etype in {'CROD', 'CTUBE', 'CBAR', 'CBEAM', 'CBUSH', 'CGAP'}:
                cell_type = cell_type_line
                property_id = element.property_id
            elif etype in {'CTRIA3', 'CTRIAR'}:
                cell_type = cell_type_tri3
                property_id = element.property_id
            elif etype in {'CQUAD4', 'CQUADR'}:
                #log.debug('cell_type_quad4')
                cell_type = cell_type_quad4
                property_id = element.property_id
            elif etype == 'CSHEAR':
                cell_type = cell_type_quad4
                property_id = element.property_id
            else:  # pragma: no cover
                raise NotImplementedError(element)

            element_id = element.element_id
            assert len(element_id) == len(property_id), etype
            nodesi = element.nodes

            if 0:  # pragma: no cover
                pids_to_filter = []
                i = filter_property_ids(etype, property_id, pids_to_filter, log)
                #nelement = len(i)
                element_id = element_id[i]
                nelement = len(element_id)
                property_id = property_id[i]
                nodesi = nodesi[i, :]

            dnode = nodesi.shape[1]
            cell_offset0, n_nodesi, cell_typei, cell_offseti = create_offset_arrays(
                grid_id, nodesi,
                nelement, cell_type, cell_offset0, dnode)

            n_nodes_.append(n_nodesi.ravel())
            cell_type_.append(cell_typei)
            cell_offset_.append(cell_offseti)
            element_ids.append(element_id)
            property_ids.append(property_id)
            #nnodes_nodes = [nnodes, node_id]
            del cell_type
            del cell_offseti
            cell_offset0 += nelement * (dnode + 1)
            #log.debug(f'end of {etype} adding')
        elif etype in midside_elements:
            #log.debug('  midside')
            midside_nodes = element.midside_nodes
            min_midside_nid = midside_nodes.min(axis=1)
            assert len(min_midside_nid) == nelement

            #print('min_midside_nid =', min_midside_nid)
            is_partial = (min_midside_nid == 0)
            is_full = ~is_partial

            element_id = element.element_id
            if etype == 'CTRIAX6':
                property_id = np.full(element_id.shape, 0, dtype='int32')
            else:
                property_id = element.property_id
            nfull = is_full.sum()
            npartial = is_partial.sum()

            if nfull and npartial:
                log.warning(f'{etype} nelement={nelement} downcasting midside elements to first order')

            #log.debug(f'{etype} nelement={nelement} nfull={nfull} npartial={npartial}')
            if nfull:
                if etype == 'CTRIA6':
                    cell_type = cell_type_tri6
                    nodes = element.nodes
                elif etype == 'CTRIAX6':
                    cell_type = cell_type_tri6
                    nodes = element.nodes[:, [0, 2, 4, 1, 3, 5]]
                elif etype == 'CQUAD8':
                    cell_type = cell_type_quad8
                    nodes = element.nodes
                else:  # pragma: no cover
                    raise NotImplementedError(f'full {etype}')
                cell_offset0 = _save_element(
                    is_full, grid_id, cell_type, cell_offset0,
                    element_id, property_id, nodes,
                    element_ids,
                    property_ids,
                    n_nodes_,
                    cell_type_,
                    cell_offset_)

            if npartial:
                #log.debug('  partial')
                if etype == 'CTRIA6':
                    cell_type = cell_type_tri3
                    nodes = element.nodes[:, :3]
                elif etype == 'CTRIAX6':
                    cell_type = cell_type_tri3
                    nodes = element.nodes[:, [0, 2, 4]]
                elif etype == 'CQUAD8':
                    cell_type = cell_type_quad8
                    nodes = element.nodes[:, :4]
                else:  # pragma: no cover
                    raise NotImplementedError(f'full {etype}')
                is_full = np.arange(element.n, dtype='int32')
                cell_offset0 = _save_element(
                    is_full, grid_id, cell_type, cell_offset0,
                    element_id, property_id, nodes,
                    element_ids, property_ids,
                    n_nodes_, cell_type_, cell_offset_)
                etypes_.append(np.full(nelement, etype))
                continue
            nodesi = element.base_nodes

        elif etype in solid_elements:
            #log.debug('  solid')
            cell_offset0, n_nodesi, cell_typei, cell_offseti = _create_solid_vtk_arrays(
                element, grid_id, cell_offset0)

            n_nodes_.append(n_nodesi.ravel())
            element_ids.append(element.element_id)
            property_ids.append(element.property_id)
            cell_type_.append(cell_typei)
            cell_offset_.append(cell_offseti)
            del n_nodesi, cell_typei, cell_offseti

        elif include_mass_in_geometry and etype == 'CONM2' and 0:
            log.debug('  CONM2 include_mass_in_geometry')
            cell_type = cell_type_point
            property_ids.append(np.full(nelement, -1))

            nodesi = element.node_id.reshape(nelement, 1)
            nodes_indexi = np.searchsorted(grid_id, nodesi)
            dnode = 1
            nnodesi = np.ones((nelement, 1), dtype='int32') * dnode
            cell_typei = np.ones(nelement, dtype='int32') * cell_type

            # (nnodes+1) = 4+1 = 5
            # [0, 5, 10, 15, 20, ... (nelements-1)*5]
            #
            # for 2 CQUAD4s elements (4 nodes; 5 columns including the node count of 4)
            # [0, 5]
            # should be length nelement
            cell_offseti = cell_offset0 + np.arange(0, nelement * (dnode + 1), dnode + 1)
            assert len(cell_offseti) == nelement

            #nnodes
            #print('nnodesi =', nnodesi)
            #print('nodes_indexi =', nodes_indexi)
            n_nodesi = np.hstack([nnodesi, nodes_indexi])
            n_nodes_.append(n_nodesi.ravel())

            cell_type_.append(cell_typei)
            cell_offset_.append(cell_offseti)
            element_ids.append(element.element_id)
            #nnodes_nodes = [nnodes, node_id]
            del cell_type
            del cell_offseti, nnodesi, nodesi
            cell_offset0 += nelement * (dnode + 1)
        else:
            if element.type in cards_to_drop_silenced:
                #log.warning(f'  drop silenced {element}')
                continue
            # more complicated element
            log.warning(f'  dropping {element}')
            continue
            #raise NotImplementedError(element.type)
        card_index[etype] = (nelement0, nelement0 + nelement)
        nelement0 += nelement
        etypes_.append(np.full(nelement, etype))
        del nelement # , dnode

        if 0:  # pragma: no cover
            #  testing to make sure things work...useful when they don't :/
            # [number of nodes, nodes]
            n_nodes = np.hstack(n_nodes_)
            element_id = np.hstack(element_ids)
            property_id = np.hstack(property_ids)
            cell_type = np.hstack(cell_type_)  # 10, 12
            cell_offset = np.hstack(cell_offset_) # 0, 5
            etypes = np.hstack(etypes_)

            assert len(element_id) == len(etypes)
            assert len(element_id) == len(property_id)
            assert len(element_id) == len(cell_type)
            assert len(element_id) == len(cell_offset)
            #log.debug('***passed checks')
    # [number of nodes, nodes]
    assert len(element_ids) == len(property_ids)

    if len(element_ids) == 0:
        element_id = np.array([], dtype='int32')
        property_id = np.array([], dtype='int32')
    else:
        n_nodes = np.hstack(n_nodes_)
        element_id = np.hstack(element_ids)
        property_id = np.hstack(property_ids)
        cell_type = np.hstack(cell_type_)  # 10, 12
        cell_offset = np.hstack(cell_offset_) # 0, 5
        etypes = np.hstack(etypes_)
        #property_id = np.hstack(property_ids)

        assert len(element_id) == len(etypes)
        assert len(element_id) == len(property_id)
        assert len(element_id) == len(cell_type)
        assert len(element_id) == len(cell_offset)

        nelement_total = len(element_id)
        #nelement_total = 5 # len(element_id)
        build_vtk_geometry(
            nelement_total, ugrid,
            n_nodes, cell_type, cell_offset)
    return etypes, element_id, property_id, element_cards, gui_elements, card_index

def _save_element(i: np.ndarray,
                  grid_id: np.ndarray,
                  cell_type: int,
                  cell_offset0: int,
                  element_id: np.ndarray,
                  property_id: list[np.ndarray],
                  nodes: list[np.ndarray],
                  element_ids: list[np.ndarray],
                  property_ids: list[np.ndarray],
                  n_nodes_: list[np.ndarray],
                  cell_type_: list[np.ndarray],
                  cell_offset_: list[np.ndarray]) -> int:
    nodesi = nodes[i, :]
    nelement, dnode = nodesi.shape
    element_idi = element_id[i]
    property_idi = property_id[i]

    cell_offset0, n_nodesi, cell_typei, cell_offseti = create_offset_arrays(
        grid_id, nodes,
        nelement, cell_type, cell_offset0, dnode)

    element_ids.append(element_idi)
    property_ids.append(property_idi)
    n_nodes_.append(n_nodesi.ravel())
    cell_type_.append(cell_typei)
    cell_offset_.append(cell_offseti)
    return cell_offset0

def gui_material_ids(model: BDF,
                     icase: int,
                     cases: Cases,
                     form: Form,
                     element_ids: np.ndarray,
                     property_ids: np.ndarray,
                     card_index: dict[str, tuple[int, int]]) -> int:
    """
    Parameters
    ----------
    model : BDF
        the vectorized model
    icase : int
        the starting index for results
    cases : dict[int, Any]
        the cases
    element_ids : (nelement,) int ndarray
        the element ids
    property_ids : (nelement,) int ndarray
        the property ids
    card_index : dict[card_type, tuple[i0, i1]]
        card_type: str
            the type of card
        value: (i0, i1): (int, int):
            the element/property index into the relevant type

    card_index = {
        'CQUAD4': (i0, i1)
    }
    """
    cards = (
        ('Rod', False, [model.conrod], []),
        ('Rod', True, [model.crod], [model.prod]),
        ('Rod', True, [model.cbar], model.bar_property_cards),
        ('Rod', True, [model.cbeam], model.beam_property_cards),
        ('Shell', True, model.shell_element_cards, model.shell_property_cards),
        ('Solid', True, model.solid_element_cards, model.solid_property_cards),
    )
    #beams = [model.cbeam if model.cbeam.n > 0]
    #solids = [card for card in model.solid_element_cards if card.n > 0]
    #shells = [card for card in model.shell_element_cards if card.n > 0]
    idtype = element_ids.dtype
    fdtype = model.fdtype
    log = model.log

    materials_dict = {}
    neid = len(element_ids)
    for (base_flag, is_pid, element_cards, property_cards) in cards:
        element_cards2 = [card for card in element_cards if card.n > 0]
        property_cards2 = [card for card in property_cards if card.n > 0]

        for card in element_cards2:
            flag = base_flag
            etype = card.type
            if etype not in card_index:
                continue

            i0, i1 = card_index[etype]
            #eids = element_ids[i0:i1]
            pids = property_ids[i0:i1]
            if not is_pid:
                continue

            if flag in materials_dict:
                material_dicti = materials_dict[flag]
            else:
                if flag == 'Rod':
                    material_id = np.full(neid, 0, dtype=idtype)
                elif flag == 'Shell':
                    #material_id = np.array([], dtype=idtype)
                    pass
                    #pshell_material_id = np.full((neid, 4), 0, dtype=idtype)
                    #pcomp_material_id = np.full((neid, 4), 0, dtype=idtype)
                elif flag == 'Solid':
                    material_id = np.full(neid, 0, dtype=idtype)
                else:  # pragma: no cover
                    raise RuntimeError(flag)

                if flag != 'Shell':
                    material_dicti = {
                        'mid': material_id,
                    }
            for card_ in property_cards2:
                card_pids = card_.property_id
                index_pid = np.array([(i, pid) for i, pid in enumerate(pids)
                                      if pid in card_pids])
                if index_pid.shape[0] == 0:
                    continue
                index = index_pid[:, 0]
                pidsi = index_pid[:, 1]
                cardi = card_.slice_card_by_id(pidsi)
                ptype = card_.type
                if flag in {'Rod'}:
                    material_idi = cardi.material_id
                    material_id[index] = material_idi
                elif flag == 'Solid' and ptype in {'PSOLID', 'PLSOLID'}:
                    material_idi = cardi.material_id
                    material_id[index] = material_idi

                elif ptype in {'PCOMP', 'PCOMPG', 'PCOMPS'}:
                    upids = np.unique(cardi.property_id)
                    ucard = card_.slice_card_by_id(upids)
                    #nlayer = ucard.nlayer.max()
                    nlayer = card_.nlayer.max()
                    #print('nlayer =', nlayer)
                    if ptype in {'PCOMP', 'PCOMPG'}:
                        flag = 'Shell - Composite'
                    elif ptype == 'PCOMPS':
                        flag = 'Solid - Composite'
                    else:
                        raise RuntimeError(ptype)

                    if flag in materials_dict:
                        material_dicti = materials_dict[flag]
                    else:
                        # hack...
                        material_id = np.full((neid, nlayer), -1, dtype=idtype)
                        thickness = np.full((neid, nlayer), np.nan, dtype=fdtype)
                        theta = np.full((neid, nlayer), np.nan, dtype=fdtype)
                        material_dicti = {
                            'mid': material_id,
                            't': thickness,
                            'theta': theta,
                        }
                    material_id = material_dicti['mid']
                    thickness = material_dicti['t']
                    theta = material_dicti['theta']

                    # Presumably the unique number of layers is small, so we
                    # iterate on that. We can then reshape the material_id once
                    # we slice off the unwanted rows
                    for ilayer in np.unique(ucard.nlayer):
                        ipid = np.where(ilayer == cardi.nlayer)[0]
                        ieid = index[ipid]
                        npid = len(ipid)
                        pids_layer = cardi.property_id[ipid]
                        cardii = card_.slice_card_by_id(pids_layer)
                        material_id_rect = cardii.material_id.reshape(npid, ilayer)
                        thickness_rect = cardii.thickness.reshape(npid, ilayer)
                        theta_rect = cardii.theta.reshape(npid, ilayer)
                        material_id[ieid, :ilayer] = material_id_rect
                        thickness[ieid, :ilayer] = thickness_rect
                        theta[ieid, :ilayer] = theta_rect
                    del material_id, cardii, ieid, ipid, npid
                    del material_id_rect, thickness_rect
                elif ptype == 'PSHELL':
                    flag = 'Shell - PSHELL'
                    if flag in materials_dict:
                        material_dicti = materials_dict[flag]
                    else:
                        # hack...
                        material_id = np.full((neid, 4), -1, dtype=idtype)
                        thickness = np.full(neid, np.nan, dtype=fdtype)
                        material_dicti = {
                            'mid': material_id,
                            't': thickness,
                        }
                    material_id = material_dicti['mid']
                    thickness = material_dicti['t']
                    material_idi = cardi.material_id
                    material_id[index] = material_idi
                    thickness[index] = cardi.total_thickness()
                    del material_id, material_idi, thickness, index

                else:  # pragma: no cover
                    raise RuntimeError(flag)
                materials_dict[flag] = material_dicti
                #if base_flag == 'Shell':
                    #del material_id
            del material_dicti
            log.debug(f'finished {etype} {flag}')

    #materials_dict2 = {}
    #for flag, material_dict in materials_dict.items():
        #mid_col_max = material_id.max(axis=0)
        #imat = (mid_col_max > 0)
        #material_id2 = material_id[:, imat]
        #materials_dict2[flag] = material_id2
    #del materials_dict

    subcase_id = 0
    pshell_result_name_map = {
        0: 'MID1 - Membrane',
        1: 'MID2 - Bending',
        2: 'MID3 - Shear',
        3: 'MID4 - Membrane/Bending',
    }

    for flag, material_dicti in materials_dict.items():
        material_id = material_dicti['mid']

        if 'Shell' in flag:
            mid_col_max = material_id.max(axis=0)
            imat = (mid_col_max > 0)
            #material_id2 = material_id[:, imat]
            #materials_dict2[flag] = material_id2

        #print(flag, material_id)
        materials_form = []
        if material_id.ndim == 1:
            material_id = material_id.reshape((neid, 1))

        nlayer = material_id.shape[1]
        for ilayer in range(nlayer):
            if flag in {'Rod', 'Solid'}:
                result_name = 'Material'
            elif 'Composite' in flag:
                result_name = f'Material Layer {ilayer+1}'

            elif 'PSHELL' in flag:
                if ilayer == 1:
                    result_namei = 'Total Thickness'
                    resulti = material_dicti['t']
                    icase = _add_finite_centroidal_gui_result(
                        icase, cases, materials_form, subcase_id,
                        result_namei, resulti)

                result_name = pshell_result_name_map[ilayer]
            else:  #  pragma: no cover
                raise RuntimeError(flag)

            material_idi = material_id[:, ilayer]
            material_idi[material_idi == 0] = -1
            if material_idi.max() == -1:
                continue
            icase = _add_integer_centroid_gui_result(
                icase, cases, materials_form, subcase_id,
                result_name, material_idi, mask_value=-1)

        if 'Composite' in flag:
            for ilayer in range(nlayer):
                result_namei = f'Thickness Layer {ilayer+1}'
                resulti = material_dicti['t'][:, ilayer]
                icase = _add_finite_centroidal_gui_result(
                    icase, cases, materials_form, subcase_id,
                    result_namei, resulti)
            for ilayer in range(nlayer):
                result_namei = f'Theta {ilayer+1}'
                resulti = material_dicti['theta'][:, ilayer]
                icase = _add_finite_centroidal_gui_result(
                    icase, cases, materials_form, subcase_id,
                    result_namei, resulti)

        if len(materials_form):
            form.append((flag, None, materials_form))
    return icase

def _mean_min_edge_length(min_edge_length: np.ndarray) -> float:
    """
    gets the mean edge length for all the elements to size the
    picker/scale factor
    """
    ifinite = np.isfinite(min_edge_length)
    mean_edge_length_finite = min_edge_length[ifinite]
    mean_edge_length = 1.
    if len(mean_edge_length_finite):
        mean_edge_length = mean_edge_length_finite.mean()
    return mean_edge_length

def _set_quality(icase: int, cases: dict[int, Any],
                 nastran_settings: NastranSettings,
                 model: BDF, subcase_id: int,
                 gui_elements: set[str],
                 element_id: np.ndarray,
                 nelements: int,
                 cards_to_read) -> tuple[float, int, Form]:
    log = model.log
    #cards_to_read = {
        #'CBAR', 'CBEAM', 'CBUSH',
        #'CQUAD4', 'CTRIA3', 'CSHEAR',
        #'CTETRA', 'CPENTA', 'CPYRAM', 'CHEXA'}
    quality_form: Form = []
    try:
        (element_id_quality, taper_ratio, area_ratio, max_skew, aspect_ratio,
         min_theta, max_theta, dideal_theta, min_edge_length, max_warp,
         #icard_type,
        ) = model.quality(cards_to_read=cards_to_read)
    except IndexError:
        mean_edge_length = 1.0
        return mean_edge_length, icase, quality_form

    mean_edge_length = _mean_min_edge_length(min_edge_length)
    if not nastran_settings.is_element_quality:
        return mean_edge_length, icase, quality_form


    if not np.array_equal(element_id, element_id_quality):
        missing_ids = np.setdiff1d(element_id_quality, element_id)
        extra_ids = np.setdiff1d(element_id, element_id_quality)
        raise RuntimeError(f'quality map error; missing_eids={missing_ids}; extra_eids={extra_ids}')

    # not included in the vtk mesh
    NO_ELEMENT = {'CONM1', 'CONM2'}
    #is_valid_element = np.ones(nelements, dtype='bool')

    if 0:
        is_valid_list = []
        for elem in model.element_cards:
            if elem.n == 0:
                continue
            if elem.type in NO_ELEMENT:
                log.info(f'n {elem.type} {elem.n}')
                is_validi = np.zeros(elem.n, dtype='bool')
            else:
                log.info('fy {elem.type} {elem.n}')
                is_validi = np.ones(elem.n, dtype='bool')
            is_valid_list.append(is_validi)
        is_valid_element = np.hstack(is_valid_list)
        nvalid_elements = is_valid_element.sum()
        assert nvalid_elements == nelements, (nvalid_elements, nelements)


    #is_active_element = ~np.isnan(icard_type)
    #nactive_element = is_active_element.sum()
    #assert nactive_element == nelements, (nactive_element, nelements)

    icase = _add_finite_centroidal_gui_result(icase, cases, quality_form, subcase_id, 'AspectRatio', aspect_ratio)
    icase = _add_finite_centroidal_gui_result(icase, cases, quality_form, subcase_id, 'TaperRatio', taper_ratio)
    icase = _add_finite_centroidal_gui_result(icase, cases, quality_form, subcase_id, 'AreaRatio', area_ratio)
    icase = _add_finite_centroidal_gui_result(icase, cases, quality_form, subcase_id, 'Max Skew', max_skew)
    icase = _add_finite_centroidal_gui_result(icase, cases, quality_form, subcase_id, 'Min Theta', min_theta)
    icase = _add_finite_centroidal_gui_result(icase, cases, quality_form, subcase_id, 'Max Theta', max_theta)
    icase = _add_finite_centroidal_gui_result(icase, cases, quality_form, subcase_id, 'dIdeal Theta', dideal_theta)
    icase = _add_finite_centroidal_gui_result(icase, cases, quality_form, subcase_id, 'Min Edge Length', min_edge_length)
    icase = _add_finite_centroidal_gui_result(icase, cases, quality_form, subcase_id, 'Max Warp', max_warp)
    #icase = _add_finite_centroidal_gui_result(icase, cases, quality_form, subcase_id, 'icard_type', icard_type[is_valid_element])

    return mean_edge_length, icase, quality_form


def _load_stress_strain(element_cards: list,
                        bdf_model: BDF,
                        op2_model: OP2,
                        name: str,
                        key: NastranKey,
                        subcase: Subcase,
                        subcase_form: Form,
                        cases: Cases, icase: int,
                        all_element_ids: np.ndarray,
                        is_stress: bool=True) -> int:
    """
    loads:
    """
    if not hasattr(bdf_model, 'celas1'):
        return icase

    is_strain = not is_stress
    # collect all the stresses
    #result_form = _get_res_form(subcase_form, 'Stress', ntimes)
    result_form = []
    if is_stress:
        stress_form = ('Stress', None, result_form)
        results = op2_model.op2_results.stress
        name_results = {
            # name, result_dicts
            'CELAS1' : (results.celas1_stress, ),
            'CELAS2' : (results.celas2_stress, ),
            'CELAS3' : (results.celas3_stress, ),
            'CELAS4' : (results.celas4_stress, ),

            #('CDAMP1', results.cdamp1_stress),
            #('CDAMP2', results.cdamp2_stress),
            #('CDAMP3', results.cdamp3_stress),
            #('CDAMP4', results.cdamp4_stress),

            'CROD' : (results.crod_stress, ),
            'CTUBE': (results.ctube_stress, ),
            'CONROD': (results.conrod_stress, ),

            'CBAR': (results.cbar_stress, ),
            'CBEAM': (results.cbeam_stress, ),

            'CTRIA3': (results.ctria3_stress, results.ctria3_composite_stress),
            'CTRIA6': (results.ctria6_stress, results.ctria6_composite_stress),
            'CTRIAR': (results.ctriar_stress, results.ctriar_composite_stress),

            'CQUAD4': (results.cquad4_stress, results.cquad4_composite_stress),
            'CQUAD8': (results.cquad8_stress, results.cquad8_composite_stress),
            'CQUADR': (results.cquadr_stress, results.cquadr_composite_stress),

            'CTETRA': (results.ctetra_stress, ),
            'CPENTA': (results.cpenta_stress, ),
            'CPYRAM': (results.cpyram_stress, ),
            'CHEXA': (results.chexa_stress, ),

            #('cplstn3', self.cplstn3_stress),
            #('cplstn4', self.cplstn4_stress),
            #('cplstn6', self.cplstn6_stress),
            #('cplstn8', self.cplstn8_stress),
            #('cplsts3', self.cplsts3_stress),
            #('cplsts4', self.cplsts4_stress),
            #('cplsts6', self.cplsts6_stress),
            #('cplsts8', self.cplsts8_stress),
        }
    else:
        results = op2_model.op2_results.strain
        stress_form = ('Strain', None, result_form)
        name_results = {
            # name, result_dicts
            'CELAS1' : (results.celas1_strain, ),
            'CELAS2' : (results.celas2_strain, ),
            'CELAS3' : (results.celas3_strain, ),
            'CELAS4' : (results.celas4_strain, ),

            #('CDAMP1', results.cdamp1_strain),
            #('CDAMP2', results.cdamp2_strain),
            #('CDAMP3', results.cdamp3_strain),
            #('CDAMP4', results.cdamp4_strain),

            'CROD' : (results.crod_strain, ),
            'CTUBE': (results.ctube_strain, ),
            'CONROD': (results.conrod_strain, ),

            'CBAR': (results.cbar_strain, ),
            'CBEAM': (results.cbeam_strain, ),

            'CTRIA3': (results.ctria3_strain, results.ctria3_composite_strain),
            'CTRIA6': (results.ctria6_strain, results.ctria6_composite_strain),
            'CTRIAR': (results.ctriar_strain, results.ctriar_composite_strain),

            'CQUAD4': (results.cquad4_strain, results.cquad4_composite_strain),
            'CQUAD8': (results.cquad8_strain, results.cquad8_composite_strain),
            'CQUADR': (results.cquadr_strain, results.cquadr_composite_strain),

            'CTETRA': (results.ctetra_strain, ),
            'CPENTA': (results.cpenta_strain, ),
            'CPYRAM': (results.cpyram_strain, ),
            'CHEXA': (results.chexa_strain, ),

            #('cplstn3', self.cplstn3_strain),
            #('cplstn4', self.cplstn4_strain),
            #('cplstn6', self.cplstn6_strain),
            #('cplstn8', self.cplstn8_strain),
            #('cplsts3', self.cplsts3_strain),
            #('cplsts4', self.cplsts4_strain),
            #('cplsts6', self.cplsts6_strain),
            #('cplsts8', self.cplsts8_strain),
        }

    name_results_dict = {}
    for name, results_dicts in name_results.items():
        results_dicts2 = []
        for result_dict in results_dicts:
            if key in result_dict:
                results_dicts2.append(result_dict[key])
        if len(results_dicts2):
            name_results_dict[name] = results_dicts2

    #assert np.array_equal(all_element_ids, np.unique(all_element_ids))
    nelement0 = 0
    nelements = len(all_element_ids)
    oxx = np.full(nelements, np.nan, dtype='float64')
    oyy = np.full(nelements, np.nan, dtype='float64')
    ozz = np.full(nelements, np.nan, dtype='float64')
    txy = np.full(nelements, np.nan, dtype='float64')
    txz = np.full(nelements, np.nan, dtype='float64')
    tyz = np.full(nelements, np.nan, dtype='float64')
    omax = np.full(nelements, np.nan, dtype='float64')
    omin = np.full(nelements, np.nan, dtype='float64')
    von_mises = np.full(nelements, np.nan, dtype='float64')

    is_von_mises = True
    model_etype_map = {
        'CELAS1': bdf_model.celas1,
        'CELAS2': bdf_model.celas2,
        'CELAS3': bdf_model.celas3,
        'CELAS4': bdf_model.celas4,

        'CDAMP1': bdf_model.cdamp1,
        'CDAMP2': bdf_model.cdamp2,
        'CDAMP3': bdf_model.cdamp3,
        'CDAMP4': bdf_model.cdamp4,

        'CROD': bdf_model.crod,
        'CTUBE': bdf_model.ctube,
        'CONROD': bdf_model.conrod,

        'CBAR': bdf_model.cbar,
        'CBEAM': bdf_model.cbeam,

        'CTETRA': bdf_model.ctetra,
        'CPENTA': bdf_model.cpenta,
        'CHEXA': bdf_model.chexa,
        'CPYRAM': bdf_model.cpyram,

        'CTRIA3': bdf_model.ctria3,
        'CTRIA6': bdf_model.ctria6,
        'CTRIAR': bdf_model.ctriar,

        'CQUAD4': bdf_model.cquad4,
        'CQUAD8': bdf_model.cquad8,
        'CQUADR': bdf_model.cquadr,
    }
    for element in element_cards:
        etype = element.type
        # not added to gui yet
        if etype in {'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
                     'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', }:
            continue

        if etype not in name_results_dict:
            nelement0 += element.n
            continue
        results_cases = name_results_dict[etype]
        model_element = model_etype_map[etype]
        model_element_id = model_element.element_id
        #print(element)
        if etype in {'CTETRA', 'CHEXA', 'CPENTA', 'CPYRAM'}:
            #oxx, oyy, ozz, txy, tyz, txz, omax, omid, omin, von_mises
            for result in results_cases:
                if result.is_complex:
                    continue
                is_von_mises = result.is_von_mises
                #if isinstance(result, RealSolidStressArray):
                elements = result.element_cid[:, 0]
                ielement0 = nelement0 + np.searchsorted(model_element_id, elements)

                # keep only the centroidal nodes
                nnodesi = result.nnodes_per_element
                stress = result.data[:, 0::nnodesi, :]

                stressi = stress[0, :, :]
                oxx[ielement0] = stressi[:, 0]
                oyy[ielement0] = stressi[:, 1]
                ozz[ielement0] = stressi[:, 2]
                txy[ielement0] = stressi[:, 3]
                tyz[ielement0] = stressi[:, 4]
                txz[ielement0] = stressi[:, 5]
                omax[ielement0] = stressi[:, 6]
                #omid[ielement0] = stressi[:, 7]
                omin[ielement0] = stressi[:, 8]
                von_mises[ielement0] = stressi[:, 9]
        #elif etype in {'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4'}:
            #for result in results_cases:
                #ielement0 = nelement0 + np.searchsorted(model_element_id, result.element)
                #stressi = result.data[0, :, :]
                #oxx[ielement0] = stressi[:, 0]
        elif etype in {'CROD', 'CTUBE', 'CONROD'}:
            for result in results_cases:
                if result.is_complex:
                    continue
                #[axial, SMa, torsion, SMt]
                ielement0 = nelement0 + np.searchsorted(model_element_id, result.element)
                stressi = result.data[0, :, :]
                oxx[ielement0] = stressi[:, 0]
                txy[ielement0] = stressi[:, 2]

        elif etype in {'CBAR'}:
            for result in results_cases:
                # if result.is_complex:
                continue
                # element does not has duplicate ids for A/B side
                ielement0 = nelement0 + np.searchsorted(model_element_id, result.element)
                # data is long
                irange = np.arange(0, len(result.element), 2)

                #type=RealBarStressArray nelements=1
                #data: [1, ntotal, 15] where
                # 15=[s1a, s2a, s3a, s4a, axial, smaxa, smina, MS_tension,
                #     s1b, s2b, s3b, s4b,        smaxb, sminb, MS_compression]
                #data.shape = (3, 1, 15)
                stressi = result.data[0, :, :]

                #[s1a, s2a, s3a, s4a, axial, smaxa, smina, MS_tension,
                # s1b, s2b, s3b, s4b,        smaxb, sminb, MS_compression]
                #oxx[ielement] = stressi[:, 0]
                #txy[ielement] = stressi[:, 2]

                #stressi[ielement, 0], # SXC
                #stressi[ielement, 1], # SXD
                #stressi[ielement, 2], # SXE
                #stressi[ielement, 3], # SXF
                oxxi = stressi[irange, 4]
                omaxi = np.column_stack([
                    stressi[irange, 5],   # omax
                    stressi[irange, 12],  # omax
                ]) #.max(axis=1)
                print(omaxi.shape)
                omaxi = omaxi.max(axis=1)
                print(omaxi.shape)
                assert len(omaxi) == len(ielement0), f'len(omaxi)={len(omaxi)}; len(ielement0)={len(ielement0)}'
                omini = np.column_stack([
                    stressi[irange, 6],   # omin
                    stressi[irange, 13],  # omin
                ]).min(axis=1)
                oxx[ielement0] = oxxi
                omax[ielement0] = omaxi
                omin[ielement0] = omini
                del oxxi, omaxi, omini
                x = 1
        elif etype in {'CBEAM'}:
            for result in results_cases:
                if result.is_complex:
                    continue
                # element has duplicate ids for A/B side
                ielement0 = nelement0 + np.searchsorted(model_element_id, result.element[::2])
                irange = np.arange(0, len(result.element), 2)

                #[sxc, sxd, sxe, sxf, smax, smin, MS_tension, MS_compression]
                stressi = result.data[0, :, :]
                #oxx[ielement] = stressi[:, 0]
                #txy[ielement] = stressi[:, 2]

                #stressi[ielement, 0], # SXC
                #stressi[ielement, 1], # SXD
                #stressi[ielement, 2], # SXE
                #stressi[ielement, 3], # SXF
                omaxi = np.column_stack([
                    stressi[irange, 4],   # omax
                    stressi[irange+1, 4], # omax
                ]).max(axis=1)
                assert len(omaxi) == len(ielement0)
                omini = np.column_stack([
                    stressi[irange, 5],   # omin
                    stressi[irange+1, 5], # omin
                ]).min(axis=1)
                omax[ielement0] = omaxi
                omin[ielement0] = omini
                del omaxi, omini
        elif etype in {'CTRIA3', 'CTRIAR', 'CQUAD4', 'CQUADR',
                       'CTRIA6', 'CQUAD8'}:
            for result in results_cases:
                if result.is_complex:
                    continue

                nnodei = result.nnodes
                if result.data.shape[2] == 8:
                    # RealPlateStressArray
                    #[fiber_distance, oxx, oyy, txy, angle, omax, omin, max_shear]
                    result_eids = result.element_node[::nnodei, 0]
                    result_eids = result_eids[::2]
                    upper = result.data[0, ::2, :]
                    lower = result.data[0, 1::2, :]
                    oxxi = np.column_stack([upper[:, 0], lower[:, 0]]).max(axis=1)
                    oyyi = np.column_stack([upper[:, 1], lower[:, 1]]).max(axis=1)
                    txyi = np.column_stack([upper[:, 2], lower[:, 2]]).max(axis=1)
                elif result.data.shape[2] == 9:
                    # RealComplatePlateStressArray
                    #[o11, o22, t12, t1z, t2z, angle, major, minor, max_shear]
                    continue
                else:  # pragma: no cover
                    raise NotImplementedError(result)
                ielement0 = nelement0 + np.searchsorted(model_element_id, result_eids)
                oxx[ielement0] = oxxi
                oxx[ielement0] = oyyi
                oxx[ielement0] = txyi
                del nnodei, oxxi, oyyi, txyi, upper, lower
        else:  # pragma: no cover
            raise NotImplementedError(etype)
        assert len(oxx) == nelements
        nelement0 += element.n
        x = 1

    vm_word = 'Von Mises' if is_von_mises else 'Max Shear'
    subcase_id = 1

    o = 'ϵ' if is_strain else 'σ'
    t = 'ϵ' if is_strain else '𝜏'
    icase = _add_finite_centroidal_gui_result(
        icase, cases, result_form, subcase_id,
        o+'xx', oxx)
    icase = _add_finite_centroidal_gui_result(
        icase, cases, result_form, subcase_id,
        o+'yy', oyy)
    icase = _add_finite_centroidal_gui_result(
        icase, cases, result_form, subcase_id,
        o+'zz', ozz)
    icase = _add_finite_centroidal_gui_result(
        icase, cases, result_form, subcase_id,
        t+'xy', txy)
    icase = _add_finite_centroidal_gui_result(
        icase, cases, result_form, subcase_id,
        t+'xz', txz)
    icase = _add_finite_centroidal_gui_result(
        icase, cases, result_form, subcase_id,
        t+'yz', tyz)
    icase = _add_finite_centroidal_gui_result(
        icase, cases, result_form, subcase_id,
        o+'max', omax)
    icase = _add_finite_centroidal_gui_result(
        icase, cases, result_form, subcase_id,
        o+'min', omax)
    icase = _add_finite_centroidal_gui_result(
        icase, cases, result_form, subcase_id,
        o+vm_word, von_mises)
    #del name_results
    # if len(result_form) == 1:
    #     stress_form.append(result_form[0])
    if len(result_form):
        formi = subcase_form
        formi.append(stress_form)
    return icase
    #x = 1

def _load_oug(model: OP2,
              name: str,
              name_results_deflects: list[tuple[str, Any, str]],
              key: NastranKey,
              subcase: Subcase,
              subcase_form: Form,
              cases: Cases, icase: int,
              node_id: np.ndarray, xyz_cid0: np.ndarray,
              dim_max: float) -> int:
    """
    loads:
     - eigenvector
     - displacement
     - velocity
     - acceleration
     - temperature

    """
    del name
    #print(f'name_results = {name_results}')
    for (res_name, results, deflects_str) in name_results_deflects:
        if subcase not in results:
            continue
        case = results[subcase]
        #print('key =', res_name, deflects_str)
        subcase_id = key # [0]

        uname = f'{res_name}Results'
        ntimes = case.data.shape[0]
        res_form = _get_res_form(subcase_form, res_name, ntimes)

        for t123_offset in [0, 3]:
            if t123_offset == 0:
                title1 = res_name + ' T_XYZ'
            else:
                assert t123_offset == 3, t123_offset
                title1 = res_name + ' R_XYZ'
            #title1 = f'v2 {title1}'

            if deflects_str == 'disp':
                nastran_res2 = DisplacementResults2(
                    subcase_id, node_id, xyz_cid0, case,
                    title=res_name,
                    t123_offset=t123_offset,
                    dim_max=dim_max,
                    data_format='%g', nlabels=None, labelsize=None,
                    ncolors=None, colormap='', set_max_min=False,
                    uname=uname)
            elif deflects_str == 'force':
                force_index_to_base_title_annotation = {
                    0: {'title': 'F_', 'corner': 'F_'},
                    3: {'title': 'M_', 'corner': 'M_'},
                }
                methods_txyz_rxyz = ['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']
                nastran_res2 = ForceResults2(
                    subcase_id, node_id, xyz_cid0, case,
                    title=res_name,
                    t123_offset=t123_offset,
                    methods_txyz_rxyz=methods_txyz_rxyz,
                    index_to_base_title_annotation=force_index_to_base_title_annotation,
                    dim_max=dim_max,
                    data_format='%g', nlabels=None, labelsize=None,
                    ncolors=None, colormap='', set_max_min=False,
                    uname=res_name)
            else:  # pragma: no cover
                raise RuntimeError((name, deflects_str))

            headers2 = []
            for itime in range(ntimes):
                # mode = 2; freq = 75.9575 Hz
                dt = case._times[itime]
                header = _get_nastran_header(case, dt, itime)
                #header_dict[(key, itime)] = header
                #keys_map[key] = KeyMap(case.subtitle, case.label,
                #                       case.superelement_adaptivity_index,
                #                       case.pval_step)
                headers2.append(header)
                cases[icase] = (nastran_res2, (itime, title1))  # do I keep this???
                formii: Form = (title1, icase, [])
                #form_dict[(key, itime)].append(formii)
                #print(f'adding formii={str(formii)} icase={icase}')
                res_form.append(formii)
                icase += 1
            nastran_res2.headers = headers2
            #form_dict[(key, itime)].append(formii)
        #form.append('')
    #print(f'res_form = {res_form}')
    return icase

def _get_res_form(subcase_form: Form, res_name: str,
                 ntimes: int) -> Form:
    # if ntimes == 1:
    #     res_form = subcase_form
    # else:
    res_form: Form = []
    subcase_form.append((res_name, None, res_form))
    return res_form

def get_case_headers(case) -> tuple[bool, list[str]]:
    """
    Returns
    -------
    scale_per_time : bool
        is this transient-like
    header_names : list[str]
        the lower left description
    """
    header_names = []
    if case.analysis_code == 1:
        scale_per_time = False
        header_names.append('Static')
    elif case.analysis_code == 2:  # modal
        scale_per_time = True
        #print(case.get_stats())
        if hasattr(case, 'eigrs'):
            for mode, eigr, freq in zip(case.modes, case.eigrs, case.model_cycles):
                header_names.append(f'mode={mode} eigr={eigr:g} freq={freq:g} Hz')
        else:
            for mode in zip(case.modes):
                header_names.append(f'mode={mode[0]}')
    elif case.analysis_code == 5:  # frequency
        scale_per_time = False
        for freq in case.freqs:
            header_names.append(f'freq={freq:g} Hz')
    elif case.analysis_code == 6:  # transient
        scale_per_time = False
        for time in case.dts:
            header_names.append(f'time={time:g} sec')
    elif case.analysis_code == 9:  # complex modes
        scale_per_time = False
        if isinstance(case.eigrs, list):
            warnings.warn('converting eigrs from a list to an array')
            case.eigrs = np.array(case.eigrs)
            case.eigis = np.array(case.eigis)
        eigr = case.eigrs
        eigi = case.eigis
        damping, abs_freqs = complex_damping_frequency(eigr, eigi)
        for mode, eigri, eigii, freq, damping, in zip(case.modes, eigr, eigi, abs_freqs, damping):
            header_names.append(f'mode={mode} eigr={eigri:g} eigi={eigii:g} f={freq:g} Hz ζ={damping:g}')
    else:  # pragma: no cover
        raise RuntimeError(case.analysis_code)
    return scale_per_time, header_names

def get_vector_scales(t123: np.ndarray, ntimes: int, dim_max: float,
                      scale_per_time: bool) -> list[float]:
    scales = []
    if scale_per_time:
        for itime in range(ntimes):
            dxyzi = t123[itime, :, :]
            normi = np.linalg.norm(dxyzi, axis=1)
            dmax = normi.max()
            scales.append(dim_max * 0.1 / dmax)
    else:
        # normalized to final time
        dxyzi = t123[-1, :, :]
        normi = np.linalg.norm(dxyzi, axis=1)
        dmax = normi.max()
        scale = dim_max * 0.1 / dmax
        scales = [scale] * ntimes
    return scales

def filter_property_ids(etype: str,
                        property_id: np.ndarray,
                        pids_to_filter: list[int],
                        log: SimpleLogger):
    #i = np.where((property_id != 2) & (property_id != 71))[0]
    if len(pids_to_filter) == 0:
        # I think None causes issues b/c it will change the shape of the array
        # this is intended to just create the nominal array
        i = slice(None, None)
        return i
    base = (property_id != pids_to_filter[0])
    for pid in pids_to_filter[1:]:
        new = (property_id != pid)
        base = np.logical_and(base, new)
    i = np.where(base)[0]
    nproperties = len(property_id)

    #print("i =", i)
    ni = len(i)
    if nproperties != ni:
        log.warning(f'filtering {etype} from {nproperties} to {ni} because pids_to_filter')
    return i

def _create_solid_vtk_arrays(element: SolidElement,
                             grid_id: np.ndarray, cell_offset0: int,
                             ) -> tuple[int, np.ndarray, np.ndarray, np.ndarray]:
    """solids"""
    cell_type_tetra4 = 10
    #cell_type_tetra10 = 24
    cell_type_pyram5 = 14   # vtkPyramid().GetCellType()
    #cell_type_pyram13 = 27  # vtkQuadraticPyramid().GetCellType()
    cell_type_hexa8 = 12
    #cell_type_hexa20 = 25
    cell_type_penta6 = 13
    #cell_type_penta15 = 26

    if element.type == 'CTETRA':
        cell_type = cell_type_tetra4
        dnode = 4
    elif element.type == 'CPYRAM':
        cell_type = cell_type_pyram5
        dnode = 5
    elif element.type == 'CPENTA':
        cell_type = cell_type_penta6
        dnode = 6
    elif element.type == 'CHEXA':
        cell_type = cell_type_hexa8
        dnode = 8
    else:  # pragma: no coer
        raise NotImplementedError(element.type)
    nelement = element.n
    nodesi = element.base_nodes[:nelement, :]
    cell_offset0, n_nodesi, cell_typei, cell_offseti = create_offset_arrays(
        grid_id, nodesi,
        nelement, cell_type, cell_offset0, dnode)
    return cell_offset0, n_nodesi, cell_typei, cell_offseti

def _add_integer_node_gui_result(icase: int,
                                 cases: Cases,
                                 form: Form,
                                 subcase_id: int,
                                 result_name: str,
                                 result_array: np.ndarray,
                                 mask_value: Optional[int]=None) -> int:
    node_res = GuiResult(subcase_id, header=result_name, title=result_name,
                         location='node', scalar=result_array, mask_value=mask_value)

    cases[icase] = (node_res, (0, result_name))
    formi: Formi = (result_name, icase, [])
    form.append(formi)

    icase += 1
    return icase

def _add_integer_centroid_gui_result(icase: int,
                                     cases: Cases,
                                     form: Form,
                                     subcase_id: int,
                                     result_name: str,
                                     result_array: np.ndarray,
                                     mask_value: Optional[int]=None) -> int:
    """
    Parameters
    ----------
    mask_value : int
        a null value for elements without the result

    """
    node_res = GuiResult(subcase_id, header=result_name, title=result_name,
                         location='centroid', scalar=result_array,
                         mask_value=mask_value)

    cases[icase] = (node_res, (0, result_name))
    formi: Formi = (result_name, icase, [])
    form.append(formi)

    icase += 1
    return icase

def _add_finite_centroidal_gui_result(icase: int,
                                      cases: Cases,
                                      form: Form,
                                      subcase_id: int,
                                      result_name: str,
                                      result_array: np.ndarray) -> int:
    if np.any(np.isfinite(result_array)):
        res = GuiResult(subcase_id, header=result_name, title=result_name,
                        location='centroid', scalar=result_array)
        form.append((result_name, icase, []))
        cases[icase] = (res, (0, result_name))
        icase += 1
    return icase

def get_property_index(property_id: np.ndarray) -> np.ndarray:
    upid = np.unique(property_id)
    uproperty_id = np.zeros(len(property_id), dtype=property_id.dtype)
    for ipid, pid in enumerate(upid):
        uproperty_id[np.where(pid == property_id)] = ipid
    return uproperty_id

def load_nastran():  # pragma: no cover
    m = Nastran3(None)
    m.load_nastran3_geometry('spike.bdf')

if __name__ == '__main__':  # pragma: no cover
    load_nastran()
