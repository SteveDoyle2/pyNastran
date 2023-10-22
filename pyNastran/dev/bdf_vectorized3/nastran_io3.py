# encoding: utf-8
from __future__ import annotations
#import os
#from itertools import count
from pathlib import PurePath
from typing import Optional, Any, TYPE_CHECKING

import numpy as np
#import vtk

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


from pyNastran.gui.gui_objects.gui_result import GuiResult# , NormalResult
#from pyNastran.gui.gui_objects.displacements import ForceTableResults, ElementalTableResults
from pyNastran.dev.bdf_vectorized3.bdf import BDF, Subcase
from .alt_actor_builder import (
    build_vtk_geometry,
    create_alt_conm2_grids, create_alt_rbe2_grids, create_alt_rbe3_grids,
    create_alt_axes,
    create_monpnt)
from pyNastran.dev.op2_vectorized3.op2_hdf5 import OP2, OP2Geom
from pyNastran.dev.op2_vectorized3.op2_hdf5 import Results
from pyNastran.utils import PathLike
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.gui.main_window import MainWindow
    from pyNastran.dev.bdf_vectorized3.bdf_interface.bdf_attributes import (
        CTETRA, CPENTA, CHEXA, CPYRAM,
    )
    #from pyNastran.dev.op2_vectorized3.bdf_interface.bdf_attributes import RBE2, CONM2, GRID


class Nastran3:
    def __init__(self, gui: MainWindow):
        self.gui = gui
        self.data_map = None
        self.model = BDF(debug=True, log=None, mode='msc')

        # the set of element types that are supported
        self.gui_elements: set[str] = set([])
        self.bar_eids = {}
        self.bar_lines = {}
        self.include_mass_in_geometry = True

    def get_nastran3_wildcard_geometry_results_functions(self):
        """
        gets the Nastran wildcard loader used in the file load menu
        """
        data = (
            'Nastran3',
            'Nastran3 (*.bdf; *.dat; *.ecd, *.h5, *.op2)', self.load_nastran3_geometry,
            #'NastranV (*.h5)', self.load_h5_results,
            'Nastran (*.op2)', self.load_op2_results,
        )
        return data

    def load_op2_results(self, op2_filename: str, plot: bool=True) -> None:
        assert isinstance(op2_filename, (str, PurePath)), op2_filename
        model = OP2(debug=True, log=None, mode='msc')
        model.read_op2(op2_filename, combine=True, build_dataframe=False,
                       skip_undefined_matrices=False, encoding=None)
        self._load_op2_results(model, plot)

    def _load_op2_results(self, model: OP2, plot: bool) -> None:
        name = 'main'
        cases: Cases = self.result_cases
        form: Form = self.get_form()
        icase = len(cases)
        self._load_results_from_model(model,
                                      form, cases, icase,
                                      name, plot)
        self.gui._finish_results_io2(name, form, cases)

    def load_h5_results(self, h5_filename: PathLike, plot: bool=True):
        model = Results(debug=True, log=None, mode='msc')
        model.read_h5(h5_filename)
        self._load_op2_results(model, plot)

    def load_nastran3_geometry(self, bdf_filename: PathLike, name: str='main', plot: bool=True):
        bdf_filename_lower = bdf_filename.lower()
        print(bdf_filename)
        if bdf_filename_lower.endswith('.op2'):
            return self.load_op2_geometry(bdf_filename)
        elif bdf_filename_lower.endswith('.h5'):
            return self.load_h5_geometry(bdf_filename)
        return self.load_bdf_geometry(bdf_filename)

    def load_nastran3_results(self, results_filename: PathLike, name: str='main', plot: bool=True):
        results_filename_lower = results_filename.lower()
        if results_filename_lower.endswith('.op2'):
            return self.load_op2_results(results_filename)
        elif results_filename_lower.endswith('.h5'):
            return self.load_h5_results(results_filename)
        raise RuntimeError(f'results_filename={results_filename!r} is not supported')

    def load_op2_geometry(self, op2_filename: PathLike, name: str='main', plot: bool=True):
        assert isinstance(op2_filename, (str, PurePath)), op2_filename
        model = OP2Geom(debug=True, log=None, mode='msc')
        model.read_op2(op2_filename)
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
                                 name: str, plot: bool):
        #self.xyz_cid0 = model.grid.xyz_cid0() # .astype('float32')
        #self.node_ids = model.grid.node_id

        xyz_cid0 = self.xyz_cid0
        node_id = self.node_id
        mmax = xyz_cid0.max(axis=0)
        mmin = xyz_cid0.min(axis=0)
        dim_max = (mmax - mmin).max()
        #nnodes = len(node_id)

        #model.transform_displacements_to_global(icd_transform, coords, xyz_cid0=None, debug=False)
        # stress
        #result_name = 'PropertyID'
        #prop_res = GuiResult(subcase_id, header=result_name, title=result_name,
                             #location='centroid', scalar=property_id, mask_value=0)

        name_results = [
            ('Displacement', model.displacements),
            ('Eigenvector', model.eigenvectors),
        ]

        subcases = set([])
        for (name, results) in name_results:
            for key, case in results.items():
                subcases.add(key)

        # SORT1
        for subcase in subcases:
            subcase_form: Form = []
            #print(subcase)
            form.append((f'Subcase {subcase}', None, subcase_form))
            icase = _load_oug(model, name, name_results,
                              key, subcase, subcase_form,
                              cases, icase,
                              node_id, xyz_cid0, dim_max)

    def load_h5_geometry(self, h5_filename: str, name: str='main', plot: bool=True):
        model = OP2Geom(debug=True, log=None, mode='msc')
        model.read_h5(h5_filename)
        ugrid, form, cases, icase = self._load_geometry_from_model(
            model, name, plot)
        self._load_results_from_model(model,
                                      form, cases, icase,
                                      name, plot)
        self.gui._finish_results_io2(name, form, cases)
        return ugrid

    def load_bdf_geometry(self, bdf_filename: str, name: str='main', plot: bool=True):
        model = BDF(debug=True, log=None, mode='msc')
        model.idtype = 'int64'
        model.read_bdf(bdf_filename)
        ugrid, form, cases, unused_icase = self._load_geometry_from_model(
            model, name, plot)
        self.gui._finish_results_io2(name, form, cases)
        return ugrid

    def _load_geometry_from_model(self, model: BDF, name: str, plot: bool):
        self.model = model
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

        #points = vtk.vtkPoints()
        #points.SetNumberOfPoints(npoints)
        #points.SetData(points_data)

        ugrid: vtkUnstructuredGrid = gui.grid
        create_alt_conm2_grids(gui, model, node_id, xyz_cid0)
        create_alt_rbe2_grids(gui, model, node_id, xyz_cid0)
        create_alt_rbe3_grids(gui, model, node_id, xyz_cid0)
        create_alt_axes(self, gui, model, node_id, xyz_cid0)
        if 0:
            create_monpnt(gui, model, node_id, xyz_cid0)

        # add alternate actors
        gui._add_alt_actors(gui.alt_grids)

        element_id, property_id = self.load_elements(ugrid, model, node_id)
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
                     property_id: np.ndarray):

        node_index = np.arange(len(node_id))
        element_index = np.arange(len(element_id))

        gui = self.gui

        # I think we specifically look up NodeID, ELementID,
        # but I think we want to look up by Index here
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
                icase, cases,
                model, subcase_id, self.gui_elements,
                element_id, nelements, self.gui_elements)
            self.mean_edge_length = mean_edge_length

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
        #style = vtk.vtkInteractorStyleImage()

        # Render and start interaction
        renderWindowInteractor.SetRenderWindow(render_window)
        renderWindowInteractor.Initialize()

        renderWindowInteractor.Start()
        #x = 1

    def load_elements(self,
                      ugrid: vtkUnstructuredGrid,
                      model: BDF,
                      grid_id: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """
        Fills the vtkUnstructuredGrid.

        Elements are added by the order in model.element_cards, not by element id.
        This makes it easier to fill the geometry/results, but harder to lookup
        a specific element id.
        """
        log = model.log
        property_ids: list[np.ndarray] = []
        element_ids: list[np.ndarray] = []
        #nodes_ids: list[np.ndarray] = []
        #nnodes: list[np.ndarray] = []
        n_nodes_ = []
        cell_type_ = []
        cell_offset_ = []

        cell_type_point = 1  # vtk.vtkVertex().GetCellType()
        cell_type_line = 3  # vtk.vtkLine().GetCellType()
        cell_type_tri3 = 5
        cell_type_tri6 = 22
        cell_type_quad4 = 9
        cell_type_quad8 = 23

        cell_offset0 = 0
        for element in model.element_cards:
            nelement = element.n
            if nelement == 0:
                continue
            #print('element')
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
            self.gui_elements = gui_elements
            etype = element.type
            #print('load', etype, nelement, element.element_id)
            if etype in basic_elements:
                #print('  basic')
                # basic elements
                # missing nodes are not allowed
                if etype in {'CONROD'}:
                    cell_type = cell_type_line
                    property_id = np.full(nelement, -1)
                elif etype in {'CROD', 'CTUBE', 'CBAR', 'CBEAM', 'CBUSH'}:
                    cell_type = cell_type_line
                    property_id = element.property_id
                elif etype in {'CTRIA3', 'CTRIAR'}:
                    cell_type = cell_type_tri3
                    property_id = element.property_id
                elif etype in {'CQUAD4', 'CQUADR'}:
                    cell_type = cell_type_quad4
                    property_id = element.property_id
                elif etype == 'CSHEAR':
                    cell_type = cell_type_quad4
                    property_id = element.property_id
                else:
                    raise NotImplementedError(element)

                element_id = element.element_id
                assert len(element_id) == len(property_id), etype
                nodesi = element.nodes

                if 0:
                    pids_to_filter = [2, 70, 71,
                                      60, 61, # skin

                                      # spar
                                      4, 5, 6, 7, 8, 9,
                                      10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                      20, 21, 22, 13, 24, 25, 26, 27, 28, 29,
                                      30, 31, 32, 13, 34, 35, 36, 37, 38, 39,
                                      #901,  # cb
                                      222]  # ribs
                    pids_to_filter = []
                    i = filter_property_ids(etype, property_id, pids_to_filter)
                    #nelement = len(i)
                    element_id = element_id[i]
                    nelement = len(element_id)
                    property_id = property_id[i]
                    nodesi = nodesi[i, :]

                nodes_indexi = np.searchsorted(grid_id, nodesi[:, :])
                dnode = nodesi.shape[1]
                nnodesi = np.ones((nelement, 1), dtype='int64') * dnode
                cell_typei = np.ones(nelement, dtype='int64') * cell_type

                cell_offseti = _cell_offset(cell_offset0, nelement, dnode)

                #nnodes
                n_nodesi = np.hstack([nnodesi, nodes_indexi])
                n_nodes_.append(n_nodesi.ravel())

                cell_type_.append(cell_typei)
                cell_offset_.append(cell_offseti)
                element_ids.append(element_id)
                property_ids.append(property_id)
                #nnodes_nodes = [nnodes, node_id]
                del cell_type
                del cell_offseti, nnodesi, nodesi
                cell_offset0 += nelement * (dnode + 1)
            elif etype in midside_elements:
                #print('  midside')
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

                log.debug(f'{etype} nelement={nelement} nfull={nfull} npartial={npartial}')
                if nfull:
                    if etype in {'CTRIA6'}:
                        cell_type = cell_type_tri6
                        nodes = element.nodes
                    elif etype == 'CTRIAX6':
                        cell_type = cell_type_tri6
                        nodes = element.nodes[:, [0, 2, 4, 1, 3, 5]]
                    elif etype == 'CQUAD8':
                        cell_type = cell_type_quad8
                        nodes = element.nodes
                    else:
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
                    if etype in {'CTRIA6'}:
                        cell_type = cell_type_tri3
                        nodes = element.nodes[:, :3]
                    elif etype == 'CTRIAX6':
                        cell_type = cell_type_tri3
                        nodes = element.nodes[:, [0, 2, 4]]
                    elif etype == 'CQUAD8':
                        cell_type = cell_type_quad8
                        nodes = element.nodes[:, :4]
                    else:
                        raise NotImplementedError(f'full {etype}')
                    is_full = np.arange(element.n, dtype='int32')
                    cell_offset0 = _save_element(
                        is_full, grid_id, cell_type, cell_offset0,
                        element_id, property_id, nodes,
                        element_ids,
                        property_ids,
                        n_nodes_,
                        cell_type_,
                        cell_offset_)
                    continue
                nodesi = element.base_nodes


            elif etype in solid_elements:
                #model.log.debug('  solid')
                cell_offset0, n_nodesi, cell_typei, cell_offseti = _create_solid_vtk_arrays(
                    element, grid_id, cell_offset0)

                n_nodes_.append(n_nodesi.ravel())
                element_ids.append(element.element_id)
                property_ids.append(element.property_id)
                cell_type_.append(cell_typei)
                cell_offset_.append(cell_offseti)
                del n_nodesi, cell_typei, cell_offseti

            elif self.include_mass_in_geometry and etype == 'CONM2' and 0:
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
                # more complicated element
                model.log.warning(f'  dropping {element}')
                continue
                #raise NotImplementedError(element.type)
            del nelement # , dnode

            if 0:
                # [number of nodes, nodes]
                n_nodes = np.hstack(n_nodes_)
                element_id = np.hstack(element_ids)
                property_id = np.hstack(property_ids)
                cell_type = np.hstack(cell_type_)  # 10, 12
                cell_offset = np.hstack(cell_offset_) # 0, 5
                #property_id = np.hstack(property_ids)

                assert len(element_id) == len(property_id)
                assert len(element_id) == len(cell_type)
                assert len(element_id) == len(cell_offset)

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
            #property_id = np.hstack(property_ids)

            assert len(element_id) == len(property_id)
            assert len(element_id) == len(cell_type)
            assert len(element_id) == len(cell_offset)

            nelement_total = len(element_id)
            #nelement_total = 5 # len(element_id)
            build_vtk_geometry(
                nelement_total, ugrid,
                n_nodes, cell_type, cell_offset)
        self.gui_elements = gui_elements
        return element_id, property_id

def _save_element(i: np.ndarray, grid_id: np.ndarray, cell_type: int, cell_offset0: int,
                  element_id,
                  property_id,
                  nodes,
                  element_ids,
                  property_ids,
                  n_nodes_,
                  cell_type_,
                  cell_offset_,
                  ):
    nodesi = nodes[i, :]
    nelement, dnode = nodesi.shape
    element_idi = element_id[i]
    property_idi = property_id[i]
    nnodesi = np.ones((nelement, 1), dtype='int32') * dnode
    nodes_indexi = np.searchsorted(grid_id, nodesi)
    n_nodesi = np.hstack([nnodesi, nodes_indexi])
    n_nodes_.append(n_nodesi.ravel())

    element_ids.append(element_idi)
    property_ids.append(property_idi)
    cell_offseti = _cell_offset(cell_offset0, nelement, dnode)

    nodes_indexi = np.searchsorted(grid_id, nodesi)
    nnodesi = np.ones((nelement, 1), dtype='int32') * dnode
    cell_typei = np.ones(nelement, dtype='int32') * cell_type
    cell_type_.append(cell_typei)
    cell_offset_.append(cell_offseti)
    cell_offset0 += nelement * (dnode + 1)
    return cell_offset0

def _cell_offset(cell_offset0: int, nelement: int, dnode: int) -> np.ndarray:
    r"""
    (nnodes+1) = 4+1 = 5
    [0, 5, 10, 15, 20, ... (nelements-1)*5]

    for 2 CQUAD4s elements (4 nodes; 5 columns including the node count of 4)
    [0, 5]
    should be length nelement
    """
    cell_offseti = cell_offset0 + np.arange(0, nelement*(dnode +1), dnode + 1)
    assert len(cell_offseti) == nelement
    return cell_offseti

def _mean_min_edge_length(min_edge_length: np.ndarray) -> float:
    ifinite = np.isfinite(min_edge_length)
    mean_edge_length_finite = min_edge_length[ifinite]
    mean_edge_length = 1.
    if len(mean_edge_length_finite):
        mean_edge_length = mean_edge_length_finite.mean()
    return mean_edge_length

def _set_quality(icase: int, cases: dict[int, Any],
                 model: BDF, subcase_id: int,
                 gui_elements: set[str],
                 element_id: np.ndarray,
                 nelements: int,
                 cards_to_read) -> tuple[float, int, Form]:

    #cards_to_read = {
        #'CBAR', 'CBEAM', 'CBUSH',
        #'CQUAD4', 'CTRIA3', 'CSHEAR',
        #'CTETRA', 'CPENTA', 'CPYRAM', 'CHEXA'}
    (element_id_quality, taper_ratio, area_ratio, max_skew, aspect_ratio,
     min_theta, max_theta, dideal_theta, min_edge_length, max_warp,
     #icard_type,
    ) = model.quality(cards_to_read=cards_to_read)

    if not np.array_equal(element_id, element_id_quality):
        raise RuntimeError('quality map error')

    # not included in the vtk mesh
    NO_ELEMENT = {'CONM1', 'CONM2'}
    #is_valid_element = np.ones(nelements, dtype='bool')

    if 0:
        is_valid_list = []
        for elem in model.element_cards:
            if elem.n == 0:
                continue
            if elem.type in NO_ELEMENT:
                model.log.info('n {elem.type} {elem.n}')
                is_validi = np.zeros(elem.n, dtype='bool')
            else:
                model.log.info('y {elem.type} {elem.n}')
                is_validi = np.ones(elem.n, dtype='bool')
            is_valid_list.append(is_validi)
        is_valid_element = np.hstack(is_valid_list)
        nvalid_elements = is_valid_element.sum()
        assert nvalid_elements == nelements, (nvalid_elements, nelements)


    #is_active_element = ~np.isnan(icard_type)
    #nactive_element = is_active_element.sum()
    #assert nactive_element == nelements, (nactive_element, nelements)

    mean_edge_length = _mean_min_edge_length(min_edge_length)

    quality_form: Form = []
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


def _load_oug(model: OP2,
              name: str, name_results,
              key,
              subcase: Subcase,
              subcase_form: Form,
              cases: Cases, icase: int,
              node_id: np.ndarray, xyz_cid0: np.ndarray,
              dim_max: float) -> int:
    for (name, results) in name_results:
        if subcase not in results:
            continue
        res_form: Form = []
        subcase_form.append((name, None, res_form))

        case = results[subcase]
        scale_per_time, header_names = get_case_headers(case)

        #print('key =', key)
        nids = case.node_gridtype[:, 0]
        inid = np.searchsorted(nids, node_id)
        #isave = (node_id[inid] == nids)
        #isave = None
        ntimes = case.data.shape[0]
        t123 = case.data[:, inid, :3]
        assert t123.shape[2] == 3, t123.shape
        #dxyz = t123[:, isave, :]
        dxyz = t123

        subcase_id = key # [0]
        titles = [name] * ntimes # legend
        headers = [name] * ntimes # sidebar
        data_formats = ['%f'] * ntimes
        unused_scalar = None
        scales = get_vector_scales(
            t123, ntimes, dim_max,
            scale_per_time)

        assert dxyz.shape[2] == 3, dxyz.shape
        res = DisplacementResults(
            subcase_id, titles, headers,
            xyz_cid0, dxyz, unused_scalar, scales,
            data_formats=data_formats, nlabels=None, labelsize=None, ncolors=None,
            colormap='jet', set_max_min=True, uname=f'{name}Results')

        for itime, header_name in enumerate(header_names):
            #assert t123.shape[1] == 3, t123.shape
            #dxyz = np.full((nnodes, 3), np.nan, dtype='float32')
            #dxyz[isave, :] = t123[isave, :]
            headers[itime] = f'Txyz {header_name}'
            formi: Form = (header_name, icase, [])
            cases[icase] = (res, (itime, name))
            res_form.append(formi)
            icase += 1

        #form_dict[(key, itime)].append(formii)
        #form.append('')
    return icase

def get_case_headers(case) -> tuple[bool, list[str]]:
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
        #cycle = freq = eigi / (2*pi)
        #radians = eigi
        scale_per_time = False
        eigrs = case.eigrs
        eigis = case.eigis
        eigrs[eigrs == -0.] = 0.
        eigis[eigis == -0.] = 0.

        denom = np.sqrt(eigrs ** 2 + eigis ** 2)
        damping = np.zeros(len(eigrs), dtype=eigrs.dtype)
        inonzero = np.where(denom != 0)[0]
        if len(inonzero):
            damping[inonzero] = -eigrs[inonzero] / denom[inonzero]

        ## not sure
        abs_freqs = np.sqrt(np.abs(eigis)) / (2 * np.pi)

        for mode, eigr, eigi, freq, damping, in zip(case.modes, eigrs, eigis, abs_freqs, damping):
            header_names.append(f'mode={mode} eigr={eigr:g} eigi={eigi:g} f={freq:g} Hz ζ={damping:g}')
    else:
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
                        pids_to_filter: list[int]):
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
        print(f'filtering {etype} from {nproperties} to {ni} because pids_to_filter')
    return i

def _create_solid_vtk_arrays(element: CTETRA | CPENTA | CHEXA | CPYRAM,
                             grid_id: np.ndarray, cell_offset0: int,
                             ) -> tuple[int, np.ndarray, np.ndarray, np.ndarray]:
    """solids"""
    cell_type_tetra4 = 10
    #cell_type_tetra10 = 24
    cell_type_pyram5 = 14   # vtk.vtkPyramid().GetCellType()
    #cell_type_pyram13 = 27  # vtk.vtkQuadraticPyramid().GetCellType()
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
    else:
        raise NotImplementedError(element.type)
    nelement = element.n
    nodesi = element.base_nodes

    nnodesi = np.ones((nelement, 1), dtype='int32') * dnode
    nodes_indexi = np.searchsorted(grid_id, nodesi[:nelement, :])
    n_nodesi = np.hstack([nnodesi, nodes_indexi])

    cell_offseti = _cell_offset(cell_offset0, nelement, dnode)

    nodes_indexi = np.searchsorted(grid_id, nodesi[:nelement, :])
    nnodesi = np.ones((nelement, 1), dtype='int32') * dnode
    cell_typei = np.ones(nelement, dtype='int32') * cell_type
    #del cell_type
    #del cell_offseti, nnodesi, nodesi
    cell_offset0 += nelement * (dnode + 1)
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
    uproperty_id = np.zeros(len(property_id), dtype='int32')
    for ipid, pid in enumerate(upid):
        uproperty_id[np.where(pid == property_id)] = ipid
    return uproperty_id

def load_nastran():  # pragma: no cover
    m = Nastran3(None)
    m.load_nastran3_geometry('spike.bdf')

if __name__ == '__main__':  # pragma: no cover
    load_nastran()
