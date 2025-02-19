"""Defines the GUI IO file for Tecplot."""
from __future__ import annotations
from typing import TYPE_CHECKING

import numpy as np
#from numpy import arange, mean, amax, amin, array
from pyNastran.gui.vtk_interface import vtkTriangle, vtkQuad, vtkTetra, vtkHexahedron

from pyNastran.converters.tecplot.tecplot import read_tecplot, Tecplot
#from pyNastran.converters.tecplot.utils import merge_tecplot_files
from pyNastran.gui.gui_objects.gui_result import GuiResult
from pyNastran.gui.utils.vtk.vtk_utils import numpy_to_vtk_points
from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid
from pyNastran.gui.utils.vtk.vectorized_geometry import (
    build_vtk_geometry, create_offset_arrays)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.gui.gui import MainWindow


class TecplotIO:
    def __init__(self, gui: MainWindow):
        self.gui = gui

    def _remove_old_cart3d_geometry(self, tecplot_filename: str) -> None:
        pass

    def get_tecplot_wildcard_geometry_results_functions(self):
        data = ('Tecplot',
                'Tecplot (*.dat; *.plt; *.tec)', self.load_tecplot_geometry,
                None, None)
        return data

    def load_tecplot_geometry(self, tecplot_filename: str, name: str='main',
                              plot: bool=True) -> None:
        model_name = name
        #key = self.case_keys[self.icase]
        #case = self.result_cases[key]

        skip_reading = self._remove_old_cart3d_geometry(tecplot_filename)
        #skip_reading = False
        if skip_reading:
            return

        #if 0:
            #fnames = os.listdir('time20000')
            #fnames = [os.path.join('time20000', fname) for fname in fnames]
            #model = merge_tecplot_files(fnames, tecplot_filename_out=None, log=self.log)
        #else:
        zones_to_exclude = None
        zones_to_include = None
        #zones_to_exclude = [0, 5, 6, 9, 10]
        #zones_to_include = np.array([54, 55, 56, 57, 58, 59, 60, 61]) - 1
        model = read_tecplot(tecplot_filename, log=self.gui.log, debug=False,
                             zones_to_exclude=zones_to_exclude,
                             zones_to_include=zones_to_include)

        self.gui.model_type = 'tecplot'
        self.gui.nnodes = sum([zone.nnodes for zone in model.zones])
        variables = None
        for zone in model.zones:
            variables = zone.variables
            break

        #self._make_tecplot_geometry(model, self.nnodes, quads_only=True) # cart3d
        is_surface, zone_ids = self._make_tecplot_geometry(model, quads_only=False)

        #self._create_cart3d_free_edegs(model, nodes, elements)


        # loadTecplotResults - regions/loads
        self.gui.scalar_bar_actor.VisibilityOn()
        self.gui.scalar_bar_actor.Modified()

        loads = []
        assert loads is not None
        if 'Mach' in loads:
            avg_mach = np.mean(loads['Mach'])
            note = ':  avg(Mach)=%g' % avg_mach
        else:
            note = ''
        self.gui.isubcase_name_map = {1: ['Tecplot%s' % note, '']}
        cases = {}
        ID = 1

        form, cases, node_ids, element_ids = self._fill_tecplot_case(
            cases, ID, model, variables, zone_ids, is_surface)
        self.gui.node_ids = node_ids
        self.gui.element_ids = element_ids
        self.gui._finish_results_io2(model_name, form, cases)

        #if 0:
            # http://www.vtk.org/Wiki/VTK/Examples/Cxx/Filtering/AppendFilter
            #points = vtkAppendFilter()
            #if VTK_MAJOR_VERSION <= 5:
                #appendFilter.AddInput(polydata)
                #appendFilter.AddInput(ug)
            #else:
            #appendFilter.AddInputData(polydata)
            #appendFilter.AddInputData()
            #appendFilter.Update()

    def _make_tecplot_geometry(self, model: Tecplot, quads_only=False):
        """
        Returns
        -------
        is_surface : bool
            the model is made up of only shells (not 100%)
        """
        grid = self.gui.grid
        #nnodes = self.gui.nnodes
        #print('nnodes=', nnodes)

        nodes, tris, quads, tets, hexas, zone_ids, names = model.stack_geometry()
        model.log.info(f'names = {names}')
        nquads = len(quads)
        ntris = len(tris)
        nhexas = len(hexas)
        ntets = len(tets)

        nshells = nquads + ntris
        nsolids = ntets + nhexas
        if nshells:
            is_surface = True
            grid = self.gui.grid
            _create_tecplot_shells(grid, nquads, quads, ntris, tris)
            self.gui.nelements = nshells

        elif nsolids:
            #if 0:
                #tris, quads = model.skin_elements()
                #is_tris = bool(len(tris))
                #is_quads = bool(len(quads))
                #self._create_tecplot_shells(is_quads, quads, is_tris, tris)
            #else:
            is_surface = False
            grid = self.gui.grid
            nelements = _create_tecplot_solids(
                grid, model,
                ntets, tets,
                nhexas, hexas,
                is_surface=is_surface)
            self.gui.nelements = nelements
        else:  # pragma: no cover
            raise NotImplementedError('shells or solids only\n'
                                      f'ntris={ntris} nquads={nquads}; ntets={ntets}; nhexas={nhexas}')

        #----------------------------------------------
        #print('nnodes', nnodes, inode)
        mmax = nodes.max(axis=0)
        mmin = nodes.min(axis=0)
        dim_max = (mmax - mmin).max()
        xmax, ymax, zmax = mmax
        xmin, ymin, zmin = mmin
        self.gui.log_info("xmin=%s xmax=%s dx=%s" % (xmin, xmax, xmax-xmin))
        self.gui.log_info("ymin=%s ymax=%s dy=%s" % (ymin, ymax, ymax-ymin))
        self.gui.log_info("zmin=%s zmax=%s dz=%s" % (zmin, zmax, zmax-zmin))
        self.gui.create_global_axes(dim_max)

        #print('nodes', nodes)
        #nnodes = nodes.shape[0]
        points = numpy_to_vtk_points(nodes)

        grid.SetPoints(points)
        grid.Modified()
        return is_surface, zone_ids

    def clear_tecplot(self):
        pass

    #def load_tecplot_results(self, cart3d_filename):
        #model = Cart3D(log=self.log, debug=False)
        #self.load_cart3d_geometry(cart3d_filename)

    def _fill_tecplot_case(self, cases, ID, model, variables, zone_ids, is_surface: bool):
        #'x', 'y', 'z',
        #result_names = ['rho', 'U', 'V', 'W', 'p']
        #result_names = zone.variables[3:]
        #nelements = elements.shape[0]
        nelements = self.gui.nelements
        nnodes = self.gui.nnodes
        #nnodes = zone.nnodes
        #nnodes = nodes.shape[0]

        #is_results = False
        is_results = True
        results_form = []

        if is_surface:
            element_id = 'FaceID'
        else:
            element_id = 'ElementID'

        geometry_form = [
            ('NodeID', 0, []),
            (element_id, 1, []),
            ('ZoneID', 2, []),
        ]
        assert isinstance(nnodes, int), 'nnodes=%s type=%s' % (nnodes, type(nnodes))
        assert isinstance(nelements, int), 'nelements=%s type=%s' % (nelements, type(nelements))
        assert nnodes > 0, nnodes
        assert nelements > 0, nelements

        nids = np.arange(1, nnodes + 1, dtype='int32')
        eids = np.arange(1, nelements + 1, dtype='int32')

        nid_res = GuiResult(ID, header='NodeID', title='NodeID',
                            location='node', scalar=nids)
        eid_res = GuiResult(ID, header=element_id, title=element_id,
                            location='centroid', scalar=eids)
        zone_res = GuiResult(ID, header='ZoneID', title='ZoneID',
                            location='centroid', scalar=zone_ids)


        icase = 0
        cases[icase] = (nid_res, (0, 'NodeID'))
        cases[icase + 1] = (eid_res, (0, element_id))
        cases[icase + 2] = (zone_res, (0, 'ZoneID'))
        icase += 3

        nvars = len(variables)
        if model.nzones == 1:
            zone0 = model.zones[0]
            results = zone0.nodal_results
        else:
            results = model.stack_results()
            #assert results.shape == (nnodes, nvars), results.shape

        if is_results and nvars:
            for iresult, result_name in enumerate(variables[3:]):
                #if results.shape[1] == 1:
                    #nodal_data = results
                    #assert len(variables) == 1, variables
                #else:
                nodal_data = results[:, iresult]

                node_res = GuiResult(ID, header=result_name, title=result_name,
                                     location='node', scalar=nodal_data)
                cases[icase] = (node_res, (0, result_name))

                results_form.append((result_name, icase, []))
                icase += 1
        form = [
            ('Geometry', None, geometry_form),
        ]
        if len(results_form):
            form.append(('Results', None, results_form))
        return form, cases, nids, eids

def _create_tecplot_shells(ugrid: vtkUnstructuredGrid,
                           nquad: int, quads: np.ndarray,
                           ntri: int, tris: np.ndarray):
    null_grid_id = np.array([])
    cell_offset0 = 0
    cell_offset_list = []
    cell_type_list = []
    n_nodes_list = []
    if nquad:
        #elem.GetCellType() = 9  # vtkQuad
        cell_type = 9
        dnode = 4
        cell_offset0, n_nodesi, cell_typei, cell_offseti = create_offset_arrays(
            null_grid_id, quads, nquad, cell_type, cell_offset0, dnode)
        n_nodes_list.append(n_nodesi)
        cell_type_list.append(cell_typei)
        cell_offset_list.append(cell_offseti)

    if ntri:
        #elem.GetCellType() = 5  # vtkTriangle
        cell_type = 5
        dnode = 3
        cell_offset0, n_nodesi, cell_typei, cell_offseti = create_offset_arrays(
            null_grid_id, tris, ntri, cell_type, cell_offset0, dnode)
        n_nodes_list.append(n_nodesi)
        cell_type_list.append(cell_typei)
        cell_offset_list.append(cell_offseti)

    n_nodes = np.hstack(n_nodes_list)
    cell_type = np.hstack(cell_type_list)
    cell_offset = np.hstack(cell_offset_list)
    nelement_total = nquad + ntri
    build_vtk_geometry(nelement_total, ugrid, n_nodes, cell_type, cell_offset)

def _create_tecplot_solids(ugrid: vtkUnstructuredGrid,
                           model: Tecplot,
                           ntet: int, tets: np.ndarray,
                           nhexa: int, hexas: np.ndarray,
                           is_surface: bool=True):
    """
    add a model with solid elements

    Parameters
    ----------
    is_surface : bool; default=True
        True : skin the model (good for large models, but doesn't load everything)
        False : load the model normally
    """
    nelement_total = ntet + nhexa
    null_grid_id = np.array([])
    cell_offset0 = 0
    cell_offset_list = []
    cell_type_list = []
    n_nodes_list = []
    if is_surface:
        nfaces = 0
        if nhexa:
            free_faces = np.array(model.get_free_faces(), dtype='int32')# + 1
            nfaces += len(free_faces)

            # elem.GetCellType() = 9  # vtkQuad
            cell_type = 9
            dnode = 4
            nquad = len(free_faces)
            quads = free_faces
            cell_offset0, n_nodesi, cell_typei, cell_offseti = create_offset_arrays(
                null_grid_id, quads, nquad, cell_type, cell_offset0, dnode)
            n_nodes_list.append(n_nodesi)
            cell_type_list.append(cell_typei)
            cell_offset_list.append(cell_offseti)
        nelement_total = nfaces
    else:
        # is_volume

        # vtkTetra = 10
        nelement_total = ntet + nhexa
        if ntet:
            # elem.GetCellType() = 10 # vtkTetra
            cell_type = 10
            dnode = 4
            cell_offset0, n_nodesi, cell_typei, cell_offseti = create_offset_arrays(
                null_grid_id, tets, ntet, cell_type, cell_offset0, dnode)
            n_nodes_list.append(n_nodesi)
            cell_type_list.append(cell_typei)
            cell_offset_list.append(cell_offseti)

        if nhexa:
            # elem.GetCellType() = 12 # vtkHexahedron
            cell_type = 12
            dnode = 8
            cell_offset0, n_nodesi, cell_typei, cell_offseti = create_offset_arrays(
                null_grid_id, hexas, nhexa, cell_type, cell_offset0, dnode)
            n_nodes_list.append(n_nodesi)
            cell_type_list.append(cell_typei)
            cell_offset_list.append(cell_offseti)
    assert nelement_total > 0, nelement_total

    n_nodes = np.hstack(n_nodes_list)
    cell_type = np.hstack(cell_type_list)
    cell_offset = np.hstack(cell_offset_list)
    build_vtk_geometry(nelement_total, ugrid,
                       n_nodes, cell_type, cell_offset)
    return nelement_total
