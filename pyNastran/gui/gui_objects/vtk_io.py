from __future__ import annotations
import os
from typing import Union, TYPE_CHECKING
import numpy as np
#from numpy import vstack, amax, amin, arange, ones, zeros, where

#VTK_TRIANGLE = 5
#from pyNastran.gui.vtk_common_core import vtkPoints
#from pyNastran.gui.vtk_interface import vtkVertex

from vtk import (
    #vtkUnstructuredGridReader,
    vtkXMLUnstructuredGridReader,
    vtkTypeFloat32Array)
from pyNastran.gui.vtk_common_core import vtkPoints, VTK_ID_TYPE
from pyNastran.gui.vtk_interface import vtkUnstructuredGrid
from pyNastran.gui.gui_objects.types import Cases, Form, Formi
from pyNastran.gui.gui_objects.gui_result import GuiResult, INT_TYPES as INT_DTYPES, REAL_TYPES as REAL_DTYPES
from pyNastran.gui.gui_objects.displacements import (
    DisplacementResults, ForceTableResults,
    #ElementalTableResults,
)
from pyNastran.gui.utils.vtk.vtk_utils import (
    #numpy_to_vtk_points,
    vtk_to_numpy, numpy_to_vtk)
#from pyNastran.gui.qt_files.colors import YELLOW_FLOAT
if TYPE_CHECKING:
    from pyNastran.gui.gui import MainWindow

cell_type_to_nnodes = {
    # celltype_int: nnodes
    1: 1, # vtkVertex
    3: 2, # vtkLine
    5: 3,# vtkTriangle
    9: 4, # vtkQuad
    10: 4, # vtkTetra
    12: 8, # vtkHexahedron
    13: 6, # vtkWedge/cpenta6
    14: 5, # vtkPyramid
    22: 6, # vtkQuadraticTriangle/TRIA6
    25: 20, # vtkQuadraticHexahedron/CHEXA20
    26: 15, # vtkQuadraticWedge/CPENTA15
    27: 13, # vtkQuadraticPyramid/PYRAM13
}
cell_type_to_offset_size = {cell_type: nnodes + 1
                            for cell_type, nnodes in cell_type_to_nnodes.items()}

#ETYPES_EXPECTED_DICT = {
    ## etype: nnodes
    #1: 1, # vertex
    #3: 2, # line
    #5: 3, # ctri3
    #9: 4, # cquad4
    #10: 4, # ctetra4
    #12: 8, # chexa8
    #13: 6, # cpenta6
    #14: 5, # cpyram5
    #22: 6, # ctria6
    #27: 13, # cpyram13
#}

Results = Union[GuiResult, DisplacementResults, ForceTableResults]
class VtkIO:
    def __init__(self, gui: MainWindow):
        self.gui = gui

    def get_vtk_wildcard_geometry_results_functions(self):
        data = (
            'VTK',
            'VTK (*.vtk, *.vtu)', self.load_vtk_geometry,
            None, None)
        return data

    def load_vtk_geometry(self, vtk_filename: str, name=None, plot: bool=True) -> None:
        model_name = name
        #skip_reading = self.remove_old_openfoam_geometry(openfoam_filename)
        #if skip_reading:
        #    return

        base, ext = os.path.splitext(vtk_filename)

        self.gui.model_type = 'vtk'
        self.gui.log.debug('vtk_filename = %s' % vtk_filename)

        geometry_form = []
        #ID = 1
        icase = 0
        cases: dict[int, Results] = {}
        subcase_id = 1

        if ext == '.vtu':
            reader = vtkXMLUnstructuredGridReader()
            reader.SetFileName(vtk_filename)
            reader.Update()  # Needed because of GetScalarRange

            ugrid: vtkUnstructuredGrid = reader.GetOutput()
            #vtkTypeInt32Array, vtkTypeFloat32Array
            self.gui.isubcase_name_map = {1: ['VTK', '']}
            icase, xyz, nnodes, node_ids = load_point_data(
                ugrid, icase, subcase_id,
                cases, geometry_form)

            icase, nelements, element_ids = load_cell_data(
                ugrid, icase, subcase_id,
                cases, geometry_form)
        else:
            raise RuntimeError(ext)

        #nnodes = model.nodes.shape[0]
        #ntris = model.tris.shape[0]
        #nquads = model.quads.shape[0]
        #nelements = ntris + nquads

        #nodes = model.nodes
        self.gui.nelements = nelements
        self.gui.nnodes = nnodes

        #print("nNodes = %s" % self.nnodes)
        #print("nElements = %s" % self.nelements)
        assert nelements > 0, nelements

        mmax = np.amax(xyz, axis=0)
        mmin = np.amin(xyz, axis=0)
        dim_max = (mmax - mmin).max()
        self.gui.create_global_axes(dim_max)
        self.gui.log.info('max = %s' % mmax)
        self.gui.log.info('min = %s' % mmin)

        #points = numpy_to_vtk_points(nodes)

        grid = self.gui.grid
        grid.SetPoints(ugrid.GetPoints())
        vtk_cells = ugrid.GetCells()
        vtk_cell_types = ugrid.GetCellTypesArray()
        cell_types = vtk_to_numpy(vtk_cell_types)
        #vtk_cell_offsets = ugrid.GetC

        ucell_types = np.unique(cell_types)
        cell_offsets = np.zeros(nelements, dtype='int32')
        for cell_type in ucell_types:
            i = np.where(cell_type == cell_types)[0]
            offset_size = cell_type_to_offset_size[cell_type]
            cell_offsets[i] = offset_size
        vtk_cell_offsets = numpy_to_vtk(cell_offsets, deep=0,
                                        array_type=VTK_ID_TYPE)
        grid.SetCells(vtk_cell_types, vtk_cell_offsets, vtk_cells)

        #grid.Allocate(self.gui.nelements, 1000)

        #elements = []
        #etypes = []
        #if ntris:
            #elements.append(tris)
            #etypes.append(5) # vtkTriangle().GetCellType()
        #if nquads:
            #elements.append(quads)
            #etypes.append(9) # vtkQuad().GetCellType()
        #create_vtk_cells_of_constant_element_types(grid, elements, etypes)

        #self.gui.nelements = nelements
        #grid.SetPoints(points)
        grid.Modified()

        # loadSurfResults - regions/loads
        self.gui.scalar_bar_actor.VisibilityOn()
        self.gui.scalar_bar_actor.Modified()

        #self.gui.isubcase_name_map = {1: ['AFLR Surface', '']}
        #cases = {}
        #ID = 1s

        #form, cases, node_ids, element_ids = self._fill_surf_case(
            #surf_filename, cases, ID, nnodes, nelements, model)
        self.gui.node_ids = node_ids
        self.gui.element_ids = element_ids
        if plot:
            self.gui._finish_results_io2(model_name, geometry_form, cases)

    def clear_vtk(self) -> None:
        pass

    #def _load_surf_results(self, openfoam_filename):
        #raise NotImplementedError()

    #def _fill_surf_case(self, surf_filename, cases, unused_ID, nnodes, nelements, model):
        #"""builds the results for the *.surf AFLR3 input file"""
        #base, ext = os.path.splitext(surf_filename)
        #assert ext == '.surf', surf_filename

        #tag_filename = base + '.tags'

        #unused_cases_new = []
        #has_tag_data = False
        #results_form = []
        #geometry_form = [
        #]
        #form = [
            #('Geometry', None, geometry_form),
        #]
        #if has_tag_data:
            #form.append(('Tag Data', None, tag_form),)

        #results_form = []
        #if len(results_form):
            #form.append(('Results', None, results_form))
        #return form, cases, nids, eids


def load_point_data(ugrid: vtkUnstructuredGrid,
                    icase: int, subcase_id: int,
                    cases: Cases,
                    geometry_form: Form) -> tuple[int, np.ndarray, int, np.ndarray]:
    vtk_points: vtkPoints = ugrid.GetPoints()
    point_data= ugrid.GetPointData()

    vtk_array_xyz: vtkTypeFloat32Array = vtk_points.GetData()
    xyz = vtk_to_numpy(vtk_array_xyz)
    nnodes = len(xyz)
    del vtk_array_xyz, vtk_points

    node_ids = np.array([], dtype='int32')
    nresults = point_data.GetNumberOfArrays()
    if nresults:
        form = []
        point_form: Formi = ('Point Data', None, form)
        geometry_form.append(point_form)

    for i in range(nresults):
        vtk_array = point_data.GetArray(i)

        array_name = vtk_array.GetName()
        dim = vtk_array.GetNumberOfComponents()
        array_type = vtk_array.GetArrayType()
        class_name = vtk_array.GetClassName()
        assert class_name in {'vtkTypeInt32Array', 'vtkTypeFloat32Array'}, (array_name, array_type, class_name)
        lower_array_name = array_name.lower()

        data = vtk_to_numpy(vtk_array)
        if lower_array_name in {'nodeid', 'nodeids'}:
            assert len(node_ids) == 0, node_ids
            node_ids = data

        if dim == 1:
            node_res = GuiResult(0, header=array_name, title=array_name,
                                location='node', scalar=data)
        elif dim == 3:
            data_type = data.dtype.str # '<c8', '<f4'
            if data_type in INT_DTYPES:
                data_format = '%i'
            elif data_type in REAL_DTYPES:
                data_format = '%i'
            #elif data_format is None:
                #data_format = '%.2f'
            else:
                raise RuntimeError(data_format)

            titles = [array_name]
            headers = [array_name]
            data_formats = [data_format]
            unused_scalar = 1
            scales = [1.]
            forces = ('force', 'loadvectors') # force, spc forces, mpc forces
            if any((key in lower_array_name for key in forces)):
                dxyz = data
                node_res = ForceTableResults(
                    subcase_id, titles, headers, dxyz, unused_scalar, scales,
                    data_formats=data_formats, nlabels=None, labelsize=None, ncolors=None,
                    colormap='jet', set_max_min=True, uname='ForceTableResults')
                node_res.validate()
            elif 'displacement' in lower_array_name:
                dxyz = data
                node_res = DisplacementResults(
                    subcase_id, titles, headers, xyz, dxyz, unused_scalar, scales,
                    data_formats=data_formats, nlabels=None, labelsize=None, ncolors=None,
                    colormap='jet', set_max_min=True, uname='DisplacementResults')
            else:
                print((array_name, array_type, class_name), headers)
                continue
            node_res.validate()
        cases[icase] = (node_res, (0, array_name))
        formi = (array_name, icase, [])
        form.append(formi)
        icase += 1

    return icase, xyz, nnodes, node_ids

def load_cell_data(ugrid: vtkUnstructuredGrid,
                   icase: int, subcase_id: int,
                   cases, geometry_form) -> tuple[int, int, np.ndarray]:
    cell_data = ugrid.GetCellData()
    nelements = 0
    element_ids = np.array([], dtype='int32')

    nresults = cell_data.GetNumberOfArrays()
    if nresults:
        form = []
        point_form = ('Cell Data', None, form)
        geometry_form.append(point_form)

    for i in range(nresults):
        vtk_array = cell_data.GetArray(i)
        array_name = vtk_array.GetName()
        dim = vtk_array.GetNumberOfComponents()
        array_type = vtk_array.GetArrayType()
        class_name = vtk_array.GetClassName()
        #expected_class_name = array_type_map[array_type]
        #assert expected_class_name == 'vtkTypeInt32Array', (array_name, array_type, class_name)
        lower_array_name = array_name.lower()

        data = vtk_to_numpy(vtk_array)

        if lower_array_name in {'elementid', 'elementids'}:
            assert len(element_ids) == 0, element_ids
            element_ids = data

        if dim == 1:
            eid_res = GuiResult(0, header=array_name, title=array_name,
                                location='centroid', scalar=data)
        else:
            print((array_name, array_type, class_name))
        cases[icase] = (eid_res, (0, array_name))
        formi = (array_name, icase, [])
        form.append(formi)
        icase += 1

    nelements = len(element_ids)
    return icase, nelements, element_ids
