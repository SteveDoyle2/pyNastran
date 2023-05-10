"""
defines:
 - create_vtk_cells_of_constant_element_type(grid, elements, etype)

"""
from __future__ import annotations
import warnings
from collections import defaultdict
from typing import Optional, TYPE_CHECKING

import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk # vtk_to_numpy

from pyNastran.gui.vtk_interface import vtkUnstructuredGrid, vtkSelectionNode
from pyNastran.gui.utils.vtk.base_utils import (
    vtkConstants, numpy_to_vtk, numpy_to_vtkIdTypeArray,
    get_numpy_idtype_for_vtk)
if TYPE_CHECKING:
    from cpylog import SimpleLogger
# // Linear cells
# VTK_EMPTY_CELL = 0
# VTK_VERTEX = 1
# VTK_POLY_VERTEX = 2
# VTK_LINE = 3
# VTK_POLY_LINE = 4
# VTK_TRIANGLE = 5
# VTK_TRIANGLE_STRIP = 6
# VTK_POLYGON = 7
# VTK_PIXEL = 8
# VTK_QUAD = 9
# VTK_TETRA = 10
# VTK_VOXEL = 11
# VTK_HEXAHEDRON = 12
# VTK_WEDGE = 13
# VTK_PYRAMID = 14
# VTK_PENTAGONAL_PRISM = 15
# VTK_HEXAGONAL_PRISM = 16
#
# // Quadratic, isoparametric cells
# VTK_QUADRATIC_EDGE = 21
# VTK_QUADRATIC_TRIANGLE = 22
# VTK_QUADRATIC_QUAD = 23
# VTK_QUADRATIC_POLYGON = 36
# VTK_QUADRATIC_TETRA = 24
# VTK_QUADRATIC_HEXAHEDRON = 25
# VTK_QUADRATIC_WEDGE = 26
# VTK_QUADRATIC_PYRAMID = 27
# VTK_BIQUADRATIC_QUAD = 28
# VTK_TRIQUADRATIC_HEXAHEDRON = 29
# VTK_QUADRATIC_LINEAR_QUAD = 30
# VTK_QUADRATIC_LINEAR_WEDGE = 31
# VTK_BIQUADRATIC_QUADRATIC_WEDGE = 32
# VTK_BIQUADRATIC_QUADRATIC_HEXAHEDRON = 33
# VTK_BIQUADRATIC_TRIANGLE = 34
#
# // Cubic, isoparametric cell
# VTK_CUBIC_LINE = 35
#
# // Special class of cells formed by convex group of points
# VTK_CONVEX_POINT_SET = 41
#
# // Polyhedron cell (consisting of polygonal faces)
# VTK_POLYHEDRON = 42
#
# // Higher order cells in parametric form
# VTK_PARAMETRIC_CURVE = 51
# VTK_PARAMETRIC_SURFACE = 52
# VTK_PARAMETRIC_TRI_SURFACE = 53
# VTK_PARAMETRIC_QUAD_SURFACE = 54
# VTK_PARAMETRIC_TETRA_REGION = 55
# VTK_PARAMETRIC_HEX_REGION = 56
#
# // Higher order cells
# VTK_HIGHER_ORDER_EDGE = 60
# VTK_HIGHER_ORDER_TRIANGLE = 61
# VTK_HIGHER_ORDER_QUAD = 62
# VTK_HIGHER_ORDER_POLYGON = 63
# VTK_HIGHER_ORDER_TETRAHEDRON = 64
# VTK_HIGHER_ORDER_WEDGE = 65
# VTK_HIGHER_ORDER_PYRAMID = 66
# VTK_HIGHER_ORDER_HEXAHEDRON = 67
#
# // Arbitrary order Lagrange elements (formulated separated from generic higher order cells)
# VTK_LAGRANGE_CURVE = 68
# VTK_LAGRANGE_TRIANGLE = 69
# VTK_LAGRANGE_QUADRILATERAL = 70
# VTK_LAGRANGE_TETRAHEDRON = 71
# VTK_LAGRANGE_HEXAHEDRON = 72
# VTK_LAGRANGE_WEDGE = 73
# VTK_LAGRANGE_PYRAMID = 74
#
# // Arbitrary order Bezier elements (formulated separated from generic higher order cells)
# VTK_BEZIER_CURVE = 75
# VTK_BEZIER_TRIANGLE = 76
# VTK_BEZIER_QUADRILATERAL = 77
# VTK_BEZIER_TETRAHEDRON = 78
# VTK_BEZIER_HEXAHEDRON = 79
# VTK_BEZIER_WEDGE = 80
# VTK_BEZIER_PYRAMID = 81

ETYPES_EXPECTED_DICT = {
    # etype: nnodes
    1: 1, # vertex
    3: 2, # line
    5: 3, # ctri3
    9: 4, # cquad4
    10: 4, # ctetra4
    14: 5, # cpyram5
    12: 8, # chexa8
    13: 6, # cpenta6
    22: 6, # ctria6
    27: 13, # cpyram13
}

def numpy_to_vtk_points(nodes, points=None, dtype='<f', deep=1):
    """common method to account for vtk endian quirks and efficiently adding points"""
    assert isinstance(nodes, np.ndarray), type(nodes)
    if points is None:
        points = vtk.vtkPoints()
        try:
            nnodes, ndim = nodes.shape
        except:
            raise RuntimeError(nodes)
        assert ndim == 3, nodes.shape
        points.SetNumberOfPoints(nnodes)

    # if we're in big endian, VTK won't work, so we byte swap
    nodes = np.asarray(nodes, dtype=np.dtype(dtype))

    points_array = numpy_to_vtk(
        num_array=nodes,
        deep=deep,
        array_type=vtk.VTK_FLOAT,
    )
    points.SetData(points_array)
    return points


def _check_shape(etype: int, elements: np.ndarray, nnodes_per_element: int) -> None:
    """
    The following table lists the supported vtk types

    ID  Class.GetCellType()
    ==  ===================
    1   vtkVertex
    3   vtkLine
    5   vtkTriangle
    12  vtkHexa
    9   vtkQuad
    10  vtkTetra
    13  vtkPenta
    14  vtkPyram
    13  vtkWedge
    25  vtkQuadraticHexahedron
    26  vtkQuadraticWedge
    27  vtkQuadraticPyramid

    """
    if isinstance(etype, int):
        try:
            nnodes = ETYPES_EXPECTED_DICT[etype]
            assert nnodes_per_element == nnodes, elements.shape
        except KeyError:
            warnings.warn(f'no recommendation for etype={etype}; nnodes_per_element={nnodes_per_element}')
    else:
        raise RuntimeError(etype)


def create_vtk_cells_of_constant_element_type(grid: vtkUnstructuredGrid,
                                              elements: np.ndarray,
                                              etype: int) -> None:
    """
    Adding constant type elements is overly complicated.

    Parameters
    ----------
    grid : vtkUnstructuredGrid()
        the unstructured grid
    elements : (nelements, nnodes_per_element) int ndarray
        the elements to add
    etype : int
        VTK cell type

    Notes
    -----
    The documentation in this method is triangle-specific as it was
    developed for a tri mesh.  It's more general than that though.

    """
    nelements, nnodes_per_element = elements.shape
    _check_shape(etype, elements, nnodes_per_element)

    # We were careful about how we defined the arrays, so the data
    # is contiguous when we ravel it.  Otherwise, you need to
    # deepcopy the arrays (deep=1).  However, numpy_to_vtk isn't so
    # good, so we could use np.copy, which is better, but it's
    # ultimately unnecessary.

    #nodes = numpy_to_vtk(elements, deep=0, array_type=vtk.VTK_ID_TYPE)
    # (nnodes_per_element + 1)  # TODO: was 4; for a tri...didn't seem to crash???
    # int8:  [-128 to 127]
    # int32: [-2_147_483_648 to 2_147_483_647]  # 2.1 billion
    # int64: [-9223372036854775808 to 9223372036854775807]
    cell_offsets = np.arange(0, nelements, dtype='int32') * (nnodes_per_element + 1)
    assert len(cell_offsets) == nelements

    # Create the array of cells
    vtk_cells = vtk.vtkCellArray()

    dtype = get_numpy_idtype_for_vtk()

    elements_vtk = np.zeros((nelements, nnodes_per_element + 1), dtype=dtype)
    elements_vtk[:, 0] = nnodes_per_element # 3 nodes/tri element
    elements_vtk[:, 1:] = elements

    cells_id_type = numpy_to_vtkIdTypeArray(elements_vtk.ravel(), deep=1)
    vtk_cells.SetCells(nelements, cells_id_type)

    # Cell types
    # 5 = vtkTriangle().GetCellType()
    cell_types = np.full(nelements, etype, dtype='int8')
    vtk_cell_types = numpy_to_vtk(
        cell_types, deep=0,
        array_type=vtk.vtkUnsignedCharArray().GetDataType())

    vtk_cell_offsets = numpy_to_vtk(cell_offsets, deep=0,
                                    array_type=vtkConstants.VTK_ID_TYPE)

    grid.SetCells(vtk_cell_types, vtk_cell_offsets, vtk_cells)

def create_vtk_cells_of_constant_element_types(grid: vtkUnstructuredGrid,
                                               elements_list, etypes_list):
    """
    Adding constant type elements is overly complicated enough as in
    ``create_vtk_cells_of_constant_element_type``.  Now we extend
    this to multiple element types.

    grid : vtkUnstructuredGrid()
        the unstructured grid
    elements_list : list[elements, ...]
        elements : (nelements, nnodes_per_element) int ndarray
            the elements to add
    etypes_list : list[etype, ...]
        etype : int
            the VTK flag as defined in
            ``create_vtk_cells_of_constant_element_type``

    """
    if isinstance(etypes_list, list) and len(etypes_list) == 1:
        create_vtk_cells_of_constant_element_type(grid, elements_list[0], etypes_list[0])
        return

    dtype = get_numpy_idtype_for_vtk()

    cell_offsets_list2 = []
    cell_types_list2 = []
    elements_list2 = []
    nelements = 0
    noffsets = 0
    for element, etype in zip(elements_list, etypes_list):
        nelement, nnodes_per_element = element.shape

        nnodesp1 = nnodes_per_element + 1  # TODO: was 4; for a tri???
        cell_offset = np.arange(0, nelement, dtype='int32') * nnodesp1 + noffsets
        noffset = nelement * nnodesp1

        cell_type = np.ones(nelement, dtype='int32') * etype
        assert len(cell_offset) == nelement

        nnodesp1 = nnodes_per_element + 1
        element_vtk = np.zeros((nelement, nnodesp1), dtype=dtype)
        element_vtk[:, 0] = nnodes_per_element # 3 nodes/tri
        element_vtk[:, 1:] = element

        cell_offsets_list2.append(cell_offset)
        cell_types_list2.append(cell_type)
        elements_list2.append(element_vtk.ravel())
        nelements += nelement
        noffsets += noffset

    cell_types_array = np.hstack(cell_types_list2)
    cell_offsets_array = np.hstack(cell_offsets_list2)
    elements_array = np.hstack(elements_list2)

    # Create the array of cells
    cells_id_type = numpy_to_vtkIdTypeArray(elements_array.ravel(), deep=1)
    vtk_cells = vtk.vtkCellArray()
    vtk_cells.SetCells(nelements, cells_id_type)

    # Cell types
    vtk_cell_types = numpy_to_vtk(
        cell_types_array, deep=0,
        array_type=vtk.vtkUnsignedCharArray().GetDataType())

    vtk_cell_offsets = numpy_to_vtk(cell_offsets_array, deep=0,
                                    array_type=vtkConstants.VTK_ID_TYPE)

    grid.SetCells(vtk_cell_types, vtk_cell_offsets, vtk_cells)

def create_unstructured_point_grid(points: vtk.vtkPoints,
                                   npoints: int) -> vtkUnstructuredGrid:
    """creates a point grid"""
    ugrid = vtkUnstructuredGrid()
    ugrid.SetPoints(points)

    cell_type_vertex = vtk.vtkVertex().GetCellType()
    etypes = [cell_type_vertex]
    elements = [np.arange(npoints, dtype='int32').reshape(npoints, 1)]
    create_vtk_cells_of_constant_element_types(ugrid, elements, etypes)
    ugrid.SetPoints(points)
    ugrid.Modified()
    return ugrid


def extract_selection_node_from_grid_to_ugrid(grid: vtkUnstructuredGrid,
                                              selection_node: vtkSelectionNode) -> vtkUnstructuredGrid:
    """
    Creates a sub-UGRID from a UGRID and a vtkSelectionNode.  In other
    words, we use a selection criteria (a definition of a subset of
    points or cells) and we create a reduced model.

    """
    selection = vtk.vtkSelection()
    selection.AddNode(selection_node)

    extract_selection = vtk.vtkExtractSelection()
    extract_selection.SetInputData(0, grid)
    extract_selection.SetInputData(1, selection)
    extract_selection.Update()

    ugrid = extract_selection.GetOutput()
    return ugrid


def create_vtk_selection_node_by_point_ids(point_ids):
    id_type_array = _convert_ids_to_vtk_idtypearray(point_ids)
    selection_node = vtk.vtkSelectionNode()
    #selection_node.SetContainingCellsOn()
    #selection_node.Initialize()
    selection_node.SetFieldType(vtk.vtkSelectionNode.POINT)
    selection_node.SetContentType(vtk.vtkSelectionNode.INDICES)
    selection_node.SetSelectionList(id_type_array)
    return selection_node

def create_vtk_selection_node_by_cell_ids(cell_ids):
    id_type_array = _convert_ids_to_vtk_idtypearray(cell_ids)
    selection_node = vtk.vtkSelectionNode()
    selection_node.SetFieldType(vtk.vtkSelectionNode.CELL)
    selection_node.SetContentType(vtk.vtkSelectionNode.INDICES)
    selection_node.SetSelectionList(id_type_array)
    return selection_node

def _convert_ids_to_vtk_idtypearray(ids):
    if isinstance(ids, int):
        ids = [ids]
    #else:
        #for idi in ids:
            #assert isinstance(idi, int), type(idi)

    id_type_array = vtk.vtkIdTypeArray()
    id_type_array.SetNumberOfComponents(1)
    for idi in ids:
        id_type_array.InsertNextValue(idi)
    return id_type_array

def find_point_id_closest_to_xyz(grid: vtkUnstructuredGrid,
                                 cell_id: int,
                                 node_xyz: np.ndarray) -> Optional[int]:
    cell = grid.GetCell(cell_id)
    if cell is None:
        return
    nnodes = cell.GetNumberOfPoints()
    points = cell.GetPoints()

    try:
        point0 = points.GetPoint(0)
    except ValueError:
        #ValueError: expects 0 <= id && id < GetNumberOfPoints()
        return None

    dist_min = vtk.vtkMath.Distance2BetweenPoints(point0, node_xyz)

    imin = 0
    #point_min = point0
    for ipoint in range(1, nnodes):
        #point = array(points.GetPoint(ipoint), dtype='float32')
        #dist = norm(point - node_xyz)
        point = points.GetPoint(ipoint)
        dist = vtk.vtkMath.Distance2BetweenPoints(point, node_xyz)
        if dist < dist_min:
            dist_min = dist
            imin = ipoint
            #point_min = point
    point_id = cell.GetPointId(imin)
    return point_id

def map_element_centroid_to_node_fringe_result(
        ugrid: vtkUnstructuredGrid,
        location: str,
        log: SimpleLogger) -> tuple[bool,
                                    tuple[int, int, float, float]]:
    """
    Maps elemental fringe results to nodal fringe results.

    If you have a 5-noded CQUAD4 (e.g., centroid + 4 nodes), only the
    centroidal value will be mapped, even though you could map the
    average the nodal values instead.  It's not wrong to do it this
    way, but it could be more accurate.

    If you have CQUAD4 with centroidal only or something like strain
    energy, this will map properly.
    """
    is_passed = False
    failed_out = (None, None, None, None)
    if location == 'node':
        log.error('Not a centroidal result.')
        return is_passed, failed_out
    assert location == 'centroid', location
    #cells = ugrid.GetCells()
    #ncells = cells.GetNumberOfCells()
    # points = ugrid.GetPoints()


    cell_data = ugrid.GetCellData()
    point_data = ugrid.GetPointData()
    # nnodes = point_data.GetNumberOfPoints()  # bad
    nnodes = ugrid.GetNumberOfPoints()

    #print('get scalars')
    vtk_cell_results = cell_data.GetScalars()
    if vtk_cell_results is None:
        log.error('Expected centroidal results, but found none.  '
                  'The result has already been mapped.')
        return is_passed, failed_out
    cell_results = vtk.util.numpy_support.vtk_to_numpy(vtk_cell_results)

    icells = np.where(np.isfinite(cell_results))[0]
    if not len(icells):
        log.error('No cells found with finite results.')
        return is_passed, failed_out
    filtered_cell_results = cell_results[icells]
    out_results_node = defaultdict(list)

    # there are no NaNs in the data
    for cell_id, res in zip(icells, filtered_cell_results):
        cell = ugrid.GetCell(int(cell_id))
        #point_ids = cell.GetPointIds()
        nnodesi = cell.GetNumberOfPoints()
        #for point_id in point_ids:
        for ipoint in range(nnodesi):
            point_id = cell.GetPointId(ipoint)
            out_results_node[point_id].append(res)

    #print('averaging')
    point_results = np.full(nnodes, np.nan, dtype='float32')
    for point_id, res in out_results_node.items():
        meani = np.average(res)
        point_results[point_id] = meani
    vtk_point_results = numpy_to_vtk(point_results, deep=0, array_type=vtk.VTK_FLOAT)
    vtk_point_results.SetName('name')

    #points.SetData(points_array)
    point_data.AddArray(vtk_point_results)

    try:
        imin = np.nanargmin(point_results)
        imax = np.nanargmax(point_results)
    except ValueError:
        # nan; just fake the result
        imin = 0
        imax = 0

    max_value = point_results[imax]
    min_value = point_results[imin]

    #out = obj.get_nlabels_labelsize_ncolors_colormap(i, name)
    #nlabels, labelsize, ncolors, colormap = out


    cell_data.SetActiveScalars(None)
    point_data.SetActiveScalars('name')
    is_passed = True
    return is_passed, (imin, imax, min_value, max_value)

def update_axis_text_size(axis: vtk.vtkAxes,
                          coord_text_scale: float,
                          width: float=1.0, height: float=0.25):
    """updates the coordinate system text size"""
    # width doesn't set the width
    # it being very large (old=0.1) makes the width constraint inactive

    texts = [
        axis.GetXAxisCaptionActor2D(),
        axis.GetYAxisCaptionActor2D(),
        axis.GetZAxisCaptionActor2D(),
    ]
    # this doesn't set the width
    # this being very large (old=0.1) makes the width constraint inactive
    for text in texts:
        text.SetWidth(coord_text_scale * width)
        text.SetHeight(coord_text_scale * height)

