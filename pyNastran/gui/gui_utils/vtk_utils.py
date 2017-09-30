"""
defines:
 - create_vtk_cells_of_constant_element_type(grid, elements, etype)
"""

import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray


def numpy_to_vtk_points(nodes, points=None, dtype='<f', deep=1):
    """
    common method to account for vtk endian quirks and efficiently adding points
    """
    assert isinstance(nodes, np.ndarray), type(nodes)
    if points is None:
        points = vtk.vtkPoints()
        nnodes = nodes.shape[0]
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


def create_vtk_cells_of_constant_element_type(grid, elements, etype):
    """
    Adding constant type elements is overly complicated.  Note that
    some of the documentation in this method is likely triangle-specific
    as it was developed for a tri mesh.

    grid : vtk.vtkUnstructuredGrid()
        the unstructured grid
    elements : (nelements, nnodes_per_element) int ndarray
        the elements to add
    etype : int
        VTK cell type

        1  = vtk.vtkVertex().GetCellType()
        3  = vtkLine().GetCellType()
        5  = vtkTriangle().GetCellType()
        9  = vtk.vtkQuad().GetCellType()
        10 = vtkTetra().GetCellType()
        #vtkPenta().GetCellType()
        #vtkHexa().GetCellType()
        #vtkPyram().GetCellType()
    """
    nelements, nnodes_per_element = elements.shape
    # We were careful about how we defined the arrays, so the data
    # is contiguous when we ravel it.  Otherwise, you need to
    # deepcopy the arrays (deep=1).  However, numpy_to_vtk isn't so
    # good, so we could use np.copy, which is better, but it's
    # ultimately unnecessary.

    #nodes = numpy_to_vtk(elements, deep=0, array_type=vtk.VTK_ID_TYPE)
    # (nnodes_per_element + 1)  # TODO: was 4; for a tri...didn't seem to crash???
    cell_offsets = np.arange(0, nelements, dtype='int32') * (nnodes_per_element + 1)
    assert len(cell_offsets) == nelements

    # Create the array of cells
    vtk_cells = vtk.vtkCellArray()

    isize = vtk.vtkIdTypeArray().GetDataTypeSize()
    if isize == 4:
        dtype = 'int32' # TODO: can we include endian?
    elif isize == 8:
        dtype = 'int64'
    else:
        msg = 'elements.dtype=%s' % str(elements.dtype)
        raise NotImplementedError(msg)

    elements_vtk = np.zeros((nelements, nnodes_per_element + 1), dtype=dtype)
    elements_vtk[:, 0] = nnodes_per_element # 3 nodes/tri element
    elements_vtk[:, 1:] = elements

    cells_id_type = numpy_to_vtkIdTypeArray(elements_vtk.ravel(), deep=1)
    vtk_cells.SetCells(nelements, cells_id_type)

    # Cell types
    # 5 = vtkTriangle().GetCellType()
    cell_types = np.ones(nelements, dtype='int32') * etype
    vtk_cell_types = numpy_to_vtk(
        cell_types, deep=0,
        array_type=vtk.vtkUnsignedCharArray().GetDataType())

    vtk_cell_offsets = numpy_to_vtk(cell_offsets, deep=0,
                                    array_type=vtk.VTK_ID_TYPE)

    grid.SetCells(vtk_cell_types, vtk_cell_offsets, vtk_cells)

def create_vtk_cells_of_constant_element_types(grid, elements_list,
                                               etypes_list):
    """
    Adding constant type elements is overly complicated enough as in
    ``create_vtk_cells_of_constant_element_type``.  Now we extend
    this to multiple element types.

    grid : vtk.vtkUnstructuredGrid()
        the unstructured grid
    elements_list : List[elements, ...]
        elements : (nelements, nnodes_per_element) int ndarray
            the elements to add
    etypes_list : List[etype, ...]
        etype : int
            the VTK flag as defined in
            ``create_vtk_cells_of_constant_element_type``
    """
    if isinstance(etypes_list, list) and len(etypes_list) == 1:
        return create_vtk_cells_of_constant_element_type(
            grid, elements_list[0], etypes_list[0])

    isize = vtk.vtkIdTypeArray().GetDataTypeSize()
    if isize == 4:
        dtype = 'int32' # TODO: can we include endian?
    elif isize == 8:
        dtype = 'int64'
    else:
        msg = 'elements_list[0].dtype=%s' % str(elements_list[0].dtype)
        raise NotImplementedError(msg)

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
                                    array_type=vtk.VTK_ID_TYPE)

    grid.SetCells(vtk_cell_types, vtk_cell_offsets, vtk_cells)
