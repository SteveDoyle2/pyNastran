"""
defines:
 - create_vtk_cells_of_constant_element_type(grid, elements, etype)
"""
import sys
import numpy as np
import vtk
from vtk.util.numpy_support import (
    create_vtk_array, get_numpy_array_type,
    get_vtk_array_type,
)

IS_TESTING = 'test' in sys.argv[0]
VTK_VERSION = [int(val) for val in vtk.VTK_VERSION.split('.')]
if VTK_VERSION[0] < 7 and not IS_TESTING:
    msg = 'VTK version=%r is no longer supported (use vtk 7 or 8)' % vtk.VTK_VERSION
    raise NotImplementedError(msg)
elif VTK_VERSION[0] in [5, 6, 7, 8]:
    # should work in 5/6
    # tested in 7.1.1
    vtkConstants = vtk
#elif VTK_VERSION[0] == vtk_9?:
    #vtkConstants = vtk.vtkConstants
else:
    msg = 'VTK version=%r is not supported (use vtk 7 or 8)' % vtk.VTK_VERSION
    raise NotImplementedError(msg)

def numpy_to_vtk_idtype(ids):
    #self.selection_node.GetProperties().Set(vtk.vtkSelectionNode.INVERSE(), 1)
    dtype = get_numpy_idtype_for_vtk()
    ids = np.asarray(ids, dtype=dtype)
    vtk_ids = numpy_to_vtkIdTypeArray(ids, deep=0)
    return vtk_ids

def get_numpy_idtype_for_vtk():
    """This gets the numpy dtype that we need to use to make vtk not crash"""
    isize = vtk.vtkIdTypeArray().GetDataTypeSize()
    if isize == 4:
        dtype = 'int32' # TODO: can we include endian?
    elif isize == 8:
        dtype = 'int64'
    else:
        msg = 'isize=%s' % str(isize)
        raise NotImplementedError(msg)
    return dtype

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
    Adding constant type elements is overly complicated.

    Parameters
    ----------
    grid : vtk.vtkUnstructuredGrid()
        the unstructured grid
    elements : (nelements, nnodes_per_element) int ndarray
        the elements to add
    etype : int
        VTK cell type

    Notes
    -----
    The documentation in this method is triangle-specific as it was
    developed for a tri mesh.  It's more general than that though.

    1 = vtk.vtkVertex().GetCellType()
    3 = vtkLine().GetCellType()
    5 = vtkTriangle().GetCellType()
    9 = vtk.vtkQuad().GetCellType()
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

    dtype = get_numpy_idtype_for_vtk()

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
                                    array_type=vtkConstants.VTK_ID_TYPE)

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

def numpy_to_vtk(num_array, deep=0, array_type=None):
    """Converts a contiguous real numpy Array to a VTK array object.

    This function only works for real arrays that are contiguous.
    Complex arrays are NOT handled.  It also works for multi-component
    arrays.  However, only 1, and 2 dimensional arrays are supported.
    This function is very efficient, so large arrays should not be a
    problem.

    If the second argument is set to 1, the array is deep-copied from
    from numpy. This is not as efficient as the default behavior
    (shallow copy) and uses more memory but detaches the two arrays
    such that the numpy array can be released.

    WARNING: You must maintain a reference to the passed numpy array, if
    the numpy data is gc'd and VTK will point to garbage which will in
    the best case give you a segfault.

    Parameters
    ----------
    - num_array :  a contiguous 1D or 2D, real numpy array.

    Notes
    -----
    This was pulled from VTK and modified to eliminate numpy 1.14 warnings.
    VTK uses a BSD license, so it's OK to do  that.
    """
    z = np.asarray(num_array)
    if not z.flags.contiguous:
        z = np.ascontiguousarray(z)

    shape = z.shape
    assert z.flags.contiguous, 'Only contiguous arrays are supported.'
    assert len(shape) < 3, \
           "Only arrays of dimensionality 2 or lower are allowed!"
    assert not np.issubdtype(z.dtype, np.complexfloating), \
           "Complex numpy arrays cannot be converted to vtk arrays."\
           "Use real() or imag() to get a component of the array before"\
           " passing it to vtk."

    # First create an array of the right type by using the typecode.
    if array_type:
        vtk_typecode = array_type
    else:
        vtk_typecode = get_vtk_array_type(z.dtype)
    result_array = create_vtk_array(vtk_typecode)

    # Fixup shape in case its empty or scalar.
    try:
        test_var = shape[0]
    except:
        shape = (0,)

    # Find the shape and set number of components.
    if len(shape) == 1:
        result_array.SetNumberOfComponents(1)
    else:
        result_array.SetNumberOfComponents(shape[1])

    result_array.SetNumberOfTuples(shape[0])

    # Ravel the array appropriately.
    arr_dtype = get_numpy_array_type(vtk_typecode)
    if np.issubdtype(z.dtype, arr_dtype) or \
       z.dtype == np.dtype(arr_dtype):
        z_flat = np.ravel(z)
    else:
        z_flat = np.ravel(z).astype(arr_dtype)
        # z_flat is now a standalone object with no references from the caller.
        # As such, it will drop out of this scope and cause memory issues if we
        # do not deep copy its data.
        deep = 1

    # Point the VTK array to the numpy data.  The last argument (1)
    # tells the array not to deallocate.
    result_array.SetVoidArray(z_flat, len(z_flat), 1)
    if deep:
        copy = result_array.NewInstance()
        copy.DeepCopy(result_array)
        result_array = copy
    else:
        result_array._numpy_reference = z
    return result_array

def numpy_to_vtkIdTypeArray(num_array, deep=0):
    """
    Notes
    -----
    This was pulled from VTK and modified to eliminate numpy 1.14 warnings.
    VTK uses a BSD license, so it's OK to do  that.
    """
    isize = vtk.vtkIdTypeArray().GetDataTypeSize()
    dtype = num_array.dtype
    if isize == 4:
        if dtype != np.int32:
            raise ValueError(
             'Expecting a numpy.int32 array, got %s instead.' % (str(dtype)))
    else:
        if dtype != np.int64:
            raise ValueError(
             'Expecting a numpy.int64 array, got %s instead.' % (str(dtype)))
    return numpy_to_vtk(num_array, deep, vtkConstants.VTK_ID_TYPE)
