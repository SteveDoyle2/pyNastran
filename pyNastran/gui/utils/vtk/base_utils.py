"""defines functions found in VTK that are overwritten for various reasons"""
import sys
import numpy as np
import vtk
from vtk.util.numpy_support import (
    create_vtk_array, get_numpy_array_type,
    get_vtk_array_type, numpy_to_vtkIdTypeArray, # numpy_to_vtk,
)


IS_TESTING = 'test' in sys.argv[0]
VTK_VERSION = [int(val) for val in vtk.VTK_VERSION.split('.')]
if VTK_VERSION[0] < 7:
    msg = f'VTK version={vtk.VTK_VERSION!r} is no longer supported (use vtk 7 or 8)'
    raise NotImplementedError(msg)
elif VTK_VERSION[0] in [7, 8, 9]:
    # tested in 7.1.1, 8.1.2, 9.0.0
    vtkConstants = vtk
#elif VTK_VERSION[0] == vtk_9?:
    #vtkConstants = vtk.vtkConstants
else:  # pragma: no cover
    msg = f'VTK version={vtk.VTK_VERSION!r} is not supported (use vtk 7, 8, or 9)'
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
    else:  # pragma: no cover
        msg = 'isize=%s' % str(isize)
        raise NotImplementedError(msg)
    return dtype


def numpy_to_vtk(num_array, deep=0, array_type=None):  # pragma: no cover
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

vtk.util.numpy_support.numpy_to_vtk = numpy_to_vtk
