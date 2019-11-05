"""defines functions found in VTK that are overwritten for various reasons"""
import sys
import numpy as np
import vtk
from vtk.util.numpy_support import (
    create_vtk_array, get_numpy_array_type,
    get_vtk_array_type, numpy_to_vtk, numpy_to_vtkIdTypeArray
)


IS_TESTING = 'test' in sys.argv[0]
VTK_VERSION = [int(val) for val in vtk.VTK_VERSION.split('.')]
if VTK_VERSION[0] < 7:
    msg = 'VTK version=%r is no longer supported (use vtk 7 or 8)' % vtk.VTK_VERSION
    raise NotImplementedError(msg)
elif VTK_VERSION[0] in [7, 8]:
    # tested in 7.1.1
    vtkConstants = vtk
#elif VTK_VERSION[0] == vtk_9?:
    #vtkConstants = vtk.vtkConstants
else:  # pragma: no cover
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
    else:  # pragma: no cover
        msg = 'isize=%s' % str(isize)
        raise NotImplementedError(msg)
    return dtype
