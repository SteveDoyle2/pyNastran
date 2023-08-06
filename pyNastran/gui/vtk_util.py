try:
    from vtkmodules.util.numpy_support import (
        numpy_to_vtk, vtk_to_numpy,
        create_vtk_array, get_numpy_array_type,
        get_vtk_array_type, numpy_to_vtkIdTypeArray,
    )
except ImportError:
    from vtk.util.numpy_support import (
        numpy_to_vtk, vtk_to_numpy,
        create_vtk_array, get_numpy_array_type,
        get_vtk_array_type, numpy_to_vtkIdTypeArray,
    )
