try:
    from vtkmodules.util.numpy_support import numpy_to_vtk, vtk_to_numpy
except ImportError:
    from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy
