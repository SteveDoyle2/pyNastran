try:
    import vtkmodules.vtkCommonCore as vtk_core
    from vtkmodules.vtkCommonCore import (
        vtkPoints, vtkArray, vtkDataArray, vtkFloatArray,
        vtkIdList, vtkIdTypeArray,
        vtkMath, vtkVersion,
        vtkTypeFloat32Array,
        vtkUnsignedCharArray,
        VTK_ID_TYPE, VTK_ID_TYPE,
        VTK_INT, VTK_FLOAT,
        VTK_FONT_FILE, VTK_VERSION)
except ModuleNotFoundError:
    import vtk as vtk_core
    from vtk import (
        vtkPoints, vtkArray, vtkDataArray, vtkFloatArray,
        vtkIdList, vtkIdTypeArray,
        vtkMath, vtkVersion,
        vtkTypeFloat32Array,
        vtkUnsignedCharArray,
        VTK_ID_TYPE, VTK_ID_TYPE,
        VTK_INT, VTK_FLOAT,
        VTK_FONT_FILE, VTK_VERSION)

if not hasattr(vtk_core, 'VTK_VERSION_FULL'):  # 9.0
    VTK_VERSION_FULL = vtk_core.VTK_VERSION
else:
    VTK_VERSION_FULL = vtk_core.VTK_VERSION_FULL
