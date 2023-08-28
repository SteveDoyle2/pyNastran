try:
    from vtkmodules.vtkCommonCore import (
        vtkPoints, vtkArray, vtkDataArray, vtkFloatArray,
        vtkIdList, vtkIdTypeArray, vtkUnsignedCharArray,
        vtkMath, vtkVersion,
        VTK_ID_TYPE, VTK_ID_TYPE,
        VTK_FLOAT,
        VTK_FONT_FILE, VTK_VERSION, VTK_VERSION_FULL)
except ImportError:
    print('error vtk_common_core')
    from vtk import (
        vtkPoints, vtkArray, vtkDataArray, vtkFloatArray,
        vtkIdList, vtkIdTypeArray, vtkUnsignedCharArray,
        vtkMath, vtkVersion,
        VTK_ID_TYPE, VTK_ID_TYPE,
        VTK_FLOAT,
        VTK_FONT_FILE, VTK_VERSION, VTK_VERSION_FULL)
