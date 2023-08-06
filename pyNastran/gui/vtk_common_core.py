try:
    from vtkmodules.vtkCommonCore import (
        vtkPoints, vtkArray, vtkDataArray, vtkFloatArray,
        vtkIdList, vtkIdTypeArray, vtkUnsignedCharArray,
        vtkMath,
        VTK_ID_TYPE, VTK_FLOAT)
except ImportError:
    print('error vtk_common_core')
    from vtk import (
        vtkPoints, vtkArray, vtkDataArray, vtkFloatArray,
        vtkIdList, vtkIdTypeArray, vtkUnsignedCharArray,
        vtkMath,
        VTK_ID_TYPE, VTK_FLOAT)
