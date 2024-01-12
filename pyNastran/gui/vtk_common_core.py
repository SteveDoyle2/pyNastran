#try:
    #import cat
    #from vtkmodules.vtkCommonCore import (
        #vtkPoints, vtkArray, vtkDataArray, vtkFloatArray,
        #vtkIdList, vtkIdTypeArray,
        #vtkMath, vtkVersion,
        #vtkTypeFloat32Array,
        #vtkUnsignedCharArray,
        #VTK_ID_TYPE, VTK_ID_TYPE,
        #VTK_INT, VTK_FLOAT,
        #VTK_FONT_FILE, VTK_VERSION, VTK_VERSION_FULL)
#except ImportError:
from vtk import (
    vtkPoints, vtkArray, vtkDataArray, vtkFloatArray,
    vtkIdList, vtkIdTypeArray,
    vtkMath, vtkVersion,
    vtkTypeFloat32Array,
    vtkUnsignedCharArray,
    VTK_ID_TYPE, VTK_ID_TYPE,
    VTK_INT, VTK_FLOAT,
    VTK_FONT_FILE, VTK_VERSION, VTK_VERSION_FULL)
