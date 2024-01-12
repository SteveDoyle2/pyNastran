#try:
    #import cat
    #from vtkmodules.vtkRenderingCore import (
        #vtkRenderer, vtkRenderWindow, vtkRenderWindowInteractor,
        #vtkActor, vtkActor2D, vtkTextActor, vtkBillboardTextActor3D,
        #vtkCamera,
        #vtkDataSetMapper, vtkPolyDataMapper,
        #vtkProp,
    #)
    #from vtkmodules.vtkRenderingUI import vtkGenericRenderWindowInteractor
    #from vtkmodules.vtkRenderingAnnotation import vtkAxesActor
#except ImportError:
from vtk import (
    vtkRenderer, vtkRenderWindow, vtkRenderWindowInteractor,
    vtkActor, vtkActor2D, vtkTextActor, vtkBillboardTextActor3D,
    vtkCamera,
    vtkDataSetMapper, vtkPolyDataMapper,
    vtkProp,
    vtkGenericRenderWindowInteractor,
    vtkAxesActor,
)
