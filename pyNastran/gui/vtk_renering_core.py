try:
    from vtkmodules.vtkRenderingCore import (
        vtkRenderer, vtkRenderWindow, vtkRenderWindowInteractor,
        vtkActor, vtkActor2D, vtkTextActor, vtkBillboardTextActor3D,
        vtkCamera,
        vtkDataSetMapper, vtkPolyDataMapper,
        vtkProp,
    )
except ImportError:
    print('error vtk_rendering_core')
    from vtk import (
        vtkRenderer, vtkRenderWindow, vtkRenderWindowInteractor,
        vtkActor, vtkActor2D, vtkTextActor, vtkBillboardTextActor3D,
        vtkCamera,
        vtkDataSetMapper, vtkPolyDataMapper,
        vtkProp,
    )
