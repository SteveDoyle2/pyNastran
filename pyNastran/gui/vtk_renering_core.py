try:
    from vtkmodules.vtkRenderingCore import (
        vtkRenderer, vtkRenderWindow, vtkRenderWindowInteractor,
        vtkActor, vtkActor2D, vtkBillboardTextActor3D,
        vtkCamera,
        vtkDataSetMapper, vtkPolyDataMapper,
    )
except ImportError:
    from vtk import (
        vtkRenderer, vtkRenderWindow, vtkRenderWindowInteractor,
        vtkActor, vtkActor2D, vtkCamera,
        vtkDataSetMapper, vtkPolyDataMapper,
    )
