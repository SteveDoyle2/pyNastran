# kills the program when you hit Cntl+C from the command line
# doesn't save the current state as presumably there's been an error
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)


import vtk

# create a rendering window and renderer
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)

# create a renderwindowinteractor
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

# create source
source = vtk.vtkConeSource()
source.SetCenter(0, 0, 0)
source.SetResolution(100)

# mapper
mapper = vtk.vtkPolyDataMapper()
if vtk.VTK_MAJOR_VERSION <= 5:
    mapper.SetInput(source.GetOutput())
else:
    mapper.SetInputConnection(source.GetOutputPort())

# actor
actor1 = vtk.vtkActor()
actor1.SetMapper(mapper)

# outline
outline = vtk.vtkOutlineFilter()
if vtk.VTK_MAJOR_VERSION <= 5:
    outline.SetInputData(source.GetOutput())
else:
    outline.SetInputConnection(source.GetOutputPort())
mapper2 = vtk.vtkPolyDataMapper()
if vtk.VTK_MAJOR_VERSION <= 5:
    mapper2.SetInput(outline.GetOutput())
else:
    mapper2.SetInputConnection(outline.GetOutputPort())

actor2 = vtk.vtkActor()
actor2.SetMapper(mapper2)

# assign actor to the renderer
ren.AddActor(actor1)
ren.AddActor(actor2)

# enable user interface interactor
iren.Initialize()
renWin.Render()
iren.Start()
