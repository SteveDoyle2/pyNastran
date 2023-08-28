import vtk
from pyNastran.gui.utils.vtk.base_utils import VTK_VERSION_SPLIT
from pyNastran.gui.vtk_interface import vtkUnstructuredGrid


# Create polyhedron (cube)

pointIds = vtk.vtkIdList()
pointIds.SetNumberOfIds(8)

for i in range(8):
    pointIds.InsertNextId(i)
    points = vtk.vtkPoints()

points.InsertNextPoint(-1.0,-1.0,-1.0)
points.InsertNextPoint( 1.0,-1.0,-1.0)
points.InsertNextPoint( 1.0, 1.0,-1.0)
points.InsertNextPoint(-1.0, 1.0,-1.0)
points.InsertNextPoint(-1.0,-1.0, 1.0)
points.InsertNextPoint( 1.0,-1.0, 1.0)
points.InsertNextPoint( 1.0, 1.0, 1.0)
points.InsertNextPoint(-1.0, 1.0, 1.0)


# list of pointIds that make up the face

faceList = [
    [0, 3, 2, 1],
    [0, 4, 7, 3],
    [4, 5, 6, 7],
    [5, 1, 2, 6],
    [0, 1, 5, 4],
    [2, 3, 7, 6],
]

faceId = vtk.vtkIdList()

faceId.InsertNextId(6) # Number faces that make up the cell.
for face in faceList: # Loop over all the faces
    faceId.InsertNextId(len(face)) # Number of points in face
    [faceId.InsertNextId(i) for i in face] # Insert the pointIds for the face


ugrid = vtkUnstructuredGrid()
ugrid.SetPoints(points)
ugrid.InsertNextCell(vtk.VTK_POLYHEDRON, faceId)

#-------------------------------------
if 1:
    ug = ugrid
    grid_mapper = vtk.vtkDataSetMapper()
    vtk_version = int(VTK_VERSION_SPLIT[0])
    grid_mapper.SetInputData(ug)

#-------------------------------------

#Create a mapper and actor
#mapper = vtk.vtkPolyDataMapper()
#mapper.SetInputConnection(text_source.GetOutputPort())

actor = vtk.vtkActor()
actor.SetMapper(grid_mapper)

#Create a renderer, render window, and interactor
renderer = vtk.vtkRenderer()
render_window = vtk.vtkRenderWindow()
render_window.AddRenderer(renderer)
render_window_interactor = vtk.vtkRenderWindowInteractor()
render_window_interactor.SetRenderWindow(render_window)

#Add the actor to the scene
renderer.AddActor(actor)
renderer.SetBackground(0, 0, 0) # Background color white

#Render and interact
render_window.Render()
render_window_interactor.Start()

#-------------------------------------

#writer = vtk.vtkXMLUnstructuredGridWriter()
#writer.SetInputData(ugrid)
#writer.SetFileName("polyhedron.vtu")
#writer.SetDataModeToAscii()
#writer.Update()
