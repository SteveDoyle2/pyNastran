import vtk


#Create a sphere
text_source = vtk.vtkTextSource()
text_source.SetText("Hello")
text_source.SetForegroundColor(1.0, 0.0, 0.0)
text_source.BackingOn()
text_source.Update()

#Create a mapper and actor
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(text_source.GetOutputPort())

actor = vtk.vtkActor()
actor.SetMapper(mapper)

#Create a renderer, render window, and interactor
renderer = vtk.vtkRenderer()
render_window = vtk.vtkRenderWindow()
render_window.AddRenderer(renderer)
render_window_interactor = vtk.vtkRenderWindowInteractor()
render_window_interactor.SetRenderWindow(render_window)

#Add the actor to the scene
renderer.AddActor(actor)
renderer.SetBackground(1, 1, 1) # Background color white

#Render and interact
render_window.Render()
render_window_interactor.Start()