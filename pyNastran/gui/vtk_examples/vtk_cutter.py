#!/usr/bin/env python

# A simple script to demonstrate the vtkCutter function
import vtk

#Create a cube
cube = vtk.vtkCubeSource()
cube.SetXLength(40)
cube.SetYLength(30)
cube.SetZLength(20)
cube_mapper = vtk.vtkPolyDataMapper()
cube_mapper.SetInputConnection(cube.GetOutputPort())

#create a plane to cut,here it cuts in the XZ direction (xz normal=(1,0,0);XY =(0,0,1),YZ =(0,1,0)
plane = vtk.vtkPlane()
plane.SetOrigin(10,0,0)
plane.SetNormal(1,0,0)

#create cutter
cutter = vtk.vtkCutter()
cutter.SetCutFunction(plane)
cutter.SetInputConnection(cube.GetOutputPort())
cutter.Update()
cutter_mapper = vtk.vtkPolyDataMapper()
cutter_mapper.SetInputConnection(cutter.GetOutputPort())

#create plane actor
plane_actor = vtk.vtkActor()
plane_actor.GetProperty().SetColor(1.0,1,0)
plane_actor.GetProperty().SetLineWidth(2)
plane_actor.SetMapper(cutter_mapper)

#create cube actor
cube_actor = vtk.vtkActor()
cube_actor.GetProperty().SetColor(0.5,1,0.5)
cube_actor.GetProperty().SetOpacity(0.5)
cube_actor.SetMapper(cube_mapper)

#create renderers and add actors of plane and cube
ren = vtk.vtkRenderer()
ren.AddActor(plane_actor)
ren.AddActor(cube_actor)

#Add renderer to renderwindow and render
ren_win = vtk.vtkRenderWindow()
ren_win.AddRenderer(ren)
ren_win.SetSize(600, 600)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(ren_win)
ren.SetBackground(0, 0, 0)
ren_win.Render()
iren.Start()
