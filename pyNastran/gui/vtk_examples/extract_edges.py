#!/usr/bin/env python

# Purdue CS530 - Introduction to Scientific Visualization
# Fall 2013

# Simple example to illustrate the typical structure of the rendering 
# pipeline in VTK. Here, we draw a sphere and enter the interactive mode.

# Our example needs the VTK Python package
import vtk

def render_demo():
    
    # Step 1. Data Source
    #   Create a sphere geometry
    sphere_src = vtk.vtkSphereSource()
    sphere_src.SetRadius(1.0)
    sphere_src.SetCenter(0.0, 0.0, 0.0)
    #   In reality, vtkSphereSource creates a polygonal approximation (a    
    # triangulation) of a sphere. Set the resolution.
    sphere_src.SetThetaResolution(20)
    sphere_src.SetPhiResolution(20)
    
    # Step 2. Data Processing
    #   We will visualize the edges of the sphere that we just created. For 
    #   that, we must first extract the edges from the sphere triangulation.
    #   VTK makes that easy.
    edge_extractor = vtk.vtkExtractEdges()
    #   Create a pipeline link between sphere source and edge extractor
    edge_extractor.SetInputConnection(sphere_src.GetOutputPort())
    #   Now our edge extractor acts as a second data source: it supplies the 
    #   geometry of the edges of the sphere
    
    # Step 3. Mappers
    #   Create a set of graphical primitives to represent the sphere
    sphere_mapper = vtk.vtkPolyDataMapper()
    #   Create a pipeline link between sphere source and mapper
    sphere_mapper.SetInputConnection(sphere_src.GetOutputPort())
    #   Similarly, create a mapper for the edges of the sphere
    edge_mapper = vtk.vtkPolyDataMapper()
    edge_mapper.SetInputConnection(edge_extractor.GetOutputPort())
    
    # Step 4. Actors
    #   Insert the graphical primitives into the scene to be rendered
    sphere_actor = vtk.vtkActor()
    #   Assign sphere mapper to actor
    sphere_actor.SetMapper(sphere_mapper)
    #   The actor now controls the graphical properties of the sphere. We can 
    #   access them to change the color of the sphere. The color is set in RGB
    #   mode using values from 0 to 1 for each of the three color channels. 
    #   Here we set the color to orange.
    sphere_actor.GetProperty().SetColor(1, 0.5, 0)
    #   Same thing for the edges of the sphere
    edge_actor = vtk.vtkActor()
    edge_actor.SetMapper(edge_mapper)
    #   We want our edges to be drawn in dark green
    edge_actor.GetProperty().SetColor(0, 0.5, 0)
    #   We also want them to be shown as thick lines
    edge_actor.GetProperty().SetLineWidth(3)
    
    # Step 5. Renderer
    #   Render the scene to form an image that can be displayed
    my_renderer = vtk.vtkRenderer()
    #   Add all our actors to the renderer
    my_renderer.AddActor(sphere_actor)
    my_renderer.AddActor(edge_actor)
    my_renderer.SetBackground(0.1, 0.2, 0.4)
    
    # Step 6. Render Window
    #   Provides a window in which the rendered image can be displayed on the 
    #   screen
    my_window = vtk.vtkRenderWindow()
    my_window.AddRenderer(my_renderer)
    #   The window size controls the resolution of the final image
    my_window.SetSize(600, 600)
    
    # Step 7. Interactor
    #   Create an interactor: the user will be able to control the 
    #   visualization through mouse and keyboard interaction
    my_interactor = vtk.vtkRenderWindowInteractor()
    my_interactor.SetRenderWindow(my_window)
    #   IMPORTANT NOTE: always initialize the interactor before doing the 
    #   rendering!
    my_interactor.Initialize()
    
    # Step 8. Actual Rendering
    #   Draw something on the screen (finally!)
    my_window.Render()
    
    # Step 9. Interaction Starts
    #   We are done: entering the interaction loop. The user is in charge
    my_interactor.Start()

if __name__=="__main__":
      render_demo()