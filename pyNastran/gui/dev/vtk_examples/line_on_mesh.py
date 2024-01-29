#!/usr/bin/env python

import vtk
import random
import numpy
from vtk import vtkPoints

from pyNastran.gui.vtk_interface import vtkTriangle, vtkPolyData
from pyNastran.gui.vtk_rendering_core import (
    vtkRenderer, vtkRenderWindow, vtkRenderWindowInteractor, vtkActor,
    vtkPolyDataMapper,
)

# Make a 32 x 32 grid
size = 32

# Define z values for the topography
topography = numpy.zeros([size, size])
for i in range(size):
    for j in range(size):
        topography[i][j] = random.randrange(0, 5)

# Define points, triangles and colors
colors = vtk.vtkUnsignedCharArray()
colors.SetNumberOfComponents(3)
points = vtkPoints()
triangles = vtk.vtkCellArray()

# Build the meshgrid manually
count = 0
for i in range(size-1):
    for j in range(size-1):

        z1 = topography[i][j]
        z2 = topography[i][j+1]
        z3 = topography[i+1][j]

        # Triangle 1
        points.InsertNextPoint(i, j, z1)
        points.InsertNextPoint(i, (j+1), z2)
        points.InsertNextPoint((i+1), j, z3)

        triangle = vtkTriangle()
        point_ids = triangle.GetPointIds()
        point_ids.SetId(0, count)
        point_ids.SetId(1, count + 1)
        point_ids.SetId(2, count + 2)

        triangles.InsertNextCell(triangle)

        z1 = topography[i][j+1]
        z2 = topography[i+1][j+1]
        z3 = topography[i+1][j]

        # Triangle 2
        points.InsertNextPoint(i, (j+1), z1)
        points.InsertNextPoint((i+1), (j+1), z2)
        points.InsertNextPoint((i+1), j, z3)

        triangle = vtkTriangle()
        point_ids = triangle.GetPointIds()
        point_ids.SetId(0, count + 3)
        point_ids.SetId(1, count + 4)
        point_ids.SetId(2, count + 5)

        count += 6

        triangles.InsertNextCell(triangle)

        # Add some color
        r = [int(i/float(size)*255),int(j/float(size)*255),0]
        if 1:
            # vtk 9.1?
            colors.InsertNextTuple(r)
            colors.InsertNextTuple(r)
            colors.InsertNextTuple(r)
            colors.InsertNextTuple(r)
            colors.InsertNextTuple(r)
            colors.InsertNextTuple(r)
        else:
            colors.InsertNextTupleValue(r)
            colors.InsertNextTupleValue(r)
            colors.InsertNextTupleValue(r)
            colors.InsertNextTupleValue(r)
            colors.InsertNextTupleValue(r)
            colors.InsertNextTupleValue(r)

# Create a polydata object
trianglePolyData = vtkPolyData()

# Add the geometry and topology to the polydata
trianglePolyData.SetPoints(points)
trianglePolyData.GetPointData().SetScalars(colors)
trianglePolyData.SetPolys(triangles)

# Clean the polydata so that the edges are shared !
cleanPolyData = vtk.vtkCleanPolyData()
cleanPolyData.SetInputData(trianglePolyData)

# Use a filter to smooth the data (will add triangles and smooth)
smooth_loop = vtk.vtkLoopSubdivisionFilter()
smooth_loop.SetNumberOfSubdivisions(3)
smooth_loop.SetInputConnection(cleanPolyData.GetOutputPort())

# Create a mapper and actor for smoothed dataset
mapper = vtkPolyDataMapper()
mapper.SetInputConnection(smooth_loop.GetOutputPort())
actor_loop = vtkActor()
actor_loop.SetMapper(mapper)
actor_loop.GetProperty().SetInterpolationToFlat()

# Update the pipeline so that vtkCellLocator finds cells !
smooth_loop.Update()

# Define a cellLocator to be able to compute intersections between lines
# and the surface
locator = vtk.vtkCellLocator()
locator.SetDataSet(smooth_loop.GetOutput())
locator.BuildLocator()

maxloop = 1000
dist = 20.0/maxloop
tolerance = 0.001

# Make a list of points. Each point is the intersection of a vertical line
# defined by p1 and p2 and the surface.
points = vtkPoints()
for i in range(maxloop):

    p1 = [2+i*dist, 16, -1]
    p2 = [2+i*dist, 16, 6]

    # Outputs (we need only pos which is the x, y, z position
    # of the intersection)
    t = vtk.mutable(0)
    pos = [0.0, 0.0, 0.0]
    pcoords = [0.0, 0.0, 0.0]
    subId = vtk.mutable(0)
    locator.IntersectWithLine(p1, p2, tolerance, t, pos, pcoords, subId)

    # Add a slight offset in z
    pos[2] += 0.01
    # Add the x, y, z position of the intersection
    points.InsertNextPoint(pos)

# Create a spline and add the points
spline = vtk.vtkParametricSpline()
spline.SetPoints(points)
functionSource = vtk.vtkParametricFunctionSource()
functionSource.SetUResolution(maxloop)
functionSource.SetParametricFunction(spline)

# Map the spline
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(functionSource.GetOutputPort())

# Define the line actor
actor = vtkActor()
actor.SetMapper(mapper)
actor.GetProperty().SetColor([1.0, 0.0, 0.0])
actor.GetProperty().SetLineWidth(3)

# Visualize
renderer = vtkRenderer()
renderWindow = vtkRenderWindow()
renderWindow.AddRenderer(renderer)
renderWindowInteractor = vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)

# Add actors and render
renderer.AddActor(actor)
renderer.AddActor(actor_loop)

renderer.SetBackground(1, 1, 1)  # Background color white
renderWindow.SetSize(800, 800)
renderWindow.Render()
renderWindowInteractor.Start()
