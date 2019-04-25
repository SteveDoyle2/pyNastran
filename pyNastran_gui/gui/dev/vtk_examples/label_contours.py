import sys
import vtk

def main():
    pointThreshold = 10

    polyData = vtk.vtkPolyData()
    contour_filter = vtk.vtkContourFilter()

    # If a file is present, read it, otherwise generate some random
    # scalars on a plane
    nargs = len(sys.argv)
    argv = sys.argv
    if nargs == 1:
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(argv[1])
        reader.Update()

        #scalar_range[2]
        scalar_range = reader.GetOutput().GetScalarRange()
        polyData = reader.GetOutput()

        print("scalar_range: ", scalar_range[0], ", ", scalar_range[1])
        contour_filter.SetValue(0, (scalar_range[1] + scalar_range[0]) / 2.0)

        contour_filter.SetInputConnection(reader.GetOutputPort())
        if nargs == 3:
            contour_filter.SetValue(0, float(argv[2]))
        elif nargs == 4:
            contour_filter.SetValue(0, float(argv[2]))
            contour_filter.SetValue(1, float(argv[3]))
        elif nargs == 5:
            contour_filter.GenerateValues(int(argv[2]), float(argv[3]), float(argv[4]))
    else:
        plane = vtk.vtkPlaneSource()
        plane.SetXResolution(10)
        plane.SetYResolution(10)
        plane.Update()

        randomScalars = vtk.vtkDoubleArray()  # was vtkArray
        randomScalars.SetNumberOfComponents(1)
        randomScalars.SetName("Isovalues")
        for i in range(plane.GetOutput().GetNumberOfPoints()):
            randomScalars.InsertNextTuple1(vtk.vtkMath.Random(-100.0, 100.0))

        plane.GetOutput().GetPointData().SetScalars(randomScalars)
        polyData = plane.GetOutput()
        contour_filter.SetInputConnection(plane.GetOutputPort())
        contour_filter.GenerateValues(5, -100, 100)
        pointThreshold = 0


    # Connect the segments of the conours into polylines
    contour_stripper = vtk.vtkStripper()
    contour_stripper.SetInputConnection(contour_filter.GetOutputPort())
    contour_stripper.Update()

    numberOfContourLines = contour_stripper.GetOutput().GetNumberOfLines()

    print("There are %s contours lines." % numberOfContourLines)

    points = contour_stripper.GetOutput().GetPoints()
    cells = contour_stripper.GetOutput().GetLines()
    scalars = contour_stripper.GetOutput().GetPointData().GetScalars()
    npoints = points.GetNumberOfPoints()

    # Create a polydata that contains point locations for the contour
    # line labels
    label_poly_data = vtk.vtkPolyData()
    label_points = vtk.vtkPoints()
    label_scalars = vtk.vtkDoubleArray() # was vtkArray
    label_scalars.SetNumberOfComponents(1)
    label_scalars.SetName("Isovalues")

    #vtkIdType *indices
    #vtkIdType npoints
    lineCount = 0
    #for(cells.InitTraversal()
    #     cells.GetNextCell(npoints, indices)
    #     npoints++)
    #numberOfPoints = vtk.vtkIdType()
    #nstart = cells.InitTraversal()
    print('npoints = %s' % npoints)
    #print('nstart=%s' % nstart)
    include_labels = False
    if include_labels:
        indices = [i for i in range(cells.GetNumberOfCells()*10 + 1)]
        print(indices)
        #nend = cells.GetNextCell(npoints, indices)
        #print('nstart=%s nend=%s' % (nstart, nend))
        ncells = cells.GetNumberOfCells()
        #for i in range(nstart, c):
        for i in range(ncells):
            #cell = cells.GetNumber
            if npoints < pointThreshold:
                continue

            print("Line ", lineCount, ": ")

            # Compute the point id to hold the label
            # Mid point or a random point
            mid_point_id = indices[npoints / 2]
            val = int(vtk.vtkMath.Random(0, npoints))
            mid_point_id = indices[val]

            # mid_point[3]
            mid_point = points.GetPoint(mid_point_id)
            print("\tmidPoint is ", mid_point_id, " with coordinate ", "(",
                mid_point[0], ", ",
                mid_point[1], ", ",
                mid_point[2], ")",
                " and value ", scalars.GetTuple1(mid_point_id)
            )
            label_points.InsertNextPoint(mid_point)
            label_scalars.InsertNextTuple1(scalars.GetTuple1(mid_point_id))
            lineCount += 1

        label_poly_data.SetPoints(label_points)
        label_poly_data.GetPointData().SetScalars(label_scalars)

    contour_mapper = vtk.vtkPolyDataMapper()
    contour_mapper.SetInputConnection(contour_stripper.GetOutputPort())
    contour_mapper.ScalarVisibilityOff()

    isolines = vtk.vtkActor()
    isolines.SetMapper(contour_mapper)

    surfaceLUT = vtk.vtkLookupTable()
    surfaceLUT.SetRange(
      polyData.GetPointData().GetScalars().GetRange())
    surfaceLUT.Build()

    surface_mapper = vtk.vtkPolyDataMapper()
    if vtk.VTK_MAJOR_VERSION <= 5:
        surface_mapper.SetInput(polyData)
    else:
        surface_mapper.SetInputData(polyData)

    surface_mapper.ScalarVisibilityOn()
    surface_mapper.SetScalarRange(
      polyData.GetPointData().GetScalars().GetRange())
    surface_mapper.SetLookupTable(surfaceLUT)

    surface = vtk.vtkActor()
    surface.SetMapper(surface_mapper)

    # The labeled data mapper will place labels at the points
    if include_labels:
        label_mapper = vtk.vtkLabeledDataMapper()
        label_mapper.SetFieldDataName("Isovalues")
        if vtk.VTK_MAJOR_VERSION <= 5:
            label_mapper.SetInput(label_poly_data)
        else:
            label_mapper.SetInputData(label_poly_data)

        label_mapper.SetLabelModeToLabelScalars()
        label_mapper.SetLabelFormat("%6.2f")

        isolabels = vtk.vtkActor2D()
        isolabels.SetMapper(label_mapper)

    # Create a renderer and render window
    renderer = vtk.vtkRenderer()

    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)

    # Create an interactor
    render_window_interactor = vtk.vtkRenderWindowInteractor()
    render_window_interactor.SetRenderWindow(render_window)

    # Add the actors to the scene
    renderer.AddActor(isolines)

    if include_labels:
        renderer.AddActor(isolabels)
    #renderer.AddActor(surface)

    # Render the scene(lights and cameras are created automatically)
    render_window.Render()
    render_window_interactor.Start()

if __name__ == '__main__':
    main()
