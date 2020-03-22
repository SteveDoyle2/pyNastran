import vtk

def main():
    sphere_source = vtk.vtkSphereSource()
    sphere_source.SetCenter(0.0, 0.0, 0.0)
    sphere_source.SetRadius(5000.0)
    sphere_source.Update()

    polydata = sphere_source.GetOutput()

    # Create a mapper
    mapper = vtk.vtkPolyDataMapper()
    #if VTK_MAJOR_VERSION <= 5
    #mapper.SetInput(polydata)
    #else
    mapper.SetInputData(polydata)
    #endif

    # Create an actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    legend = vtk.vtkLegendBoxActor()
    legend.SetNumberOfEntries(2)

    colors = vtk.vtkNamedColors()

    legend_box = vtk.vtkCubeSource()
    legend_box.Update()

    color = [0., 0., 0., 0.]
    colors.GetColor('tomato', color)
    legend.SetEntry(0, legend_box.GetOutput(), 'Box', color[:3])

    color = [0., 0., 0., 0.]
    colors.GetColor('banana', color)
    legend.SetEntry(1, sphere_source.GetOutput(), 'Ball', color[:3])

    # place legend in lower right
    legend.GetPositionCoordinate().SetCoordinateSystemToView()
    legend.GetPositionCoordinate().SetValue(.5, -1.0)

    legend.GetPosition2Coordinate().SetCoordinateSystemToView()
    legend.GetPosition2Coordinate().SetValue(1.0, -0.5)

    legend.UseBackgroundOn()
    background = [0., 0., 0., 0.]
    colors.GetColor('warm_grey', background)
    legend.SetBackgroundColor(background[:3])

    # A renderer and render window
    renderer = vtk.vtkRenderer()
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)

    # Add the actors to the scene
    renderer.AddActor(actor)
    renderer.AddActor(legend)
    renderer.SetBackground(0,1,1) # Background color cyan

    # Render an image (lights and cameras are created automatically)
    render_window.Render()

    # An interactor
    render_window_interactor = vtk.vtkRenderWindowInteractor()
    render_window_interactor.SetRenderWindow(render_window)

    # Begin mouse interaction
    render_window_interactor.Start()

main()
