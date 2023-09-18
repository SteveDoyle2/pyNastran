"""
converted from:
 - http://www.vtk.org/Wiki/VTK/Examples/Cxx/PolyData/ExtractSelectionCells
"""

import vtk
from pyNastran.gui.vtk_interface import vtkUnstructuredGrid

def main():
    sphere_source = vtk.vtkSphereSource()
    sphere_source.Update()

    print("There are %s input points" % sphere_source.GetOutput().GetNumberOfPoints())
    print("There are %s input cells" % sphere_source.GetOutput().GetNumberOfCells())

    ids = vtk.vtkIdTypeArray()
    ids.SetNumberOfComponents(1)

    # Specify that we want to extract cells 10 through 19
    i = 10
    while i < 20:
        ids.InsertNextValue(i)
        i += 1

    selection_node = vtk.vtkSelectionNode()
    selection_node.SetFieldType(vtk.vtkSelectionNode.CELL)
    selection_node.SetContentType(vtk.vtkSelectionNode.INDICES)
    selection_node.SetSelectionList(ids)


    selection = vtk.vtkSelection()
    selection.AddNode(selection_node)

    extract_selection = vtk.vtkExtractSelection()
    if vtk.VTK_MAJOR_VERSION <= 5:
        extract_selection.SetInputConnection(0, sphere_source.GetOutputPort())
        extract_selection.SetInput(1, selection)
    else:
        extract_selection.SetInputConnection(0, sphere_source.GetOutputPort())
        extract_selection.SetInputData(1, selection)
    extract_selection.Update()

    # In selection
    grid_selected = vtkUnstructuredGrid()
    grid_selected.ShallowCopy(extract_selection.GetOutput())

    print("There are %s points in the selection" % grid_selected.GetNumberOfPoints())
    print("There are %s cells in the selection" % grid_selected.GetNumberOfCells())

    # Get points that are NOT in the selection
    # invert the selection
    selection_node.GetProperties().Set(vtk.vtkSelectionNode.INVERSE(), 1)
    extract_selection.Update()

    not_selected = vtkUnstructuredGrid()
    not_selected.ShallowCopy(extract_selection.GetOutput())

    print("There are %s points NOT in the selection" % not_selected.GetNumberOfPoints())
    print("There are %s cells NOT in the selection" % not_selected.GetNumberOfCells())


    backfaces = vtk.vtkProperty()
    backfaces.SetColor(1, 0, 0)

    input_mapper = vtk.vtkDataSetMapper()
    input_mapper.SetInputConnection(sphere_source.GetOutputPort())
    input_actor = vtk.vtkActor()
    input_actor.SetMapper(input_mapper)
    input_actor.SetBackfaceProperty(backfaces)

    selected_mapper = vtk.vtkDataSetMapper()
    if vtk.VTK_MAJOR_VERSION <= 5:
        selected_mapper.SetInputConnection(grid_selected.GetProducerPort())
    else:
        selected_mapper.SetInputData(grid_selected)

    selected_actor = vtk.vtkActor()
    selected_actor.SetMapper(selected_mapper)
    selected_actor.SetBackfaceProperty(backfaces)

    not_selected_mapper = vtk.vtkDataSetMapper()
    if vtk.VTK_MAJOR_VERSION <= 5:
        not_selected_mapper.SetInputConnection(not_selected.GetProducerPort())
    else:
        not_selected_mapper.SetInputData(not_selected)

    not_selected_actor = vtk.vtkActor()
    not_selected_actor.SetMapper(not_selected_mapper)
    not_selected_actor.SetBackfaceProperty(backfaces)

    # There will be one render window
    render_window = vtk.vtkRenderWindow()
    render_window.SetSize(900, 300)

    # And one interactor
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)

    # Define viewport ranges
    # (xmin, ymin, xmax, ymax)
    left_viewport = [0.0, 0.0, 0.5, 1.0]
    center_viewport = [0.33, 0.0, .66, 1.0]
    right_viewport = [0.66, 0.0, 1.0, 1.0]

    # Create a camera for all renderers
    camera = vtk.vtkCamera()

    # Setup the renderers
    left_renderer = vtk.vtkRenderer()
    render_window.AddRenderer(left_renderer)
    left_renderer.SetViewport(left_viewport)
    left_renderer.SetBackground(.6, .5, .4)
    left_renderer.SetActiveCamera(camera)

    center_renderer = vtk.vtkRenderer()
    render_window.AddRenderer(center_renderer)
    center_renderer.SetViewport(center_viewport)
    center_renderer.SetBackground(.3, .1, .4)
    center_renderer.SetActiveCamera(camera)

    right_renderer = vtk.vtkRenderer()
    render_window.AddRenderer(right_renderer)
    right_renderer.SetViewport(right_viewport)
    right_renderer.SetBackground(.4, .5, .6)
    right_renderer.SetActiveCamera(camera)

    left_renderer.AddActor(input_actor)
    center_renderer.AddActor(selected_actor)
    right_renderer.AddActor(not_selected_actor)

    left_renderer.ResetCamera()

    render_window.Render()
    interactor.Start()

if __name__ == '__main__':  # pragma: no cover
    main()
