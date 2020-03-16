"""
based on:
 - http://www.vtk.org/Wiki/VTK/Examples/Python/PolyData/ExtractSelectionCells
"""

import vtk
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

def main():
    """
    ^y
    |
    8----9----10---11
    |  / |  / |  / |
    4----5----6----7
    |  / |  / |  / |
    0----1----2----3 --> x
    """
    ugrid = vtk.vtkUnstructuredGrid()
    xyzs = [
        [0., 0., 0.],
        [1., 0., 0.],
        [2., 0., 0.],
        [3., 0., 0.],

        [0., 1., 0.],
        [1., 1., 0.],
        [2., 1., 0.],
        [3., 1., 0.],

        [0., 2., 0.],
        [1., 2., 0.],
        [2., 2., 0.],
        [3., 2., 0.],
    ]
    # we make the lower triangle first, then the upper one to finish off the quad
    # go accross each row, left to right
    tris = [
        [0, 1, 5],
        [0, 5, 4],
        [1, 2, 6],
        [1, 6, 5],
        [2, 3, 7],
        [2, 7, 6],

        [4, 5, 9],
        [4, 9, 8],
        [5, 6, 10],
        [5, 10, 9],
        [6, 7, 11],
        [6, 11, 10],
    ]
    ids_to_show = [0, 1, #2, 3,
                   4, 5, 6, 7, 8, 9, 10, 11, 12]
    ids_to_show_updated = [
        0, 1, #2, 3,
        #4, 5,
        6, 7, 8, 9, 10, 11, 12]

    nnodes = len(xyzs)
    points = vtk.vtkPoints()
    points.SetNumberOfPoints(nnodes)

    nid = 0
    for i, xyz in enumerate(xyzs):
        points.InsertPoint(nid, *xyz)
        nid += 1

    for tri in tris:
        elem = vtk.vtkTriangle()
        (n1, n2, n3) = tri
        elem.GetPointIds().SetId(0, n1)
        elem.GetPointIds().SetId(1, n2)
        elem.GetPointIds().SetId(2, n3)
        ugrid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())

    ugrid.SetPoints(points)

    grid_mapper = vtk.vtkDataSetMapper()
    if vtk.VTK_MAJOR_VERSION <= 5:
        grid_mapper.SetInputConnection(ugrid.GetProducerPort())
    else:
        grid_mapper.SetInputData(ugrid)
    input_actor = vtk.vtkActor()
    input_actor.SetMapper(grid_mapper)

    render_window = vtk.vtkRenderWindow()
    camera = vtk.vtkCamera()
    interactor = vtk.vtkRenderWindowInteractor()

    interactor.SetRenderWindow(render_window)



    print("There are %s input points" % ugrid.GetNumberOfPoints())
    print("There are %s input cells" % ugrid.GetNumberOfCells())

    ids = vtk.vtkIdTypeArray()
    ids.SetNumberOfComponents(1)
    ids.Allocate(len(ids_to_show))

    for id_to_show in ids_to_show:
        ids.InsertNextValue(id_to_show)


    selection_node = vtk.vtkSelectionNode()
    selection_node.SetFieldType(vtk.vtkSelectionNode.CELL)
    selection_node.SetContentType(vtk.vtkSelectionNode.INDICES)
    selection_node.SetSelectionList(ids)


    selection = vtk.vtkSelection()
    selection.AddNode(selection_node)

    extract_selection = vtk.vtkExtractSelection()
    if vtk.VTK_MAJOR_VERSION <= 5:
        extract_selection.SetInput(0, ugrid)
        extract_selection.SetInput(1, selection)
    else:
        extract_selection.SetInputData(0, ugrid)
        extract_selection.SetInputData(1, selection)
    extract_selection.Update()

    # In selection
    grid_selected = vtk.vtkUnstructuredGrid()
    grid_selected.ShallowCopy(extract_selection.GetOutput())

    print("There are %s points in the selection" % grid_selected.GetNumberOfPoints())
    print("There are %s cells in the selection" % grid_selected.GetNumberOfCells())

    # Get points that are NOT in the selection
    # invert the selection
    selection_node.GetProperties().Set(vtk.vtkSelectionNode.INVERSE(), 1)
    extract_selection.Update()

    not_selected = vtk.vtkUnstructuredGrid()
    not_selected.ShallowCopy(extract_selection.GetOutput())

    print("There are %s points NOT in the selection" % not_selected.GetNumberOfPoints())
    print("There are %s cells NOT in the selection" % not_selected.GetNumberOfCells())


    selected_mapper = vtk.vtkDataSetMapper()
    if vtk.VTK_MAJOR_VERSION <= 5:
        selected_mapper.SetInputConnection(grid_selected.GetProducerPort())
    else:
        selected_mapper.SetInputData(grid_selected)

    selected_actor = vtk.vtkActor()
    selected_actor.SetMapper(selected_mapper)

    not_selected_mapper = vtk.vtkDataSetMapper()
    if vtk.VTK_MAJOR_VERSION <= 5:
        not_selected_mapper.SetInputConnection(not_selected.GetProducerPort())
    else:
        not_selected_mapper.SetInputData(not_selected)

    not_selected_actor = vtk.vtkActor()
    not_selected_actor.SetMapper(not_selected_mapper)

    # There will be one render window
    render_window = vtk.vtkRenderWindow()
    render_window.SetSize(900, 300)

    # And one interactor
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)

    # Define viewport ranges
    # (xmin, ymin, xmax, ymax)
    left_viewport = [0.0, 0.0, 0.5, 1.0]
    right_viewport = [0.5, 0.0, 1.0, 1.0]

    # Create a camera for all renderers
    camera = vtk.vtkCamera()

    # Setup the renderers
    left_renderer = vtk.vtkRenderer()
    render_window.AddRenderer(left_renderer)
    left_renderer.SetViewport(left_viewport)
    left_renderer.SetBackground(.6, .5, .4)
    left_renderer.SetActiveCamera(camera)

    right_renderer = vtk.vtkRenderer()
    render_window.AddRenderer(right_renderer)
    right_renderer.SetViewport(right_viewport)
    right_renderer.SetBackground(.3, .1, .4)
    right_renderer.SetActiveCamera(camera)
    right_renderer.AddActor(selected_actor)

    left_renderer.AddActor(input_actor)
    right_renderer.AddActor(not_selected_actor)

    right_renderer.ResetCamera()

    interactor.Start()
    #-----------------
    ids.Reset()
    ids.Allocate(len(ids_to_show_updated))
    for id_to_show in ids_to_show_updated:
        ids.InsertNextValue(id_to_show)
    ids.Modified()
    grid_selected.Modified()

    #ids.Update()
    ids.Modified()
    #selection_node.Update()
    selection_node.Modified()
    #selection.Update()
    selection.Modified()
    extract_selection.Update()
    extract_selection.Modified()
    #grid_selected.Update()
    grid_selected.Modified()
    selected_mapper.Update()
    selected_mapper.Modified()
    #selected_actor.Update()
    selected_actor.Modified()

    right_renderer.Modified()
    #right_renderer.Update()

    interactor.Modified()
    #interactor.Update()
    #-----------------
    render_window.Render()

    render_window.Modified()
    #render_window.Update()


if __name__ == '__main__':  # pragma: no cover
    main()
