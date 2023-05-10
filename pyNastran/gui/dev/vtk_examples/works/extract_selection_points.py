"""
converted from:
 - http://www.vtk.org/Wiki/VTK/Examples/Cxx/PolyData/extract_selectionCells
"""

import vtk
from pyNastran.gui.vtk_interface import vtkUnstructuredGrid

def main():
    point_source = vtk.vtkPointSource()
    point_source.SetNumberOfPoints(50)
    point_source.Update()

    print("There are %s input points\n" % point_source.GetOutput().GetNumberOfPoints())

    ids = vtk.vtkIdTypeArray()
    ids.SetNumberOfComponents(1)

    # Set values
    i = 10
    while i < 20:
        ids.InsertNextValue(i)
        i += 1

    selection_node = vtk.vtkSelectionNode()
    selection_node.SetFieldType(1) # POINT
    #  CELL_DATA = 0
    #  POINT_DATA = 1
    #  FIELD_DATA = 2
    #  VERTEX_DATA = 3
    #  EDGE_DATA = 4

    selection_node.SetContentType(4) # INDICES
    #SELECTIONS = 0
    #GLOBALIDS = 1
    #PEDIGREEIDS = 2
    #VALUES = 3
    #INDICES = 4
    #FRUSTUM = 5
    #LOCATIONS = 6
    #THRESHOLDS = 7
    #BLOCKS = 8
    selection_node.SetSelectionList(ids)

    selection = vtk.vtkSelection()
    selection.AddNode(selection_node)

    extract_selection = vtk.vtkExtractSelection()

    extract_selection.SetInputConnection(0, point_source.GetOutputPort())
    if vtk.VTK_MAJOR_VERSION <= 5:
        extract_selection.SetInput(1, selection)
    else:
        extract_selection.SetInputData(1, selection)
    extract_selection.Update()

    # In selection
    selected = vtkUnstructuredGrid()
    selected.ShallowCopy(extract_selection.GetOutput())

    print("There are %s points in the selection" % selected.GetNumberOfPoints())
    print("There are %s cells in the selection" % selected.GetNumberOfCells())

    # Get points that are NOT in the selection
    # invert the selection
    selection_node.GetProperties().Set(vtk.vtkSelectionNode.INVERSE(), 1)
    extract_selection.Update()

    not_selected = vtkUnstructuredGrid()
    not_selected.ShallowCopy(extract_selection.GetOutput())

    print("There are %s points NOT in the selection" % not_selected.GetNumberOfPoints())
    print("There are %s cells NOT in the selection" % not_selected.GetNumberOfCells())

    input_mapper = vtk.vtkDataSetMapper()
    input_mapper.SetInputConnection(point_source.GetOutputPort())
    input_actor = vtk.vtkActor()
    input_actor.SetMapper(input_mapper)

    selected_mapper = vtk.vtkDataSetMapper()
    if vtk.VTK_MAJOR_VERSION <= 5:
        selected_mapper.SetInputConnection(selected.GetProducerPort())
    else:
        selected_mapper.SetInputData(selected)
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
    left_viewport = [0.0, 0.0, 0.33, 1.0]
    center_viewport = [0.33, 0.0, .66, 1.0]
    right_viewport = [0.66, 0.0, 1.0, 1.0]

    # Setup the renderers
    left_renderer = vtk.vtkRenderer()
    render_window.AddRenderer(left_renderer)
    left_renderer.SetViewport(left_viewport)
    left_renderer.SetBackground(.6, .5, .4)

    center_renderer = vtk.vtkRenderer()
    render_window.AddRenderer(center_renderer)
    center_renderer.SetViewport(center_viewport)
    center_renderer.SetBackground(.3, .1, .4)

    right_renderer = vtk.vtkRenderer()
    render_window.AddRenderer(right_renderer)
    right_renderer.SetViewport(right_viewport)
    right_renderer.SetBackground(.4, .5, .6)

    left_renderer.AddActor(input_actor)
    center_renderer.AddActor(selected_actor)
    right_renderer.AddActor(not_selected_actor)

    left_renderer.ResetCamera()
    center_renderer.ResetCamera()
    right_renderer.ResetCamera()

    render_window.Render()
    interactor.Start()

if __name__ == '__main__':  # pragma: no cover
    main()
