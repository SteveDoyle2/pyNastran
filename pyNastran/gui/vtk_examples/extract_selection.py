"""
converted from:
 - http://www.vtk.org/Wiki/VTK/Examples/Cxx/PolyData/ExtractSelectionCells
"""

import vtk

def main():
    pointSource = vtk.vtkPointSource()
    pointSource.SetNumberOfPoints(50)
    pointSource.Update()

    print("There are %s input points\n" % pointSource.GetOutput().GetNumberOfPoints())

    ids = vtk.vtkIdTypeArray()
    ids.SetNumberOfComponents(1)

    # Set values
    i = 10
    while i < 20:
        ids.InsertNextValue(i)
        i += 1

    selectionNode = vtk.vtkSelectionNode()
    selectionNode.SetFieldType(1) # POINT
    #  CELL_DATA = 0
    #  POINT_DATA = 1
    #  FIELD_DATA = 2
    #  VERTEX_DATA = 3
    #  EDGE_DATA = 4

    selectionNode.SetContentType(4) # INDICES
    #SELECTIONS = 0
    #GLOBALIDS = 1
    #PEDIGREEIDS = 2
    #VALUES = 3
    #INDICES = 4
    #FRUSTUM = 5
    #LOCATIONS = 6
    #THRESHOLDS = 7
    #BLOCKS = 8
    selectionNode.SetSelectionList(ids)

    selection = vtk.vtkSelection()
    selection.AddNode(selectionNode)

    extractSelection = vtk.vtkExtractSelection()

    extractSelection.SetInputConnection(0, pointSource.GetOutputPort())
    if vtk.VTK_MAJOR_VERSION <= 5:
        extractSelection.SetInput(1, selection)
    else:
        extractSelection.SetInputData(1, selection)
    extractSelection.Update()

    # In selection
    selected = vtk.vtkUnstructuredGrid()
    selected.ShallowCopy(extractSelection.GetOutput())

    print("There are %s points in the selection" % selected.GetNumberOfPoints())
    print("There are %s cells in the selection" %selected.GetNumberOfCells())

    # Get points that are NOT in the selection
    # invert the selection
    selectionNode.GetProperties().Set(vtk.vtkSelectionNode.INVERSE(), 1)
    extractSelection.Update()

    notSelected = vtk.vtkUnstructuredGrid()
    notSelected.ShallowCopy(extractSelection.GetOutput())

    print("There are %s points NOT in the selection" % notSelected.GetNumberOfPoints())
    print("There are %s cells NOT in the selection" % notSelected.GetNumberOfCells())

    inputMapper = vtk.vtkDataSetMapper()
    inputMapper.SetInputConnection(pointSource.GetOutputPort())
    inputActor = vtk.vtkActor()
    inputActor.SetMapper(inputMapper)

    selectedMapper = vtk.vtkDataSetMapper()
    if vtk.VTK_MAJOR_VERSION <= 5:
        selectedMapper.SetInputConnection(selected.GetProducerPort())
    else:
        selectedMapper.SetInputData(selected)
    selectedActor = vtk.vtkActor()
    selectedActor.SetMapper(selectedMapper)

    notSelectedMapper = vtk.vtkDataSetMapper()
    if vtk.VTK_MAJOR_VERSION <= 5:
        notSelectedMapper.SetInputConnection(notSelected.GetProducerPort())
    else:
        notSelectedMapper.SetInputData(notSelected)
    notSelectedActor = vtk.vtkActor()
    notSelectedActor.SetMapper(notSelectedMapper)


    # There will be one render window
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.SetSize(900, 300)

    # And one interactor
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(renderWindow)

    # Define viewport ranges
    # (xmin, ymin, xmax, ymax)
    leftViewport = [0.0, 0.0, 0.33, 1.0]
    centerViewport = [0.33, 0.0, .66, 1.0]
    rightViewport = [0.66, 0.0, 1.0, 1.0]

    # Setup the renderers
    leftRenderer = vtk.vtkRenderer()
    renderWindow.AddRenderer(leftRenderer)
    leftRenderer.SetViewport(leftViewport)
    leftRenderer.SetBackground(.6, .5, .4)

    centerRenderer = vtk.vtkRenderer()
    renderWindow.AddRenderer(centerRenderer)
    centerRenderer.SetViewport(centerViewport)
    centerRenderer.SetBackground(.3, .1, .4)

    rightRenderer = vtk.vtkRenderer()
    renderWindow.AddRenderer(rightRenderer)
    rightRenderer.SetViewport(rightViewport)
    rightRenderer.SetBackground(.4, .5, .6)

    leftRenderer.AddActor(inputActor)
    centerRenderer.AddActor(selectedActor)
    rightRenderer.AddActor(notSelectedActor)

    leftRenderer.ResetCamera()
    centerRenderer.ResetCamera()
    rightRenderer.ResetCamera()

    renderWindow.Render()
    interactor.Start()

if __name__ == '__main__':
    main()
