# kills the program when you hit Cntl+C from the command line
# doesn't save the current state as presumably there's been an error
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import vtk
from vtk import vtkActor2D

apply_colors = False

def main():
    inputFilename = 'element_labels_input.vtk'

    # read file
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(inputFilename)
    reader.ReadAllScalarsOn()
    reader.SetScalarsName(reader.GetScalarsNameInFile(0))
    reader.Update()

    ncells = reader.GetOutput().GetNumberOfCells()
    npoints = reader.GetOutput().GetNumberOfPoints()

    # get attributes
    ugrid = reader.GetOutput()
    cell_data = ugrid.GetCellData()
    point_data = ugrid.GetPointData()
    cell_scalars = cell_data.GetScalars(reader.GetScalarsNameInFile(0))

    # validate that attributes are read correctly
    print('name0:')
    for i in range(ncells):
        print(i, ": ", cell_scalars.GetComponent(i, 0))

    cell_scalars = cell_data.GetScalars(reader.GetScalarsNameInFile(1))
    print('\nname1:')
    for i in range(ncells):
        print(i, ": ", cell_scalars.GetComponent(i, 0))

    point_scalars = point_data.GetScalars(reader.GetScalarsNameInFile(0))
    print(point_scalars)
    print('\nname2:')
    for i in range(npoints):
        print(i, ": ", point_scalars.GetComponent(i, 0))
    #point_scalars = cell_data.GetScalars(reader.GetScalarsNameInFile(2))
    #print(point_scalars)

    # geometry filter
    geometry_filter = vtk.vtkUnstructuredGridGeometryFilter()
    geometry_filter.SetInputConnection(reader.GetOutputPort())
    geometry_filter.Update()

    # Generate data arrays containing point and cell ids
    ids = vtk.vtkIdFilter()
    ids.SetInputConnection(geometry_filter.GetOutputPort())
    ids.PointIdsOff()
    ids.CellIdsOff()
    ids.FieldDataOn()

    # Create labels for cells
    cell_centers = vtk.vtkCellCenters()
    cell_centers.SetInputConnection(ids.GetOutputPort())

    node_points = vtk.vtkVertexGlyphFilter()
    node_points.SetInputConnection(ids.GetOutputPort())
    node_points.Update()
    #node_points = vtk.vtkPointCentered()
    #node_points.SetInputConnection(ids.GetOutputPort())

    # lut
    lut = vtk.vtkLookupTable()
    lut.SetNumberOfTableValues(ncells)
    lut.Build()
    lut.SetTableValue(0, 1, 0, 0, 1) # red.
    lut.SetTableValue(1, 0, 1, 0, 1) # green.
    lut.SetTableValue(2, 0, 0, 1, 1) # blue.

    # mapper
    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputConnection(geometry_filter.GetOutputPort())
    mapper.SetScalarVisibility(1)
    if apply_colors:
        mapper.SetLookupTable(lut)
        mapper.SetScalarModeToUseCellData()
        mapper.SetScalarRange(11, 13)
    mapper.GetInput().GetCellData().SetActiveScalars("cell_tag")

    # label mapper
    cell_label_mapper = vtk.vtkLabeledDataMapper()
    cell_label_mapper.SetInputConnection(cell_centers.GetOutputPort())
    cell_label_mapper.SetLabelModeToLabelScalars()

    #point_label_mapper = vtk.vtkLabeledDataMapper()
    #point_label_mapper.SetInputConnection(node_points.GetOutputPort())
    #point_label_mapper.SetLabelModeToLabelScalars()

    # actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetRepresentationToWireframe()

    # label actor
    cell_label_actor = vtkActor2D()
    cell_label_actor.SetMapper(cell_label_mapper)

    #point_label_actor = vtk.vtkActor2D()
    #point_label_actor.SetMapper(point_label_mapper)

    # renderer
    from pyNastran.gui.vtk_rendering_core import vtkRenderer, vtkRenderWindow, vtkActor2D
    renderer = vtkRenderer()
    renderWindow = vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)
    renderer.AddActor(actor)
    renderer.AddActor(cell_label_actor)
    #renderer.AddActor(point_label_actor)
    renderWindow.Render()
    renderWindowInteractor.Start()

if __name__ == '__main__':   # pragma: no cover
    main()
