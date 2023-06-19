import os
import vtk
from pyNastran.gui.vtk_renering_core import (
    vtkRenderer, vtkRenderWindow, vtkRenderWindowInteractor,
    vtkActor, vtkActor2D)
#from pyNastran.gui.vtk_interface import vtkUnstructuredGrid

def main():
    inputFilename = '3_cells.vtk'
    assert os.path.exists(inputFilename)

    # read file.
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(inputFilename)
    reader.ReadAllScalarsOn()
    reader.SetScalarsName(reader.GetScalarsNameInFile(0))
    reader.Update()

    ncell = reader.GetOutput().GetNumberOfCells()
    assert ncell > 0

    # get attributes.
    ugrid = reader.GetOutput()
    cellData = ugrid.GetCellData()
    data = cellData.GetScalars(reader.GetScalarsNameInFile(0))

    # validate that attributes are read correctly.
    #for (int i=0 i<ncell i++)
    #{
        #std::cout<< i << ": " << data.GetComponent(i,0)<< std::endl
    #}
    data = cellData.GetScalars(reader.GetScalarsNameInFile(1))
    #for (int i=0 i<ncell i++)
    #{
        #std::cout<< i << ": " << data.GetComponent(i,0)<< std::endl
    #}

    data = cellData.GetScalars(reader.GetScalarsNameInFile(0))

    # geometry filter.
    geometryFilter = vtk.vtkUnstructuredGridGeometryFilter()
    geometryFilter.SetInputConnection(reader.GetOutputPort())
    geometryFilter.Update()

    # Generate data arrays containing point and cell ids
    ids = vtk.vtkIdFilter()
    ids.SetInputConnection(geometryFilter.GetOutputPort())
    ids.PointIdsOff()
    ids.CellIdsOff()
    ids.FieldDataOn()

    # Create labels for cells
    cc = vtk.vtkCellCenters()
    cc.SetInputConnection(ids.GetOutputPort())

    # lut
    lut = vtk.vtkLookupTable()
    lut.SetNumberOfTableValues(ncell)
    lut.Build()
    lut.SetTableValue(0, 1, 0, 0, 1) # red.
    lut.SetTableValue(1, 0, 1, 0, 1) # green.
    lut.SetTableValue(2, 0, 0, 1, 1) # blue.

    # mapper.
    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputConnection(geometryFilter.GetOutputPort())
    mapper.SetLookupTable(lut)
    mapper.SetScalarVisibility(1)
    mapper.SetScalarModeToUseCellData()
    mapper.SetScalarRange(11, 13)
    mapper.GetInput().GetCellData().SetActiveScalars("cell_tag")

    # label mapper.
    label_mapper = vtk.vtkLabeledDataMapper()
    label_mapper.SetInputConnection(cc.GetOutputPort())
    label_mapper.SetLabelModeToLabelScalars()

    # actor.
    actor = vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetRepresentationToWireframe()

    # label actor.
    label_actor = vtkActor2D()
    label_actor.SetMapper(label_mapper)

    # renderer.
    renderer = vtkRenderer()
    renderWindow = vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindowInteractor = vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)
    renderer.AddActor(actor)
    renderer.AddActor(label_actor)
    renderWindow.Render()
    renderWindowInteractor.Start()

main()
