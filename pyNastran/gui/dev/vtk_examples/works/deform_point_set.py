import sys
import vtk
from pyNastran.gui.vtk_renering_core import (
    vtkRenderer, vtkRenderWindow, vtkRenderWindowInteractor,
    vtkActor,
    vtkPolyDataMapper)

def main():
    input_data = vtk.vtkPolyData()
    #double bounds[6]

    argc = 1
    if argc == 1:

        # Create a sphere to warp
        sphere = vtk.vtkSphereSource()
        sphere.SetThetaResolution(51)
        sphere.SetPhiResolution(17)
        sphere.Update()
        bounds = sphere.GetOutput().GetBounds()

        # Generate some scalars on the polydata
        ele = vtk.vtkElevationFilter()
        ele.SetInputConnection(sphere.GetOutputPort())
        ele.SetLowPoint(0, 0, -0.5)
        ele.SetHighPoint(0, 0, 0.5)
        ele.SetLowPoint((bounds[1] + bounds[0]) / 2.0,
                        (bounds[3] + bounds[2]) / 2.0,
                        -bounds[5])
        ele.SetHighPoint((bounds[1] + bounds[0]) / 2.0,
                         (bounds[3] + bounds[2]) / 2.0,
                         bounds[5])

        ele.Update()
        input_data.ShallowCopy(ele.GetOutput())

    else:
        input_filename = sys.argv[1]

        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(input_filename)
        reader.Update()

        input_data.ShallowCopy(reader.GetOutput())
        bounds = input_data.GetBounds()


    # Now create a control mesh, in this case a octagon that encloses
    # the point set

    points = vtk.vtkPoints()
    points.SetNumberOfPoints(6)
    points.SetPoint(0,
                    bounds[0] - .1 * (bounds[1] - bounds[0]),
                    (bounds[3] + bounds[2]) / 2.0,
                    (bounds[5] + bounds[4]) / 2.0)
    points.SetPoint(1,
                    bounds[1] + .1 * (bounds[1] - bounds[0]),
                    (bounds[3] + bounds[2]) / 2.0,
                    (bounds[5] + bounds[4]) / 2.0)
    points.SetPoint(2,
                    (bounds[1] + bounds[0]) / 2.0,
                    bounds[2] - .1 * (bounds[3] - bounds[2]),
                    (bounds[5] + bounds[4]) / 2.0)
    points.SetPoint(3,
                    (bounds[1] + bounds[0]) / 2.0,
                    bounds[3] + .1 * (bounds[3] - bounds[2]),
                    (bounds[5] + bounds[4]) / 2.0)
    points.SetPoint(4,
                    (bounds[1] + bounds[0]) / 2.0,
                    (bounds[3] + bounds[2]) / 2.0,
                    bounds[4] - .1 * (bounds[5] - bounds[4]))
    points.SetPoint(5,
                    (bounds[1] + bounds[0]) / 2.0,
                    (bounds[3] + bounds[2]) / 2.0,
                    bounds[5] + .1 * (bounds[5] - bounds[4]))

    tris = vtk.vtkCellArray()
    tris.InsertNextCell(3)
    tris.InsertCellPoint(2)
    tris.InsertCellPoint(0)
    tris.InsertCellPoint(4)

    tris.InsertNextCell(3)
    tris.InsertCellPoint(1)
    tris.InsertCellPoint(2)
    tris.InsertCellPoint(4)

    tris.InsertNextCell(3)
    tris.InsertCellPoint(3)
    tris.InsertCellPoint(1)
    tris.InsertCellPoint(4)

    tris.InsertNextCell(3)
    tris.InsertCellPoint(0)
    tris.InsertCellPoint(3)
    tris.InsertCellPoint(4)

    tris.InsertNextCell(3)
    tris.InsertCellPoint(0)
    tris.InsertCellPoint(2)
    tris.InsertCellPoint(5)

    tris.InsertNextCell(3)
    tris.InsertCellPoint(2)
    tris.InsertCellPoint(1)
    tris.InsertCellPoint(5)

    tris.InsertNextCell(3)
    tris.InsertCellPoint(1)
    tris.InsertCellPoint(3)
    tris.InsertCellPoint(5)

    tris.InsertNextCell(3)
    tris.InsertCellPoint(3)
    tris.InsertCellPoint(0)
    tris.InsertCellPoint(5)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetPolys(tris)

    # Display the control mesh
    mesh_mapper = vtkPolyDataMapper()
    mesh_mapper.SetInputData(polydata)
    mesh_actor = vtkActor()
    mesh_actor.SetMapper(mesh_mapper)
    mesh_actor.GetProperty().SetRepresentationToWireframe()
    mesh_actor.GetProperty().SetColor(0, 0, 0)

    # Do the initial weight generation
    deform = vtk.vtkDeformPointSet()
    deform.SetInputData(input_data)
    deform.SetControlMeshData(polydata)
    deform.Update() # this creates the initial weights

    # Now move one point and deform
    control_point = points.GetPoint(5)
    points.SetPoint(5,
                    control_point[0],
                    control_point[1],
                    bounds[5] + .8 * (bounds[5] - bounds[4]))
    points.Modified()

    # Display the warped polydata
    poly_mapper = vtkPolyDataMapper()
    poly_mapper.SetInputConnection(deform.GetOutputPort())
    poly_actor = vtkActor()
    poly_actor.SetMapper(poly_mapper)

    renderer = vtkRenderer()
    ren_win = vtkRenderWindow()
    ren_win.AddRenderer(renderer)
    iren = vtkRenderWindowInteractor()
    iren.SetRenderWindow(ren_win)

    renderer.AddActor(poly_actor)
    renderer.AddActor(mesh_actor)

    renderer.GetActiveCamera().SetPosition(1, 1, 1)
    renderer.ResetCamera()
    renderer.SetBackground(.2, .3, .4)

    ren_win.SetSize(300, 300)
    ren_win.Render()

    iren.Start()

if __name__ == '__main__':   # pragma: no cover
    main()
