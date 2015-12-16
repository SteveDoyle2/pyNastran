#include "vtkSmartPointer.h"

#include "vtkActor.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkRenderer.h"
#include "vtkSphereSource.h"
#include "vtkElevationFilter.h"
#include "vtkProperty.h"
#include "vtkCellArray.h"
#include "vtkPoints.h"
#include "vtkDeformPointSet.h"
#include "vtkCamera.h"
#include "vtkXMLPolyDataReader.h"

import vtk

def main():
    input = vtk.vtkPolyData()
    #double bounds[6]

    argc = 1
    if(argc == 1):

        # Create a sphere to warp
        sphere = vtk.vtkSphereSource()
        sphere.SetThetaResolution(51)
        sphere.SetPhiResolution(17)
        sphere.Update()
        bounds = sphere.GetOutput().GetBounds()

        # Generate some scalars on the polydata
        ele = vtk.vtkElevationFilter()
        ele.SetInputConnection(sphere.GetOutputPort())
        ele.SetLowPoint(0,0,-0.5)
        ele.SetHighPoint(0,0,0.5)
        ele.SetLowPoint((bounds[1] + bounds[0]) / 2.0,
                         (bounds[3] + bounds[2]) / 2.0,
                         -bounds[5])
        ele.SetHighPoint((bounds[1] + bounds[0]) / 2.0,
                          (bounds[3] + bounds[2]) / 2.0,
                          bounds[5])

        ele.Update()
        input.ShallowCopy(ele.GetOutput())

    else:
        inputFilename = argv[1]

        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(inputFilename)
        reader.Update()

        input.ShallowCopy(reader.GetOutput())
        bounds = input.GetBounds()


    # Now create a control mesh, in this case a octagon that encloses
    # the point set

    pts = vtk.vtkPoints()
    pts.SetNumberOfPoints(6)
    pts.SetPoint(0,
                  bounds[0] - .1 * (bounds[1] - bounds[0]),
                  (bounds[3] + bounds[2]) / 2.0,
                  (bounds[5] + bounds[4]) / 2.0)
    pts.SetPoint(1,
                  bounds[1] + .1 * (bounds[1] - bounds[0]),
                  (bounds[3] + bounds[2]) / 2.0,
                  (bounds[5] + bounds[4]) / 2.0)
    pts.SetPoint(2,
                  (bounds[1] + bounds[0]) / 2.0,
                  bounds[2] - .1 * (bounds[3] - bounds[2]),
                  (bounds[5] + bounds[4]) / 2.0)
    pts.SetPoint(3,
                  (bounds[1] + bounds[0]) / 2.0,
                  bounds[3] + .1 * (bounds[3] - bounds[2]),
                  (bounds[5] + bounds[4]) / 2.0)
    pts.SetPoint(4,
                  (bounds[1] + bounds[0]) / 2.0,
                  (bounds[3] + bounds[2]) / 2.0,
                  bounds[4] - .1 * (bounds[5] - bounds[4]))
    pts.SetPoint(5,
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

    pd = vtk.vtkPolyData()
    pd.SetPoints(pts)
    pd.SetPolys(tris)

    # Display the control mesh
    meshMapper = vtk.vtkPolyDataMapper()
    meshMapper.SetInputData(pd)
    meshActor = vtk.vtkActor()
    meshActor.SetMapper(meshMapper)
    meshActor.GetProperty().SetRepresentationToWireframe()
    meshActor.GetProperty().SetColor(0,0,0)

    # Do the intitial weight generation
    deform = vtk.vtkDeformPointSet()
    deform.SetInputData(input)
    deform.SetControlMeshData(pd)
    deform.Update() # this creates the initial weights

    # Now move one point and deform
    #double controlPoint[3]
    controlPoint = pts.GetPoint(5)
    pts.SetPoint(5, controlPoint[0],
                  controlPoint[1],
                  bounds[5] + .8 * (bounds[5] - bounds[4]))
    pts.Modified()

    # Display the warped polydata
    polyMapper = vtk.vtkPolyDataMapper()
    polyMapper.SetInputConnection(deform.GetOutputPort())
    polyActor = vtk.vtkActor()
    polyActor.SetMapper(polyMapper)

    renderer = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(renderer)
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)

    renderer.AddActor(polyActor)
    renderer.AddActor(meshActor)

    renderer.GetActiveCamera().SetPosition(1,1,1)
    renderer.ResetCamera()
    renderer.SetBackground(.2, .3, .4)

    renWin.SetSize(300,300)
    renWin.Render()

    iren.Start()

    return

main()