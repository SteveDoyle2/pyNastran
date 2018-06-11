#!/usr/bin/env python
from __future__ import print_function
import vtk


def main():
    input_filename, number_of_cuts = get_program_parameters()

    colors = vtk.vtkNamedColors()

    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(input_filename)
    reader.Update()

    bounds = reader.GetOutput().GetBounds()
    print(bounds)

    plane = vtk.vtkPlane()
    plane.SetOrigin((bounds[1] + bounds[0]) / 2.0,
                    (bounds[3] + bounds[2]) / 2.0,
                    (bounds[5] + bounds[4]) / 2.0)
    plane.SetNormal(0., 0., 1.)

    # Create Scalars.
    scalars = vtk.vtkDoubleArray()
    number_of_points = reader.GetOutput().GetNumberOfPoints()
    scalars.SetNumberOfTuples(number_of_points)
    pts = reader.GetOutput().GetPoints()
    for i in range(number_of_points):
        point = pts.GetPoint(i)
        scalars.SetTuple1(i, plane.EvaluateFunction(point))
    reader.GetOutput().GetPointData().SetScalars(scalars)
    reader.GetOutput().GetPointData().GetScalars().GetRange()

    # Create the cutter.
    cutter = vtk.vtkContourFilter()
    cutter.SetInputConnection(reader.GetOutputPort())
    cutter.ComputeScalarsOff()
    cutter.ComputeNormalsOff()
    cutter.GenerateValues(
        number_of_cuts,
        0.99 * reader.GetOutput().GetPointData().GetScalars().GetRange()[0],
        0.99 * reader.GetOutput().GetPointData().GetScalars().GetRange()[1])

    cutter_mapper = vtk.vtkPolyDataMapper()
    cutter_mapper.SetInputConnection(cutter.GetOutputPort())
    cutter_mapper.ScalarVisibilityOff()

    # Create the cut actor.
    cutter_actor = vtk.vtkActor()
    cutter_actor.GetProperty().SetColor(colors.GetColor3d("Banana"))
    cutter_actor.GetProperty().SetLineWidth(2)
    cutter_actor.SetMapper(cutter_mapper)

    # Create the model actor
    model_mapper = vtk.vtkPolyDataMapper()
    model_mapper.SetInputConnection(reader.GetOutputPort())
    model_mapper.ScalarVisibilityOff()

    model_actor = vtk.vtkActor()
    model_actor.GetProperty().SetColor(colors.GetColor3d("Flesh"))
    model_actor.SetMapper(model_mapper)

    # Create renderers and add the plane and model actors.
    renderer = vtk.vtkRenderer()
    renderer.AddActor(cutter_actor)
    renderer.AddActor(model_actor)

    # Add renderer to renderwindow and render
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    render_window.SetSize(600, 600)

    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)

    renderer.SetBackground(colors.GetColor3d("Burlywood"))
    renderer.GetActiveCamera().SetPosition(0, -1, 0)
    renderer.GetActiveCamera().SetFocalPoint(0, 0, 0)
    renderer.GetActiveCamera().SetViewUp(0, 0, 1)
    renderer.GetActiveCamera().Azimuth(30)
    renderer.GetActiveCamera().Elevation(30)

    renderer.ResetCamera()
    render_window.Render()

    interactor.Start()

def get_program_parameters():
    import argparse
    description = 'Cutting a surface model of the skin with a series of planes produces contour lines.'
    epilogue = '''
    Lines are wrapped with tubes for visual clarity.
    '''
    parser = argparse.ArgumentParser(description=description, epilog=epilogue,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('filename1', help='Torso.vtp.')
    parser.add_argument('-n', type=int, default=20, help='Number of cuts.')
    args = parser.parse_args()
    return args.filename1, args.n


if __name__ == '__main__':
    main()
