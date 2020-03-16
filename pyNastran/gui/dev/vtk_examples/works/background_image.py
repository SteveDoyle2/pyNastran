"""
based on:
 - http://www.vtk.org/Wiki/VTK/Examples/Cxx/Images/BackgroundImage

tested with:
 - VTK 6.3, Python 2.7
 - VTK 7.0, Python 3.5

>>> python background_image.py image_filename.jpg

"""
import sys
from vtk import (
    vtkJPEGReader, vtkImageCanvasSource2D, vtkImageActor, vtkPolyDataMapper,
    vtkRenderer, vtkRenderWindow, vtkRenderWindowInteractor, vtkSuperquadricSource,
    vtkActor, VTK_MAJOR_VERSION
)

def main(argv):
    #  Verify input arguments
    if len(argv) > 1:
        # Read the image
        jpeg_reader = vtkJPEGReader()
        if not jpeg_reader.CanReadFile(argv[1]):
            print("Error reading file %s" % argv[1])
            return

        jpeg_reader.SetFileName(argv[1])
        jpeg_reader.Update()
        image_data = jpeg_reader.GetOutput()
    else:
        canvas_source = vtkImageCanvasSource2D()
        canvas_source.SetExtent(0, 100, 0, 100, 0, 0)
        canvas_source.SetScalarTypeToUnsignedChar()
        canvas_source.SetNumberOfScalarComponents(3)
        canvas_source.SetDrawColor(127, 127, 100)
        canvas_source.FillBox(0, 100, 0, 100)
        canvas_source.SetDrawColor(100, 255, 255)
        canvas_source.FillTriangle(10, 10, 25, 10, 25, 25)
        canvas_source.SetDrawColor(255, 100, 255)
        canvas_source.FillTube(75, 75, 0, 75, 5.0)
        canvas_source.Update()
        image_data = canvas_source.GetOutput()

    # Create an image actor to display the image
    image_actor = vtkImageActor()

    if VTK_MAJOR_VERSION <= 5:
        image_actor.SetInput(image_data)
    else:
        image_actor.SetInputData(image_data)

    # Create a superquadric
    superquadric_source = vtkSuperquadricSource()
    superquadric_source.SetPhiRoundness(1.1)
    superquadric_source.SetThetaRoundness(.2)

    # Create a mapper and actor
    superquadric_mapper = vtkPolyDataMapper()
    superquadric_mapper.SetInputConnection(superquadric_source.GetOutputPort())

    superquadric_actor = vtkActor()
    superquadric_actor.SetMapper(superquadric_mapper)

    vtk_renderer = vtkRenderer()
    vtk_renderer.SetLayer(1)

    # Set up the render window and renderers such that there is
    # a background layer and a foreground layer
    background_renderer = vtkRenderer()
    background_renderer.SetLayer(0)
    background_renderer.InteractiveOff()
    background_renderer.AddActor(image_actor)

    render_window = vtkRenderWindow()
    render_window.SetNumberOfLayers(2)
    render_window.AddRenderer(background_renderer)
    render_window.AddRenderer(vtk_renderer)

    render_window_interactor = vtkRenderWindowInteractor()
    render_window_interactor.SetRenderWindow(render_window)

    # Add actors to the renderers
    vtk_renderer.AddActor(superquadric_actor)

    # Render once to figure out where the background camera will be
    render_window.Render()

    # Set up the background camera to fill the renderer with the image
    origin = image_data.GetOrigin()
    spacing = image_data.GetSpacing()
    extent = image_data.GetExtent()

    camera = background_renderer.GetActiveCamera()
    camera.ParallelProjectionOn()

    xc = origin[0] + 0.5*(extent[0] + extent[1]) * spacing[0]
    yc = origin[1] + 0.5*(extent[2] + extent[3]) * spacing[1]
    # xd = (extent[1] - extent[0] + 1) * spacing[0]
    yd = (extent[3] - extent[2] + 1) * spacing[1]
    d = camera.GetDistance()
    camera.SetParallelScale(0.5 * yd)
    camera.SetFocalPoint(xc, yc, 0.0)
    camera.SetPosition(xc, yc, d)

    # Render again to set the correct view
    render_window.Render()

    # Interact with the window
    render_window_interactor.Start()


if __name__ == '__main__':
    main(sys.argv)
