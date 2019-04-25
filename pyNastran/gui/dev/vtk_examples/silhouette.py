"""per https://lorensen.github.io/VTKExamples/site/Cxx/PolyData/Silhouette/"""
import sys
import vtk
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

def main():
  nargs = len(sys.argv)

  if nargs < 2:
      sphere_source = vtk.vtkSphereSource()
      sphere_source.Update()
      polyData = sphere_source.GetOutput()
  else:
      reader = vtkXMLPolyDataReader()
      reader.SetFileName(sys.argv[1])

      clean = vtk.vtkCleanPolyData()
      clean.SetInputConnection(reader.GetOutputPort())
      clean.Update()
      poly_data = clean.GetOutput()

  colors = vtk.vtkNamedColors()

  #create mapper and actor for original model
  mapper = vtk.vtkPolyDataMapper()
  mapper.SetInputData(poly_data)
  mapper.ScalarVisibilityOff()

  actor = vtk.vtkActor()
  actor.SetMapper(mapper)
  actor.GetProperty().SetInterpolationToFlat()
  #actor.GetProperty().SetColor(colors.GetColor3d("Banana").GetData())
  actor.GetProperty().SetColor(1., 0., 1.)

  #create renderer and renderWindow
  renderer = vtk.vtkRenderer()
  render_window = vtk.vtkRenderWindow()
  render_window.AddRenderer(renderer)

  renderer.AddActor(actor) #view the original model

  #Compute the silhouette
  silhouette = vtk.vtkPolyDataSilhouette()
  silhouette.SetInputData(poly_data)
  silhouette.SetCamera(renderer.GetActiveCamera())
  silhouette.SetEnableFeatureAngle(0)

  #create mapper and actor for silouette
  mapper2 = vtk.vtkPolyDataMapper()
  mapper2.SetInputConnection(silhouette.GetOutputPort())

  actor2 = vtk.vtkActor()
  actor2.SetMapper(mapper2)
  #actor2.GetProperty().SetColor(colors.GetColor3d("Tomato").GetData())
  actor2.GetProperty().SetColor(1., 1., 0.)
  actor2.GetProperty().SetLineWidth(5)

  renderer.AddActor(actor2)
  #renderer.SetBackground(colors.GetColor3d("Silver").GetData())
  renderer.SetBackground(0.1, 0.2, 0.5)
  renderer.ResetCamera()
  renderer.GetActiveCamera().Azimuth(30)
  renderer.GetActiveCamera().Elevation(30)
  renderer.GetActiveCamera().Dolly(1.5)
  renderer.ResetCameraClippingRange()

  iren = vtk.vtkRenderWindowInteractor()
  iren.SetRenderWindow(render_window)

  #render and interact
  render_window.SetSize(640, 480)
  render_window.Render()
  iren.Start()

main()
