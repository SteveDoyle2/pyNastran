import vtk
import math

vtk_major_version = vtk.VTK_MAJOR_VERSION
"""
vtk.vtkPoints(): A set of points cam be collected in a vtkPoints
object by first setting the number of points with the
vtk_points.SetNumberOfPoints(npts) method and then enumerating
through the points and using the vtk_points.SetPoint(i, point) method.
"""
npts = 100
vtk_points = vtk.vtkPoints()
vtk_points.SetNumberOfPoints(100)
for i in range (npts):
    x = math.sin(math.pi*i/20.)
    y = math.cos(math.pi*i/20.)
    z = 2*i/float(npts)
    vtk_points.SetPoint(i, (x,y,z))
"""
vtk.vtkCellArray(): the ordering of the points or the connectivity of
the points can be captured with a vtkCellArray object. This is done by
first declaring that you will be inserting a cell with a set number of
points using the vtk_cell_array.InsertNextCell(npts) and then looping
over the range of the number of points with the following loop,
for i in range(npts):
 vtk_cell_array.InsertCellPoint(i)
"""
vtk_cell_array = vtk.vtkCellArray()
vtk_cell_array.InsertNextCell(npts)
for i in range(npts):
    vtk_cell_array.InsertCellPoint(i)
"""
vtk.vtkFLoatArray(): A set of values associated with the initial
points, be it force magnitudes at the point or bending moments, can be
captured as a scalar array of floating point values with the
vtkFloatArray object. this is done by first declaring the number of
scalars with the vtk_float_array.SetNumberOfValues(npts) method and then
enumerating through the points and using the
vtk_float_array.SetValue(i, value) method.
"""
value = lambda i: math.fabs(math.sin(math.pi*i/30.))
vtk_float_array = vtk.vtkFloatArray()
vtk_float_array.SetNumberOfValues(npts)
for i in range(npts):
    vtk_float_array.SetValue(i, value(i))
"""
vtk.vtkPolyData(): The point connectivity and scalar values can be
encapsulated in the vtkPolyData object this is done by using the,
vtk_poly_data.SetPoints(vtk_points), vtk_poly_data.SetLines(vtk_cell_array) and
vtk_poly_data.GetPointData().SetScalars(vtkFLoatArray) methods
"""
vtk_poly_data = vtk.vtkPolyData()
vtk_poly_data.SetPoints(vtk_points)
vtk_poly_data.SetLines(vtk_cell_array)
vtk_poly_data.GetPointData().SetScalars(vtk_float_array)
"""
vtk.vtkSplineFilter(): The data can be smoothly interpolated across a
number number of subdivisions using the vtkSplineFilter object. This
is accomplished by first setting the input of the spline filter to be
the previously constructed vtkPolyData object using the
vtk_spline_filter.SetInput(vtk_poly_data) method, and then setting the
number of desired subdivisions to create a smooth curve with the
vtk_spline_filter.SetNumberOfSubdivisions(int), and finally call the
vtk_spline_filter.Update() method.
"""

vtk_spline_filter = vtk.vtkSplineFilter()
if vtk_major_version == 6:
    vtk_spline_filter.SetInputData(vtk_poly_data)
elif vtk_major_version in [7, 8]:
    vtk_spline_filter.SetInputData(vtk_poly_data)
else:
    # this was written in ~2009, so we had to hack it
    # http://www.patinaed.org/sdk/2009/09/vtk-vtktubefilter-example-in-python.html
    raise NotImplementedError(vtk_major_version)
    #vtk_spline_filter.SetInput(vtk_poly_data)
vtk_spline_filter.SetNumberOfSubdivisions(5*npts)
vtk_spline_filter.Update()

"""
vtk.vtkTubeFilter(): The data can be visualized as colored tube using
the vtkTubeFilter. This is accomplished by first setting the
vtkplineFilter output port as the input connection to the
vtkTubleFilter with the,
vtk_tube_filter.SetInputConnection(vtk_spline_filter.GetOutputPort())
statement. The tubes radius of the can be set with the
vtk_tube_filter.SetRadius(radius) method. The smoothness around the tube
can be adjusted with the vtk_tube_filter.SetNumberOfSides(nSides)
method. If it is desired to have the end caps on the tube this can be
done with the vtk_tube_filter.CappingOn() method.
"""
vtk_tube_filter = vtk.vtkTubeFilter()
vtk_tube_filter.SetInputConnection(vtk_spline_filter.GetOutputPort())
vtk_tube_filter.SetRadius(0.15)
vtk_tube_filter.SetNumberOfSides(20)
vtk_tube_filter.CappingOn()
"""
vtk.vtkPolyDataMapper(): The out put of the vtkTubeFilter can be
visulized by creating a vtkPolyDataMapper object and then setting the
vtkTubeFilters output port as the input connection of the
vtkPolyDataMapper with the
vtk_poly_data_mapper.SetInputConnection(vtk_tube_filter.GetOutputPort())
method.
"""
vtk_poly_data_mapper = vtk.vtkPolyDataMapper()
vtk_poly_data_mapper.SetInputConnection(vtk_tube_filter.GetOutputPort())

"""
vtk.vtkActor(): the vtkPolyDataMapper construction can be pipped into
a vtkActor object using the vtk_actor.SetMapper(vtk_poly_data_mapper)
method.
"""
vtk_actor = vtk.vtkActor()
vtk_actor.SetMapper(vtk_poly_data_mapper)
"""
vtk.vtkRenderer(): The vtkActor can be rendered for the screen
by using the vtkRenderer.AddActor(vtk_actor) method.
"""
vtk_renderer = vtk.vtkRenderer()
vtk_renderer.AddActor(vtk_actor)
"""
vtk.vtkRenderWindow(): the vtkRenderer can be added to a render window
with the vtk_render_window.AddRenderer(vtk_renderer) method. The final
rendering can be exicuted with the vtk_render_window.Render() method
"""
vtk_render_window = vtk.vtkRenderWindow()
vtk_render_window.AddRenderer(vtk_renderer)
vtk_render_window.Render()
"""
vtk.vtkRenderWindowInteractor(): The render window can be held on
the screen by adding set it with the
vtk_render_window_interactor.SetRenderWindow(vtk_render_window) method. Then you
must start and initialize the vtkRenderWindowInteractor with the
vtk_render_window_interactor.Initialize() and
vtk_render_window_interactor.Start() methods.
"""
vtk_render_window_interactor = vtk.vtkRenderWindowInteractor()
vtk_render_window_interactor.SetRenderWindow(vtk_render_window)
vtk_render_window_interactor.Initialize()
vtk_render_window_interactor.Start()
