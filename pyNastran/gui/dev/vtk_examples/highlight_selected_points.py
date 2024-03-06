import vtk
from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from pyNastran.gui.vtk_rendering_core import vtkRenderer, vtkRenderWindow, vtkRenderWindowInteractor, vtkActor

class myInteractorStyle(vtk.vtkInteractorStyleRubberBandPick):
    #def __init__(self, parent=None
                 #):
        #super(QtGui.QDockWidget, self).__init__('Python Console', parent=parent)
        #vtk.vtkInteractorStyleRubberBandPick.__init__(self)

    def SetPoints(self, points):
        self.Points = points

# Define interaction style
class InteractorStyle(vtk.vtkInteractorStyleRubberBandPick):
    def __init__(self):
        self.__init__(vtk.vtkInteractorStyleRubberBandPick)
    #public:
    #static InteractorStyle* New()
    #vtkTypeMacro(InteractorStyle,vtkInteractorStyleRubberBandPick)

    def InteractorStyle(self):
        self.selected_mapper = vtk.vtkDataSetMapper()
        self.selected_actor = vtkActor()
        self.selected_actor.SetMapper(selected_mapper)

    def OnLeftButtonUp(self):
        # Forward events
        vtk.vtkInteractorStyleRubberBandPick.OnLeftButtonUp()

        frustum = vtk.vtkAreaPicker(self.GetInteractor().GetPicker()).GetFrustum()

        extract_geometry = vtk.vtkExtractGeometry()
        extract_geometry.SetImplicitFunction(frustum)

        if vtk.VTK_MAJOR_VERSION <= 5:
            extract_geometry.SetInput(self.Points)
         else:
            extract_geometry.SetInputData(self.Points)
        extract_geometry.Update()

        glyph_filter = vtk.vtkVertexGlyphFilter()
        glyph_filter.SetInputConnection(extract_geometry.GetOutputPort())
        glyph_filter.Update()

        selected = glyph_filter.GetOutput()
        print("Selected %s points" % selected.GetNumberOfPoints())
        print("Selected %s cells" % selected.GetNumberOfCells())
        if vtk.VTK_MAJOR_VERSION <= 5:
            self.SelectedMapper.SetInput(selected)
        else:
            self.SelectedMapper.SetInputData(selected)

        self.SelectedMapper.ScalarVisibilityOff()

        ids = vtkIdTypeArray.SafeDownCast(selected.GetPointData().GetArray("OriginalIds"))
        for i in range(ids.GetNumberOfTuples()):
            print("Id %s : %s" % (i, ids.GetValue(i)))

        prop = self.SelectedActor.GetProperty()
        prop.SetColor(1.0, 0.0, 0.0) #(R,G,B)
        prop.SetPointSize(3)

        self.CurrentRenderer.AddActor(selected_actor)
        self.GetInteractor().GetRenderWindow().Render()
        #self.HighlightProp(NULL)

    def SetPoints(self, points):
        self.Points = points

    #private:
    #vtkPolyData Points
    #vtkActor SelectedActor
    #vtkDataSetMapper Selected_mapper


#vtkStandardNewMacro(InteractorStyle)

#def main():
if 1:
    point_source = vtk.vtkPointSource()
    point_source.SetNumberOfPoints(20)
    point_source.Update()

    id_filter = vtk.vtkIdFilter()
    id_filter.SetInputConnection(point_source.GetOutputPort())
    id_filter.SetIdsArrayName("OriginalIds")
    id_filter.Update()

    surface_filter = vtk.vtkDataSetSurfaceFilter()
    surface_filter.SetInputConnection(id_filter.GetOutputPort())
    surface_filter.Update()

    poly_input = surface_filter.GetOutput()

    # Create a mapper and actor
    mapper = vtk.vtkPolyDataMapper()

    if vtk.VTK_MAJOR_VERSION <= 5:
        mapper.SetInputConnection(poly_input.GetProducerPort())
    else:
        mapper.SetInputData(poly_input)
    mapper.ScalarVisibilityOff()

    actor = vtkActor()
    actor.SetMapper(mapper)

    # Visualize
    renderer = vtkRenderer()
    render_window = vtkRenderWindow()
    render_window.AddRenderer(renderer)

    area_picker = vtk.vtkAreaPicker()
    render_window_interactor = vtk.vtkRenderWindowInteractor()
    render_window_interactor.SetPicker(area_picker)
    render_window_interactor.SetRenderWindow(render_window)

    renderer.AddActor(actor)
    #renderer.SetBackground(1,1,1) # Background color white

    render_window.Render()

    style = vtkRenderWindowInteractor()
    #style = myInteractorStyle()
    #style = InteractorStyle()
    #style = QVTKRenderWindowInteractor()
    #style.SetPoints(poly_input)
    render_window_interactor.SetInteractorStyle(style)

    render_window_interactor.Start()

#if __name__ == '__main__':  # pragma: no cover
#    main()
