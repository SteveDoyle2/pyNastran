#include <vtkVersion.h
#include <vtkSmartPointer.h
#include <vtkPointData.h
#include <vtkIdTypeArray.h
#include <vtkDataSetSurfaceFilter.h
#include <vtkRendererCollection.h
#include <vtkProperty.h
#include <vtkPlanes.h
#include <vtkObjectFactory.h
#include <vtkPolyDataMapper.h
#include <vtkActor.h
#include <vtkRenderWindow.h
#include <vtkRenderer.h
#include <vtkRenderWindowInteractor.h
#include <vtkPolyData.h
#include <vtkPointSource.h
#include <vtkInteractorStyleRubberBandPick.h
#include <vtkAreaPicker.h
#include <vtkExtractGeometry.h
#include <vtkDataSetMapper.h
#include <vtkUnstructuredGrid.h
#include <vtkVertexGlyphFilter.h
#include <vtkIdFilter.h

import vtk

from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

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
        self.SelectedMapper = vtk.vtkDataSetMapper()
        self.SelectedActor = vtk.vtkActor()
        self.SelectedActor.SetMapper(SelectedMapper)

    def OnLeftButtonUp(self):
        # Forward events
        vtkInteractorStyleRubberBandPick.OnLeftButtonUp()

        frustum = vtk.vtkAreaPicker(self.GetInteractor().GetPicker()).GetFrustum()

        extractGeometry = vtk.vtkExtractGeometry.New()
        extractGeometry.SetImplicitFunction(frustum)

        if VTK_MAJOR_VERSION <= 5:
            extractGeometry.SetInput(this.Points)
        else:
            extractGeometry.SetInputData(this.Points)
        extractGeometry.Update()

        glyphFilter = vtk.vtkVertexGlyphFilter()
        glyphFilter.SetInputConnection(extractGeometry.GetOutputPort())
        glyphFilter.Update()

        selected = glyphFilter.GetOutput()
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

        self.SelectedActor.GetProperty().SetColor(1.0, 0.0, 0.0) #(R,G,B)
        self.SelectedActor.GetProperty().SetPointSize(3)

        self.CurrentRenderer.AddActor(SelectedActor)
        self.GetInteractor().GetRenderWindow().Render()
        #self.HighlightProp(NULL)

    def SetPoints(self, points):
        self.Points = points

    #private:
    #vtkPolyData Points
    #vtkActor SelectedActor
    #vtkDataSetMapper SelectedMapper


#vtkStandardNewMacro(InteractorStyle)

def main():
    pointSource = vtk.vtkPointSource()
    pointSource.SetNumberOfPoints(20)
    pointSource.Update()

    idFilter = vtk.vtkIdFilter()
    idFilter.SetInputConnection(pointSource.GetOutputPort())
    idFilter.SetIdsArrayName("OriginalIds")
    idFilter.Update()

    surfaceFilter = vtk.vtkDataSetSurfaceFilter()
    surfaceFilter.SetInputConnection(idFilter.GetOutputPort())
    surfaceFilter.Update()

    poly_input = surfaceFilter.GetOutput()

    # Create a mapper and actor
    mapper = vtk.vtkPolyDataMapper()

    if vtk.VTK_MAJOR_VERSION <= 5:
        mapper.SetInputConnection(poly_input.GetProducerPort())
    else:
        mapper.SetInputData(poly_input)
    mapper.ScalarVisibilityOff()

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    # Visualize
    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)

    areaPicker = vtk.vtkAreaPicker()
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetPicker(areaPicker)
    renderWindowInteractor.SetRenderWindow(renderWindow)

    renderer.AddActor(actor)
    #renderer.SetBackground(1,1,1) # Background color white

    renderWindow.Render()

    style = vtk.vtkRenderWindowInteractor()
    #style = myInteractorStyle()
    #style = InteractorStyle()
    #style = QVTKRenderWindowInteractor()
    #style.SetPoints(poly_input)
    renderWindowInteractor.SetInteractorStyle(style)

    renderWindowInteractor.Start()

if __name__ == '__main__':  # pragma: no cover
    main()
