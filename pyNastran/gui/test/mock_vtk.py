#import sys
from vtk import vtkCamera
#import vtk.util

#VTK_FLOAT = 1
#VTK_UNSIGNED_CHAR = 2
#VTK_UNSIGNED_SHORT = 3
#VTK_UNSIGNED_INT = 4
#VTK_VERSION = '7.1.1'
#class vtkLine:
    #def __init__(self):
        #pass
#class vtkTriangle:
    #def __init__(self):
        #pass
#class vtkQuad:
    #def __init__(self):
        #pass
#class vtkTetra:
    #def __init__(self):
        #pass
#class vtkWedge:
    #def __init__(self):
        #pass
#class vtkHexahedron:
    #def __init__(self):
        #pass
#class vtkQuadraticTriangle:
    #def __init__(self):
        #pass
#class vtkQuadraticQuad:
    #def __init__(self):
        #pass
#class vtkQuadraticTetra:
    #def __init__(self):
        #pass
#class vtkQuadraticWedge:
    #def __init__(self):
        #pass
#class vtkQuadraticHexahedron:
    #def __init__(self):
        #pass
#class vtkPyramid:
    #def __init__(self):
        #pass


#class vtkPoints:
    #def __init__(self):
        #self.npoints = 0
    #def SetNumberOfPoints(self, npoints):
        #assert isinstance(npoints, integer_types), 'npoints=%s type=%s' % (npoints, type(npoints))
        #self.npoints = npoints

#class Arrays:
    #def __init__(self):
        #pass
    #def AddArray(self, grid):
        #pass
    #def SetActiveScalars(self, name):
        #pass
    #def GetNumberOfArrays(self):
        #return 4
    #def GetArrayName(self, i):
        #return 'fakename'
    #def RemoveArray(self, name):
        #pass
    #def SetActiveVectors(self, name):
        #pass

#class vtkArray:
    #def __init__(self):
        #pass
    #def SetNumberOfComponents(self, ncomp):
        #assert isinstance(ncomp, int), ncomp
#class vtkLongArray(vtkArray):
    #def __init__(self):
        #Arrays.__init__(self)
    #def GetDataTypeSize(self):
        #return 8
#class vtkIdTypeArray:
    #def __init__(self):
        #pass
    #def GetDataTypeSize(self):
        #return 4


#class vtkDataArray:
    #def __init__(self):
        #pass
    #def CreateDataArray(self, vtk_array_type):
        #pass

#class vtkGenericRenderWindowInteractor:
    #def __init__(self):
        #pass
#class vtkInteractorStyleRubberBandZoom:
    #def __init__(self):
        #pass
#class vtkInteractorStyleTrackballCamera:
    #def __init__(self):
        #pass

#class vtkColorTransferFunction:
    #def __init__(self):
        #pass
    #def SetNanColor(self, red, green, blue):
        #assert isinstance(red, float), red
        #assert isinstance(green, float), green
        #assert isinstance(blue, float), blue
    #def SetColorSpaceToHSV(self):
        #pass
    ##def SetColorSpaceToRGB(self):
        ##pass
    #def AddRGBPoint(self, value, red, green, blue):
        #assert isinstance(value, float), value
        #assert isinstance(red, float), red
        #assert isinstance(green, float), green
        #assert isinstance(blue, float), blue
    #def HSVWrapOff(self):
        #pass
    #def SetRange(self, min_value, max_value):
        #assert isinstance(min_value, float), min_value
        #assert isinstance(max_value, float), max_value

#class vtkScalarBarActor:
    #def __init__(self):
        #pass
    #def SetTitle(self, title):
        #assert isinstance(title, str), 'title=%r' % title
    #def SetLookupTable(self, color_function):
        #assert isinstance(color_function, vtkColorTransferFunction), 'color_function=%r' % color_function
    #def SetOrientationToVertical(self):
        #pass
    #def SetPosition(self, x, y):
        #assert isinstance(x, float), x
        #assert isinstance(y, float), y
    #def SetHeight(self, height):
        #assert isinstance(height, float), height
    #def SetWidth(self, width):
        #assert isinstance(width, float), width
    #def SetLabelFormat(self, label_format):
        #assert isinstance(label_format, str), 'label_format=%r' % label_format
    #def SetNumberOfLabels(self, nlabels):
        #assert isinstance(nlabels, int), nlabels
    #def SetMaximumNumberOfColors(self, ncolors):
        #assert isinstance(ncolors, int), ncolors
    #def VisibilityOff(self):
        #pass
    #def Modified(self):
        #pass

class GeometryProperty:
    def __init__(self):
        pass
    def SetRepresentationToPoints(self):
        pass
    def SetColor(self, rgb_floats):
        for val in rgb_floats:
            assert isinstance(val, float), rgb_floats
    def SetBackfaceColor(self, rgb_floats):
        for val in rgb_floats:
            assert isinstance(val, float), rgb_floats

    def SetPointSize(self, size):
        assert isinstance(size, int), type(size)

class GridMapper:
    def __init__(self):
        pass
    def InterpolateScalarsBeforeMappingOff(self):
        pass
    def InterpolateScalarsBeforeMappingOn(self):
        pass


class vtkRenderer:
    def __init__(self):
        pass
    def AddActor(self, actor):
        pass
    def GetActiveCamera(self):
        return vtkCamera()
    def ResetCameraClippingRange(self):
        pass
    def ResetCamera(self):
        pass
    def Render(self):
        pass
    def SetBackground(self, a, b, c):
        pass
    def SetBackground2(self, a, b, c):
        pass
    def SetGradientBackground(self, a):
        pass

class vtkLODActor:
    def __init__(self):
        pass
    def SetVisibility(self, is_visible):
        pass

#class vtkCamera:
    #def __init__(self):
        #pass
    #def GetPosition(self):
        #return (0., 0., 0.)
    #def SetPosition(self, xyz):
        #pass
    #def GetFocalPoint(self):
        #return (0., 0., 0.)
    #def SetFocalPoint(self, xyz):
        #pass
    #def GetViewUp(self):
        #return (0., 0., 1.)

class VTKInteractor:
    def __init__(self):
        pass
    def Render(self):
        pass
    def GetRenderWindow(self):
        return vtkRenderer()
