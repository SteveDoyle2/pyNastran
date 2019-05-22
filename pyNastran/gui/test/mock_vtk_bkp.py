from pyNastran.gui.test.mock_vtk import GeometryProperty


#class Grid:
    #def Reset(self):
        #pass
    #def Allocate(self, nelements, delta):
        #pass
    #def InsertNextCell(self, *cell):
        #pass
    #def SetPoints(self, *cell):
        #pass
    #def Modified(self):
        #pass
    #def Update(self):
        #pass
    #def SetCells(self, vtk_cell_types, vtk_cell_locations, vtk_cells):
        #pass
    #def GetCellData(self):
        #return Arrays()
    #def GetPointData(self):
        #return Arrays()

#class vtkTextProperty:
    #def __init__(self):
        #pass
    #def SetFontFamilyToArial(self):
        #pass
    #def BoldOn(self):
        #pass
    #def BoldOff(self):
        #pass
    #def ShadowOn(self):
        #pass
    #def ShadowOff(self):
        #pass

#class vtkTextActor:
    #def __init__(self):
        #pass
    #def SetInput(self, string):
        #assert isinstance(string, str), 'type(string)=%s' % type(string)
    #def VisibilityOn(self):
        #pass
    #def VisibilityOff(self):
        #pass


class vtkActor:
    def __init__(self):
        self._prop = GeometryProperty()
    def GetProperty(self):
        return self._prop
    def SetBackfaceProperty(self, prop):
        pass
    def VisibilityOn(self):
        pass
    def VisibilityOff(self):
        pass
    def Modified(self):
        pass


class vtkArrowSource:
    def __init__(self):
        pass

class vtkGlyph3D:
    def SetScaleFactor(self, value):
        pass

class vtkPolyDataMapper:
    def __init__(self):
        pass

