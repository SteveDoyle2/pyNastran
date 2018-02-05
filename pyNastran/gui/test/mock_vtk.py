class MockArrays(object):
    def __init__(self):
        pass
    def AddArray(self, grid):
        pass
    def SetActiveScalars(self, name):
        pass
    def GetNumberOfArrays(self):
        return 4
    def GetArrayName(self, i):
        return 'fakename'
    def RemoveArray(self, name):
        pass
    def SetActiveVectors(self, name):
        pass

class MockGridMapper(object):
    def __init__(self):
        pass
    def InterpolateScalarsBeforeMappingOff(self):
        pass

class MockGeometryProperty(object):
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


class MockGeometryActor(object):
    def __init__(self):
        self._prop = MockGeometryProperty()
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


class MockGrid(object):
    def Reset(self):
        pass
    def Allocate(self, nelements, delta):
        pass
    def InsertNextCell(self, *cell):
        pass
    def SetPoints(self, *cell):
        pass
    def Modified(self):
        pass
    def Update(self):
        pass
    def SetCells(self, vtk_cell_types, vtk_cell_locations, vtk_cells):
        pass
    def GetCellData(self):
        return MockArrays()
    def GetPointData(self):
        return MockArrays()

class MockArrowSource(object):
    def __init__(self):
        pass

class MockGlyph3D(object):
    def SetScaleFactor(self, value):
        pass

class MockPolyDataMapper(object):
    def __init__(self):
        pass

class MockLODActor(object):
    def __init__(self):
        pass
    def SetVisibility(self, is_visible):
        pass

class MockVTKInteractor(object):
    def __init__(self):
        pass
    def Render(self):
        pass
