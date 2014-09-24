from __future__ import print_function
from pyNastran.utils.log import get_logger

def add_dummy_gui_functions(test):
    def removeOldGeometry(self):
        pass
    def cycleResults(self):
        pass
    def TurnTextOn():
        pass
    def TurnTextOff():
        pass
    def update_axes_length(value):
        pass
    def passer():
        pass

    class Grid(object):
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

    class ScalarBar(object):
        def VisibilityOff(self):
            pass
        def VisibilityOn(self):
            pass
        def Modified(self):
            pass

    test.grid = Grid()
    test.grid2 = Grid()
    test.scalarBar = ScalarBar()
    test.log = get_logger(log=None, level='debug')
    test.log_info = print
    test.removeOldGeometry = removeOldGeometry
    test.cycleResults = cycleResults
    test.TurnTextOn = TurnTextOn
    test.TurnTextOff = TurnTextOff
    test.update_axes_length = update_axes_length
    test.cycleResults_explicit = passer

