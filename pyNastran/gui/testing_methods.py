from __future__ import print_function
from six import iteritems
from pyNastran.utils.log import get_logger

class GeometryActor(object):
    def VisibilityOn(self):
        pass
    def VisibilityOff(self):
        pass
    def Modified(self):
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


def add_dummy_gui_functions(test, debug=True):
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
    def passer1(a):
        pass
    def passer2(a, b):
        pass
    def create_alternate_vtk_grid(name):
        pass
    def _add_alt_actors(alt_grids):
        for name, grid in iteritems(alt_grids):
            test.geometry_actors[name] = GeometryActor()
    def log_info(msg):
        if debug:
            print('INFO:  ', msg)

    def log_error(msg):
        if debug:
            print('ERROR:  ', msg)


    test.form = []
    test.result_cases = {}
    test._finish_results_io = passer1
    test._finish_results_io2 = passer2
    test.geometry_actors = {
        'main' : GeometryActor(),
    }
    test.debug = True
    test.grid = Grid()
    test.grid2 = Grid()
    test.scalarBar = ScalarBar()
    test.alt_geometry_actor = ScalarBar()
    test.alt_grids = {
        'main' : test.grid,
        'caero' : Grid(),
        'caero_sub' : Grid(),
    }
    test.geometry_properties = {
        #'main' : None,
        #'caero' : None,
        #'caero_sub' : None,
    }
    test.create_alternate_vtk_grid = create_alternate_vtk_grid
    test._add_alt_actors = _add_alt_actors

    test.log = get_logger(log=None, level='debug')
    test.log_error = log_error
    #test.log_info = print
    test.log_info = log_info
    test.removeOldGeometry = removeOldGeometry
    test.cycleResults = cycleResults
    test.TurnTextOn = TurnTextOn
    test.TurnTextOff = TurnTextOff
    test.update_axes_length = update_axes_length
    test.cycleResults_explicit = passer

