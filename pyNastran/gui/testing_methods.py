from __future__ import print_function
#import os
from collections import OrderedDict

from six import iteritems, integer_types

from pyNastran.gui.qt_files.gui_qt_common import GuiCommon
from pyNastran.gui.qt_files.scalar_bar import ScalarBar
from pyNastran.gui.gui_objects.alt_geometry_storage import AltGeometry
from pyNastran.gui.formats import CLASS_MAP

from pyNastran.bdf.cards.base_card import deprecated
from pyNastran.utils.log import get_logger

from pyNastran.gui.test.mock_vtk import (
    GeometryProperty, GridMapper, # Grid, vtkActor,
    ArrowSource, Glyph3D, PolyDataMapper, VTKInteractor, vtkRenderer,
)
from vtk import vtkTextActor, vtkLODActor, vtkActor
import vtk

#class ScalarBar(object):
    #def VisibilityOff(self):
        #pass
    #def VisibilityOn(self):
        #pass
    #def Modified(self):
        #pass
    #@property
    #def is_shown(self):
        #return True


class Button(object):
    def __init__(self):
        pass
    def setChecked(self, is_checked):
        pass

class MockTreeView(object):
    def __init__(self):
        self.fringe = Button()
        self.vector = Button()
        self.disp = Button()

class MockResultCaseWindow(object):
    def __init__(self):
        self.treeView = MockTreeView()

class MockResWidget(object):
    def __init__(self):
        self.result_case_window = MockResultCaseWindow()
    def update_results(self, form, name):
        """fake method"""
        pass
    def update_method(self, methods):
        """fake method"""
        pass

class FakeGUIMethods(GuiCommon):
    """all the methods in here are faked"""
    def __init__(self, inputs=None):
        if inputs is None:
            inputs = {
                'magnify' : 1,
                'debug' : False,
                'console' : True,
                'is_groups' : False,
            }
        #GuiCommon.__init__(self, inputs=inputs)

        res_widget = MockResWidget()
        kwds = {
            'inputs' : inputs,
            'res_widget' : res_widget
        }
        #GuiAttributes.__init__(self, **kwds)
        GuiCommon.__init__(self, **kwds)
        self.res_widget = res_widget
        self.vtk_interactor = VTKInteractor()
        self.debug = False
        self._form = []
        self.result_cases = OrderedDict()
        #self.geometry_actors = {
            #'main' : vtkActor(),
        #}
        self.main_grid_mappers = {'main' : GridMapper()}
        self.grid = vtk.vtkUnstructuredGrid()
        #self.scalarBar = ScalarBar()
        self.scalar_bar = ScalarBar()
        self.alt_geometry_actor = ScalarBar()
        self.alt_grids = {
            'main' : self.grid,
        }
        self.main_geometry_actors = {
            'main' :  vtkActor(),
        }

        self.glyph_source = ArrowSource()
        self.glyphs = Glyph3D()
        self.glyph_mapper = PolyDataMapper()
        self.arrow_actor = vtkLODActor()
        self.arrow_actor_centroid = vtkLODActor()

        #self.geometry_properties = {
            #'main' : None,
            #'caero' : None,
            #'caero_sub' : None,
        #}
        #self._add_alt_actors = _add_alt_actors

        level = 'debug' if self.debug else 'info'
        self.log = get_logger(log=None, level=level)
        self.rend = vtkRenderer()

        self.text_actors[0] = vtkTextActor()
        self.text_actors[1] = vtkTextActor()
        self.text_actors[2] = vtkTextActor()
        self.text_actors[3] = vtkTextActor()
        self.format_class_map = CLASS_MAP

    def setup_fake_text_actors(self):
        for icase in self.result_cases:
            self.label_actors[icase] = []

    @property
    def scalarBar(self):
        return self.scalar_bar.scalar_bar

    @property
    def grid_selected(self):
        return self.grid

    #def hide_legend(self):
        #pass
    #def show_legend(self):
        #pass

    #def update_scalar_bar(self, title, min_value, max_value, norm_value,
                          #data_format,
                          #nlabels=None, labelsize=None,
                          #ncolors=None, colormap='jet',
                          #is_shown=True):
        #pass

    def update_legend(self, icase_fringe, icase_disp, icase_vector,
                      name, min_value, max_value, data_format, scale, phase,
                      arrow_scale,
                      nlabels, labelsize, ncolors, colormap,
                      use_fringe_internal=False, use_disp_internal=False,
                      use_vector_internal=False, external_call=True):
        pass

    def update_menu_bar(self):
        pass

    def _finish_results_io2(self, form, cases, reset_labels=True):
        """
        This is not quite the same as the main one.
        It's more or less just _set_results
        """
        if self.node_ids is None:  # pragma: no cover
            raise RuntimeError('implement self.node_ids for this format')
        if self.element_ids is None:  # pragma: no cover
            raise RuntimeError('implement self.element_ids for this format')
        #assert hasattr(self, 'gui'), 'gui does not exist for this format'
        assert hasattr(self, 'isubcase_name_map'), 'isubcase_name_map does not exist for this format'
        assert isinstance(self.nnodes, integer_types), 'nnodes=%r must be an integer' % self.nnodes
        assert isinstance(self.nelements, integer_types), 'nelements=%r must be an integer' % self.nelements

        assert len(cases) > 0, cases
        if isinstance(cases, OrderedDict):
            self.case_keys = list(cases.keys())
        else:
            self.case_keys = sorted(cases.keys())
            assert isinstance(cases, dict), type(cases)

        #print('self.case_keys = ', self.case_keys)
        for key in self.case_keys:
            assert isinstance(key, integer_types), key
            obj, (i, name) = cases[key]
            value = cases[key]
            if isinstance(value[0], int):
                raise RuntimeError('old style key is being used.\n key=%s\n type=%s value=%s' % (
                    key, type(value[0]), value))
            #assert len(value) == 2, 'value=%s; len=%s' % (str(value), len(value))

            subcase_id = obj.subcase_id
            case = obj.get_result(i, name)
            unused_result_type = obj.get_title(i, name)
            vector_size = obj.get_vector_size(i, name)
            #location = obj.get_location(i, name)
            methods = obj.get_methods(i)
            data_format = obj.get_data_format(i, name)
            scale = obj.get_scale(i, name)
            phase = obj.get_phase(i, name)
            unused_label2 = obj.get_header(i, name)
            unused_flag = obj.is_normal_result(i, name)
            #scalar_result = obj.get_scalar(i, name)
            nlabels, labelsize, ncolors, colormap = obj.get_nlabels_labelsize_ncolors_colormap(i, name)
            if vector_size == 3:
                plot_value = obj.get_plot_value(i, name) # vector
                scale = 1.0
                phase = 2.0
                obj.set_scale(i, name, scale)
                obj.set_phase(i, name, phase)
                assert obj.deflects(i, name) in [True, False], obj.deflects(i, name)
                xyz, deflected_xyz = obj.get_vector_result(i, name)
            else:
                scalar_result = obj.get_scalar(i, name)


            default_data_format = obj.get_default_data_format(i, name)
            default_min, default_max = obj.get_default_min_max(i, name)
            default_scale = obj.get_default_scale(i, name)
            default_title = obj.get_default_title(i, name)
            default_phase = obj.get_default_phase(i, name)
            out_labels = obj.get_default_nlabels_labelsize_ncolors_colormap(i, name)
            nlabels = 4
            labelsize = 10
            ncolors = 20
            colormap = 'jet'
            obj.set_nlabels_labelsize_ncolors_colormap(
                i, name, nlabels, labelsize, ncolors, colormap)
            default_nlabels, default_labelsize, default_ncolors, default_colormap = out_labels

            #default_max, default_min = obj.get_default_min_max(i, name)
            min_value, max_value = obj.get_min_max(i, name)

        self.result_cases = cases

        if len(self.case_keys) > 1:
            self.icase = -1
            self.ncases = len(self.result_cases)  # number of keys in dictionary
        elif len(self.case_keys) == 1:
            self.icase = -1
            self.ncases = 1
        else:
            self.icase = -1
            self.ncases = 0

    #def cycle_results(self, icase=None, show_msg=True):
        #"""fake method"""
        #pass

    #def cycle_results_explicit(self):
        #"""fake method"""
        #pass

    #def create_annotation(self, label, slot, x, y, z):
        #"""fake method"""
        #pass

    def update_axes_length(self, value):
        self.settings.dim_max = value

    #def passer(self):
        #"""fake method"""
        #pass

    #def passer1(self, a):
        #"""fake method"""
        #pass

    #def passer2(self, a, b):
        #"""fake method"""
        #pass

    @property
    def displacement_scale_factor(self):
        return 1 * self.settings.dim_max

    def _add_alt_actors(self, alt_grids):
        for name, grid in iteritems(alt_grids):
            self.geometry_actors[name] = vtkActor()

    #test.log_error = log_error
    #test.log_info = print
    #test.log_info = log_info
    #test.cycle_results = cycle_results
    #test.turn_text_on =  turn_text_on
    #test.turn_text_off = turn_text_off
    #test.cycle_results_explicit = passer

    def setWindowTitle(self, msg):
        assert isinstance(msg, str), 'msg=%r type=%r' % (msg, type(msg))
        return

    def get_edges(self):
        pass

    def getWindowTitle(self):  # pragma: no cover
        """fake QMainWindow method"""
        return 'title'

    def resize(self, height, width):  # pragma: no cover
        """fake QMainWindow method"""
        assert isinstance(height, int), 'height=%r' % height
        assert isinstance(width, int), 'width=%r' % width

    def setFont(self, font):
        """fake QMainWindow method"""
        pass

