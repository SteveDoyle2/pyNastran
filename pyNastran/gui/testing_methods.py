from typing import Any
import numpy as np
from cpylog import get_logger

from pyNastran.gui.qt_version import qt_version
from vtk import (
    vtkTextActor, vtkLODActor, vtkActor,
    #GeometryProperty,
    #GridMapper,
    #Grid,
    vtkArrowSource,
    vtkGlyph3D,
    vtkPolyDataMapper,
)

from pyNastran.gui.gui_objects.settings import Settings
from pyNastran.gui.vtk_interface import vtkUnstructuredGrid
from pyNastran.gui.test.mock_vtk import (
    #GeometryProperty,
    GridMapper,
    VTKInteractor, vtkRenderer,
)
from pyNastran.gui.menus.groups_modify.groups_modify import Group

#from pyNastran.gui.test.mock_vtk_bkp import (
    #Grid,
    #vtkActor,
    #vtkArrowSource,
    #vtkGlyph3D,
    #vtkPolyDataMapper,
#)
import vtkmodules

from pyNastran.gui.gui_common import GuiVTKCommon
from pyNastran.gui.qt_files.scalar_bar import ScalarBar
#from pyNastran.gui.gui_objects.alt_geometry_storage import AltGeometry
from pyNastran.gui.formats import CLASS_MAP

#from pyNastran.bdf.cards.base_card import deprecated

#class ScalarBar:
    #def VisibilityOff(self):
        #pass
    #def VisibilityOn(self):
        #pass
    #def Modified(self):
        #pass
    #@property
    #def is_shown(self):
        #return True


class Button:
    def __init__(self):
        pass
    def setChecked(self, is_checked):
        pass

class MockTreeView:
    def __init__(self):
        self.fringe = Button()
        self.vector = Button()
        self.disp = Button()

class MockResultCaseWindow:
    def __init__(self):
        self.tree_view = MockTreeView()

class MockResWidget:
    def __init__(self):
        self.result_case_window = MockResultCaseWindow()
    def update_results(self, form, name):
        """fake method"""
        pass
    def update_method(self, methods):
        """fake method"""
        pass
    def update_icase(self, icase):
        pass
    def set_methods_table_visible(self, is_visible: bool) -> None:
        assert isinstance(is_visible, bool), is_visible
    def set_coord_transform_visible(self, is_visible: bool, methods: list[str]) -> None:
        assert isinstance(is_visible, bool), is_visible
        assert isinstance(methods, list), methods
        for method in methods:
            assert isinstance(method, str), method
    def set_derivation_visible(self, is_visible: bool, methods_dict: dict[str, Any]) -> None:
        assert isinstance(is_visible, bool), is_visible
        assert isinstance(methods_dict, dict), methods_dict
        unused_tooltip = methods_dict.get('tooltip', '')
        unused_label =  methods_dict.get('label', '')
        default = []
        methods =  methods_dict.get('derivation', default)
        for method in methods:
            assert isinstance(method, str), method

    def set_nodal_combine_visible(self, is_visible: bool, methods: list[str]) -> None:
        assert isinstance(is_visible, bool), is_visible
        assert isinstance(methods, list), methods
        for method in methods:
            assert isinstance(method, str), method
    def set_output_checkbox(self,
                            is_enabled_fringe: bool, is_checked_fringe: bool,
                            is_enabled_disp: bool, is_checked_disp: bool,
                            is_enabled_vector: bool, is_checked_vector: bool) -> None:
        pass

class FakeGUIMethods(GuiVTKCommon):
    """all the methods in here are faked"""
    def __init__(self, inputs=None):
        if inputs is None:
            inputs = {
                'magnify' : 1,
                'debug' : False,
                'console' : True,
                'is_groups' : True,
            }
        self.rend = vtkRenderer()
        GuiVTKCommon.__init__(self, inputs=inputs)
        self.fake_init()

        # the gui is not actually running
        self.is_gui = False


        res_widget = MockResWidget()
        kwds = {
            'inputs' : inputs,
            'res_widget' : res_widget
        }
        #GuiAttributes.__init__(self, **kwds)
        #GuiVTKCommon.__init__(self, **kwds)
        self.res_widget = res_widget
        self.vtk_interactor = VTKInteractor()
        self.debug = False
        self._form = []
        #self.geometry_actors = {
            #'main' : vtkActor(),
        #}
        settings = Settings(self)
        self.scalar_bar = ScalarBar(settings)
        self.grid = vtkUnstructuredGrid()
        self.main_grid_mappers = {'main' : GridMapper()}
        if 1:  # pragma: no cover
            #self.scalar_bar_actor = ScalarBar()
            self.alt_geometry_actor = ScalarBar(self)
            self.alt_grids = {
                'main' : self.grid,
            }
            self.main_geometry_actors = {
                'main' :  vtkActor(),
            }

            self.glyph_source = vtkArrowSource()
            self.glyphs = vtkGlyph3D()
            self.glyph_mapper = vtkPolyDataMapper()
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

        self.corner_text_actors[0] = vtkTextActor()
        self.corner_text_actors[1] = vtkTextActor()
        self.corner_text_actors[2] = vtkTextActor()
        self.corner_text_actors[3] = vtkTextActor()
        self.format_class_map = CLASS_MAP

    def cell_centroid(self, cell_id, dtype='float32'):
        return np.zeros(3, dtype=dtype)

    def setup_fake_text_actors(self):
        for icase in self.result_cases:
            self.label_actors[icase] = []

    #@property
    #def scalar_bar_actor(self):
        #return self.scalar_bar.scalar_bar

    @property
    def grid_selected(self):
        return self.grid

    #def hide_legend(self):
        #pass
    #def show_legend(self):
        #pass

    #def update_scalar_bar(self, title, min_value, max_value,
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

    def _finish_results_io2(self, model_name: str, form, cases, reset_labels: bool=True):
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
        assert isinstance(self.nnodes, int), 'nnodes=%r must be an integer' % self.nnodes
        assert isinstance(self.nelements, int), 'nelements=%r must be an integer' % self.nelements

        assert len(cases) > 0, cases
        self.case_keys = list(cases.keys())
        #self.case_keys = sorted(cases.keys())
        assert isinstance(cases, dict), type(cases)

        #print('self.case_keys = ', self.case_keys)
        for key in self.case_keys:
            assert isinstance(key, int), key
            obj, (i, name) = cases[key]
            value = cases[key]
            if isinstance(value[0], int):
                raise RuntimeError('old style key is being used.\n key=%s\n type=%s value=%s' % (
                    key, type(value[0]), value))
            #assert len(value) == 2, 'value=%s; len=%s' % (str(value), len(value))

            unused_subcase_id = obj.subcase_id
            unused_methods = obj.get_methods(i, name)
            #method = methods[0]
            unused_fringe = obj.get_fringe_result(i, name)
            unused_fringe, unused_disp = obj.get_fringe_vector_result(i, name)
            unused_legend_title = obj.get_legend_title(i, name)
            vector_size = obj.get_vector_size(i, name)
            #location = obj.get_location(i, name)
            unused_data_format = obj.get_data_format(i, name)
            scale = obj.get_scale(i, name)
            #arrow_scale = obj.get_arrow_scale(i, name)
            phase = obj.get_phase(i, name)
            unused_label2 = obj.get_annotation(i, name)
            unused_flag = obj.is_normal_result(i, name)
            #scalar_result = obj.get_scalar(i, name)
            is_method_array = obj.is_method_array
            outi = obj.get_nlabels_labelsize_ncolors_colormap(i, name)
            nlabels, labelsize, ncolors, colormap = outi
            if colormap == '':
                colormap = 'jet'

            if vector_size == 3:
                #unused_plot_value = obj.get_plot_value(i, name) # vector
                scale = 1.1
                phase = 2.2
                obj.set_scale(i, name, scale)
                #obj.set_arrow_scale(i, name, scale)
                obj.set_phase(i, name, phase)
                assert obj.deflects(i, name) in {True, False}, obj.deflects(i, name)
                unused_xyz, unused_deflected_xyz = obj.get_vector_result(i, name)
            #else:
                #unused_scalar_result = obj.get_scalar(i, name)

            unused_default_data_format = obj.get_default_data_format(i, name)
            default_min, unused_default_max = obj.get_default_min_max(i, name)
            unused_default_scale = obj.get_default_scale(i, name)
            #unused_default_arrow_scale = obj.get_default_arrow_scale(i, name)
            unused_default_title = obj.get_default_legend_title(i, name)
            unused_default_phase = obj.get_default_phase(i, name)
            out_labels = obj.get_default_nlabels_labelsize_ncolors_colormap(i, name)
            nlabels = 4
            labelsize = 10
            ncolors = 20
            obj.set_nlabels_labelsize_ncolors_colormap(
                i, name, nlabels, labelsize, ncolors, colormap)
            (unused_default_nlabels, unused_default_labelsize,
             unused_default_ncolors, unused_default_colormap) = out_labels

            unused_min_value, unused_max_value = obj.get_min_max(i, name)
            unused_imin_value, unused_imax_value = obj.get_imin_imax(i, name)

        self.model_data.result_cases = cases

        if len(self.case_keys) > 1:
            self.icase = -1
            self.ncases = len(self.result_cases)  # number of keys in dictionary
        elif len(self.case_keys) == 1:
            self.icase = -1
            self.ncases = 1
        else:
            self.icase = -1
            self.ncases = 0

        if self.is_groups and len(self.element_ids):
            #eids = np.arange(172)
            #eids = []
            #self.hide_elements_mask(eids)
            elements_pound = self.element_ids[-1]
            main_group = Group(
                'main', '', elements_pound,
                editable=False)
            main_group.element_ids = self.element_ids
            self.groups['main'] = main_group
            self.post_group(main_group)
            #self.show_elements_mask(np.arange(self.nelements))

        for unused_module_name, module in self.modules.items():
            module.post_load_geometry()

    #def cycle_results(self, icase=None, show_msg=True):
        #"""fake method"""
        #pass

    #def cycle_results_explicit(self):
        #"""fake method"""
        #pass

    #def create_annotation(self, label, x, y, z):
        #"""fake method"""
        #return None

    #def update_axes_length(self, value):
        #self.settings.dim_max = value

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
        for name, unused_grid in alt_grids.items():
            self.geometry_actors[name] = vtkActor()

    #test.log_error = log_error
    #test.log_info = print
    #test.log_info = log_info
    #test.cycle_results = cycle_results
    #test.turn_corner_text_on =  turn_corner_text_on
    #test.turn_corner_text_off = turn_corner_text_off
    #test.cycle_results_explicit = passer

    def setWindowTitle(self, msg):
        assert isinstance(msg, str), 'msg=%r type=%r' % (msg, type(msg))
        return

    def get_edges(self):
        pass

    def getWindowTitle(self):  # pragma: no cover
        """fake QMainWindow method"""
        return 'title'

    def resize(self, height: int, width: int):  # pragma: no cover
        """fake QMainWindow method"""
        assert isinstance(height, int), 'height=%r' % height
        assert isinstance(width, int), 'width=%r' % width

    def setFont(self, font):
        """fake QMainWindow method"""
        pass
