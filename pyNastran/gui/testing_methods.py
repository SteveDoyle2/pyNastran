from __future__ import print_function
from collections import OrderedDict

from six import iteritems

from pyNastran.utils.log import get_logger
from pyNastran.gui.qt_files.alt_geometry_storage import AltGeometry


class GuiAttributes(object):
    """All methods in this class must not require VTK"""
    def __init__(self, inputs, res_widget):
        """
        These variables are common between the GUI and
        the batch mode testing that fakes the GUI
        """
        self.case_keys = {}
        self.res_widget = res_widget
        self._show_flag = True

        self.is_testing = False
        self.is_groups = False
        self._logo = None
        self._script_path = None
        self._icon_path = ''

        self.title = None
        self.min_value = None
        self.max_value = None
        self.blue_to_red = False
        self._is_axes_shown = True
        self.nvalues = 9
        self.is_wireframe = False
        #-------------

        # window variables
        self._picker_window_shown = False
        self._legend_window_shown = False
        self._clipping_window_shown = False
        self._edit_geometry_properties_window_shown = False
        self._modify_groups_window_shown = False
        self._label_window_shown = False
        #self._label_window = None
        #-------------
        # inputs dict
        self.is_edges = False
        self.is_edges_black = self.is_edges
        # self.is_nodal = inputs['is_nodal']
        # self.is_centroidal = inputs['is_centroidal']
        self.magnify = inputs['magnify']

        #self.format = ''
        debug = inputs['debug']
        self.debug = debug
        assert debug in [True, False], 'debug=%s' % debug

        #-------------
        # file
        self.menu_bar_format = None
        self.format = None
        self.infile_name = None
        self.out_filename = None
        self.dirname = ''
        self.last_dir = '' # last visited directory while opening file
        self._default_python_file = None

        #-------------
        # internal params
        self.show_info = True
        self.show_debug = True
        self.show_gui = True
        self.show_command = True
        self.coord_id = 0

        self.ncases = 0
        self.icase = 0
        self.nNodes = 0
        self.nElements = 0

        self.supported_formats = []
        self.model_type = None

        self.tools = []
        self.checkables = []
        self.actions = {}

        # actor_slots
        self.text_actors = {}
        self.geometry_actors = OrderedDict()
        self.alt_grids = {} #additional grids

        # coords
        self.transform = {}
        self.axes = {}

        #geom = Geom(color, line_thickness, etc.)
        #self.geometry_properties = {
        #    'name' : Geom(),
        #}
        self.geometry_properties = OrderedDict()

        self.itext = 0

        self.pick_state = 'node/centroid' # if self.is_centroidal else 'nodal'
        self.label_actors = {}
        self.label_ids = {}
        self.cameras = {}
        self.label_scale = 1.0 # in percent

        self.is_horizontal_scalar_bar = False

        self.result_cases = {}
        self.num_user_points = 0

        self._is_displaced = False
        self._xyz_nominal = None

        self.nvalues = 9
        self.dim_max = 1.0
        self.nid_maps = {}
        self.eid_maps = {}
        self.name = 'main'

        self.groups = {}
        self.group_active = 'main'

    def set_quad_grid(self, name, nodes, elements, color, line_width=5, opacity=1.):
        """
        Makes a CQUAD4 grid
        """
        self.create_alternate_vtk_grid(name, color=color, line_width=line_width,
                                       opacity=opacity, representation='wire')

        nnodes = nodes.shape[0]
        nquads = elements.shape[0]
        #print(nodes)
        if nnodes == 0:
            return
        if nquads == 0:
            return

        #print('adding quad_grid %s; nnodes=%s nquads=%s' % (name, nnodes, nquads))
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(nnodes)
        for nid, node in enumerate(nodes):
            #print(nid, node)
            points.InsertPoint(nid, *list(node))

        #assert vtkQuad().GetCellType() == 9, elem.GetCellType()
        self.alt_grids[name].Allocate(nquads, 1000)
        for element in elements:
            elem = vtkQuad()
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, element[0])
            point_ids.SetId(1, element[1])
            point_ids.SetId(2, element[2])
            point_ids.SetId(3, element[3])
            self.alt_grids[name].InsertNextCell(9, elem.GetPointIds())
        self.alt_grids[name].SetPoints(points)

        self._add_alt_actors({name : self.alt_grids[name]})

        #if name in self.geometry_actors:
        self.geometry_actors[name].Modified()

    def create_coordinate_system(self, dim_max, label='', origin=None, matrix_3x3=None, Type='xyz'):
        pass

    @property
    def nid_map(self):
        return self.nid_maps[self.name]

    @nid_map.setter
    def nid_map(self, nid_map):
        self.nid_maps[self.name] = nid_map

    @property
    def eid_map(self):
        return self.eid_maps[self.name]

    @eid_map.setter
    def eid_map(self, eid_map):
        self.eid_maps[self.name] = eid_map

    @property
    def displacement_scale_factor(self):
        """
        # dim_max = max_val * scale
        # scale = dim_max / max_val
        # 0.25 added just cause

        scale = self.displacement_scale_factor / tnorm_abs_max
        """
        #scale = self.dim_max / tnorm_abs_max * 0.25
        scale = self.dim_max * 0.25
        return scale

    def set_script_path(self, script_path):
        """Sets the path to the custom script directory"""
        self._script_path = script_path

    def set_icon_path(self, icon_path):
        """Sets the path to the icon directory where custom icons are found"""
        self._icon_path = icon_path

    def form(self):
        formi = self.res_widget.get_form()
        return formi

    def get_form(self):
        return self._form

    def set_form(self, formi):
        self._form = formi
        data = []
        for key in self.case_keys:
            #print(key)
            if isinstance(key, int):
                obj, (i, name) = self.result_cases[key]
                t = (i, [])
            else:
                t = (key[1], [])
            data.append(t)

        self.res_widget.update_results(formi)

        key = list(self.case_keys)[0]
        location = self.get_case_location(key)
        method = 'centroid' if location else 'nodal'

        data2 = [(method, None, [])]
        self.res_widget.update_methods(data2)


class CoordProperties(object):
    def __init__(self, label, Type, is_visible, scale):
        self.label = label
        #self.axes = axes
        self.Type = Type
        self.is_visible = is_visible
        self.representation = 'coord'
        self.scale = scale


class GeometryProperty(object):
    def __init__(self):
        pass
    def SetRepresentationToPoints(self):
        pass
    def SetPointSize(self, size):
        assert isinstance(size, int), type(size)


class GeometryActor(object):
    def __init__(self):
        self._prop = GeometryProperty()
    def GetProperty(self):
        return self._prop
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


class MockResWidget(object):
    def __init__(self):
        pass


class GUIMethods(GuiAttributes):
    def __init__(self, inputs=None):
        if inputs is None:
            inputs = {
                'magnify' : 1,
                'debug' : False,
                'console' : True,
            }

        res_widget = MockResWidget()
        GuiAttributes.__init__(self, inputs, res_widget)
        self.is_testing = True
        self.debug = False
        self._form = []
        self.result_cases = {}
        self._finish_results_io = self.passer1
        #self.geometry_actors = {
            #'main' : GeometryActor(),
        #}
        self.grid = Grid()
        self.scalarBar = ScalarBar()
        self.alt_geometry_actor = ScalarBar()
        self.alt_grids = {
            'main' : self.grid,
        }
        #self.geometry_properties = {
            ##'main' : None,
            ##'caero' : None,
            ##'caero_sub' : None,
        #}
        #self._add_alt_actors = _add_alt_actors

        level = 'debug' if self.debug else 'info'
        self.log = get_logger(log=None, level=level)

    def _finish_results_io2(self, form, cases):
        """
        This is not quite the same as the main one.
        It's more or less just _set_results
        """
        assert len(cases) > 0, cases
        if isinstance(cases, OrderedDict):
            self.case_keys = list(cases.keys())
        else:
            self.case_keys = sorted(cases.keys())
            assert isinstance(cases, dict), type(cases)

        for key in self.case_keys:
            value = cases[key]
            #print('value[0] = %s' % value[0])
            assert not isinstance(value[0], int), 'key=%s\n type=%s value=%s' % (key, type(value[0]), value)

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

    def _remove_old_geometry(self, filename):
        skip_reading = False
        return skip_reading

    def cycle_results(self):
        pass

    def  turn_text_on(self):
        pass

    def turn_text_off(self):
        pass

    def create_global_axes(self, dim_max):
        pass

    def update_axes_length(self, value):
        self.dim_max = value

    def passer(self):
        pass

    def passer1(self, a):
        pass

    def passer2(self, a, b):
        pass

    @property
    def displacement_scale_factor(self):
        return 1 * self.dim_max

    def create_alternate_vtk_grid(self, name, color=None, line_width=None, opacity=None,
                                  point_size=None, bar_scale=None,
                                  representation=None, is_visible=True):
        self.alt_grids[name] = Grid()
        geom = AltGeometry(self, name, color=color, line_width=line_width,
                           point_size=point_size, bar_scale=bar_scale,
                           opacity=opacity, representation=representation,
                           is_visible=is_visible)
        self.geometry_properties[name] = geom

    def _add_alt_actors(self, alt_grids):
        for name, grid in iteritems(alt_grids):
            self.geometry_actors[name] = GeometryActor()

    def log_debug(self, msg):
        if self.debug:
            print('DEBUG:  ', msg)

    def log_info(self, msg):
        if self.debug:
            print('INFO:  ', msg)

    def log_error(self, msg):
        if self.debug:
            print('ERROR:  ', msg)

    def log_warning(self, msg):
        if self.debug:
            print('WARNING:  ', msg)

    #test.log_error = log_error
    #test.log_info = print
    #test.log_info = log_info
    #test.cycle_results = cycle_results
    #test.turn_text_on =  turn_text_on
    #test.turn_text_off = turn_text_off
    #test.cycle_results_explicit = passer

