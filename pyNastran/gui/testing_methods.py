from __future__ import print_function
import os
from collections import OrderedDict

from six import iteritems, integer_types
import numpy as np

import vtk
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
from pyNastran.utils.log import get_logger
from pyNastran.gui.qt_files.alt_geometry_storage import AltGeometry
#from pyNastran.gui.gui_objects.gui_result import GuiResult
#from pyNastran.converters.nastran.displacements import DisplacementResults
from pyNastran.bdf.cards.base_card import deprecated
from pyNastran.gui.gui_utils.utils import create_res_obj
from pyNastran.gui.gui_utils.vtk_utils import (
    create_vtk_cells_of_constant_element_type, numpy_to_vtk_points)


#class GuiAttributes(QMainWindow):
class GuiAttributes(object):
    """All methods in this class must not require VTK"""
    def __init__(self, **kwds):
        """
        These variables are common between the GUI and
        the batch mode testing that fakes the GUI
        """
        inputs = kwds['inputs']
        res_widget = kwds['res_widget']
        self.dev = False

        # the result type being currently shown
        # for a Nastran NodeID/displacement, this is 'node'
        # for a Nastran ElementID/PropertyID, this is 'element'
        self.result_location = None

        self.case_keys = {}
        self.res_widget = res_widget
        self._show_flag = True
        self._camera_mode = None
        self.observers = {}

        if 'test' in inputs:
            self.is_testing_flag = inputs['test']
        else:
            self.is_testing_flag = False
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
        self._legend_window_shown = False
        self._preferences_window_shown = False
        self._clipping_window_shown = False
        self._edit_geometry_properties_window_shown = False
        self._modify_groups_window_shown = False
        #self._label_window = None
        #-------------
        # inputs dict
        self.is_edges = False
        self.is_edges_black = self.is_edges
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
        self.nnodes = 0
        self.nelements = 0

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
        self.follower_nodes = {}

        self.itext = 0

        self.pick_state = 'node/centroid'
        self.label_actors = {-1 : []}
        self.label_ids = {}
        self.cameras = {}
        self.label_scale = 1.0 # in percent

        self.is_horizontal_scalar_bar = False

        self.result_cases = {}
        self.num_user_points = 0

        self._is_displaced = False
        self._is_forces = False
        self._is_normals = False

        self._xyz_nominal = None

        self.nvalues = 9
        self.dim_max = 1.0
        self.nid_maps = {}
        self.eid_maps = {}
        self.name = 'main'

        self.groups = {}
        self.group_active = 'main'

        #if not isinstance(res_widget, MockResWidget):
            #if qt_version == 4:
                #QMainWindow.__init__(self)
            #elif qt_version == 5:
                #super(QMainWindow, self).__init__()

        self.main_grids = {}
        self.main_grid_mappers = {}
        self.main_geometry_actors = {}

        self.main_edge_mappers = {}
        self.main_edge_actors = {}

    #-------------------------------------------------------------------
    # deprecated attributes
    def deprecated(self, old_name, new_name, deprecated_version):
        # type: (str, str, str, Optional[List[int]]) -> None
        """
        Throws a deprecation message and crashes if past a specific version.

        Parameters
        ----------
        old_name : str
            the old function name
        new_name : str
            the new function name
        deprecated_version : float
            the version the method was first deprecated in
        """
        deprecated(old_name, new_name, deprecated_version, levels=[0])

    @property
    def nNodes(self):
        """gets the number of nodes"""
        self.deprecated('self.nNodes', 'self.nnodes', '1.1')
        return self.nnodes

    @nNodes.setter
    def nNodes(self, nnodes):
        """sets the number of nodes"""
        self.deprecated('self.nNodes', 'self.nnodes', '1.1')
        self.nnodes = nnodes

    @property
    def nElements(self):
        """gets the number of elements"""
        self.deprecated('self.nElements', 'self.nelements', '1.1')
        return self.nelements

    @nElements.setter
    def nElements(self, nelements):
        """sets the number of elements"""
        self.deprecated('self.nElements', 'self.nelements', '1.1')
        self.nelements = nelements

    #-------------------------------------------------------------------
    # geom
    @property
    def grid(self):
        #print('get grid; %r' % self.name)
        return self.main_grids[self.name]

    @grid.setter
    def grid(self, grid):
        #print('set grid; %r' % self.name)
        self.main_grids[self.name] = grid

    @property
    def grid_mapper(self):
        return self.main_grid_mappers[self.name]

    @grid_mapper.setter
    def grid_mapper(self, grid_mapper):
        self.main_grid_mappers[self.name] = grid_mapper

    @property
    def geom_actor(self):
        return self.main_geometry_actors[self.name]

    @geom_actor.setter
    def geom_actor(self, geom_actor):
        self.main_geometry_actors[self.name] = geom_actor

    #-------------------------------------------------------------------
    # edges
    @property
    def edge_mapper(self):
        return self.main_edge_mappers[self.name]

    @edge_mapper.setter
    def edge_mapper(self, edge_mapper):
        self.main_edge_mappers[self.name] = edge_mapper

    @property
    def edge_actor(self):
        return self.main_edge_actors[self.name]

    @edge_actor.setter
    def edge_actor(self, edge_actor):
        self.main_edge_actors[self.name] = edge_actor

    def set_glyph_scale_factor(self, scale):
        """sets the glyph scale factor"""
        self.glyph_scale_factor = scale
        self.glyphs.SetScaleFactor(scale)

    #-------------------------------------------------------------------
    def _load_patran_nod(self, nod_filename):
        """reads a Patran formatted *.nod file"""
        from pyNastran.bdf.patran_utils.read_patran import read_patran
        data_dict = read_patran(nod_filename, fdtype='float32', idtype='int32')
        nids = data_dict['nids']
        data = data_dict['data']
        data_headers = data_dict['headers']
        ndata = data.shape[0]
        if len(data.shape) == 1:
            shape = (ndata, 1)
            data = data.reshape(shape)

        if ndata != self.node_ids.shape[0]:
            inids = np.searchsorted(self.node_ids, nids)
            assert np.array_equal(nids, self.node_ids[inids]), 'the node ids are invalid'
            data2 = np.full(data.shape, np.nan, data.dtype)
            data2[inids, :] = data
        else:
            data2 = data

        A = {}
        fmt_dict = {}
        headers = data_headers['SEC']
        for i, header in enumerate(headers):
            A[header] = data2[:, i]
            fmt_dict[header] = '%f'

        out_filename_short = os.path.relpath(nod_filename)
        result_type = 'node'
        self._add_cases_to_form(A, fmt_dict, headers, result_type,
                                out_filename_short, update=True,
                                is_scalar=True)

    def _add_cases_to_form(self, A, fmt_dict, headers, result_type,
                           out_filename_short, update=True, is_scalar=True,
                           is_deflection=False, is_force=False):
        """
        common method between:
          - _load_csv
          - _load_deflection_csv

        Parameters
        ----------
        A : dict[key] = (n, m) array
            the numpy arrays
            key : str
                the name
            n : int
                number of nodes/elements
            m : int
                secondary dimension
                N/A : 1D array
                3 : deflection
        fmt_dict : dict[header] = fmt
            the format of the arrays
            header : str
                the name
            fmt : str
                '%i', '%f'
        headers : List[str]???
            the titles???
        result_type : str
            'node', 'centroid'
        out_filename_short : str
            the display name
        update : bool; default=True
            ???

        # A = np.loadtxt('loadtxt_spike.txt', dtype=('float,int'))
        # dtype=[('f0', '<f8'), ('f1', '<i4')])
        # A['f0']
        # A['f1']
        """
        #print('A =', A)
        formi = []
        form = self.get_form()
        icase = len(self.case_keys)
        islot = 0
        for case_key in self.case_keys:
            if isinstance(case_key, tuple):
                islot = case_key[0]
                break

            if is_scalar:
                out = create_res_obj(islot, headers, A, fmt_dict, result_type)
            else:
                out = create_res_obj(islot, headers, A, fmt_dict, result_type,
                                     self.dim_max, self.xyz_cid0)
            res_obj, title, header = out

            #cases[icase] = (stress_res, (subcase_id, 'Stress - isElementOn'))
            #form_dict[(key, itime)].append(('Stress - IsElementOn', icase, []))
            #key = (res_obj, (0, title))
            self.case_keys.append(icase)
            self.result_cases[icase] = (res_obj, (islot, title))
            formi.append((header, icase, []))

            # TODO: double check this should be a string instead of an int
            self.label_actors[icase] = []
            self.label_ids[icase] = set([])
            icase += 1
        form.append((out_filename_short, None, formi))

        self.ncases += len(headers)
        #cases[(ID, 2, 'Region', 1, 'centroid', '%i')] = regions
        self.res_widget.update_results(form, 'main')


    #-------------------------------------------------------------------
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
        assert isinstance(nodes, np.ndarray), type(nodes)

        points = numpy_to_vtk_points(nodes)
        grid = self.alt_grids[name]
        grid.SetPoints(points)

        etype = 9  # vtk.vtkQuad().GetCellType()
        create_vtk_cells_of_constant_element_type(grid, elements, etype)

        self._add_alt_actors({name : self.alt_grids[name]})

        #if name in self.geometry_actors:
        self.geometry_actors[name].Modified()

    def create_coordinate_system(self, dim_max, label='', origin=None, matrix_3x3=None,
                                 Type='xyz'):
        """
        Creates a coordinate system

        Parameters
        ----------
        dim_max : float
            the max model dimension; 10% of the max will be used for
            the coord length
        label : str
            the coord id or other unique label (default is empty to
            indicate the global frame)
        origin : (3, ) ndarray/list/tuple
            the origin
        matrix_3x3 : (3, 3) ndarray
            a standard Nastran-style coordinate system
        Type : str
            a string of 'xyz', 'Rtz', 'Rtp' (xyz, cylindrical, spherical)
            that changes the axis names
        """
        pass

    @property
    def nid_map(self):
        return self.nid_maps[self.name]

    @nid_map.setter
    def nid_map(self, nid_map):
        self.nid_maps[self.name] = nid_map

    @property
    def eid_map(self):
        try:
            return self.eid_maps[self.name]
        except:
            msg = 'KeyError: key=%r; keys=%s' % (self.name, list(self.eid_maps.keys()))
            raise KeyError(msg)

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
        """
        Sets the path to the icon directory where custom icons are found
        """
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
            assert isinstance(key, int), key
            obj, (i, name) = self.result_cases[key]
            t = (i, [])
            data.append(t)

        self.res_widget.update_results(formi, self.name)

        key = list(self.case_keys)[0]
        location = self.get_case_location(key)
        method = 'centroid' if location else 'nodal'

        data2 = [(method, None, [])]
        self.res_widget.update_methods(data2)

    def _remove_old_geometry(self, geom_filename):
        skip_reading = False
        if self.dev:
            return skip_reading

        self.eid_map = {}
        self.nid_map = {}
        params_to_delete = (
            'case_keys', 'icase', 'isubcase_name_map',
            'result_cases', 'eid_map', 'nid_map',
        )
        if geom_filename is None or geom_filename is '':
            skip_reading = True
            return skip_reading
        else:
            self.turn_text_off()
            self.grid.Reset()

            self.result_cases = {}
            self.ncases = 0
            for param in params_to_delete:
                if hasattr(self, param):  # TODO: is this correct???
                    try:
                        delattr(self, param)
                    except AttributeError:
                        msg = 'cannot delete %r; hasattr=%r' % (param, hasattr(self, param))
                        self.log.warning(msg)

            skip_reading = False
        #self.scalarBar.VisibilityOff()
        self.scalarBar.Modified()
        return skip_reading


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
    def SetCells(self, vtk_cell_types, vtk_cell_locations, vtk_cells):
        pass


class ScalarBar(object):
    def VisibilityOff(self):
        pass
    def VisibilityOn(self):
        pass
    def Modified(self):
        pass


class ArrowSource(object):
    def __init__(self):
        pass
class Glyph3D(object):
    def SetScaleFactor(self, value):
        pass
class PolyDataMapper(object):
    def __init__(self):
        pass
class LODActor(object):
    def __init__(self):
        pass


class MockResWidget(object):
    def __init__(self):
        pass
    def update_results(self, form, name):
        """fake method"""
        pass

class FakeGUIMethods(GuiAttributes):
    """all the methods in here are faked"""
    def __init__(self, inputs=None):
        if inputs is None:
            inputs = {
                'magnify' : 1,
                'debug' : False,
                'console' : True,
            }

        res_widget = MockResWidget()
        kwds = {
            'inputs' : inputs,
            'res_widget' : res_widget
        }
        GuiAttributes.__init__(self, **kwds)
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

        self.glyph_source = ArrowSource()
        self.glyphs = Glyph3D()
        self.glyph_mapper = PolyDataMapper()
        self.arrow_actor = LODActor()

        #self.geometry_properties = {
            #'main' : None,
            #'caero' : None,
            #'caero_sub' : None,
        #}
        #self._add_alt_actors = _add_alt_actors

        level = 'debug' if self.debug else 'info'
        self.log = get_logger(log=None, level=level)

    def _finish_results_io2(self, form, cases, reset_labels=True):
        """
        This is not quite the same as the main one.
        It's more or less just _set_results
        """
        #assert self.node_ids is not None
        #assert self.element_ids is not None

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
            result_type = obj.get_title(i, name)
            vector_size = obj.get_vector_size(i, name)
            #location = obj.get_location(i, name)
            methods = obj.get_methods(i)
            data_format = obj.get_data_format(i, name)
            scale = obj.get_scale(i, name)
            phase = obj.get_phase(i, name)
            label2 = obj.get_header(i, name)
            flag = obj.is_normal_result(i, name)
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

    def cycle_results(self):
        """fake method"""
        pass

    def cycle_results_explicit(self):
        """fake method"""
        pass

    def _create_annotation(self, label, icase, x, y, z):
        """fake method"""
        pass

    def  turn_text_on(self):
        """fake method"""
        pass

    def turn_text_off(self):
        """fake method"""
        pass

    def create_global_axes(self, dim_max):
        """fake method"""
        pass

    def update_axes_length(self, value):
        self.dim_max = value

    def passer(self):
        """fake method"""
        pass

    def passer1(self, a):
        """fake method"""
        pass

    def passer2(self, a, b):
        """fake method"""
        pass

    @property
    def displacement_scale_factor(self):
        return 1 * self.dim_max

    def create_alternate_vtk_grid(self, name, color=None, line_width=None,
                                  opacity=None, point_size=None, bar_scale=None,
                                  representation=None, is_visible=True,
                                  follower_nodes=None, is_pickable=False):
        """Fake creates an AltGeometry object"""
        self.alt_grids[name] = Grid()
        geom = AltGeometry(self, name, color=color, line_width=line_width,
                           point_size=point_size, bar_scale=bar_scale,
                           opacity=opacity, representation=representation,
                           is_visible=is_visible, is_pickable=is_pickable)
        self.geometry_properties[name] = geom
        if follower_nodes is not None:
            self.follower_nodes[name] = follower_nodes

    def duplicate_alternate_vtk_grid(self, name, name_duplicate_from, color=None, line_width=5,
                                     opacity=1.0, point_size=1, bar_scale=0.0, is_visible=True,
                                     follower_nodes=None, is_pickable=False):
        """Fake copies the VTK actor"""
        if name_duplicate_from == 'main':
            grid_copy_from = self.grid
            representation = 'toggle'
        else:
            grid_copy_from = self.alt_grids[name_duplicate_from]
            props = self.geometry_properties[name_duplicate_from]
            representation = props.representation

        self.alt_grids[name] = Grid()
        geom = AltGeometry(self, name, color=color, line_width=line_width,
                           point_size=point_size, bar_scale=bar_scale,
                           opacity=opacity, representation=representation,
                           is_visible=is_visible, is_pickable=is_pickable)
        self.geometry_properties[name] = geom
        if follower_nodes is not None:
            self.follower_nodes[name] = follower_nodes

    def _add_alt_actors(self, alt_grids):
        for name, grid in iteritems(alt_grids):
            self.geometry_actors[name] = GeometryActor()

    def log_debug(self, msg):
        """turns logs into prints to aide testing debug"""
        if self.debug:
            print('DEBUG:  ', msg)

    def log_info(self, msg):
        """turns logs into prints to aide testing debug"""
        if self.debug:
            print('INFO:  ', msg)

    def log_error(self, msg):
        """turns logs into prints to aide testing debug"""
        if self.debug:
            print('ERROR:  ', msg)

    def log_warning(self, msg):
        """turns logs into prints to aide testing debug"""
        if self.debug:
            print('WARNING:  ', msg)

    #test.log_error = log_error
    #test.log_info = print
    #test.log_info = log_info
    #test.cycle_results = cycle_results
    #test.turn_text_on =  turn_text_on
    #test.turn_text_off = turn_text_off
    #test.cycle_results_explicit = passer

