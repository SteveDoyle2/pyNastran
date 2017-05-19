from __future__ import print_function
from collections import OrderedDict

from six import iteritems, integer_types
import numpy as np

import vtk
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
from pyNastran.utils.log import get_logger
from pyNastran.gui.qt_files.alt_geometry_storage import AltGeometry

#from pyNastran.gui.qt_version import qt_version
#if qt_version == 4:
#    from PyQt4.QtGui import QMainWindow
#elif qt_version == 5:
#    from PyQt5.QtWidgets import  QMainWindow

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
        self.follower_nodes = {}

        self.itext = 0

        self.pick_state = 'node/centroid'
        self.label_actors = {}
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
    def numpy_to_vtk_points(self, nodes, points=None, dtype='<f'):
        """
        common method to account for vtk endian quirks and
        efficiently adding points
        """
        assert isinstance(nodes, np.ndarray), type(nodes)
        if points is None:
            points = vtk.vtkPoints()
            nnodes = nodes.shape[0]
            points.SetNumberOfPoints(nnodes)

            # if we're in big endian, VTK won't work, so we byte swap
            nodes = np.asarray(nodes, dtype=np.dtype(dtype))

        points_array = numpy_to_vtk(
            num_array=nodes,
            deep=True,
            array_type=vtk.VTK_FLOAT,
        )
        points.SetData(points_array)
        return points

    def create_vtk_cells_of_constant_element_type(self, grid, elements, etype):
        """
        Adding constant type elements is overly complicated.  Note that
        some of the documentation in this method is likely triangle-specific
        as it was developed for a tri mesh.

        grid : vtk.vtkUnstructuredGrid()
            the unstructured grid
        elements : (nelements, nnodes_per_element) int ndarray
            the elements to add
        etype : int
            VTK cell type

            1  = vtk.vtkVertex().GetCellType()
            3  = vtkLine().GetCellType()
            5  = vtkTriangle().GetCellType()
            9  = vtk.vtkQuad().GetCellType()
            10 = vtkTetra().GetCellType()
            #vtkPenta().GetCellType()
            #vtkHexa().GetCellType()
            #vtkPyram().GetCellType()
        """
        nelements, nnodes_per_element = elements.shape
        # We were careful about how we defined the arrays, so the data
        # is contiguous when we ravel it.  Otherwise, you need to
        # deepcopy the arrays (deep=1).  However, numpy_to_vtk isn't so
        # good, so we could use np.copy, which is better, but it's
        # ultimately unnecessary.

        #nodes = numpy_to_vtk(elements, deep=0, array_type=vtk.VTK_ID_TYPE)
        # (nnodes_per_element + 1)  # TODO: was 4; for a tri...didn't seem to crash???
        cell_offsets = np.arange(0, nelements, dtype='int32') * (nnodes_per_element + 1)
        assert len(cell_offsets) == nelements

        # Create the array of cells
        vtk_cells = vtk.vtkCellArray()

        isize = vtk.vtkIdTypeArray().GetDataTypeSize()
        if isize == 4:
            dtype = 'int32' # TODO: can we include endian?
        elif isize == 8:
            dtype = 'int64'
        else:
            msg = 'elements.dtype=%s' % str(elements.dtype)
            raise NotImplementedError(msg)

        elements_vtk = np.zeros((nelements, nnodes_per_element + 1), dtype=dtype)
        elements_vtk[:, 0] = nnodes_per_element # 3 nodes/tri element
        elements_vtk[:, 1:] = elements

        cells_id_type = numpy_to_vtkIdTypeArray(elements_vtk.ravel(), deep=1)
        vtk_cells.SetCells(nelements, cells_id_type)

        # Cell types
        # 5 = vtkTriangle().GetCellType()
        cell_types = np.ones(nelements, dtype='int32') * etype
        vtk_cell_types = numpy_to_vtk(
            cell_types, deep=0,
            array_type=vtk.vtkUnsignedCharArray().GetDataType())

        vtk_cell_offsets = numpy_to_vtk(cell_offsets, deep=0,
                                        array_type=vtk.VTK_ID_TYPE)

        grid.SetCells(vtk_cell_types, vtk_cell_offsets, vtk_cells)

    def create_vtk_cells_of_constant_element_types(self, grid, elements_list,
                                                   etypes_list):
        """
        Adding constant type elements is overly complicated enough as in
        ``create_vtk_cells_of_constant_element_type``.  Now we extend
        this to multiple element types.

        grid : vtk.vtkUnstructuredGrid()
            the unstructured grid
        elements_list : List[elements, ...]
            elements : (nelements, nnodes_per_element) int ndarray
                the elements to add
        etypes_list : List[etype, ...]
            etype : int
                the VTK flag as defined in
                ``create_vtk_cells_of_constant_element_type``
        """
        if isinstance(etypes_list, list) and len(etypes_list) == 1:
            return self.create_vtk_cells_of_constant_element_type(
                grid, elements_list[0], etypes_list[0])

        isize = vtk.vtkIdTypeArray().GetDataTypeSize()
        if isize == 4:
            dtype = 'int32' # TODO: can we include endian?
        elif isize == 8:
            dtype = 'int64'
        else:
            msg = 'elements_list[0].dtype=%s' % str(elements_list[0].dtype)
            raise NotImplementedError(msg)

        cell_offsets_list2 = []
        cell_types_list2 = []
        elements_list2 = []
        nelements = 0
        noffsets = 0
        for element, etype in zip(elements_list, etypes_list):
            nelement, nnodes_per_element = element.shape

            nnodesp1 = nnodes_per_element + 1  # TODO: was 4; for a tri???
            cell_offset = np.arange(0, nelement, dtype='int32') * nnodesp1 + noffsets
            noffset = nelement * nnodesp1

            cell_type = np.ones(nelement, dtype='int32') * etype
            assert len(cell_offset) == nelement

            nnodesp1 = nnodes_per_element + 1
            element_vtk = np.zeros((nelement, nnodesp1), dtype=dtype)
            element_vtk[:, 0] = nnodes_per_element # 3 nodes/tri
            element_vtk[:, 1:] = element

            cell_offsets_list2.append(cell_offset)
            cell_types_list2.append(cell_type)
            elements_list2.append(element_vtk.ravel())
            nelements += nelement
            noffsets += noffset

        cell_types_array = np.hstack(cell_types_list2)
        cell_offsets_array = np.hstack(cell_offsets_list2)
        elements_array = np.hstack(elements_list2)

        # Create the array of cells
        cells_id_type = numpy_to_vtkIdTypeArray(elements_array.ravel(), deep=1)
        vtk_cells = vtk.vtkCellArray()
        vtk_cells.SetCells(nelements, cells_id_type)

        # Cell types
        vtk_cell_types = numpy_to_vtk(
            cell_types_array, deep=0,
            array_type=vtk.vtkUnsignedCharArray().GetDataType())

        vtk_cell_offsets = numpy_to_vtk(cell_offsets_array, deep=0,
                                        array_type=vtk.VTK_ID_TYPE)

        grid.SetCells(vtk_cell_types, vtk_cell_offsets, vtk_cells)


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

        points = self.numpy_to_vtk_points(nodes)
        grid = self.alt_grids[name]
        grid.SetPoints(points)

        etype = 9  # vtk.vtkQuad().GetCellType()
        self.create_vtk_cells_of_constant_element_type(grid, elements, etype)

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

        params_to_delete = (
            'case_keys', 'icase', 'iSubcaseNameMap',
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


class FakeGUIMethods(GuiAttributes):
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

        #print('self.case_keys = ', self.case_keys)
        for key in self.case_keys:
            assert isinstance(key, integer_types), key
            obj, (i, name) = cases[key]
            value = cases[key]
            assert not isinstance(value[0], int), 'key=%s\n type=%s value=%s' % (key, type(value[0]), value)
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

    def create_alternate_vtk_grid(self, name, color=None, line_width=None,
                                  opacity=None, point_size=None, bar_scale=None,
                                  representation=None, is_visible=True,
                                  follower_nodes=None, is_pickable=False):
        self.alt_grids[name] = Grid()
        geom = AltGeometry(self, name, color=color, line_width=line_width,
                           point_size=point_size, bar_scale=bar_scale,
                           opacity=opacity, representation=representation,
                           is_visible=is_visible, is_pickable=is_pickable)
        self.geometry_properties[name] = geom
        self.follower_nodes[name] = follower_nodes

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

