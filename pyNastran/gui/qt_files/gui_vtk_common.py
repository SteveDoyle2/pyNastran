# coding: utf-8
from collections import OrderedDict

import numpy as np
import vtk

from pyNastran.gui.qt_version import qt_version
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.gui.qt_files.gui_qt_common import GuiQtCommon
from pyNastran.gui.qt_files.mark_actions import create_annotation
from pyNastran.gui.utils.vtk.vtk_utils import map_element_centroid_to_node_fringe_result
from pyNastran.gui.menus.menus import Group


#class Interactor(vtk.vtkGenericRenderWindowInteractor):
    #def __init__(self):
        ##vtk.vtkGenericRenderWindowInteractor()
        #pass

    #def HighlightProp(self):
        #print('highlight')


#class PyNastranRenderWindowInteractor(QVTKRenderWindowInteractor):
    #def __init__(self, parent=None):

        #render_window = vtk.vtkRenderWindow()
        #iren = Interactor()
        #iren.SetRenderWindow(render_window)
        #kwargs = {
            #'iren' : iren,
            #'rw' : render_window,
        #}
        #QVTKRenderWindowInteractor.__init__(self, parent=parent,
                                            #iren=iren, rw=render_window)
        #self.Highlight

class GuiVTKCommon(GuiQtCommon):
    """this class has VTK functionality, but no interactive/menu capability"""
    def __init__(self, **kwds):
        if qt_version == 'pyqt5':
            super(GuiVTKCommon, self).__init__(**kwds)
        elif qt_version == 'pyside2':
            GuiQtCommon.__init__(self, **kwds)
        else:  #: pragma: no cover
            raise NotImplementedError(qt_version)

    def fake_init(self):
        """vtk setup"""
        self.create_vtk_actors(create_rend=False)
        #self._create_vtk_objects()

        #self.build_vtk_frame()
        #self.add_geometry()
        self._build_vtk_frame_post(build_lookup_table=False)
    #---------------------------------------------------------------------------
    # basic init
    def create_vtk_actors(self, create_rend=True):
        """creates the vtk actors used by the GUI"""
        if create_rend:
            self.rend = vtk.vtkRenderer()

        xyz = [0., 0., 0.]
        self.min_max_actors.append(create_annotation(self, 'Min', *xyz))
        self.min_max_actors.append(create_annotation(self, 'Max', *xyz))
        for actor in self.min_max_actors:
            actor.SetVisibility(False)

        # vtk actors
        self.grid = vtk.vtkUnstructuredGrid()

        # edges
        self.edge_actor = vtk.vtkLODActor()
        self.edge_actor.DragableOff()
        self.edge_mapper = vtk.vtkPolyDataMapper()

        self.create_cell_picker()

    def create_cell_picker(self):
        """creates the vtk picker objects"""
        self.cell_picker = vtk.vtkCellPicker()
        self.node_picker = vtk.vtkPointPicker()

        self.area_picker = vtk.vtkAreaPicker()  # vtkRenderedAreaPicker?
        self.rubber_band_style = vtk.vtkInteractorStyleRubberBandPick()
        #vtk.vtkInteractorStyleRubberBand2D
        #vtk.vtkInteractorStyleRubberBand3D
        #vtk.vtkInteractorStyleRubberBandZoom
        #vtk.vtkInteractorStyleAreaSelectHover
        #vtk.vtkInteractorStyleDrawPolygon

        #vtk.vtkAngleWidget
        #vtk.vtkAngleRepresentation2D
        #vtk.vtkAngleRepresentation3D
        #vtk.vtkAnnotation
        #vtk.vtkArrowSource
        #vtk.vtkGlyph2D
        #vtk.vtkGlyph3D
        #vtk.vtkHedgeHog
        #vtk.vtkLegendBoxActor
        #vtk.vtkLegendScaleActor
        #vtk.vtkLabelPlacer

        self.cell_picker.SetTolerance(0.001)
        self.node_picker.SetTolerance(0.001)

    def add_geometry(self):
        """
        #(N,)  for stress, x-disp
        #(N,3) for warp vectors/glyphs
        grid_result = vtk.vtkFloatArray()

        point_data = self.grid.GetPointData()
        cell_data = self.grid.GetCellData()

        self.grid.GetCellData().SetScalars(grid_result)
        self.grid.GetPointData().SetScalars(grid_result)


        self.grid_mapper   <-input-> self.grid
        vtkDataSetMapper() <-input-> vtkUnstructuredGrid()

        self.grid_mapper   <--map--> self.geom_actor <-add-> self.rend
        vtkDataSetMapper() <--map--> vtkActor()      <-add-> vtkRenderer()
        """
        if self.is_groups:
            # solid_bending: eids 1-182
            self._setup_element_mask()
            #eids = np.arange(172)
            #eids = arange(172)
            #self.update_element_mask(eids)
        else:
            self.grid_selected = self.grid
        #print('grid_selected =', self.grid_selected)

        self.grid_mapper = vtk.vtkDataSetMapper()
        self.grid_mapper.SetInputData(self.grid_selected)

        self._make_contour_filter()

        #if 0:
            #self.warp_filter = vtk.vtkWarpVector()
            #self.warp_filter.SetScaleFactor(50.0)
            #self.warp_filter.SetInput(self.grid_mapper.GetUnstructuredGridOutput())

            #self.geom_filter = vtk.vtkGeometryFilter()
            #self.geom_filter.SetInput(self.warp_filter.GetUnstructuredGridOutput())

            #self.geom_mapper = vtk.vtkPolyDataMapper()
            #self.geom_actor.setMapper(self.geom_mapper)

        #if 0:
            #from vtk.numpy_interface import algorithms
            #arrow = vtk.vtkArrowSource()
            #arrow.PickableOff()

            #self.glyph_transform = vtk.vtkTransform()
            #self.glyph_transform_filter = vtk.vtkTransformPolyDataFilter()
            #self.glyph_transform_filter.SetInputConnection(arrow.GetOutputPort())
            #self.glyph_transform_filter.SetTransform(self.glyph_transform)

            #self.glyph = vtk.vtkGlyph3D()
            #self.glyph.setInput(xxx)
            #self.glyph.SetSource(self.glyph_transform_filter.GetOutput())

            #self.glyph.SetVectorModeToUseVector()
            #self.glyph.SetColorModeToColorByVector()
            #self.glyph.SetScaleModeToScaleByVector()
            #self.glyph.SetScaleFactor(1.0)

            #self.append_filter = vtk.vtkAppendFilter()
            #self.append_filter.AddInputConnection(self.grid.GetOutput())


        #self.warpVector = vtk.vtkWarpVector()
        #self.warpVector.SetInput(self.grid_mapper.GetUnstructuredGridOutput())
        #grid_mapper.SetInput(Filter.GetOutput())

        self.geom_actor = vtk.vtkLODActor()
        self.geom_actor.DragableOff()
        self.geom_actor.SetMapper(self.grid_mapper)
        #geometryActor.AddPosition(2, 0, 2)
        #geometryActor.GetProperty().SetDiffuseColor(0, 0, 1) # blue
        #self.geom_actor.GetProperty().SetDiffuseColor(1, 0, 0)  # red

        #if 0:
            #id_filter = vtk.vtkIdFilter()

            #ids = np.array([1, 2, 3], dtype='int32')
            #id_array = numpy_to_vtk(
                #num_array=ids,
                #deep=True,
                #array_type=vtk.VTK_INT,
            #)

            #id_filter.SetCellIds(id_array.GetOutputPort())
            #id_filter.CellIdsOff()
            #self.grid_mapper.SetInputConnection(id_filter.GetOutputPort())
        self.rend.AddActor(self.geom_actor)
        self.build_glyph()

    def build_glyph(self):
        """builds the glyph actor"""
        glyph_source, glyphs, glyph_mapper, arrow_actor = build_glyph(self.grid)
        self.rend.AddActor(arrow_actor)

        self.glyph_source = glyph_source
        self.glyphs = glyphs
        self.glyph_mapper = glyph_mapper
        self.arrow_actor = arrow_actor
        #-----------------------------------------
        glyphs_centroid = vtk.vtkGlyph3D()
        glyphs_centroid.SetVectorModeToUseVector()
        glyphs_centroid.SetScaleModeToScaleByVector()
        glyphs_centroid.SetColorModeToColorByScale()
        glyphs_centroid.ScalingOn()
        glyphs_centroid.ClampingOn()
        glyphs_centroid.SetSourceConnection(glyph_source.GetOutputPort())

        glyph_mapper_centroid = vtk.vtkPolyDataMapper()
        glyph_mapper_centroid.SetInputConnection(glyphs_centroid.GetOutputPort())
        glyph_mapper_centroid.ScalarVisibilityOff()

        arrow_actor_centroid = vtk.vtkLODActor()
        arrow_actor_centroid.SetMapper(glyph_mapper_centroid)

        self.glyphs_centroid = glyphs_centroid
        self.glyph_mapper_centroid = glyph_mapper_centroid
        self.arrow_actor_centroid = arrow_actor_centroid

    def _build_vtk_frame_post(self, build_lookup_table=True):
        if build_lookup_table:
            self.build_lookup_table()

        text_size = self.settings.text_size # was 14
        dtext_size = text_size + 1
        self.create_text([5, 5 + 3 * dtext_size], 'Max  ', text_size)  # text actor 0
        self.create_text([5, 5 + 2 * dtext_size], 'Min  ', text_size)  # text actor 1
        self.create_text([5, 5 + 1 * dtext_size], 'Word1', text_size)  # text actor 2
        self.create_text([5, 5], 'Word2', text_size)  # text actor 3

        self.get_edges()
        if self.is_edges:
            prop = self.edge_actor.GetProperty()
            prop.EdgeVisibilityOn()
        else:
            prop = self.edge_actor.GetProperty()
            prop.EdgeVisibilityOff()

    def get_edges(self):
        """Create the edge actor"""
        edges = vtk.vtkExtractEdges()
        edge_mapper = self.edge_mapper
        edge_actor = self.edge_actor

        edges.SetInputData(self.grid_selected)
        edge_mapper.SetInputConnection(edges.GetOutputPort())

        edge_actor.SetMapper(edge_mapper)
        edge_actor.GetProperty().SetColor(0., 0., 0.)
        edge_mapper.SetLookupTable(self.color_function)
        edge_mapper.SetResolveCoincidentTopologyToPolygonOffset()

        prop = edge_actor.GetProperty()
        prop.SetColor(0., 0., 0.)
        edge_actor.SetVisibility(self.is_edges)
        self.rend.AddActor(edge_actor)

    #---------------------------------------------------------------------------
    # properties

    @property
    def logo(self):
        """Gets the pyNastran icon path, which can be overwritten"""
        return self._logo

    @logo.setter
    def logo(self, logo):
        """Sets the pyNastran icon path, which can be overwritten"""
        self._logo = logo

    @property
    def legend_shown(self):
        """determines if the legend is shown"""
        return self.scalar_bar.is_shown

    @property
    def scalar_bar_actor(self):
        """gets the scalar bar actor"""
        return self.scalar_bar.scalar_bar

    @property
    def color_function(self):
        """gets the scalar bar's color function"""
        return self.scalar_bar.color_function

    @property
    def result_name(self):
        """
        creates the self.result_name variable
        """
        # case_key = (1, 'ElementID', 1, 'centroid', '%.0f')
        case_key = self.case_keys[self.icase]
        assert isinstance(case_key, integer_types), case_key
        unused_obj, (unused_i, res_name) = self.result_cases[case_key]
        return res_name

    #---------------------------------------------------------------------------
    # basic interaction
    def on_show_debug(self):
        """sets a flag for showing/hiding DEBUG messages"""
        self.settings.show_debug = not self.settings.show_debug

    def on_show_info(self):
        """sets a flag for showing/hiding INFO messages"""
        self.settings.show_info = not self.settings.show_info

    def on_show_command(self):
        """sets a flag for showing/hiding COMMAND messages"""
        self.settings.show_command = not self.settings.show_command

    def on_show_warning(self):
        """sets a flag for showing/hiding WARNING messages"""
        self.settings.show_warning = not self.settings.show_warning

    def on_show_error(self):
        """sets a flag for showing/hiding ERROR messages"""
        self.settings.show_error = not self.settings.show_error

    def hide_axes(self, cids=None):
        """
        Show a set of coordinate systems

        ..todo :: fix the coords
        """
        if cids is None:
            cids = self.axes.keys()
        for cid in cids:
            axis = self.axes[cid]
            axis.VisibilityOff()
        self.corner_axis.EnabledOff()

    def show_axes(self, cids=None):
        """
        Show a set of coordinate systems

        ..todo :: fix the coords
        """
        if cids is None:
            cids = self.axes.keys()
        for cid in cids:
            axis = self.axes[cid]
            axis.VisibilityOn()
        self.corner_axis.EnabledOn()

    def delete_actor(self, name):
        """deletes an actor and associated properties"""
        if name != 'main':
            if name in self.geometry_actors:
                actor = self.geometry_actors[name]
                self.rend.RemoveActor(actor)
                del self.geometry_actors[name]

            if name in self.geometry_properties:
                unused_prop = self.geometry_properties[name]
                del self.geometry_properties[name]
            self.Render()

    def show_only(self, names):
        """
        Show these actors only

        names : str, List[str]
            names to show
            If they're hidden, show them.
            If they're shown and shouldn't be, hide them.

        ..todo :: update the GeomeryProperties
        """
        raise NotImplementedError('show_only')

    def hide_actors(self, except_names=None):
        """
        Hide all the actors

        except_names : str, List[str], None
            list of names to exclude
            None : hide all

        ..note :: If an actor is hidden and in the except_names, it will still be hidden.
        ..todo :: update the GeomeryProperties
        """
        if except_names is None:
            except_names = []
        elif isinstance(except_names, str):
            except_names = [except_names]

        # hide everything but the main grid
        for key, actor in self.geometry_actors.items():
            if key not in except_names:
                actor.VisibilityOff()

        self.hide_axes()
        self.hide_legend()
        #self.settings.set_background_color_to_white()

    #---------------------------------------------------------------------------
    def Render(self):
        """Renders the GUI"""
        #self.vtk_interactor.Render()
        self.vtk_interactor.GetRenderWindow().Render()

    def cell_centroid_grid(self, grid, cell_id, dtype='float32'):
        """gets the cell centroid"""
        if not self.is_gui:
            centroid = np.zeros(3, dtype=dtype)
            return centroid
        cell = grid.GetCell(cell_id)
        nnodes = cell.GetNumberOfPoints()
        points = cell.GetPoints()
        assert nnodes > 0, 'nnodes=%s cell_id=%s cell=%s' % (nnodes, cell_id, cell)
        centroid = np.zeros(3, dtype=dtype)
        for ipoint in range(nnodes):
            point = np.array(points.GetPoint(ipoint), dtype=dtype)
            centroid += point
        centroid /= nnodes
        return centroid

    def cell_centroid(self, cell_id, dtype='float32'):
        """gets the cell centroid"""
        return self.cell_centroid_grid(self.grid_selected, cell_id, dtype=dtype)

    def build_lookup_table(self):
        """build the nominal lookup table for the scalar bar"""
        scalar_range = self.grid_selected.GetScalarRange()
        self.grid_mapper.SetScalarRange(scalar_range)
        self.grid_mapper.SetLookupTable(self.color_function)
        self.rend.AddActor(self.scalar_bar_actor)

    def _set_results(self, form, cases):
        assert len(cases) > 0, cases
        if isinstance(cases, OrderedDict):
            self.case_keys = list(cases.keys())
        else:
            self.case_keys = sorted(cases.keys())
            assert isinstance(cases, dict), type(cases)

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
        self.icase_disp = None
        self.icase_vector = None
        self.icase_fringe = None
        self.set_form(form)

    #===========================================================================
    # groups
    # fake functions
    def show_eids(self, eids):
        pass

    # fake functions
    #---------------
    # real functions
    def find_result_by_name(self, desired_name):
        for icase in range(self.ncases):
            name, result = self.get_name_result_data(icase)
            if name == desired_name:
                return result
        raise RuntimeError('cannot find name=%r' % desired_name)

    def post_group_by_name(self, name):
        """posts a group with a specific name"""
        assert isinstance(name, str), name
        group = self.groups[name]
        self.post_group(group, update_groups=False)

    def post_group(self, group, update_groups=False):
        """posts a group object"""
        name = group.name
        if update_groups and name not in self.groups:
            self.groups[name] = group
        eids = group.element_ids
        self.show_eids(eids)
        assert isinstance(name, str), name
        self.group_active = name

    def create_group_with_name(self, name, eids, model_name='main'):
        """
        Creates a group from the root model (model_name)

        Parameters
        ----------
        name : str
            the name of the group
        eids : (neids, ) int ndarray
            the elements in the group
        model_name : str; default='main'
            the name of the parent model
        """
        elements_pound = self.groups[model_name].elements_pound
        element_str = ''
        group = Group(
            name, element_str, elements_pound,
            editable=True)

        # TODO: make sure all the eids exist
        group.element_ids = eids
        self.log_command('create_group_with_name(%r, %r)' % (name, eids))
        self.groups[name] = group

    def map_element_centroid_to_node_fringe_result(self, update_limits=True, show_msg=True):
        """
        Maps elemental fringe results to nodal fringe results.

        If you have a 5-noded CQUAD4 (e.g., centroid + 4 nodes), only the
        centroidal value will be mapped, even though you could map the
        average the nodal values instead.  It's not wrong to do it this
        way, but it could be more accurate.

        If you have CQUAD4 with centroidal only or something like strain
        energy, this will map properly.

        >>> is_passed = self.map_element_centroid_to_node_fringe_result(update_limits=True)
        """
        is_passed = False
        icase_fringe = self.icase_fringe
        if icase_fringe is None:
            self.log.error('No fringe result shown')
            return is_passed

        (obj, (i, name)) = self.result_cases[icase_fringe]
        location = obj.get_location(i, name)
        is_passed, out_data = map_element_centroid_to_node_fringe_result(
            self.grid, location, self.log)
        if not is_passed:
            return is_passed
        imin, imax, min_value_actual, max_value_actual = out_data

        title = obj.get_title(i, name)
        min_value, max_value = obj.get_min_max(i, name)
        if update_limits:
            min_value = min_value_actual
            max_value = max_value_actual

        data_format = obj.get_data_format(i, name)
        nlabels, labelsize, ncolors, colormap = obj.get_nlabels_labelsize_ncolors_colormap(i, name)
        is_legend_shown = self.scalar_bar.is_shown

        self.update_scalar_bar(title, min_value, max_value,
                               data_format,
                               nlabels=nlabels, labelsize=labelsize,
                               ncolors=ncolors, colormap=colormap,
                               is_shown=is_legend_shown)

        location_nodal = 'node'
        self._update_min_max_actors(location_nodal, icase_fringe,
                                    imin, min_value_actual,
                                    imax, max_value_actual)

        subcase_id = obj.subcase_id
        subtitle, label = self.get_subtitle_label(subcase_id)
        label2 = obj.get_header(i, name)
        if label2:
            label += '; ' + label2

        self.update_text_actors(subcase_id, subtitle,
                                imin, min_value_actual,
                                imax, max_value_actual, label, location_nodal)

        self.vtk_interactor.Render()
        self.log_command(f'map_element_centroid_to_node_fringe_result('
                         f'update_limits={update_limits}, show_msg={show_msg})')
        return is_passed



def build_glyph(grid):
    """builds the glyph actor"""
    glyphs = vtk.vtkGlyph3D()
    #if filter_small_forces:
        #glyphs.SetRange(0.5, 1.)

    glyphs.SetVectorModeToUseVector()
    #apply_color_to_glyph = False
    #if apply_color_to_glyph:
    #glyphs.SetScaleModeToScaleByScalar()
    glyphs.SetScaleModeToScaleByVector()
    glyphs.SetColorModeToColorByScale()
    #glyphs.SetColorModeToColorByScalar()  # super tiny
    #glyphs.SetColorModeToColorByVector()  # super tiny

    glyphs.ScalingOn()
    glyphs.ClampingOn()
    #glyphs.Update()

    glyph_source = vtk.vtkArrowSource()
    #glyph_source.InvertOn()  # flip this arrow direction
    glyphs.SetInputData(grid)


    glyphs.SetSourceConnection(glyph_source.GetOutputPort())
    #glyphs.SetScaleModeToDataScalingOff()
    #glyphs.SetScaleFactor(10.0)  # bwb
    #glyphs.SetScaleFactor(1.0)  # solid-bending
    glyph_mapper = vtk.vtkPolyDataMapper()
    glyph_mapper.SetInputConnection(glyphs.GetOutputPort())
    glyph_mapper.ScalarVisibilityOff()

    arrow_actor = vtk.vtkLODActor()
    arrow_actor.SetMapper(glyph_mapper)

    prop = arrow_actor.GetProperty()
    prop.SetColor(1., 0., 0.)
    #self.grid.GetPointData().SetActiveVectors(None)
    arrow_actor.SetVisibility(False)
    return glyph_source, glyphs, glyph_mapper, arrow_actor
